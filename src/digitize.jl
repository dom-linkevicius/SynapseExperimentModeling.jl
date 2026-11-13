# Digitizing data from published figure images; a port of the previous project's Python
# digitizer (scikit-image) to Julia.
#
# A figure is described by its plot panels - where each plot frame sits in the image and
# how its tick marks map pixels to data - and by the series drawn in them, told apart by
# ink colour. Series come in two kinds:
#   * MarkerSeries: separate markers, found by template matching (normalised
#     cross-correlation);
#   * TraceSeries: densely sampled data whose markers merge into a band, read column by
#     column as the centre of the band.
#
# Pixel coordinates are 0-based, x to the right and y down, as image viewers show them.

using FileIO: load, save
using ImageCore: RGB, N0f8, red, green, blue
using ImageFiltering: imfilter, centered, mapwindow, Fill, Algorithm
using ImageMorphology: feature_transform, distance_transform, label_components, strel_diamond
using DataStructures: BinaryMinHeap
using CairoMakie: Makie, Figure, Axis, Rect, Circle, image!, scatter!, lines!, text!, hidedecorations!, hidespines!,
                  resize_to_layout!
import XLSX

"""
    AxisMap(ticks; log = false)

Linear map between pixel position and data value along one axis, least-squares fitted to
`ticks`, a list of `pixel => value` pairs read off the figure's tick marks. On a log axis
the map is linear in `log10(value)`.
"""
struct AxisMap
    offset::Float64     # pixel = offset + slope * (log ? log10(value) : value)
    slope::Float64
    log::Bool
end

function AxisMap(ticks; log::Bool = false)
    px = Float64[first(t) for t in ticks]
    v = Float64[log ? log10(last(t)) : last(t) for t in ticks]
    pm, vm = sum(px) / length(px), sum(v) / length(v)
    slope = sum((v .- vm) .* (px .- pm)) / sum(abs2, v .- vm)
    return AxisMap(pm - slope * vm, slope, log)
end

"""
    pixel_to_value(m::AxisMap, px) → value
Data value at pixel position `px` along the axis mapped by `m`, rounded to the precision a
tenth of a pixel resolves.
"""
function pixel_to_value(m::AxisMap, px)
    t = (px - m.offset) / m.slope
    m.log && return round(exp10(t); sigdigits = 4)
    return round(t; digits = ceil(Int, -log10(0.1 / abs(m.slope))))
end

"""
    PlotPanel(; name, frame, x, y, skip = [], title = name, crop, grid, conditions = (;),
              xvar = nothing, yvar = nothing, margin = 3, fontsize = nothing, axis = (;),
              frame_width = nothing, read = nothing)

One plot box of a figure. `frame` is `(left, top, right, bottom)` of the plot frame (or axis)
lines, `x` and `y` are the `AxisMap`s from its tick marks, and `skip` lists boxes, in the
same `(left, top, right, bottom)` form, whose ink is not data (legends, annotations, an
inset). Data are searched for inside the frame lines, `margin` pixels in from them, or
inside `read` (same form) where only part of the plot is to be read, e.g. the time range of
one trace; the digitized plot still spans the whole frame, as the original does. `crop`
is the region of the image shown above the digitized plot in the calibration figure, and
`grid` the panel's `(row, column)` there. `conditions` are the experimental conditions
shared by the panel's series, e.g. `(CaN_peptide = true,)`. `xvar`, `yvar` (`quantity =>
unit`), `fontsize` and `axis` (Makie `Axis` attributes) are for panels whose quantities,
text size or axes differ from the rest of the figure's, such as an inset. `frame_width` is
the width of the digitized frame in the calibration figure (default `FRAME_WIDTH`), e.g. to
draw the segments of a broken axis at one scale.
"""
Base.@kwdef struct PlotPanel
    name::String
    frame::NTuple{4,Int}
    x::AxisMap
    y::AxisMap
    skip::Vector{NTuple{4,Int}} = NTuple{4,Int}[]
    title::String = name
    crop::NTuple{4,Int}
    grid::Tuple{Int,Int}
    conditions::NamedTuple = (;)
    xvar::Union{Nothing,Pair{String,String}} = nothing
    yvar::Union{Nothing,Pair{String,String}} = nothing
    margin::Int = 3
    fontsize::Union{Nothing,Float64} = nothing
    axis::NamedTuple = (;)
    frame_width::Union{Nothing,Float64} = nothing
    read::Union{Nothing,NTuple{4,Int}} = nothing
end

"""
    MarkerSeries(; name, panel, ink, shape = :circle, radius, hollow = false, stroke = 2, threshold = 0.45,
                 locate = :centroid, near = [], fixed = [], conditions = (;), colour, label = name, saved = true)

A series of separate markers in panel `panel` with experimental `conditions` of its own
(on top of the panel's), found in the colour test `ink` (a key of
`INK`) by a template of `shape` (`:circle`, `:square`, `:diamond`, `:triangle` or
`:triangle_down`) with half-width `radius` pixels; open markers (`hollow`) by its outline,
`stroke` pixels wide. Template-match peaks below
`threshold` times the best match in the panel are dropped. A marker's centre is the
centroid of its watershed basin (`locate = :centroid`), or the template-match peak itself
(`:peak`) where markers sit on lines that would pull a centroid off. In figures too small
or crowded for markers to stand out from lines, `near` lists their approximate centres
(pixels) instead, each snapped to the best template match within `radius`; the series then
competes with others in its ink only at those markers. `fixed` lists centres read by eye
and kept as given, for markers no template can match, e.g. ones on a frame line. Error bars
are not read. `colour` and `label` are how the series is drawn in the digitized plot. A
series with `saved = false` only competes for markers, so that those of the others are told
apart, and is neither written nor drawn.
"""
Base.@kwdef struct MarkerSeries
    name::String
    panel::String
    ink::Symbol
    shape::Symbol = :circle
    radius::Int
    hollow::Bool = false
    stroke::Int = 2
    threshold::Float64 = 0.45
    locate::Symbol = :centroid
    near::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    fixed::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    conditions::NamedTuple = (;)
    colour::String
    label::String = name
    saved::Bool = true
end

"""
    TraceSeries(; name, panel, ink, curve = 1, hollow = false, gap = 1, conditions = (;), shape = :circle,
                radius, colour, label = "", label_below = false, saved = true, span = nothing)

A densely sampled series in panel `panel`, with experimental `conditions` of its own, whose
markers merge into a band of ink `ink`, read as one point per pixel column: the centre of
the band. These points are not the individual measurements. Several traces in one ink are
told apart by `curve`, their order from the top (1) where all bands are separate. A trace
that covers only part of the panel's x range, e.g. one that ends early, is given that range
as `span` (x data values, `(from, to)`); outside it the trace is not tracked. Open
markers (`hollow`) are filled in first, and gaps of up to `gap` pixels within a band are
bridged. `shape`, `radius`, `colour` and `label` are how the trace is drawn in the
digitized plot; the label goes above the trace's right end, or below it (`label_below`).
A trace with `saved = false` is only tracked, so that the bands of the others are told
apart, and is neither written nor drawn.
"""
Base.@kwdef struct TraceSeries
    name::String
    panel::String
    ink::Symbol
    curve::Int = 1
    hollow::Bool = false
    gap::Int = 1
    conditions::NamedTuple = (;)
    shape::Symbol = :circle
    radius::Int
    colour::String
    label::String = ""
    label_below::Bool = false
    saved::Bool = true
    span::Union{Nothing,NTuple{2,Float64}} = nothing
end

const Series = Union{MarkerSeries,TraceSeries}

"""
    PaperFigure(; paper, image, doi, description, panels, series, xvar, yvar, fontsize = 40,
                linewidth = 5, ticksize = 14, shared_ylabel = false, axis = (;))

A published figure to digitize, the image `papers/<paper>/<image>`. Its data go to
`data/exp_raw/<paper>.xlsx`, in the tab named after the image (`Fig10.jpg` → `FIG10`), with
the same columns for every figure: `panel`, `series`, `x`, `y` and `merged` (where traces
ran together). What they mean goes to the workbook's `META` tab, one row per panel and
series: the figure's `doi` and `description`, how the series was read, the x and y
quantities and units (`xvar`, `yvar` as `quantity => unit`, unless the panel sets its
own), and the panel's and series' conditions, one column each. `load_paper_data` joins the
two back together.

In the calibration figure, `fontsize` and `ticksize` are the size of the tick labels and
marks drawn over the original, in image pixels, and the `xticks` and `yticks` of `axis`
(Makie `Axis` attributes, as in the original) are the ticks drawn. `linewidth`,
`shared_ylabel` and the rest of `axis` are no longer used.
"""
Base.@kwdef struct PaperFigure
    paper::String
    image::String
    doi::String
    description::String
    panels::Vector{PlotPanel}
    series::Vector{Series}
    xvar::Union{Nothing,Pair{String,String}} = nothing
    yvar::Union{Nothing,Pair{String,String}} = nothing
    fontsize::Float64 = 40
    linewidth::Float64 = 5
    ticksize::Float64 = 14
    shared_ylabel::Bool = false
    axis::NamedTuple = (;)
end

figure_stem(fig::PaperFigure) = first(splitext(fig.image))
figure_tab(fig::PaperFigure) = uppercase(figure_stem(fig))
paper_dir(paper) = pkgdir(@__MODULE__, "papers", paper)
exp_raw_path(paper) = pkgdir(@__MODULE__, "data", "exp_raw", paper * ".xlsx")

# Colour tests on 0-255 RGB channel matrices, each returning the mask of that ink.
# Thresholds come from sampling marker pixels of the figures that use them.
const INK = Dict{Symbol,Function}(
    :blue => (r, g, b) -> (b .- r .> 45) .& (b .> 85) .& (r .< 160) .& (b .- g .> 30),
    :orange => (r, g, b) -> (r .- b .> 90) .& (r .> 180) .& (r .- g .> 40),
    # thin black lines and open markers come out mid-grey in small JPEGs
    :black => (r, g, b) -> (max.(r, g, b) .< 150) .& (max.(r, g, b) .- min.(r, g, b) .< 60),
    # anything from black to light grey, e.g. grey markers under a black fit line
    :grey => (r, g, b) -> (max.(r, g, b) .< 215) .& (max.(r, g, b) .- min.(r, g, b) .< 40),
)

function solid_template(shape::Symbol, r::Integer)
    # a square needs a ring of background around it: an all-ones template has no variance,
    # so its normalised cross-correlation is undefined
    shape === :square && return [float(max(abs(i), abs(j)) <= r) for i in -(r + 1):(r + 1), j in -(r + 1):(r + 1)]
    shape === :circle && return [float(i^2 + j^2 <= r^2) for i in -r:r, j in -r:r]
    shape === :diamond && return [float(abs(i) + abs(j) <= r) for i in -(r + 1):(r + 1), j in -(r + 1):(r + 1)]
    # equilateral, apex up: half-base r, height √3 r, widening linearly from the apex
    if shape === :triangle
        h = round(Int, sqrt(3) * r)
        apex = -(h ÷ 2)
        return [float(apex <= i < apex + h && abs(j) <= r * (i - apex) / max(h - 1, 1))
                for i in -(r + 1):(r + 1), j in -(r + 1):(r + 1)]
    end
    shape === :triangle_down && return reverse(solid_template(:triangle, r); dims = 1)
    error("unknown marker shape $shape; use :circle, :square, :diamond, :triangle or :triangle_down")
end

# An open marker is matched by its outline: error bars, caps and fit lines that run through
# it don't look like an outline, whereas once filled they would match a solid shape. The
# outline is `stroke` px wide: 2 by default, to tolerate the half-pixel misfits of thin
# JPEG strokes, wider for markers drawn with thick strokes.
function marker_template(shape::Symbol, r::Integer; hollow::Bool = false, stroke::Integer = 2)
    t = solid_template(shape, r)
    hollow || return t
    # pad so the shape has background all round, then keep what is within `stroke` px of it
    p = zeros(size(t) .+ 2)
    p[2:(end - 1), 2:(end - 1)] .= t
    edge = distance_transform(feature_transform(p .== 0)) .<= stroke
    return (p .* edge)[2:(end - 1), 2:(end - 1)]
end

"""
    match_template(img, template) → Matrix
Normalised cross-correlation of `template` centred on every pixel of `img`, zero-padded
at the edges: scikit-image's `match_template(img, template, pad_input=True)`.
"""
function match_template(img::AbstractMatrix, template::AbstractMatrix)
    n = length(template)
    t0 = template .- sum(template) / n
    box = centered(ones(size(template)))
    # FIR keeps the window sums of a 0/1 image exact, so flat windows give exactly zero
    # variance instead of FFT round-off
    num = imfilter(img, centered(t0), Fill(0), Algorithm.FIR())
    s1 = imfilter(img, box, Fill(0), Algorithm.FIR())
    s2 = imfilter(img .^ 2, box, Fill(0), Algorithm.FIR())
    den = sqrt.(max.((s2 .- s1 .^ 2 ./ n) .* sum(abs2, t0), 0))
    return [d > eps() ? c / d : 0.0 for (c, d) in zip(num, den)]
end

"""
    peak_local_max(score; min_distance, threshold_rel) → Vector{CartesianIndex{2}}
Local maxima of `score` within a `(2min_distance + 1)²` window that exceed
`threshold_rel * maximum(score)`, thinned so no two are closer than `min_distance`
(Chebyshev), strongest first: scikit-image's `peak_local_max(..., exclude_border=False)`.
"""
function peak_local_max(score::AbstractMatrix; min_distance::Integer, threshold_rel::Real)
    w = 2min_distance + 1
    localmax = mapwindow(maximum, score, (w, w))
    cut = max(minimum(score), threshold_rel * maximum(score))
    peaks = [I for I in CartesianIndices(score) if score[I] == localmax[I] && score[I] > cut]
    # ties in row-major order, as scikit-image sorts them
    sort!(peaks; by = I -> (-score[I], I[1], I[2]))
    kept = CartesianIndex{2}[]
    for I in peaks
        all(J -> maximum(abs, Tuple(I - J)) >= min_distance, kept) && push!(kept, I)
    end
    return kept
end

"""
    watershed(img, seeds, mask) → Matrix{Int}
Flood `img` from the labelled `seeds` (0 = none), lowest values first, through the
4-connected pixels of `mask`: scikit-image's `watershed(img, seeds, mask=mask)`, which
also drops seeds outside the mask.
"""
function watershed(img::AbstractMatrix, seeds::AbstractMatrix{<:Integer}, mask::AbstractMatrix{Bool})
    labels = seeds .* mask
    lin = LinearIndices(labels)
    cart = CartesianIndices(labels)
    # (value, age, linear index): equal values flood first-come first-served
    heap = BinaryMinHeap{Tuple{Float64,Int,Int}}()
    for I in cart
        labels[I] > 0 && push!(heap, (img[I], 0, lin[I]))
    end
    age = 0
    neighbours = (CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(0, 1))
    while !isempty(heap)
        I = cart[last(pop!(heap))]
        for d in neighbours
            J = I + d
            (checkbounds(Bool, labels, J) && mask[J] && labels[J] == 0) || continue
            labels[J] = labels[I]
            age += 1
            push!(heap, (img[J], age, lin[J]))
        end
    end
    return labels
end

"""
    fill_holes(mask) → BitMatrix
`mask` with every enclosed (4-connected) background region filled, so open markers become
solid: scipy's `binary_fill_holes`.
"""
function fill_holes(mask::AbstractMatrix{Bool})
    labels = label_components(.!mask, strel_diamond((3, 3)))
    h, w = size(labels)
    outside = Set(labels[I] for I in CartesianIndices(labels) if I[1] in (1, h) || I[2] in (1, w))
    return mask .| [l > 0 && l ∉ outside for l in labels]
end

# 1-based rows and columns of the image searched for data: inside the frame lines
function search_region(p::PlotPanel)
    left, top, right, bottom = something(p.read, p.frame)
    return (top + p.margin + 1):(bottom - p.margin + 1), (left + p.margin + 1):(right - p.margin + 1)
end

# template-match score of `s`'s marker over `panel`'s search region of `mask`
function marker_score(mask, panel::PlotPanel, s::MarkerSeries)
    rows, cols = search_region(panel)
    crop = mask[rows, cols]
    return match_template(Float64.(crop), marker_template(s.shape, s.radius; s.hollow, s.stroke)), crop
end

# `s`'s template score as a rival of other series in its ink; a series whose markers are
# given (`near`) competes only within `radius` of them
function rival_score(mask, panel::PlotPanel, s::MarkerSeries)
    score = first(marker_score(mask, panel, s))
    isempty(s.near) && return score
    rows, cols = search_region(panel)
    keep = falses(size(score))
    for (x, y) in s.near
        i, j = round(Int, y) - rows[1] + 2, round(Int, x) - cols[1] + 2
        keep[max(1, i - s.radius):min(end, i + s.radius), max(1, j - s.radius):min(end, j + s.radius)] .= true
    end
    return score .* keep
end

# centroids (row, col) of the watershed basins of `peaks` in `crop` that have a plausible
# marker area; watershed splits blobs where several markers (or a marker plus a fit line)
# fused, and seeding it with the template peaks keeps one basin per marker
function basin_centroids(crop, peaks, template)
    seeds = zeros(Int, size(crop))
    for (k, I) in enumerate(peaks)
        seeds[I] = k
    end
    distance = distance_transform(feature_transform(.!crop))
    basins = watershed(-distance, seeds, crop)
    n = length(peaks)
    area, sum_row, sum_col = zeros(Int, n), zeros(n), zeros(n)
    for I in CartesianIndices(basins)
        k = basins[I]
        k == 0 && continue
        area[k] += 1
        sum_row[k] += I[1]
        sum_col[k] += I[2]
    end
    return [(sum_row[k] / area[k], sum_col[k] / area[k]) for k in 1:n
            if 0.25 * sum(template) <= area[k] <= 6 * sum(template)]
end

# template-match peak refined to sub-pixel by a parabola through it and its neighbours
function refine_peak(score, I)
    function shift(a, b, c)
        den = a - 2b + c
        return den < 0 ? clamp((a - c) / 2den, -0.5, 0.5) : 0.0
    end
    i, j = Tuple(I)
    di = 1 < i < size(score, 1) ? shift(score[i - 1, j], score[i, j], score[i + 1, j]) : 0.0
    dj = 1 < j < size(score, 2) ? shift(score[i, j - 1], score[i, j], score[i, j + 1]) : 0.0
    return (i + di, j + dj)
end

# best score within `r` pixels of `I`: different templates peak at slightly different
# offsets on the same marker, so rivals are compared over the marker, not at one pixel
local_max(score, I, r) = maximum(view(score, max(1, I[1] - r):min(size(score, 1), I[1] + r),
                                      max(1, I[2] - r):min(size(score, 2), I[2] + r)))

"""
    find_markers(mask, panel, series; rivals = ()) → Vector{NTuple{2,Float64}}
Pixel centres `(x, y)` of the markers of `series` in the ink `mask`, inside `panel`'s frame.
`rivals` are the template scores of other series in the same ink and panel: a marker goes
to whichever series' template matches it best within its radius.
"""
function find_markers(mask::AbstractMatrix{Bool}, panel::PlotPanel, s::MarkerSeries; rivals = ())
    rows, cols = search_region(panel)
    score, crop = marker_score(mask, panel, s)
    template = marker_template(s.shape, s.radius; s.hollow, s.stroke)

    if isempty(s.near)
        # Marker shapes are degraded where they overlap, so the cut is relative to the best
        # match in this panel rather than an absolute correlation.
        peaks = peak_local_max(score; min_distance = max(2, s.radius - 1), threshold_rel = s.threshold)
        filter!(I -> all(rival -> score[I] >= local_max(rival, I, s.radius), rivals), peaks)
    else
        # the best match within `radius` of each given approximate centre
        peaks = map(s.near) do (x, y)
            i, j = round(Int, y) - rows[1] + 2, round(Int, x) - cols[1] + 2
            window = CartesianIndices((max(1, i - s.radius):min(size(score, 1), i + s.radius),
                                       max(1, j - s.radius):min(size(score, 2), j + s.radius)))
            window[argmax(score[window])]
        end
    end
    centres = s.locate === :peak || !isempty(s.near) ? [refine_peak(score, I) for I in peaks] :
              basin_centroids(crop, peaks, template)

    pts = NTuple{2,Float64}[]
    for (i, j) in centres
        # from 1-based crop indices back to 0-based image pixels
        x, y = cols[1] - 2 + j, rows[1] - 2 + i
        any(((l, t, r, b),) -> l <= x <= r && t <= y <= b, panel.skip) && continue
        push!(pts, (x, y))
    end
    return append!(pts, s.fixed)
end

"""
    ink_runs(column, gap) → Vector{UnitRange{Int}}
Runs of `true` in `column`, joining runs separated by at most `gap` false pixels.
"""
function ink_runs(column, gap)
    runs = UnitRange{Int}[]
    for i in findall(column)
        if !isempty(runs) && i - last(runs[end]) <= gap + 1
            runs[end] = first(runs[end]):i
        else
            push!(runs, i:i)
        end
    end
    return runs
end

# Largest jump, in pixels, a trace may make from one column to the next; ink farther from
# a trace than this is not part of it.
const MAX_JUMP = 25

"""
    trace_bands(mask, panel, series) → Dict{String, Vector{NTuple{3}}}
`(x, y, merged)` for every pixel column of each `TraceSeries` in `series` (all in the same
ink and panel): the centre of its band in that column. Tracking starts at a column where
every band is separate - the `curve`-th band from the top belongs to the series with that
`curve` - and runs outward column by column, each series taking the band nearest its last
position, also where a column shows as many bands as series. Where bands have merged, each series sits its usual thickness in from its own
edge of the merged band, and is flagged `merged`. A series is tracked only in the columns
within its `span`.
"""
function trace_bands(mask::AbstractMatrix{Bool}, panel::PlotPanel, series)
    rows, cols = search_region(panel)
    crop = mask[rows, cols]
    # annotations and insets are in skip boxes; blank them so they can't pass for a band
    for (l, t, r, b) in panel.skip
        crop[max(1, t + 2 - rows[1]):min(end, b + 2 - rows[1]), max(1, l + 2 - cols[1]):min(end, r + 2 - cols[1])] .= false
    end
    any(s -> s.hollow, series) && (crop = fill_holes(crop))
    series = sort(series; by = s -> s.curve)
    n = length(series)

    # bands in every column, ignoring single-pixel specks
    bands = [filter(r -> length(r) > 1, ink_runs(view(crop, :, j), maximum(s -> s.gap, series)))
             for j in axes(crop, 2)]
    # series present in each column: those whose span covers it
    within(s, j) = isnothing(s.span) || s.span[1] <= pixel_to_value(panel.x, cols[1] - 2 + j) <= s.span[2]
    present = [[k for k in 1:n if within(series[k], j)] for j in axes(crop, 2)]
    resolved = findall(j -> length(present[j]) == n && length(bands[j]) == n, axes(crop, 2))
    isempty(resolved) && error("panel $(panel.name): no column shows all $n bands apart")
    thickness = map(1:n) do k
        sort([length(bands[j][k]) for j in resolved])[cld(length(resolved), 2)]
    end

    mid(r) = (first(r) + last(r)) / 2
    distance(r, y) = y < first(r) ? first(r) - y : y > last(r) ? y - last(r) : 0.0
    pos, merged = fill(NaN, n, size(crop, 2)), falses(n, size(crop, 2))
    function place!(j, prev)
        b, ks = bands[j], present[j]
        # bands go to the series in order only at the start; after it, each series follows its nearest band
        if !isempty(ks) && length(b) == length(ks) && all(k -> isnan(prev[k]), ks)
            pos[ks, j] .= mid.(b)
            return [isnan(pos[k, j]) ? prev[k] : pos[k, j] for k in 1:n]
        end
        (isempty(b) || isempty(ks)) && return prev
        owner = Dict(k => argmin(r -> distance(r, prev[k]), b) for k in ks)
        # a series farther than MAX_JUMP from every band (e.g. one that has ended) owns none
        reach = Dict(k => distance(owner[k], prev[k]) <= MAX_JUMP for k in ks)
        for k in ks
            r = owner[k]
            reach[k] || continue
            shared = filter(i -> reach[i] && owner[i] == r, ks)
            if length(shared) == 1
                pos[k, j] = mid(r)
            else
                merged[k, j] = true
                pos[k, j] = k == first(shared) ? first(r) + (thickness[k] - 1) / 2 :
                            k == last(shared) ? last(r) - (thickness[k] - 1) / 2 : mid(r)
            end
        end
        return [isnan(pos[k, j]) ? prev[k] : pos[k, j] for k in 1:n]
    end
    seed = resolved[cld(length(resolved), 2)]
    start = place!(seed, fill(NaN, n))
    prev = start
    for j in (seed + 1):size(crop, 2)
        prev = place!(j, prev)
    end
    prev = start
    for j in (seed - 1):-1:1
        prev = place!(j, prev)
    end

    # from 1-based crop indices back to 0-based image pixels
    return Dict(s.name => [(cols[1] - 2 + j, rows[1] - 2 + pos[k, j], merged[k, j])
                           for j in axes(crop, 2) if !isnan(pos[k, j])]
                for (k, s) in enumerate(series))
end

"""
    digitize(fig::PaperFigure) → DataFrame
Read every series of `fig`, write the data to the figure's tab of
`data/exp_raw/<paper>.xlsx` and what they are to its `META` tab, and save the calibration
figure - the digitized data drawn over the original (see `calibration_figure`) - as
`papers/<paper>/<image stem>_digitized.png`.
"""
function digitize(fig::PaperFigure)
    dir = paper_dir(fig.paper)
    img = RGB{N0f8}.(load(joinpath(dir, fig.image)))
    r, g, b = (Int.(reinterpret.(channel.(img))) for channel in (red, green, blue))  # 0-255
    masks = Dict(ink => INK[ink](r, g, b) for ink in unique(s.ink for s in fig.series))

    tables = DataFrame[]
    for p in fig.panels
        markers = [s for s in fig.series if s isa MarkerSeries && s.panel == p.name]
        for s in markers
            rivals = [rival_score(masks[s.ink], p, o) for o in markers if o !== s && o.ink == s.ink]
            pts = find_markers(masks[s.ink], p, s; rivals)
            if !s.saved
                println("  $(rpad(p.name, 10)) $(rpad(s.name, 10)) n = $(length(pts)), rival only, not saved")
                continue
            end
            rows = [(panel = p.name, series = s.name, x = pixel_to_value(p.x, x), y = pixel_to_value(p.y, y),
                     merged = false) for (x, y) in pts]
            isempty(rows) || push!(tables, DataFrame(rows))
            println("  $(rpad(p.name, 10)) $(rpad(s.name, 10)) n = $(length(pts))")
        end
        traces = [s for s in fig.series if s isa TraceSeries && s.panel == p.name]
        for ink in unique(s.ink for s in traces)
            group = filter(s -> s.ink == ink, traces)
            for (name, pts) in trace_bands(masks[ink], p, group)
                isempty(pts) && continue
                if !only(s for s in group if s.name == name).saved
                    println("  $(rpad(p.name, 10)) $(rpad(name, 10)) tracked only, not saved")
                    continue
                end
                rows = [(panel = p.name, series = name, x = pixel_to_value(p.x, x), y = pixel_to_value(p.y, y),
                         merged = m) for (x, y, m) in pts]
                push!(tables, DataFrame(rows))
                println("  $(rpad(p.name, 10)) $(rpad(name, 10)) n = $(length(pts)) columns",
                        length(group) > 1 ? ", $(count(last, pts)) merged" : "")
            end
        end
    end
    data = sort!(vcat(tables...), [:panel, :series, :x])

    path = write_paper_data(exp_raw_path(fig.paper), figure_tab(fig), data, metadata(fig))
    println("  data:        ", relpath(path, pkgdir(@__MODULE__)))
    path = calibration_figure(fig, img, data)
    println("  calibration: ", relpath(path, pkgdir(@__MODULE__)))
    return data
end

# META columns that every workbook has, in this order; condition columns go in between
const META_FIRST = ["figure", "panel", "series", "method", "x_quantity", "x_unit", "y_quantity", "y_unit"]
const META_LAST = ["description", "doi"]

# One META row per series of `fig`: what its x and y are, how they were read, and the
# conditions of its panel and its own.
saved(s::Series) = s.saved

function metadata(fig::PaperFigure)
    rows = map(filter(saved, fig.series)) do s
        p = only(filter(p -> p.name == s.panel, fig.panels))
        xvar, yvar = something(p.xvar, fig.xvar, missing), something(p.yvar, fig.yvar, missing)
        (ismissing(xvar) || ismissing(yvar)) && error("$(fig.image) panel $(p.name): set xvar and yvar")
        DataFrame([(; figure = figure_tab(fig), panel = p.name, series = s.name,
                    method = s isa MarkerSeries ? "markers" : "trace",
                    x_quantity = first(xvar), x_unit = last(xvar), y_quantity = first(yvar), y_unit = last(yvar),
                    p.conditions..., s.conditions..., description = fig.description, doi = fig.doi)])
    end
    return vcat(rows...; cols = :union)
end

# figure tabs in figure order: FIG2 before FIG10
figure_order(tab) = (something(tryparse(Int, filter(isdigit, tab)), typemax(Int)), tab)

"""
    write_paper_data(path, tab, data, meta) → path
Write `data` (columns `panel, series, x, y, merged`) to sheet `tab` of a paper's workbook at
`path`, and `meta` (one row per panel and series, `figure` = `tab`, with the `META_FIRST`
and `META_LAST` columns and any conditions) in place of that tab's rows of its `META`
sheet, keeping all other sheets and rows. `META` comes first, then figure tabs in figure
order, then any other tabs. Used both for digitized figures and for data shared as tables.
"""
function write_paper_data(path, tab, data, meta)
    # the whole workbook is rewritten, as XLSX.jl's read-write mode loses track of sheets
    # when one is replaced among several
    sheets = Dict{String,DataFrame}()
    if isfile(path)
        for name in XLSX.sheetnames(XLSX.readxlsx(path))
            sheets[name] = DataFrame(XLSX.readtable(path, name))
        end
    end
    sheets[tab] = data
    old = get(sheets, "META", nothing)
    # new rows first, so their column order wins over the old file's
    meta = isnothing(old) ? meta : vcat(meta, old[old.figure .!= tab, :]; cols = :union)
    # conditions of figures no longer in the workbook leave empty columns behind
    meta = select(meta, [c for c in names(meta) if !all(ismissing, meta[!, c])])
    meta = select(meta, META_FIRST, Not([META_FIRST; META_LAST]), META_LAST)
    sort!(meta, :figure; by = figure_order, alg = Base.Sort.DEFAULT_STABLE)
    tabs = sort(filter(!=("META"), collect(keys(sheets))); by = figure_order)
    mkpath(dirname(path))
    XLSX.writetable(path, "META" => meta, (t => sheets[t] for t in tabs)...; overwrite = true)
    return path
end

"""
    load_paper_data(paper, figure = nothing) → DataFrame
Data of `paper` (e.g. `"Quintana_etal_2005"`) in `data/exp_raw/<paper>.xlsx`, from the tab
`figure` (e.g. `"FIG3"`, or a non-figure tab such as `"TRACES"`) or from all of its tabs,
joined with their `META` rows: one row per point with `paper`, `figure`, `panel`, `series`,
`x`, `y`, `merged`, how the series was obtained (`method`), the x and y quantities and
units, the conditions (`missing` where one doesn't apply), the description and DOI.
"""
function load_paper_data(paper, figure = nothing)
    path = exp_raw_path(paper)
    meta = DataFrame(XLSX.readtable(path, "META"))
    figures = isnothing(figure) ? unique(meta.figure) : [figure]
    data = vcat([insertcols!(DataFrame(XLSX.readtable(path, f)), 1, :figure => f) for f in figures]...)
    # whole numbers come back from Excel as integers
    data.x, data.y = Float64.(data.x), Float64.(data.y)
    out = leftjoin(data, meta; on = [:figure, :panel, :series], order = :left)
    out = select(out, [c for c in names(out) if !all(ismissing, out[!, c])])
    return insertcols!(out, 1, :paper => paper)
end

# Width of each plot frame in the calibration figure; the crop of the original around it,
# the fonts and the markers are drawn at this frame's scale.
const FRAME_WIDTH = 600

# Rect and Circle fill the marker size exactly; the :rect and :circle symbols are drawn
# smaller than it
plot_marker(shape) = shape === :square ? Rect : shape === :triangle ? :utriangle :
                     shape === :triangle_down ? :dtriangle : shape === :diamond ? :diamond : Circle

# pixel position of data value `v` along the axis mapped by `m`, the inverse of pixel_to_value
value_to_pixel(m::AxisMap, v) = m.offset + m.slope * (m.log ? log10(v) : v)

# tick values and labels from an `xticks`/`yticks` axis attribute (values, or (values,
# labels)); without one, Makie's automatic ticks, or the decades on a log axis
function tick_spec(spec, lo, hi, log::Bool)
    if isnothing(spec)
        values = log ? exp10.(ceil(Int, log10(lo)):floor(Int, log10(hi))) :
                 Makie.get_tickvalues(Makie.WilkinsonTicks(5), lo, hi)
        return values, string.(values)
    end
    spec isa Tuple && return collect(first(spec)), collect(last(spec))
    return collect(spec), string.(collect(spec))
end

"""
    calibration_figure(fig::PaperFigure[, img, data]; colour = :red, alpha = 0.5) → path
Draw `fig`'s digitized `data` over the crop of the original image around each panel: the
points in `colour` with opacity `alpha`, mapped back to pixels through the panel's axis
maps, and the frame and the tick marks and labels of those maps (labels inside the frame),
so they can be checked against the original's. The figure goes to
`papers/<paper>/<image stem>_digitized.png`. Without `img` and `data`, the image and the
saved data are loaded, so the figure can be redrawn without digitizing again.
"""
function calibration_figure(fig::PaperFigure, img, df; colour = :red, alpha = 0.5)
    f = Figure()
    ink = (colour, alpha)
    X(px) = px + 0.5      # 0-based pixel to its centre, x to the right
    Y(py) = -(py + 0.5)   # y down in the image, up in the axis
    for p in fig.panels
        row, col = p.grid
        l, t, r, b = p.crop
        fl, ft, fr, fb = p.frame
        s = something(p.frame_width, FRAME_WIDTH) / (fr - fl)
        fs = something(p.fontsize, fig.fontsize) * s
        ax = Axis(f[row, col]; width = (r - l + 1) * s, height = (b - t + 1) * s, limits = (l, r + 1, -(b + 1), -t))
        hidedecorations!(ax)
        hidespines!(ax)
        image!(ax, (l, r + 1), (-(b + 1), -t), rotr90(img[(t + 1):(b + 1), (l + 1):(r + 1)]))

        lines!(ax, X.([fl, fr, fr, fl, fl]), Y.([fb, fb, ft, ft, fb]); color = ink, linewidth = 1.5)
        attrs = merge(fig.axis, p.axis)
        tl = 1.5 * fig.ticksize
        xlo, xhi = extrema(pixel_to_value.(Ref(p.x), (fl, fr)))
        for (v, label) in zip(tick_spec(get(attrs, :xticks, nothing), xlo, xhi, p.x.log)...)
            xlo <= v <= xhi || continue
            px = X(value_to_pixel(p.x, v))
            lines!(ax, [px, px], Y.([fb - tl, fb + tl]); color = ink, linewidth = 2)
            text!(ax, px, Y(fb - tl - 2); text = label, align = (:center, :bottom), color = ink, fontsize = fs)
        end
        ylo, yhi = extrema(pixel_to_value.(Ref(p.y), (fb, ft)))
        for (v, label) in zip(tick_spec(get(attrs, :yticks, nothing), ylo, yhi, p.y.log)...)
            ylo <= v <= yhi || continue
            py = Y(value_to_pixel(p.y, v))
            lines!(ax, X.([fl - tl, fl + tl]), [py, py]; color = ink, linewidth = 2)
            text!(ax, X(fl + tl + 3), py; text = label, align = (:left, :center), color = ink, fontsize = fs)
        end

        for se in filter(se -> se.panel == p.name && saved(se), fig.series)
            d = df[(df.panel .== p.name) .& (df.series .== se.name), :]
            nrow(d) == 0 && continue
            style = se.hollow ? (color = :transparent, strokecolor = ink, strokewidth = 1.5) : (color = ink,)
            # trace points at least 3 px across, to stay visible over dense data
            size = se isa TraceSeries ? max((2se.radius + 1) * s, 3) : (2se.radius + 1) * s
            scatter!(ax, X.(value_to_pixel.(Ref(p.x), d.x)), Y.(value_to_pixel.(Ref(p.y), d.y));
                     marker = plot_marker(se.shape), markersize = size, style...)
            if se isa TraceSeries && !isempty(se.label)
                text!(ax, X(value_to_pixel(p.x, d.x[end])), Y(value_to_pixel(p.y, d.y[end])); text = se.label, color = ink,
                      fontsize = fs, align = (:right, se.label_below ? :top : :bottom), offset = (0, se.label_below ? -8 : 8))
            end
        end
    end
    resize_to_layout!(f)
    path = joinpath(paper_dir(fig.paper), figure_stem(fig) * "_digitized.png")
    save(path, f)
    return path
end

calibration_figure(fig::PaperFigure; kw...) =
    calibration_figure(fig, RGB{N0f8}.(load(joinpath(paper_dir(fig.paper), fig.image))),
                       load_paper_data(fig.paper, figure_tab(fig)); kw...)
