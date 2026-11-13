# Digitize scatter data from published figures, one section per figure. Each writes its
# data to data/exp_raw/<paper>.xlsx (tab named after the image, columns panel, series, x, y,
# merged), what they are to that workbook's META tab (quantities, units, conditions, how
# they were read), and a calibration figure, the original above the digitized plot, next
# to the image in papers/<paper>/. load_paper_data(paper, figure) reads them back joined.
#
# Pixel positions are 0-based (x right, y down), as image viewers show them.

using Pkg
Pkg.activate(dirname(@__DIR__))  # repo root; this file is in scripts/
using Revise
using SynapseExperimentModeling
using CairoMakie: rich, superscript, subscript

# Figures to digitize when this file is run, by the names they are defined under below,
# e.g. [:quintana_fig1, :quintana_fig3]; :all runs every figure.
RUN = Symbol[]

# tick pixel positions paired with their axis values, for AxisMap
ticks(px, values) = [p => v for (p, v) in zip(px, values)]

## Williams et al. 2023, Cell Calcium, doi:10.1016/j.ceca.2023.102752
# Fig 10a,c: Ca2+ binding of the CaM C- and N-lobes, WT vs N53I CaM, without (a) and with
# (c) CaN peptide. Squares in a, circles in c; the black lines are fits, not data. Methods:
# 3 µM full-length CaM, intrinsic Tyr (C-lobe) and Phe (N-lobe) fluorescence, 50 mM HEPES
# pH 7.40, 100 mM KCl, 5 mM NTA, 1 mM MgCl2, 50 µM EGTA, 22 °C; free Ca2+ from Fluo-5N
# (Kd 85 µM, a) or Fluo-5F (Kd 2.3 µM, c).
# Panels: frame (left, top, right, bottom), x tick pixels at 1E-9 ... 1E-2, y tick pixels at
# 1.0 ... 0.0, legend text box, and the crop shown in the calibration figure.
williams_panel(name, grid, frame, xticks, yticks, legend, crop) = PlotPanel(;
    name, grid, frame, crop, title = name[3:end], skip = [legend],
    conditions = (CaM_lobe = string(last(name)), CaN_peptide = startswith(name, "c"), CaM_uM = 3.0, temperature_C = 22.0,
                  indicator = startswith(name, "c") ? "Fluo-5F" : "Fluo-5N", indicator_Kd_uM = startswith(name, "c") ? 2.3 : 85.0),
    x = AxisMap([p => exp10(e) for (p, e) in zip(xticks, -9:-2)]; log = true),
    y = AxisMap([p => v for (p, v) in zip(yticks, 1.0:-0.2:0.0)]))

williams_fig10 = PaperFigure(;
    paper = "Williams_etal_2023", image = "Fig10.jpg", doi = "10.1016/j.ceca.2023.102752",
    description = "Ca2+ binding of the CaM C- and N-lobes (fluorescence), WT vs N53I CaM, without (a) and with (c) CaN peptide",
    xvar = "Ca_free" => "M", yvar = "normalized_fluorescence" => "1", fontsize = 47, shared_ylabel = true,
    panels = [
        williams_panel("a_CaM-C", (1, 1), (198, 16, 1093, 912),
                       [298.5, 398, 498, 597.5, 695, 794.5, 895, 995],
                       [91, 241, 390.5, 538, 688, 837], (755, 672, 1050, 810), (0, 0, 1130, 1070)),
        williams_panel("a_CaM-N", (1, 2), (1345, 11, 2240, 907),
                       [1445, 1544.5, 1643.5, 1743, 1842, 1941, 2040, 2139.5],
                       [87, 236, 385, 533, 683, 831.5], (1885, 648, 2180, 790), (1240, 0, 2296, 1070)),
        williams_panel("c_CaM-C", (2, 1), (242, 1914, 1136, 2809),
                       [341, 441, 540.5, 639.5, 738.5, 837.5, 936.5, 1037],
                       [1989, 2139, 2287, 2437, 2586.5, 2734], (700, 2385, 1090, 2765), (0, 1885, 1160, 2995)),
        williams_panel("c_CaM-N", (2, 2), (1380, 1912, 2273, 2804),
                       [1480, 1579.5, 1679, 1778, 1876, 1976, 2075, 2174.5],
                       [1987, 2135.5, 2284, 2432, 2581.5, 2730.5], (1842, 2352, 2232, 2732), (1240, 1885, 2296, 2995)),
    ],
    series = [MarkerSeries(; name, panel, ink, shape, radius = 7, conditions = (CaM = name,), colour,
                           label = label * suffix)
              for (panel, shape, suffix) in (("a_CaM-C", :square, ""), ("a_CaM-N", :square, ""),
                                             ("c_CaM-C", :circle, " +\nCaN peptide"),
                                             ("c_CaM-N", :circle, " +\nCaN peptide"))
              for (name, ink, colour, label) in (("WT", :blue, "#2B459E", "WT CaM"),
                                                 ("N53I", :orange, "#F68B25", "N53I CaM"))],
    axis = (; xlabel = rich("[Ca", superscript("2+"), "]", subscript("free")),
            ylabel = "Normalized Fluorescence",
            xticks = (exp10.(-10:2:-2), ["1E-10", "1E-8", "1E-6", "1E-4", "1E-2"]),
            xminorticks = exp10.(-10:-1), xminorticksvisible = true,
            yticks = (0:0.2:1, ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])),
)

# Fig 7a,b (left): Ca2+ dissociation from isolated N53I CaM at 20 °C, by Tyr fluorescence after
# rapid mixing with EGTA (a) and by quin-2 fluorescence (b). Dense orange data points with a
# black single-exponential fit on top; only the orange ink is read, as one band per pixel column
# (the centre of the noise cloud), bridging the fit line. The rate summaries on the right are
# fitted values and panels c, d contain CaN peptide; none are read.
williams_fig7 = PaperFigure(;
    paper = "Williams_etal_2023", image = "Fig7.jpg", doi = "10.1016/j.ceca.2023.102752",
    description = "Stopped-flow Ca2+ dissociation from isolated N53I CaM at 20 C (Applied Photophysics SX.18MV, dead " *
                  "time 1.4 ms): 3 uM CaM with 30 uM CaCl2 in 10 mM MOPS pH 7.0, 150 mM KCl (Methods, given for the " *
                  "CaN peptide experiments), rapidly mixed with 10 mM EGTA (a; Tyr fluorescence, ex 275 nm, UG1 " *
                  "filter) or 150 uM quin-2 (b; ex 330 nm, em 510 nm band-pass). Averages of at least three traces. " *
                  "The black lines are exponential fits, not data. Points are the centre of the band the dense data " *
                  "points form, one per pixel column, not the individual measurements",
    xvar = "time" => "s", fontsize = 40, linewidth = 3, ticksize = 12,
    panels = [PlotPanel(; name = "a", title = "", grid = (1, 1), frame = (204, 100, 1240, 759), crop = (0, 0, 1300, 930),
                        conditions = (CaM = "N53I", CaM_uM = 3.0, Ca_total_uM = 30.0, EGTA_uM = 10000.0, temperature_C = 20.0),
                        yvar = "Tyr_fluorescence" => "mV",
                        x = AxisMap(ticks([203.5, 546.5, 890.5, 1230.5], 0:3)),
                        y = AxisMap(ticks([119.0, 278.5, 438.5, 600.0, 759.0], 9.0:-0.5:7.0)),
                        axis = (; yticks = 7.0:0.5:9.0)),
              PlotPanel(; name = "b", title = "", grid = (2, 1), frame = (213, 1040, 1250, 1626), crop = (0, 930, 1300, 1803),
                        conditions = (CaM = "N53I", CaM_uM = 3.0, Ca_total_uM = 30.0, quin2_uM = 150.0, temperature_C = 20.0),
                        yvar = "quin2_fluorescence" => "mV",
                        x = AxisMap(ticks([212.0, 555.0, 897.0, 1240.5], 0:3)),
                        y = AxisMap(ticks([1055.5, 1149.5, 1244.5, 1340.0, 1435.0, 1530.5, 1626.0], 8.2:-0.2:7.0)),
                        axis = (; yticks = 7.0:0.2:8.2))],
    series = [TraceSeries(; name = "N53I CaM", panel, ink = :orange, gap = 30, radius = 1, colour = "#F68B25")
              for panel in ("a", "b")],
    axis = (; xlabel = "Time (sec)", ylabel = "Fluorescence (mV)", xticks = (0:3, ["0.00", "1.00", "2.00", "3.00"]),
            topspinevisible = false, rightspinevisible = false, xtickalign = 0, ytickalign = 0))

## Quintana et al. 2005, BBRC, doi:10.1016/j.bbrc.2005.06.152
# Small JPEGs of Origin plots: 1 px axes with outward ticks, black open markers. Dense
# kinetic traces are read as bands (TraceSeries), sparse points as markers.
quintana_axis =(; xtickalign = 0, ytickalign = 0, xminortickalign = 0, yminortickalign = 0,
                 xtickcolor = :black, ytickcolor = :black, xminortickcolor = :black, yminortickcolor = :black,
                 xticklabelfont = :regular, yticklabelfont = :regular, xlabelfont = :regular, ylabelfont = :regular,
                 topspinevisible = false, rightspinevisible = false)
quintana(image; kw...) = PaperFigure(; paper = "Quintana_etal_2005", image, doi = "10.1016/j.bbrc.2005.06.152",
                                     fontsize = 10, linewidth = 1, ticksize = 4, kw...)

# Fig 1B: CaN activity vs CaM (means; the SD bars are not read), with the Hill fit;
# triangles are CaNA only.
# The ~8 px open markers barely stand out from the fit line, so their approximate centres
# (read off the image) are given and snapped to the best template match.
quintana_fig1 = quintana("Fig1.jpg";
    description = "CaN activity (RII-P dephosphorylation) vs CaM; means, SD bars not read",
    xvar = "CaM" => "nM", yvar = "CaN_activity" => "umol/min/mg",
    panels = [PlotPanel(; name = "B", grid = (1, 1), frame = (198, 12, 443, 211), margin = 2, crop = (138, 0, 443, 246),
                        x = AxisMap(ticks([198, 248.5, 299, 350, 400], exp10.(-1:3)); log = true),
                        y = AxisMap(ticks([26, 54.5, 83, 111, 140, 168, 197], 0.30:-0.05:0.0)))],
    series = [MarkerSeries(; name, panel = "B", ink = :black, shape, radius = 4, hollow = true, near,
                           conditions = (enzyme = name,), colour = "#000000", label = "")
              for (name, shape, near) in (
                  ("CaN", :circle, [(222, 191), (259, 191), (274, 178), (294, 151), (309, 127), (314, 112),
                                    (325, 106), (360, 61), (375, 45), (411, 48), (426, 27)]),
                  ("CaNA", :triangle, [(224, 191), (426, 192)]))],
    axis = (; quintana_axis..., xlabel = "CaM (nM)", ylabel = "Calcineurin Activity (µmol/min/mg)",
            xticks = (exp10.(-1:3), ["0.1", "1", "10", "100", "1000"]),
            xminorticks = [k * exp10(e) for e in -1:3 for k in 2:9], xminorticksvisible = true,
            yticks = (0:0.05:0.3, ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"])))

# Fig 2B,C: steady-state dissociation of CaM from CaN (B) and CaNA only (C), open circles
# merging into a band along the single-exponential fit.
quintana_fig2 = quintana("Fig2.jpg";
    description = "Steady-state dissociation of Ca2+/CaM(C75)ACR from CaN (B) or CaNA only (C) on adding excess unlabelled CaM",
    xvar = "time" => "s", yvar = "fluorescence" => "a.u.",
    panels = [PlotPanel(; name = "B", grid = (1, 1), frame = (439, 0, 728, 235), margin = 2, crop = (360, 0, 729, 290),
                        axis = (; xticks = 0:1000:6000, yticks = (60000:20000:140000, string.(60000:20000:140000))),
                        skip = [(565, 95, 675, 118)],
                        x = AxisMap(ticks([439, 483.5, 528, 573, 617.5, 662, 706.5], 0:1000:6000)),
                        y = AxisMap(ticks([23.5, 70.5, 118, 164.5, 211.5], 140000:-20000:60000))),
              PlotPanel(; name = "C", grid = (1, 2), frame = (253, 300, 544, 536), margin = 2, crop = (190, 285, 600, 577),
                        axis = (; xticks = 0:2000:8000, yticks = (60000:10000:100000, string.(60000:10000:100000))),
                        skip = [(375, 370, 485, 395)],
                        x = AxisMap(ticks([282, 340.5, 398.5, 456, 514], 0:2000:8000)),
                        y = AxisMap(ticks([325, 375, 424, 474, 523], 100000:-10000:60000)))],
    series = [TraceSeries(; name, panel, ink = :black, hollow = true, gap = 6, conditions = (enzyme = name,),
                          radius = 3, colour = "#000000")
              for (name, panel) in (("CaN", "B"), ("CaNA", "C"))],
    axis = (; quintana_axis..., xlabel = "Time (sec)", ylabel = "Fluorescence (Arbitrary Units)"))

# Fig 3: stopped-flow association at 4.0, 1.0 and 0.5 µM CaN (open circles merging into
# three bands; the axis says ms but the data are recorded in s), and the inset of the fitted rate against CaN
# concentration (labelled "CaM (M^-1)" in the figure).
quintana_fig3 = quintana("Fig3.jpg";
    description = "Stopped-flow association of CaM(C75)ACR with CaN: traces at three CaN concentrations (main) " *
                  "and fitted rate vs CaN concentration (inset). The main time axis is labelled ms in the paper " *
                  "but recorded here in s: the inset rates (~4-104 1/s) and the 1.7 ms stopped-flow dead time " *
                  "only fit seconds",
    panels = [PlotPanel(; name = "main", title = "", grid = (1, 1), frame = (51, 4, 344, 242), margin = 2,
                        crop = (0, 0, 354, 278), xvar = "time" => "s", yvar = "fluorescence" => "a.u.",
                        axis = (; xlabel = "Time (ms)", ylabel = "Fluorescence (Arbitraty Units)",
                                xticks = (0:0.05:0.25, ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"]),
                                yticks = (-1.0:0.05:-0.6, ["-1.00", "-0.95", "-0.90", "-0.85", "-0.80",
                                                           "-0.75", "-0.70", "-0.65", "-0.60"])),
                        skip = [(176, 112, 335, 240), (285, 3, 335, 19), (290, 56, 335, 71), (288, 99, 335, 114)],
                        x = AxisMap(ticks([51, 109, 168, 227, 285, 344], 0:0.05:0.25)),
                        y = AxisMap(ticks([4, 34, 64, 93, 123, 182, 212, 242], [-0.60, -0.65, -0.70, -0.75, -0.80, -0.90, -0.95, -1.00]))),
              PlotPanel(; name = "inset", title = "", grid = (1, 2), frame = (199, 119, 321, 217), margin = 2,
                        crop = (176, 108, 334, 240), xvar = "CaN" => "M", yvar = "rate" => "1/s", fontsize = 5,
                        axis = (; xlabel = rich("CaM (M", superscript("-1"), ")"), ylabel = rich("Rate (s", superscript("-1"), ")"),
                                xticks = (0:5e-7:2e-6, ["0.0", "5.0x10⁻⁷", "1.0x10⁻⁶", "1.5x10⁻⁶", "2.0x10⁻⁶"]),
                                yticks = 0:20:120),
                        skip = [(200, 126, 285, 143)],
                        x = AxisMap(ticks([211, 235.5, 260, 284, 308.5], 0:5e-7:2e-6)),
                        y = AxisMap(ticks([118.5, 133.75, 149, 164, 179, 194, 210], 120:-20:0)))],
    series = [[TraceSeries(; name, panel = "main", ink = :black, curve, hollow = true, conditions = (CaN_uM = conc,),
                           radius = 3, colour = "#000000", label, label_below = curve == 3)
               for (curve, conc, name, label) in ((1, 4.0, "4.0 uM", "4.0 µM"), (2, 1.0, "1.0 uM", "1.0 µM"),
                                                  (3, 0.5, "0.5 uM", "0.5 µM"))];
              MarkerSeries(; name = "rate", panel = "inset", ink = :black, radius = 2, hollow = true,
                           near = [(214, 207), (218, 204), (224, 199), (236, 189), (260, 170), (308, 131)],
                           colour = "#000000", label = "")],
    axis = quintana_axis)

# Fig 4: stopped-flow dissociation after EGTA, from CaN (A) and CaNA only (B); grey markers
# merged into a band under the black fit.
quintana_fig4 = quintana("Fig4.jpg";
    description = "Stopped-flow dissociation of Ca2+/CaM(C75)ACR from CaN (A) or CaNA only (B) on chelating Ca2+ with EGTA",
    xvar = "time" => "s", yvar = "fluorescence" => "a.u.",
    panels = [PlotPanel(; name = "A", grid = (1, 1), frame = (43, 5, 301, 216), margin = 2, crop = (0, 0, 308, 250),
                        axis = (; yticks = (1.0:0.5:3.5, ["1.0", "1.5", "2.0", "2.5", "3.0", "3.5"])),
                        skip = [(118, 75, 250, 113)],
                        x = AxisMap(ticks([63, 102.5, 142, 182, 222, 262, 301], 0:2:12)),
                        y = AxisMap(ticks([5.5, 44, 82, 120, 158, 197], 3.5:-0.5:1.0))),
              PlotPanel(; name = "B", grid = (2, 1), frame = (42, 274, 301, 485), margin = 2, crop = (0, 255, 308, 518),
                        axis = (; yticks = (0.8:0.2:2.4, ["0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"])),
                        skip = [(145, 385, 230, 408)],
                        x = AxisMap(ticks([62, 102, 141.5, 181, 221, 261, 301], 0:2:12)),
                        y = AxisMap(ticks([287, 311.5, 336, 361, 386, 411, 435, 460, 485], 2.4:-0.2:0.8)))],
    series = [TraceSeries(; name, panel, ink = :grey, gap = 4, conditions = (enzyme = name,), radius = 4,
                          colour = "#808080")
              for (name, panel) in (("CaN", "A"), ("CaNA", "B"))],
    axis = (; quintana_axis..., xlabel = "Time (sec)", ylabel = "Fluorescence (Arbitrary Units)", xticks = 0:2:12))

## Shifman et al. 2006, PNAS, doi:10.1073/pnas.0606433103
# Kaleidagraph plots: 3 px frames with inward ticks on all four sides, bold tick labels,
# filled black markers (open ones where noted); lines are fits. Legend symbols are decoded
# from the PDF's symbol-font glyphs. Error bars are not read.
shifman_axis = (; xtickalign = 1, ytickalign = 1, xtickcolor = :black, ytickcolor = :black,
                xlabelfont = :regular, ylabelfont = :regular)
shifman(image; kw...) = PaperFigure(; paper = "Shifman_etal_2006", image, doi = "10.1073/pnas.0606433103",
                                    linewidth = 3, ticksize = 12, kw...)

# Fig 1B: Ca2+ binding to 5 µM WT CaM in a competition assay with 5 µM Fluo4FF at 25 °C,
# alone (triangles), with 5 µM CaMKII (squares) or 5 µM CaMKII-cbp (diamonds).
# Squares and diamonds overlap each other and pass for triangle clusters, so their centres
# are given; the triangles are found automatically, centred on their template peaks as
# they stack into merged blobs, with a threshold above the fit lines' dashes.
shifman_fig1_squares = [(272, 1310), (280, 1246), (288, 1179), (303, 1128), (305, 1047), (312, 1097), (314, 982),
                        (324, 916), (326, 1026), (341, 868), (344, 974), (356, 899), (358, 822), (373, 839),
                        (375, 772), (386, 766.5), (395.5, 731.5), (410, 731.5), (429, 676), (459, 643.5),
                        (484.5, 609.5), (497, 538), (539, 609), (564, 604), (576, 528), (758, 632), (779, 512),
                        (975, 636), (1107, 464)]
shifman_fig1_diamonds = [(235, 1594), (235, 1297), (235, 1250), (236, 1204), (236, 1157), (236, 1109), (237, 1063),
                         (237, 1017), (238, 968), (239, 922), (240, 878), (241, 836), (243, 794), (248, 756),
                         (252, 725), (262, 704), (286, 725), (316, 690), (358, 672), (400, 659), (443, 642),
                         (491.5, 638), (542, 638), (579.5, 603.5), (628, 596.5)]
shifman_fig1 = shifman("Fig1.jpeg";
    description = "Ca2+ binding to 5 uM WT CaM alone or with 5 uM CaMKII or CaMKII-cbp, measured by competition " *
                  "with 5 uM Fluo4FF at 25 C; dotted lines are equilibrium fits, not data",
    xvar = "Ca_free" => "uM", yvar = "bound_Ca_per_CaM" => "1", fontsize = 42,
    panels = [PlotPanel(; name = "B", grid = (1, 1), frame = (223, 400, 1437, 1628), crop = (0, 370, 1478, 1799),
                        conditions = (CaM = "WT", CaM_uM = 5.0, Fluo4FF_uM = 5.0, temperature_C = 25.0),
                        skip = [(765, 1320, 1400, 1570), (223, 400, 1437, 414), (223, 1614, 1437, 1628)],
                        x = AxisMap(ticks([441.5, 664, 884, 1107, 1326], 10:10:50)),
                        y = AxisMap(ticks([641, 883, 1124, 1364, 1605], 4:-1:0)))],
    series = [MarkerSeries(; name, panel = "B", ink = :black, shape, radius, near, conditions = (; target, target_uM),
                           locate = :peak, threshold = 0.54, colour = "#000000", label = name)
              for (name, shape, radius, near, target, target_uM) in (
                  ("CaM", :triangle, 13, NTuple{2,Float64}[], "none", 0.0),
                  ("CaM + CaMKII", :square, 11, shifman_fig1_squares, "CaMKII", 5.0),
                  ("CaM + CaMKII-cbp", :diamond, 17, shifman_fig1_diamonds, "CaMKII-cbp", 5.0))],
    axis = (; shifman_axis..., xlabel = rich("[Ca", superscript("2+"), subscript("free"), "], µM"),
            ylabel = rich("Bound Ca", superscript("2+"), " per CaM"), xticks = 0:10:50, yticks = 0:5))

# Fig 3D: Fluo4FF fraction bound as Ca2+ is titrated into 10 µM Fluo4FF alone (filled
# squares) or with 10 µM CaM-C^WT (open down triangles), CaM-N^WT (filled triangles) or WT
# CaM (filled circles); the thick lines are fits. The legend gives circles for Fluo4FF alone
# and squares for WT CaM, but the squares are what Fluo4FF alone gives, so the two are
# swapped here (see the description). The x axis is broken, so its two segments are read as
# panels D and D_right. Markers sit on the fit lines and on each other, so their centres are
# given; the filled triangles are mostly hidden under the open ones and only those in view
# are read.
shifman_fig3_near = Dict(
    ("D", "WT CaM") => [(1122.7, 1515), (1134.3, 1505), (1146, 1481.7), (1169.3, 1445), (1194.3, 1409.3),
                         (1218.3, 1371.5), (1241.7, 1317.3), (1266.7, 1288.3), (1291.7, 1250), (1316, 1222.3),
                         (1340, 1197.3), (1365, 1170), (1389, 1147), (1412.3, 1123.3), (1460.7, 1076.7),
                         (1509.3, 1071.7)],
    ("D", "CaM-C^WT") => [(1119.3, 1493.3), (1131, 1476), (1145.3, 1456.7), (1172, 1418.3), (1195.5, 1384),
                          (1219.3, 1346), (1242.7, 1284), (1266.7, 1249.3), (1291.7, 1221.7), (1316, 1198.3),
                          (1340, 1161.7), (1365, 1139.3), (1389, 1115), (1411.7, 1097.3), (1461.7, 1065),
                          (1510, 1047.3)],
    ("D", "CaM-N^WT") => [(1117.7, 1501.7), (1132.7, 1466.7), (1147.7, 1446), (1167.7, 1401.7), (1191.7, 1365),
                          (1215, 1323.3), (1243.3, 1276.7), (1265.7, 1240), (1291.7, 1212.7), (1341.7, 1156.7),
                          (1390, 1123.3), (1433.3, 1093.3)],
    ("D", "Fluo4FF") => [(1130.7, 1381), (1143.3, 1341.7), (1155, 1318.3), (1179.3, 1254.3), (1204, 1205),
                         (1228.3, 1194), (1252.3, 1153.7), (1276.7, 1127.7), (1301, 1119), (1325.7, 1093.3),
                         (1350, 1075.7), (1374, 1073.3), (1399.3, 1061.7), (1422.3, 1054.3), (1471.7, 1028.3),
                         (1520, 1021)],
    ("D_right", "WT CaM") => NTuple{2,Float64}[],
    ("D_right", "CaM-C^WT") => [(1602.5, 965.5), (1700, 964)],
    ("D_right", "CaM-N^WT") => [(1676, 964)],
    ("D_right", "Fluo4FF") => [(1610.5, 955.5), (1705.5, 950)])
shifman_fig3_series = (("WT CaM", :circle, 8, false, "WT", 10.0), ("CaM-C^WT", :triangle_down, 9, true, "CaM-C^WT", 10.0),
                       ("CaM-N^WT", :triangle, 8, false, "CaM-N^WT", 10.0), ("Fluo4FF", :square, 8, false, "none", 0.0))
shifman_fig3_y = AxisMap(ticks([953.5, 1083.5, 1212.5, 1342.5, 1472.5], 1.0:-0.2:0.2))
shifman_fig3 = shifman("Fig3.jpeg";
    description = "Fluo4FF fraction bound vs total Ca2+ titrated into 10 uM Fluo4FF alone or with 10 uM WT or mutant " *
                  "CaM (axis labels: Fraction Bound, [Ca2+]). Lines are fits, not data. The x axis is broken: panel D " *
                  "is 0-90 uM, D_right 220-240 uM. Legend symbols were decoded from the PDF; its filled circle " *
                  "(Fluo4FF alone) and filled square (WT CaM) are swapped here: the squares match 10 uM Fluo4FF " *
                  "alone with its KD of 13 uM against total Ca2+ and exceed the most WT CaM could bind at 10-20 uM, " *
                  "and with the paper's KDs (Table 2, as stepwise constants) every curve matches Fluo4FF's fraction " *
                  "bound with its competing CaM, not CaM's",
    xvar = "Ca_total" => "uM", yvar = "Fluo4FF_fraction_bound" => "1", fontsize = 34, shared_ylabel = true,
    panels = [PlotPanel(; name = "D", grid = (1, 1), frame = (1086, 920, 1535, 1602), crop = (950, 880, 1545, 1713),
                        conditions = (Fluo4FF_uM = 10.0, temperature_C = 25.0),
                        skip = [(1086, 920, 1535, 932), (1086, 1590, 1535, 1602)],
                        x = AxisMap(ticks([1184, 1280.5, 1378.5, 1475.5], 20:20:80)), y = shifman_fig3_y,
                        axis = (; xticks = 0:20:80)),
              PlotPanel(; name = "D_right", title = "", grid = (1, 2), frame = (1588, 920, 1763, 1602),
                        crop = (1580, 880, 1799, 1713), conditions = (Fluo4FF_uM = 10.0, temperature_C = 25.0),
                        frame_width = 600 * (1763 - 1588) / (1535 - 1086),
                        skip = [(1588, 920, 1763, 932), (1588, 1590, 1763, 1602)],
                        x = AxisMap(ticks([1665.5, 1762.5], 220:20:240)), y = shifman_fig3_y,
                        axis = (; xticks = 220:20:240))],
    series = [MarkerSeries(; name, panel, ink = hollow ? :grey : :black, shape, radius, hollow,
                           near = shifman_fig3_near[(panel, name)], conditions = (; CaM, CaM_uM),
                           colour = "#000000", label = "")
              for panel in ("D", "D_right") for (name, shape, radius, hollow, CaM, CaM_uM) in shifman_fig3_series
              if !isempty(shifman_fig3_near[(panel, name)])],
    axis = (; shifman_axis..., xlabel = rich("[Ca", superscript("2+"), "], µM"), ylabel = "Fraction Bound",
            yticks = 0:0.2:1))

# Fig 4B-D: CaMKII (0.28 µM) autophosphorylation rate, normalized, with WT CaM (filled down
# triangles), CaM-C^WT (circles) or CaM-N^WT (squares): against CaM at 300 µM Ca2+ in 30 s
# reactions (B), against Ca2+ with 5 µM CaM (C) and against CaM at 12 µM Ca2+ (D), both in
# 1 min reactions. Means of four experiments; the error bars (SD) are not read. Markers sit
# on error bars, fits and labels, so their centres are given; the points at 0 lie on the
# frame and are not read.
shifman_fig4_series = (("WT", :triangle_down), ("CaM-C^WT", :circle), ("CaM-N^WT", :square))
shifman_fig4_near = Dict(
    ("B", "CaM-C^WT") => [(1119, 576.7), (1153.5, 385), (1208.5, 247.5), (1265, 132.5), (1321, 133.5), (1376, 94),
                          (1433.5, 97.5)],
    ("B", "CaM-N^WT") => [(1107, 609), (1119.7, 607.5), (1151.8, 599.5), (1208.5, 533.5), (1321, 391.5), (1432, 350),
                          (1543.5, 346), (1655, 309)],
    ("C", "WT") => [(309.5, 1212.5), (527.5, 1016), (801, 924)],
    ("C", "CaM-C^WT") => [(312, 1342.5), (527.5, 1262.5), (802.5, 1101)],
    ("C", "CaM-N^WT") => [(309.5, 1369), (527.5, 1345), (802.5, 1304)],
    ("D", "WT") => [(1211, 1135), (1332, 959), (1468.5, 842.5), (1743.5, 992.5)],
    ("D", "CaM-C^WT") => [(1213, 1390), (1332, 1270), (1468.5, 1165), (1743.5, 1161.5)],
    ("D", "CaM-N^WT") => [(1205, 1389), (1332, 1377.5), (1468.5, 1332.5), (1743.5, 1294)])
shifman_fig4 = shifman("Fig4.jpeg";
    description = "Rate of CaMKII (0.28 uM) autophosphorylation, normalized, with WT, CaM-C^WT or CaM-N^WT CaM: vs " *
                  "CaM at 300 uM Ca2+, 30 s reactions, normalized to saturating Ca2+/CaM for 1 min (B); vs Ca2+ with " *
                  "5 uM CaM (C) and vs CaM at 12 uM Ca2+ (D), 1 min reactions. Means of four experiments, SD bars " *
                  "not read; lines are fits; the points at 0 lie on the frame and are not read",
    xvar = "CaM" => "uM", yvar = "autophosphorylation_rate" => "1", fontsize = 28, linewidth = 5,
    panels = [PlotPanel(; name = "B", grid = (1, 1), frame = (1096, 5, 1711, 621), crop = (900, 0, 1799, 740),
                        conditions = (CaMKII_uM = 0.28, Ca_uM = 300.0, reaction_s = 30.0),
                        x = AxisMap(ticks([1207.5, 1319.5, 1431.5, 1543.5, 1655], 10:10:50)),
                        y = AxisMap(ticks([72.5, 208.5, 343.5, 479, 614], 0.8:-0.2:0.0)),
                        axis = (; xlabel = "[CaM], µM", xticks = 0:10:50, yticks = 0:0.2:0.8)),
              PlotPanel(; name = "C", grid = (2, 1), frame = (256, 833, 855, 1433), crop = (0, 780, 870, 1536),
                        conditions = (CaMKII_uM = 0.28, CaM_uM = 5.0, reaction_s = 60.0), xvar = "Ca" => "uM",
                        x = AxisMap(ticks([364.5, 473.5, 582, 691, 799.5], 2:2:10)),
                        y = AxisMap(ticks([926, 1018, 1110, 1202.5, 1294.5, 1386.5], 1.0:-0.2:0.0)),
                        axis = (; xlabel = rich("[Ca", superscript("2+"), "], µM"), xticks = 0:2:10, yticks = 0:0.2:1.2)),
              PlotPanel(; name = "D", grid = (2, 2), frame = (1196, 832, 1797, 1433), crop = (900, 780, 1799, 1536),
                        conditions = (CaMKII_uM = 0.28, Ca_uM = 12.0, reaction_s = 60.0),
                        x = AxisMap(ticks([1332, 1469, 1605.5, 1741.5], 5:5:20)),
                        y = AxisMap(ticks([913, 993, 1072.5, 1153.5, 1233, 1313.5, 1393.5], 1.2:-0.2:0.0)),
                        axis = (; xlabel = "[CaM], µM", xticks = 0:5:20, yticks = 0:0.2:1.4))],
    series = [MarkerSeries(; name, panel, ink = :black, shape, radius = 8, near = shifman_fig4_near[(panel, name)],
                           conditions = (CaM = name,), colour = "#000000", label = "")
              for panel in ("B", "C", "D") for (name, shape) in shifman_fig4_series if haskey(shifman_fig4_near, (panel, name))],
    axis = (; shifman_axis..., ylabel = "Rate of Autophosphorylation\n(normalized)"))

# Fig 5: time course of CaMKII autophosphorylation in saturating Ca2+/CaM (quench flow),
# normalized, with WT CaM (circles), CaM-N^WT (filled triangles) or CaM-C^WT (open down
# triangles); error bars not read. Centres are given (markers on error bars and fits); the
# 0.1 s markers lie on the frame line, so they are read by eye.
shifman_fig5 = shifman("Fig5.jpeg";
    description = "Time course of CaMKII autophosphorylation (quench flow, 30 C) with WT, CaM-C^WT or CaM-N^WT CaM, " *
                  "normalized to fully autophosphorylated CaMKII (30 s WT, 100 s mutants). Concentrations after 1:1 " *
                  "mixing of the two solutions (0.7 uM CaMKII subunits, 2 mM Ca2+, 100 uM ATP, and 6 uM WT, 15 uM " *
                  "CaM-C^WT or 30 uM CaM-N^WT). Means with error bars not read; lines are fits; the 0.1 s markers lie " *
                  "on the frame line and were read by eye",
    xvar = "time" => "s", yvar = "autophosphorylated_CaMKII" => "1", fontsize = 45, linewidth = 11, ticksize = 25,
    panels = [PlotPanel(; name = "main", title = "", grid = (1, 1), frame = (350, 13, 1786, 1450), margin = 7,
                        crop = (0, 0, 1799, 1627),
                        conditions = (CaMKII_uM = 0.7, Ca_uM = 2000.0, ATP_uM = 100.0, temperature_C = 30.0),
                        x = AxisMap(ticks([581.5, 812.5, 1044.5, 1276, 1508.5, 1739.5], 5:5:30)),
                        y = AxisMap(ticks([144, 405, 666, 927.5, 1188.5], 1.0:-0.2:0.2)))],
    series = [MarkerSeries(; name = "WT", panel = "main", ink = :black, shape = :circle, radius = 17,
                           near = [(368.5, 1024), (397.25, 621.25), (485, 245), (813, 198)], fixed = [(356, 1244)],
                           conditions = (CaM = "WT", CaM_uM = 6.0), colour = "#000000", label = ""),
              MarkerSeries(; name = "CaM-N^WT", panel = "main", ink = :black, shape = :triangle, radius = 17,
                           near = [(366, 1270), (396, 1127.5), (485, 948), (813, 580), (1738, 412)], fixed = [(353.5, 1349)],
                           conditions = (CaM = "CaM-N^WT", CaM_uM = 30.0), colour = "#000000", label = ""),
              MarkerSeries(; name = "CaM-C^WT", panel = "main", ink = :grey, shape = :triangle_down, radius = 18, hollow = true,
                           near = [(363.5, 1337.5), (397.25, 1300), (485, 1165), (813, 764), (1738, 388)],
                           fixed = [(356, 1391)], conditions = (CaM = "CaM-C^WT", CaM_uM = 15.0), colour = "#000000",
                           label = "")],
    axis = (; shifman_axis..., xlabel = "Time (sec)", ylabel = "Autophosphorylated CaMKII\n(normalized)",
            xticks = 0:5:30, yticks = (0:0.2:1, ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])))

## VanScyoc et al. 2002, Biophys J, doi:10.1016/S0006-3495(02)75286-7
# Fig 3: Ca2+ titrations of 6 µM CaM(1-75) (A), CaM(76-148) (B) and CaM(1-148) (C) at 22 °C,
# monitored by Phe fluorescence (ex/em 250/280 nm, filled triangles; sites I and II) and Tyr
# fluorescence (277/320 nm, open triangles; sites III and IV). Image: the journal web page's
# GIF (3335 x 4239, black and white), converted to PNG. 8 px frames on all four sides with
# inward ticks; lines are fits. In A and B the Ca2+-independent series sit on a dotted line
# with their origin offset to -0.1. Open triangles are drawn with ~8 px strokes, and both
# kinds sit on fit lines, so centres are the template peaks. Each titration point is read at
# both wavelength pairs, so the two series share their x values (within 1-3 px in A and C).
# In B's dense mid-titration cluster the filled triangles of the offset row overlap into one
# band and the open ones crowd each other; the markers no template peak separates are given
# as fixed centres: the filled ones 29 px below their apexes (the offset of the separate
# markers' peaks), the open ones at their template peaks below the match threshold.
vanscyoc_fixed = Dict(
    ("B", "Phe") => [(1393, 2476), (1427, 2476), (1463, 2476), (1487, 2476), (1517, 2476), (1537, 2482),
                     (1557, 2476), (1676, 2476), (1702, 2476), (1730, 2476), (1752, 2476)],
    ("B", "Tyr") => [(1702, 1648), (1730, 1621), (1751, 1608)])
vanscyoc_x = AxisMap(ticks([444.5, 1145.5, 1846.5, 2546.5], exp10.(-7:-4)); log = true)
vanscyoc_panel(name, grid, frame, yticks, letter, crop, CaM; axis = (;)) = PlotPanel(;
    name, grid, frame, crop, margin = 8, frame_width = 900, skip = [letter],
    conditions = (CaM = CaM, CaM_uM = 6.0, temperature_C = 22.0, OGB5N_uM = 0.1, OGB5N_Kd_uM = 29.6),
    x = vanscyoc_x, y = AxisMap(ticks(yticks, 1.0:-0.5:0.0)), axis)
vanscyoc_fig3 = PaperFigure(;
    paper = "VanScyoc_etal_2002", image = "Fig3.png", doi = "10.1016/S0006-3495(02)75286-7",
    description = "Ca2+ titrations of 6 uM CaM(1-75) (A), CaM(76-148) (B) and CaM(1-148) (C) at 22 C, monitored by " *
                  "Phe fluorescence (ex/em 250/280 nm; sites I and II) and Tyr fluorescence (277/320 nm; sites III and " *
                  "IV), normalized as (f - fmin)/(fmax - fmin); free Ca2+ from 0.1 uM Oregon Green 488 BAPTA-5N " *
                  "(OGB-5N), as Kd (f_high - f)/(f - f_low) with their Kd of 29.6 uM (50 mM HEPES pH 7.4, 100 mM " *
                  "KCl, 22 C). Lines are fits, " *
                  "not data. In A and B the Ca2+-independent series (Tyr in A, Phe in B) are plotted, and recorded here, " *
                  "with their origin offset to -0.1",
    xvar = "Ca_free" => "M", yvar = "normalized_fluorescence" => "1", fontsize = 95, linewidth = 8, ticksize = 40,
    panels = [vanscyoc_panel("A", (1, 1), (378, 74, 3249, 1265), [161.5, 623.5, 1077.5], (412, 222, 513, 329),
                             (0, 30, 3334, 1290), "CaM(1-75)"; axis = (; xticklabelsvisible = false)),
              vanscyoc_panel("B", (2, 1), (378, 1380, 3249, 2571), [1471.5, 1931.5, 2386.5], (418, 1530, 509, 1639),
                             (0, 1340, 3334, 2600), "CaM(76-148)";
                             axis = (; xticklabelsvisible = false, ylabel = "Normalized Fluorescence Intensity")),
              vanscyoc_panel("C", (3, 1), (378, 2711, 3249, 3899), [2806.5, 3307.5, 3803.5], (402, 2877, 497, 2988),
                             (0, 2670, 3334, 4238), "CaM(1-148)";
                             axis = (; xlabel = rich("[Calcium] ", subscript("free"))))],
    series = [MarkerSeries(; name, panel, ink = :black, shape = :triangle, radius = 29, hollow, stroke = 8, locate = :peak,
                           threshold = 0.75, fixed = get(vanscyoc_fixed, (panel, name), NTuple{2,Float64}[]),
                           conditions = (ex_nm = ex, em_nm = em), colour = "#000000", label = "")
              for panel in ("A", "B", "C")
              for (name, hollow, ex, em) in (("Phe", false, 250.0, 280.0), ("Tyr", true, 277.0, 320.0))],
    axis = (; ylabel = "", xticks = (exp10.(-7:-3), [rich("10", superscript(string(e))) for e in -7:-3]),
            xminorticks = [k * exp10(e) for e in -8:-4 for k in 2:9], xminorticksvisible = true,
            yticks = (0:0.5:1, ["0.0", "0.5", "1.0"])))

# Fig 4: Ca2+ titrations of the full-length mutants R74A (A), R90A (B) and R90G (C), read as in
# Fig 3 (same image source and markers): filled triangles at 250/280 nm (Phe), open at 277/320
# nm (Tyr); solid lines are fits, dashed lines the CaM(1-148) fit for comparison. In B the
# filled triangles rise and the open ones fall, the reverse of A, C and Fig 3C. The x axis runs
# from 10^-7 (left frame) to 10^-3 M (right frame). Panel letters are skipped. Markers no
# template peak separates are fixed centres, 29 px below their apexes as in Fig 3: the pair
# where the curves cross in A, B's open triangle at 4.2 µM (wholly under the filled one, so
# given its centre), and C's two points at 10^-7 M on the frame line.
vanscyoc_x4 = AxisMap(ticks([357.5, 1077.5, 1806.5, 2524.5, 3253.5], exp10.(-7:-3)); log = true)
vanscyoc4_fixed = Dict(("A", "Phe") => [(1506.7, 437.3)], ("A", "Tyr") => [(1506.7, 419.0)],
                       ("B", "Tyr") => [(1528.3, 1665.7)], ("C", "Phe") => [(364.5, 2796.0)],
                       ("C", "Tyr") => [(364.5, 3774.0)])
vanscyoc4_panel(name, grid, frame, yticks, letter, crop, CaM; axis = (;)) = PlotPanel(;
    name, grid, frame, crop, margin = 8, frame_width = 900, skip = [letter],
    conditions = (CaM = CaM, CaM_uM = 6.0, temperature_C = 22.0, OGB5N_uM = 0.1, OGB5N_Kd_uM = 29.6),
    x = vanscyoc_x4, y = AxisMap(ticks(yticks, 0.8:-0.2:(1.0 - 0.2 * length(yticks)))), axis)
vanscyoc_fig4 = PaperFigure(;
    paper = "VanScyoc_etal_2002", image = "Fig4.png", doi = "10.1016/S0006-3495(02)75286-7",
    description = "Ca2+ titrations of 6 uM R74A CaM(1-148) (A), R90A CaM(1-148) (B) and R90G CaM(1-148) (C) at 22 C, " *
                  "monitored at ex/em 250/280 nm (filled triangles; Phe, sites I and II) and 277/320 nm (open " *
                  "triangles; Tyr, sites III and IV), each normalized as (f - fmin)/(fmax - fmin); free Ca2+ from 0.1 " *
                  "uM Oregon Green 488 BAPTA-5N (OGB-5N), Kd 29.6 uM (50 mM HEPES pH 7.4, 100 mM KCl, 0.05 mM EGTA, " *
                  "5 mM NTA, 22 C). Solid lines are fits and dashed lines the CaM(1-148) fit, not data. In B the 250/280 " *
                  "series rises and the 277/320 series falls with Ca2+, the reverse of A and C. Five overlapping or " *
                  "frame-line markers were read by eye; B's 277/320 point at 4.2 uM is wholly hidden under the 250/280 " *
                  "one and is given its centre. Image: the journal web page's GIF, converted to PNG",
    xvar = "Ca_free" => "M", yvar = "normalized_fluorescence" => "1", fontsize = 95, linewidth = 8, ticksize = 40,
    panels = [vanscyoc4_panel("A", (1, 1), (361, 45, 3249, 1236), [339.5, 541.5, 741.5, 941.5, 1141.5],
                              (395, 200, 495, 300), (0, 0, 3334, 1290), "R74A(1-148)";
                              axis = (; xticklabelsvisible = false)),
              vanscyoc4_panel("B", (2, 1), (361, 1361, 3249, 2555), [1657.5, 1856.5, 2058.5, 2260.5, 2462.5],
                              (395, 1520, 495, 1630), (0, 1300, 3334, 2600), "R90A(1-148)";
                              axis = (; xticklabelsvisible = false, ylabel = "Normalized Fluorescence Intensity")),
              vanscyoc4_panel("C", (3, 1), (361, 2702, 3249, 3893), [2996.5, 3198.5, 3396.5, 3598.5],
                              (385, 2855, 490, 2965), (0, 2650, 3334, 4238), "R90G(1-148)";
                              axis = (; xlabel = rich("[Calcium] ", subscript("free"))))],
    series = [MarkerSeries(; name, panel, ink = :black, shape = :triangle, radius = 29, hollow, stroke = 8, locate = :peak,
                           threshold = 0.75, fixed = get(vanscyoc4_fixed, (panel, name), NTuple{2,Float64}[]),
                           conditions = (ex_nm = ex, em_nm = em), colour = "#000000", label = "")
              for panel in ("A", "B", "C")
              for (name, hollow, ex, em) in (("Phe", false, 250.0, 280.0), ("Tyr", true, 277.0, 320.0))],
    axis = (; ylabel = "", xticks = (exp10.(-7:-3), [rich("10", superscript(string(e))) for e in -7:-3]),
            xminorticks = [k * exp10(e) for e in -7:-4 for k in 2:9], xminorticksvisible = true,
            yticks = (0:0.2:1, ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])))

## Gaertner et al. 2004, J Biol Chem, doi:10.1074/jbc.M405352200
# Fig 1A: Ca2+ dissociation from CaM (open circles) and CaM + RC3 (small filled dots) at
# 3 °C, monitored by quin-2 fluorescence and converted by the authors to mol Ca2+ released
# per mol CaM; lines are exponential fits. Only the CaM circles are read: the CaM + RC3 dots
# are skipped above 3.3 mol/mol and before 0.004 s, where no circle is, and so is a match on
# the CaM + RC3 fit line beside the circles at 0.004 s, 1.7 mol/mol. The x axis is broken
# after 0.25 s, and the segment past the break has a single tick (1.00 s), so only the circle
# centred on that tick is read, in panel A_break (its x map is that tick at the scale before
# the break; the other circle there, at an unknown time, is skipped).
# Fig 1B: Ca2+ dissociation from CaM (open circles) and CaM + RC3 (small filled dots) at
# 22 °C, monitored by Tyr fluorescence after rapid mixing with EGTA; lines are exponential
# fits. SigmaPlot axes: left and bottom spines, outward ticks. The x axis is broken after
# 0.06 s, and the points past the break (at 0, around 1 s) are not read: that segment has a
# single tick. Both series are black, so CaM + RC3, whose dots merge into a band, is read in
# a second panel over the same axes; each panel skips the other series (CaM circles all
# above 0.5, CaM + RC3 dots below). The CaM circles touch their neighbours and sit on the
# fit line, so centres are the template peaks.
gaertner_x = AxisMap(ticks([140.5, 323.5, 506.5, 689.5], 0:0.02:0.06))
gaertner_y = AxisMap(ticks([746.5, 826.5, 907.5, 987.5, 1068.5, 1148.5], 1.0:-0.2:0.0))
gaertner_A_y = AxisMap(ticks([81.5, 194.5, 307.5, 420.5, 533.5], 4:-1:0))
gaertner_A_conditions = (CaM_uM = 2.0, Ca_total_uM = 100.0, quin2_uM = 150.0, MgCl2_uM = 1000.0, temperature_C = 3.0)
gaertner_A_axis = (; ylabel = rich("Moles Ca", superscript("2+"), " released/Mole CaM"), yticks = 0:4)
gaertner_fig1 = PaperFigure(;
    paper = "Gaertner_etal_2004", image = "Fig1.jpg", doi = "10.1074/jbc.M405352200",
    description = "A: Ca2+ dissociation from CaM alone or with RC3 at 3 C, measured by quin-2 fluorescence (ex " *
                  "334.5 nm, emission through a 435 nm cut-on filter): 2 uM CaM and 100 uM CaCl2 in 20 mM MOPS pH " *
                  "7.0, 100 mM KCl, 1 mM MgCl2, with or without 50 uM RC3, rapidly mixed with 150 uM quin-2 " *
                  "(syringe concentrations, before mixing); averages of five runs. The authors converted the signal " *
                  "to mol Ca2+ released per mol CaM by setting the amplitude of the slow release from CaM alone to " *
                  "2 mol. Lines are double (CaM) or triple (CaM + RC3) exponential fits, not data. Only CaM alone " *
                  "is read; past the x-axis break only the point at the 1.00 s tick (panel A_break). " *
                  "B: Ca2+ dissociation from CaM alone or with RC3 (neurogranin) at 22 C: 2 uM CaM and 100 uM CaCl2, " *
                  "with or without 50 uM RC3, in 20 mM MOPS pH 7.0, 100 mM KCl (no MgCl2 listed), rapidly mixed " *
                  "with 5 mM EGTA so that free Ca2+ after mixing was < 10 nM (figure legend; the Methods give 10 uM " *
                  "CaM and pH 7.5). Concentrations here are the syringe concentrations, before mixing; using them " *
                  "as the concentrations after mixing is a simplification, since diluting all of them alike leaves " *
                  "the scaled C-lobe trace unchanged. " *
                  "Tyr fluorescence (ex 276 nm, emission bandpass at 355 nm; C-lobe), scaled; averages of five runs. " *
                  "Lines are single (CaM) or double (CaM + RC3) exponential fits, not data. The points after the " *
                  "x-axis break (at 0, around 1 s) are not read",
    xvar = "time" => "s", yvar = "Tyr_fluorescence_scaled" => "1", fontsize = 22, linewidth = 2, ticksize = 7,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (138, 70, 736, 533), crop = (0, 0, 737, 635),
                        skip = [(138, 60, 736, 165), (138, 165, 146, 533), (147, 334, 149, 341)],
                        conditions = gaertner_A_conditions,
                        yvar = "Ca_released_per_CaM" => "mol/mol",
                        x = AxisMap(ticks([138.5, 257.5, 377.0, 496.5, 615.5], 0:0.05:0.20)), y = gaertner_A_y,
                        axis = (; gaertner_A_axis...,
                                xticks = (0:0.05:0.25, ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"]))),
              PlotPanel(; name = "A_break", title = "", grid = (1, 2), frame = (744, 70, 806, 533), crop = (737, 0, 832, 635),
                        skip = [(744, 60, 775, 533)], conditions = gaertner_A_conditions,
                        yvar = "Ca_released_per_CaM" => "mol/mol", frame_width = 600 * (806 - 744) / (736 - 138),
                        x = AxisMap(ticks([798.5 - 119.25, 798.5], [0.95, 1.0])), y = gaertner_A_y,
                        axis = (; gaertner_A_axis..., ylabel = "", yticklabelsvisible = false,
                                xticks = ([1.0], ["1.00"]))),
              [PlotPanel(; name, title = "", grid = (2, col), frame = (140, 700, 736, 1169), crop = (0, 640, 832, 1279),
                         skip, x = gaertner_x, y = gaertner_y,
                         conditions = (CaM_uM = 2.0, Ca_total_uM = 100.0, EGTA_uM = 5000.0, temperature_C = 22.0))
               for (name, col, skip) in (("B", 1, [(140, 950, 736, 1169)]), ("B_RC3", 2, [(140, 700, 736, 945)]))]...],
    series = [[MarkerSeries(; name = "CaM", panel, ink = :black, radius = 6, hollow = true, locate = :peak,
                            threshold = 0.5, conditions = (RC3_uM = 0.0,), colour = "#000000", label = "")
               for panel in ("A", "A_break", "B")]...,
              TraceSeries(; name = "CaM + RC3", panel = "B_RC3", ink = :black, radius = 2, conditions = (RC3_uM = 50.0,),
                          colour = "#000000")],
    axis = (; quintana_axis..., xlabel = "Time (s)", ylabel = "Tyrosine Fluorescence, scaled",
            xticks = (0:0.02:0.06, ["0.00", "0.02", "0.04", "0.06"]),
            yticks = (0:0.2:1, ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])))

# Fig 2A: steady-state Ca2+ titration of 30 µM CaM with 30 µM 5,5'-Br2BAPTA, alone (filled
# circles) or with 60 µM RC3 (filled squares), read by Br2BAPTA absorbance; the lines join
# the points. Same SigmaPlot axes as Fig 1, the x axis offset below y = 0. Only CaM alone is
# read. At ~10 px the two filled shapes are too alike for the templates to tell apart, so
# the circles' centres are given: snapped to the circle template where they stand alone,
# kept as read by eye where they overlap a square (the first one and four near saturation).
# Panel B (CaM + RC3 association traces and fitted rates) is not read.
gaertner_fig2_circles = [(125.2, 172.4), (133.5, 197.8), (141.2, 218.9), (149.7, 236.8), (157.7, 252.9), (165.8, 266.2),
                         (173.7, 279.1), (181.5, 290.9), (193.2, 306.9), (204.9, 322.1), (216.2, 335.9), (227.9, 349.3),
                         (239.0, 361.4), (250.7, 373.0), (261.9, 384.9), (272.9, 394.4), (283.9, 403.3), (294.8, 411.3),
                         (313.0, 421.7), (330.5, 430.0), (348.0, 437.7), (365.3, 442.9), (382.5, 447.5), (399.2, 451.0),
                         (414.0, 453.8), (603.1, 462.5)]
gaertner_fig2_overlapped = [(117.2, 139.6), (447.3, 456.0), (479.0, 458.8), (509.0, 460.6), (534.0, 461.5)]
gaertner_fig2 = PaperFigure(;
    paper = "Gaertner_etal_2004", image = "Fig2.jpg", doi = "10.1074/jbc.M405352200",
    description = "Steady-state Ca2+ binding: CaCl2 titrated into 30 uM decalcified CaM and 30 uM 5,5'-Br2BAPTA " *
                  "(Kd 1.6 uM, stated in the paper), with or without 60 uM RC3, in 20 mM MOPS pH 7.5, 100 mM KCl, " *
                  "at room temperature (recorded as 22 C), read by Br2BAPTA absorbance at 263 nm (264 nm in the " *
                  "legend), scaled by the authors. x is the total Ca2+ added (axis label [Ca2+]); free Ca2+ could " *
                  "not leave 30 uM Br2BAPTA unsaturated up to ~150 uM. Lines join the points. Only CaM alone is " *
                  "read; its first point and four near saturation overlap CaM + RC3 squares and were read by eye",
    xvar = "Ca_total" => "uM", yvar = "Br2BAPTA_absorbance_scaled" => "1", fontsize = 22, linewidth = 2, ticksize = 7,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (110, 90, 660, 482), margin = 2,
                        crop = (0, 0, 711, 595), skip = [(345, 196, 575, 286)],
                        conditions = (CaM_uM = 30.0, Br2BAPTA_uM = 30.0, Br2BAPTA_Kd_uM = 1.6, temperature_C = 22.0),
                        x = AxisMap(ticks([109.5, 201.5, 292.5, 384.5, 476.5, 568.0, 659.5], 0:50:300)),
                        y = AxisMap(ticks([97.0, 171.5, 245.5, 319.5, 393.5, 468.0], 1.0:-0.2:0.0)))],
    series = [MarkerSeries(; name = "CaM", panel = "A", ink = :black, shape = :circle, radius = 5,
                           near = gaertner_fig2_circles, fixed = gaertner_fig2_overlapped,
                           conditions = (RC3_uM = 0.0,), colour = "#000000", label = "CaM alone")],
    axis = (; quintana_axis..., xlabel = rich("[Ca", superscript("2+"), "] (µM)"),
            ylabel = rich("Br", subscript("2"), "BAPTA absorbance, scaled"), xticks = 0:50:300,
            yticks = (0:0.2:1, ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])))

# Fig 3A: Ca2+ dissociation at 22 °C into quin-2 from CaM alone, with CaM kinase II, or with
# CaM kinase II autophosphorylated in ATPγS, converted by the authors to mol Ca2+ released
# per mol CaM. Small open circles merging along the fit lines, read as three bands (top to
# bottom: CaM, + CaMKII, + phospho-CaMKII) over 0-1.03 s only, the extent of the CaM trace:
# beyond it the tracker loses CaM and takes stray label ink. The two CaM kinase II traces are
# tracked so the bands are told apart, but not saved. A lone circle at t = 0, 1.75 mol/mol, on
# none of the curves, is skipped. Panels B and C (fit residuals) are not read.
gaertner_fig3 = PaperFigure(;
    paper = "Gaertner_etal_2004", image = "Fig3.jpg", doi = "10.1074/jbc.M405352200",
    description = "Ca2+ dissociation measured by quin-2 fluorescence (ex 334.5 nm, emission through a 435 nm " *
                  "cut-on filter) at 22 C: 2 uM CaM and 20 uM CaCl2 in 20 mM MOPS pH 7.0, 100 mM KCl, 1 mM MgCl2, " *
                  "alone or with 2 uM CaM kinase II or CaM kinase II + 1 mM ATPgammaS, rapidly mixed with 150 uM " *
                  "quin-2 (syringe concentrations, before mixing); averages of five runs. The authors converted the " *
                  "signal to mol Ca2+ released per mol CaM by setting the amplitude of Ca2+ release from CaM alone " *
                  "to 2 mol. Lines are single (CaM) or triple (CaM kinase II) exponential fits, not data. Traces are " *
                  "the centres of the bands the circles form, one point per pixel column, not the individual " *
                  "measurements. Only CaM alone is saved, over 0-1.03 s, the extent of its trace. A lone circle at " *
                  "t = 0 (1.75 mol/mol) on none of the curves is not read",
    xvar = "time" => "s", yvar = "Ca_released_per_CaM" => "mol/mol", fontsize = 22, linewidth = 2, ticksize = 7,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (157, 82, 828, 554), read = (185, 82, 356, 552),
                        crop = (0, 0, 851, 668), skip = [(182, 348, 192, 365)],
                        conditions = (CaM_uM = 2.0, Ca_total_uM = 20.0, quin2_uM = 150.0, MgCl2_uM = 1000.0,
                                      temperature_C = 22.0),
                        x = AxisMap(ticks([187.5, 348.0, 507.5, 667.5, 827.5], 0:4)),
                        y = AxisMap(ticks([104.5, 216.5, 328.5, 440.5, 552.5], 4:-1:0)))],
    series = [TraceSeries(; name, panel = "A", ink = :black, curve, hollow = true, gap = 4, radius = 3, saved = curve == 1,
                          conditions = (; target, target_uM), colour = "#000000", label = name)
              for (curve, name, target, target_uM) in ((1, "CaM", "none", 0.0), (2, "CaM + CaMKII", "CaMKII", 2.0),
                                                       (3, "CaM + phospho-CaMKII", "phospho-CaMKII", 2.0))],
    axis = (; quintana_axis..., xlabel = "Time (s)", ylabel = rich("Moles Ca", superscript("2+"), " released/Mole CaM"),
            xticks = 0:4, yticks = 0:4))

## Johnson et al. 1996, J Biol Chem, doi:10.1074/jbc.271.2.761
# Fig 1A: Ca2+ dissociation from CaM into quin-2 at 10 °C, over 0-2 s (main) and 0-20 ms
# (inset). Boxed frames with inward ticks, bold labels. Data and fit are both black lines
# lying on each other, so each trace is read as the centre of their band; in the inset the
# fit separates only before ~2.5 ms, where it drops out of the frame below the data. The
# first 1.6 ms of each trace is premixing (the dead time). Labels inside the frames are
# skipped, and so is the inset in the main panel.
johnson_fig1 = PaperFigure(;
    paper = "Johnson_etal_1996", image = "Fig1.jpg", doi = "10.1074/jbc.271.2.761",
    description = "Ca2+ dissociation from bovine brain CaM measured by quin-2 fluorescence (ex 330 nm, em 510 nm): " *
                  "8 uM CaM + 60 uM Ca2+ in 20 mM HEPES pH 7.0 (no KCl or Mg2+ listed) rapidly mixed with an equal " *
                  "volume of 150 uM quin-2 in the same buffer at 10 C. Concentrations here are the syringe " *
                  "concentrations, before mixing. Main panel: 0-2 s, average of four traces; inset: 0-20 ms, average " *
                  "of six traces. Traces were triggered at time zero and the first 1.6 ms (premixing, the dead time) " *
                  "is shown; the authors fitted only after mixing was complete (fast 405/s, slow 2.4/s). Data and fit " *
                  "lines are both black and overlap, so the points are the centre of their joint band. Panel B " *
                  "(CaM with MLCK or peptides) is not read",
    xvar = "time" => "s", yvar = "quin2_fluorescence_change" => "a.u.", fontsize = 32, linewidth = 4, ticksize = 6,
    panels = [PlotPanel(; name = "A", grid = (1, 1), frame = (124, 11, 1190, 770), margin = 8, crop = (0, 0, 1225, 880),
                        skip = [(130, 12, 190, 70), (140, 500, 290, 600), (730, 128, 985, 180), (515, 262, 1186, 766)],
                        x = AxisMap(ticks([341.5, 555.0, 768.0, 983.0], 0.4:0.4:1.6)),
                        y = AxisMap(ticks([117.0, 231.0, 350.5, 469.5, 589.0, 704.0], 4.2:-0.2:3.2)),
                        axis = (; xticks = (0.4:0.4:2.0, ["0.4", "0.8", "1.2", "1.6", "2.0"]),
                                yticks = (3.2:0.2:4.2, ["3.2", "3.4", "3.6", "3.8", "4.0", "4.2"]))),
              PlotPanel(; name = "A_inset", title = "", grid = (1, 2), frame = (527, 268, 1171, 718), margin = 8,
                        crop = (515, 260, 1180, 760), fontsize = 22, skip = [(785, 400, 940, 495)],
                        x = AxisMap(ticks([656.0, 786.0, 917.0, 1045.0], 0.004:0.004:0.016)),
                        y = AxisMap(ticks([346.0, 460.5, 573.0, 687.5], 3.75:-0.15:3.30)),
                        axis = (; xlabel = "", ylabel = "",
                                xticks = (0.004:0.004:0.016, ["4.0", "8.0", "12.0", "16.0msec"]),
                                yticks = (3.30:0.15:3.75, ["3.30", "3.45", "3.60", "3.75"])))],
    series = [TraceSeries(; name = "CaM", panel, ink = :black, gap = 3, radius = 1,
                          conditions = (CaM_uM = 8.0, Ca_total_uM = 60.0, quin2_uM = 150.0, temperature_C = 10.0),
                          colour = "#000000")
              for panel in ("A", "A_inset")],
    axis = (; xlabel = "Time (s)", ylabel = "ΔQuin Fluorescence", xticklabelfont = :bold, yticklabelfont = :bold,
            xlabelfont = :bold, ylabelfont = :bold, xtickcolor = :black, ytickcolor = :black))

# Fig 2: Ca2+ dissociation at 10 °C. Above y = 0, quin-2 fluorescence as Ca2+ leaves the
# C-terminal sites of CaM41/75, alone or with MARCKS peptide; below, CaM Tyr fluorescence as
# EGTA strips CaM, alone or with MARCKS. Dense noisy traces, data and fit lines overlapping,
# each shifted to start at 0. The quin-2 and Tyr traces are different quantities, so they
# are read as two panels over the same axes, split at y = 0, over 0-5.2 s only (the traces
# without MARCKS end at 5 s). The + MARCKS traces in that range are tracked so the bands are
# told apart, but not saved. Labels are skipped.
johnson_x2 = AxisMap(ticks([506.0, 818.5, 1128.5, 1441.5], 4:4:16))
johnson_y2 = AxisMap(ticks([233.0, 407.5, 580.5, 753.5, 934.0, 1107.5], 0.4:-0.2:-0.6))
johnson_fig2 = PaperFigure(;
    paper = "Johnson_etal_1996", image = "Fig2.jpg", doi = "10.1074/jbc.271.2.761",
    description = "Ca2+ dissociation at 10 C in 20 mM HEPES pH 7.0. quin2: 4 uM CaM41/75 + 60 uM Ca2+ rapidly mixed " *
                  "with an equal volume of 150 uM quin-2; quin-2 fluorescence as in Fig 1. Tyr: 4 uM CaM + 100 uM " *
                  "Ca2+ rapidly mixed with an equal volume of 5 mM EGTA; CaM Tyr fluorescence (ex 275 nm, UG1 " *
                  "filter). Concentrations are the syringe concentrations, before mixing. Each trace is an average of " *
                  "five, shifted by the authors to start at 0; data and fit lines overlap, so the points are the centre " *
                  "of their joint band, one per pixel column. Only 0-5.2 s is read (these traces end at 5 s), and the " *
                  "first ~0.1 s lies under the frame line. The traces with 10 uM MARCKS peptide are not saved",
    xvar = "time" => "s", fontsize = 40, linewidth = 7, ticksize = 10,
    panels = [PlotPanel(; name = "quin2", title = "", grid = (1, 1), frame = (199, 63, 1753, 1197), read = (203, 66, 598, 590),
                        crop = (0, 0, 1799, 1348), skip = [(225, 68, 560, 132), (416, 368, 805, 415)],
                        yvar = "quin2_fluorescence_change" => "a.u.",
                        conditions = (Ca_total_uM = 60.0, quin2_uM = 150.0, temperature_C = 10.0),
                        x = johnson_x2, y = johnson_y2),
              PlotPanel(; name = "Tyr", title = "", grid = (1, 2), frame = (199, 63, 1753, 1197), read = (203, 590, 598, 1195),
                        crop = (0, 0, 1799, 1348), skip = [(416, 820, 805, 866), (238, 1108, 550, 1170)],
                        yvar = "Tyr_fluorescence_change" => "a.u.",
                        conditions = (Ca_total_uM = 100.0, EGTA_uM = 5000.0, temperature_C = 10.0),
                        x = johnson_x2, y = johnson_y2)],
    series = [TraceSeries(; name, panel, ink = :black, curve, gap = 6, radius = 1, saved = MARCKS_uM == 0,
                          conditions = (; CaM, CaM_uM = 4.0, MARCKS_uM), colour = "#000000")
              for (panel, curve, name, CaM, MARCKS_uM) in (("quin2", 1, "CaM41/75", "CaM41/75", 0.0),
                                                           ("quin2", 2, "CaM41/75 + MARCKS", "CaM41/75", 10.0),
                                                           ("Tyr", 1, "CaM + MARCKS", "CaM", 10.0),
                                                           ("Tyr", 2, "CaM", "CaM", 0.0))],
    axis = (; xlabel = "Time (s)", ylabel = "ΔQuin2 or ΔCaM Tyrosine Fluorescence",
            xticks = (4:4:20, ["4.0", "8.0", "12.0", "16.0", "20.0"]),
            yticks = (-0.6:0.2:0.4, ["-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4"]),
            xticklabelfont = :bold, yticklabelfont = :bold, xlabelfont = :bold, ylabelfont = :bold,
            xtickcolor = :black, ytickcolor = :black))

# Fig 3: C-lobe Tyr fluorescence of CaM (open circles) and CaM41/75 (asterisks) against log
# free [Ca2+] in EGTA buffer at 10 °C, normalized to its full increase; the RS-20, C2K and
# MARCKS series are not read, nor are the S.E. bars. Every symbol overlaps others (the CaM
# and CaM41/75 curves coincide, and all five series meet at both plateaus), and there is no
# asterisk template, so all centres were read by eye from magnified crops, at the titration
# columns; template detection is skipped over the whole frame. Plateau points are ±2%.
johnson_fig3_circles = [(791.5, 1056.0), (972.3, 1050.0), (1063.0, 1035.0), (1153.5, 904.0), (1244.5, 642.5),
                        (1335.3, 289.3), (1426.0, 113.0), (1516.5, 97.0), (1697.5, 101.0)]
johnson_fig3_asterisks = [(972.3, 1018.5), (1063.0, 1005.0), (1153.5, 928.0), (1244.5, 656.0), (1335.3, 246.7),
                          (1426.0, 82.0), (1516.5, 97.0), (1697.5, 101.0)]
johnson_fig3 = PaperFigure(;
    paper = "Johnson_etal_1996", image = "Fig3.jpg", doi = "10.1074/jbc.271.2.761",
    description = "Ca2+ dependence of the C-lobe Tyr fluorescence (ex 275 nm, em 305 nm) of 1 uM CaM or CaM41/75 " *
                  "titrated with Ca2+ in 200 mM HEPES, 2 mM EGTA, pH 7.0, at 10 C, against log free [Ca2+] (M) as " *
                  "set by the EGTA buffer; 100% is the full increase (~2.4-fold). Means of three titrations; S.E. " *
                  "bars not read; lines are sigmoidal fits. The RS-20, C2K and MARCKS series are not read. All " *
                  "centres were read by eye (symbols overlap; asterisks have no template); plateau points are +-2%",
    xvar = "log10_Ca_free" => "log10(M)", yvar = "Tyr_fluorescence_increase" => "%", fontsize = 40, linewidth = 6,
    ticksize = 14,
    panels = [PlotPanel(; name = "main", title = "", grid = (1, 1), frame = (247, 16, 1786, 1138), crop = (0, 0, 1799, 1345),
                        skip = [(247, 16, 1786, 1138)], conditions = (CaM_uM = 1.0, EGTA_uM = 2000.0, temperature_C = 10.0),
                        x = AxisMap(ticks([336.5, 791.5, 1244.5, 1697.5], -8.0:1.0:-5.0)),
                        y = AxisMap(ticks([93.0, 576.0, 1059.5], [100.0, 50.0, 0.0])))],
    series = [MarkerSeries(; name = "CaM", panel = "main", ink = :black, shape = :circle, radius = 12, hollow = true,
                           fixed = johnson_fig3_circles, conditions = (CaM = "CaM",), colour = "#000000", label = "CaM"),
              MarkerSeries(; name = "CaM41/75", panel = "main", ink = :black, shape = :diamond, radius = 12,
                           fixed = johnson_fig3_asterisks, conditions = (CaM = "CaM41/75",), colour = "#808080",
                           label = "CaM41/75 (asterisks)")],
    axis = (; xlabel = "Log [Ca]", ylabel = "% Increase Fluorescence", xticks = (-8:1:-5, ["-8.0", "-7.0", "-6.0", "-5.0"]),
            yticks = 0:50:100,
            xticklabelfont = :bold, yticklabelfont = :bold, xlabelfont = :bold, ylabelfont = :bold,
            xtickcolor = :black, ytickcolor = :black))

# Fig 6: CaM Tyr fluorescence at 10 °C, three noisy traces with fits on top: Ca2+-saturated CaM
# mixed with Ca2+ (the saturated level), Ca2+-saturated CaM mixed with EGTA, and CaM in EGTA
# mixed with Ca2+ (a rapid Ca2+ transient that only partly fills the C-lobe). Read as three
# bands (top to bottom at 0.5 s); the last two merge after ~1 s and are flagged there. Labels
# are skipped, and the frame is read from 0.09 s and above y = 4.0 to keep the inward tick
# marks out (a tick next to the merged band would pass for a third band). Gaps of up to 70 px
# are bridged: the noisy data line strays up to ~0.03 (~65 px) from its fit line, which would
# otherwise split a trace into two bands.
johnson_fig6 = PaperFigure(;
    paper = "Johnson_etal_1996", image = "Fig6.jpg", doi = "10.1074/jbc.271.2.761",
    description = "CaM Tyr fluorescence (as in Fig 3) after rapid mixing in 200 mM HEPES pH 7.0 at 10 C: 4 uM CaM + " *
                  "100 uM Ca2+ mixed with 100 uM Ca2+ (Ca2+-saturated level); 4 uM CaM + 200 uM Ca2+ mixed with 2 mM " *
                  "EGTA (dissociation from saturated C-lobe sites); 4 uM CaM + 2 mM EGTA mixed with 200 uM Ca2+ " *
                  "(transient occupancy of the C-lobe sites by the rapid Ca2+ transient of Fig 4). Concentrations " *
                  "are the syringe concentrations, before mixing (equal volumes). Averages of seven traces; data and " *
                  "fit lines overlap, so the points are the centre of their joint band, one per pixel column; the two " *
                  "EGTA traces merge after ~1 s (flagged merged). Read from 0.09 s (inward tick marks before it)",
    xvar = "time" => "s", yvar = "Tyr_fluorescence" => "a.u.", fontsize = 40, linewidth = 6, ticksize = 12,
    panels = [PlotPanel(; name = "main", title = "", grid = (1, 1), frame = (196, 23, 1750, 1050), read = (218, 23, 1740, 1015),
                        crop = (0, 0, 1799, 1212),
                        skip = [(680, 50, 1170, 112), (640, 612, 1230, 748), (635, 902, 1225, 1018)],
                        conditions = (CaM_uM = 4.0, temperature_C = 10.0),
                        x = AxisMap(ticks([503.5, 816.0, 1125.0, 1437.5], 1.0:1.0:4.0)),
                        y = AxisMap(ticks([127.5, 354.0, 571.0, 790.5, 1010.5], 4.4:-0.1:4.0)))],
    series = [TraceSeries(; name, panel = "main", ink = :black, curve, gap = 70, radius = 1,
                          conditions = (; Ca_uM, EGTA_uM, mixed_with_Ca_uM, mixed_with_EGTA_uM), colour = "#000000")
              for (curve, name, Ca_uM, EGTA_uM, mixed_with_Ca_uM, mixed_with_EGTA_uM) in (
                  (1, "CaM + Ca vs Ca", 100.0, 0.0, 100.0, 0.0), (2, "CaM + Ca vs EGTA", 200.0, 0.0, 0.0, 2000.0),
                  (3, "CaM + EGTA vs Ca", 0.0, 2000.0, 200.0, 0.0))],
    axis = (; xlabel = "Time (s)", ylabel = "CaM Tyr Fluorescence", xticks = (1:5, ["1.0", "2.0", "3.0", "4.0", "5.0"]),
            yticks = (4.0:0.1:4.4, ["4.0", "4.1", "4.2", "4.3", "4.4"]),
            xticklabelfont = :bold, yticklabelfont = :bold, xlabelfont = :bold, ylabelfont = :bold,
            xtickcolor = :black, ytickcolor = :black))

## Putkey et al. 2008, J Biol Chem, doi:10.1074/jbc.M703831200
# Fig 5: Ca2+ binding to recombinant CaM alone and with PEP-19 polypeptides; only CaM alone is
# saved. Image: the journal's high-resolution figure (1351 x 1800). Boxed frames, 3-4 px lines,
# inward ticks. A: Br2BAPTA titration, filled circles (CaM), open squares (+ PEP-19), open
# triangles (+ PEP(28-62)), filled diamonds (+ PEP(39-62)); lines are fits. The other three
# series are matched only as rivals, so each marker goes to the template it fits best, up to
# ~105 uM. The 0.4 tick and one diamond the rivals lost are skipped. Beyond ~105 uM the CaM
# circles lie under the + PEP-19 squares and + PEP(28-62) triangles (and from ~135 uM the
# + PEP(39-62) diamonds too), so there the centre of each marker cluster is given as the
# CaM point: the density peak of the black ink blurred over a marker width (Gaussian, 4 px),
# with the frame lines and bottom tick marks blanked; the cluster at the 200 uM tick takes
# its x from the image, as blanking the tick shifts its peak. B:
# quin-2 dissociation traces; the + PEP-19 and + PEP(28-62) traces (top, to 0.55 s) and
# + PEP(39-62) (bottom) are tracked so the bands are told apart, but not saved; the top ones
# only up to 0.545 s, where they end and CaM rises into their band. The inward tick marks on
# the left and bottom frame lines are skipped (they would pass for bands among the early
# data). D: C-lobe association traces of apo-CaM in Ca2+/Br2BAPTA buffers, top to bottom 2.9,
# 2.3, 1.6 and 0.9 uM free Ca2+ (larger and faster at higher Ca2+). E: as D at 3 uM free Ca2+,
# CaM alone (bottom) and the + PEP-19 / + PEP(28-62) band (top, tracked, not saved), read up
# to 0.18 s before CaM merges into it. C (EGTA controls, no CaM) and the insets (fitted rates)
# are not read; labels and panel letters are skipped.
putkey08_axis = (; xtickalign = 1, ytickalign = 1, xtickcolor = :black, ytickcolor = :black)
putkey08_B_x = AxisMap(ticks([223.5, 323.0, 423.5, 526.0], 0.2:0.2:0.8))
putkey08_B_y = AxisMap(ticks([717.0, 851.0, 918.0, 987.0, 1054.0], [0.30, 0.20, 0.15, 0.10, 0.05]))
putkey08_B_skip = [(130, 676, 270, 708), (160, 736, 255, 768), (250, 815, 305, 848), (465, 858, 605, 890),
                   (398, 852, 622, 1110), (185, 625, 235, 672),
                   [(121, y - 3, 135, y + 4) for y in (717, 784, 851, 918, 987, 1054)]...,
                   [(x - 3, 1108, x + 3, 1121) for x in (223, 323, 423, 526)]...]
putkey08_B_conditions = (CaM_uM = 2.0, Ca_total_uM = 20.0, quin2_uM = 300.0, temperature_C = 23.0)
putkey08_A_conditions = (CaM_uM = 30.0, Br2BAPTA_uM = 30.0, Br2BAPTA_Kd_uM = 1.59, temperature_C = 23.0)
putkey08_A_x = AxisMap(ticks([224.0, 325.5, 427.5], 50:50:150))
putkey08_A_y = AxisMap(ticks([112.5, 214.0, 314.5, 416.5], 0.8:-0.2:0.2))
putkey08_A_clusters = [(341.0, 449.0), (360.0, 461.0), (378.0, 471.0), (396.0, 478.0), (412.0, 485.0), (427.0, 489.0),
                       (460.0, 497.0), (492.0, 504.0), (524.0, 507.0), (553.0, 509.0), (582.0, 510.0), (609.0, 510.0)]
putkey08_fig5 = PaperFigure(;
    paper = "Putkey_etal_2008", image = "Fig5.jpg", doi = "10.1074/jbc.M703831200",
    description = "Recombinant CaM in 20 mM MOPS pH 7.5, 100 mM KCl. A: equilibrium Ca2+ titration of 30 uM CaM " *
                  "with 30 uM 5,5'-Br2BAPTA (Kd 1.59 uM), read by Br2BAPTA absorbance at 263 nm relative to its " *
                  "Ca2+-free value, against total Ca2+ (axis labelled mM, but the 30 uM concentrations make it uM); " *
                  "temperature not stated (recorded as 23 C, the stopped-flow temperature). Lines are fits. B: " *
                  "stopped-flow Ca2+ dissociation at 23 C, 2 uM CaM + 20 uM Ca2+ rapidly mixed with 300 uM quin-2 " *
                  "(syringe concentrations, 1:1 mixing; dead time 1.7 ms), quin-2 fluorescence. D, E: stopped-flow " *
                  "association at 20 C: 2 uM apo-CaM with 5 uM BAPTA (syringe A) rapidly mixed with 250 uM Br2BAPTA " *
                  "plus Ca2+ (syringe B; 1:1 mixing), free Ca2+ in the optical chamber 0.9, 1.6, 2.3 and 2.9 uM (D; " *
                  "as stated in the D legend) or 3 uM (E), read by the change in Br2BAPTA absorbance as Ca2+ moves " *
                  "to CaM. Traces are the centres of the bands their dots form, one point per pixel column, not the " *
                  "individual measurements. Only CaM alone is saved. In A, beyond ~105 uM total Ca2+, " *
                  "the CaM circles are hidden under the other series' markers, so each point is the centre " *
                  "of a marker cluster (CaM with + PEP-19 and + PEP(28-62), and from ~135 uM also + PEP(39-62)), " *
                  "not a CaM circle itself",
    xvar = "time" => "s", fontsize = 30, linewidth = 3, ticksize = 8,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (120, 13, 628, 520), read = (120, 13, 335, 520),
                        crop = (0, 0, 650, 590),
                        skip = [(250, 20, 600, 205), (188, 30, 235, 75), (120, 305, 137, 324), (134, 68, 148, 83)],
                        conditions = putkey08_A_conditions,
                        xvar = "Ca_total" => "uM", yvar = "Br2BAPTA_absorbance_relative" => "1",
                        x = putkey08_A_x, y = putkey08_A_y, axis = (; xticks = 50:50:250, yticks = 0.2:0.2:1.0)),
              PlotPanel(; name = "B", title = "", grid = (2, 1), frame = (121, 609, 625, 1123),
                        crop = (0, 590, 650, 1230), skip = putkey08_B_skip, conditions = putkey08_B_conditions,
                        yvar = "quin2_fluorescence_change" => "a.u.", x = putkey08_B_x, y = putkey08_B_y,
                        axis = (; xticks = 0:0.2:1, yticks = 0.05:0.05:0.30)),
              PlotPanel(; name = "D", title = "", grid = (3, 1), frame = (811, 606, 1324, 1121), crop = (680, 590, 1350, 1170),
                        skip = [(825, 620, 875, 670)],
                        conditions = (CaM_uM = 2.0, BAPTA_uM = 5.0, Br2BAPTA_uM = 250.0, temperature_C = 20.0),
                        yvar = "absorbance_change" => "AU",
                        x = AxisMap(ticks([897.5, 982.5, 1066.5, 1153.5, 1237.5], 0.05:0.05:0.25)),
                        y = AxisMap(ticks([709.0, 812.5, 914.0, 1018.5], 0.020:-0.005:0.005)),
                        axis = (; xticks = 0.05:0.05:0.25, yticks = 0:0.005:0.025)),
              PlotPanel(; name = "E", title = "", grid = (3, 2), frame = (810, 1200, 1322, 1716), read = (810, 1200, 1270, 1716),
                        crop = (680, 1180, 1350, 1799),
                        skip = [(825, 1215, 870, 1265), (905, 1212, 995, 1245), (880, 1300, 1010, 1335),
                                (885, 1425, 940, 1455), (1035, 1410, 1310, 1695)],
                        conditions = (CaM_uM = 2.0, BAPTA_uM = 5.0, Br2BAPTA_uM = 250.0, Ca_free_uM = 3.0,
                                      temperature_C = 20.0),
                        yvar = "absorbance_change" => "AU",
                        x = AxisMap(ticks([936.5, 1067.5, 1197.5], 0.05:0.05:0.15)),
                        y = AxisMap(ticks([1246.5, 1364.0, 1482.0, 1599.0], 0.020:-0.005:0.005)),
                        axis = (; xticks = 0:0.05:0.20, yticks = 0:0.005:0.020))],
    series = [MarkerSeries(; name = "CaM", panel = "A", ink = :black, shape = :circle, radius = 6, locate = :peak,
                           fixed = putkey08_A_clusters, conditions = (target = "none",), colour = "#000000", label = "CaM"),
              [MarkerSeries(; name, panel = "A", ink = :black, shape, radius, hollow, locate = :peak, saved = false,
                            colour = "#000000")
               for (name, shape, radius, hollow) in (("CaM + PEP-19", :square, 7, true),
                                                     ("CaM + PEP(28-62)", :triangle, 8, true),
                                                     ("CaM + PEP(39-62)", :diamond, 8, false))]...,
              [TraceSeries(; name, panel = "B", ink = :black, curve, gap = 8, radius = 1, saved = name == "CaM", span,
                           conditions = (target = "none",), colour = "#000000")
               for (curve, name, span) in ((1, "CaM + PEP-19/PEP(28-62)", (-0.01, 0.545)), (2, "CaM", nothing),
                                           (3, "CaM + PEP(39-62)", nothing))]...,
              [TraceSeries(; name, panel = "D", ink = :black, curve, gap = 12, radius = 1,
                           conditions = (target = "none", Ca_free_uM), colour = "#000000",
                           label = replace(name, "uM" => "µM") * " free Ca²⁺", label_below = curve <= 2)
               for (curve, name, Ca_free_uM) in ((1, "2.9 uM", 2.9), (2, "2.3 uM", 2.3), (3, "1.6 uM", 1.6),
                                                 (4, "0.9 uM", 0.9))]...,
              [TraceSeries(; name, panel = "E", ink = :black, curve, gap = 8, radius = 1, saved = name == "CaM",
                           conditions = (target = "none",), colour = "#000000")
               for (curve, name) in ((1, "CaM + PEP-19/PEP(28-62)"), (2, "CaM"))]...],
    axis = putkey08_axis)

## Putkey et al. 2003, J Biol Chem, doi:10.1074/jbc.C300372200
# Fig 1A,B: CaM alone (open circles) and with PEP-19 (filled circles); only CaM alone is saved.
# Image: the journal's figure (459 x 1280), thin grey frames with inward ticks. A: quin-2
# dissociation converted by the authors to mol Ca2+ released per mol CaM, x axis broken after
# 0.02 s and resumed at 0.45 s: the dense open circles before the break are read as a band
# (A, below the + PEP-19 trace), the separate ones after it as markers (A_break), with the
# filled + PEP-19 circles just above as rivals only; the inward x ticks are skipped. B:
# competitive Br2BAPTA titration; conditions and the y quantity (relative Br2BAPTA absorbance,
# unlabelled in the figure; its 0.2 tick is printed "2") are from the online supplement. The
# open circles touch or sit on the filled + PEP-19 ones, so their centres are given, each
# snapped to the best open-circle match. The insets (Tyr fluorescence), C (CaM-CaMKII),
# labels and legends are not read.
putkey03_B_circles = [(53.0, 483.0), (61.25, 512.5), (70.0, 535.75), (77.5, 554.25), (86.75, 569.5), (95.5, 584.0),
                      (104.25, 595.0), (112.0, 605.75), (120.5, 615.25), (133.25, 628.25), (145.0, 640.5),
                      (158.25, 651.25), (170.0, 662.5), (182.5, 672.0), (193.75, 681.25), (205.5, 690.5), (218.0, 698.75),
                      (228.75, 706.25), (241.25, 715.0), (259.5, 726.0), (278.0, 734.5), (297.0, 743.25), (315.5, 750.0),
                      (333.25, 754.25), (351.75, 759.5), (369.5, 762.5), (403.75, 768.25)]
putkey03_A_y = AxisMap(ticks([58.0, 132.5, 207.0, 281.5], 2.0:-0.5:0.5))
putkey03_A_conditions = (CaM_uM = 2.0, Ca_total_uM = 20.0, quin2_uM = 300.0, temperature_C = 23.0)
putkey03_fig1 = PaperFigure(;
    paper = "Putkey_etal_2003", image = "Fig1.jpg", doi = "10.1074/jbc.C300372200",
    description = "Recombinant CaM in 20 mM MOPS pH 7.5, 100 mM KCl, alone (open circles) or with PEP-19. A: " *
                  "stopped-flow Ca2+ dissociation at 23 C (Applied Photophysics SX.18MV, dead time 1.7 ms), 2 uM CaM, " *
                  "20 uM Ca2+ and 300 uM quin-2 (not stated whether syringe or final concentrations; the same lab's " *
                  "2008 paper gives these as syringe concentrations, 1:1 mixing), quin-2 fluorescence converted by " *
                  "the authors to mol Ca2+ released per mol CaM. The x axis is broken: A is 0-0.02 s, A_break " *
                  "0.45-0.55 s. Before the break the dense open circles are read as the centre of their band, one " *
                  "point per pixel column, not the individual measurements. B: equilibrium Ca2+ binding by the " *
                  "competitive assay of Linse et al. (online supplement): 30 uM CaM and 30 uM 5,5'-Br2BAPTA (Kd 1.59 " *
                  "uM, measured at the assay's pH and ionic strength), titrated with a 1.58 mM Ca2+ stock made 30 uM " *
                  "in both CaM and Br2BAPTA, so only Ca2+ varies; Br2BAPTA absorbance at 262 nm relative to its " *
                  "Ca2+-free value (y axis unlabelled in the figure, its 0.2 tick printed '2'), against total Ca2+. " *
                  "The authors' fits include a factor F (typically > 0.85) for errors in the CaM concentration; the " *
                  "stated 30 uM is recorded. Temperature not stated (recorded as 23 C, the lab's room temperature). " *
                  "Lines are fits. Only CaM alone is saved",
    xvar = "time" => "s", yvar = "Ca_released_per_CaM" => "mol/mol", fontsize = 12, linewidth = 1, ticksize = 3,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (57, 8, 250, 356), read = (57, 228, 250, 356),
                        crop = (0, 0, 255, 370), skip = [(118, 268, 202, 292), (147, 345, 158, 357), (241, 345, 251, 357)],
                        conditions = putkey03_A_conditions,
                        x = AxisMap(ticks([152.5, 246.5], [0.01, 0.02])), y = putkey03_A_y,
                        axis = (; xticks = 0:0.01:0.02, yticks = 0:0.5:2.0)),
              PlotPanel(; name = "A_break", title = "", grid = (1, 2), frame = (270, 8, 456, 356), read = (270, 12, 456, 64),
                        crop = (262, 0, 458, 370), conditions = putkey03_A_conditions,
                        frame_width = 600 * (456 - 270) / (250 - 57),
                        x = AxisMap(ticks([299.0, 377.0], [0.45, 0.50])), y = putkey03_A_y,
                        axis = (; xticks = [0.45, 0.50], yticks = 0:0.5:2.0)),
              PlotPanel(; name = "B", title = "", grid = (2, 1), frame = (43, 439, 439, 784), crop = (0, 430, 458, 830),
                        skip = [(63, 712, 178, 768), (160, 455, 418, 640), (195, 640, 418, 670), (78, 446, 112, 480)],
                        conditions = (CaM_uM = 30.0, Br2BAPTA_uM = 30.0, Br2BAPTA_Kd_uM = 1.59, temperature_C = 23.0),
                        xvar = "Ca_total" => "uM", yvar = "Br2BAPTA_absorbance_relative" => "1",
                        x = AxisMap(ticks([143.0, 242.0, 340.0], 50:50:150)),
                        y = AxisMap(ticks([508.0, 578.0, 646.0, 715.0], 0.8:-0.2:0.2)),
                        axis = (; xticks = 0:50:200, yticks = 0.2:0.2:0.8))],
    series = [TraceSeries(; name = "CaM", panel = "A", ink = :grey, hollow = true, gap = 4, radius = 1,
                          conditions = (target = "none",), colour = "#000000"),
              MarkerSeries(; name = "CaM", panel = "A_break", ink = :grey, radius = 3, hollow = true, locate = :peak,
                           conditions = (target = "none",), colour = "#000000", label = ""),
              MarkerSeries(; name = "CaM + PEP-19", panel = "A_break", ink = :grey, radius = 3, locate = :peak,
                           saved = false, colour = "#000000"),
              MarkerSeries(; name = "CaM", panel = "B", ink = :grey, radius = 4, hollow = true, near = putkey03_B_circles,
                           conditions = (target = "none",), colour = "#000000", label = "")],
    axis = putkey08_axis)

## Hoffman et al. 2014, J Biol Chem, doi:10.1074/jbc.M114.560656
# Fig 2C: Ca2+ dissociation from CaM (black filled circles) and CaM/Ng (red) by quin-2
# fluorescence; only CaM is read. Image: the vector figure rendered from the PDF at 600 dpi
# (page 5, Figure 2 region). Origin axes: left and bottom spines, outward ticks. The black line
# is a single-exponential fit; the legend is skipped.
hoffman_fig2 = PaperFigure(;
    paper = "Hoffman_etal_2014", image = "Fig2.png", doi = "10.1074/jbc.M114.560656",
    description = "Stopped-flow Ca2+ dissociation from recombinant mammalian CaM at room temperature (24 C; Applied " *
                  "Photophysics SX.17MV, dead time 1.7 ms) in 50 mM MOPS pH 7.2, 100 mM KCl: Ca2+-loaded CaM rapidly " *
                  "mixed 1:1 with quin-2, final concentrations 5 uM CaM, 25 uM Ca2+ and 150 uM quin-2; quin-2 " *
                  "fluorescence (ex 334.5 nm, 450 nm cut-on filter). Every 50th data point is shown. The line is a " *
                  "single-exponential fit. Only CaM alone is read (CaM/Ng, red, is not). Image rendered from the " *
                  "PDF's vector figure at 600 dpi",
    xvar = "time" => "s", yvar = "quin2_fluorescence" => "a.u.", fontsize = 60, linewidth = 4, ticksize = 15,
    panels = [PlotPanel(; name = "C", title = "", grid = (1, 1), frame = (3256, 290, 4365, 1196), crop = (3150, 250, 4450, 1400),
                        skip = [(4000, 1005, 4310, 1145)],
                        conditions = (CaM_uM = 5.0, Ca_total_uM = 25.0, quin2_uM = 150.0, temperature_C = 24.0),
                        x = AxisMap(ticks([3341.5, 3512.5, 3682.5, 3852.5, 4023.5, 4193.5, 4364.5], 0.0:0.1:0.6)),
                        y = AxisMap(ticks([303.5, 415.5, 526.5, 637.5, 749.5, 860.5, 972.5, 1083.5], 8.6:-0.2:7.2)),
                        axis = (; xticks = 0:0.1:0.6, yticks = 7.0:0.2:8.6))],
    series = [MarkerSeries(; name = "CaM", panel = "C", ink = :black, radius = 11, locate = :peak,
                           conditions = (target = "none",), colour = "#000000", label = "")],
    axis = (; xtickalign = 0, ytickalign = 0))

## Masada et al. 2012, Biochemistry, doi:10.1021/bi300646y
# Fig 5: Ca2+ dissociation from CaM alone by quin-2 fluorescence, WT (A), CaM12 (B, N-lobe
# sites dead) and CaM34 (C, C-lobe sites dead). Image: the 300 ppi greyscale figure from the
# PDF (875 x 2067). Grey data under a black single-exponential fit line, read as the centre
# of the grey-plus-black band per pixel column. C is read from its 0-0.04 s inset only; the
# main C panel is a flat plateau after 0.04 s. Panel labels are skipped.
masada_x = AxisMap(ticks([143.0, 284.5, 427.0, 568.0, 709.5, 852.0], 0.0:0.2:1.0))
masada_conditions = (CaM_uM = 2.0, Ca_total_uM = 20.0, quin2_uM = 150.0, MgCl2_uM = 1000.0, temperature_C = 22.0)
masada_fig5 = PaperFigure(;
    paper = "Masada_etal_2012", image = "Fig5.png", doi = "10.1021/bi300646y",
    description = "Stopped-flow Ca2+ dissociation at 22 C (Applied Photophysics SX.20MV, dead time 5 ms) from " *
                  "recombinant N-terminally His-tagged rat CaM in 20 mM MOPS pH 7.4, 100 mM KCl, 1 mM MgCl2: 2 uM " *
                  "CaM with 20 uM CaCl2 rapidly mixed 1:1 with 150 uM quin-2 (syringe concentrations), quin-2 " *
                  "fluorescence relative to its final value; averages of 5-8 traces. CaM12 and CaM34 carry Asp to " *
                  "Ala in the Ca2+-binding loops of EF-hands 1-2 and 3-4. The black lines are single-exponential " *
                  "fits, not data; points are the centre of the band the grey data and the fit line form, one per " *
                  "pixel column, not the individual measurements. In the C inset the rising phase of the grey " *
                  "data is mostly hidden under the fit line, so there the points largely follow the fit",
    xvar = "time" => "s", yvar = "quin2_fluorescence_relative" => "1", fontsize = 24, linewidth = 3, ticksize = 7,
    panels = [PlotPanel(; name = "A", title = "", grid = (1, 1), frame = (143, 86, 855, 563), margin = 5,
                        crop = (0, 0, 874, 690), skip = [(705, 478, 870, 535)],
                        conditions = (masada_conditions..., CaM = "WT"), x = masada_x,
                        y = AxisMap(ticks([112.0, 255.5, 398.0, 541.5], 1.00:-0.02:0.94)),
                        axis = (; xticks = 0:0.2:1, yticks = 0.94:0.02:1.0)),
              PlotPanel(; name = "B", title = "", grid = (2, 1), frame = (143, 779, 855, 1255), margin = 5,
                        crop = (0, 690, 874, 1380), skip = [(705, 1162, 870, 1220)],
                        conditions = (masada_conditions..., CaM = "CaM12"), x = masada_x,
                        y = AxisMap(ticks([803.0, 947.0, 1090.0, 1233.5], 1.00:-0.02:0.94)),
                        axis = (; xticks = 0:0.2:1, yticks = 0.94:0.02:1.0)),
              PlotPanel(; name = "C_inset", title = "", grid = (3, 1), frame = (389, 1622, 556, 1872),
                        read = (393, 1624, 556, 1870),
                        crop = (310, 1570, 630, 1915), conditions = (masada_conditions..., CaM = "CaM34"),
                        x = AxisMap(ticks([389.0, 472.0, 554.5], 0.00:0.02:0.04)),
                        y = AxisMap(ticks([1623.0, 1748.0, 1872.0], 1.00:-0.01:0.98)),
                        axis = (; xticks = 0:0.02:0.04, yticks = 0.98:0.01:1.0))],
    series = [TraceSeries(; name, panel, ink = :grey, gap = 4, radius = 1, conditions = (target = "none",),
                          colour = "#000000")
              for (panel, name) in (("A", "CaM WT"), ("B", "CaM12"), ("C_inset", "CaM34"))],
    axis = (; xtickalign = 0, ytickalign = 0))

## Digitize the figures named in RUN
figures = RUN === :all ? [n for n in names(@__MODULE__; all = true)
                          if isdefined(@__MODULE__, n) && getfield(@__MODULE__, n) isa PaperFigure] : RUN
for name in figures
    println(name)
    digitize(getfield(@__MODULE__, name))
end
