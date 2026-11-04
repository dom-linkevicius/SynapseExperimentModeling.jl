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

## Williams et al. 2023, Cell Calcium, doi:10.1016/j.ceca.2023.102752
# Fig 10a,c: Ca2+ binding of the CaM C- and N-lobes, WT vs N53I CaM, without (a) and with
# (c) CaN peptide. Squares in a, circles in c; the black lines are fits, not data.
# Panels: frame (left, top, right, bottom), x tick pixels at 1E-9 ... 1E-2, y tick pixels at
# 1.0 ... 0.0, legend text box, and the crop shown in the calibration figure.
williams_panel(name, grid, frame, xticks, yticks, legend, crop) = PlotPanel(;
    name, grid, frame, crop, title = name[3:end], skip = [legend],
    conditions = (CaM_lobe = string(last(name)), CaN_peptide = startswith(name, "c")),
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

## Quintana et al. 2005, BBRC, doi:10.1016/j.bbrc.2005.06.152
# Small JPEGs of Origin plots: 1 px axes with outward ticks, black open markers. Dense
# kinetic traces are read as bands (TraceSeries), sparse points as markers.
ticks(px, values) = [p => v for (p, v) in zip(px, values)]
quintana_axis = (; xtickalign = 0, ytickalign = 0, xminortickalign = 0, yminortickalign = 0,
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
                        conditions = (Fluo4FF_uM = 10.0,), skip = [(1086, 920, 1535, 932), (1086, 1590, 1535, 1602)],
                        x = AxisMap(ticks([1184, 1280.5, 1378.5, 1475.5], 20:20:80)), y = shifman_fig3_y,
                        axis = (; xticks = 0:20:80)),
              PlotPanel(; name = "D_right", title = "", grid = (1, 2), frame = (1588, 920, 1763, 1602),
                        crop = (1580, 880, 1799, 1713), conditions = (Fluo4FF_uM = 10.0,),
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

## Digitize the figures named in RUN
figures = RUN === :all ? [n for n in names(@__MODULE__; all = true)
                          if isdefined(@__MODULE__, n) && getfield(@__MODULE__, n) isa PaperFigure] : RUN
for name in figures
    println(name)
    digitize(getfield(@__MODULE__, name))
end
