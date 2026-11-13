# Import data that authors shared as tables (not digitized from figures) into the same
# layout as the digitized figures: data/exp_raw/<paper>.xlsx with a data tab (columns panel,
# series, x, y, merged) and META rows saying what they are, read back with
# load_paper_data(paper). One section per dataset; the original files stay untouched in
# papers/<paper>/.

using Pkg
Pkg.activate(dirname(@__DIR__))  # repo root; this file is in scripts/
using Revise
using SynapseExperimentModeling, DataFrames
import XLSX

## Faas et al. 2011, Nat Neurosci, doi:10.1038/nn.2746
# Flash-photolysis traces shared privately by the authors, in
# papers/Faas_etal_2011/Faas_etal_2011_original.xlsx. Sheet1 has one column per trace, named
# in row 1: the trace's parameters in the rows above "TIME", then F/F_t=0 at the times in
# column E (ms, not equidistant); the authors' notes are in column A. Sheet2 holds their model
# equations and stays in the original only. All experiments were at 35.0 ± 0.1 °C (Online
# Methods).
faas_parameters = [  # Sheet1 label in column D => condition column
    "fast DM fraction" => "DMn_fast_fraction",
    "tau fast" => "DMn_tau_fast_ms",
    "tau slow" => "DMn_tau_slow_ms",
    "DM total" => "DMn_total_M",
    "Kd_DM" => "DMn_Kd_M",
    "K_on DM Ca2+" => "DMn_kon_per_M_per_ms",
    "K_off DM (calculated)" => "DMn_koff_per_ms",
    "Kd_PP (calculated)" => "photoproduct_Kd_M",
    "K_on PP (same as DM)" => "photoproduct_kon_per_M_per_ms",
    "K_off PP" => "photoproduct_koff_per_ms",
    "D total" => "OGB5N_total_M",
    "Kd_D" => "OGB5N_Kd_M",
    "K_on D (calculated)" => "OGB5N_kon_per_M_per_ms",
    "K_off D" => "OGB5N_koff_per_ms",
    "Ratio D (F_max/F_min)" => "OGB5N_Fmax_over_Fmin",
    "Ca2+ free_t=0" => "Ca_free_t0_M",
    "B total" => "CaM_total_M",
    "pockes cell delay" => "pockels_cell_delay_us",
    "first approximation of uncaged fraction (calculated)" => "uncaged_fraction_approx",
]
let paper = "Faas_etal_2011"
    m = XLSX.readxlsx(joinpath(dirname(@__DIR__), "papers", paper, "Faas_etal_2011_original.xlsx"))["Sheet1"][:]
    header = findfirst(isequal("TIME"), m[:, 5])
    row(label) = only(findall(isequal(label), m[1:header, 4]))
    traces = [(j, m[1, j]) for j in 6:size(m, 2) if !ismissing(m[1, j])]
    batch(trace) = first(rsplit(trace, "_"; limit = 2))  # "#1021_WT_360" → "#1021_WT"

    data = DataFrame([(panel = batch(t), series = t, x = Float64(m[i, 5]), y = Float64(m[i, j]), merged = false)
                      for (j, t) in traces for i in (header + 1):size(m, 1) if !ismissing(m[i, j])])
    notes = join([string(v) for v in m[:, 1] if !ismissing(v)], " ")
    description = "Flash photolysis of DM-nitrophen (DMn) with WT CaM and the Ca2+ indicator OGB-5N: " *
                  "fluorescence relative to t = 0 (F/F_t=0) vs time, at non-equidistant time points. " *
                  "Original data shared privately by the authors, not digitized from the paper; the " *
                  "original file, with the authors' model equations, is " *
                  "papers/Faas_etal_2011/Faas_etal_2011_original.xlsx. Each trace's DMn, OGB-5N, Ca2+, " *
                  "CaM and laser parameters are its condition columns. Authors' notes: " * notes
    meta = DataFrame([(; figure = "TRACES", panel = batch(t), series = t, method = "shared privately",
                       x_quantity = "time", x_unit = "ms", y_quantity = "fluorescence_ratio", y_unit = "1",
                       CaM = "WT", (Symbol(name) => m[row(label), j] for (label, name) in faas_parameters)...,
                       temperature_C = 35.0,
                       description, doi = "10.1038/nn.2746")
                      for (j, t) in traces])
    path = write_paper_data(SynapseExperimentModeling.exp_raw_path(paper), "TRACES",
                            sort!(data, [:panel, :series, :x]), meta)
    println("$paper: $(length(traces)) traces, $(nrow(data)) points → ", relpath(path, dirname(@__DIR__)))
end
