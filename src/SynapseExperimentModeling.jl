module SynapseExperimentModeling

using DataFrames

include("settings.jl")
include("digitize.jl")
include("kinetics.jl")

export SOLVER, ABSTOL, RELTOL
export AxisMap, PlotPanel, MarkerSeries, TraceSeries, PaperFigure, digitize, load_paper_data, write_paper_data
export equilibrium

end
