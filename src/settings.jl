# Settings shared by the models in model_zoo/.

using NoLimits: AutoVern7, Rodas5P, AutoFiniteDiff

# ODE solver: Vern7, switching to Rodas5P where the problem turns stiff
const SOLVER = AutoVern7(Rodas5P(autodiff = AutoFiniteDiff()))
# solver tolerances, for concentrations in µM
const ABSTOL = 1e-9
const RELTOL = 1e-8
