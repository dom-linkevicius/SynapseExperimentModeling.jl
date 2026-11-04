using Pkg
Pkg.activate(dirname(@__DIR__))
using Revise
using SynapseExperimentModeling, NoLimits, JLD2, Random
import Turing

include(joinpath(dirname(@__DIR__), "model_zoo", "calmodulin.jl"))

const MODELS = ("calmodulin_faas", "calmodulin_faas_shifman_fig1", "calmodulin_full")

model_name, chain_id = ARGS[1], parse(Int, ARGS[2])
version = findfirst(==(model_name), MODELS)
version === nothing && error("unknown model $(model_name), expected one of $(MODELS)")
dm = datamodel(version)

res = fit_model(
    dm,
    NoLimits.MCMC(;
        sampler = Turing.NUTS(0.75),
        turing_kwargs = (
            n_samples = 1000,
            n_adapt = 500,
            progress = true
        )
    );
    rng = Xoshiro(chain_id)
)

dir = mkpath(joinpath(dirname(@__DIR__), "results", model_name))
jldsave(joinpath(dir, "chain_$(chain_id).jld2"); chain = get_chain(res))
save_fit(joinpath(dir, "fit_$(chain_id).jld2"), res; include_data = true)
