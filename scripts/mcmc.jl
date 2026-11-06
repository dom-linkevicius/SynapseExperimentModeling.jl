using Pkg
const ROOT = "/home/d/d-linkevicius/SynapseExperimentModeling.jl"
Pkg.activate(ROOT)
using Revise
using SynapseExperimentModeling, NoLimits, JLD2, Random
import Turing

include(joinpath(ROOT, "model_zoo", "calmodulin.jl"))

const MODELS = ("calmodulin_faas", "calmodulin_faas_shifman_fig1", "calmodulin_full")
const METRICS = Dict("diag" => Turing.Inference.AHMC.DiagEuclideanMetric,
                     "dense" => Turing.Inference.AHMC.DenseEuclideanMetric)
const N_ADAPT = 500
const N_SAMPLES = 1000

model_name, metric, chain_id = ARGS[1], ARGS[2], parse(Int, ARGS[3])
version = findfirst(==(model_name), MODELS)
version === nothing && error("unknown model $(model_name), expected one of $(MODELS)")
haskey(METRICS, metric) || error("unknown metric $(metric), expected one of $(keys(METRICS))")
dm = datamodel(version)
run_name = "$(model_name), $(metric) metric, chain $(chain_id)"

println(run_name, ": started ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
flush(stdout)
flusher = Timer(_ -> (flush(stdout); flush(stderr)), 60; interval = 60)

t0 = time()
function log_iteration(rng, model, sampler, sample, state, i; kwargs...)
    println("iteration ", i, "/", N_ADAPT + N_SAMPLES, i <= N_ADAPT ? " (adapt)" : "", ", ",
            round(time() - t0; digits = 1), " s, tree depth ", get(sample.stats, :tree_depth, missing),
            ", steps ", get(sample.stats, :n_steps, missing), ", step size ", get(sample.stats, :step_size, missing))
    flush(stdout)
end

res = fit_model(
    dm,
    NoLimits.MCMC(;
        sampler = Turing.NUTS(0.75; metricT = METRICS[metric]),
        turing_kwargs = (
            n_samples = N_ADAPT + N_SAMPLES,
            n_adapt = N_ADAPT,
            discard_initial = 0,
            callback = log_iteration,
            progress = true
        )
    );
    rng = Xoshiro(chain_id)
)

dir = mkpath(joinpath(ROOT, "results", model_name))
jldsave(joinpath(dir, "$(metric)_chain_$(chain_id).jld2"); chain = get_chain(res))
save_fit(joinpath(dir, "$(metric)_fit_$(chain_id).jld2"), res; include_data = true)

println(run_name, ": finished ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
close(flusher)
