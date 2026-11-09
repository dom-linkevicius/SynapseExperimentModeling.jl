using Pkg
const ROOT = "/home/d/d-linkevicius/SynapseExperimentModeling.jl"
Pkg.activate(ROOT)
using Revise
using SynapseExperimentModeling, NoLimits, JLD2, Random
using LinearAlgebra: diag, eigvals, Symmetric
import Turing

include(joinpath(ROOT, "model_zoo", "calmodulin.jl"))

const SOLVE_FAILURES = Dict{String, Int}()
const SOLVE_FAILURES_LOCK = ReentrantLock()
function record_solve_failure(retcode)
    lock(SOLVE_FAILURES_LOCK) do
        SOLVE_FAILURES[string(retcode)] = get(SOLVE_FAILURES, string(retcode), 0) + 1
    end
    return nothing
end
@eval NoLimits _ll_drop_solve(retcode) = $record_solve_failure(retcode)

const METRICS = Dict("diag" => Turing.Inference.AHMC.DiagEuclideanMetric,
                     "dense" => Turing.Inference.AHMC.DenseEuclideanMetric)
const N_ADAPT = 1000
const N_SAMPLES = 1000

model_name, metric, chain_id = ARGS[1], ARGS[2], parse(Int, ARGS[3])
haskey(MODEL_EXPERIMENTS, Symbol(model_name)) || error("unknown model $(model_name), expected one of $(keys(MODEL_EXPERIMENTS))")
haskey(METRICS, metric) || error("unknown metric $(metric), expected one of $(keys(METRICS))")
dm = datamodel(model_name)
run_name = "$(model_name), $(metric) metric, chain $(chain_id)"

println(run_name, ": started ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
flush(stdout)
flusher = Timer(_ -> (flush(stdout); flush(stderr)), 60; interval = 60)

t0 = time()
last_metric = Ref{Any}(nothing)
last_failures = Ref(Dict{String, Int}())
function log_iteration(rng, model, sampler, sample, state, i; kwargs...)
    failures = lock(() -> copy(SOLVE_FAILURES), SOLVE_FAILURES_LOCK)
    new_failures = Dict(k => v - get(last_failures[], k, 0) for (k, v) in failures if v > get(last_failures[], k, 0))
    last_failures[] = failures
    println("iteration ", i, "/", N_ADAPT + N_SAMPLES, i <= N_ADAPT ? " (adapt)" : "", ", ",
            round(time() - t0; digits = 1), " s, tree depth ", get(sample.stats, :tree_depth, missing),
            ", steps ", get(sample.stats, :n_steps, missing), ", step size ", get(sample.stats, :step_size, missing),
            ", divergent ", get(sample.stats, :numerical_error, missing),
            ", acceptance ", round(get(sample.stats, :acceptance_rate, NaN); digits = 3),
            ", log density ", round(get(sample.stats, :log_density, NaN); digits = 1),
            ", failed ODE solves ", sum(values(new_failures); init = 0),
            isempty(new_failures) ? "" : " ($(join(("$k: $v" for (k, v) in new_failures), ", ")))")
    M⁻¹ = state.hamiltonian.metric.M⁻¹
    if M⁻¹ != last_metric[]
        last_metric[] = copy(M⁻¹)
        println("metric changed at iteration ", i, ", diag(M⁻¹) = ",
                round.(M⁻¹ isa AbstractMatrix ? diag(M⁻¹) : M⁻¹; sigdigits = 3))
        M⁻¹ isa AbstractMatrix && println("  eigvals(M⁻¹) = ", round.(eigvals(Symmetric(M⁻¹)); sigdigits = 3))
    end
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

println("failed ODE solves in total: ", isempty(SOLVE_FAILURES) ? 0 : SOLVE_FAILURES)
println(run_name, ": finished ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
close(flusher)
