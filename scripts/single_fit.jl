using Pkg
Pkg.activate(dirname(@__DIR__))
using Revise
using SynapseExperimentModeling, NoLimits, CairoMakie, AlgebraOfGraphics
import Turing
using TerminalLoggers, Logging

include(joinpath(dirname(@__DIR__), "model_zoo", "calmodulin.jl"))
dm = datamodel()

###res = fit_model(
###    dm,
###    NoLimits.FOCEI(;
###        interaction = true,
###        optim_kwargs = (;
###            maxiters = 1,
###            local_maxiters = 20,
###            show_trace = true,
###            extended_trace = false,
###        )
###    )
###)
###
###### continue fit
###res = fit_model(
###    dm, 
###    NoLimits.FOCEI(;
###        interaction = true,
###        optim_kwargs = (;
###            maxiters = 1,
###            local_maxiters = 20,
###            show_trace = true,
###            extended_trace = false,
###        )
###    ),
###    theta_0_untransformed = get_params(res; scale = :untransformed)
###)

res = with_logger(TerminalLogger()) do
    fit_model(
        dm, 
        NoLimits.MCMC(; 
            sampler = Turing.NUTS(0.75), 
            turing_kwargs = (
                n_samples = 1000, 
                n_adapt = 500, 
                progress = true
            )
        )
    )
end

plot_calmodulin(res)


###
###
###begin
###    chain = get_chain(res)
###    ps = names(chain, :parameters)
###    fig = Figure(size = (1600, 1000));
###    for (i, p) in enumerate(ps[1:17])
###        ax = Axis(fig[fldmod1(i, 4)...]; title = string(p))
###        lines!(ax, vec(chain[p]))
###    end
###    fig
###end


###println(get_params(res; scale = :untransformed))
###println("objective ", get_objective(res), ", converged ", get_converged(res))
###f = plot_fits(
###    res, 
###    ncols = 10,
###    kwargs_layout = (; size = (2500, 1300)),
###    #kwargs_subplot = (; yticklabelsvisible = false),
###) # ;individuals_idx = 1:9)
###foreach(delete!, [c for c in f.content if c isa Legend])
###
###df = DataFrame(res)
###
###pred = predict(res, faas_data(); re_mode = :ebe)
###dd = faas_data()
###
###df_joint = innerjoin(dd, pred; on = [:series => :id, :t => :time])
###
###spec = AlgebraOfGraphics.data(df_joint) * (
###    mapping(:t, :OGB5_F_F0, layout=:series) * visual(Scatter, color = :lightblue, label="Data") +
###    mapping(:t, :prediction, layout=:series) * visual(Lines, color = :orange, label="Prediction")
###)
###draw(
###    spec,
###    figure = (;
###        size = (2000, 1000)
###    ),
###    axis = (;
###        ylabel = "OGB5_F_F0",
###        xlabel = "Time",
###    ),
###    legend = (;
###        position = :top,
###    )
###)