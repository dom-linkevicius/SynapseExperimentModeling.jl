### DL: Byrne et al (2008) scheme with unique CaM binding sites within lobes.
### DL: Had to re-write from Linkevicius et al. (2025) who had independent lobes,
### DL: since downstream models will need tracking of CaM molecules and 
### DL: model written with independent lobes cannot track that

using NoLimits, Distributions, DataFrames
using LinearAlgebra: kron
using BracketingNonlinearSolve: IntervalNonlinearProblem, ITP, solve
using SynapseExperimentModeling: load_paper_data, SOLVER, ABSTOL, RELTOL
using CairoMakie: Figure, Label, Scatter, Lines, colsize!, Relative
using AlgebraOfGraphics: AlgebraOfGraphics, mapping, visual, draw!, sorter

lobe_weights(k01, k10, k02, k20, k13, k31, Ca) = [one(Ca), Ca * k01 / k10, Ca * k02 / k20, Ca^2 * k01 / k10 * k13 / k31]

lobe_bound(w) = (w[2] + w[3] + 2 * w[4]) / sum(w)

function calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total)
    w_N = lobe_weights(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, Ca)
    w_C = lobe_weights(k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, Ca)
    return total * kron(w_C, w_N) / (sum(w_N) * sum(w_C))
end

const CALMODULIN_CALCIUM = [n_N + n_C for n_C in (0, 1, 1, 2) for n_N in (0, 1, 1, 2)]

bound_calcium(CaM) = sum(CALMODULIN_CALCIUM .* CaM)

calcium_excess(Ca, p) = Ca + p[15] * Ca / (Ca + p[16]) +
                        p[14] * (lobe_bound(lobe_weights(p[1:6]..., Ca)) + lobe_bound(lobe_weights(p[7:12]..., Ca))) - p[13]

function free_calcium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                      k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                      Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
    p = [k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k01_C, k10_C, k02_C, k20_C, k13_C, k31_C,
         Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF]
    return solve(IntervalNonlinearProblem(calcium_excess, (zero(Ca_total), Ca_total), p), ITP()).u
end

const EXPERIMENTS = ("Faas_etal_2011", "Shifman_etal_2006_Fig1B", "Shifman_etal_2006_Fig3D")

experiment_id(name) = findfirst(==(name), EXPERIMENTS)

function initial_calcium(experiment, k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                         k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                         Ca_free, Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
    EXPERIMENTS[experiment] == "Shifman_etal_2006_Fig3D" || return Ca_free
    return free_calcium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                        k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                        Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
end

calmodulin_faas = @Model begin
    @helpers begin
        calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                               k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total) =
            calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                   k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total)
    end

    @fixedEffects begin
        tv_k01_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_N = RealNumber(log10(0.5072), prior = Normal(log10(0.5072), 3))    # log10(µM^-1 ms^-1)
        tv_k23_N = RealNumber(log10(0.5), prior = Normal(log10(0.5), 3))          # log10(µM^-1 ms^-1)
        tv_K01d_N = RealNumber(log10(33.0), prior = Normal(log10(33.0), 3))       # log10(µM)
        tv_K13d_N = RealNumber(log10(230.0), prior = Normal(log10(230.0), 3))     # log10(µM)
        tv_K02d_N = RealNumber(log10(3.45), prior = Normal(log10(3.45), 3))       # log10(µM)
        tv_k01_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_C = RealNumber(log10(0.00371), prior = Normal(log10(0.00371), 3))  # log10(µM^-1 ms^-1)
        tv_k23_C = RealNumber(log10(0.118), prior = Normal(log10(0.118), 3))      # log10(µM^-1 ms^-1)
        tv_K01d_C = RealNumber(log10(18.5), prior = Normal(log10(18.5), 3))       # log10(µM)
        tv_K02d_C = RealNumber(log10(116.0), prior = Normal(log10(116.0), 3))     # log10(µM)
        tv_K13d_C = RealNumber(log10(0.38), prior = Normal(log10(0.38), 3))       # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0,
                        prior = truncated(Normal(-1, 1); lower = -5, upper = 5))                 # 1
        Ω = RealNumber(1.5, scale = :log,
                        prior = truncated(Normal(0, 1); lower = 0)) # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0,
                                 prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
    end

    @covariates begin
        t = Covariate()                                         # ms
        f_frac = ConstantCovariate(; constant_on = :series)     # 1
        τ_f = ConstantCovariate(; constant_on = :series)        # ms
        τ_s = ConstantCovariate(; constant_on = :series)        # ms
        kon_DMn = ConstantCovariate(; constant_on = :series)    # µM^-1 ms^-1
        koff_DMn = ConstantCovariate(; constant_on = :series)   # ms^-1
        Kd_DMn = ConstantCovariate(; constant_on = :series)     # µM
        koff_PP = ConstantCovariate(; constant_on = :series)    # ms^-1
        DMn_t = ConstantCovariate(; constant_on = :series)      # µM
        CaM_t = ConstantCovariate(; constant_on = :series)      # µM
        OGB5_t = ConstantCovariate(; constant_on = :series)     # µM
        Ca_free = ConstantCovariate(; constant_on = :series)    # µM
        kon_D = ConstantCovariate(; constant_on = :series)      # µM^-1 ms^-1
        koff_D = ConstantCovariate(; constant_on = :series)     # ms^-1
        Fmax_Fmin = ConstantCovariate(; constant_on = :series)  # 1
    end

    @randomEffects begin
        U = RandomEffect(LogitNormal(μ, Ω); column = :U_group)  # 1
    end

    @preDifferentialEquation begin
        kon_PP = kon_DMn
        Kd_D = koff_D / kon_D

        K01d_N = 10^tv_K01d_N
        k01_N = 10^tv_k01_N
        k10_N = K01d_N * 10^tv_k01_N

        K02d_N = 10^tv_K02d_N
        k02_N = 10^tv_k02_N
        k20_N = K02d_N * 10^tv_k02_N

        K13d_N = 10^tv_K13d_N
        k13_N = 10^tv_k13_N
        k31_N = K13d_N * 10^tv_k13_N

        K23d_N = K01d_N * K13d_N / K02d_N
        k23_N = 10^tv_k23_N
        k32_N = K23d_N * 10^tv_k23_N

        K01d_C = 10^tv_K01d_C
        k01_C = 10^tv_k01_C
        k10_C = K01d_C * 10^tv_k01_C

        K02d_C = 10^tv_K02d_C
        k02_C = 10^tv_k02_C
        k20_C = K02d_C * 10^tv_k02_C

        K13d_C = 10^tv_K13d_C
        k13_C = 10^tv_k13_C
        k31_C = K13d_C * 10^tv_k13_C

        K23d_C = K01d_C * K13d_C / K02d_C
        k23_C = 10^tv_k23_C
        k32_C = K23d_C * 10^tv_k23_C

        DMn₀ = (Kd_DMn * DMn_t) / (Ca_free + Kd_DMn)
        CaDMn₀ = DMn_t - DMn₀
        OGB5₀ = (Kd_D * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀ = OGB5_t - OGB5₀
        CaM₀ = calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                      k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                                      Ca_free, CaM_t)

        DMn_s₀_u = (1 - f_frac) * DMn₀ * U
        CaDMn_s₀_u = (1 - f_frac) * CaDMn₀ * U
        DMn_f₀_u = f_frac * DMn₀ * U
        CaDMn_f₀_u = f_frac * CaDMn₀ * U
    end

    @DifferentialEquation begin
        JN01_C0(t) = k01_N * Ca * CaM_N0C0 - k10_N * CaM_N1C0
        JN02_C0(t) = k02_N * Ca * CaM_N0C0 - k20_N * CaM_N2C0
        JN13_C0(t) = k13_N * Ca * CaM_N1C0 - k31_N * CaM_N3C0
        JN23_C0(t) = k23_N * Ca * CaM_N2C0 - k32_N * CaM_N3C0
        JN01_C1(t) = k01_N * Ca * CaM_N0C1 - k10_N * CaM_N1C1
        JN02_C1(t) = k02_N * Ca * CaM_N0C1 - k20_N * CaM_N2C1
        JN13_C1(t) = k13_N * Ca * CaM_N1C1 - k31_N * CaM_N3C1
        JN23_C1(t) = k23_N * Ca * CaM_N2C1 - k32_N * CaM_N3C1
        JN01_C2(t) = k01_N * Ca * CaM_N0C2 - k10_N * CaM_N1C2
        JN02_C2(t) = k02_N * Ca * CaM_N0C2 - k20_N * CaM_N2C2
        JN13_C2(t) = k13_N * Ca * CaM_N1C2 - k31_N * CaM_N3C2
        JN23_C2(t) = k23_N * Ca * CaM_N2C2 - k32_N * CaM_N3C2
        JN01_C3(t) = k01_N * Ca * CaM_N0C3 - k10_N * CaM_N1C3
        JN02_C3(t) = k02_N * Ca * CaM_N0C3 - k20_N * CaM_N2C3
        JN13_C3(t) = k13_N * Ca * CaM_N1C3 - k31_N * CaM_N3C3
        JN23_C3(t) = k23_N * Ca * CaM_N2C3 - k32_N * CaM_N3C3
        JC01_N0(t) = k01_C * Ca * CaM_N0C0 - k10_C * CaM_N0C1
        JC02_N0(t) = k02_C * Ca * CaM_N0C0 - k20_C * CaM_N0C2
        JC13_N0(t) = k13_C * Ca * CaM_N0C1 - k31_C * CaM_N0C3
        JC23_N0(t) = k23_C * Ca * CaM_N0C2 - k32_C * CaM_N0C3
        JC01_N1(t) = k01_C * Ca * CaM_N1C0 - k10_C * CaM_N1C1
        JC02_N1(t) = k02_C * Ca * CaM_N1C0 - k20_C * CaM_N1C2
        JC13_N1(t) = k13_C * Ca * CaM_N1C1 - k31_C * CaM_N1C3
        JC23_N1(t) = k23_C * Ca * CaM_N1C2 - k32_C * CaM_N1C3
        JC01_N2(t) = k01_C * Ca * CaM_N2C0 - k10_C * CaM_N2C1
        JC02_N2(t) = k02_C * Ca * CaM_N2C0 - k20_C * CaM_N2C2
        JC13_N2(t) = k13_C * Ca * CaM_N2C1 - k31_C * CaM_N2C3
        JC23_N2(t) = k23_C * Ca * CaM_N2C2 - k32_C * CaM_N2C3
        JC01_N3(t) = k01_C * Ca * CaM_N3C0 - k10_C * CaM_N3C1
        JC02_N3(t) = k02_C * Ca * CaM_N3C0 - k20_C * CaM_N3C2
        JC13_N3(t) = k13_C * Ca * CaM_N3C1 - k31_C * CaM_N3C3
        JC23_N3(t) = k23_C * Ca * CaM_N3C2 - k32_C * CaM_N3C3

        D(Ca) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                -kon_PP * PP * Ca + koff_PP * CaPP +
                -kon_D * OGB5 * Ca + koff_D * CaOGB5 +
                -JN01_C0(t) - JN02_C0(t) - JN13_C0(t) - JN23_C0(t) +
                -JN01_C1(t) - JN02_C1(t) - JN13_C1(t) - JN23_C1(t) +
                -JN01_C2(t) - JN02_C2(t) - JN13_C2(t) - JN23_C2(t) +
                -JN01_C3(t) - JN02_C3(t) - JN13_C3(t) - JN23_C3(t) +
                -JC01_N0(t) - JC02_N0(t) - JC13_N0(t) - JC23_N0(t) +
                -JC01_N1(t) - JC02_N1(t) - JC13_N1(t) - JC23_N1(t) +
                -JC01_N2(t) - JC02_N2(t) - JC13_N2(t) - JC23_N2(t) +
                -JC01_N3(t) - JC02_N3(t) - JC13_N3(t) - JC23_N3(t)
        D(DMn_s) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s - DMn_s / τ_s
        D(CaDMn_s) ~ kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s / τ_s
        D(DMn_f) ~ -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f - DMn_f / τ_f
        D(CaDMn_f) ~ kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f / τ_f
        D(PP) ~ -kon_PP * PP * Ca + koff_PP * CaPP +
                2 * DMn_s / τ_s + 2 * DMn_f / τ_f + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(CaPP) ~ kon_PP * PP * Ca - koff_PP * CaPP + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(OGB5) ~ -kon_D * OGB5 * Ca + koff_D * CaOGB5
        D(CaOGB5) ~ kon_D * OGB5 * Ca - koff_D * CaOGB5
        D(CaM_N0C0) ~ -JN01_C0(t) - JN02_C0(t) - JC01_N0(t) - JC02_N0(t)
        D(CaM_N1C0) ~ JN01_C0(t) - JN13_C0(t) - JC01_N1(t) - JC02_N1(t)
        D(CaM_N2C0) ~ JN02_C0(t) - JN23_C0(t) - JC01_N2(t) - JC02_N2(t)
        D(CaM_N3C0) ~ JN13_C0(t) + JN23_C0(t) - JC01_N3(t) - JC02_N3(t)
        D(CaM_N0C1) ~ -JN01_C1(t) - JN02_C1(t) + JC01_N0(t) - JC13_N0(t)
        D(CaM_N1C1) ~ JN01_C1(t) - JN13_C1(t) + JC01_N1(t) - JC13_N1(t)
        D(CaM_N2C1) ~ JN02_C1(t) - JN23_C1(t) + JC01_N2(t) - JC13_N2(t)
        D(CaM_N3C1) ~ JN13_C1(t) + JN23_C1(t) + JC01_N3(t) - JC13_N3(t)
        D(CaM_N0C2) ~ -JN01_C2(t) - JN02_C2(t) + JC02_N0(t) - JC23_N0(t)
        D(CaM_N1C2) ~ JN01_C2(t) - JN13_C2(t) + JC02_N1(t) - JC23_N1(t)
        D(CaM_N2C2) ~ JN02_C2(t) - JN23_C2(t) + JC02_N2(t) - JC23_N2(t)
        D(CaM_N3C2) ~ JN13_C2(t) + JN23_C2(t) + JC02_N3(t) - JC23_N3(t)
        D(CaM_N0C3) ~ -JN01_C3(t) - JN02_C3(t) + JC13_N0(t) + JC23_N0(t)
        D(CaM_N1C3) ~ JN01_C3(t) - JN13_C3(t) + JC13_N1(t) + JC23_N1(t)
        D(CaM_N2C3) ~ JN02_C3(t) - JN23_C3(t) + JC13_N2(t) + JC23_N2(t)
        D(CaM_N3C3) ~ JN13_C3(t) + JN23_C3(t) + JC13_N3(t) + JC23_N3(t)
    end

    @initialDE begin
        DMn_s = DMn_s₀_u
        CaDMn_s = CaDMn_s₀_u
        DMn_f = DMn_f₀_u
        CaDMn_f = CaDMn_f₀_u
        PP = 0.0
        CaPP = 0.0
        OGB5 = OGB5₀
        CaOGB5 = CaOGB5₀
        CaM_N0C0 = CaM₀[1]
        CaM_N1C0 = CaM₀[2]
        CaM_N2C0 = CaM₀[3]
        CaM_N3C0 = CaM₀[4]
        CaM_N0C1 = CaM₀[5]
        CaM_N1C1 = CaM₀[6]
        CaM_N2C1 = CaM₀[7]
        CaM_N3C1 = CaM₀[8]
        CaM_N0C2 = CaM₀[9]
        CaM_N1C2 = CaM₀[10]
        CaM_N2C2 = CaM₀[11]
        CaM_N3C2 = CaM₀[12]
        CaM_N0C3 = CaM₀[13]
        CaM_N1C3 = CaM₀[14]
        CaM_N2C3 = CaM₀[15]
        CaM_N3C3 = CaM₀[16]
        Ca = Ca_free
    end

    @formulas begin
        OGB5_F_F0 ~ Normal((OGB5(t) + Fmax_Fmin * CaOGB5(t)) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_OGB5_F_F0)
    end
end

calmodulin_faas = set_solver_config(calmodulin_faas; alg = SOLVER, kwargs = (abstol = ABSTOL, reltol = RELTOL))

calmodulin_faas_shifman_fig1 = @Model begin
    @helpers begin
        calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                               k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total) =
            calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                   k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total)
        bound_calcium(CaM) = bound_calcium(CaM)
    end

    @fixedEffects begin
        tv_k01_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_N = RealNumber(log10(0.5072), prior = Normal(log10(0.5072), 3))    # log10(µM^-1 ms^-1)
        tv_k23_N = RealNumber(log10(0.5), prior = Normal(log10(0.5), 3))          # log10(µM^-1 ms^-1)
        tv_K01d_N = RealNumber(log10(33.0), prior = Normal(log10(33.0), 3))       # log10(µM)
        tv_K13d_N = RealNumber(log10(230.0), prior = Normal(log10(230.0), 3))     # log10(µM)
        tv_K02d_N = RealNumber(log10(3.45), prior = Normal(log10(3.45), 3))       # log10(µM)
        tv_k01_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_C = RealNumber(log10(0.00371), prior = Normal(log10(0.00371), 3))  # log10(µM^-1 ms^-1)
        tv_k23_C = RealNumber(log10(0.118), prior = Normal(log10(0.118), 3))      # log10(µM^-1 ms^-1)
        tv_K01d_C = RealNumber(log10(18.5), prior = Normal(log10(18.5), 3))       # log10(µM)
        tv_K02d_C = RealNumber(log10(116.0), prior = Normal(log10(116.0), 3))     # log10(µM)
        tv_K13d_C = RealNumber(log10(0.38), prior = Normal(log10(0.38), 3))       # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0,
                        prior = truncated(Normal(-1, 1); lower = -5, upper = 5))                 # 1
        Ω = RealNumber(1.5, scale = :log,
                        prior = truncated(Normal(0, 1); lower = 0)) # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0,
                                 prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
        σ_bound_Ca_per_CaM = RealNumber(0.3, scale = :log, lower = 1e-10, upper = 1.0,
                                        prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
    end

    @covariates begin
        t = Covariate()                                         # ms
        f_frac = ConstantCovariate(; constant_on = :series)     # 1
        τ_f = ConstantCovariate(; constant_on = :series)        # ms
        τ_s = ConstantCovariate(; constant_on = :series)        # ms
        kon_DMn = ConstantCovariate(; constant_on = :series)    # µM^-1 ms^-1
        koff_DMn = ConstantCovariate(; constant_on = :series)   # ms^-1
        Kd_DMn = ConstantCovariate(; constant_on = :series)     # µM
        koff_PP = ConstantCovariate(; constant_on = :series)    # ms^-1
        DMn_t = ConstantCovariate(; constant_on = :series)      # µM
        CaM_t = ConstantCovariate(; constant_on = :series)      # µM
        OGB5_t = ConstantCovariate(; constant_on = :series)     # µM
        Ca_free = ConstantCovariate(; constant_on = :series)    # µM
        kon_D = ConstantCovariate(; constant_on = :series)      # µM^-1 ms^-1
        koff_D = ConstantCovariate(; constant_on = :series)     # ms^-1
        Fmax_Fmin = ConstantCovariate(; constant_on = :series)  # 1
    end

    @randomEffects begin
        U = RandomEffect(LogitNormal(μ, Ω); column = :U_group)  # 1
    end

    @preDifferentialEquation begin
        kon_PP = kon_DMn
        Kd_D = koff_D / kon_D

        K01d_N = 10^tv_K01d_N
        k01_N = 10^tv_k01_N
        k10_N = K01d_N * 10^tv_k01_N

        K02d_N = 10^tv_K02d_N
        k02_N = 10^tv_k02_N
        k20_N = K02d_N * 10^tv_k02_N

        K13d_N = 10^tv_K13d_N
        k13_N = 10^tv_k13_N
        k31_N = K13d_N * 10^tv_k13_N

        K23d_N = K01d_N * K13d_N / K02d_N
        k23_N = 10^tv_k23_N
        k32_N = K23d_N * 10^tv_k23_N

        K01d_C = 10^tv_K01d_C
        k01_C = 10^tv_k01_C
        k10_C = K01d_C * 10^tv_k01_C

        K02d_C = 10^tv_K02d_C
        k02_C = 10^tv_k02_C
        k20_C = K02d_C * 10^tv_k02_C

        K13d_C = 10^tv_K13d_C
        k13_C = 10^tv_k13_C
        k31_C = K13d_C * 10^tv_k13_C

        K23d_C = K01d_C * K13d_C / K02d_C
        k23_C = 10^tv_k23_C
        k32_C = K23d_C * 10^tv_k23_C

        DMn₀ = (Kd_DMn * DMn_t) / (Ca_free + Kd_DMn)
        CaDMn₀ = DMn_t - DMn₀
        OGB5₀ = (Kd_D * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀ = OGB5_t - OGB5₀
        CaM₀ = calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                      k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                                      Ca_free, CaM_t)

        DMn_s₀_u = (1 - f_frac) * DMn₀ * U
        CaDMn_s₀_u = (1 - f_frac) * CaDMn₀ * U
        DMn_f₀_u = f_frac * DMn₀ * U
        CaDMn_f₀_u = f_frac * CaDMn₀ * U
    end

    @DifferentialEquation begin
        JN01_C0(t) = k01_N * Ca * CaM_N0C0 - k10_N * CaM_N1C0
        JN02_C0(t) = k02_N * Ca * CaM_N0C0 - k20_N * CaM_N2C0
        JN13_C0(t) = k13_N * Ca * CaM_N1C0 - k31_N * CaM_N3C0
        JN23_C0(t) = k23_N * Ca * CaM_N2C0 - k32_N * CaM_N3C0
        JN01_C1(t) = k01_N * Ca * CaM_N0C1 - k10_N * CaM_N1C1
        JN02_C1(t) = k02_N * Ca * CaM_N0C1 - k20_N * CaM_N2C1
        JN13_C1(t) = k13_N * Ca * CaM_N1C1 - k31_N * CaM_N3C1
        JN23_C1(t) = k23_N * Ca * CaM_N2C1 - k32_N * CaM_N3C1
        JN01_C2(t) = k01_N * Ca * CaM_N0C2 - k10_N * CaM_N1C2
        JN02_C2(t) = k02_N * Ca * CaM_N0C2 - k20_N * CaM_N2C2
        JN13_C2(t) = k13_N * Ca * CaM_N1C2 - k31_N * CaM_N3C2
        JN23_C2(t) = k23_N * Ca * CaM_N2C2 - k32_N * CaM_N3C2
        JN01_C3(t) = k01_N * Ca * CaM_N0C3 - k10_N * CaM_N1C3
        JN02_C3(t) = k02_N * Ca * CaM_N0C3 - k20_N * CaM_N2C3
        JN13_C3(t) = k13_N * Ca * CaM_N1C3 - k31_N * CaM_N3C3
        JN23_C3(t) = k23_N * Ca * CaM_N2C3 - k32_N * CaM_N3C3
        JC01_N0(t) = k01_C * Ca * CaM_N0C0 - k10_C * CaM_N0C1
        JC02_N0(t) = k02_C * Ca * CaM_N0C0 - k20_C * CaM_N0C2
        JC13_N0(t) = k13_C * Ca * CaM_N0C1 - k31_C * CaM_N0C3
        JC23_N0(t) = k23_C * Ca * CaM_N0C2 - k32_C * CaM_N0C3
        JC01_N1(t) = k01_C * Ca * CaM_N1C0 - k10_C * CaM_N1C1
        JC02_N1(t) = k02_C * Ca * CaM_N1C0 - k20_C * CaM_N1C2
        JC13_N1(t) = k13_C * Ca * CaM_N1C1 - k31_C * CaM_N1C3
        JC23_N1(t) = k23_C * Ca * CaM_N1C2 - k32_C * CaM_N1C3
        JC01_N2(t) = k01_C * Ca * CaM_N2C0 - k10_C * CaM_N2C1
        JC02_N2(t) = k02_C * Ca * CaM_N2C0 - k20_C * CaM_N2C2
        JC13_N2(t) = k13_C * Ca * CaM_N2C1 - k31_C * CaM_N2C3
        JC23_N2(t) = k23_C * Ca * CaM_N2C2 - k32_C * CaM_N2C3
        JC01_N3(t) = k01_C * Ca * CaM_N3C0 - k10_C * CaM_N3C1
        JC02_N3(t) = k02_C * Ca * CaM_N3C0 - k20_C * CaM_N3C2
        JC13_N3(t) = k13_C * Ca * CaM_N3C1 - k31_C * CaM_N3C3
        JC23_N3(t) = k23_C * Ca * CaM_N3C2 - k32_C * CaM_N3C3

        D(Ca) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                -kon_PP * PP * Ca + koff_PP * CaPP +
                -kon_D * OGB5 * Ca + koff_D * CaOGB5 +
                -JN01_C0(t) - JN02_C0(t) - JN13_C0(t) - JN23_C0(t) +
                -JN01_C1(t) - JN02_C1(t) - JN13_C1(t) - JN23_C1(t) +
                -JN01_C2(t) - JN02_C2(t) - JN13_C2(t) - JN23_C2(t) +
                -JN01_C3(t) - JN02_C3(t) - JN13_C3(t) - JN23_C3(t) +
                -JC01_N0(t) - JC02_N0(t) - JC13_N0(t) - JC23_N0(t) +
                -JC01_N1(t) - JC02_N1(t) - JC13_N1(t) - JC23_N1(t) +
                -JC01_N2(t) - JC02_N2(t) - JC13_N2(t) - JC23_N2(t) +
                -JC01_N3(t) - JC02_N3(t) - JC13_N3(t) - JC23_N3(t)
        D(DMn_s) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s - DMn_s / τ_s
        D(CaDMn_s) ~ kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s / τ_s
        D(DMn_f) ~ -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f - DMn_f / τ_f
        D(CaDMn_f) ~ kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f / τ_f
        D(PP) ~ -kon_PP * PP * Ca + koff_PP * CaPP +
                2 * DMn_s / τ_s + 2 * DMn_f / τ_f + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(CaPP) ~ kon_PP * PP * Ca - koff_PP * CaPP + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(OGB5) ~ -kon_D * OGB5 * Ca + koff_D * CaOGB5
        D(CaOGB5) ~ kon_D * OGB5 * Ca - koff_D * CaOGB5
        D(CaM_N0C0) ~ -JN01_C0(t) - JN02_C0(t) - JC01_N0(t) - JC02_N0(t)
        D(CaM_N1C0) ~ JN01_C0(t) - JN13_C0(t) - JC01_N1(t) - JC02_N1(t)
        D(CaM_N2C0) ~ JN02_C0(t) - JN23_C0(t) - JC01_N2(t) - JC02_N2(t)
        D(CaM_N3C0) ~ JN13_C0(t) + JN23_C0(t) - JC01_N3(t) - JC02_N3(t)
        D(CaM_N0C1) ~ -JN01_C1(t) - JN02_C1(t) + JC01_N0(t) - JC13_N0(t)
        D(CaM_N1C1) ~ JN01_C1(t) - JN13_C1(t) + JC01_N1(t) - JC13_N1(t)
        D(CaM_N2C1) ~ JN02_C1(t) - JN23_C1(t) + JC01_N2(t) - JC13_N2(t)
        D(CaM_N3C1) ~ JN13_C1(t) + JN23_C1(t) + JC01_N3(t) - JC13_N3(t)
        D(CaM_N0C2) ~ -JN01_C2(t) - JN02_C2(t) + JC02_N0(t) - JC23_N0(t)
        D(CaM_N1C2) ~ JN01_C2(t) - JN13_C2(t) + JC02_N1(t) - JC23_N1(t)
        D(CaM_N2C2) ~ JN02_C2(t) - JN23_C2(t) + JC02_N2(t) - JC23_N2(t)
        D(CaM_N3C2) ~ JN13_C2(t) + JN23_C2(t) + JC02_N3(t) - JC23_N3(t)
        D(CaM_N0C3) ~ -JN01_C3(t) - JN02_C3(t) + JC13_N0(t) + JC23_N0(t)
        D(CaM_N1C3) ~ JN01_C3(t) - JN13_C3(t) + JC13_N1(t) + JC23_N1(t)
        D(CaM_N2C3) ~ JN02_C3(t) - JN23_C3(t) + JC13_N2(t) + JC23_N2(t)
        D(CaM_N3C3) ~ JN13_C3(t) + JN23_C3(t) + JC13_N3(t) + JC23_N3(t)
    end

    @initialDE begin
        DMn_s = DMn_s₀_u
        CaDMn_s = CaDMn_s₀_u
        DMn_f = DMn_f₀_u
        CaDMn_f = CaDMn_f₀_u
        PP = 0.0
        CaPP = 0.0
        OGB5 = OGB5₀
        CaOGB5 = CaOGB5₀
        CaM_N0C0 = CaM₀[1]
        CaM_N1C0 = CaM₀[2]
        CaM_N2C0 = CaM₀[3]
        CaM_N3C0 = CaM₀[4]
        CaM_N0C1 = CaM₀[5]
        CaM_N1C1 = CaM₀[6]
        CaM_N2C1 = CaM₀[7]
        CaM_N3C1 = CaM₀[8]
        CaM_N0C2 = CaM₀[9]
        CaM_N1C2 = CaM₀[10]
        CaM_N2C2 = CaM₀[11]
        CaM_N3C2 = CaM₀[12]
        CaM_N0C3 = CaM₀[13]
        CaM_N1C3 = CaM₀[14]
        CaM_N2C3 = CaM₀[15]
        CaM_N3C3 = CaM₀[16]
        Ca = Ca_free
    end

    @formulas begin
        OGB5_F_F0 ~ Normal((OGB5(t) + Fmax_Fmin * CaOGB5(t)) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_OGB5_F_F0)
        bound_Ca_per_CaM ~ Normal(bound_calcium([CaM_N0C0(t), CaM_N1C0(t), CaM_N2C0(t), CaM_N3C0(t),
                                                 CaM_N0C1(t), CaM_N1C1(t), CaM_N2C1(t), CaM_N3C1(t),
                                                 CaM_N0C2(t), CaM_N1C2(t), CaM_N2C2(t), CaM_N3C2(t),
                                                 CaM_N0C3(t), CaM_N1C3(t), CaM_N2C3(t), CaM_N3C3(t)]) / CaM_t,
                                  σ_bound_Ca_per_CaM)
    end
end

calmodulin_faas_shifman_fig1 = set_solver_config(calmodulin_faas_shifman_fig1; alg = SOLVER, kwargs = (abstol = ABSTOL, reltol = RELTOL))

calmodulin_full = @Model begin
    @helpers begin
        calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                               k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total) =
            calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                   k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total)
        initial_calcium(experiment, k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                        k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                        Ca_free, Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF) =
            initial_calcium(experiment, k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                            k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                            Ca_free, Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
        bound_calcium(CaM) = bound_calcium(CaM)
    end

    @fixedEffects begin
        tv_k01_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_N = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_N = RealNumber(log10(0.5072), prior = Normal(log10(0.5072), 3))    # log10(µM^-1 ms^-1)
        tv_k23_N = RealNumber(log10(0.5), prior = Normal(log10(0.5), 3))          # log10(µM^-1 ms^-1)
        tv_K01d_N = RealNumber(log10(33.0), prior = Normal(log10(33.0), 3))       # log10(µM)
        tv_K13d_N = RealNumber(log10(230.0), prior = Normal(log10(230.0), 3))     # log10(µM)
        tv_K02d_N = RealNumber(log10(3.45), prior = Normal(log10(3.45), 3))       # log10(µM)
        tv_k01_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k02_C = RealNumber(log10(0.275), prior = Normal(log10(0.275), 3))      # log10(µM^-1 ms^-1)
        tv_k13_C = RealNumber(log10(0.00371), prior = Normal(log10(0.00371), 3))  # log10(µM^-1 ms^-1)
        tv_k23_C = RealNumber(log10(0.118), prior = Normal(log10(0.118), 3))      # log10(µM^-1 ms^-1)
        tv_K01d_C = RealNumber(log10(18.5), prior = Normal(log10(18.5), 3))       # log10(µM)
        tv_K02d_C = RealNumber(log10(116.0), prior = Normal(log10(116.0), 3))     # log10(µM)
        tv_K13d_C = RealNumber(log10(0.38), prior = Normal(log10(0.38), 3))       # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0,
                        prior = truncated(Normal(-1, 1); lower = -5, upper = 5))                 # 1
        Ω = RealNumber(1.5, scale = :log,
                        prior = truncated(Normal(0, 1); lower = 0)) # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0,
                                 prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
        σ_bound_Ca_per_CaM = RealNumber(0.3, scale = :log, lower = 1e-10, upper = 1.0,
                                        prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
        σ_Fluo4FF_fraction_bound = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0,
                                              prior = truncated(Normal(0, 1); lower = 1e-10, upper = 1))  # 1
    end

    @covariates begin
        t = Covariate()                                         # ms
        f_frac = ConstantCovariate(; constant_on = :series)     # 1
        τ_f = ConstantCovariate(; constant_on = :series)        # ms
        τ_s = ConstantCovariate(; constant_on = :series)        # ms
        kon_DMn = ConstantCovariate(; constant_on = :series)    # µM^-1 ms^-1
        koff_DMn = ConstantCovariate(; constant_on = :series)   # ms^-1
        Kd_DMn = ConstantCovariate(; constant_on = :series)     # µM
        koff_PP = ConstantCovariate(; constant_on = :series)    # ms^-1
        DMn_t = ConstantCovariate(; constant_on = :series)      # µM
        CaM_t = ConstantCovariate(; constant_on = :series)      # µM
        OGB5_t = ConstantCovariate(; constant_on = :series)     # µM
        Ca_free = ConstantCovariate(; constant_on = :series)    # µM
        kon_D = ConstantCovariate(; constant_on = :series)      # µM^-1 ms^-1
        koff_D = ConstantCovariate(; constant_on = :series)     # ms^-1
        Fmax_Fmin = ConstantCovariate(; constant_on = :series)  # 1
        Fluo4FF_t = ConstantCovariate(; constant_on = :series)  # µM
        Kd_Fluo4FF = ConstantCovariate(; constant_on = :series) # µM
        kon_Fluo4FF = ConstantCovariate(; constant_on = :series)    # µM^-1 ms^-1
        Ca_total = ConstantCovariate(; constant_on = :series)   # µM
        N_lobe = ConstantCovariate(; constant_on = :series)     # 1
        C_lobe = ConstantCovariate(; constant_on = :series)     # 1
        experiment = ConstantCovariate(; constant_on = :series) # 1
    end

    @randomEffects begin
        U = RandomEffect(LogitNormal(μ, Ω); column = :U_group)  # 1
    end

    @preDifferentialEquation begin
        kon_PP = kon_DMn
        Kd_D = koff_D / kon_D
        koff_Fluo4FF = Kd_Fluo4FF * kon_Fluo4FF

        K01d_N = 10^tv_K01d_N
        k01_N = N_lobe * 10^tv_k01_N
        k10_N = K01d_N * 10^tv_k01_N

        K02d_N = 10^tv_K02d_N
        k02_N = N_lobe * 10^tv_k02_N
        k20_N = K02d_N * 10^tv_k02_N

        K13d_N = 10^tv_K13d_N
        k13_N = N_lobe * 10^tv_k13_N
        k31_N = K13d_N * 10^tv_k13_N

        K23d_N = K01d_N * K13d_N / K02d_N
        k23_N = N_lobe * 10^tv_k23_N
        k32_N = K23d_N * 10^tv_k23_N

        K01d_C = 10^tv_K01d_C
        k01_C = C_lobe * 10^tv_k01_C
        k10_C = K01d_C * 10^tv_k01_C

        K02d_C = 10^tv_K02d_C
        k02_C = C_lobe * 10^tv_k02_C
        k20_C = K02d_C * 10^tv_k02_C

        K13d_C = 10^tv_K13d_C
        k13_C = C_lobe * 10^tv_k13_C
        k31_C = K13d_C * 10^tv_k13_C

        K23d_C = K01d_C * K13d_C / K02d_C
        k23_C = C_lobe * 10^tv_k23_C
        k32_C = K23d_C * 10^tv_k23_C

        Ca₀ = initial_calcium(experiment, k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                              k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                              Ca_free, Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
        DMn₀ = (Kd_DMn * DMn_t) / (Ca₀ + Kd_DMn)
        CaDMn₀ = DMn_t - DMn₀
        OGB5₀ = (Kd_D * OGB5_t) / (Ca₀ + Kd_D)
        CaOGB5₀ = OGB5_t - OGB5₀
        Fluo4FF₀ = (Kd_Fluo4FF * Fluo4FF_t) / (Ca₀ + Kd_Fluo4FF)
        CaFluo4FF₀ = Fluo4FF_t - Fluo4FF₀
        CaM₀ = calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                      k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                                      Ca₀, CaM_t)

        DMn_s₀_u = (1 - f_frac) * DMn₀ * U
        CaDMn_s₀_u = (1 - f_frac) * CaDMn₀ * U
        DMn_f₀_u = f_frac * DMn₀ * U
        CaDMn_f₀_u = f_frac * CaDMn₀ * U
    end

    @DifferentialEquation begin
        JN01_C0(t) = k01_N * Ca * CaM_N0C0 - k10_N * CaM_N1C0
        JN02_C0(t) = k02_N * Ca * CaM_N0C0 - k20_N * CaM_N2C0
        JN13_C0(t) = k13_N * Ca * CaM_N1C0 - k31_N * CaM_N3C0
        JN23_C0(t) = k23_N * Ca * CaM_N2C0 - k32_N * CaM_N3C0
        JN01_C1(t) = k01_N * Ca * CaM_N0C1 - k10_N * CaM_N1C1
        JN02_C1(t) = k02_N * Ca * CaM_N0C1 - k20_N * CaM_N2C1
        JN13_C1(t) = k13_N * Ca * CaM_N1C1 - k31_N * CaM_N3C1
        JN23_C1(t) = k23_N * Ca * CaM_N2C1 - k32_N * CaM_N3C1
        JN01_C2(t) = k01_N * Ca * CaM_N0C2 - k10_N * CaM_N1C2
        JN02_C2(t) = k02_N * Ca * CaM_N0C2 - k20_N * CaM_N2C2
        JN13_C2(t) = k13_N * Ca * CaM_N1C2 - k31_N * CaM_N3C2
        JN23_C2(t) = k23_N * Ca * CaM_N2C2 - k32_N * CaM_N3C2
        JN01_C3(t) = k01_N * Ca * CaM_N0C3 - k10_N * CaM_N1C3
        JN02_C3(t) = k02_N * Ca * CaM_N0C3 - k20_N * CaM_N2C3
        JN13_C3(t) = k13_N * Ca * CaM_N1C3 - k31_N * CaM_N3C3
        JN23_C3(t) = k23_N * Ca * CaM_N2C3 - k32_N * CaM_N3C3
        JC01_N0(t) = k01_C * Ca * CaM_N0C0 - k10_C * CaM_N0C1
        JC02_N0(t) = k02_C * Ca * CaM_N0C0 - k20_C * CaM_N0C2
        JC13_N0(t) = k13_C * Ca * CaM_N0C1 - k31_C * CaM_N0C3
        JC23_N0(t) = k23_C * Ca * CaM_N0C2 - k32_C * CaM_N0C3
        JC01_N1(t) = k01_C * Ca * CaM_N1C0 - k10_C * CaM_N1C1
        JC02_N1(t) = k02_C * Ca * CaM_N1C0 - k20_C * CaM_N1C2
        JC13_N1(t) = k13_C * Ca * CaM_N1C1 - k31_C * CaM_N1C3
        JC23_N1(t) = k23_C * Ca * CaM_N1C2 - k32_C * CaM_N1C3
        JC01_N2(t) = k01_C * Ca * CaM_N2C0 - k10_C * CaM_N2C1
        JC02_N2(t) = k02_C * Ca * CaM_N2C0 - k20_C * CaM_N2C2
        JC13_N2(t) = k13_C * Ca * CaM_N2C1 - k31_C * CaM_N2C3
        JC23_N2(t) = k23_C * Ca * CaM_N2C2 - k32_C * CaM_N2C3
        JC01_N3(t) = k01_C * Ca * CaM_N3C0 - k10_C * CaM_N3C1
        JC02_N3(t) = k02_C * Ca * CaM_N3C0 - k20_C * CaM_N3C2
        JC13_N3(t) = k13_C * Ca * CaM_N3C1 - k31_C * CaM_N3C3
        JC23_N3(t) = k23_C * Ca * CaM_N3C2 - k32_C * CaM_N3C3

        D(Ca) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                -kon_PP * PP * Ca + koff_PP * CaPP +
                -kon_D * OGB5 * Ca + koff_D * CaOGB5 +
                -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF +
                -JN01_C0(t) - JN02_C0(t) - JN13_C0(t) - JN23_C0(t) +
                -JN01_C1(t) - JN02_C1(t) - JN13_C1(t) - JN23_C1(t) +
                -JN01_C2(t) - JN02_C2(t) - JN13_C2(t) - JN23_C2(t) +
                -JN01_C3(t) - JN02_C3(t) - JN13_C3(t) - JN23_C3(t) +
                -JC01_N0(t) - JC02_N0(t) - JC13_N0(t) - JC23_N0(t) +
                -JC01_N1(t) - JC02_N1(t) - JC13_N1(t) - JC23_N1(t) +
                -JC01_N2(t) - JC02_N2(t) - JC13_N2(t) - JC23_N2(t) +
                -JC01_N3(t) - JC02_N3(t) - JC13_N3(t) - JC23_N3(t)
        D(DMn_s) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s - DMn_s / τ_s
        D(CaDMn_s) ~ kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s / τ_s
        D(DMn_f) ~ -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f - DMn_f / τ_f
        D(CaDMn_f) ~ kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f / τ_f
        D(PP) ~ -kon_PP * PP * Ca + koff_PP * CaPP +
                2 * DMn_s / τ_s + 2 * DMn_f / τ_f + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(CaPP) ~ kon_PP * PP * Ca - koff_PP * CaPP + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(OGB5) ~ -kon_D * OGB5 * Ca + koff_D * CaOGB5
        D(CaOGB5) ~ kon_D * OGB5 * Ca - koff_D * CaOGB5
        D(Fluo4FF) ~ -kon_Fluo4FF * Fluo4FF * Ca + koff_Fluo4FF * CaFluo4FF
        D(CaFluo4FF) ~ kon_Fluo4FF * Fluo4FF * Ca - koff_Fluo4FF * CaFluo4FF
        D(CaM_N0C0) ~ -JN01_C0(t) - JN02_C0(t) - JC01_N0(t) - JC02_N0(t)
        D(CaM_N1C0) ~ JN01_C0(t) - JN13_C0(t) - JC01_N1(t) - JC02_N1(t)
        D(CaM_N2C0) ~ JN02_C0(t) - JN23_C0(t) - JC01_N2(t) - JC02_N2(t)
        D(CaM_N3C0) ~ JN13_C0(t) + JN23_C0(t) - JC01_N3(t) - JC02_N3(t)
        D(CaM_N0C1) ~ -JN01_C1(t) - JN02_C1(t) + JC01_N0(t) - JC13_N0(t)
        D(CaM_N1C1) ~ JN01_C1(t) - JN13_C1(t) + JC01_N1(t) - JC13_N1(t)
        D(CaM_N2C1) ~ JN02_C1(t) - JN23_C1(t) + JC01_N2(t) - JC13_N2(t)
        D(CaM_N3C1) ~ JN13_C1(t) + JN23_C1(t) + JC01_N3(t) - JC13_N3(t)
        D(CaM_N0C2) ~ -JN01_C2(t) - JN02_C2(t) + JC02_N0(t) - JC23_N0(t)
        D(CaM_N1C2) ~ JN01_C2(t) - JN13_C2(t) + JC02_N1(t) - JC23_N1(t)
        D(CaM_N2C2) ~ JN02_C2(t) - JN23_C2(t) + JC02_N2(t) - JC23_N2(t)
        D(CaM_N3C2) ~ JN13_C2(t) + JN23_C2(t) + JC02_N3(t) - JC23_N3(t)
        D(CaM_N0C3) ~ -JN01_C3(t) - JN02_C3(t) + JC13_N0(t) + JC23_N0(t)
        D(CaM_N1C3) ~ JN01_C3(t) - JN13_C3(t) + JC13_N1(t) + JC23_N1(t)
        D(CaM_N2C3) ~ JN02_C3(t) - JN23_C3(t) + JC13_N2(t) + JC23_N2(t)
        D(CaM_N3C3) ~ JN13_C3(t) + JN23_C3(t) + JC13_N3(t) + JC23_N3(t)
    end

    @initialDE begin
        DMn_s = DMn_s₀_u
        CaDMn_s = CaDMn_s₀_u
        DMn_f = DMn_f₀_u
        CaDMn_f = CaDMn_f₀_u
        PP = 0.0
        CaPP = 0.0
        OGB5 = OGB5₀
        CaOGB5 = CaOGB5₀
        Fluo4FF = Fluo4FF₀
        CaFluo4FF = CaFluo4FF₀
        CaM_N0C0 = CaM₀[1]
        CaM_N1C0 = CaM₀[2]
        CaM_N2C0 = CaM₀[3]
        CaM_N3C0 = CaM₀[4]
        CaM_N0C1 = CaM₀[5]
        CaM_N1C1 = CaM₀[6]
        CaM_N2C1 = CaM₀[7]
        CaM_N3C1 = CaM₀[8]
        CaM_N0C2 = CaM₀[9]
        CaM_N1C2 = CaM₀[10]
        CaM_N2C2 = CaM₀[11]
        CaM_N3C2 = CaM₀[12]
        CaM_N0C3 = CaM₀[13]
        CaM_N1C3 = CaM₀[14]
        CaM_N2C3 = CaM₀[15]
        CaM_N3C3 = CaM₀[16]
        Ca = Ca₀
    end

    @formulas begin
        OGB5_F_F0 ~ Normal((OGB5(t) + Fmax_Fmin * CaOGB5(t)) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_OGB5_F_F0)
        bound_Ca_per_CaM ~ Normal(bound_calcium([CaM_N0C0(t), CaM_N1C0(t), CaM_N2C0(t), CaM_N3C0(t),
                                                 CaM_N0C1(t), CaM_N1C1(t), CaM_N2C1(t), CaM_N3C1(t),
                                                 CaM_N0C2(t), CaM_N1C2(t), CaM_N2C2(t), CaM_N3C2(t),
                                                 CaM_N0C3(t), CaM_N1C3(t), CaM_N2C3(t), CaM_N3C3(t)]) / CaM_t,
                                  σ_bound_Ca_per_CaM)
        Fluo4FF_fraction_bound ~ Normal(CaFluo4FF(t) / Fluo4FF_t, σ_Fluo4FF_fraction_bound)
    end
end

calmodulin_full = set_solver_config(calmodulin_full; alg = SOLVER, kwargs = (abstol = ABSTOL, reltol = RELTOL))

function faas_data()
    df = select(load_paper_data("Faas_etal_2011"), :series, :x => :t, :y => :OGB5_F_F0,
                :DMn_fast_fraction => :f_frac,
                :DMn_tau_fast_ms => :τ_f,
                :DMn_tau_slow_ms => :τ_s,
                :DMn_kon_per_M_per_ms => ByRow(x -> 1e-6 * x) => :kon_DMn,
                :DMn_koff_per_ms => :koff_DMn,
                :DMn_Kd_M => ByRow(x -> 1e6 * x) => :Kd_DMn,
                :photoproduct_koff_per_ms => :koff_PP,
                :DMn_total_M => ByRow(x -> 1e6 * x) => :DMn_t,
                :CaM_total_M => ByRow(x -> 1e6 * x) => :CaM_t,
                :OGB5N_total_M => ByRow(x -> 1e6 * x) => :OGB5_t,
                :Ca_free_t0_M => ByRow(x -> 1e6 * x) => :Ca_free,
                :OGB5N_kon_per_M_per_ms => ByRow(x -> 1e-6 * x) => :kon_D,
                :OGB5N_koff_per_ms => :koff_D,
                :OGB5N_Fmax_over_Fmin => :Fmax_Fmin)
    return insertcols!(df, :U_group => string.(df.series), :bound_Ca_per_CaM => missing, :Fluo4FF_fraction_bound => missing,
                       :Fluo4FF_t => 0.0, :Kd_Fluo4FF => 13.0, :kon_Fluo4FF => 1.0,
                       :Ca_total => 0.0, :N_lobe => 1.0, :C_lobe => 1.0,
                       :experiment => experiment_id("Faas_etal_2011"))
end

function shifman_data()
    fig1 = load_paper_data("Shifman_etal_2006", "FIG1")
    fig3 = load_paper_data("Shifman_etal_2006", "FIG3")
    d = vcat(fig1[fig1.series .== "CaM", :], fig3[fig3.series .!= "Fluo4FF", :]; cols = :union)
    free = d.figure .== "FIG1"
    return DataFrame(series = string.("Shifman_", d.figure, "_", d.series, "_", 1:nrow(d)), U_group = "Shifman", t = 0.001,
                     OGB5_F_F0 = missing,
                     bound_Ca_per_CaM = ifelse.(free, d.y, missing),
                     Fluo4FF_fraction_bound = ifelse.(free, missing, d.y),
                     f_frac = 0.0, τ_f = 1.0, τ_s = 1.0, kon_DMn = 0.0, koff_DMn = 0.0, Kd_DMn = 1.0, koff_PP = 0.0,
                     DMn_t = 0.0, CaM_t = d.CaM_uM, OGB5_t = 0.0, Ca_free = ifelse.(free, d.x, 0.0),
                     kon_D = 1.0, koff_D = 1.0, Fmax_Fmin = 1.0,
                     Fluo4FF_t = d.Fluo4FF_uM, Kd_Fluo4FF = 13.0, kon_Fluo4FF = 1.0,
                     Ca_total = ifelse.(free, 0.0, d.x),
                     N_lobe = Float64.(d.CaM .!= "CaM-C^WT"), C_lobe = Float64.(d.CaM .!= "CaM-N^WT"),
                     experiment = experiment_id.(ifelse.(free, "Shifman_etal_2006_Fig1B", "Shifman_etal_2006_Fig3D")))
end

function calmodulin_data(version = 3)
    df = vcat(faas_data(), shifman_data())
    return df[df.experiment .<= version, :]
end

datamodel(version = 3, df = calmodulin_data(version)) =
    DataModel((calmodulin_faas, calmodulin_faas_shifman_fig1, calmodulin_full)[version], df; primary_id = :series, time_col = :t)

function shifman_grid(df, n)
    grid = mapreduce(vcat, groupby(df, [:experiment, :N_lobe, :C_lobe, :CaM_t, :Fluo4FF_t]); init = df[1:0, :]) do g
        x = EXPERIMENTS[g.experiment[1]] == "Shifman_etal_2006_Fig3D" ? :Ca_total : :Ca_free
        r = repeat(g[1:1, :], n)
        r[!, x] = collect(range(extrema(g[!, x])..., n))
        return r
    end
    grid.series = string.("Shifman_grid_", 1:nrow(grid))
    return grid
end

with_prediction(df, pred, outcome) =
    innerjoin(dropmissing(df, outcome),
              select(pred[pred.observable .== outcome, :], :id => :series, :time => :t, :prediction);
              on = [:series, :t], order = :left)

calmodulin_variant(N_lobe, C_lobe) = N_lobe == 0 ? "CaM-C^WT" : C_lobe == 0 ? "CaM-N^WT" : "WT CaM"

function plot_calmodulin(res; re_mode = :ebe, n = 100)
    version = maximum(get_df(get_data_model(res)).experiment)
    shifman = filter(:experiment => <=(version), shifman_data())
    df = vcat(faas_data(), shifman_grid(shifman, n))
    pred = predict(res, df; re_mode)
    faas = with_prediction(df, pred, :OGB5_F_F0)
    faas.prefix = first.(rsplit.(faas.series, "_"; limit = 2))
    curve_1b = with_prediction(df, pred, :bound_Ca_per_CaM)
    curve_3d = with_prediction(df, pred, :Fluo4FF_fraction_bound)
    data_1b = dropmissing(shifman, :bound_Ca_per_CaM)
    data_3d = dropmissing(shifman, :Fluo4FF_fraction_bound)
    for d in (curve_3d, data_3d)
        d.variant = calmodulin_variant.(d.N_lobe, d.C_lobe)
    end

    data_visual = visual(Scatter, color = :black)
    prediction_visual = visual(Lines, color = :orange)
    fig = Figure(size = (2400, 1200))
    draw!(fig[1, 1],
          AlgebraOfGraphics.data(faas) * mapping(layout = :prefix, group = :series) *
          (mapping(:t, :OGB5_F_F0) * visual(Lines, color = :black) + mapping(:t, :prediction) * prediction_visual);
          axis = (; xlabel = "t (ms)", ylabel = "OGB-5N F/F₀"), facet = (; linkyaxes = :none))
    if version >= 2
        draw!(fig[1, 2],
              AlgebraOfGraphics.data(data_1b) * mapping(:Ca_free, :bound_Ca_per_CaM) * data_visual +
              AlgebraOfGraphics.data(curve_1b) * mapping(:Ca_free, :prediction) * prediction_visual;
              axis = (; xlabel = "[Ca²⁺] free (µM)", ylabel = "bound Ca²⁺ per CaM"))
    end
    if version >= 3
        draw!(fig[1, 3],
              (AlgebraOfGraphics.data(data_3d) * mapping(:Ca_total, :Fluo4FF_fraction_bound) * data_visual +
               AlgebraOfGraphics.data(curve_3d) * mapping(:Ca_total, :prediction) * prediction_visual) *
              mapping(row = :variant => sorter("WT CaM", "CaM-N^WT", "CaM-C^WT"));
              axis = (; xlabel = "[Ca²⁺] total (µM)", ylabel = "Fluo4FF fraction bound"))
    end
    for (i, title) in enumerate(("Faas et al. 2011", "Shifman et al. 2006, Fig 1B", "Shifman et al. 2006, Fig 3D")[1:version])
        Label(fig[0, i], title; tellwidth = false)
    end
    version > 1 && colsize!(fig.layout, 1, Relative(0.6))
    return fig
end
