### DL: Faas et al (2011) scheme with identical CaM binding sites within lobes.
### DL: Had to re-write from Linkevicius et al. (2025) who had independent lobes,
### DL: since downstream models will need tracking of CaM molecules and 
### DL: model written with independent lobes cannot track that

using NoLimits, Distributions, DataFrames
using LinearAlgebra: kron, Diagonal
using BracketingNonlinearSolve: IntervalNonlinearProblem, ITP, solve
using SynapseExperimentModeling: load_paper_data, SOLVER, ABSTOL, RELTOL
using CairoMakie: Figure, Axis, Label, Scatter, Lines, GridLayout, scatter!, lines!, colsize!, colgap!, Relative,
                  rich, superscript, subscript, Bottom
using AlgebraOfGraphics: AlgebraOfGraphics, mapping, visual, draw!, legend!, sorter
using CairoMakie: Makie
using PairPlots: PairPlots, pairplot

lobe_weights(kon1, koff1, kon2, koff2, Ca) = [one(Ca), 2 * Ca * kon1 / koff1, Ca^2 * kon1 / koff1 * kon2 / koff2]

lobe_bound(w) = (w[2] + 2 * w[3]) / sum(w)

function calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca, total)
    w_N = lobe_weights(kon1_N, koff1_N, kon2_N, koff2_N, Ca)
    w_C = lobe_weights(kon1_C, koff1_C, kon2_C, koff2_C, Ca)
    return total * kron(w_C, w_N) / (sum(w_N) * sum(w_C))
end

const CALMODULIN_CALCIUM = [n_N + n_C for n_C in 0:2 for n_N in 0:2]

bound_calcium(CaM) = sum(CALMODULIN_CALCIUM .* CaM)

const N_LOBE_CALCIUM = [n_N for n_C in 0:2 for n_N in 0:2]

const C_LOBE_CALCIUM = [n_C for n_C in 0:2 for n_N in 0:2]

N_lobe_calcium(CaM) = sum(N_LOBE_CALCIUM .* CaM)

C_lobe_calcium(CaM) = sum(C_LOBE_CALCIUM .* CaM)

const Q10 = 2.2  # 1; Bayley et al. 1984, Biochem Biophys Res Commun 120:185, doi:10.1016/0006-291x(84)91431-1; Martin et al. 1985, Eur J Biochem 151:543, doi:10.1111/j.1432-1033.1985.tb09137.x; Gaertner et al. 2004, J Biol Chem 279:39374, doi:10.1074/jbc.M405352200

const T_REF = 35.0  # °C; Faas et al. 2011, Nat Neurosci 14:301, doi:10.1038/nn.2746

temperature_factor(temperature) = Q10^((temperature - T_REF) / 10)

upper_truncated_normal(μ, σ, upper) = truncated(Normal(μ, σ); upper = upper)
lower_truncated_normal(μ, σ, lower) = truncated(Normal(μ, σ); lower = lower)
truncated_normal(μ, σ, lower, upper) = truncated(Normal(μ, σ); lower = lower, upper = upper)

calcium_excess(Ca, p) = Ca + p[9] * (lobe_bound(lobe_weights(p[1:4]..., Ca)) + lobe_bound(lobe_weights(p[5:8]..., Ca))) +
                        p[11] * Ca / (Ca + p[12]) - p[10]

function free_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF)
    p = [kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF]
    return solve(IntervalNonlinearProblem(calcium_excess, (zero(Ca_total), Ca_total), p), ITP()).u
end

initial_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF) =
    Ca_total > 0 ?
    free_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF) :
    Ca_free

function equilibrium_readout(readout, kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                             Ca_free, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF)
    Ca = initial_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t, Ca_total,
                         Fluo4FF_t, Kd_Fluo4FF)
    N = lobe_bound(lobe_weights(kon1_N, koff1_N, kon2_N, koff2_N, Ca))
    C = lobe_bound(lobe_weights(kon1_C, koff1_C, kon2_C, koff2_C, Ca))
    readout == 2 && return N + C
    return Ca / (Ca + Kd_Fluo4FF)
end

const KD_FLUO4FF = 13.0  # µM; Shifman et al. 2006, PNAS 103:13968, doi:10.1073/pnas.0606433103

const KON_EGTA = 0.0015  # µM^-1 ms^-1; Smith et al. 1984, Anal Biochem 143:188, doi:10.1016/0003-2697(84)90575-X

const KD_EGTA = exp10(-0.69)  # µM; Hellam & Podolsky 1969, J Physiol 200:807, doi:10.1113/jphysiol.1969.sp008723

const KD_QUIN2 = 0.115  # µM; Tsien et al. 1982, J Cell Biol 94:325, doi:10.1083/jcb.94.2.325

const KOFF_QUIN2 = 0.083  # ms^-1; Quast et al. 1984, Biochem Biophys Res Commun 123:604, doi:10.1016/0006-291x(84)90272-9

const T_QUIN2 = 37.0  # °C; Quast et al. 1984, Biochem Biophys Res Commun 123:604, doi:10.1016/0006-291x(84)90272-9

const EA_KOFF_QUIN2 = 22.0  # kJ mol^-1; Quast et al. 1984, Biochem Biophys Res Commun 123:604, doi:10.1016/0006-291x(84)90272-9

const DH_BAPTA = 13.9  # kJ mol^-1; Harrison & Bers 1987, Biochim Biophys Acta 925:133, doi:10.1016/0304-4165(87)90102-4

const KON_BAPTA = 0.45  # µM^-1 ms^-1; Naraghi 1997, Cell Calcium 22:255, doi:10.1016/s0143-4160(97)90064-6

const KOFF_BAPTA = 0.079  # ms^-1; Naraghi 1997, Cell Calcium 22:255, doi:10.1016/s0143-4160(97)90064-6

const T_BAPTA = 22.0  # °C; Naraghi 1997, Cell Calcium 22:255, doi:10.1016/s0143-4160(97)90064-6

const KD_BR2BAPTA = 1.59  # µM; Putkey et al. 2008, J Biol Chem 283:1401, doi:10.1074/jbc.M703831200

const T_BR2BAPTA = 23.0  # °C; Putkey et al. 2003, J Biol Chem 278:49667, doi:10.1074/jbc.C300372200

const DH_BR2BAPTA = 16.9  # kJ mol^-1; Harrison & Bers 1987, Biochim Biophys Acta 925:133, doi:10.1016/0304-4165(87)90102-4

const GAS_CONSTANT = 8.314e-3  # kJ mol^-1 K^-1

arrhenius(k_ref, E, T_ref, temperature) = k_ref * exp(E / GAS_CONSTANT * (1 / (T_ref + 273.15) - 1 / (temperature + 273.15)))

quin2_koff(temperature) = arrhenius(KOFF_QUIN2, EA_KOFF_QUIN2, T_QUIN2, temperature)

quin2_Kd(temperature) = arrhenius(KD_QUIN2, -DH_BAPTA, T_QUIN2, temperature)

quin2_kon(temperature) = quin2_koff(temperature) / quin2_Kd(temperature)

bapta_koff(temperature) = arrhenius(KOFF_BAPTA, EA_KOFF_QUIN2, T_BAPTA, temperature)

bapta_Kd(temperature) = arrhenius(KOFF_BAPTA / KON_BAPTA, -DH_BAPTA, T_BAPTA, temperature)

bapta_kon(temperature) = bapta_koff(temperature) / bapta_Kd(temperature)

br2bapta_kon(temperature) = bapta_kon(temperature)

br2bapta_Kd(temperature) = arrhenius(KD_BR2BAPTA, -DH_BR2BAPTA, T_BR2BAPTA, temperature)

br2bapta_koff(temperature) = br2bapta_kon(temperature) * br2bapta_Kd(temperature)
const OMEGA_PRIOR_SD = 0.2  # log10

function normalized_indicator_bound(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                    CaM_t, Ca_total, Ca_total_first, Ca_total_last, indicator_t, Kd_indicator)
    bound(x) = equilibrium_readout(3, kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                   zero(x), CaM_t, x, indicator_t, Kd_indicator)
    return (bound(Ca_total) - bound(Ca_total_first)) / (bound(Ca_total_last) - bound(Ca_total_first))
end

preloaded_chelator(chelator_t, kon, koff, Ca, apo_start) = apo_start > 0 ? chelator_t * kon * Ca / (kon * Ca + koff) : zero(Ca)


#=
calcium_excess(Ca, p) = Ca + p[15] * Ca / (Ca + p[16]) +
                        p[14] * (lobe_bound(lobe_weights(p[1:6]..., Ca)) + lobe_bound(lobe_weights(p[7:12]..., Ca))) - p[13]

function free_calcium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                      k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                      Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
    p = [k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k01_C, k10_C, k02_C, k20_C, k13_C, k31_C,
         Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF]
    return solve(IntervalNonlinearProblem(calcium_excess, (zero(Ca_total), Ca_total), p), ITP()).u
end
=#

const EXPERIMENTS = ("Faas_etal_2011", "Shifman_etal_2006_Fig1B", "Shifman_etal_2006_Fig3D", "Gaertner_etal_2004_Fig1A",
                     "Gaertner_etal_2004_Fig2A", "Gaertner_etal_2004_Fig3A", "Williams_etal_2023_Fig7b",
                     "Putkey_etal_2008_Fig5A", "Putkey_etal_2008_Fig5B", "Putkey_etal_2008_Fig5D", "Putkey_etal_2008_Fig5E",
                     "Putkey_etal_2003_Fig1A", "Putkey_etal_2003_Fig1B", "Hoffman_etal_2014_Fig2C",
                     "Masada_etal_2012_Fig5A", "Masada_etal_2012_Fig5B", "Masada_etal_2012_Fig5C")

const MODEL_EXPERIMENTS = (calmodulin_full = EXPERIMENTS,)

experiment_id(name) = findfirst(==(name), EXPERIMENTS)

#=
function initial_calcium(experiment, k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                         k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                         Ca_free, Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
    EXPERIMENTS[experiment] == "Shifman_etal_2006_Fig3D" || return Ca_free
    return free_calcium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                        k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C,
                        Ca_total, CaM_t, Fluo4FF_t, Kd_Fluo4FF)
end
=#

#=
calmodulin_faas = @Model begin
    @helpers begin
        calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca, total) =
            calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca, total)
        temperature_factor(temperature) = temperature_factor(temperature)
    end

    @fixedEffects begin
        tv_kon1_N = RealNumber(log10(0.77), upper = log10(60.0), prior = truncated(Normal(log10(0.77), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_kon2_N = RealNumber(log10(32.0), upper = log10(60.0), prior = truncated(Normal(log10(32.0), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_Kd1_N = RealNumber(log10(193.0), lower = log10(0.001), prior = truncated(Normal(log10(193.0), 3); lower = log10(0.001)))    # log10(µM)
        tv_Kd_app_N = RealNumber(log10(sqrt(193.0 * 0.788)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(193.0 * 0.788)), 3); lower = log10(1.0), upper = log10(200.0)))  # log10(µM)
        tv_kon1_C = RealNumber(log10(0.084), upper = log10(60.0), prior = truncated(Normal(log10(0.084), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_kon2_C = RealNumber(log10(0.025), upper = log10(60.0), prior = truncated(Normal(log10(0.025), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_Kd1_C = RealNumber(log10(27.8), lower = log10(0.001), prior = truncated(Normal(log10(27.8), 3); lower = log10(0.001)))      # log10(µM)
        tv_Kd_app_C = RealNumber(log10(sqrt(27.8 * 0.264)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(27.8 * 0.264)), 3); lower = log10(1.0), upper = log10(200.0)))   # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0)                               # 1
        Ω = RealNumber(1.5, scale = :log, lower = 1e-10, upper = 5.0)                 # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0)      # 1
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
        temperature = ConstantCovariate(; constant_on = :series)  # °C
    end

    @randomEffects begin
        U = RandomEffect(LogitNormal(μ, Ω); column = :U_group)  # 1
    end

    @preDifferentialEquation begin
        kon_PP = kon_DMn
        Kd_D = koff_D / kon_D
        Q_T = temperature_factor(temperature)

        Kd_app_N = 10^tv_Kd_app_N

        Kd1_N = 10^tv_Kd1_N
        kon1_N = Q_T * 10^tv_kon1_N
        koff1_N = Q_T * Kd1_N * 10^tv_kon1_N

        Kd2_N = Kd_app_N^2 / Kd1_N
        kon2_N = Q_T * 10^tv_kon2_N
        koff2_N = Q_T * Kd2_N * 10^tv_kon2_N

        Kd_app_C = 10^tv_Kd_app_C

        Kd1_C = 10^tv_Kd1_C
        kon1_C = Q_T * 10^tv_kon1_C
        koff1_C = Q_T * Kd1_C * 10^tv_kon1_C

        Kd2_C = Kd_app_C^2 / Kd1_C
        kon2_C = Q_T * 10^tv_kon2_C
        koff2_C = Q_T * Kd2_C * 10^tv_kon2_C

        DMn₀ = (Kd_DMn * DMn_t) / (Ca_free + Kd_DMn)
        CaDMn₀ = DMn_t - DMn₀
        OGB5₀ = (Kd_D * OGB5_t) / (Ca_free + Kd_D)
        CaOGB5₀ = OGB5_t - OGB5₀
        CaM₀ = calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t)

        DMn_s₀_u = (1 - f_frac) * DMn₀ * U
        CaDMn_s₀_u = (1 - f_frac) * CaDMn₀ * U
        DMn_f₀_u = f_frac * DMn₀ * U
        CaDMn_f₀_u = f_frac * CaDMn₀ * U
    end

    @DifferentialEquation begin
        JN01_C0(t) = 2 * kon1_N * Ca * CaM_N0C0 - koff1_N * CaM_N1C0
        JN12_C0(t) = kon2_N * Ca * CaM_N1C0 - 2 * koff2_N * CaM_N2C0
        JN01_C1(t) = 2 * kon1_N * Ca * CaM_N0C1 - koff1_N * CaM_N1C1
        JN12_C1(t) = kon2_N * Ca * CaM_N1C1 - 2 * koff2_N * CaM_N2C1
        JN01_C2(t) = 2 * kon1_N * Ca * CaM_N0C2 - koff1_N * CaM_N1C2
        JN12_C2(t) = kon2_N * Ca * CaM_N1C2 - 2 * koff2_N * CaM_N2C2
        JC01_N0(t) = 2 * kon1_C * Ca * CaM_N0C0 - koff1_C * CaM_N0C1
        JC12_N0(t) = kon2_C * Ca * CaM_N0C1 - 2 * koff2_C * CaM_N0C2
        JC01_N1(t) = 2 * kon1_C * Ca * CaM_N1C0 - koff1_C * CaM_N1C1
        JC12_N1(t) = kon2_C * Ca * CaM_N1C1 - 2 * koff2_C * CaM_N1C2
        JC01_N2(t) = 2 * kon1_C * Ca * CaM_N2C0 - koff1_C * CaM_N2C1
        JC12_N2(t) = kon2_C * Ca * CaM_N2C1 - 2 * koff2_C * CaM_N2C2

        D(Ca) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                -kon_PP * PP * Ca + koff_PP * CaPP +
                -kon_D * OGB5 * Ca + koff_D * CaOGB5 +
                -JN01_C0(t) - JN12_C0(t) - JN01_C1(t) - JN12_C1(t) - JN01_C2(t) - JN12_C2(t) +
                -JC01_N0(t) - JC12_N0(t) - JC01_N1(t) - JC12_N1(t) - JC01_N2(t) - JC12_N2(t)
        D(DMn_s) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s - DMn_s / τ_s
        D(CaDMn_s) ~ kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s / τ_s
        D(DMn_f) ~ -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f - DMn_f / τ_f
        D(CaDMn_f) ~ kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f / τ_f
        D(PP) ~ -kon_PP * PP * Ca + koff_PP * CaPP +
                2 * DMn_s / τ_s + 2 * DMn_f / τ_f + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(CaPP) ~ kon_PP * PP * Ca - koff_PP * CaPP + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(OGB5) ~ -kon_D * OGB5 * Ca + koff_D * CaOGB5
        D(CaOGB5) ~ kon_D * OGB5 * Ca - koff_D * CaOGB5
        D(CaM_N0C0) ~ -JN01_C0(t) - JC01_N0(t)
        D(CaM_N1C0) ~ JN01_C0(t) - JN12_C0(t) - JC01_N1(t)
        D(CaM_N2C0) ~ JN12_C0(t) - JC01_N2(t)
        D(CaM_N0C1) ~ -JN01_C1(t) + JC01_N0(t) - JC12_N0(t)
        D(CaM_N1C1) ~ JN01_C1(t) - JN12_C1(t) + JC01_N1(t) - JC12_N1(t)
        D(CaM_N2C1) ~ JN12_C1(t) + JC01_N2(t) - JC12_N2(t)
        D(CaM_N0C2) ~ -JN01_C2(t) + JC12_N0(t)
        D(CaM_N1C2) ~ JN01_C2(t) - JN12_C2(t) + JC12_N1(t)
        D(CaM_N2C2) ~ JN12_C2(t) + JC12_N2(t)
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
        CaM_N0C1 = CaM₀[4]
        CaM_N1C1 = CaM₀[5]
        CaM_N2C1 = CaM₀[6]
        CaM_N0C2 = CaM₀[7]
        CaM_N1C2 = CaM₀[8]
        CaM_N2C2 = CaM₀[9]
        Ca = Ca_free
    end

    @formulas begin
        OGB5_F_F0 ~ Normal((OGB5(t) + Fmax_Fmin * CaOGB5(t)) / (OGB5₀ + Fmax_Fmin * CaOGB5₀), σ_OGB5_F_F0)
    end
end

calmodulin_faas = set_solver_config(calmodulin_faas; alg = SOLVER, kwargs = (abstol = ABSTOL, reltol = RELTOL))
=#

#=
calmodulin_faas_shifman_fig1 = @Model begin
    @helpers begin
        calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                               k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total) =
            calmodulin_equilibrium(k01_N, k10_N, k02_N, k20_N, k13_N, k31_N, k23_N, k32_N,
                                   k01_C, k10_C, k02_C, k20_C, k13_C, k31_C, k23_C, k32_C, Ca, total)
        bound_calcium(CaM) = bound_calcium(CaM)
    end

    @fixedEffects begin
        tv_k01_N = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k02_N = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k13_N = RealNumber(log10(0.5072), upper = log10(60.0), prior = truncated(Normal(log10(0.5072), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_k23_N = RealNumber(log10(0.5), upper = log10(60.0), prior = truncated(Normal(log10(0.5), 3); upper = log10(60.0)))          # log10(µM^-1 ms^-1)
        tv_K01d_N = RealNumber(log10(33.0), lower = log10(0.001), prior = truncated(Normal(log10(33.0), 3); lower = log10(0.001)))       # log10(µM)
        tv_Kd_app_N = RealNumber(log10(sqrt(33.0 * 230.0)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(33.0 * 230.0)), 3); lower = log10(1.0), upper = log10(200.0)))  # log10(µM)
        tv_K02d_N = RealNumber(log10(3.45), lower = log10(0.001), prior = truncated(Normal(log10(3.45), 3); lower = log10(0.001)))       # log10(µM)
        tv_k01_C = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k02_C = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k13_C = RealNumber(log10(0.00371), upper = log10(60.0), prior = truncated(Normal(log10(0.00371), 3); upper = log10(60.0)))  # log10(µM^-1 ms^-1)
        tv_k23_C = RealNumber(log10(0.118), upper = log10(60.0), prior = truncated(Normal(log10(0.118), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_K01d_C = RealNumber(log10(18.5), lower = log10(0.001), prior = truncated(Normal(log10(18.5), 3); lower = log10(0.001)))       # log10(µM)
        tv_K02d_C = RealNumber(log10(116.0), lower = log10(0.001), prior = truncated(Normal(log10(116.0), 3); lower = log10(0.001)))     # log10(µM)
        tv_Kd_app_C = RealNumber(log10(sqrt(18.5 * 0.38)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(18.5 * 0.38)), 3); lower = log10(1.0), upper = log10(200.0)))   # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0)                               # 1
        Ω = RealNumber(1.5, scale = :log, lower = 1e-10, upper = 5.0)                 # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0)      # 1
        σ_bound_Ca_per_CaM = RealNumber(0.3, scale = :log, lower = 1e-10, upper = 1.0)  # 1
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

        Kd_app_N = 10^tv_Kd_app_N

        K01d_N = 10^tv_K01d_N
        k01_N = 10^tv_k01_N
        k10_N = K01d_N * 10^tv_k01_N

        K02d_N = 10^tv_K02d_N
        k02_N = 10^tv_k02_N
        k20_N = K02d_N * 10^tv_k02_N

        K13d_N = Kd_app_N^2 / K01d_N
        k13_N = 10^tv_k13_N
        k31_N = K13d_N * 10^tv_k13_N

        K23d_N = Kd_app_N^2 / K02d_N
        k23_N = 10^tv_k23_N
        k32_N = K23d_N * 10^tv_k23_N

        Kd_app_C = 10^tv_Kd_app_C

        K01d_C = 10^tv_K01d_C
        k01_C = 10^tv_k01_C
        k10_C = K01d_C * 10^tv_k01_C

        K02d_C = 10^tv_K02d_C
        k02_C = 10^tv_k02_C
        k20_C = K02d_C * 10^tv_k02_C

        K13d_C = Kd_app_C^2 / K01d_C
        k13_C = 10^tv_k13_C
        k31_C = K13d_C * 10^tv_k13_C

        K23d_C = Kd_app_C^2 / K02d_C
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
=#

#=
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
        tv_k01_N = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k02_N = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k13_N = RealNumber(log10(0.5072), upper = log10(60.0), prior = truncated(Normal(log10(0.5072), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_k23_N = RealNumber(log10(0.5), upper = log10(60.0), prior = truncated(Normal(log10(0.5), 3); upper = log10(60.0)))          # log10(µM^-1 ms^-1)
        tv_K01d_N = RealNumber(log10(33.0), lower = log10(0.001), prior = truncated(Normal(log10(33.0), 3); lower = log10(0.001)))       # log10(µM)
        tv_Kd_app_N = RealNumber(log10(sqrt(33.0 * 230.0)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(33.0 * 230.0)), 3); lower = log10(1.0), upper = log10(200.0)))  # log10(µM)
        tv_K02d_N = RealNumber(log10(3.45), lower = log10(0.001), prior = truncated(Normal(log10(3.45), 3); lower = log10(0.001)))       # log10(µM)
        tv_k01_C = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k02_C = RealNumber(log10(0.275), upper = log10(60.0), prior = truncated(Normal(log10(0.275), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_k13_C = RealNumber(log10(0.00371), upper = log10(60.0), prior = truncated(Normal(log10(0.00371), 3); upper = log10(60.0)))  # log10(µM^-1 ms^-1)
        tv_k23_C = RealNumber(log10(0.118), upper = log10(60.0), prior = truncated(Normal(log10(0.118), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_K01d_C = RealNumber(log10(18.5), lower = log10(0.001), prior = truncated(Normal(log10(18.5), 3); lower = log10(0.001)))       # log10(µM)
        tv_K02d_C = RealNumber(log10(116.0), lower = log10(0.001), prior = truncated(Normal(log10(116.0), 3); lower = log10(0.001)))     # log10(µM)
        tv_Kd_app_C = RealNumber(log10(sqrt(18.5 * 0.38)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(18.5 * 0.38)), 3); lower = log10(1.0), upper = log10(200.0)))   # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0)                               # 1
        Ω = RealNumber(1.5, scale = :log, lower = 1e-10, upper = 5.0)                 # 1
        σ_OGB5_F_F0 = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0)      # 1
        σ_bound_Ca_per_CaM = RealNumber(0.3, scale = :log, lower = 1e-10, upper = 1.0)  # 1
        σ_Fluo4FF_fraction_bound = RealNumber(0.05, scale = :log, lower = 1e-10, upper = 1.0)  # 1
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

        Kd_app_N = 10^tv_Kd_app_N

        K01d_N = 10^tv_K01d_N
        k01_N = N_lobe * 10^tv_k01_N
        k10_N = K01d_N * 10^tv_k01_N

        K02d_N = 10^tv_K02d_N
        k02_N = N_lobe * 10^tv_k02_N
        k20_N = K02d_N * 10^tv_k02_N

        K13d_N = Kd_app_N^2 / K01d_N
        k13_N = N_lobe * 10^tv_k13_N
        k31_N = K13d_N * 10^tv_k13_N

        K23d_N = Kd_app_N^2 / K02d_N
        k23_N = N_lobe * 10^tv_k23_N
        k32_N = K23d_N * 10^tv_k23_N

        Kd_app_C = 10^tv_Kd_app_C

        K01d_C = 10^tv_K01d_C
        k01_C = C_lobe * 10^tv_k01_C
        k10_C = K01d_C * 10^tv_k01_C

        K02d_C = 10^tv_K02d_C
        k02_C = C_lobe * 10^tv_k02_C
        k20_C = K02d_C * 10^tv_k02_C

        K13d_C = Kd_app_C^2 / K01d_C
        k13_C = C_lobe * 10^tv_k13_C
        k31_C = K13d_C * 10^tv_k13_C

        K23d_C = Kd_app_C^2 / K02d_C
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
=#

calmodulin_full = @Model begin
    @helpers begin
        calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca, total) =
            calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca, total)
        temperature_factor(temperature) = temperature_factor(temperature)
        upper_truncated_normal(μ, σ, upper) = upper_truncated_normal(μ, σ, upper)
        lower_truncated_normal(μ, σ, lower) = lower_truncated_normal(μ, σ, lower)
        truncated_normal(μ, σ, lower, upper) = truncated_normal(μ, σ, lower, upper)
        initial_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t, Ca_total,
                        Fluo4FF_t, Kd_Fluo4FF) =
            initial_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t, Ca_total,
                            Fluo4FF_t, Kd_Fluo4FF)
        equilibrium_readout(readout, kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                            Ca_free, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF) =
            equilibrium_readout(readout, kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                Ca_free, CaM_t, Ca_total, Fluo4FF_t, Kd_Fluo4FF)
        normalized_indicator_bound(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                   CaM_t, Ca_total, Ca_total_first, Ca_total_last, indicator_t, Kd_indicator) =
            normalized_indicator_bound(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                       CaM_t, Ca_total, Ca_total_first, Ca_total_last, indicator_t, Kd_indicator)
        preloaded_chelator(chelator_t, kon, koff, Ca, apo_start) = preloaded_chelator(chelator_t, kon, koff, Ca, apo_start)
    end

    @fixedEffects begin
        tv_kon1_N = RealNumber(log10(0.77), upper = log10(60.0), prior = truncated(Normal(log10(0.77), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_kon2_N = RealNumber(log10(32.0), upper = log10(60.0), prior = truncated(Normal(log10(32.0), 3); upper = log10(60.0)))      # log10(µM^-1 ms^-1)
        tv_Kd1_N = RealNumber(log10(193.0), lower = log10(0.001), prior = truncated(Normal(log10(193.0), 3); lower = log10(0.001)))    # log10(µM)
        tv_Kd_app_N = RealNumber(log10(sqrt(193.0 * 0.788)), lower = log10(1.0), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(193.0 * 0.788)), 3); lower = log10(1.0), upper = log10(200.0)))  # log10(µM)
        tv_kon1_C = RealNumber(log10(0.084), upper = log10(60.0), prior = truncated(Normal(log10(0.084), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_kon2_C = RealNumber(log10(0.025), upper = log10(60.0), prior = truncated(Normal(log10(0.025), 3); upper = log10(60.0)))    # log10(µM^-1 ms^-1)
        tv_Kd1_C = RealNumber(log10(27.8), lower = log10(0.001), prior = truncated(Normal(log10(27.8), 3); lower = log10(0.001)))      # log10(µM)
        tv_Kd_app_C = RealNumber(log10(sqrt(27.8 * 0.264)), lower = log10(0.1), upper = log10(200.0),
                                 prior = truncated(Normal(log10(sqrt(27.8 * 0.264)), 3); lower = log10(0.1), upper = log10(200.0)))   # log10(µM)
        μ = RealNumber(-1.0, lower = -5.0, upper = 5.0)                               # 1
        Ω = RealNumber(1.5, scale = :log, lower = 1e-10, upper = 5.0)                 # 1
        σ = RealVector(fill(0.9, 17), scale = :log, lower = fill(1e-10, 17), upper = fill(1.0, 17))  # units of y
        ω = RealVector(fill(0.5, 8), scale = :log, lower = fill(1e-10, 8), upper = fill(1.0, 8),
                       prior = product_distribution(fill(truncated(Normal(0.0, OMEGA_PRIOR_SD); lower = 0.0, upper = 1.0), 8)))  # log10
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
        Kd_Ca_indicator = ConstantCovariate(; constant_on = :series)     # µM
        Ca_indicator_ratio = ConstantCovariate(; constant_on = :series)  # 1
        kon_D = ConstantCovariate(; constant_on = :series)      # µM^-1 ms^-1
        koff_D = ConstantCovariate(; constant_on = :series)     # ms^-1
        Fmax_Fmin = ConstantCovariate(; constant_on = :series)  # 1
        Ca_total = ConstantCovariate(; constant_on = :series)   # µM
        chelator_t = ConstantCovariate(; constant_on = :series)     # µM
        kon_chelator = ConstantCovariate(; constant_on = :series)   # µM^-1 ms^-1
        koff_chelator = ConstantCovariate(; constant_on = :series)  # ms^-1
        chelator2_t = ConstantCovariate(; constant_on = :series)     # µM
        kon_chelator2 = ConstantCovariate(; constant_on = :series)   # µM^-1 ms^-1
        koff_chelator2 = ConstantCovariate(; constant_on = :series)  # ms^-1
        apo_start = ConstantCovariate(; constant_on = :series)       # 1
        Fluo4FF_t = ConstantCovariate(; constant_on = :series)  # µM
        Kd_Fluo4FF = ConstantCovariate(; constant_on = :series) # µM
        N_lobe = ConstantCovariate(; constant_on = :series)     # 1
        C_lobe = ConstantCovariate(; constant_on = :series)     # 1
        temperature = ConstantCovariate(; constant_on = :series)  # °C
        noise = ConstantCovariate(; constant_on = :series)      # 1
        t_first = Covariate()                                   # ms
        t_last = Covariate()                                    # ms
        readout = Covariate()                                   # 1
        Ca_ratio_row = Covariate()                              # 1
        Ca_total_row = Covariate()                              # µM
        Ca_total_first = Covariate()                            # µM
        Ca_total_last = Covariate()                             # µM
        y_ref0 = Covariate()                                    # units of y
        y_ref1 = Covariate()                                    # units of y
    end

    @randomEffects begin
        U = RandomEffect(LogitNormal(μ, Ω); column = :U_group)  # 1
        z_kon1_N = RandomEffect(upper_truncated_normal(0.0, ω[1], log10(60.0) - tv_kon1_N); column = :setup)  # 1
        z_kon2_N = RandomEffect(upper_truncated_normal(0.0, ω[2], log10(60.0) - tv_kon2_N); column = :setup)  # 1
        z_Kd1_N = RandomEffect(lower_truncated_normal(0.0, ω[3], log10(0.001) - tv_Kd1_N); column = :setup)  # 1
        z_Kd_app_N = RandomEffect(truncated_normal(0.0, ω[4], log10(1.0) - tv_Kd_app_N, log10(200.0) - tv_Kd_app_N); column = :setup)  # 1
        z_kon1_C = RandomEffect(upper_truncated_normal(0.0, ω[5], log10(60.0) - tv_kon1_C); column = :setup)  # 1
        z_kon2_C = RandomEffect(upper_truncated_normal(0.0, ω[6], log10(60.0) - tv_kon2_C); column = :setup)  # 1
        z_Kd1_C = RandomEffect(lower_truncated_normal(0.0, ω[7], log10(0.001) - tv_Kd1_C); column = :setup)  # 1
        z_Kd_app_C = RandomEffect(truncated_normal(0.0, ω[8], log10(0.1) - tv_Kd_app_C, log10(200.0) - tv_Kd_app_C); column = :setup)  # 1
    end

    @preDifferentialEquation begin
        kon_PP = kon_DMn
        Kd_D = koff_D / kon_D
        Q_T = temperature_factor(temperature)

        Kd_app_N = 10^(tv_Kd_app_N + z_Kd_app_N)

        Kd1_N = 10^(tv_Kd1_N + z_Kd1_N)
        k1_N = Q_T * 10^(tv_kon1_N + z_kon1_N)
        kon1_N = N_lobe * k1_N
        koff1_N = Kd1_N * k1_N

        Kd2_N = Kd_app_N^2 / Kd1_N
        k2_N = Q_T * 10^(tv_kon2_N + z_kon2_N)
        kon2_N = N_lobe * k2_N
        koff2_N = Kd2_N * k2_N

        Kd_app_C = 10^(tv_Kd_app_C + z_Kd_app_C)

        Kd1_C = 10^(tv_Kd1_C + z_Kd1_C)
        k1_C = Q_T * 10^(tv_kon1_C + z_kon1_C)
        kon1_C = C_lobe * k1_C
        koff1_C = Kd1_C * k1_C

        Kd2_C = Kd_app_C^2 / Kd1_C
        k2_C = Q_T * 10^(tv_kon2_C + z_kon2_C)
        kon2_C = C_lobe * k2_C
        koff2_C = Kd2_C * k2_C

        Ca_free = Kd_Ca_indicator * Ca_indicator_ratio
        Ca₀ = initial_calcium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C, Ca_free, CaM_t, Ca_total,
                              Fluo4FF_t, Kd_Fluo4FF)
        DMn₀ = (Kd_DMn * DMn_t) / (Ca₀ + Kd_DMn)
        CaDMn₀ = DMn_t - DMn₀
        OGB5₀ = (Kd_D * OGB5_t) / (Ca₀ + Kd_D)
        CaOGB5₀ = OGB5_t - OGB5₀
        CaM₀ = calmodulin_equilibrium(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                      (1 - apo_start) * Ca₀, CaM_t)
        CaChelator₀ = preloaded_chelator(chelator_t, kon_chelator, koff_chelator, Ca₀, apo_start)

        DMn_s₀_u = (1 - f_frac) * DMn₀ * U
        CaDMn_s₀_u = (1 - f_frac) * CaDMn₀ * U
        DMn_f₀_u = f_frac * DMn₀ * U
        CaDMn_f₀_u = f_frac * CaDMn₀ * U
    end

    @DifferentialEquation begin
        JN01_C0(t) = 2 * kon1_N * Ca * CaM_N0C0 - koff1_N * CaM_N1C0
        JN12_C0(t) = kon2_N * Ca * CaM_N1C0 - 2 * koff2_N * CaM_N2C0
        JN01_C1(t) = 2 * kon1_N * Ca * CaM_N0C1 - koff1_N * CaM_N1C1
        JN12_C1(t) = kon2_N * Ca * CaM_N1C1 - 2 * koff2_N * CaM_N2C1
        JN01_C2(t) = 2 * kon1_N * Ca * CaM_N0C2 - koff1_N * CaM_N1C2
        JN12_C2(t) = kon2_N * Ca * CaM_N1C2 - 2 * koff2_N * CaM_N2C2
        JC01_N0(t) = 2 * kon1_C * Ca * CaM_N0C0 - koff1_C * CaM_N0C1
        JC12_N0(t) = kon2_C * Ca * CaM_N0C1 - 2 * koff2_C * CaM_N0C2
        JC01_N1(t) = 2 * kon1_C * Ca * CaM_N1C0 - koff1_C * CaM_N1C1
        JC12_N1(t) = kon2_C * Ca * CaM_N1C1 - 2 * koff2_C * CaM_N1C2
        JC01_N2(t) = 2 * kon1_C * Ca * CaM_N2C0 - koff1_C * CaM_N2C1
        JC12_N2(t) = kon2_C * Ca * CaM_N2C1 - 2 * koff2_C * CaM_N2C2

        D(Ca) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s +
                -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f +
                -kon_PP * PP * Ca + koff_PP * CaPP +
                -kon_D * OGB5 * Ca + koff_D * CaOGB5 +
                -kon_chelator * Chelator * Ca + koff_chelator * CaChelator +
                -kon_chelator2 * Chelator2 * Ca + koff_chelator2 * CaChelator2 +
                -JN01_C0(t) - JN12_C0(t) - JN01_C1(t) - JN12_C1(t) - JN01_C2(t) - JN12_C2(t) +
                -JC01_N0(t) - JC12_N0(t) - JC01_N1(t) - JC12_N1(t) - JC01_N2(t) - JC12_N2(t)
        D(DMn_s) ~ -kon_DMn * DMn_s * Ca + koff_DMn * CaDMn_s - DMn_s / τ_s
        D(CaDMn_s) ~ kon_DMn * DMn_s * Ca - koff_DMn * CaDMn_s - CaDMn_s / τ_s
        D(DMn_f) ~ -kon_DMn * DMn_f * Ca + koff_DMn * CaDMn_f - DMn_f / τ_f
        D(CaDMn_f) ~ kon_DMn * DMn_f * Ca - koff_DMn * CaDMn_f - CaDMn_f / τ_f
        D(PP) ~ -kon_PP * PP * Ca + koff_PP * CaPP +
                2 * DMn_s / τ_s + 2 * DMn_f / τ_f + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(CaPP) ~ kon_PP * PP * Ca - koff_PP * CaPP + CaDMn_s / τ_s + CaDMn_f / τ_f
        D(OGB5) ~ -kon_D * OGB5 * Ca + koff_D * CaOGB5
        D(CaOGB5) ~ kon_D * OGB5 * Ca - koff_D * CaOGB5
        D(Chelator) ~ -kon_chelator * Chelator * Ca + koff_chelator * CaChelator
        D(CaChelator) ~ kon_chelator * Chelator * Ca - koff_chelator * CaChelator
        D(Chelator2) ~ -kon_chelator2 * Chelator2 * Ca + koff_chelator2 * CaChelator2
        D(CaChelator2) ~ kon_chelator2 * Chelator2 * Ca - koff_chelator2 * CaChelator2
        D(CaM_N0C0) ~ -JN01_C0(t) - JC01_N0(t)
        D(CaM_N1C0) ~ JN01_C0(t) - JN12_C0(t) - JC01_N1(t)
        D(CaM_N2C0) ~ JN12_C0(t) - JC01_N2(t)
        D(CaM_N0C1) ~ -JN01_C1(t) + JC01_N0(t) - JC12_N0(t)
        D(CaM_N1C1) ~ JN01_C1(t) - JN12_C1(t) + JC01_N1(t) - JC12_N1(t)
        D(CaM_N2C1) ~ JN12_C1(t) + JC01_N2(t) - JC12_N2(t)
        D(CaM_N0C2) ~ -JN01_C2(t) + JC12_N0(t)
        D(CaM_N1C2) ~ JN01_C2(t) - JN12_C2(t) + JC12_N1(t)
        D(CaM_N2C2) ~ JN12_C2(t) + JC12_N2(t)
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
        Chelator = chelator_t - CaChelator₀
        CaChelator = CaChelator₀
        Chelator2 = chelator2_t
        CaChelator2 = 0.0
        CaM_N0C0 = CaM₀[1]
        CaM_N1C0 = CaM₀[2]
        CaM_N2C0 = CaM₀[3]
        CaM_N0C1 = CaM₀[4]
        CaM_N1C1 = CaM₀[5]
        CaM_N2C1 = CaM₀[6]
        CaM_N0C2 = CaM₀[7]
        CaM_N1C2 = CaM₀[8]
        CaM_N2C2 = CaM₀[9]
        Ca = Ca₀
    end

    @formulas begin
        μ_r = readout == 1 ? (OGB5(t) + Fmax_Fmin * CaOGB5(t)) / (OGB5₀ + Fmax_Fmin * CaOGB5₀) :
              readout == 4 ? (CaChelator(t) - CaChelator(t_first)) / (CaChelator(t_last) - CaChelator(t_first) + 1e-6) :  # µM
              readout == 5 ? normalized_indicator_bound(kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                                         CaM_t, Ca_total_row, Ca_total_first, Ca_total_last,
                                                         Fluo4FF_t, Kd_Fluo4FF) :
              equilibrium_readout(readout, kon1_N, koff1_N, kon2_N, koff2_N, kon1_C, koff1_C, kon2_C, koff2_C,
                                  Kd_Ca_indicator * Ca_ratio_row, CaM_t, Ca_total_row, Fluo4FF_t, Kd_Fluo4FF)
        μ_y = y_ref0 + (y_ref1 - y_ref0) * μ_r
        y ~ Normal(μ_y, σ[Int(noise)])
    end
end

calmodulin_full = set_solver_config(calmodulin_full; alg = SOLVER, kwargs = (abstol = ABSTOL, reltol = RELTOL))

const READOUT = (OGB5_F_F0 = 1, bound_Ca_per_CaM = 2, Fluo4FF_fraction_bound = 3, chelator_bound_normalized = 4,
                  indicator_fraction_normalized = 5)

const KINETIC_READOUTS = (READOUT.chelator_bound_normalized,)

const CONSTANTS_RE = (; U = (; none = 0.5))  # 1

const NOISE = (Faas = 1, Shifman_Fig1B = 2, Shifman_Fig3D = 3, Gaertner_Fig1A = 4, Gaertner_Fig2A = 5, Gaertner_Fig3A = 6,
               Williams_Fig7b = 7, Putkey08_Fig5A = 8, Putkey08_Fig5B = 9, Putkey08_Fig5D = 10, Putkey08_Fig5E = 11,
               Putkey03_Fig1A = 12, Putkey03_Fig1B = 13, Hoffman_Fig2C = 14, Masada_Fig5A = 15, Masada_Fig5B = 16,
               Masada_Fig5C = 17)

const DEFAULT_COVARIATES =(f_frac = 0.0, τ_f = 1.0, τ_s = 1.0, kon_DMn = 0.0, koff_DMn = 0.0, Kd_DMn = 1.0, koff_PP = 0.0,
                            DMn_t = 0.0, OGB5_t = 0.0, kon_D = 1.0, koff_D = 1.0, Fmax_Fmin = 1.0,
                            Ca_free = 0.0, Kd_Ca_indicator = 1.0, Ca_indicator_ratio = 0.0, Ca_total = 0.0,
                            chelator_t = 0.0, kon_chelator = 0.0, koff_chelator = 0.0, chelator2_t = 0.0,
                            kon_chelator2 = 0.0, koff_chelator2 = 0.0, apo_start = 0.0, Fluo4FF_t = 0.0, Kd_Fluo4FF = KD_FLUO4FF, N_lobe = 1.0, C_lobe = 1.0,
                            t_first = 0.0, t_last = 1.0, Ca_ratio_row = 0.0, Ca_total_row = 0.0,
                            Ca_total_first = 0.0, Ca_total_last = 1.0, y_ref0 = 0.0, y_ref1 = 1.0)

function experiment_rows(n; columns...)
    df = DataFrame()
    for (name, value) in pairs(merge(DEFAULT_COVARIATES, (; columns...)))
        df[!, name] = value isa AbstractVector ? value : fill(value, n)
    end
    return df
end

function index_rows!(df)
    sort!(df, :series; alg = Base.Sort.DEFAULT_STABLE)
    return transform!(groupby(df, :series), :series => (s -> Float64.(eachindex(s))) => :t)
end

function faas_data()
    d = load_paper_data("Faas_etal_2011")
    Ca_free = 1e6 .* d.Ca_free_t0_M
    Kd = 1e6 .* d.OGB5N_Kd_M
    return experiment_rows(nrow(d); series = d.series, U_group = string.(d.series), setup = "Faas",
                           experiment = experiment_id("Faas_etal_2011"), noise = NOISE.Faas,
                           readout = READOUT.OGB5_F_F0, t = d.x, y = d.y,
                           f_frac = d.DMn_fast_fraction, τ_f = d.DMn_tau_fast_ms, τ_s = d.DMn_tau_slow_ms,
                           kon_DMn = 1e-6 .* d.DMn_kon_per_M_per_ms, koff_DMn = d.DMn_koff_per_ms,
                           Kd_DMn = 1e6 .* d.DMn_Kd_M, koff_PP = d.photoproduct_koff_per_ms,
                           DMn_t = 1e6 .* d.DMn_total_M, CaM_t = 1e6 .* d.CaM_total_M, OGB5_t = 1e6 .* d.OGB5N_total_M,
                           Ca_free = Ca_free, Kd_Ca_indicator = Kd, Ca_indicator_ratio = Ca_free ./ Kd,
                           kon_D = 1e-6 .* d.OGB5N_kon_per_M_per_ms, koff_D = d.OGB5N_koff_per_ms,
                           Fmax_Fmin = d.OGB5N_Fmax_over_Fmin, temperature = d.temperature_C)
end

function shifman_data()
    fig1 = load_paper_data("Shifman_etal_2006", "FIG1")
    fig3 = load_paper_data("Shifman_etal_2006", "FIG3")
    d = vcat(fig1[fig1.series .== "CaM", :], fig3[fig3.series .!= "Fluo4FF", :]; cols = :union)
    free = d.figure .== "FIG1"
    Ca_free = ifelse.(free, d.x, 0.0)
    return index_rows!(experiment_rows(nrow(d); series = string.("Shifman_", d.figure, "_", d.CaM),
                           U_group = "none", setup = string.("Shifman_", d.CaM),
                           experiment = experiment_id.(ifelse.(free, "Shifman_etal_2006_Fig1B", "Shifman_etal_2006_Fig3D")),
                           noise = ifelse.(free, NOISE.Shifman_Fig1B, NOISE.Shifman_Fig3D),
                           readout = ifelse.(free, READOUT.bound_Ca_per_CaM, READOUT.Fluo4FF_fraction_bound), y = d.y,
                           CaM_t = d.CaM_uM, Ca_free = Ca_free, Kd_Ca_indicator = KD_FLUO4FF,
                           Ca_ratio_row = Ca_free ./ KD_FLUO4FF, Ca_total_row = ifelse.(free, 0.0, d.x),
                           Fluo4FF_t = d.Fluo4FF_uM, N_lobe = Float64.(d.CaM .!= "CaM-C^WT"),
                           C_lobe = Float64.(d.CaM .!= "CaM-N^WT"), temperature = d.temperature_C))
end

function gaertner_data()
    fig1 = load_paper_data("Gaertner_etal_2004", "FIG1")
    fig2 = load_paper_data("Gaertner_etal_2004", "FIG2")
    fig3 = load_paper_data("Gaertner_etal_2004", "FIG3")
    quin2(d, figure, noise) =
        experiment_rows(nrow(d); series = "Gaertner_FIG" * figure * "_CaM", U_group = "none", setup = "Gaertner_Mg",
                        experiment = experiment_id("Gaertner_etal_2004_Fig" * figure), noise,
                        readout = READOUT.chelator_bound_normalized, t = 1000 .* d.x, y = d.y,
                        y_ref0 = d.y[1], y_ref1 = d.y[end], t_first = 1000 * d.x[1], t_last = 1000 * d.x[end], CaM_t = d.CaM_uM, Ca_total = d.Ca_total_uM,
                        chelator_t = d.quin2_uM, kon_chelator = quin2_kon.(d.temperature_C),
                        koff_chelator = quin2_koff.(d.temperature_C),
                        temperature = d.temperature_C)
    fig1a = quin2(sort(fig1[in(("A", "A_break")).(fig1.panel) .& (fig1.series .== "CaM"), :], :x), "1A",
                  NOISE.Gaertner_Fig1A)
    fig3a = quin2(sort(fig3[fig3.series .== "CaM", :], :x), "3A", NOISE.Gaertner_Fig3A)
    e = sort(fig2[fig2.series .== "CaM", :], :x)
    fig2a = index_rows!(experiment_rows(nrow(e); series = "Gaertner_FIG2A_CaM", U_group = "none", setup = "Putkey_Gaertner",
                                        experiment = experiment_id("Gaertner_etal_2004_Fig2A"), noise = NOISE.Gaertner_Fig2A,
                                        readout = READOUT.indicator_fraction_normalized, y = e.y,
                                        y_ref0 = e.y[1], y_ref1 = e.y[end], CaM_t = e.CaM_uM, Ca_total_row = e.x, Ca_total_first = e.x[1],
                                        Ca_total_last = e.x[end], Fluo4FF_t = e.Br2BAPTA_uM,
                                        Kd_Fluo4FF = e.Br2BAPTA_Kd_uM, temperature = e.temperature_C))
    return vcat(fig1a, fig3a, fig2a)
end


function williams_data()
    fig7 = load_paper_data("Williams_etal_2023", "FIG7")
    d = sort(fig7[fig7.panel .== "b", :], :x)
    return experiment_rows(nrow(d); series = "Williams_FIG7b_N53I", U_group = "none", setup = "Williams_N53I_stopped_flow",
                           experiment = experiment_id("Williams_etal_2023_Fig7b"), noise = NOISE.Williams_Fig7b,
                           readout = READOUT.chelator_bound_normalized, t = 1000 .* d.x, y = d.y, y_ref0 = d.y[1],
                           y_ref1 = d.y[end], t_first = 1000 * d.x[1], t_last = 1000 * d.x[end], CaM_t = d.CaM_uM,
                           Ca_total = d.Ca_total_uM, chelator_t = d.quin2_uM, kon_chelator = quin2_kon.(d.temperature_C),
                           koff_chelator = quin2_koff.(d.temperature_C), temperature = d.temperature_C)
end

function titration_rows(d, series, setup, experiment, noise)
    d = sort(d, :x)
    return index_rows!(experiment_rows(nrow(d); series, U_group = "none", setup, experiment = experiment_id(experiment), noise,
                                       readout = READOUT.indicator_fraction_normalized, y = d.y, y_ref0 = d.y[1],
                                       y_ref1 = d.y[end], CaM_t = d.CaM_uM, Ca_total_row = d.x, Ca_total_first = d.x[1],
                                       Ca_total_last = d.x[end], Fluo4FF_t = d.Br2BAPTA_uM, Kd_Fluo4FF = d.Br2BAPTA_Kd_uM,
                                       temperature = d.temperature_C))
end

function trace_rows(d, series, setup, experiment, noise; columns...)
    d = sort(d[.!d.merged, :], :x)
    return experiment_rows(nrow(d); series, U_group = "none", setup, experiment = experiment_id(experiment), noise,
                           readout = READOUT.chelator_bound_normalized, t = 1000 .* d.x, y = d.y, y_ref0 = d.y[1],
                           y_ref1 = d.y[end], t_first = 1000 * d.x[1], t_last = 1000 * d.x[end], temperature = d.temperature_C,
                           columns...)
end

quin2_rows(d, series, setup, experiment, noise; columns...) =
    trace_rows(d, series, setup, experiment, noise; CaM_t = d.CaM_uM[1], Ca_total = d.Ca_total_uM[1], chelator_t = d.quin2_uM[1],
               kon_chelator = quin2_kon(d.temperature_C[1]), koff_chelator = quin2_koff(d.temperature_C[1]), columns...)

association_rows(d, series, setup, experiment, noise) =
    trace_rows(d, series, setup, experiment, noise; CaM_t = d.CaM_uM[1] / 2, Kd_Ca_indicator = 1.0,
               Ca_indicator_ratio = d.Ca_free_uM[1], apo_start = 1.0, chelator_t = d.Br2BAPTA_uM[1] / 2,
               kon_chelator = br2bapta_kon(d.temperature_C[1]), koff_chelator = br2bapta_koff(d.temperature_C[1]),
               chelator2_t = d.BAPTA_uM[1] / 2, kon_chelator2 = bapta_kon(d.temperature_C[1]),
               koff_chelator2 = bapta_koff(d.temperature_C[1]))

function putkey08_data()
    fig5 = load_paper_data("Putkey_etal_2008", "FIG5")
    panel(names...) = fig5[in(names).(fig5.panel), :]
    fig5a = titration_rows(panel("A", "A_tail"), "Putkey08_FIG5A_CaM", "Putkey_Gaertner", "Putkey_etal_2008_Fig5A",
                           NOISE.Putkey08_Fig5A)
    fig5b = quin2_rows(panel("B", "B_late"), "Putkey08_FIG5B_CaM", "Putkey_Gaertner", "Putkey_etal_2008_Fig5B",
                       NOISE.Putkey08_Fig5B)
    fig5d = [association_rows(g, "Putkey08_FIG5D_" * g.series[1], "Putkey_Gaertner", "Putkey_etal_2008_Fig5D",
                              NOISE.Putkey08_Fig5D)
             for g in groupby(panel("D"), :series)]
    fig5e = association_rows(panel("E"), "Putkey08_FIG5E_CaM", "Putkey_Gaertner", "Putkey_etal_2008_Fig5E",
                             NOISE.Putkey08_Fig5E)
    return vcat(fig5a, fig5b, fig5d..., fig5e)
end

function putkey03_data()
    fig1 = load_paper_data("Putkey_etal_2003", "FIG1")
    panel(names...) = fig1[in(names).(fig1.panel), :]
    a = panel("A", "A_break")
    fig1a = quin2_rows(a[a.x .>= 0.002, :], "Putkey03_FIG1A_CaM", "Putkey_Gaertner", "Putkey_etal_2003_Fig1A",
                       NOISE.Putkey03_Fig1A)
    fig1b = titration_rows(panel("B"), "Putkey03_FIG1B_CaM", "Putkey_Gaertner", "Putkey_etal_2003_Fig1B", NOISE.Putkey03_Fig1B)
    return vcat(fig1a, fig1b)
end

hoffman_data() = quin2_rows(load_paper_data("Hoffman_etal_2014", "FIG2"), "Hoffman_FIG2C_CaM", "Hoffman",
                            "Hoffman_etal_2014_Fig2C", NOISE.Hoffman_Fig2C)

function masada_data()
    fig5 = load_paper_data("Masada_etal_2012", "FIG5")
    fig5 = fig5[fig5.x .>= 0.005, :]
    rows(panel, figure, variant; lobes...) =
        quin2_rows(fig5[fig5.panel .== panel, :], "Masada_FIG5" * figure * "_" * variant, "Masada_" * variant,
                   "Masada_etal_2012_Fig5" * figure, NOISE[Symbol("Masada_Fig5", figure)]; lobes...)
    return vcat(rows("A", "A", "WT"), rows("B", "B", "CaM12"; N_lobe = 0.0), rows("C_inset", "C", "CaM34"; C_lobe = 0.0))
end

function calmodulin_data(model = "calmodulin_full")
    df = vcat(faas_data(), shifman_data(), gaertner_data(), williams_data(), putkey08_data(), putkey03_data(), hoffman_data(),
              masada_data())
    return df[in(experiment_id.(MODEL_EXPERIMENTS[Symbol(model)])).(df.experiment), :]
end

datamodel(model = "calmodulin_full", df = calmodulin_data(model)) =
    DataModel(getfield(@__MODULE__, Symbol(model)), df; primary_id = :series, time_col = :t)

function curve_grid(df, n)
    grid = mapreduce(vcat, groupby(df, [:experiment, :setup, :readout]); init = df[1:0, :]) do g
        x = g.readout[1] in (READOUT.Fluo4FF_fraction_bound, READOUT.indicator_fraction_normalized) ? :Ca_total_row : :Ca_free
        lo, hi = extrema(g[!, x])
        r = repeat(g[1:1, :], n)
        r[!, x] = collect(range(lo, hi, n))
        r.series .= string("grid_", g.setup[1], "_", g.experiment[1], "_", g.readout[1])
        r.t = Float64.(1:n)
        return r
    end
    grid.Ca_ratio_row = grid.Ca_free ./ grid.Kd_Ca_indicator
    return grid
end

function time_grid(df, n)
    nrow(df) == 0 && return df
    return mapreduce(vcat, groupby(df, :series)) do g
        r = repeat(g[1:1, :], n)
        r.t = exp10.(range(log10(minimum(g.t)), log10(maximum(g.t)), n))
        r.series .= "grid_" * g.series[1]
        return r
    end
end

with_prediction(df, pred) =
    innerjoin(df, select(pred, :id => :series, :time => :t, :prediction); on = [:series, :t], order = :left)

rows_of(df, experiments...) = df[in(experiment_id.(experiments)).(df.experiment), :]

calmodulin_variant(N_lobe, C_lobe) = N_lobe == 0 ? "CaM-C^WT" : C_lobe == 0 ? "CaM-N^WT" : "WT CaM"

function plot_calmodulin(res; re_mode = :ebe, n = 100)
    experiments = EXPERIMENTS[unique(get_df(get_data_model(res)).experiment)]
    fitted(e) = EXPERIMENTS[e] in experiments
    shifman = filter(:experiment => fitted, shifman_data())
    gaertner = filter(:experiment => fitted, gaertner_data())
    putkey08 = filter(:experiment => fitted, putkey08_data())
    putkey03 = filter(:experiment => fitted, putkey03_data())
    hoffman = filter(:experiment => fitted, hoffman_data())
    masada = filter(:experiment => fitted, masada_data())
    kinetic(d) = d[in(KINETIC_READOUTS).(d.readout), :]
    equilibrium(d) = d[.!in(KINETIC_READOUTS).(d.readout), :]
    df = vcat(filter(:experiment => fitted, faas_data()), curve_grid(shifman, n), curve_grid(equilibrium(gaertner), n),
              kinetic(gaertner),
              filter(:experiment => fitted, williams_data()), time_grid(rows_of(gaertner, "Gaertner_etal_2004_Fig1A"), n),
              curve_grid(equilibrium(putkey08), n), kinetic(putkey08), curve_grid(equilibrium(putkey03), n),
              kinetic(putkey03), time_grid(rows_of(putkey03, "Putkey_etal_2003_Fig1A"), 2n), hoffman, time_grid(hoffman, n),
              masada)
    curves = with_prediction(df, predict(res, df; re_mode))
    faas = rows_of(curves, "Faas_etal_2011")
    faas.prefix = first.(rsplit.(faas.series, "_"; limit = 2))
    data_1b = rows_of(shifman, "Shifman_etal_2006_Fig1B")
    curve_1b = rows_of(curves, "Shifman_etal_2006_Fig1B")
    data_3d = rows_of(shifman, "Shifman_etal_2006_Fig3D")
    curve_3d = rows_of(curves, "Shifman_etal_2006_Fig3D")
    for d in (curve_3d, data_3d)
        d.variant = calmodulin_variant.(d.N_lobe, d.C_lobe)
    end
    seconds = :t => (t -> t / 1000)
    moles_released = rich("Moles Ca", superscript("2+"), " released/Mole CaM")

    marker_style = (; markersize = 10, alpha = 0.6, strokecolor = :black, strokewidth = 1)
    trace_style = (; color = :black, alpha = 0.5)
    data_visual = visual(Scatter; color = :black, marker_style...)
    trace_visual = visual(Lines; trace_style...)
    prediction_visual = visual(Lines, color = :orange)
    overlaid(d_data, d_curve, x, y_data, y_curve, data_layer = visual(Scatter; marker_style...)) =
        AlgebraOfGraphics.data(d_data) * mapping(x, y_data) * data_layer +
        AlgebraOfGraphics.data(d_curve) * mapping(x, y_curve) * visual(Lines)
    legend_top!(grid) = legend!(fig[0, length(titles), Bottom()], grid; orientation = :horizontal, tellwidth = false,
                                tellheight = true, framevisible = false)
    fig = Figure(size = (2940, 1440))
    titles = ["Faas et al. 2011"]
    draw!(fig[1:2, 1],
          AlgebraOfGraphics.data(faas) * mapping(layout = :prefix, group = :series) *
          (mapping(:t, :y) * visual(Lines, color = :black) + mapping(:t, :prediction) * prediction_visual);
          axis = (; xlabel = "t (ms)", ylabel = "OGB-5N F/F₀"), facet = (; linkyaxes = :none))
    shifman_1b = "Shifman_etal_2006_Fig1B" in experiments
    shifman_3d = "Shifman_etal_2006_Fig3D" in experiments
    if shifman_1b || shifman_3d
        push!(titles, "Shifman et al. 2006")
        if shifman_1b
            draw!(fig[shifman_3d ? 1 : (1:2), length(titles)],
                  AlgebraOfGraphics.data(data_1b) * mapping(:Ca_free, :y) * data_visual +
                  AlgebraOfGraphics.data(curve_1b) * mapping(:Ca_free, :prediction) * prediction_visual;
                  axis = (; title = "Fig 1B", xlabel = rich("[Ca", superscript("2+"), subscript("free"), "], µM"),
                          ylabel = rich("Bound Ca", superscript("2+"), " per CaM"), limits = (0.095, 54.98, -0.095, 5.003)))
        end
        if shifman_3d
            panel = fig[shifman_1b ? 2 : (1:2), length(titles)] = GridLayout()
            layer = overlaid(data_3d, curve_3d, :Ca_total_row, :y, :prediction) *
                    mapping(color = :variant => sorter("WT CaM", "CaM-N^WT", "CaM-C^WT") => "CaM")
            grid = draw!(panel[1, 1], layer;
                         axis = (; title = "Fig 3D", xlabel = rich("[Ca", superscript("2+"), "], µM"), ylabel = "Fraction Bound",
                                 limits = (-0.1, 92.24, 0.0, 1.052), xticks = 0:20:80))
            draw!(panel[1, 2], layer; axis = (; xlabel = "", ylabel = "", yticklabelsvisible = false,
                                              limits = (204.0, 240.1, 0.0, 1.052), xticks = 220:20:240))
            colsize!(panel, 2, Relative(0.28))
            colgap!(panel, 6)
            legend_top!(grid)
        end
    end
    trace_panel!(position, experiment; axis...) =
        draw!(position, AlgebraOfGraphics.data(rows_of(curves, experiment)) *
                        (mapping(seconds, :y) * trace_visual + mapping(seconds, :prediction) * prediction_visual);
              axis = (; title = "Fig " * last(split(experiment, "_Fig")), axis...))
    if "Williams_etal_2023_Fig7b" in experiments
        push!(titles, "Williams et al. 2023")
        trace_panel!(fig[1:2, length(titles)], "Williams_etal_2023_Fig7b"; xlabel = "Time (sec)",
                     ylabel = "Fluorescence (mV)", limits = (0.003, 3.028, 6.999, 8.231))
    end
    gaertner_figures = filter(in(experiments), ["Gaertner_etal_2004_Fig1A", "Gaertner_etal_2004_Fig2A",
                                                "Gaertner_etal_2004_Fig3A"])
    if !isempty(gaertner_figures)
        push!(titles, "Gaertner et al. 2004")
        column = fig[1:2, length(titles)] = GridLayout()
        for (k, experiment) in enumerate(gaertner_figures)
            if experiment == "Gaertner_etal_2004_Fig2A"
                draw!(column[k, 1],
                      AlgebraOfGraphics.data(rows_of(gaertner, experiment)) * mapping(:Ca_total_row, :y) * data_visual +
                      AlgebraOfGraphics.data(rows_of(curves, experiment)) * mapping(:Ca_total_row, :prediction) *
                      prediction_visual;
                      axis = (; title = "Fig 2A", xlabel = rich("[Ca", superscript("2+"), "] (µM)"),
                              ylabel = rich("Br", subscript("2"), "BAPTA absorbance, scaled"),
                              limits = (0.25, 300.2, -0.038, 1.019)))
            elseif experiment == "Gaertner_etal_2004_Fig1A"
                d = filter(:series => !startswith("grid_"), rows_of(curves, experiment))
                g = filter(:series => startswith("grid_"), rows_of(curves, experiment))
                panel = column[k, 1] = GridLayout()
                broken = [Axis(panel[1, 1]; title = "Fig 1A", xlabel = "Time (s)",
                               ylabel = moles_released, limits = (-0.005, 0.255, 0.0, 4.1),
                               xticks = (0:0.05:0.25, ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"]), yticks = 0:4),
                          Axis(panel[1, 2]; limits = (0.98, 1.005, 0.0, 4.1), xticks = ([1.0], ["1.00"]), yticks = 0:4,
                               yticklabelsvisible = false)]
                colsize!(panel, 2, Relative(0.1))
                colgap!(panel, 6)
                for ax in broken
                    lines!(ax, d.t ./ 1000, d.y; trace_style...)
                    lines!(ax, g.t ./ 1000, g.prediction; color = :orange)
                end
            else
                trace_panel!(column[k, 1], experiment; xlabel = "Time (s)", ylabel = moles_released,
                             limits = (-0.192, 4.003, -0.0134, 4.201))
            end
        end
    end
    titration_panel!(position, data, experiment; axis...) =
        draw!(position, AlgebraOfGraphics.data(rows_of(data, experiment)) * mapping(:Ca_total_row, :y) * data_visual +
                        AlgebraOfGraphics.data(rows_of(curves, experiment)) * mapping(:Ca_total_row, :prediction) *
                        prediction_visual;
              axis = (; title = "Fig " * last(split(experiment, "_Fig")), axis...))
    putkey08_figures = filter(in(experiments), ["Putkey_etal_2008_Fig5A", "Putkey_etal_2008_Fig5B", "Putkey_etal_2008_Fig5D",
                                                "Putkey_etal_2008_Fig5E"])
    if !isempty(putkey08_figures)
        push!(titles, "Putkey et al. 2008")
        column = fig[1:2, length(titles)] = GridLayout()
        for (k, experiment) in enumerate(putkey08_figures)
            if experiment == "Putkey_etal_2008_Fig5A"
                titration_panel!(column[k, 1], putkey08, experiment; xlabel = "Total Calcium (µM)",
                                 ylabel = rich("Relative Br", subscript("2"), "BAPTA Absorbance"),
                                 limits = (-1.06, 248.57, -0.0049, 0.9965), xticks = 50:50:250, yticks = 0.2:0.2:1.0)
            elseif experiment == "Putkey_etal_2008_Fig5D"
                d = rows_of(curves, experiment)
                d.trace = replace.(d.series, "Putkey08_FIG5D_" => "")
                grid = draw!(column[k, 1],
                             overlaid(d, d, seconds, :y, :prediction, visual(Lines, alpha = 0.5)) *
                             mapping(color = :trace => "free Ca²⁺");
                             axis = (; title = "Fig 5D", xlabel = "Time (sec)", ylabel = "Change in Absorbance",
                                     limits = (-0.00071, 0.30071, 0.0, 0.025), xticks = 0.05:0.05:0.25,
                                     yticks = 0:0.005:0.025))
                legend_top!(grid)
            elseif experiment == "Putkey_etal_2008_Fig5B"
                trace_panel!(column[k, 1], experiment; xlabel = "Time (sec)", ylabel = "Change in Fluorescence",
                             limits = (-0.002, 0.998, -0.00128, 0.37969), xticks = 0:0.2:1, yticks = 0.05:0.05:0.30)
            else
                trace_panel!(column[k, 1], experiment; xlabel = "Time (sec)", ylabel = "Change in Absorbance",
                             limits = (0.00147, 0.19764, 3.2e-5, 0.02198), xticks = 0:0.05:0.20, yticks = 0:0.005:0.020)
            end
        end
    end
    putkey03_figures = filter(in(experiments), ["Putkey_etal_2003_Fig1A", "Putkey_etal_2003_Fig1B"])
    if !isempty(putkey03_figures)
        push!(titles, "Putkey et al. 2003")
        column = fig[1:2, length(titles)] = GridLayout()
        for (k, experiment) in enumerate(putkey03_figures)
            if experiment == "Putkey_etal_2003_Fig1A"
                d = filter(:series => !startswith("grid_"), rows_of(curves, experiment))
                g = filter(:series => startswith("grid_"), rows_of(curves, experiment))
                panel = column[k, 1] = GridLayout()
                broken = [Axis(panel[1, 1]; title = "Fig 1A", xlabel = "Time (sec)", ylabel = "Calcium Release (moles/mole CaM)",
                               limits = (-0.00016, 0.02037, 0.0, 2.3356), xticks = (0:0.01:0.02, ["0", "0.01", "0.02"]),
                               yticks = 0:0.5:2.0),
                          Axis(panel[1, 2]; limits = (0.43141, 0.55064, 0.0, 2.3356), xticks = ([0.45, 0.50], ["0.45", "0.50"]),
                               yticks = 0:0.5:2.0, yticklabelsvisible = false)]
                colsize!(panel, 2, Relative(0.49))
                colgap!(panel, 6)
                early, late = d[d.t .< 100, :], d[d.t .>= 100, :]
                lines!(broken[1], early.t ./ 1000, early.y; trace_style...)
                scatter!(broken[2], late.t ./ 1000, late.y; color = :black, marker_style...)
                for ax in broken
                    lines!(ax, g.t ./ 1000, g.prediction; color = :orange)
                end
            else
                titration_panel!(column[k, 1], putkey03, experiment; xlabel = "Total Calcium Concentration (µM)",
                                 ylabel = "Relative Absorbance", limits = (-0.85, 200.17, 0.0, 1.0015), xticks = 0:50:200,
                                 yticks = 0.2:0.2:0.8)
            end
        end
    end
    if "Hoffman_etal_2014_Fig2C" in experiments
        push!(titles, "Hoffman et al. 2014")
        d = filter(:series => !startswith("grid_"), rows_of(curves, "Hoffman_etal_2014_Fig2C"))
        g = filter(:series => startswith("grid_"), rows_of(curves, "Hoffman_etal_2014_Fig2C"))
        draw!(fig[1:2, length(titles)],
              AlgebraOfGraphics.data(d) * mapping(seconds, :y) * data_visual +
              AlgebraOfGraphics.data(g) * mapping(seconds, :prediction) * prediction_visual;
              axis = (; title = "Fig 2C", xlabel = "Time (sec)", ylabel = "Fluorescence Intensity",
                      limits = (-0.05025, 0.60046, 6.9982, 8.6245), xticks = 0:0.1:0.6, yticks = 7.0:0.2:8.6))
    end
    masada_figures = filter(in(experiments), ["Masada_etal_2012_Fig5A", "Masada_etal_2012_Fig5B", "Masada_etal_2012_Fig5C"])
    if !isempty(masada_figures)
        push!(titles, "Masada et al. 2012")
        column = fig[1:2, length(titles)] = GridLayout()
        for (k, experiment) in enumerate(masada_figures)
            if experiment == "Masada_etal_2012_Fig5C"
                trace_panel!(column[k, 1], experiment; xlabel = "Time (sec)", ylabel = "Relative Quin-2 fluorescence",
                             limits = (0.0, 0.04036, 0.98, 1.00008), xticks = 0:0.02:0.04, yticks = 0.98:0.01:1.0)
            else
                trace_panel!(column[k, 1], experiment; xlabel = "Time (sec)", ylabel = "Relative Quin-2 fluorescence",
                             limits = (0.0, 1.00423, 0.937, experiment == "Masada_etal_2012_Fig5A" ? 1.00363 : 1.00334),
                             xticks = 0:0.2:1.0, yticks = 0.94:0.02:1.0)
            end
        end
    end
    for (i, title) in enumerate(titles)
        Label(fig[0, i], title; tellwidth = false)
    end
    length(titles) > 1 && colsize!(fig.layout, 1, Relative(0.3))
    return fig
end

const SETUP_PARAMETER_LABELS = (kon1_N = "z_kon1_N (log₁₀)", kon2_N = "z_kon2_N (log₁₀)", Kd1_N = "z_Kd1_N (log₁₀)",
                                Kd_app_N = "z_Kd_app_N (log₁₀)", kon1_C = "z_kon1_C (log₁₀)", kon2_C = "z_kon2_C (log₁₀)",
                                Kd1_C = "z_Kd1_C (log₁₀)", Kd_app_C = "z_Kd_app_C (log₁₀)")

function setup_random_effects(res)
    re = get_random_effects(res; flatten = false)
    tables = [DataFrame(:setup => re[Symbol("z_", p)].setup, p => re[Symbol("z_", p)].value) for p in keys(SETUP_PARAMETER_LABELS)]
    return reduce((a, b) -> innerjoin(a, b; on = :setup), tables)
end

function padded_limits(values; pad = 0.15)
    lo, hi = extrema(values)
    margin = max(pad * (hi - lo), 0.1)
    return (; low = lo - margin, high = hi + margin)
end

function plot_setup_parameters(df::DataFrame; bins = 6, common_limits = false)
    ps = collect(keys(SETUP_PARAMETER_LABELS))
    shared = padded_limits([vec(Matrix(df[!, ps])); 0.0])
    limits = NamedTuple{Tuple(ps)}(Tuple((; lims = common_limits ? shared : padded_limits([df[!, p]; 0.0])) for p in ps))
    colors = [Makie.to_colormap(:tab10)[mod1(i, 10)] for i in eachindex(df.setup)]
    marker_style = (; markersize = 14, strokecolor = :black, strokewidth = 1)
    setups = [PairPlots.Series(df[df.setup .== s, ps]; color = colors[i]) => (PairPlots.Scatter(; marker_style...),)
              for (i, s) in enumerate(df.setup)]
    population = NamedTuple{Tuple(ps)}(Tuple(0.0 for _ in ps))
    fig = Figure(size = (2940, 1440), fontsize = 18)
    grid = pairplot(fig[1, 1], PairPlots.Series(df[!, ps]; color = (:grey, 0.6), bins = Dict(p => bins for p in ps)) =>
                               (PairPlots.MarginHist(),),
                    setups..., PairPlots.Truth(population; color = :black);
                    labels = Dict(p => string(SETUP_PARAMETER_LABELS[p]) for p in ps), axis = limits)
    entries = [[Makie.MarkerElement(; marker = :circle, color = c, marker_style...) for c in colors];
               Makie.LineElement(; color = :black)]
    Makie.Legend(grid[begin, end + 1], entries, [df.setup; "population (z = 0)"]; tellheight = false, valign = :top,
                 halign = :right, labelsize = 28, patchsize = (40, 30), rowgap = 8)
    return fig
end

plot_setup_parameters(res; kwargs...) = plot_setup_parameters(setup_random_effects(res); kwargs...)
