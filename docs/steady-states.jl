#===
# Steady-state solutions

across a range of glucose levels, from 3mM to 30 mM.
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using PythonCall
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14
plt.matplotlib.rcParams["font.family"] = "sans-serif"

# Extract values from ensemble simulations with a symbol
extract(sols, k, scale=1) = getindex.(sols.u, k * scale)

# Default model
@named sys = make_model()
prob = SteadyStateProblem(sys, []) ## Use default u0
alg = DynamicSS(Rodas5())
sol = solve(prob, alg) ## Warm up

# Galactose model: glycolysis produces zero net ATP
# By increasing the ATP consumed in the first part of glycolysis from 2 to 4
@named sys_gal = make_model(gk_atp_stoich=4)
prob_gal = SteadyStateProblem(sys_gal, [])

# FFA model: Additional flux reducing mitochondrial NAD/NADH couple
# A 50% increase w.r.t baseline CAC flux
@named sys_ffa = make_model()
@unpack J_CAC, J_FFA = sys_ffa
prob_ffa = SteadyStateProblem(sys_ffa, [], [J_FFA => sol[J_CAC] * 0.5])

# ATP ssynthase-limiting (Oligomycin) model
@unpack VmaxF1 = sys
VmaxF1_val = sys.defaults[VmaxF1]
prob_f1 = SteadyStateProblem(sys, [], [VmaxF1 => VmaxF1_val * 0.1])

# ETC-limiting (Rotenone) model
@unpack VmaxETC = sys
VmaxETC_val = sys.defaults[VmaxETC]
prob_etc = SteadyStateProblem(sys, [], [VmaxETC => VmaxETC_val * 0.1])

# Proton leak-increasing (FCCP) model
@unpack pHleak = sys
pHleak_val = sys.defaults[pHleak]
prob_hl = SteadyStateProblem(sys, [], [pHleak => pHleak_val * 5])

# Simulating on a range of glucose
# Test on a range of glucose (3 mM to 30 mM)
idxGlc = indexof(sys.GlcConst, parameters(sys))
glc = range(3.0, 30.0, length=101)

prob_func = function (prob, i, repeat)
    prob.p[idxGlc] = glc[i]
    return prob
end

alg = DynamicSS(Rodas5())
trajectories = length(glc)

# Run the simulations
sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories)
sim_gal = solve(EnsembleProblem(prob_gal; prob_func), alg; trajectories)
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; trajectories)
sim_f1 = solve(EnsembleProblem(prob_f1; prob_func), alg; trajectories)
sim_etc = solve(EnsembleProblem(prob_etc; prob_func), alg; trajectories)
sim_hl = solve(EnsembleProblem(prob_hl; prob_func), alg; trajectories)

@unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, x, degavg = sys

## Plot results
function plot_comparisions(k)
    plt.figure()
    plt.plot(glc, extract(sim, k), color="black", label="Default")
    plt.plot(glc, extract(sim_gal,k), label="Galactose")
    plt.plot(glc, extract(sim_ffa, k), label="FFA")
    plt.plot(glc, extract(sim_f1, k), label="Oligo")
    plt.plot(glc[5:end], extract(sim_hl, k)[5:end], label="FCCP")
    plt.plot(glc, extract(sim_etc, k), label="Rote")
    plt.legend()
    plt.gcf()
end

# Mitochondrial membrane potential
plot_comparisions(ΔΨm * 1000)

# ATP:ADP ratio
plot_comparisions(ATP_c / ADP_c)

# Terminal node / branching node ratio
plot_comparisions(x[1] / x[3])

# Mitochondrial NADH:NAD ratio
plot_comparisions(NADH_m / NAD_m)

# Cytosolic NADH:NAD ratio
plot_comparisions(NADH_c / NAD_c)

# Average node degree
plot_comparisions(degavg)

# Mitochondrial calcium
plot_comparisions(Ca_m)
