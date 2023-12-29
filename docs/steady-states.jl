#===
# Steady-state solutions

across a range of glucose levels, from 3mM to 30 mM.
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using PythonCall
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 12
plt.matplotlib.rcParams["font.family"] = "sans-serif"

# Default model
@named sys = make_model()
prob = SteadyStateProblem(sys, []) ## Use default u0
alg = DynamicSS(TRBDF2())
## Warm up
sol = solve(prob, alg)

# Galactose model: glycolysis produces zero net ATP
# By increasing the ATP consumed in the first part of glycolysis from 2 to 4
@named sys_gal = make_model(gk_atp_stoich=4)
prob_gal = SteadyStateProblem(sys_gal, [])

# FFA model: Additional flux reducing the mitochondrial NAD/NADH couple
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
# increase proton leak by 5X
@unpack pHleak = sys
pHleak_val = sys.defaults[pHleak]
prob_hl = SteadyStateProblem(sys, [], [pHleak => pHleak_val * 5])

# Simulating on a range of glucose
# Test on a range of glucose (3 mM to 30 mM)
glc = range(3.0, 30.0, length=101)
prob_func = function (prob, i, repeat)
    sys = prob.f.sys
    idxGlc = indexof(sys.GlcConst, parameters(sys))
    prob.p[idxGlc] = glc[i]
    return prob
end

trajectories = length(glc)
alg = DynamicSS(TRBDF2())

# Run the simulations
sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories)
sim_gal = solve(EnsembleProblem(prob_gal; prob_func), alg; trajectories)
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; trajectories)
sim_f1 = solve(EnsembleProblem(prob_f1; prob_func), alg; trajectories)
sim_etc = solve(EnsembleProblem(prob_etc; prob_func), alg; trajectories)
sim_hl = solve(EnsembleProblem(prob_hl; prob_func), alg; trajectories)

# Plot results
function plot_comparisions(k; figsize=(6, 6), title="", ylabel="", legend_loc="best")
    fig, ax = plt.subplots(; figsize)
    ax.plot(glc, extract(sim, k), color="black", label="Default")
    ax.plot(glc, extract(sim_gal,k), label="Galactose")
    ax.plot(glc, extract(sim_ffa, k), label="FFA")
    ax.plot(glc, extract(sim_f1, k), label="Oligo")
    ax.plot(glc[5:end], extract(sim_hl, k)[5:end], label="FCCP")
    ax.plot(glc, extract(sim_etc, k), label="Rote")
    ax.set_title(title)
    ax.set_xlabel("Glucose (mM)")
    ax.set_ylabel(ylabel)
    ax.legend(loc=legend_loc)
    ax.grid()
    fig.tight_layout()
    return fig
end

@unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, x, degavg, J_O2 = sys

# Mitochondrial membrane potential
fig = plot_comparisions(ΔΨm * 1000, title="(A) Mitochondrial membrane potential", ylabel="ΔΨm (mV)")

# ATP:ADP ratio
fig = plot_comparisions(ATP_c / ADP_c, title="(B) ATP-to-ADP ratio", ylabel="ATP:ADP")

# Mitochondrial NADH:NAD ratio
fig = plot_comparisions(NADH_m / NAD_m, title="(C) Mito. NADH-to-NAD ratio", ylabel="NADH:NAD (mito)", legend_loc="upper left")

# Cytosolic NADH:NAD ratio
fig = plot_comparisions(NADH_c / NAD_c, title="(D) Cyto. NADH-to-NAD ratio", ylabel="NADH:NAD (cyto)", legend_loc="upper left")

# Mitochondrial calcium
fig = plot_comparisions(Ca_m * 1000, title="(E) Mito. calcium", ylabel="[Ca]m (μM)")

# Oxygen consumption rate
fig = plot_comparisions(J_O2, title="(F) Oxygen consumption", ylabel="VO2 (mM/s)")

# Average node degree
fig = plot_comparisions(degavg, title="(G) Average node degree")

# Node ratio
fig = plot_comparisions(x[3] / x[1], title="(H) Degree-3 to drgree-1 ratio")
