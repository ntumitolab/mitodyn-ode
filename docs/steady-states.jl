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

# Default model
@named sys = make_model()
prob = SteadyStateProblem(sys, []) ## Use default u0
alg = DynamicSS(Rodas5())
## Warm up
sol = solve(prob, alg)

# Galactose model: glycolysis produces zero net ATP
# By increasing the ATP consumed in the first part of glycolysis from 2 to 4
prob_gal = SteadyStateProblem(sys, [], [sys.ATPstiochGK => 4])

# NAD reducing model: an additional reaction to reduce mitochondrial NAD
prob_ffa = SteadyStateProblem(sys, [], [sys.kFFA => sol[0.10 * sys.J_DH / sys.NAD_m]])

# ATP synthase-limiting model (Oligomycin)
prob_f1 = SteadyStateProblem(sys, [], [sys.VmaxF1 => 0.1 * sys.defaults[sys.VmaxF1]])

# ETC-limiting model (Rotenone ?)
prob_etc = SteadyStateProblem(sys, [], [sys.VmaxETC => 0.1 * sys.defaults[sys.VmaxETC]])

# Proton leak-increasing (FCCP) model
# increase proton leak by 7X
prob_hl = SteadyStateProblem(sys, [], [sys.pHleak => 10 * sys.defaults[sys.pHleak]])

# Simulating on a range of glucose
# Test on a range of glucose (3 mM to 30 mM)
glc = range(3.0, 30.0, step=0.3)
idxGlc = indexof(sys.GlcConst, parameters(sys))

prob_func = function (prob, i, repeat)
    prob.p[idxGlc] = glc[i]
    return prob
end

trajectories = length(glc)
alg = DynamicSS(Rodas5())

# Run the simulations
sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories)
sim_gal = solve(EnsembleProblem(prob_gal; prob_func), alg; trajectories)
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; trajectories)
sim_f1 = solve(EnsembleProblem(prob_f1; prob_func), alg; trajectories)
sim_etc = solve(EnsembleProblem(prob_etc; prob_func), alg; trajectories)
sim_hl = solve(EnsembleProblem(prob_hl; prob_func), alg; trajectories);

# Plot results
function plot_comparisions(k; figsize=(6, 6), title="", ylabel="", legend_loc="best")
    fig, ax = plt.subplots(; figsize)
    ax.plot(glc, extract(sim, k), color="black", label="Default")
    ax.plot(glc, extract(sim_gal,k), label="Galactose")
    ax.plot(glc, extract(sim_ffa, k), label="NAD reducing")
    ax.plot(glc, extract(sim_f1, k), label="F1 block")
    ax.plot(glc[9:end], extract(sim_hl, k)[9:end], label="H leak")
    ax.plot(glc, extract(sim_etc, k), label="ETC block")
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
