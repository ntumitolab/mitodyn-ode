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
prob = ODEProblem(sys, [], Inf)
alg = Rodas5()

sol = solve(prob, alg, save_everystep=false, callback=TerminateSteadyState())

# Galactose model: glycolysis produces zero net ATP
# By increasing the ATP consumed in the first part of glycolysis from 2 to 4
prob_gal = ODEProblem(sys, [], Inf, [sys.ATPstiochGK => 4])

# NAD reducing model: an additional reaction to reduce mitochondrial NAD
prob_ffa = ODEProblem(sys, [], Inf, [sys.kFFA => sol[0.10 * sys.J_DH / sys.NAD_m][end]])

# ATP synthase-limiting model (Oligomycin)
prob_f1 = ODEProblem(sys, [], Inf, [sys.VmaxF1 => 0.1 * sys.defaults[sys.VmaxF1]])

# ETC-limiting model (Rotenone ?)
prob_etc = ODEProblem(sys, [], Inf, [sys.VmaxETC => 0.1 * sys.defaults[sys.VmaxETC]])

# Proton leak-increasing (FCCP) model
# increase proton leak by 7X
prob_hl = ODEProblem(sys, [], Inf, [sys.pHleak => 10 * sys.defaults[sys.pHleak]])

# Simulating on a range of glucose
# Test on a range of glucose (3 mM to 30 mM)
glc = range(3.0, 30.0, step=0.3)
idxGlc = indexof(sys.GlcConst, parameters(sys))

prob_func = function (prob, i, repeat)
    prob.p[idxGlc] = glc[i]
    return prob
end

trajectories = length(glc)
alg = Rodas5()
callback=TerminateSteadyState()

# Run the simulations
sim = solve(EnsembleProblem(prob; prob_func), alg; save_everystep=false, trajectories, callback)
sim_gal = solve(EnsembleProblem(prob_gal; prob_func), alg; save_everystep=false, trajectories, callback)
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; save_everystep=false, trajectories, callback)
sim_f1 = solve(EnsembleProblem(prob_f1; prob_func), alg; save_everystep=false, trajectories, callback)
sim_etc = solve(EnsembleProblem(prob_etc; prob_func), alg; save_everystep=false, trajectories, callback)
sim_hl = solve(EnsembleProblem(prob_hl; prob_func), alg; save_everystep=false, trajectories, callback);

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
fig = plot_comparisions(ΔΨm * 1000, title="(A) Mitochondrial membrane potential", ylabel="ΔΨm (mV)");
fig |> PNG

# ATP:ADP ratio
fig = plot_comparisions(ATP_c / ADP_c, title="(B) ATP-to-ADP ratio", ylabel="ATP:ADP");
fig |> PNG

# Mitochondrial NADH:NAD ratio
fig = plot_comparisions(NADH_m / NAD_m, title="(C) Mito. NADH-to-NAD ratio", ylabel="NADH:NAD (mito)", legend_loc="upper left");
fig |> PNG

# Cytosolic NADH:NAD ratio
fig = plot_comparisions(NADH_c / NAD_c, title="(D) Cyto. NADH-to-NAD ratio", ylabel="NADH:NAD (cyto)", legend_loc="upper left");
fig |> PNG

# Mitochondrial calcium
fig = plot_comparisions(Ca_m * 1000, title="(E) Mito. calcium", ylabel="[Ca]m (μM)");
fig |> PNG

# Oxygen consumption rate
fig = plot_comparisions(J_O2, title="(F) Oxygen consumption", ylabel="VO2 (mM/s)");
fig |> PNG

# Average node degree
fig = plot_comparisions(degavg, title="(G) Average node degree");
fig |> PNG

# Node ratio
fig = plot_comparisions(x[3] / x[1], title="(H) Degree-3 to drgree-1 ratio");
fig |> PNG
