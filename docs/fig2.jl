#===
# Figure 2

Steady-state solutions across a range of glucose levels.
===#
using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using PythonCall
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

# PNG output in Literate.jl
PNG(fig) = display("image/png", fig)

# Default model
@named sys = make_model()
prob = ODEProblem(sys, [], Inf)
alg = Rodas5()
sol = solve(prob, alg, save_everystep=false, callback=TerminateSteadyState())

# Galactose model: glycolysis produces zero net ATP
# By increasing the ATP consumed in the first part of glycolysis from 2 to 4
prob_gal = ODEProblem(sys, [], Inf, [sys.ATPstiochGK => 4])

# FFA model: Additional flux reducing mitochondrial NAD/NADH couple
# A 50% increase w.r.t baseline CAC flux
prob_ffa = ODEProblem(sys, [], Inf, [sys.kFFA => sol[0.10 * sys.J_DH / sys.NAD_m][end]])

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
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; save_everystep=false, trajectories, callback);

# ## Steady states for a range of glucose
function plot_steady_state(glc, sols, sys; figsize=(10, 10), title="")

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, x, degavg = sys

    glc5 = glc ./ 5
    g3p = extract(sols, G3P * 1000)
    pyr = extract(sols, Pyr * 1000)
    ca_c = extract(sols, Ca_c * 1000)
    ca_m = extract(sols, Ca_m * 1000)
    nad_ratio_c = extract(sols, NADH_c / NAD_c)
    nad_ratio_m = extract(sols, NADH_m / NAD_m)
    atp_c = extract(sols, ATP_c * 1000)
    adp_c = extract(sols, ADP_c * 1000)
    amp_c = extract(sols, AMP_c * 1000)
    td = extract(sols, ATP_c / ADP_c)
    dpsi = extract(sols, ΔΨm * 1000)
    x1 = extract(sols, x[1])
    x2 = extract(sols, x[2])
    x3 = extract(sols, x[3])
    deg = extract(sols, degavg)

    numrows = 3
    numcols = 3
    fig, axs = plt.subplots(numrows, numcols; figsize)

    axs[0, 0].plot(glc5, g3p)
    axs[0, 0].set(title="(A) G3P (μM)", ylim=(0.0, 10.0))
    axs[0, 1].plot(glc5, pyr)
    axs[0, 1].set(title="(B) Pyruvate (μM)")
    axs[0, 2].plot(glc5, ca_c, label="cyto")
    axs[0, 2].plot(glc5, ca_m, label="mito")
    axs[0, 2].legend()
    axs[0, 2].set(title="(C) Calcium (μM)", ylim=(0.0, 1.5))
    axs[1, 0].plot(glc5, nad_ratio_c, label="cyto")
    axs[1, 0].plot(glc5, nad_ratio_m, label="mito")
    axs[1, 0].legend()
    axs[1, 0].set(title="(D) NADH:NAD")
    axs[1, 1].plot(glc5, atp_c, label="ATP")
    axs[1, 1].plot(glc5, adp_c, label="ADP")
    axs[1, 1].plot(glc5, amp_c, label="AMP")
    axs[1, 1].legend()
    axs[1, 1].set(title="(E) Adenylates (μM)")
    axs[1, 2].plot(glc5, td)
    axs[1, 2].set(title="(F) ATP:ADP")
    axs[2, 0].plot(glc5, dpsi, label="cyto")
    axs[2, 0].set(title="(G) ΔΨ (mV)", xlabel="Glucose (X)")
    axs[2, 1].plot(glc5, x1, label="X1")
    axs[2, 1].plot(glc5, x2, label="X2")
    axs[2, 1].plot(glc5, x3, label="X3")
    axs[2, 1].legend()
    axs[2, 1].set(title="(H) Mitochondrial nodes", xlabel="Glucose (X)")
    axs[2, 2].plot(glc5, deg)
    axs[2, 2].set(title="(I) Average Node Degree", xlabel="Glucose (X)")

    for i in 0:numrows-1, j in 0:numcols-1
        axs[i, j].set_xticks(1:6)
        axs[i, j].grid()
    end
    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

#---
fig_glc_default = plot_steady_state(glc, sim, sys, title="")

# Default parameters
fig_glc_default.savefig("Fig2.tif", dpi=300, pil_kwargs=pydict(Dict("compression" => "tiff_lzw")))

# Adding free fatty acids
fig_glc_ffa = plot_steady_state(glc, sim_ffa, sys, title="FFA parameters")

# Using galactose instead of glucose as the hydrocarbon source
fig_glc_gal = plot_steady_state(glc, sim_gal, sys, title="Galactose parameters")

# ## Comparing default, FFA, and galactose models

function plot_ffa_gal(glc, sim, sim_gal, sim_ffa, sys; figsize=(10, 10), title="", labels=["Default", "Gal", "FFA"])

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, degavg, J_O2 = sys

    glc5 = glc ./ 5
    numcol = 2
    numrow = 3
    fig, axs = plt.subplots(numrow, numcol; figsize)

    axs[0, 0].set(title="(A) Cytosolic NADH:NAD")
    k = NADH_c / NAD_c
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[0, 0].plot(glc5, yy)
    axs[0, 0].legend(lines, labels)

    axs[0, 1].set(title="(B) Mitochondrial NADH:NAD")
    k = NADH_m / NAD_m
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[0, 1].plot(glc5, yy)
    axs[0, 1].legend(lines, labels)

    axs[1, 0].set(title="(C) ATP:ADP")
    k = ATP_c / ADP_c
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[1, 0].plot(glc5, yy)
    axs[1, 0].legend(lines, labels)

    axs[1, 1].set(title="(D) ΔΨm (mV)")
    k = ΔΨm * 1000
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[1, 1].plot(glc5, yy)
    axs[1, 1].legend(lines, labels)

    axs[2, 0].set(title="(E) Average node degree")
    k = degavg
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[2, 0].plot(glc5, yy)
    axs[2, 0].legend(lines, labels)
    axs[2, 0].set(xlabel="Glucose (X)")

    axs[2, 1].set(title="(F) Oxygen consumption")
    k = J_O2
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = axs[2, 1].plot(glc5, yy)
    axs[2, 1].legend(lines, labels)
    axs[2, 1].set(xlabel="Glucose (X)", ylabel="mM/s")

    for i in 0:numrow-1, j in 0:numcol-1
        axs[i, j].set_xticks(1:6)
        axs[i, j].grid()
    end
    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

figFFAGal = plot_ffa_gal(glc, sim, sim_gal, sim_ffa, sys)

# Export figure
figFFAGal.savefig("s1-fig1.tif", dpi=300, pil_kwargs=pydict(Dict("compression" => "tiff_lzw")))
