#===
# Figure 2

Steady-state solutions across a range of glucose levels.
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
import MitochondrialDynamics: second, μM, mV, mM, Hz
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14
## plt.matplotlib.rcParams["font.sans-serif"] = "Arial"
## plt.matplotlib.rcParams["font.family"] = "sans-serif"

# Default model
@named sys = make_model()
prob = SteadyStateProblem(sys, []) ## Use default u0
sol = solve(prob)

# Galactose model: glycolysis produces zero net ATP
@named sys_gal = make_model(gk_atp_stoich=4)

# FFA model: Additional flux producing mitochondrial NADH
@unpack J_CAC = sys
@named sys_ffa = make_model(j_ffa=sol[J_CAC] * 0.5)

# Simulating on a range of glucose
@unpack GlcConst = sys
idxGlc = findfirst(isequal(GlcConst), parameters(sys))

#---

prob = SteadyStateProblem(sys, [])
prob_gal = SteadyStateProblem(sys_gal, [])
prob_ffa = SteadyStateProblem(sys_ffa, [])
glc = range(3.0mM, 30.0mM, length=51)  # Range of glucose

function prob_func_glc(prob, i, repeat)
    p = copy(prob.p)
    p[idxGlc] = glc[i]
    remake(prob; p=p)
end

alg = DynamicSS(Rodas5())
prob_func=prob_func_glc
trajectories=length(glc)
sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories)
sim_gal = solve(EnsembleProblem(prob_gal; prob_func), alg; trajectories)
sim_ffa = solve(EnsembleProblem(prob_ffa; prob_func), alg; trajectories)

#---

extract(sols, k, scale=1) = map(s->s[k] * scale, sols)

function plot_steady_state(glc, sols, sys; figsize=(10, 10), title="")
    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, x, degavg = sys
    glc5 = glc ./ 5
    g3p = extract(sols, G3P, 1000)
    pyr = extract(sols, Pyr, 1000)
    ca_c = extract(sols, Ca_c, 1000)
    ca_m = extract(sols, Ca_m, 1000)
    nad_ratio_c = extract(sols, NADH_c/NAD_c)
    nad_ratio_m = extract(sols, NADH_m/NAD_m)
    atp_c = extract(sols, ATP_c, 1000)
    adp_c = extract(sols, ADP_c, 1000)
    amp_c = extract(sols, AMP_c, 1000)
    td = extract(sols, ATP_c / ADP_c)
    dpsi = extract(sols, ΔΨm, 1000)
    x1 = extract(sols, x[1])
    x2 = extract(sols, x[2])
    x3 = extract(sols, x[3])
    deg = extract(sols, degavg)

    fig, ax = plt.subplots(3, 3; figsize)

    ax[0, 0].plot(glc5, g3p)
    ax[0, 0].set(title="(A) G3P (μM)", ylim=(0.0, 10.0))
    ax[0, 1].plot(glc5, pyr)
    ax[0, 1].set(title="(B) Pyruvate (μM)")
    ax[0, 2].plot(glc5, ca_c, label="cyto")
    ax[0, 2].plot(glc5, ca_m, label="mito")
    ax[0, 2].legend()
    ax[0, 2].set(title="(C) Calcium (μM)", ylim=(0.0, 1.5))
    ax[1, 0].plot(glc5, nad_ratio_c, label="cyto")
    ax[1, 0].plot(glc5, nad_ratio_m, label="mito")
    ax[1, 0].legend()
    ax[1, 0].set(title="(D) NADH:NAD")
    ax[1, 1].plot(glc5, atp_c, label="ATP")
    ax[1, 1].plot(glc5, adp_c, label="ADP")
    ax[1, 1].plot(glc5, amp_c, label="AMP")
    ax[1, 1].legend()
    ax[1, 1].set(title="(E) Adenylates (μM)")
    ax[1, 2].plot(glc5, td)
    ax[1, 2].set(title="(F) ATP:ADP")
    ax[2, 0].plot(glc5, dpsi, label="cyto")
    ax[2, 0].set(title="(G) ΔΨ (mV)", xlabel="Glucose (X)")
    ax[2, 1].plot(glc5, x1, label="X1")
    ax[2, 1].plot(glc5, x2, label="X2")
    ax[2, 1].plot(glc5, x3, label="X3")
    ax[2, 1].set(title="(H) Mitochondrial nodes", xlabel="Glucose (X)")
    ax[2, 2].plot(glc5, deg)
    ax[2, 2].set(title="(I) Average Node Degree", xlabel="Glucose (X)")

    for i in 0:2, j in 0:2
        ax[i, j].set_xticks(1:6)
        ax[i, j].grid()
    end

    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

#---
fig_glc_default = plot_steady_state(glc, sim, sys, title="Default parameters")
#---
fig_glc_ffa = plot_steady_state(glc, sim_ffa, sys_ffa, title="FFA parameters")
#---
fig_glc_gal = plot_steady_state(glc, sim_gal, sys_gal, title="Galactose parameters")

# Compare default, FFA, and Gal parameter sets

fig, ax = plt.subplots()
lineObjects = ax.plot(1:10, rand(10, 3))

ax.legend(lineObjects, ["foo", "bar", "baz"])

fig

function plot_fig2(glc, sim, sim_gal, sim_ffa, sys; figsize=(8, 8), title="", labels=["Default", "Gal", "FFA"])

    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, degavg = sys

    glc5 = glc ./ 5
    fig, ax = plt.subplots(2, 2; figsize)

    ax[0, 0].set(title="(A) Mitochondrial NADH:NAD")
    k = NADH_m/NAD_m
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = ax[0, 0].plot(glc5, yy)
    ax[0, 0].legend(lines, labels)

    ax[0, 1].set(title="(B) ATP:ADP")
    k = ATP_c/ADP_c
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = ax[0, 1].plot(glc5, yy)
    ax[0, 1].legend(lines, labels)

    ax[1, 0].set(title="(C) ΔΨm (mV)")
    k = ΔΨm
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)] .* 1000
    lines = ax[1, 0].plot(glc5, yy)
    ax[1, 0].legend(lines, labels)

    ax[1, 1].set(title="(D) Average node degree")
    k = degavg
    yy = [extract(sim, k) extract(sim_gal, k) extract(sim_ffa, k)]
    lines = ax[1, 1].plot(glc5, yy)
    ax[1, 1].legend(lines, labels)

    for i in 0:1, j in 0:1
        ax[i, j].set_xticks(1:6)
        ax[i, j].grid()
    end

    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

fig2 = plot_fig2(glc, sim, sim_gal, sim_ffa, sys)

# Tiff figure
fig2.savefig("Fig2.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
