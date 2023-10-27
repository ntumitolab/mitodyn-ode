#===
# Calcium overload

Steady-state solutions across a range of glucose levels.
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

using MitochondrialDynamics: μM

# Default model
@named sys = make_model()
prob = SteadyStateProblem(sys, [], jac=true) ## Use default u0
alg = DynamicSS(Rodas5())
sol = solve(prob, alg)

# High calcium model
@unpack RestingCa, ActivatedCa = sys
prob_ca5 = SteadyStateProblem(sys, [], [RestingCa=>0.45μM, ActivatedCa=>1.25μM], jac=true)
prob_ca10 = SteadyStateProblem(sys, [], [RestingCa=>0.9μM, ActivatedCa=>2.5μM], jac=true)

# Simulating on a range of glucose
@unpack GlcConst = sys
idxGlc = findfirst(isequal(GlcConst), parameters(sys))

# Test on a range of glucose
glc = 3.5:0.5:30.0
prob_func = function (prob, i, repeat)
    prob.p[idxGlc] = glc[i]
    prob
end

alg = DynamicSS(Rodas5())
trajectories=length(glc)

sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories)
sim_ca5 = solve(EnsembleProblem(prob_ca5; prob_func), alg; trajectories)
sim_ca10 = solve(EnsembleProblem(prob_ca10; prob_func), alg; trajectories)

# ## Steady states for a range of glucose

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

    ax[1, 1].plot(glc5, g3p)
    ax[1, 1].set(title="(A) G3P (μM)")
    ax[1, 2].plot(glc5, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)")
    ax[1, 3].plot(glc5, ca_c, label="cyto")
    ax[1, 3].plot(glc5, ca_m, label="mito")
    ax[1, 3].legend()
    ax[1, 3].set(title="(C) Calcium (μM)")
    ax[2, 1].plot(glc5, nad_ratio_c, label="cyto")
    ax[2, 1].plot(glc5, nad_ratio_m, label="mito")
    ax[2, 1].legend()
    ax[2, 1].set(title="(D) NADH:NAD")
    ax[2, 2].plot(glc5, atp_c, label="ATP")
    ax[2, 2].plot(glc5, adp_c, label="ADP")
    ax[2, 2].plot(glc5, amp_c, label="AMP")
    ax[2, 2].legend()
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 3].plot(glc5, td)
    ax[2, 3].set(title="(F) ATP:ADP")
    ax[3, 1].plot(glc5, dpsi, label="cyto")
    ax[3, 1].set(title="(G) ΔΨ (mV)", xlabel="Glucose (X)")
    ax[3, 2].plot(glc5, x1, label="X1")
    ax[3, 2].plot(glc5, x2, label="X2")
    ax[3, 2].plot(glc5, x3, label="X3")
    ax[3, 2].set(title="(H) Mitochondrial nodes", xlabel="Glucose (X)")
    ax[3, 3].plot(glc5, deg)
    ax[3, 3].set(title="(I) Average Node Degree", xlabel="Glucose (X)")

    for a in ax
        a.set_xticks(1:6)
        a.grid()
    end
    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

# Default model
fig_glc_default = plot_steady_state(glc, sim, sys, title="Default parameters")

# High calcium
fig_ca5 = plot_steady_state(glc, sim_ca5, sys, title="Calcium 5X")
fig_ca10 = plot_steady_state(glc, sim_ca10, sys, title="Calcium 10X")

# ## Comparing default and high calcium models

function plot_comparision(glc, sim, sim_ca5, sim_ca10, sys;
    figsize=(8, 8), title="", labels=["Default", "Ca 5X", "Ca 10X"]
)

    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, degavg, J_O2 = sys

    glc5 = glc ./ 5
    fig, ax = plt.subplots(3, 2; figsize)

    ax[1, 1].set(title="(A) Cytosolic NADH:NAD")
    k = NADH_c/NAD_c
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)]
    lines = ax[1, 1].plot(glc5, yy)
    ax[1, 1].legend(lines, labels)

    ax[1, 2].set(title="(B) Mitochondrial NADH:NAD")
    k = NADH_m/NAD_m
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)]
    lines = ax[1, 2].plot(glc5, yy)
    ax[1, 2].legend(lines, labels)

    ax[2, 1].set(title="(C) ATP:ADP")
    k = ATP_c/ADP_c
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)]
    lines = ax[2, 1].plot(glc5, yy)
    ax[2, 1].legend(lines, labels)

    ax[2, 2].set(title="(D) ΔΨm (mV)")
    k = ΔΨm
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)] .* 1000
    lines = ax[2, 2].plot(glc5, yy)
    ax[2, 2].legend(lines, labels)

    ax[3, 1].set(title="(E) Average node degree")
    k = degavg
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)]
    lines = ax[3, 1].plot(glc5, yy)
    ax[3, 1].legend(lines, labels, loc="lower right")
    ax[3, 1].set(xlabel="Glucose (X)")

    ax[3, 2].set(title="(F) Oxygen consumption")
    k = J_O2
    yy = [extract(sim, k) extract(sim_ca5, k) extract(sim_ca10, k)]
    lines = ax[3, 2].plot(glc5, yy)
    ax[3, 2].legend(lines, labels)
    ax[3, 2].set(xlabel="Glucose (X)", ylabel="mM/s")

    for a in ax
        a.set_xticks(1:6)
        a.grid()
    end

    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

figcomp = plot_comparision(glc, sim, sim_ca5, sim_ca10, sys)

# Export figure
## figcomp.savefig("S1_HighCa.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))

# ## mitochondria membrane potential vs average node degree

function plot_dpsi_k(sim, sim_ca5, sim_ca10, sys; figsize=(6,6), title="", labels=["Default", "Ca 5X", "Ca 10X"])
    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)
    @unpack ΔΨm, degavg = sys

    fig, ax = plt.subplots(1, 1; figsize)

    ax.plot(extract(sim, ΔΨm, 1000), extract(sim, degavg), "v", label=labels[1])
    ax.plot(extract(sim_ca5, ΔΨm, 1000), extract(sim_ca5, degavg), "o", label=labels[2])
    ax.plot(extract(sim_ca10, ΔΨm, 1000), extract(sim_ca10, degavg), "x", label=labels[3])
    ax.set(xlabel="ΔΨm (mV)", ylabel="Average node degree", title=title)
    ax.legend()
    ax.grid()

    return fig
end

fig = plot_dpsi_k(sim, sim_ca5, sim_ca10, sys)

# ## x-axis as Ca2+ and y-axis as average node

function plot_ca_k(sim, sim_ca5, sim_ca10, sys; figsize=(6,6), title="", labels=["Default", "Ca 5X", "Ca 10X"])
    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)
    @unpack Ca_m, degavg = sys

    fig, ax = plt.subplots(1, 1; figsize)

    ax.plot(extract(sim, Ca_m, 1000), extract(sim, degavg), "v", label=labels[1])
    ax.plot(extract(sim_ca5, Ca_m, 1000), extract(sim_ca5, degavg), "o", label=labels[2])
    ax.plot(extract(sim_ca10, Ca_m, 1000), extract(sim_ca10, degavg), "x", label=labels[3])
    ax.set(xlabel="Mitochondrial Ca (mM)", ylabel="Average node degree", title=title)
    ax.legend()
    ax.grid()

    return fig
end

fig = plot_ca_k(sim, sim_ca5, sim_ca10, sys)