#===
# Calcium overload

Steady-state solutions across a range of glucose levels.
===#
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using ModelingToolkit
using Startup
using Startup: μM
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

# `sol[k]` is a time series for an ODE solution; steady states are its last value.
extract_ss(sim, k) = map(s -> last(s[k]), sim)

# Default model
@time "Build system" @named sys = make_model()
# Integrate the physical ODE to steady state. `SteadyStateProblem` is lowered
# to an SCC nonlinear problem whose block-wise pseudo-transient flow is
# unstable for the fission/fusion block and gives nonphysical concentrations.
@time "Build problem" prob = ODEProblem(sys, [], (0.0, 1e5))
alg = TRBDF2()
@time "Solve problem" sol = solve(prob, alg)

# High calcium model
@unpack RestingCa, ActivatedCa = sys
prob_ca5 = ODEProblem(sys, [RestingCa=>0.45μM, ActivatedCa=>1.25μM], (0.0, 1e5))
prob_ca10 = ODEProblem(sys, [RestingCa=>0.9μM, ActivatedCa=>2.5μM], (0.0, 1e5))

# Simulating on a range of glucose
@unpack Glc = sys

# Test on a range of glucose
glc = 3.5:0.5:30.0
prob_func = (prob, ctx) -> begin
    remake(prob, p=[Glc => glc[ctx.sim_id]])
end

trajectories=length(glc)

@time sim = map(glc) do g
    _prob = remake(prob, p=[Glc => g])
    solve(_prob, alg)
end;

@time sim_ca5 = map(glc) do g
    _prob = remake(prob_ca5, p=[Glc => g])
    solve(_prob, alg)
end;

@time sim_ca10 = map(glc) do g
    _prob = remake(prob_ca10, p=[Glc => g])
    solve(_prob, alg)
end;

# ## Steady states for a range of glucose

function plot_steady_state(glc, sols, sys; figsize=(10, 10), title="")

    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, x1, x2, x3, degavg = sys

    glc5 = glc ./ 5
    g3p = extract_ss(sols, G3P * 1000)
    pyr = extract_ss(sols, Pyr * 1000)
    ca_c = extract_ss(sols, Ca_c * 1000)
    ca_m = extract_ss(sols, Ca_m * 1000)
    nad_ratio_c = extract_ss(sols, NADH_c/NAD_c)
    nad_ratio_m = extract_ss(sols, NADH_m/NAD_m)
    atp_c = extract_ss(sols, ATP_c * 1000)
    adp_c = extract_ss(sols, ADP_c * 1000)
    amp_c = extract_ss(sols, AMP_c * 1000)
    td = extract_ss(sols, ATP_c / ADP_c)
    dpsi = extract_ss(sols, ΔΨm * 1000)
    x1 = extract_ss(sols, x1)
    x2 = extract_ss(sols, x2)
    x3 = extract_ss(sols, x3)
    deg = extract_ss(sols, degavg)

    numrows = 3
    numcols = 3

    fig, ax = plt.subplots(numrows, numcols; figsize)

    ax[0, 0].plot(glc5, g3p)
    ax[0, 0].set(ylabel="G3P (μM)")
    ax[0, 0].set_title("a", loc="left")
    ax[0, 1].plot(glc5, pyr)
    ax[0, 1].set(ylabel="Pyruvate (μM)")
    ax[0, 1].set_title("b", loc="left")
    ax[0, 2].plot(glc5, ca_c, label="cyto")
    ax[0, 2].plot(glc5, ca_m, label="mito")
    ax[0, 2].legend()
    ax[0, 2].set(ylabel="Calcium (μM)")
    ax[0, 2].set_title("c", loc="left")
    ax[1, 0].plot(glc5, nad_ratio_c, label="cyto")
    ax[1, 0].plot(glc5, nad_ratio_m, label="mito")
    ax[1, 0].legend()
    ax[1, 0].set(ylabel="NADH:NAD (ratio)")
    ax[1, 0].set_title("d", loc="left")
    ax[1, 1].plot(glc5, atp_c, label="ATP")
    ax[1, 1].plot(glc5, adp_c, label="ADP")
    ax[1, 1].plot(glc5, amp_c, label="AMP")
    ax[1, 1].legend()
    ax[1, 1].set(ylabel="Adenylates (μM)")
    ax[1, 1].set_title("e", loc="left")
    ax[1, 2].plot(glc5, td)
    ax[1, 2].set(ylabel="ATP:ADP (ratio)")
    ax[1, 2].set_title("f", loc="left")
    ax[2, 0].plot(glc5, dpsi, label="cyto")
    ax[2, 0].set(ylabel="ΔΨ (mV)", xlabel="Glucose (X)")
    ax[2, 0].set_title("g", loc="left")
    ax[2, 1].plot(glc5, x1, label="X1")
    ax[2, 1].plot(glc5, x2, label="X2")
    ax[2, 1].plot(glc5, x3, label="X3")
    ax[2, 1].set(ylabel="Mitochondrial nodes (a.u.)", xlabel="Glucose (X)")
    ax[2, 1].set_title("h", loc="left")
    ax[2, 2].plot(glc5, deg)
    ax[2, 2].set(ylabel="Avg. Node Degree (a.u.)", xlabel="Glucose (X)")
    ax[2, 2].set_title("i", loc="left")

    for i in 0:numrows-1, j in 0:numcols-1
        ax[i, j].set_xticks(1:6)
        ax[i, j].grid()
    end
    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

# Default model
fig_glc_default = plot_steady_state(glc, sim, sys, title="Calcium 1X")
# High calcium (5X)
fig_ca5 = plot_steady_state(glc, sim_ca5, sys, title="Calcium 5X")
# High calcium (10X)
fig_ca10 = plot_steady_state(glc, sim_ca10, sys, title="Calcium 10X")

# ## Comparing default and high calcium models
function plot_comparision(glc, sim, sim_ca5, sim_ca10, sys;
    figsize=(8, 10), title="", labels=["Ca 1X", "Ca 5X", "Ca 10X"]
)
    @unpack G3P, Pyr, Ca_c, Ca_m, NADH_c, NADH_m, NAD_c, NAD_m, ATP_c, ADP_c, AMP_c, ΔΨm, degavg, J_O2 = sys

    glc5 = glc ./ 5

    numrows = 3
    numcols = 2
    fig, ax = plt.subplots(numrows, numcols; figsize)

    ax[0, 0].set_title("a", loc="left")
    ax[0, 0].set_ylabel("Cyto. NADH:NAD (ratio)")
    k = NADH_c/NAD_c
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[0, 0].plot(glc5, yy)
    ax[0, 0].legend(lines, labels)

    ax[0, 1].set_title("b", loc="left")
    ax[0, 1].set_ylabel("Mito. NADH:NAD (ratio)")
    k = NADH_m/NAD_m
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[0, 1].plot(glc5, yy)
    ax[0, 1].legend(lines, labels)

    ax[1, 0].set_title("c", loc="left")
    ax[1, 0].set_ylabel("ATP:ADP (ratio)")
    k = ATP_c/ADP_c
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[1, 0].plot(glc5, yy)
    ax[1, 0].legend(lines, labels)

    ax[1, 1].set_title("d", loc="left")
    ax[1, 1].set_ylabel("ΔΨm (mV)")
    k = ΔΨm * 1000
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[1, 1].plot(glc5, yy)
    ax[1, 1].legend(lines, labels)

    ax[2, 0].set_title("e", loc="left")
    ax[2, 0].set_ylabel("Avg. node degree (ratio)")
    k = degavg
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[2, 0].plot(glc5, yy)
    ax[2, 0].legend(lines, labels, loc="lower right")
    ax[2, 0].set(xlabel="Glucose (X)")

    ax[2, 1].set_title("f", loc="left")
    ax[2, 1].set_ylabel("VO2 (mM/s)")
    k = J_O2
    yy = [extract_ss(sim, k) extract_ss(sim_ca5, k) extract_ss(sim_ca10, k)]
    lines = ax[2, 1].plot(glc5, yy)
    ax[2, 1].legend(lines, labels)
    ax[2, 1].set(xlabel="Glucose (X)")

    for i in 0:numrows-1, j in 0:numcols-1
        ax[i, j].set_xticks(1:6)
        ax[i, j].grid()
    end

    fig.suptitle(title)
    fig.tight_layout()
    return fig
end

figcomp = plot_comparision(glc, sim, sim_ca5, sim_ca10, sys)

# Export figure
exportTIF(figcomp, "S1_HighCa.tif")

# ## MMP vs <k>
function plot_dpsi_k(sim, sim_ca5, sim_ca10, sys; figsize=(6,6), title="", labels=["Ca 1X", "Ca 5X", "Ca 10X"])
    @unpack ΔΨm, degavg = sys

    fig, ax = plt.subplots(1, 1; figsize)

    ax.plot(extract_ss(sim, ΔΨm * 1000), extract_ss(sim, degavg), "v", label=labels[1])
    ax.plot(extract_ss(sim_ca5, ΔΨm * 1000), extract_ss(sim_ca5, degavg), "o", label=labels[2])
    ax.plot(extract_ss(sim_ca10, ΔΨm * 1000), extract_ss(sim_ca10, degavg), "x", label=labels[3])
    ax.set(xlabel="ΔΨm (mV)", ylabel="Average node degree", title=title)
    ax.legend()
    ax.grid()

    return fig
end

fig = plot_dpsi_k(sim, sim_ca5, sim_ca10, sys)

#---
## exportTIF(fig, "S1_HighCa_dpsi_k.tif")

# ## x-axis as Ca2+ and y-axis as average node degree

function plot_ca_k(sim, sim_ca5, sim_ca10, sys; figsize=(6,6), title="", labels=["Ca 1X", "Ca 5X", "Ca 10X"])
    @unpack Ca_m, degavg = sys

    fig, ax = plt.subplots(1, 1; figsize)

    ax.plot(extract_ss(sim, Ca_m * 1000), extract_ss(sim, degavg), "v", label=labels[1])
    ax.plot(extract_ss(sim_ca5, Ca_m * 1000), extract_ss(sim_ca5, degavg), "o", label=labels[2])
    ax.plot(extract_ss(sim_ca10, Ca_m * 1000), extract_ss(sim_ca10, degavg), "x", label=labels[3])
    ax.set(xlabel="Mitochondrial Ca (mM)", ylabel="Average node degree", title=title)
    ax.legend()
    ax.grid()
    return fig
end

fig = plot_ca_k(sim, sim_ca5, sim_ca10, sys)

#---
exportTIF(fig, "S1_HighCa_ca_k.tif")

# ## x-axis as ATP and y-axis as average node degree
function plot_atp_k(sim, sim_ca5, sim_ca10, sys; figsize=(6,6), title="", labels=["Ca 1X", "Ca 5X", "Ca 10X"])
    @unpack ATP_c, ADP_c, degavg = sys

    k = ATP_c / ADP_c

    fig, ax = plt.subplots(1, 1; figsize)

    ax.plot(extract_ss(sim, k), extract_ss(sim, degavg), "v", label=labels[1])
    ax.plot(extract_ss(sim_ca5, k), extract_ss(sim_ca5, degavg), "o", label=labels[2])
    ax.plot(extract_ss(sim_ca10, k), extract_ss(sim_ca10, degavg), "x", label=labels[3])
    ax.set(xlabel="ATP:ADP ratio", ylabel="Average node degree", title=title)
    ax.legend()
    ax.grid()

    return fig
end

fig = plot_atp_k(sim, sim_ca5, sim_ca10, sys)

#---
exportTIF(fig, "S1_HighCa_atp_k.tif")
