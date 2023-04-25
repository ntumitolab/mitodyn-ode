# # Supplimentary figures
# ## FigS1
# Response to elevated glucose concentrations in steps

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

@variables t, Glc(t)
glc_step(t) = 5.0mM * (1 + (t >= 20minute) + (t >= 40minute))
@named sys = make_model(; glceq=Glc~glc_step(t))

ts = range(0, 60minute; step=0.5minute)
prob = ODEProblem(sys, [], ts[end])
sol = solve(prob, tstops=[20minute, 40minute], saveat=ts);

#---

function plot_figs1(
    sol;
    figsize=(10, 10),
    tight=true,
    grid=true
)
    @unpack G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ATP_c, ADP_c, AMP_c, ΔΨm, degavg, x = sol.prob.f.sys
    ts = sol.t
    tsm = ts ./ 60
    g3p = sol[G3P * 1000]
    pyr = sol[Pyr * 1000]
    nadh_c = sol[NADH_c * 1000]
    nadh_m = sol[NADH_m * 1000]
    ca_c = sol[Ca_c * 1000]
    ca_m = sol[Ca_m * 1000]
    atp_c = sol[ATP_c * 1000]
    adp_c = sol[ADP_c * 1000]
    amp_c = sol[AMP_c * 1000]
    dpsi = sol[ΔΨm * 1000]
    k = sol[degavg]
    x1 = sol[x[1]]
    x2 = sol[x[2]]
    x3 = sol[x[3]]

    fig, ax = plt.subplots(3, 3; figsize)

    ax[1, 1].plot(tsm, g3p)
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0.0, 8.0))

    ax[1, 2].plot(tsm, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0.0, 100.0))

    ax[1, 3].plot(tsm, ca_c, label="cyto")
    ax[1, 3].plot(tsm, ca_m, label="mito")
    ax[1, 3].set(title="(C) Calcium (μM)")
    ax[1, 3].legend()

    ax[2, 1].plot(tsm, nadh_c, label="cyto")
    ax[2, 1].plot(tsm, nadh_m, label="mito")
    ax[2, 1].set(title="(D) NADH (μM)")
    ax[2, 1].legend()

    ax[2, 2].plot(tsm, atp_c, label="ATP")
    ax[2, 2].plot(tsm, adp_c, label="ADP")
    ax[2, 2].plot(tsm, amp_c, label="AMP")
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 2].legend()

    ax[2, 3].plot(tsm, atp_c ./ adp_c)
    ax[2, 3].set(title="(F) ATP:ADP", ylim=(0, 45))

    ax[3, 1].plot(tsm, dpsi)
    ax[3, 1].set(title="(G) ΔΨ (mV)", ylim=(80, 150), xlabel="Time (minutes)")

    ax[3, 2].plot(tsm, x1, label="X1")
    ax[3, 2].plot(tsm, x2, label="X2")
    ax[3, 2].plot(tsm, x3, label="X3")
    ax[3, 2].set(title="(H) Mitochondrial nodes", xlabel="Time (minutes)")
    ax[3, 2].legend(loc="upper left")

    ax[3, 3].plot(tsm, k)
    ax[3, 3].set(title="(I) Average Node Degree", xlabel="Time (minutes)")


    for a in ax
        a.grid(grid)
    end

    fig.set_tight_layout(tight)
    return fig
end

#---

figs1 = plot_figs1(sol)
figs1

# TIF file
# `figs1.savefig("FigS1.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

#===
## Figure S2

Response to adding glucose in healthy and DM cells.

In the diabetic parameter set. We adjusted

- PDH capacity
- ETC capacity
- F1 synthase capacity
- Proton leak rate
===#

@unpack VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxVmaxPDH = findfirst(isequal(VmaxPDH), parameters(sys))
idxpHleak = findfirst(isequal(pHleak), parameters(sys))
idxVmaxF1 =  findfirst(isequal(VmaxF1), parameters(sys))
idxVmaxETC =  findfirst(isequal(VmaxETC), parameters(sys))

function remake_dm(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxETC] *= rETC
    p[idxVmaxF1] *= rF1
    p[idxpHleak] *= rHL
    p[idxVmaxPDH] *= rPDH
    return remake(prob, p=p)
end

prob_dm = remake_dm(prob)

#---

solDM = solve(prob_dm, tstops=[20minute, 40minute], saveat=ts);

#---
function plot_figs2(
    sol, solDM;
    figsize=(12, 12),
    labels=["Baseline", "Diabetic"],
    tight=true,
    grid=true
)
    @unpack G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ATP_c, ADP_c, ΔΨm, degavg = sol.prob.f.sys
    ts = sol.t
    tsm = ts ./ 60

    g3p = sol[G3P * 1000]
    pyr = sol[Pyr * 1000]
    nadh_c = sol[NADH_c * 1000]
    nadh_m = sol[NADH_m * 1000]
    ca_c = sol[Ca_c * 1000]
    ca_m = sol[Ca_m * 1000]
    td = sol[ATP_c/ADP_c]
    dpsi = sol[ΔΨm * 1000]
    k = sol[degavg]

    g3pDM = solDM[G3P * 1000]
    pyrDM = solDM[Pyr * 1000]
    nadh_cDM = solDM[NADH_c * 1000]
    nadh_mDM = solDM[NADH_m * 1000]
    ca_cDM = solDM[Ca_c * 1000]
    ca_mDM = solDM[Ca_m * 1000]
    tdDM = solDM[ATP_c/ADP_c]
    dpsiDM = solDM[ΔΨm * 1000]
    kDM = solDM[degavg]

    fig, ax = plt.subplots(3, 3; figsize)

    ax[1, 1].plot(tsm, g3p, label=labels[1])
    ax[1, 1].plot(tsm, g3pDM, label=labels[2])
    ax[1, 1].set(ylabel="G3P (μM)")
    ax[1, 1].set_title("(A)", loc="left")

    ax[1, 2].plot(tsm, pyr, label=labels[1])
    ax[1, 2].plot(tsm, pyrDM, label=labels[2])
    ax[1, 2].set(ylabel="Pyruvate (μM)")
    ax[1, 2].set_title("(B)", loc="left")

    ax[1, 3].plot(tsm, nadh_c, label=labels[1])
    ax[1, 3].plot(tsm, nadh_cDM, label=labels[2])
    ax[1, 3].set(ylabel="Cytosolic NADH (μM)")
    ax[1, 3].set_title("(C)", loc="left")

    ax[2, 1].plot(tsm, nadh_m, label=labels[1])
    ax[2, 1].plot(tsm, nadh_mDM, label=labels[2])
    ax[2, 1].set(ylabel="Mitochondrial NADH (μM)")
    ax[2, 1].set_title("(D)", loc="left")

    ax[2, 2].plot(tsm, ca_c, label=labels[1])
    ax[2, 2].plot(tsm, ca_cDM, label=labels[2])
    ax[2, 2].set(ylabel="Cytosolic Calcium (μM)")
    ax[2, 2].set_title("(E)", loc="left")

    ax[2, 3].plot(tsm, ca_m, label=labels[1])
    ax[2, 3].plot(tsm, ca_mDM, label=labels[2])
    ax[2, 3].set(ylabel="Mitochondrial Calcium (μM)")
    ax[2, 3].set_title("(F)", loc="left")

    ax[3, 1].plot(tsm, td, label=labels[1])
    ax[3, 1].plot(tsm, tdDM, label=labels[2])
    ax[3, 1].set(ylabel="ATP:ADP")
    ax[3, 1].set_title("(G)", loc="left")

    ax[3, 2].plot(tsm, dpsi, label=labels[1])
    ax[3, 2].plot(tsm, dpsiDM, label=labels[2])
    ax[3, 2].set(ylabel="ΔΨm (mV)")
    ax[3, 2].set_title("(H)", loc="left")

    ax[3, 3].plot(tsm, k, label=labels[1])
    ax[3, 3].plot(tsm, kDM, label=labels[2])
    ax[3, 3].set(ylabel="Average Node Degree")
    ax[3, 3].set_title("(I)", loc="left")

    for a in ax
        a.grid(grid)
        a.legend()
    end

    fig.set_tight_layout(tight)
    return fig
end

figs2 = plot_figs2(sol, solDM)
figs2

#===
## Figure S3

Changes in response to both glucose stimulation and chemical agents.

Using the Glucose-Oligomycin-FCCP protocol.
===#

@named sys = make_model()

@unpack GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxGlc = findfirst(isequal(GlcConst), parameters(sys))
idxVmaxPDH = findfirst(isequal(VmaxPDH), parameters(sys))
idxpHleak = findfirst(isequal(pHleak), parameters(sys))
idxVmaxF1 =  findfirst(isequal(VmaxF1), parameters(sys))
idxVmaxETC =  findfirst(isequal(VmaxETC), parameters(sys))

tend = 80minute
ts = range(0, tend, 401)

prob = ODEProblem(sys, [], ts[end])

function remake_dm(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxETC] *= rETC
    p[idxVmaxF1] *= rF1
    p[idxpHleak] *= rHL
    p[idxVmaxPDH] *= rPDH
    return remake(prob, p=p)
end

prob_dm = remake_dm(prob)

# Define events
function add_glucose!(i)
    i.p[idxGlc] += 20mM
    set_proposed_dt!(i, 0.1)
end

add_glucose_cb = PresetTimeCallback(20minute, add_glucose!)

function add_oligomycin!(i)
    i.p[idxVmaxF1] *= 0.05
    set_proposed_dt!(i, 0.1)
end

add_oligomycin_cb = PresetTimeCallback(40minute, add_oligomycin!)

function add_rotenone!(i)
    i.p[idxVmaxETC] *= 0.05
    set_proposed_dt!(i, 0.1)
end

add_rotenone_cb = PresetTimeCallback(60minute, add_rotenone!)

function add_fccp!(i)
    i.p[idxpHleak] *= 10
    set_proposed_dt!(i, 0.1)
end

add_fccp_cb = PresetTimeCallback(60minute, add_fccp!)

cbs = CallbackSet(add_glucose_cb, add_oligomycin_cb, add_fccp_cb)
sols3 = solve(prob; callback=cbs, saveat=ts)
solDMs3 = solve(prob_dm; callback=cbs, saveat=ts)

figs3 = plot_figs2(sols3, solDMs3)
figs3

#===
## Figure S4
Oxygen consumption in response to both glucose stimulation and chemical agents.
===#

function plot_jo2(sol, solDM;
    labels=["Baseline", "Diabetic"],
    figsize=(6, 6),
    tight=true,
    grid=true
)
    ts = sol.t
    tsm = ts ./ 60

    @unpack J_O2 = sol.prob.f.sys
    jo2 = sol[J_O2]
    jo2DM = solDM[J_O2]

    fig, ax = plt.subplots(; figsize)

    ax.plot(tsm, jo2, label=labels[1])
    ax.plot(tsm, jo2DM, label=labels[2])
    ax.set(xlabel="Time (minute)", ylabel="Rate (mM/s)", title="Oxygen consumption rate")

    ax.grid(grid)
    ax.legend()
    fig.set_tight_layout(tight)
    return fig
end

#---

figs4 = plot_jo2(sols3, solDMs3)
figs4

# ## Figure S5
# Baseline vs. Diabetic models using Glucose-Oligomycin-Rotenone protocol.
cbs = CallbackSet(add_glucose_cb, add_oligomycin_cb, add_rotenone_cb)
sols5 = solve(prob; callback=cbs, saveat=ts)
sols5DM = solve(prob_dm; callback=cbs, saveat=ts)
figs5 = plot_figs2(sols5, sols5DM)
figs5

# TIF file
# `figs5.savefig("FigS5.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure S6
# Oxygen consumption of Baseline and Diabetic models using the protocol from Figure S5.

figs6 = plot_jo2(sols5, sols5DM)
figs6

# TIF file
# `figs6.savefig("FigS6.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`
