md"""
# Figure S2 : DM step response

Response to adding glucose in healthy and DM cells.

In the diabetic parameter set. We adjusted

- PDH capacity
- ETC capacity
- F1 synthase capacity
- Proton leak rate
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC, J_ANT, J_O2
using MitochondrialDynamics: G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ΔΨm, ATP_c, ADP_c, AMP_c, degavg, t, x
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

glc_step(t) = 5.0mM * (1 + (t >= 20minute) + (t >= 40minute))
@named sys = make_model(; glcrhs=glc_step(t))

ts = range(0, 60minute, 201)
prob = ODEProblem(sys, [], ts[end])
sol = solve(prob, tstops=[20minute, 40minute], saveat=ts);

#---

pidx = Dict(k => i for (i, k) in enumerate(parameters(sys)))

function make_dm_prob(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    idxVmaxPDH = pidx[VmaxPDH]
    idxpHleak = pidx[pHleak]
    idxVmaxF1 = pidx[VmaxF1]
    idxVmaxETC = pidx[VmaxETC]
    p[idxVmaxPDH] *= rPDH
    p[idxpHleak] *= rHL
    p[idxVmaxF1] *= rF1
    p[idxVmaxETC] *= rETC
    return remake(prob, p=p)
end

#---

probDM = make_dm_prob(prob)
solDM = solve(probDM, tstops=[20minute, 40minute], saveat=ts);

#---

function plot_figs2(
    sol, solDM;
    figsize=(12, 12),
    labels=["Baseline", "Diabetic"],
    tight=true
)
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
        a.grid()
        a.legend()
    end

    fig.set_tight_layout(tight)
    return fig
end

#---

figs2 = plot_figs2(sol, solDM)
figs2

# TIF file
# `figs2.savefig("FigS2.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

#---

#===
## Figure S3

Changes in response to both glucose stimulation and chemical agents.

Using the Glucose-Oligomycin-FCCP protocol.
===#

add_glucose(t) = 5mM + 15mM * (t >= 20minute)
add_oligomycin(t) = 1.0 - 0.95 * (t >= 40minute)
add_rotenone(t) = 1.0 - 0.95 * (t >= 60minute)
add_fccp(t) = 1.0 + 9.0 * (t >= 60minute)

#---

tend = 80minute
ts = range(0, tend, 201)

@named syss3 = make_model(;
    glcrhs=add_glucose(t),
    rf1=add_oligomycin(t),
    rhleak=add_fccp(t)
)

probs3 = ODEProblem(syss3, [], tend)
probs3DM = make_dm_prob(probs3; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)

tstops = [20minute, 40minute, 60minute]

sols3 = solve(probs3; tstops, saveat=ts)
solDMs3 = solve(probs3DM; tstops, saveat=ts)

figs3 = plot_figs2(sols3, solDMs3)
figs3

# TIF file
# `figs3.savefig("FigS3.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure S4
# Oxygen consumption in response to both glucose stimulation and chemical agents.

function plot_jo2(sol, solDM;
    labels=["Baseline", "Diabetic"],
    figsize=(6, 6),
    tight=true
)
    ts = sol.t
    tsm = ts ./ 60

    jo2 = sol[J_O2]
    jo2DM = solDM[J_O2]

    fig, ax = plt.subplots(; figsize)

    ax.plot(tsm, jo2, label=labels[1])
    ax.plot(tsm, jo2DM, label=labels[2])
    ax.set(xlabel="Time (minute)", ylabel="Rate (mM/s)", title="Oxygen consumption rate")

    ax.grid()
    ax.legend()
    fig.set_tight_layout(tight)
    return fig
end

#---

figs4 = plot_jo2(sols3, solDMs3)
figs4

# TIF file
# `figs4.savefig("FigS4.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure S5
# Baseline vs. Diabetic models using Glucose-Oligomycin-Rotenone protocol.

@named sys2 = make_model(;
    glcrhs = add_glucose(t),
    rf1 = add_oligomycin(t),
    retc = add_rotenone(t)
)

#---

prob2 = ODEProblem(sys2, [], tend)
prob2DM = make_dm_prob(prob2; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
sols5 = solve(prob2, saveat=ts)
sols5DM = solve(prob2DM, saveat=ts)

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
