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

import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14
## plt.matplotlib.rcParams["font.sans-serif"] = "Arial"
## plt.matplotlib.rcParams["font.family"] = "sans-serif"

#---

glc_step(t) = 5.0mM * (1 + (t >= 20minute) + (t >= 40minute))
@named sys = make_model(; glcrhs=glc_step(t))

tend = 60minute
prob = ODEProblem(sys, [], tend)
sol = solve(prob, tstops=[20minute, 40minute]);

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
solDM = solve(probDM, tstops=[20minute, 40minute]);

#---

function plot_figs2(sol, solDM;
    tspan=(0.0, 60minute),
    density=301,
    figsize=(12, 12),
    labels=["Baseline", "Diabetic"],
    tight=true
)
    ts = range(tspan[1], tspan[2], length=density)
    tsm = ts ./ 60

    g3p = sol.(ts, idxs=G3P) .* 1000
    pyr = sol.(ts, idxs=Pyr) .* 1000
    nadh_c = sol.(ts, idxs=NADH_c) .* 1000
    nadh_m = sol.(ts, idxs=NADH_m) .* 1000
    ca_c = sol.(ts, idxs=Ca_c) .* 1000
    ca_m = sol.(ts, idxs=Ca_m) .* 1000
    atp_c = sol.(ts, idxs=ATP_c) .* 1000
    adp_c = sol.(ts, idxs=ADP_c) .* 1000
    td = atp_c ./ adp_c
    dpsi = sol.(ts, idxs=ΔΨm) .* 1000
    k = sol.(ts, idxs=degavg)

    g3pDM = solDM.(ts, idxs=G3P) .* 1000
    pyrDM = solDM.(ts, idxs=Pyr) .* 1000
    nadh_cDM = solDM.(ts, idxs=NADH_c) .* 1000
    nadh_mDM = solDM.(ts, idxs=NADH_m) .* 1000
    ca_cDM = solDM.(ts, idxs=Ca_c) .* 1000
    ca_mDM = solDM.(ts, idxs=Ca_m) .* 1000
    atp_cDM = solDM.(ts, idxs=ATP_c) .* 1000
    adp_cDM = solDM.(ts, idxs=ADP_c) .* 1000
    tdDM = atp_cDM ./ adp_cDM
    dpsiDM = solDM.(ts, idxs=ΔΨm) .* 1000
    kDM = solDM.(ts, idxs=degavg)

    fig, ax = plt.subplots(3, 3; figsize)

    ax[0, 0].plot(tsm, g3p, label=labels[1])
    ax[0, 0].plot(tsm, g3pDM, label=labels[2])
    ax[0, 0].set(ylabel="G3P (μM)")
    ax[0, 0].set_title("(A)", loc="left")

    ax[0, 1].plot(tsm, pyr, label=labels[1])
    ax[0, 1].plot(tsm, pyrDM, label=labels[2])
    ax[0, 1].set(ylabel="Pyruvate (μM)")
    ax[0, 1].set_title("(B)", loc="left")

    ax[0, 2].plot(tsm, nadh_c, label=labels[1])
    ax[0, 2].plot(tsm, nadh_cDM, label=labels[2])
    ax[0, 2].set(ylabel="Cytosolic NADH (μM)")
    ax[0, 2].set_title("(C)", loc="left")

    ax[1, 0].plot(tsm, nadh_m, label=labels[1])
    ax[1, 0].plot(tsm, nadh_mDM, label=labels[2])
    ax[1, 0].set(ylabel="Mitochondrial NADH (μM)")
    ax[1, 0].set_title("(D)", loc="left")

    ax[1, 1].plot(tsm, ca_c, label=labels[1])
    ax[1, 1].plot(tsm, ca_cDM, label=labels[2])
    ax[1, 1].set(ylabel="Cytosolic Calcium (μM)")
    ax[1, 1].set_title("(E)", loc="left")

    ax[1, 2].plot(tsm, ca_m, label=labels[1])
    ax[1, 2].plot(tsm, ca_mDM, label=labels[2])
    ax[1, 2].set(ylabel="Mitochondrial Calcium (μM)")
    ax[1, 2].set_title("(F)", loc="left")

    ax[2, 0].plot(tsm, td, label=labels[1])
    ax[2, 0].plot(tsm, tdDM, label=labels[2])
    ax[2, 0].set(ylabel="ATP:ADP")
    ax[2, 0].set_title("(G)", loc="left")

    ax[2, 1].plot(tsm, dpsi, label=labels[1])
    ax[2, 1].plot(tsm, dpsiDM, label=labels[2])
    ax[2, 1].set(ylabel="ΔΨm (mV)")
    ax[2, 1].set_title("(H)", loc="left")

    ax[2, 2].plot(tsm, k, label=labels[1])
    ax[2, 2].plot(tsm, kDM, label=labels[2])
    ax[2, 2].set(ylabel="Average Node Degree")
    ax[2, 2].set_title("(I)", loc="left")

    for i in 0:2, j in 0:2
        ax[i, j].grid()
        ax[i, j].legend()
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

md"""
## Figure S3

Changes in response to both glucose stimulation and chemical agents.

Using the Glucose-Oligomycin-FCCP protocol.
"""

add_glucose(t) = 5mM + 15mM * (t >= 20minute)
add_oligomycin(t) = 1.0 - 0.95 * (t >= 40minute)
add_rotenone(t) = 1.0 - 0.95 * (t >= 60minute)
add_fccp(t) = 1.0 + 9.0 * (t >= 60minute)

#---

tend = 80minute

@named syss3 = make_model(;
    glcrhs=add_glucose(t),
    rf1=add_oligomycin(t),
    rhleak=add_fccp(t)
)

probs3 = ODEProblem(syss3, [], tend)
probs3DM = make_dm_prob(probs3; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)

tstops = [20minute, 40minute, 60minute]

sols3 = solve(probs3; tstops)
solDMs3 = solve(probs3DM; tstops)

figs3 = plot_figs2(sols3, solDMs3; tspan=(0.0, tend))
figs3

# TIF file
# `figs3.savefig("FigS3.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure S4
# Oxygen consumption in response to both glucose stimulation and chemical agents.

function plot_jo2(sol, solDM;
    tspan=(0.0, 80minute),
    labels=["Baseline", "Diabetic"],
    density=301,
    figsize=(6, 6),
    tight=true
)
    ts = range(tspan[1], tspan[2], length=density)
    tsm = ts ./ 60

    jo2 = sol.(ts, idxs=J_O2)
    jo2DM = solDM.(ts, idxs=J_O2)

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

figs4 = plot_jo2(sols3, solDMs3, tspan=(0.0, tend))
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
sols5 = solve(prob2)
sols5DM = solve(prob2DM)

figs5 = plot_figs2(sols5, sols5DM; tspan=(0.0, tend))
figs5

# TIF file
# `figs5.savefig("FigS5.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure S6
# Oxygen consumption of Baseline and Diabetic models using the protocol from Figure S5.

figs6 = plot_jo2(sols5, sols5DM)
figs6

# TIF file
# `figs6.savefig("FigS6.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`
