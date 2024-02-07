#===
# Glucose-Oligomycin-Rotenone
## Step responses to both glucose stimulation and chemical agents
===#
using OrdinaryDiffEq
using DiffEqCallbacks
using ModelingToolkit
using DisplayAs: PNG
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

@named sys = make_model()

@unpack GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxGlc = findfirst(isequal(GlcConst), parameters(sys))
idxVmaxPDH = findfirst(isequal(VmaxPDH), parameters(sys))
idxpHleak = findfirst(isequal(pHleak), parameters(sys))
idxVmaxF1 =  findfirst(isequal(VmaxF1), parameters(sys))
idxVmaxETC =  findfirst(isequal(VmaxETC), parameters(sys))

tend = 100minute
ts = range(0, tend, 401)

probs5 = ODEProblem(sys, [], ts[end])

function remake_dm(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxETC] *= rETC
    p[idxVmaxF1] *= rF1
    p[idxpHleak] *= rHL
    p[idxVmaxPDH] *= rPDH
    return remake(prob, p=p)
end

prob_dmS5 = remake_dm(probs5)

# Define events
function add_glucose!(i)
    i.p[idxGlc] += 15mM
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

#---
function plot_figs2(sol, solDM; figsize=(12, 12), labels=["Baseline", "Diabetic"])
    @unpack G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ATP_c, ADP_c, ΔΨm, degavg, J_O2 = sol.prob.f.sys
    ts = sol.t
    tsm = ts ./ 60

    g3p = sol[G3P * 1000]
    jo2 = sol[J_O2]
    nadh_c = sol[NADH_c * 1000]
    nadh_m = sol[NADH_m * 1000]
    ca_c = sol[Ca_c * 1000]
    ca_m = sol[Ca_m * 1000]
    td = sol[ATP_c/ADP_c]
    dpsi = sol[ΔΨm * 1000]
    k = sol[degavg]

    g3pDM = solDM[G3P * 1000]
    jo2DM = solDM[J_O2]
    nadh_cDM = solDM[NADH_c * 1000]
    nadh_mDM = solDM[NADH_m * 1000]
    ca_cDM = solDM[Ca_c * 1000]
    ca_mDM = solDM[Ca_m * 1000]
    tdDM = solDM[ATP_c/ADP_c]
    dpsiDM = solDM[ΔΨm * 1000]
    kDM = solDM[degavg]

    numrows = 3
    numcols = 3
    fig, ax = plt.subplots(numrows, numcols; figsize)

    ax[0, 0].plot(tsm, g3p, label=labels[1])
    ax[0, 0].plot(tsm, g3pDM, label=labels[2])
    ax[0, 0].set(ylabel="G3P (μM)")
    ax[0, 0].set_title("(A)", loc="left")

    ax[0, 1].plot(tsm, jo2, label=labels[1])
    ax[0, 1].plot(tsm, jo2DM, label=labels[2])
    ax[0, 1].set(ylabel="OCR (mM/s)")
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

    for i in 0:numrows-1, j in 0:numcols-1
        ax[i, j].legend()
        ax[i, j].grid()
    end

    fig.tight_layout()
    return fig
end

# Start simulations
cbs = CallbackSet(add_glucose_cb, add_oligomycin_cb, add_rotenone_cb)
sols5 = solve(probs5, Rodas5(); callback=cbs, saveat=ts)
sols5DM = solve(prob_dmS5, Rodas5(); callback=cbs, saveat=ts)

figs5 = plot_figs2(sols5, sols5DM);
figs5 |> PNG

# TIFF file
exportTIF(figs5, "FigS1-Glucose-Oligomycin-Rotenone.tif")
