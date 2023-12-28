# # Step Response to elevated glucose concentrations

using DifferentialEquations
using ModelingToolkit
using MKL
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
using PythonCall
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14
## plt.matplotlib.rcParams["font.sans-serif"] = "Arial"
## plt.matplotlib.rcParams["font.family"] = "sans-serif"

# ## Step Response to elevated glucose concentrations

@variables t, Glc(t)
glc_step(t) = 5.0mM * (1 + (t >= 20minute) + (t >= 40minute))
@named sys = make_model(; glceq=Glc~glc_step(t))

ts = range(0, 60minute; step=0.5minute)
prob = ODEProblem(sys, [], ts[end])
sol = solve(prob, TRBDF2(), tstops=[20minute, 40minute], saveat=ts);

#---

function plot_figs1( sol; figsize=(10, 10))
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

    numrows = 3
    numcols = 3
    fig, ax = plt.subplots(numrows, numcols; figsize)

    ax[0, 0].plot(tsm, g3p)
    ax[0, 0].set(title="(A) G3P (μM)", ylim=(0.0, 8.0))

    ax[0, 1].plot(tsm, pyr)
    ax[0, 1].set(title="(B) Pyruvate (μM)", ylim=(0.0, 100.0))

    ax[0, 2].plot(tsm, ca_c, label="cyto")
    ax[0, 2].plot(tsm, ca_m, label="mito")
    ax[0, 2].set(title="(C) Calcium (μM)")
    ax[0, 2].legend()

    ax[1, 0].plot(tsm, nadh_c, label="cyto")
    ax[1, 0].plot(tsm, nadh_m, label="mito")
    ax[1, 0].set(title="(D) NADH (μM)")
    ax[1, 0].legend()

    ax[1, 1].plot(tsm, atp_c, label="ATP")
    ax[1, 1].plot(tsm, adp_c, label="ADP")
    ax[1, 1].plot(tsm, amp_c, label="AMP")
    ax[1, 1].set(title="(E) Adenylates (μM)")
    ax[1, 1].legend()

    ax[1, 2].plot(tsm, atp_c ./ adp_c)
    ax[1, 2].set(title="(F) ATP:ADP", ylim=(0, 45))

    ax[2, 0].plot(tsm, dpsi)
    ax[2, 0].set(title="(G) ΔΨ (mV)", ylim=(80, 150), xlabel="Time (minutes)")

    ax[2, 1].plot(tsm, x1, label="X1")
    ax[2, 1].plot(tsm, x2, label="X2")
    ax[2, 1].plot(tsm, x3, label="X3")
    ax[2, 1].set(title="(H) Mitochondrial nodes", xlabel="Time (minutes)")
    ax[2, 1].legend(loc="upper left")

    ax[2, 2].plot(tsm, k)
    ax[2, 2].set(title="(I) Average Node Degree", xlabel="Time (minutes)")

    for i in 0:numrows-1, j in 0:numcols-1
        ax[i, j].grid()
    end

    fig.tight_layout()
    return fig
end

#---

fig = plot_figs1(sol)

# Export figure
## figs1.savefig("FigS_glc_step.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))

#===
## Step sesponse to glucose addition in healthy and DM cells.

In the diabetic parameter set. We

- decreased PDH capacity
- decreased ETC capacity
- decreased ATP synthase capacity
- increases proton leak rate
===#

@named sys = make_model(; glceq=Glc~glc_step(t))

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
ts = range(0, 60minute; step=0.5minute)
solDM = solve(prob_dm, TRBDF2(), tstops=[20minute, 40minute], saveat=ts);

#---
function plot_figs2(sol, solDM; figsize=(10, 10), labels=["Baseline", "Diabetic"])
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
        ax[i, j].grid()
    end

    fig.tight_layout()
    return fig
end

figs2 = plot_figs2(sol, solDM)
