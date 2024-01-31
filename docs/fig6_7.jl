#===
# Figure 6 and 7
===#
using OrdinaryDiffEq
using DiffEqCallbacks
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

#---
glc = 3.0:0.5:30.0
@named sys = make_model()
prob = ODEProblem(sys, [], Inf)

# Parameters
@unpack GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxGlc = indexof(GlcConst, parameters(sys))
idxVmaxPDH = indexof(VmaxPDH, parameters(sys))
idxpHleak = indexof(pHleak, parameters(sys))
idxVmaxF1 =  indexof(VmaxF1, parameters(sys))
idxVmaxETC =  indexof(VmaxETC, parameters(sys))

# ## Fig 6
function remake_rotenone(prob; rETC=0.1)
    p = copy(prob.p)
    p[idxVmaxETC] = prob.p[idxVmaxETC] * rETC
    return remake(prob, p=p)
end

function remake_oligomycin(prob; rF1=0.1)
    p = copy(prob.p)
    p[idxVmaxF1] = prob.p[idxVmaxF1] * rF1
    return remake(prob, p=p)
end

function remake_fccp(prob; rHL=5)
    p = copy(prob.p)
    p[idxpHleak] = prob.p[idxpHleak] * rHL
    return remake(prob, p=p)
end

function remake_dm(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxETC] = rETC * prob.p[idxVmaxETC]
    p[idxVmaxF1] = rF1 * prob.p[idxVmaxF1]
    p[idxpHleak] = rHL * prob.p[idxpHleak]
    p[idxVmaxPDH] = rPDH * prob.p[idxVmaxPDH]
    return remake(prob, p=p)
end

function prob_func_glc(prob, i, repeat)
    prob.p[idxGlc] = glc[i]
    prob
end

# DM cells

prob_dm = remake_dm(prob)
alg = TRBDF2()
prob_func=prob_func_glc
opt = (
    trajectories=length(glc),
    callback=TerminateSteadyState(),
    save_everystep=false,
    save_start = false,
)

sols = solve(EnsembleProblem(prob; prob_func), alg; opt...)
solsDM = solve(EnsembleProblem(prob_dm; prob_func), alg; opt...);

#---
function plot_fig6(sols, solsDM, glc; figsize=(10, 10), labels=["Baseline", "Diabetic"])
    glc5 = glc ./ 5
    numrows = 3
    numcols = 3

    fig, ax = plt.subplots(numrows, numcols; figsize)

    @unpack G3P = sys
    ax[0, 0].plot(glc5, extract(sols, G3P * 1000), label=labels[1])
    ax[0, 0].plot(glc5, extract(solsDM, G3P * 1000), label=labels[2])
    ax[0, 0].set_title("(A) G3P", loc="left")
    ax[0, 0].set(ylabel="Conc. (μM)")

    @unpack Pyr = sys
    ax[0, 1].plot(glc5, extract(sols, Pyr * 1000), label=labels[1])
    ax[0, 1].plot(glc5, extract(solsDM, Pyr * 1000), label=labels[2])
    ax[0, 1].set_title("(B) Pyr", loc="left")
    ax[0, 1].set(ylabel="Conc. (μM)")

    @unpack NADH_c, NAD_c = sys
    ax[0, 2].plot(glc5, extract(sols, NADH_c/NAD_c), label=labels[1])
    ax[0, 2].plot(glc5, extract(solsDM, NADH_c/NAD_c), label=labels[2])
    ax[0, 2].set_title("(C) NADH:NAD (cyto)", loc="left")

    @unpack NADH_m, NAD_m = sys
    ax[1, 0].plot(glc5, extract(sols, NADH_m/NAD_m), label=labels[1])
    ax[1, 0].plot(glc5, extract(solsDM, NADH_m/NAD_m), label=labels[2])
    ax[1, 0].set_title("(D) NADH:NAD (mito)", loc="left")

    @unpack Ca_c = sys
    ax[1, 1].plot(glc5, extract(sols, Ca_c * 1000), label=labels[1])
    ax[1, 1].plot(glc5, extract(solsDM, Ca_c * 1000), label=labels[2])
    ax[1, 1].set_title("(E) Calcium (cyto)", loc="left")
    ax[1, 1].set(ylabel="Conc. (μM)")

    @unpack Ca_m = sys
    ax[1, 2].plot(glc5, extract(sols, Ca_m * 1000), label=labels[1])
    ax[1, 2].plot(glc5, extract(solsDM, Ca_m * 1000), label=labels[2])
    ax[1, 2].set_title("(F) Calcium (mito)", loc="left")
    ax[1, 2].set(ylabel="Conc. (μM)")

    @unpack ΔΨm = sys
    ax[2, 0].plot(glc5, extract(sols, ΔΨm * 1000), label=labels[1])
    ax[2, 0].plot(glc5, extract(solsDM, ΔΨm * 1000), label=labels[2])
    ax[2, 0].set_title("(G) ΔΨ", loc="left")
    ax[2, 0].set(xlabel="Glucose (X)", ylabel="mV")

    @unpack ATP_c, ADP_c = sys
    ax[2, 1].plot(glc5, extract(sols, ATP_c/ADP_c), label=labels[1])
    ax[2, 1].plot(glc5, extract(solsDM, ATP_c/ADP_c), label=labels[2])
    ax[2, 1].set_title("(H) ATP:ADP", loc="left")
    ax[2, 1].set(xlabel="Glucose (X)")

    @unpack degavg = sys
    ax[2, 2].plot(glc5, extract(sols, degavg), label=labels[1])
    ax[2, 2].plot(glc5, extract(solsDM, degavg), label=labels[2])
    ax[2, 2].set_title("(I) Avg, Node degree", loc="left")
    ax[2, 2].set(xlabel="Glucose (X)")

    for i in 0:numrows-1, j in 0:numcols-1
        ax[i, j].grid()
        ax[i, j].legend()
    end

    fig.tight_layout()
    return fig
end

#---
fig6 = plot_fig6(sols, solsDM, glc);
fig6 |> PNG

# Export figure
exportTIF(fig6, "Fig6.tif")

# ## Figure 7
prob_fccp = remake_fccp(prob)
prob_rotenone = remake_rotenone(prob)
prob_oligomycin = remake_oligomycin(prob)

#---
sols = solve(EnsembleProblem(prob; prob_func), alg; opt...)
solsDM = solve(EnsembleProblem(prob_dm; prob_func), alg; opt...)
solsFCCP = solve(EnsembleProblem(prob_fccp; prob_func), alg; opt...)
solsRot = solve(EnsembleProblem(prob_oligomycin; prob_func), alg; opt...)
solsOligo = solve(EnsembleProblem(prob_rotenone; prob_func), alg; opt...);

#---
function plot_fig7(sols, solsDM, solsFCCP, solsRot, solsOligo, glc; figsize=(12, 6))
    sys = sols[begin].prob.f.sys
    @unpack J_HL, J_ANT = sys
    ## Gather ATP synthesis rate (fusion) and proton leak rate (fission)
    jHL_baseline = extract(sols, J_HL)
    jANT_baseline = extract(sols, J_ANT)
    ff_baseline = jANT_baseline ./ jHL_baseline
    jHL_dm = extract(solsDM, J_HL)
    jANT_dm = extract(solsDM, J_ANT)
    ff_dm = jANT_dm ./ jHL_dm
    jHL_fccp = extract(solsFCCP, J_HL)
    jANT_fccp = extract(solsFCCP, J_ANT)
    ff_fccp = jANT_fccp ./ jHL_fccp
    jHL_rot = extract(solsRot, J_HL)
    jANT_rot = extract(solsRot, J_ANT)
    ff_rot = jANT_rot ./ jHL_rot
    jHL_oligo = extract(solsOligo, J_HL)
    jANT_oligo = extract(solsOligo, J_ANT)
    ff_oligo = jANT_oligo ./ jHL_oligo

    glc5 = glc ./ 5

    fig, ax = plt.subplots(1, 2; figsize)

    ax[0].plot(glc5, ff_baseline, "b-", label="Baseline")
    ax[0].plot(glc5, ff_dm, "r--", label="Diabetic")
    ax[0].plot(glc5, ff_rot, "g--", label="Rotenone")
    ax[0].plot(glc5, ff_oligo, "c--", label="Oligomycin")
    ax[0].plot(glc5, ff_fccp, "k--", label="Uncoupler")
    ax[0].set(xlabel="Glucose (X)", ylabel="Fusion rate / Fission rate", xlim=(0.0, 6.0), ylim=(0.0, 2.5))
    ax[0].set_title("(A)", loc="left")
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(jHL_baseline, jANT_baseline, "bo-", label="Baseline")
    ax[1].plot(jHL_dm, jANT_dm, "ro-", label="Diabetic")
    ax[1].plot(jHL_rot, jANT_rot, "go-", label="Rotenone")
    ax[1].plot(jHL_oligo, jANT_oligo, "co-", label="Oligomycin")
    ax[1].plot(jHL_fccp, jANT_fccp, "ko-", label="Uncoupler")
    ax[1].set(xlabel="Proton leak rate (mM/s)", ylabel="ATP synthase rate (mM/s)", xlim=(0.0, 0.45), ylim=(0.0, 0.15))
    ax[1].set_title("(B)", loc="left")
    ax[1].grid()
    ax[1].legend()

    fig.tight_layout()
    return fig
end

#---

fig7 = plot_fig7(sols, solsDM, solsFCCP, solsRot, solsOligo, glc);
fig7 |> PNG

# Export figure
exportTIF(fig7, "Fig7.tif")
