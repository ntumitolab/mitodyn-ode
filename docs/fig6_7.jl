#===
# Figure 6 and 7
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
# using MitochondrialDynamics: GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC, J_ANT, J_HL
# using MitochondrialDynamics: G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ΔΨm, ATP_c, ADP_c, degavg
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

glc = range(3.0mM, 30.0mM, length=201)
tend = 50minute
@named sys = make_model()
prob = SteadyStateProblem(sys, [])

#---

@unpack GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxGlc = findfirst(isequal(GlcConst), parameters(sys))
idxVmaxPDH = findfirst(isequal(VmaxPDH), parameters(sys))
idxpHleak = findfirst(isequal(pHleak), parameters(sys))
idxVmaxF1 =  findfirst(isequal(VmaxF1), parameters(sys))
idxVmaxETC =  findfirst(isequal(VmaxETC), parameters(sys))

# ## Fig 6

function make_dm_prob(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxPDH] *= rPDH
    p[idxpHleak] *= rHL
    p[idxVmaxF1] *= rF1
    p[idxVmaxETC] *= rETC
    return remake(prob, p=p)
end

function make_rotenone_prob(prob; rETC=0.1)
    p = copy(prob.p)
    prob.p[idxVmaxETC] *= rETC
    return remake(prob, p=p)
end

function make_oligomycin_prob(prob; rF1=0.1)
    p = copy(prob.p)
    p[idxVmaxF1] *= rF1
    return remake(prob, p=p)
end

function make_fccp_prob(prob; rHL=5)
    p = copy(prob.p)
    p[idxpHleak] *= rHL
    return remake(prob, p=p)
end

function remake_glc(prob, g)
    p = copy(prob.p)
    p[idxGlc] = g
    remake(prob; p=p)
end

#---

prob_dm = make_dm_prob(prob)

sols = map(glc) do g
    solve(remake_glc(prob, g), DynamicSS(Rodas5()))
end

solsDM = map(glc) do g
    solve(remake_glc(prob_dm, g), DynamicSS(Rodas5()))
end;

# TODO: less information

function plot_fig6(sols, solsDM, glc; figsize=(12, 12), tight=true, labels=["Baseline", "Diabetic"])
    extract(sols, k, scale=1) = map(s->s[k] * scale, sols)

    @unpack ATP_c, ADP_c = sol.prob.f.sys
    glc5 = glc ./ 5
    td = extract(sols, ATP_c/ADP_c)
    tdDM = extract(solsDM, ATP_c/ADP_c)

    fig, ax = plt.subplots(3, 3; figsize)

    sys = sol.prob.f.sys
    @unpack G3P = sys
    ax[1, 1].plot(glc5, extract(sols, G3P, 1000), label=labels[1])
    ax[1, 1].plot(glc5, extract(solsDM, G3P, 1000), label=labels[2])
    ax[1, 1].set_title("(A) G3P", loc="left")
    ax[1, 1].set(ylabel="Conc. (μM)")

    @unpack Pyr = sys
    ax[1, 2].plot(glc5, extract(sols, Pyr, 1000), label=labels[1])
    ax[1, 2].plot(glc5, extract(solsDM, Pyr, 1000), label=labels[2])
    ax[1, 2].set_title("(B) Pyr", loc="left")
    ax[1, 2].set(ylabel="Conc. (μM)")

    @unpack NADH_c, NAD_c = sys
    ax[1, 3].plot(glc5, extract(sols, NADH_c/NAD_c), label=labels[1])
    ax[1, 3].plot(glc5, extract(solsDM, NADH_c/NAD_c), label=labels[2])
    ax[1, 3].set_title("(C) NADH:NAD (cyto)", loc="left")

    @unpack NADH_m, NAD_m = sys
    ax[2, 1].plot(glc5, extract(sols, NADH_m/NAD_m), label=labels[1])
    ax[2, 1].plot(glc5, extract(solsDM, NADH_m/NAD_m), label=labels[2])
    ax[2, 1].set_title("(D) NADH:NAD (mito)", loc="left")

    @unpack Ca_c = sys
    ax[2, 2].plot(glc5, extract(sols, Ca_c, 1000), label=labels[1])
    ax[2, 2].plot(glc5, extract(solsDM, Ca_c, 1000), label=labels[2])
    ax[2, 2].set_title("(E) Calcium (cyto)", loc="left")
    ax[2, 2].set(ylabel="Conc. (μM)")

    @unpack Ca_m = sys
    ax[2, 3].plot(glc5, extract(sols, Ca_m, 1000), label=labels[1])
    ax[2, 3].plot(glc5, extract(solsDM, Ca_m, 1000), label=labels[2])
    ax[2, 3].set_title("(F) Calcium (mito)", loc="left")
    ax[2, 3].set(ylabel="Conc. (μM)")

    @unpack ΔΨm = sys
    ax[3, 1].plot(glc5, extract(sols, ΔΨm, 1000), label=labels[1])
    ax[3, 1].plot(glc5, extract(solsDM, ΔΨm, 1000), label=labels[2])
    ax[3, 1].set_title("(G) ΔΨ", loc="left")
    ax[3, 1].set(xlabel="Glucose (X)", ylabel="mV")

    @unpack ATP_c, ADP_c = sys
    ax[3, 2].plot(glc5, extract(sols, ATP_c/ADP_c), label=labels[1])
    ax[3, 2].plot(glc5, extract(solsDM, ATP_c/ADP_c), label=labels[2])
    ax[3, 2].set_title("(H) ATP:ADP", loc="left")
    ax[3, 2].set(xlabel="Glucose (X)")

    @unpack degavg = sys
    ax[3, 3].plot(glc5, extract(sols, degavg), label=labels[1])
    ax[3, 3].plot(glc5, extract(solsDM, degavg), label=labels[2])
    ax[3, 3].set_title("(I) Average node degree", loc="left")
    ax[3, 3].set(xlabel="Glucose (X)")

    for a in ax
        a.grid()
    end

    fig.set_tight_layout(tight)
    return fig
end

#---

fig6 = plot_fig6(sols, solsDM, glc)
fig6

# Generating tiff file
# `fig6.savefig("Fig6.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`

# ## Figure 7

prob_fccp = make_fccp_prob(prob)
prob_rotenone = make_rotenone_prob(prob)
prob_oligomycin = make_oligomycin_prob(prob)

solsFCCP = map(glc) do g
    solve(remake_glc(prob_fccp, g), DynamicSS(Rodas5()))
end

solsRot = map(glc) do g
    solve(remake_glc(prob_rotenone, g), DynamicSS(Rodas5()))
end

solsOligo = map(glc) do g
    solve(remake_glc(prob_oligomycin, g), DynamicSS(Rodas5()))
end

#---

function plot_fig7(sols, solsDM, solsFCCP, solsRot, solsOligo, glc;
    figsize=(12, 6), tight=true
)
    ## Gather ATP synthesis rate (fusion) and proton leak rate (fission)
    jHL_baseline = getindex.(sols, J_HL)
    jANT_baseline = getindex.(sols, J_ANT)
    ff_baseline = jANT_baseline ./ jHL_baseline
    jHL_dm = getindex.(solsDM, J_HL)
    jANT_dm = getindex.(solsDM, J_ANT)
    ff_dm = jANT_dm ./ jHL_dm
    jHL_fccp = getindex.(solsFCCP, J_HL)
    jANT_fccp = getindex.(solsFCCP, J_ANT)
    ff_fccp = jANT_fccp ./ jHL_fccp
    jHL_rot = getindex.(solsRot, J_HL)
    jANT_rot = getindex.(solsRot, J_ANT)
    ff_rot = jANT_rot ./ jHL_rot
    jHL_oligo = getindex.(solsOligo, J_HL)
    jANT_oligo = getindex.(solsOligo, J_ANT)
    ff_oligo = jANT_oligo ./ jHL_oligo

    glc5 = glc ./ 5

    fig, ax = plt.subplots(1, 2; figsize)

    ax[1].plot(glc5, ff_baseline, "b-", label="Baseline")
    ax[1].plot(glc5, ff_dm, "r--", label="Diabetic")
    ax[1].plot(glc5, ff_rot, "g--", label="Rotenone")
    ax[1].plot(glc5, ff_oligo, "c--", label="Oligomycin")
    ax[1].plot(glc5, ff_fccp, "k--", label="Uncoupler")
    ax[1].set(xlabel="Glucose (X)", ylabel="Fusion rate / Fission rate", xlim=(0.0, 6.0), ylim=(0.0, 2.5))
    ax[1].set_title("(A)", loc="left")
    ax[1].grid()
    ax[1].legend()

    ax[2].plot(jHL_baseline, jANT_baseline, "bo-", label="Baseline")
    ax[2].plot(jHL_dm, jANT_dm, "ro-", label="Diabetic")
    ax[2].plot(jHL_rot, jANT_rot, "go-", label="Rotenone")
    ax[2].plot(jHL_oligo, jANT_oligo, "co-", label="Oligomycin")
    ax[2].plot(jHL_fccp, jANT_fccp, "ko-", label="Uncoupler")
    ax[2].set(xlabel="Proton leak rate (mM/s)", ylabel="ATP synthase rate (mM/s)", xlim=(0.0, 0.45), ylim=(0.0, 0.15))
    ax[2].set_title("(B)", loc="left")
    ax[2].grid()
    ax[2].legend()

    fig.set_tight_layout(tight)
    return fig
end

#---

fig7 = plot_fig7(sols, solsDM, solsFCCP, solsRot, solsOligo, glc)
fig7

# Generating tiff file
# `fig7.savefig("Fig7.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`
