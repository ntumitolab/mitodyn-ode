md"""
# Figure 2

Steady-state solutions across a range of glucose levels.
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: GlcConst, G3P, Pyr, NADH_c, NADH_m, ATP_c, ADP_c, AMP_c, Ca_m, Ca_c, x, ΔΨm, degavg
import MitochondrialDynamics: second, μM, mV, mM, Hz
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14
## plt.matplotlib.rcParams["font.sans-serif"] = "Arial"
## plt.matplotlib.rcParams["font.family"] = "sans-serif"

#---

@named sys = make_model()

#---

pidx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
idxGlc = pidx[GlcConst]

function remake_glc(prob, g)
    p = copy(prob.p)
    p[idxGlc] = g
    remake(prob; p=p)
end

#---

##  using default u0
prob = SteadyStateProblem(sys, [])

glc = range(3.0mM, 30.0mM, length=101)  # Range of glucose

sols = map(glc) do g
    solve(remake_glc(prob, g), alg=DynamicSS(Rodas5()))
end;

#---

extract(sols, k) = map(s->s[k], sols)

function plot_fig2(glc, sols; figsize=(12, 12))
    glc5 = glc ./ 5
    g3p = extract(sols, G3P) .* 1000
    pyr = extract(sols, Pyr) .* 1000
    ca_c = extract(sols, Ca_c) .* 1000
    ca_m = extract(sols, Ca_m) .* 1000
    nadh_c = extract(sols, NADH_c) .* 1000
    nadh_m = extract(sols, NADH_m) .* 1000
    atp_c = extract(sols, ATP_c) .* 1000
    adp_c = extract(sols, ADP_c) .* 1000
    amp_c = extract(sols, AMP_c) .* 1000
    dpsi = extract(sols, ΔΨm) .* 1000
    x1 = extract(sols, x[1])
    x2 = extract(sols, x[2])
    x3 = extract(sols, x[3])
    deg = extract(sols, degavg)

    fig, ax = plt.subplots(3, 3; figsize)

    ax[1, 1].plot(glc5, g3p)
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0.0, 10.0))
    ax[1, 2].plot(glc5, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0.0, 120.0))
    ax[1, 3].plot(glc5, ca_c, label="cyto")
    ax[1, 3].plot(glc5, ca_m, label="mito")
    ax[1, 3].legend()
    ax[1, 3].set(title="(C) Calcium (μM)", ylim=(0.0, 1.5))
    ax[2, 1].plot(glc5, nadh_c, label="cyto")
    ax[2, 1].plot(glc5, nadh_m, label="mito")
    ax[2, 1].legend()
    ax[2, 1].set(title="(D) NADH (μM)")
    ax[2, 2].plot(glc5, atp_c, label="ATP")
    ax[2, 2].plot(glc5, adp_c, label="ADP")
    ax[2, 2].plot(glc5, amp_c, label="AMP")
    ax[2, 2].legend()
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 3].plot(glc5, atp_c ./ adp_c)
    ax[2, 3].set(title="(F) ATP/ADP ratio")
    ax[3, 1].plot(glc5, dpsi, label="cyto")
    ax[3, 1].set(title="(G) ΔΨ (mV)", ylim=(80, 160), xlabel="Glucose (X)")
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

    plt.tight_layout()
    return fig
end

#---

fig2 = plot_fig2(glc, sols)
plt.gcf()

# Tiff figure
fig2.savefig("Fig2.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
