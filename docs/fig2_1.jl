md"""
# Free fatty acid and galactose treatment

- Free fatty acid: adding more citric acid cycle (CAC) flux
- Galactose: Tweak the stoichiometry of glucokinase (GK) so that glycolysis yields no ATP.
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: GlcConst, G3P, Pyr, NADH_c, NADH_m, ATP_c, ADP_c, AMP_c, Ca_m, Ca_c, x, ΔΨm, degavg
import MitochondrialDynamics: second, μM, mV, mM, Hz
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"
using MitochondrialDynamics: GlcConst, G3P, Pyr, NADH_c, NADH_m, ATP_c, ADP_c, AMP_c, Ca_m, Ca_c, x, ΔΨm, degavg, J_CAC

#---

@named sys = make_model()
prob = SteadyStateProblem(sys, []) ## Use default u0
sol = solve(prob)

#---

pidx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
idxGlc = pidx[GlcConst]

function remake_glc(prob, g)
    p = copy(prob.p)
    p[idxGlc] = g
    remake(prob; p=p)
end

sols = map(glc) do g
    solve(remake_glc(prob, g), DynamicSS(Rodas5()))
end

# Additional CAC flux: 3.4 μM/s
@named sysffa = make_model(j_ffa=sol[J_CAC] * 0.5)
probffa = SteadyStateProblem(sysffa, [])
solsffa = map(glc) do g
    solve(remake_glc(probffa, g), DynamicSS(Rodas5()))
end

@named sys_gal = make_model(gk_atp_stoich=4)
prob_gal = SteadyStateProblem(sys_gal, [])
sols_gal = map(glc) do g
    solve(remake_glc(prob_gal, g), DynamicSS(Rodas5()))
end;

#---

function plot_ffa(
    sols, solsffa, glc;
    figsize=(12, 12),
    tight=true,
    labels=["Default", "FFA"]
)
    glc5 = glc ./ 5
    fig, ax = plt.subplots(3, 3; figsize)

    ax[1, 1].plot(glc5, extract(sols, G3P) .* 1000, label=labels[1])
    ax[1, 1].plot(glc5, extract(solsffa, G3P) .* 1000, label=labels[2])
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0, 10))

    ax[1, 2].plot(glc5, extract(sols, Pyr) .* 1000, label=labels[1])
    ax[1, 2].plot(glc5, extract(solsffa, Pyr) .* 1000, label=labels[2])
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0, 160))

    ax[1, 3].plot(glc5, extract(sols, Ca_c) .* 1000, label=labels[1])
    ax[1, 3].plot(glc5, extract(solsffa, Ca_c) .* 1000, label=labels[2])
    ax[1, 3].set(title="(C) Cytosolic Calcium (μM)", ylim=(0.0, 0.4))

    ax[2, 1].plot(glc5, extract(sols, Ca_m) .* 1000, label=labels[1])
    ax[2, 1].plot(glc5, extract(solsffa, Ca_m) .* 1000, label=labels[2])
    ax[2, 1].set(title="(D) Mitochondrial Calcium (μM)", ylim=(0.0, 1.5))

    ax[2, 2].plot(glc5, extract(sols, NADH_c) .* 1000, label=labels[1])
    ax[2, 2].plot(glc5, extract(solsffa, NADH_c) .* 1000, label=labels[2])
    ax[2, 2].set(title="(E) Cytosolic NADH (μM)")

    ax[2, 3].plot(glc5, extract(sols, NADH_m) .* 1000, label=labels[1])
    ax[2, 3].plot(glc5, extract(solsffa, NADH_m) .* 1000, label=labels[2])
    ax[2, 3].set(title="(F) Mitochondrial NADH (μM)")

    ax[3, 1].plot(glc5, extract(sols, ATP_c) ./ extract(sols, ADP_c), label=labels[1])
    ax[3, 1].plot(glc5, extract(solsffa, ATP_c) ./ extract(solsffa, ADP_c), label=labels[2])
    ax[3, 1].set(title="(G) ATP:ADP")

    ax[3, 2].plot(glc5, extract(sols, ΔΨm) .* 1000, label=labels[1])
    ax[3, 2].plot(glc5, extract(solsffa, ΔΨm) .* 1000, label=labels[2])
    ax[3, 2].set(title="(H) ΔΨ (mV)")

    ax[3, 3].plot(glc5, extract(sols, degavg), label=labels[1])
    ax[3, 3].plot(glc5, extract(solsffa, degavg), label=labels[2])
    ax[3, 3].set(title="(I) Average Degree")

    for a in ax
        a.set_xticks(1:6)
        a.grid()
        a.legend()
    end

    fig.set_tight_layout(tight)
    return fig
end

#---

fig_ffa = plot_ffa(sols, solsffa, glc; labels=["Default", "FFA"])
plt.gcf()

# Tiff file
fig_ffa.savefig("Fig_FFA.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))

#---

fig_gal = plot_ffa(sols, sols_gal, glc; labels=["Default", "Galactose"])
plt.gcf()

# Tiff file
fig_gal.savefig("Fig_GAL.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
