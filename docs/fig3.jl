md"""
# Figure 3

Steady-state solutions for a range of glucose concentrations and OXPHOS capacities.
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
import PythonPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

@named sys = make_model()
@unpack GlcConst = sys
iGlc = findfirst(isequal(GlcConst), parameters(sys))

prob = SteadyStateProblem(sys, [])

function solve_fig3(glc, r, protein, prob; alg=DynamicSS(Rodas5()))
    idx = findfirst(isequal(protein), parameters(sys))
    p = copy(prob.p)
    p[iGlc] = glc
    p[idx] = prob.p[idx] * r
    return solve(remake(prob, p=p), alg)
end

#---

## TODO: use ensmeble simulation

rGlc1 = range(3.0, 30.0, 50)
rGlc2 = range(4.0, 30.0, 50)
rF1 = range(0.1, 2.0, 50)
rETC = range(0.1, 2.0, 50)
rHL = range(0.1, 5.0, 50)

@unpack VmaxF1, VmaxETC, pHleak = sys
uInf_f1 = [solve_fig3(glc, r, VmaxF1, prob) for r in rF1, glc in rGlc1];
uInf_etc = [solve_fig3(glc, r, VmaxETC, prob) for r in rETC, glc in rGlc1];
uInf_hl = [solve_fig3(glc, r, pHleak, prob) for r in rHL, glc in rGlc2];

#---

function plot_fig3(;
    figsize=(10, 10),
    levels=40,
    cmaps=["bwr", "magma", "viridis"],
    ylabels=[
        "ATP synthase capacity (X)",
        "ETC capacity (X)",
        "Proton leak rate (X)"
    ],
    cbarlabels=["<k>", "ΔΨ", "ATP/ADP"],
    xxs=(rGlc1, rGlc1, rGlc2),
    xscale=5.0,
    yys=(rF1, rETC, rHL),
    zs=(uInf_f1, uInf_etc, uInf_hl),
    extremes=((1.0, 2.0), (80.0, 180.0), (0.0, 60.0))
)
    ## mapping functions
    @unpack degavg, ΔΨm, ATP_c, ADP_c = sys
    fs = (s -> s[degavg], s -> s[ΔΨm] * 1000, s -> s[ATP_c] / s[ADP_c])

    fig, axes = plt.subplots(3, 3; figsize)

    for col in 1:3
        f = fs[col]
        cm = cmaps[col]
        cbl = cbarlabels[col]
        vmin, vmax = extremes[col]

        ## lvls = LinRange(vmin, vmax, levels)
        for row in 1:3
            xx = xxs[row] ./ xscale
            yy = yys[row]
            z = zs[row]
            ax = axes[row, col]

            ylabel = ylabels[row]

            mesh = ax.pcolormesh(
                xx, yy, map(f, z);
                shading="gouraud",
                rasterized=true,
                ## levels=levels,
                cmap=cm,
                vmin=vmin,
                vmax=vmax
            )

            ax.set(ylabel=ylabel, xlabel="Glucose (X)")

            ## Arrow annotation: https://matplotlib.org/stable/tutorials/text/annotations.html#plotting-guide-annotation
            if row == 1
                ax.text(5.5, 1, "Oligomycin", ha="center", va="center", rotation=-90, size=16, bbox=Dict("boxstyle" => "rarrow", "fc" => "w", "ec" => "k", "lw" => 2, "alpha" => 0.5))
            elseif row == 2
                ax.text(5.5, 1, "Rotenone", ha="center", va="center", rotation=-90, size=16, bbox=Dict("boxstyle" => "rarrow", "fc" => "w", "ec" => "k", "lw" => 2, "alpha" => 0.5))
            elseif row == 3
                ax.text(5.5, 2.5, "FCCP", ha="center", va="center", rotation=90, size=16, bbox=Dict("boxstyle" => "rarrow", "fc" => "w", "ec" => "k", "lw" => 2, "alpha" => 0.5))
            end
            cbar = fig.colorbar(mesh, ax=ax)
            cbar.ax.set_title(cbl)
        end
    end

    plt.tight_layout()
    return fig
end

#---

fig3 = plot_fig3(figsize=(13, 10))
plt.gcf()

# TIFF file
fig3.savefig("Fig3.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
