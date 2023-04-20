md"""
# Figure 3

Steady-state solutions for a range of glucose concentrations and OXPHOS capacities.
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

@named sys = make_model()
@unpack GlcConst, VmaxF1, VmaxETC, pHleak = sys
iGlc = findfirst(isequal(GlcConst), parameters(sys))
iVmaxF1 = findfirst(isequal(VmaxF1), parameters(sys))
iVmaxETC = findfirst(isequal(VmaxETC), parameters(sys))
ipHleak = findfirst(isequal(pHleak), parameters(sys))
prob = SteadyStateProblem(sys, [])

# TODO: use ensmeble simulation
# function solve_fig3(glc, r, protein, prob; alg=DynamicSS(Rodas5()))
#     idx = findfirst(isequal(protein), parameters(sys))
#     p = copy(prob.p)
#     p[iGlc] = glc
#     p[idx] = prob.p[idx] * r
#     return solve(remake(prob, p=p), alg)
# end

#---


## TODO: use ensmeble simulation

rGlcF1 = range(3.0, 30.0, 51)
rGlcETC = range(3.0, 30.0, 51)
rGlcHL = range(4.0, 30.0, 51)
rF1 = range(0.1, 2.0, 51)
rETC = range(0.1, 2.0, 51)
rHL = range(0.1, 5.0, 51)

#---

Cf1 = CartesianIndices((length(rF1), length(rGlcF1)))
Cetc = CartesianIndices((length(rGlcETC), length(rETC)))
Chl = CartesianIndices((length(rGlcHL), length(rHL)))

# https://discourse.julialang.org/t/julia-usage-how-to-get-2d-indexes-from-1d-index-when-accessing-a-2d-array/61440/2
probf_f1 = function (prob, i, repeat)
    glc = Cf1[i][2]
    rf1 = Cf1[i][1]
    prob.p[iGlc] = rGlcF1[glc]
    prob.p[iVmaxF1] *= rF1[rf1]
    return prob
end

probf_etc = function (prob, i, repeat)
    glc = Cetc[i][2]
    retc = Cetc[i][1]
    prob.p[iGlc] = rGlcETC[glc]
    prob.p[iVmaxETC] *= rETC[retc]
    return prob
end

probf_hl = function (prob, i, repeat)
    glc = Chl[i][2]
    rhl = Chl[i][1]
    prob.p[iGlc] = rGlcHL[glc]
    prob.p[ipHleak] *= rHL[rhl]
    return prob
end

#---

alg = DynamicSS(Rodas5())
eprobf1 = EnsembleProblem(prob; prob_func = probf_f1)
simf1 = solve(eprobf1, alg; trajectories=length(Cf1))
eprobetc = EnsembleProblem(prob; prob_func = probf_etc)
simetc =solve(eprobetc, alg; trajectories=length(Cetc))
eprobhl = EnsembleProblem(prob; prob_func = probf_hl)
simhl =solve(eprobhl, alg; trajectories=length(Chl))

#---

function plot_fig3(;
    figsize=(10, 10),
    cmaps=["bwr", "magma", "viridis"],
    ylabels=[
        "ATP synthase capacity (X)",
        "ETC capacity (X)",
        "Proton leak rate (X)"
    ],
    cbarlabels=["<k>", "ΔΨ", "ATP/ADP"],
    xxs=(rGlcF1, rGlcETC, rGlcHL),
    xscale=5.0,
    yys=(rF1, rETC, rHL),
    zs=(simf1, simetc, simhl),
    extremes=((1.0, 2.0), (80.0, 180.0), (0.0, 60.0))
)
    ## mapping functions
    @unpack degavg, ΔΨm, ATP_c, ADP_c = sys
    fs = (s -> s[degavg], s -> s[ΔΨm * 1000] , s -> s[ATP_c / ADP_c])

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
                xx, yy, reshape(map(f, z), length(xx), length(yy));
                shading="gouraud",
                rasterized=true,
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

    fig.tight_layout()
    return fig
end

#---

fig3 = plot_fig3(figsize=(13, 10))
fig3

# TIFF file
# `fig3.savefig("Fig3.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`
