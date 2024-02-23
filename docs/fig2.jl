#===
# Figure 2

Steady-state solutions for a range of glucose concentrations and OXPHOS capacities by chemicals.
===#
using DifferentialEquations
using ModelingToolkit
using DisplayAs: PNG
using MitochondrialDynamics
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

#---
@named sys = make_model()
@unpack GlcConst, VmaxF1, VmaxETC, pHleak = sys

iGlc = indexof(sys.GlcConst, parameters(sys))
iVmaxF1 = indexof(sys.VmaxF1, parameters(sys))
iVmaxETC = indexof(sys.VmaxETC, parameters(sys))
ipHleak = indexof(sys.pHleak, parameters(sys))

prob = SteadyStateProblem(sys, [])

# Range for two parameters

rGlcF1 = range(3.0, 30.0, 51)
rGlcETC = range(3.0, 30.0, 51)
rGlcHL = range(4.0, 30.0, 51)
rF1 = range(0.1, 2.0, 51)
rETC = range(0.1, 2.0, 51)
rHL = range(0.1, 5.0, 51)

# Testing
prob.p
p = copy(prob.p)
newp = ModelingToolkit.MTKParameters(sys, [GlcConst => 10.0])
p

newprob = remake(prob, ps=Dict(GlcConst => 10.0))
newprob.ps[GlcConst] = 10.0
prob.ps[GlcConst]


function solve_fig3(glc, r, k, prob)
    newprob = remake(prob)
    newprob.ps[GlcConst] = glc
    newprob.ps[k] = r * prob.ps[k]
    return solve(newprob, DynamicSS(Rodas5()))
end

solsf1 = [solve_fig3(glc, r, VmaxF1, prob) for r in rF1, glc in rGlcF1];
solsetc = [solve_fig3(glc, r, VmaxETC, prob) for r in rETC, glc in rGlcETC];
solshl = [solve_fig3(glc, r, pHleak, prob) for r in rHL, glc in rGlcHL];

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
    zs=(solsf1, solsetc, solshl),
    extremes=((1.0, 1.8), (80.0, 180.0), (0.0, 60.0))
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
            ax = axes[row-1, col-1]

            ylabel = ylabels[row]

            mesh = ax.pcolormesh(
                xx, yy, map(f, z);
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

fig = plot_fig3(figsize=(13, 10));
fig |> PNG

# Export figure
## exportTIF(fig, "Fig2-2Dsteadystate.tif")
