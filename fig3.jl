using DifferentialEquations
using LabelledArrays
using Parameters
using Setfield
using FromFile
@from "Model/Model.jl" using Model
@from "Model/utils.jl" import second, μM, mV, mM, Hz

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"


u0 = LVector(g3p = 2.8μM,
             pyr = 8.5μM,
             nadh_c = 1μM,
             nadh_m = 60μM,
             atp_c = 4000μM,
             adp_c = 500μM,
             ca_m = 0.250μM,
             dpsi = 100mV,
             x2=0.2, x3=0.05)

param0 = MitoDynNode()
ssprob = SteadyStateProblem(model!, u0, param0)

function _solve(glc, p::MitoDynNode; tspan=3000.0)
    ssalg= DynamicSS(Rodas5(), tspan=tspan)
    solve(SteadyStateProblem(model!, u0, setglc(p, glc)), ssalg)
end

function solve_f1(glc, r, param0)
    p = @set param0.f1fo.VMAX *= r
    _solve(glc, p)
end

function solve_etc(glc, r, param0)
    p = @set param0.etc.VMAX *= r
    _solve(glc, p)
end

function solve_hl(glc, r, param0)
    p = @set param0.hleak.P_H *= r
    _solve(glc, p)
end

rGlc1 = LinRange(3.0, 30.0, 50)
rGlc2 = LinRange(4.0, 30.0, 50)
rF1 = LinRange(0.1, 2.0, 50)
rETC = LinRange(0.1, 2.0, 50)
rHL = LinRange(0.1, 5.0, 50)

uInf_f1 = [solve_f1(glc, r, param0) for r in rF1, glc in rGlc1]
uInf_etc = [solve_etc(glc, r, param0) for r in rETC, glc in rGlc1]
uInf_hl = [solve_hl(glc, r, param0) for r in rHL, glc in rGlc2]

function plot_fig3(; figsize=(13,10),
                     levels = 40,
                     cmaps = ["bwr", "magma", "viridis"],
                     ylabels = ["ATP synthase capacity (X)",
                                "ETC capacity (X)",
                                "Proton leak rate (X)"],
                     cbarlabels = ["Average node degree",
			                       "Mitochondrial membrane potential",
		                           "ATP/ADP ratio"],
                     xxs = (rGlc1, rGlc1, rGlc2),
                     xscale = 5.0,
                     yys = (rF1, rETC, rHL),
                     zs  = (uInf_f1, uInf_etc, uInf_hl),
                     extremes=((1.0, 2.0),(80.0, 180.0),(0.0, 60.0)))


    # mapping functions
    fs = (s->avgdeg(s.u), s->s.u.dpsi*1000, s->s.u[:atp_c] / s.u[:adp_c])

    fig, axes = plt.subplots(3, 3, figsize=figsize)

    for col in 1:3
        f = fs[col]
        cm = cmaps[col]
        cbl = cbarlabels[col]
        vmin, vmax = extremes[col]

        # lvls = LinRange(vmin, vmax, levels)
        for row in 1:3
            xx = xxs[row] ./ xscale
            yy = yys[row]
            z  = zs[row]
            ax = axes[row, col]

            ylabel = ylabels[row]

            mesh = ax.pcolormesh(xx, yy, map(f, z),
                shading="gouraud", rasterized=true,
                # levels=levels,
                cmap=cm, vmin=vmin, vmax=vmax)

            ax.set(ylabel=ylabel, xlabel="Glucose (X)")

            # Arrow annotation
            # https://matplotlib.org/stable/tutorials/text/annotations.html#plotting-guide-annotation
            if row == 1
                ax.text(5.5, 1, "Oligomycin", ha="center", va="center", rotation=-90, size=16, bbox=Dict("boxstyle"=>"rarrow", "fc"=>"w", "ec"=>"k", "lw"=>2, "alpha"=>0.5))

                cbar = plt.colorbar(mesh, ax=ax, location="top")
                cbar.ax.set_title(cbl)
            elseif row == 2
                ax.text(5.5, 1, "Rotenone", ha="center", va="center", rotation=-90, size=16, bbox=Dict("boxstyle"=>"rarrow", "fc"=>"w", "ec"=>"k", "lw"=>2, "alpha"=>0.5))
            elseif row == 3
                ax.text(5.5, 2.5, "FCCP", ha="center", va="center", rotation=90, size=16, bbox=Dict("boxstyle"=>"rarrow", "fc"=>"w", "ec"=>"k", "lw"=>2, "alpha"=>0.5))
            end

            # cbar = fig.colorbar(cnt, ax=ax)
            # cbar.ax.set_ylabel(cbl)
        end
        # combine the colorbar but no success
        # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/colorbar_placement.html
        # fig.colorbar(cnt, ax=axes[:, col])
    end

    plt.tight_layout()

    return fig
end

fig3 = plot_fig3()

fig3.savefig("figures/Fig3.pdf")
