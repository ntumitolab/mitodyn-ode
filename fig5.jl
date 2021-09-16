using DifferentialEquations
using LabelledArrays
using Parameters

using FromFile
@from "Model/Model.jl" using Model
@from "Model/Model.jl" using Model: CaOsciSmooth
@from "Model/utils.jl" import second, μM, mV, mM, Hz

using Setfield

# Plotting config
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"


tend = 2000.0

u0 = LVector(g3p = 2.8μM,
             pyr = 8.5μM,
             nadh_c = 1μM,
             nadh_m = 60μM,
             atp_c = 4000μM,
             adp_c = 500μM,
             ca_m = 0.250μM,
             dpsi = 100mV,
             x2 = 0.2,
             x3 = 0.05)
param = setglc(MitoDynNode(), 10mM)
sssol = solve(SteadyStateProblem(model!, u0, param), DynamicSS(Rodas5(), tspan=tend))
caavg = cacyto(sssol.u[:adp_c], sssol.u[:atp_c], param, nothing)

pOsc = @set param.cai = CaOsciSmooth(ka_ca=(caavg - 0.09μM) * 1.5)

sol = solve(ODEProblem(model!, sssol.u, tend, pOsc))

# Same legend box for twin plots
# https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend

function plottwin!(ax, ts, ys, cac; title="", ylabel="", xlabel="Time (minute)", ylim=(), label="")
    ax.set_title(title, loc="left")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    isempty(ylim) || ax.set_ylim(ylim...)
    lx = ax.plot(ts, ys, "k-")
    axca = ax.twinx()
    lca = axca.plot(ts, cac, "b--")
    axca.set_ylabel("Calcium (μM)")
    axca.set_ylim(0.0, 1.0)
    ax.legend([first(lx), first(lca)], [label, "Ca (c)"], loc="upper left")
    return ax
end

function plot_fig5(sol; ts = LinRange(1520, tend, 200),
                        figsize=(10,12))

    g3p = sol.(ts, idxs=1)
    pyr = sol.(ts, idxs=2)
    nadh_c = sol.(ts, idxs=3)
    nadh_m = sol.(ts, idxs=4)
    atp_c = sol.(ts, idxs=5)
    adp_c = sol.(ts, idxs=6)
    # amp_c = ampcyto.(adp_c, atp_c, Ref(param))
    ca_m = sol.(ts, idxs=7)
    dpsi = sol.(ts, idxs=8)
    x2 = sol.(ts, idxs=9)
    x3 = sol.(ts, idxs=10)
    x1 = getx1.(x2, x3)
    ca_c = pOsc.cai.(ts)

    td = atp_c ./ adp_c
    avgDeg = avgdeg.(x2, x3, x1)

    for arr in (g3p, pyr, nadh_c, nadh_m, ca_m, ca_c, atp_c, adp_c, dpsi)
        arr .*= 1000
    end

    tsm = ts ./ 60

    fig, ax = plt.subplots(6, 1, figsize=figsize)

    ax[1].plot(tsm, [ca_c ca_m], label=["Ca (c)", "Ca (m)"])
    ax[1].set_title("A", loc="left")
    ax[1].set_ylabel("Conc. (μM)")
    # ax[1].set_xlabel("Time (minute)")
    ax[1].set_ylim(0.0, 1.0)
    ax[1].set_xlim(tsm[1], tsm[end])
    ax[1].legend(loc="right")

    ax[2].plot(tsm, [g3p nadh_c], label=["G3P", "NADH (c)"])
    ax[2].set_title("B", loc="left")
    ax[2].set_ylabel("Conc. (μM)")
    # ax[2].set_ylim()
    ax[2].set_xlim(tsm[1], tsm[end])
    ax[2].legend(loc="right")


    ax[3].plot(tsm, [pyr nadh_m], label=["Pyr", "NADH (m)"])
    ax[3].set_title("C", loc="left")
    ax[3].set_ylabel("Conc. (μM)")
    ax[3].set_xlim(tsm[1], tsm[end])
    ax[3].legend(loc="right")

    ax[4].plot(tsm, td, label="ATP:ADP")
    ax[4].set_title("D", loc="left")
    ax[4].set_ylabel("a.u.")
    ax[4].set_xlim(tsm[1], tsm[end])
    ax[4].legend(loc="right")

    ax[5].plot(tsm, dpsi, label="ΔΨ (mV)")
    ax[5].set_title("E", loc="left")
    ax[5].set_ylabel("mV")
    ax[5].set_xlim(tsm[1], tsm[end])
    ax[5].legend(loc="right")

    ax[6].plot(tsm, avgDeg, label="Average node degree")
    ax[6].set_title("F", loc="left")
    ax[6].set_ylabel("a.u.")
    ax[6].set_xlim(tsm[1], tsm[end])
    ax[6].legend(loc="right")
    ax[6].set_xlabel("Time (minute)")
    plt.tight_layout()

    fig
end

fig5 = plot_fig5(sol)

fig5.savefig("Fig5.pdf")
