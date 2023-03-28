md"""
# Figure 5

Calcium oscillation
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: GlcConst, degavg, ΔΨm, Ca_c, Ca_m, t, G3P
using MitochondrialDynamics: Pyr, NADH_c, NADH_m, ATP_c, ADP_c, degavg
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

tend = 2000.0
@named sys = make_model(; glcrhs=10mM)
sssol = solve(SteadyStateProblem(sys, []), DynamicSS(Rodas5(), tspan=tend))
caavg = sssol[Ca_c]

#---

function cac_wave(t)
    ca_r = 0.09μM
    period = 2minute
    ka_ca = (caavg - ca_r) * 1.5
    x = 5 * ((t / period) % 1.0)
    return ca_r + ka_ca * (x * exp(1 - x))^4
end

#---

@register_symbolic cac_wave(t)
@named sysosci = make_model(; cacrhs=cac_wave(t), glcrhs=10mM)

#---

prob = ODEProblem(sysosci, sssol.u, tend)
sol = solve(prob)

#---

function plot_fig5(
    sol;
    tspan=(1520.0, 2000.0),
    npoints=200,
    figsize=(10, 12)
)
    ts = LinRange(tspan[1], tspan[2], npoints)
    tsm = ts ./ 60
    ca_c = sol(ts, idxs=Ca_c) .* 1000
    ca_m = sol(ts, idxs=Ca_m) .* 1000
    g3p = sol(ts, idxs=G3P) .* 1000
    pyr = sol(ts, idxs=Pyr) .* 1000
    nadh_c = sol(ts, idxs=NADH_c) .* 1000
    nadh_m = sol(ts, idxs=NADH_m) .* 1000
    atp_c = sol(ts, idxs=ATP_c) .* 1000
    adp_c = sol(ts, idxs=ADP_c) .* 1000
    td = atp_c.u ./ adp_c
    dpsi = sol(ts, idxs=ΔΨm) .* 1000
    k = sol(ts, idxs=degavg)

    fig, ax = plt.subplots(6, 1; figsize)

    ax[1].plot(tsm, ca_c, label="Ca(cyto)")
    ax[1].plot(tsm, ca_m, label="Ca(mito)")
    ax[1].set_title("A", loc="left")
    ax[1].set(ylabel="Conc. (μM)")

    ax[2].plot(tsm, g3p, label="G3P")
    ax[2].plot(tsm, nadh_c, label="NADH (cyto)")
    ax[2].set_title("B", loc="left")
    ax[2].set(ylabel="Conc. (μM)")

    ax[3].plot(tsm, pyr, label="Pyr")
    ax[3].plot(tsm, nadh_m, label="NADH (mito)")
    ax[3].set_title("C", loc="left")
    ax[3].set(ylabel="Conc. (μM)")

    ax[4].plot(tsm, td, label="ATP:ADP")
    ax[4].set_title("D", loc="left")

    ax[5].plot(tsm, dpsi, label="ΔΨm")
    ax[5].set_title("E", loc="left")
    ax[5].set(ylabel="mV")

    ax[6].plot(tsm, k, label="<k>")
    ax[6].set_title("F", loc="left")
    ax[6].set(xlabel="Time (minute)")

    for a in ax
        a.grid()
        a.legend()
        a.set_xlim(tsm[begin], tsm[end])
    end

    plt.tight_layout()
    return fig
end

#---

fig5 = plot_fig5(sol)
plt.gcf()

# Generate Tiff file

fig5.savefig("Fig5.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
