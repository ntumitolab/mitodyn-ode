md"""
# Figure 5

Calcium oscillation
"""

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

tend = 2000.0
@named sys = make_model()
@unpack GlcConst, Ca_c = sys
prob = SteadyStateProblem(sys, [], [GlcConst => 10mM])
sssol = solve(prob, DynamicSS(Rodas5(), tspan=tend))
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

@variables t
@register_symbolic cac_wave(t)
@named sysosci = make_model(; cacrhs=cac_wave(t))

#---
npoints=201
ts = range(1520.0, tend, npoints)
prob = ODEProblem(sysosci, sssol.u, tend, [GlcConst => 10mM])
sol = solve(prob, saveat=ts)

#---

function plot_fig5(sol, figsize=(10, 12))
    ts = sol.t
    tsm = ts ./ 60
    @unpack Ca_c, Ca_m, G3P, Pyr, NADH_c, NADH_m, ATP_c, ADP_c, ΔΨm, degavg = sys
    ## TODO: panel B and C not needed
    fig, ax = plt.subplots(6, 1; figsize)

    ax[1].plot(tsm, sol[Ca_c * 1000], label="Ca(cyto)")
    ax[1].plot(tsm, sol[Ca_m * 1000], label="Ca(mito)")
    ax[1].set_title("A", loc="left")
    ax[1].set(ylabel="Conc. (μM)")

    ax[2].plot(tsm, sol[G3P * 1000], label="G3P")
    ax[2].plot(tsm, sol[NADH_c * 1000], label="NADH (cyto)")
    ax[2].set_title("B", loc="left")
    ax[2].set(ylabel="Conc. (μM)")

    ax[3].plot(tsm, sol[Pyr * 1000], label="Pyr")
    ax[3].plot(tsm, sol[NADH_m * 1000], label="NADH (mito)")
    ax[3].set_title("C", loc="left")
    ax[3].set(ylabel="Conc. (μM)")

    ## TODO: describe stage II behavior
    ax[4].plot(tsm, sol[ATP_c / ADP_c], label="ATP:ADP")
    ax[4].set_title("D", loc="left")

    ax[5].plot(tsm, sol[ΔΨm * 1000], label="ΔΨm")
    ax[5].set_title("E", loc="left")
    ax[5].set(ylabel="mV")

    ax[6].plot(tsm, sol[degavg], label="<k>")
    ax[6].set_title("F", loc="left")
    ax[6].set(xlabel="Time (minute)")

    for a in ax
        a.grid()
        a.legend(loc="center left")
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
