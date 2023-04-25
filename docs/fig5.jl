#===
# Figure 5

Calcium oscillation
===#

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
prob = SteadyStateProblem(sys, [], [GlcConst => 10])
sssol = solve(prob, DynamicSS(Rodas5(), tspan=tend))
caavg = sssol[Ca_c]

# Calcium wave independent to ATP:ADP ratio

function cac_wave(t)
    ca_r = 0.09μM
    period = 2minute
    ka_ca = (caavg - ca_r) * 1.5
    x = 5 * ((t / period) % 1.0)
    return ca_r + ka_ca * (x * exp(1 - x))^4
end

@variables t
@register_symbolic cac_wave(t)
@named sysosci = make_model(; caceq=Ca_c~cac_wave(t))

#---
ts = range(1520.0, tend; step=2.0)
prob = ODEProblem(sysosci, sssol.u, tend, [GlcConst => 10])
sol = solve(prob, saveat=ts)

#---

function plot_fig5(sol, figsize=(10, 10))
    ts = sol.t
    tsm = ts ./ 60
    @unpack Ca_c, Ca_m, G3P, Pyr, NADH_c, NADH_m, ATP_c, ADP_c, ΔΨm, degavg = sys
    ## TODO: panel B and C not needed
    fig, ax = plt.subplots(4, 1; figsize)

    ax[1].plot(tsm, sol[Ca_c * 1000], label="Cyto. Ca (μM)")
    ax[1].plot(tsm, sol[Ca_m * 1000], label="Mito. Ca (μM)")
    ax[1].set_title("A", loc="left")

    ## TODO: describe state III behavior
    ax[2].plot(tsm, sol[ATP_c / ADP_c], label="ATP:ADP")
    ax[2].set_title("B", loc="left")

    ax[3].plot(tsm, sol[ΔΨm * 1000], label="ΔΨm (mV)")
    ax[3].set_title("C", loc="left")

    ax[4].plot(tsm, sol[degavg], label="<k>")
    ax[4].set_title("D", loc="left")
    ax[4].set(xlabel="Time (minute)")

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
fig5

# Generate Tiff file
# `fig5.savefig("Fig5.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))`
