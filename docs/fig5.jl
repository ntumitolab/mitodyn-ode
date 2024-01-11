#===
# Figure 5

Calcium oscillation
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute
using PythonCall
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

# PNG output in Literate.jl
PNG(fig) = display("image/png", fig)

#---
alg = TRBDF2()
tend = 2000.0
@named sys = make_model()
@unpack GlcConst, Ca_c = sys
prob = ODEProblem(sys, [], Inf, [GlcConst => 10])
sssol = solve(prob, alg, save_everystep=false, callback=TerminateSteadyState())
caavg = sssol[Ca_c][end]

# Calcium wave independent to ATP:ADP ratio

function cac_wave(t, amplitude=1.5)
    ca_r = 0.09μM
    period = 2minute
    ka_ca = (caavg - ca_r) * amplitude
    x = 5 * ((t / period) % 1.0)
    return ca_r + ka_ca * (x * exp(1 - x))^4
end

@variables t
@register_symbolic cac_wave(t)
@named sysosci = make_model(; caceq=Ca_c~cac_wave(t))

#---
ts = range(1520.0, tend; step=2.0)
prob = ODEProblem(sysosci, [], tend, [GlcConst => 10])
sol = solve(prob, TRBDF2(), saveat=ts)

#---
function plot_fig5(sol, figsize=(10, 10))
    ts = sol.t
    tsm = ts ./ 60
    @unpack Ca_c, Ca_m, ATP_c, ADP_c, ΔΨm, degavg, J_ANT, J_HL = sys
    fig, ax = plt.subplots(5, 1; figsize)

    ax[0].plot(tsm, sol[Ca_c * 1000], label="Cyto. Ca (μM)")
    ax[0].plot(tsm, sol[Ca_m * 1000], label="Mito. Ca (μM)")
    ax[0].set_title("A", loc="left")

    ax[1].plot(tsm, sol[ATP_c / ADP_c], label="ATP:ADP")
    ax[1].set_title("B", loc="left")

    ax[2].plot(tsm, sol[ΔΨm * 1000], label="ΔΨm (mV)")
    ax[2].set_title("C", loc="left")

    ax[3].plot(tsm, sol[degavg], label="<k>")
    ax[3].set_title("D", loc="left")

    ax[4].plot(tsm, sol[J_ANT], label="ANT (mM/s)")
    ax[4].plot(tsm, sol[J_HL], label="H leak (mM/s)")
    ax[4].set_title("E", loc="left")
    ax[4].set(xlabel="Time (minute)")

    for i in 0:4
        ax[i].grid()
        ax[i].legend(loc="center left")
        ax[i].set_xlim(tsm[begin], tsm[end])
    end

    plt.tight_layout()
    return fig
end

#---

fig5 = plot_fig5(sol)

# Export figure
exportTIF(fig5, "Fig5.tif")
