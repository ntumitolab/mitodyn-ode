#===
# Figure 6: Calcium oscillation
===#
using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute, hil
import PythonPlot as plt
plt.matplotlib.rcParams["font.size"] = 14

#---
@named sys = make_model()
@unpack Ca_c, GlcConst, kATPCa, kATP = sys
alg = TRBDF2()

# Calcium oscillation function
function cac_wave(; ca_base = 0.09μM, ca_act = 0.25μM, n=4, katp=25, amplitude=0.5, period=2minute)
    @variables t Ca_c(t) ATP_c(t) ADP_c(t)
    @parameters (RestingCa=ca_base, ActivatedCa=ca_act, NCac=n, KatpCac=katp)
    x = 5 * ((t / period) % 1.0) ## An oscillating function
    w = (x * exp(1 - x))^4  ## Scale from 0 to 1
    caceq = Ca_c ~ RestingCa + ActivatedCa * hil(ATP_c, KatpCac * ADP_c, NCac) * (1 + amplitude * (2w-1))
    return caceq
end

@named sysosci = make_model(; caceq=cac_wave(amplitude=0.8))

#---
tend = 4000.0
ts = range(tend-480, tend; length=201)
prob = ODEProblem(sysosci, [], tend, [GlcConst => 10mM])
sol = solve(prob, alg, saveat=ts)

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

    ax[3].plot(tsm, sol[degavg], label="Average node degree")
    ax[3].set_title("D", loc="left")

    ax[4].plot(tsm, sol[J_ANT], label="ATP export (mM/s)")
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
exportTIF(fig5, "Fig6-ca-oscillation.tif")

# Tuning ca-dependent ATP consumption rate (kATPCa)
# kATPCa : 90 -> 10
prob2 = ODEProblem(sysosci, [], tend, [GlcConst => 10mM, kATPCa=>10Hz/mM, kATP=>0.055Hz])
sol2 = solve(prob2, alg, saveat=ts)
plot_fig5(sol2)

# kATPCa : 90 -> 0.1
prob4 = ODEProblem(sysosci, [], tend, [GlcConst => 10mM, kATPCa=>0.1Hz/mM, kATP=>0.06Hz])
sol4 = solve(prob4, alg, saveat=ts)
plot_fig5(sol4)
