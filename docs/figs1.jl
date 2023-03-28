# # FigS1
# Response to elevated glucose concentrations in steps

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC, J_ANT, J_O2
using MitochondrialDynamics: G3P, Pyr, NADH_c, NADH_m, Ca_c, Ca_m, ΔΨm, ATP_c, ADP_c, AMP_c, degavg, t, x
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

#---

glc_step(t) = 5.0mM * (1 + (t >= 20minute) + (t >= 40minute))
@named sys = make_model(; glcrhs=glc_step(t))

tend = 60minute
prob = ODEProblem(sys, [], tend)
sol = solve(prob, tstops=[20minute, 40minute]);

#---

function plot_figs1(
    sol;
    figsize=(10, 10),
    tspan=(0.0, 60minute),
    tight=true,
    grid=true
)
    ts = sol.t
    tsm = ts ./ 60
    g3p = sol.(ts, idxs=G3P) .* 1000
    pyr = sol.(ts, idxs=Pyr) .* 1000
    nadh_c = sol.(ts, idxs=NADH_c) .* 1000
    nadh_m = sol.(ts, idxs=NADH_m) .* 1000
    ca_c = sol.(ts, idxs=Ca_c) .* 1000
    ca_m = sol.(ts, idxs=Ca_m) .* 1000
    atp_c = sol.(ts, idxs=ATP_c) .* 1000
    adp_c = sol.(ts, idxs=ADP_c) .* 1000
    amp_c = sol.(ts, idxs=AMP_c) .* 1000
    dpsi = sol.(ts, idxs=ΔΨm) .* 1000
    k = sol.(ts, idxs=degavg)
    x1 = sol.(ts, idxs=x[1])
    x2 = sol.(ts, idxs=x[2])
    x3 = sol.(ts, idxs=x[3])

    fig, ax = plt.subplots(3, 3; figsize)

    ax[1, 1].plot(tsm, g3p)
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0.0, 8.0))

    ax[1, 2].plot(tsm, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0.0, 80.0))

    ax[1, 3].plot(tsm, ca_c, label="cyto")
    ax[1, 3].plot(tsm, ca_m, label="mito")
    ax[1, 3].set(title="(C) Calcium (μM)")
    ax[1, 3].legend()

    ax[2, 1].plot(tsm, nadh_c, label="cyto")
    ax[2, 1].plot(tsm, nadh_m, label="mito")
    ax[2, 1].set(title="(D) NADH (μM)")
    ax[2, 1].legend()

    ax[2, 2].plot(tsm, atp_c, label="ATP")
    ax[2, 2].plot(tsm, adp_c, label="ADP")
    ax[2, 2].plot(tsm, amp_c, label="AMP")
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 2].legend()

    ax[2, 3].plot(tsm, atp_c ./ adp_c)
    ax[2, 3].set(title="(F) ATP:ADP", ylim=(0, 40))

    ax[3, 1].plot(tsm, dpsi)
    ax[3, 1].set(title="(G) ΔΨ (mV)", ylim=(80, 150), xlabel="Time (seconds)")

    ax[3, 2].plot(tsm, x1, label="X1")
    ax[3, 2].plot(tsm, x2, label="X2")
    ax[3, 2].plot(tsm, x3, label="X3")
    ax[3, 2].set(title="(H) Mitochondrial nodes", xlabel="Time (seconds)")
    ax[3, 2].legend()

    ax[3, 3].plot(tsm, k)
    ax[3, 3].set(title="(I) Average Node Degree", xlabel="Time (seconds)")

    if grid
        for a in ax
            a.grid()
        end
    end

    fig.set_tight_layout(tight)
    return fig
end

#---

figs1 = plot_figs1(sol)
plt.gcf()

# TIF file
figs1.savefig("FigS1.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
