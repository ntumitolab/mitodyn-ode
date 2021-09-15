# Supplementary figure 1: Response to elevated glucose concentrations in steps.
using DifferentialEquations
using LabelledArrays
using Parameters
using FromFile
@from "Model/Model.jl" using Model
@from "Model/utils.jl" import second, μM, mV, mM, Hz

# Plotting
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
             x2=0.25, x3=0.05)
tend = 3000.0
param = MitoDynNode()
prob = ODEProblem(model!, u0, tend, param)

sol = solve(prob)
uInf = sol.u[end]

@show getx1(uInf) avgdeg(uInf)

"Create a function to change glucose level for callbacks"
function make_glc_event(glc, newdt=0.01)
    return i -> begin
        i.p = setglc(i.p, glc)
        set_proposed_dt!(i, newdt)
    end
end

cb1 = PresetTimeCallback([1000.0second], make_glc_event(10mM))
cb2 = PresetTimeCallback([2000.0second], make_glc_event(15mM))
cbs = CallbackSet(cb1, cb2)

sol1s = solve(ODEProblem(model!, uInf, tend, param), callback=cbs)

function plot_fig1s(sol1s, tend; figsize=(15,15))
    ts = LinRange(0.0, tend, 201)
    g3p = sol1s.(ts, idxs=1)
    pyr = sol1s.(ts, idxs=2)
    nadh_c = sol1s.(ts, idxs=3)
    nadh_m = sol1s.(ts, idxs=4)
    atp_c = sol1s.(ts, idxs=5)
    adp_c = sol1s.(ts, idxs=6)
    amp_c = ampcyto.(adp_c, atp_c, Ref(param))
    ca_m = sol1s.(ts, idxs=7)
    ca_c = cacyto.(adp_c, atp_c, Ref(param), ts)
    ΔΨ = sol1s.(ts, idxs=8)
    x2 = sol1s.(ts, idxs=9)
    x3 = sol1s.(ts, idxs=10)
    x1 = Model.getx1.(x2, x3)
    degree = Model.avgdeg.(x2, x3, x1)

    # convert unit to μM
    for arr in (g3p, pyr, nadh_c, nadh_m, atp_c, adp_c, amp_c, ca_m, ca_c, ΔΨ)
        arr .*= 1000
    end

    fig, ax = plt.subplots(3, 3, figsize=figsize)
    ax[1, 1].plot(ts, g3p, label=nothing)
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0.0, 10.0))

    ax[1, 2].plot(ts, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0.0, 80.0))

    ax[1, 3].plot(ts, [ca_c ca_m], label=["cyto", "mito"])
    ax[1, 3].set(ylim=(0.0, 1.5), title="(C) Calcium (μM)")
    ax[1, 3].legend(loc="upper left")

    ax[2, 1].plot(ts, [nadh_c nadh_m], label=["cyto", "mito"])
    ax[2, 1].set(title="(D) NADH (μM)")
    ax[2, 1].legend()

    ax[2, 2].plot(ts, [atp_c adp_c amp_c], label=["ATP", "ADP", "AMP"])
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 2].legend()

    ax[2, 3].plot(ts, atp_c ./ adp_c)
    ax[2, 3].set(title="(F) ATP/ADP ratio", ylim=(0.0, 45.0) )

    ax[3, 1].plot(ts, ΔΨ)
    ax[3, 1].set(title="(G) ΔΨ (mV)", ylim=(80, 150), xlabel="Time (seconds)")

    ax[3, 2].plot(ts, [x1 x2 x3], label=["X1", "X2", "X3"])
    ax[3, 2].set(title="(H) Mitochondrial nodes", xlabel="Time (seconds)", ylim=(0.0, 0.6))
    ax[3, 2].legend()

    ax[3, 3].plot(ts, degree)
    ax[3, 3].set(title="(I) Average Node Degree", xlabel="Time (seconds)")

    plt.tight_layout()
    return fig
end

fig1s = plot_fig1s(sol1s, tend, figsize=(15,10))

fig1s.savefig("S1_Fig1.pdf")
