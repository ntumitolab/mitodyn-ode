# # Fig 7-1 Oxygen consumption in response to glucose stimulation and chemical agents.

# Load packages
using DifferentialEquations
using LabelledArrays
using Parameters
using FromFile
@from "Model/Model.jl" using Model
@from "Model/utils.jl" import second, μM, mV, mM, Hz
using Setfield

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"

tend = 80minute

u0 = LVector(g3p = 2.8μM,
             pyr = 8.5μM,
             nadh_c = 1μM,
             nadh_m = 60μM,
             atp_c = 4000μM,
             adp_c = 500μM,
             ca_m = 0.250μM,
             dpsi = 100mV,
             x2 = 0.20,
             x3 = 0.05)

param = MitoDynNode()

paramDM = let rPDH = 0.5, rETC = 0.75, rHL  = 1.4, rF1  = 0.5
    p = @set param.pdh.VMAX *= rPDH
    p = @set p.hleak.P_H *= rHL
    p = @set p.f1fo.VMAX *= rF1
    p = @set p.etc.VMAX *= rETC
end


# Add glucose to 20mM
function add_glc!(integrator, glc=20mM, newdt = 0.01)
    integrator.p = setglc(integrator.p, glc)
    set_proposed_dt!(integrator, newdt)
end

# Oligomycin: cut ATP synthase to 10%
function oligomycin!(integrator, f1Cut = 0.05, newdt = 0.01)
    p = integrator.p
    integrator.p = @set p.f1fo.VMAX *= f1Cut
    set_proposed_dt!(integrator, newdt)
end

# Rotenone: cut ETC
function rotenone!(integrator, etcCut = 0.05, newdt = 0.01)
    p = integrator.p
    integrator.p = @set p.etc.VMAX *= etcCut
    set_proposed_dt!(integrator, newdt)
end

# FCCP: increase proton leak
function fccp!(integrator, leakRate = 10.0, newdt = 0.01)
    p = integrator.p
    integrator.p = @set p.hleak.P_H *= leakRate
    set_proposed_dt!(integrator, newdt)
end

ssalg=DynamicSS(Rodas5(), tspan=tend)

u∞ = solve(SteadyStateProblem(model!, u0, param), ssalg).u
u∞DM = solve(SteadyStateProblem(model!, u0, paramDM), ssalg).u

cbs = CallbackSet( 	PresetTimeCallback([20minute], add_glc!),
                    PresetTimeCallback([40minute], oligomycin!),
                    PresetTimeCallback([60minute], fccp!))

sol = solve(ODEProblem(model!, u∞, tend, param), callback=cbs)
solDM = solve(ODEProblem(model!, u∞DM, tend, paramDM), callback=cbs)

function plot_fig71(sol, solDM, tend = 80minute;
    figsize=(20,20), lables=["Baseline","Diabetic"])

    ts = range(0.0, tend, length=301)

    g3p = [sol.(ts, idxs=1) solDM.(ts, idxs=1)]
    pyr = [sol.(ts, idxs=2) solDM.(ts, idxs=2)]
    nadh_c = [sol.(ts, idxs=3) solDM.(ts, idxs=3)]
    nadh_m = [sol.(ts, idxs=4) solDM.(ts, idxs=4)]
    atp_c = [sol.(ts, idxs=5) solDM.(ts, idxs=5)]
    adp_c = [sol.(ts, idxs=6) solDM.(ts, idxs=6)]
    td = atp_c ./ adp_c
    ca_m = [sol.(ts, idxs=7) solDM.(ts, idxs=7)]
    ca_c = cacyto.(adp_c, atp_c, Ref(param), nothing)
    dpsi = [sol.(ts, idxs=8) solDM.(ts, idxs=8)]
    x2 = [sol.(ts, idxs=9) solDM.(ts, idxs=9)]
    x3 = [sol.(ts, idxs=10) solDM.(ts, idxs=10)]
    x1 = getx1.(x2, x3)
    avgDeg = avgdeg.(x2, x3, x1)

    for arr in (g3p, pyr, nadh_c, nadh_m, atp_c, adp_c, ca_m, ca_c, dpsi)
        arr .*= 1000
    end

    ll = lables
    tsm = ts ./ 60
    fig, ax = plt.subplots(3, 3, figsize=figsize)

    ax[1].plot(tsm, g3p)
    ax[1].legend(ll)
    ax[1].set_ylabel("G3P (μM)")
    ax[1].set_xlabel("Time (minute)")
    ax[1].set_title("A", loc="left")
    ax[1].set_ylim(0.0, 9.0)

    ax[4].plot(tsm, pyr)
    ax[4].legend(ll)
    ax[4].set_ylabel("Pyruvate (μM)")
    ax[4].set_xlabel("Time (minute)")
    ax[4].set_title("B", loc="left")
    ax[4].set_ylim(0.0, 700.0)

    ax[7].plot(tsm, nadh_c)
    ax[7].legend(ll)
    ax[7].set_ylabel("Cytosolic NADH (μM)")
    ax[7].set_xlabel("Time (minute)")
    ax[7].set_title("C", loc="left")
    ax[7].set_ylim(0.0, 8.0)

    ax[2].plot(tsm, nadh_m)
    ax[2].legend(ll)
    ax[2].set_ylabel("Mitochondrial NADH (μM)")
    ax[2].set_xlabel("Time (minute)")
    ax[2].set_title("D", loc="left")
    ax[2].set_ylim(0.0, 200.0)

    ax[5].plot(tsm, ca_c)
    ax[5].legend(ll)
    ax[5].set_ylabel("Cytosolic Calcium (μM)")
    ax[5].set_xlabel("Time (minute)")
    ax[5].set_title("E", loc="left")
    ax[5].set_ylim(0.0, 0.4)

    ax[8].plot(tsm, ca_m)
    ax[8].legend(ll)
    ax[8].set_ylabel("Mitochondrial Calcium (μM)")
    ax[8].set_xlabel("Time (minute)")
    ax[8].set_title("F", loc="left")
    ax[8].set_ylim(0.0, 1.2)

    ax[3].plot(tsm, td)
    ax[3].legend(ll)
    ax[3].set_ylabel("ATP/ADP ratio")
    ax[3].set_xlabel("Time (minute)")
    ax[3].set_title("G", loc="left")
    ax[3].set_ylim()

    ax[6].plot(tsm, dpsi)
    ax[6].legend(ll, loc="upper right")
    ax[6].set_ylabel("Mitochodrial membrane potential (mV)")
    ax[6].set_xlabel("Time (minute)")
    ax[6].set_title("H", loc="left")
    ax[6].set_ylim()

    ax[9].plot(tsm, avgDeg)
    ax[9].legend(ll)
    ax[9].set_ylabel("Average Node Degree")
    ax[9].set_xlabel("Time (minute)")
    ax[9].set_title("I", loc="left")
    ax[9].set_ylim()

    fig
end

fig71 = plot_fig71(sol, solDM, tend, figsize=(20, 16))

fig71.savefig("figures/S1Fig3-1.pdf")

# Fig 7-2

ssalg=DynamicSS(Rodas5(), tspan=tend)

u∞ = solve(SteadyStateProblem(model!, u0, param), ssalg).u
u∞DM = solve(SteadyStateProblem(model!, u0, paramDM), ssalg).u

cbs = CallbackSet( 	PresetTimeCallback([20minute], add_glc!),
                    PresetTimeCallback([40minute], oligomycin!),
                    PresetTimeCallback([60minute], rotenone!))

sol = solve(ODEProblem(model!, u∞, tend, param), callback=cbs)
solDM = solve(ODEProblem(model!, u∞DM, tend, paramDM), callback=cbs)

function plot_fig72(sol, solDM, tend = 80minute;
    figsize=(20,20), lables=["Baseline","Diabetic"], etcCut = 0.05)

    ts = range(0.0, tend, length=301)

    g3p = [sol.(ts, idxs=1) solDM.(ts, idxs=1)]
    pyr = [sol.(ts, idxs=2) solDM.(ts, idxs=2)]
    nadh_c = [sol.(ts, idxs=3) solDM.(ts, idxs=3)]
    nadh_m = [sol.(ts, idxs=4) solDM.(ts, idxs=4)]
    atp_c = [sol.(ts, idxs=5) solDM.(ts, idxs=5)]
    adp_c = [sol.(ts, idxs=6) solDM.(ts, idxs=6)]
    td = atp_c ./ adp_c
    ca_m = [sol.(ts, idxs=7) solDM.(ts, idxs=7)]
    ca_c = cacyto.(adp_c, atp_c, Ref(param), nothing)
    dpsi = [sol.(ts, idxs=8) solDM.(ts, idxs=8)]
    x2 = [sol.(ts, idxs=9) solDM.(ts, idxs=9)]
    x3 = [sol.(ts, idxs=10) solDM.(ts, idxs=10)]
    x1 = getx1.(x2, x3)
    avgDeg = avgdeg.(x2, x3, x1)
    x13 = x1 ./ x3

    jo2 = [param.etc.(nadh_m[:, 1], dpsi[:, 1]) .* 1//10  paramDM.etc.(nadh_m[:, 2], dpsi[:, 2]) .* 1//10]
    jo2[ts .>= 60minute, :] .*= etcCut

    for arr in (g3p, pyr, nadh_c, nadh_m, atp_c, adp_c, ca_m, ca_c, dpsi)
        arr .*= 1000
    end

    ll = lables

    tsm = ts ./ 60

    fig, ax = plt.subplots(3, 3, figsize=figsize)

    ax[1].plot(tsm, jo2)
    ax[1].legend(ll)
    ax[1].set_ylabel("O2 consumption")
    ax[1].set_xlabel("Time (minute)")
    ax[1].set_title("A", loc="left")
    ax[1].set_ylim()

    ax[4].plot(tsm, pyr)
    ax[4].legend(ll)
    ax[4].set_ylabel("Pyruvate (μM)")
    ax[4].set_xlabel("Time (minute)")
    ax[4].set_title("B", loc="left")
    ax[4].set_ylim(0.0, 700.0)

    ax[7].plot(tsm, nadh_c)
    ax[7].legend(ll)
    ax[7].set_ylabel("Cytosolic NADH (μM)")
    ax[7].set_xlabel("Time (minute)")
    ax[7].set_title("C", loc="left")
    ax[7].set_ylim(0.0, 20.0)

    ax[2].plot(tsm, nadh_m)
    ax[2].legend(ll)
    ax[2].set_ylabel("Mitochondrial NADH (μM)")
    ax[2].set_xlabel("Time (minute)")
    ax[2].set_title("D", loc="left")
    ax[2].set_ylim()

    ax[5].plot(tsm, ca_c)
    ax[5].legend(ll)
    ax[5].set_ylabel("Cytosolic Calcium (μM)")
    ax[5].set_xlabel("Time (minute)")
    ax[5].set_title("E", loc="left")
    ax[5].set_ylim(0.0, 0.4)

    ax[8].plot(tsm, ca_m)
    ax[8].legend(ll)
    ax[8].set_ylabel("Mitochondrial Calcium (μM)")
    ax[8].set_xlabel("Time (minute)")
    ax[8].set_title("F", loc="left")
    ax[8].set_ylim(0.0, 1.2)

    ax[3].plot(tsm, td)
    ax[3].legend(ll)
    ax[3].set_ylabel("ATP/ADP ratio")
    ax[3].set_xlabel("Time (minute)")
    ax[3].set_title("G", loc="left")
    ax[3].set_ylim()

    ax[6].plot(tsm, dpsi)
    ax[6].legend(ll, loc="upper right")
    ax[6].set_ylabel("Mitochodrial membrane potential (mV)")
    ax[6].set_xlabel("Time (minute)")
    ax[6].set_title("H", loc="left")
    ax[6].set_ylim()

    ax[9].plot(tsm, avgDeg)
    ax[9].legend(ll)
    ax[9].set_ylabel("Average Node Degree")
    ax[9].set_xlabel("Time (minute)")
    ax[9].set_title("I", loc="left")
    ax[9].set_ylim()

    return fig
end

fig72 = plot_fig72(sol, solDM, tend, figsize=(20, 16))

fig72.savefig("figures/S1Fig3-2.pdf")
