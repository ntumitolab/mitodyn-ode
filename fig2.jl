using DifferentialEquations
using LabelledArrays
using Parameters

using FromFile
@from "Model/Model.jl" using Model
@from "Model/utils.jl" import second, μM, mV, mM, Hz

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"

glc = range(3.0mM, 30.0mM, length=101)  # Range of glucose
tend = 3000.0  # Time span
u0 = LVector(g3p = 2.8μM,
             pyr = 8.5μM,
             nadh_c = 1μM,
             nadh_m = 60μM,
             atp_c = 4000μM,
             adp_c = 500μM,
             ca_m = 0.250μM,
             dpsi = 100mV,
             x2=0.25,
             x3=0.05)

param0 = MitoDynNode()

sols = map(glc) do g
    param = setglc(param0, g)
    prob = SteadyStateProblem(model!, u0, param)
    sol = solve(prob, DynamicSS(Rodas5(), tspan=tend))
end

"Plot figure 2"
function plot_fig2(sols, glc; size=(20,15),
                   glucoseScale = 5mM,
                   xlabelfontsize=16)

    # Collect data

    g3p = map(s->s.u.g3p, sols)
    pyr = map(s->s.u.pyr, sols)
    nadh_c = map(s->s.u.nadh_c, sols)
    nadh_m = map(s->s.u.nadh_m, sols)
    atp_c = map(s->s.u.atp_c, sols)
    adp_c = map(s->s.u.adp_c, sols)
    amp_c = ampcyto.(adp_c, atp_c, Ref(param0))
    ca_m = map(s->s.u.ca_m, sols)
    dpsi = map(s->s.u.dpsi, sols)
    ca_c  = cacyto.(adp_c, atp_c, Ref(param0), 0.0)
    x2 = map(s->s.u.x2, sols)
    x3 = map(s->s.u.x3, sols)
    x1 = getx1.(x2, x3)
    degree = avgdeg.(x2, x3, x1)

    for x in (g3p, pyr, nadh_c, nadh_m, atp_c, adp_c, amp_c, ca_m, dpsi, ca_c)
        x .*= 1000
    end


    # Plot
    glc5 = glc ./ glucoseScale

    fig, ax = plt.subplots(3, 3, figsize=size)

    ax[1, 1].plot(glc5, g3p)
    ax[1, 1].set(title="(A) G3P (μM)", ylim=(0.0, 10.0))

    ax[1, 2].plot(glc5, pyr)
    ax[1, 2].set(title="(B) Pyruvate (μM)", ylim=(0.0, 125.0))

    ax[1, 3].plot(glc5, [ca_c ca_m], label=["cyto", "mito"])
    ax[1, 3].set(ylim=(0.0, 1.5), title="(C) Calcium (μM)")
    ax[1, 3].legend(loc="right")

    ax[2, 1].plot(glc5, [nadh_c nadh_m], label=["cyto", "mito"])
    ax[2, 1].set(title="(D) NADH (μM)")
    ax[2, 1].legend(loc="right")

    ax[2, 2].plot(glc5, [atp_c adp_c amp_c], label=["ATP", "ADP", "AMP"])
    ax[2, 2].set(title="(E) Adenylates (μM)")
    ax[2, 2].legend(loc="right")

    ax[2, 3].plot(glc5, atp_c ./ adp_c)
    ax[2, 3].set(title="(F) ATP/ADP ratio" , ylim=(0.0, 45.0))

    ax[3, 1].plot(glc5, dpsi)
    ax[3, 1].set(title="(G) ΔΨ (mV)", ylim=(80, 160))
    ax[3, 1].set_xlabel("Glucose (X)", fontsize=xlabelfontsize)
    # ax[3, 1].legend(loc="best")

    ax[3, 2].plot(glc5, [x1 x2 x3], label=["X1", "X2", "X3"])
    ax[3, 2].set_title("(H) Mitochondrial nodes")
    ax[3, 2].set_xlabel("Glucose (X)", fontsize=xlabelfontsize)
    ax[3, 2].legend()

    ax[3, 3].plot(glc5, degree)
    ax[3, 3].set_title("(I) Average Node Degree")
    ax[3, 3].set_xlabel("Glucose (X)", fontsize=xlabelfontsize)

    plt.tight_layout()
    return fig
end

fig2 = plot_fig2(sols, glc, size=(12, 8))

fig2.savefig("figures/Fig2.pdf")
