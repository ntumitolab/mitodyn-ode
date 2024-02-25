using DifferentialEquations
using ModelingToolkit

# Define some variables
@parameters σ ρ β
@variables t x(t) y(t) z(t)
D = Differential(t)

eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]

noiseeqs = [
    0.1 * x,
    0.1 * y,
    0.1 * z
]

@named ode = ODESystem(eqs, t, [x, y, z], [σ, ρ, β])
@named sde = SDESystem(ode, noiseeqs)

u0map = [
    x => 1.0,
    y => 0.0,
    z => 0.0,
]

parammap = [
    σ => 10.0,
    β => 26.0,
    ρ => 2.33,
]

tspan = (0.0, 100.0)

odeprob = ODEProblem(ode |> complete, u0map, tspan, parammap)
sdeprob = SDEProblem(sde |> complete, u0map, tspan, parammap)

odesol = solve(odeprob, Tsit5())
sdesol = solve(sdeprob, SOSRI())

odesol(0.0:2.0:100.0, idxs=x)

sdesol(0.0:2.0:100.0, idxs=x)
