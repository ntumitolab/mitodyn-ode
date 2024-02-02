using OrdinaryDiffEq
using SteadyStateDiffEq
using ModelingToolkit

@variables t x(t) RHS(t)
@parameters τ
D = Differential(t)

eqs = [
    RHS ~ (1 - x) / τ
    D(x) ~ RHS
]

@named sys = ODESystem(eqs, defaults=[x=>0.0, τ=>3.0])

sys = structural_simplify(sys)

prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
sol[sys.x]
sol[sys.RHS]

ssprob = SteadyStateProblem(sys, [])
sssol = solve(ssprob, DynamicSS(Rodas5()))
sssol[sys.x]
sssol[sys.RHS] # MethodError: No method matching... (details below)
