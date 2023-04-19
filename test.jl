using ModelingToolkit
@variables t

@variables x(t) y(t)

D = Differential(t)
eqs = [
    2 * D(x) ~ 2 * x - 2 * x * y
    D(y) ~ x * y - y
]

@named sys = ODESystem(eqs, t)

equations(sys)

sys = structural_simplify(sys)

equations(sys)
