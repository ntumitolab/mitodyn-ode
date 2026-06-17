# # Sensitiviy analysis
# Influence of parameters on the average degree of mitochondrial node
using DifferentialEquations
using OrdinaryDiffEqSDIRK
using SteadyStateDiffEq
using ModelingToolkit
using MitochondrialDynamics
using Tables
using MarkdownTables

# ## Setup the ODE system
@time "Build system" @named sys = make_model()
@time "Build problem" prob = SteadyStateProblem(sys, [])

params = parameters(sys)
sol0 = solve(prob, DynamicSS(FBDF()); reltol=1e-8, abstol=1e-8)
## Average degree of mitochondrial node
@unpack degavg = sys
d0 = sol0[degavg]

println("The default average degree of mitochondrial node is: ", d0)

function _calc_sens(k)
    original_value = prob.ps[k]
    _prob = remake(prob, p=[k => original_value * 1.01]) ## Increase 1% of the parameter value
    sol = solve(_prob, DynamicSS(FBDF()); reltol = 1e-8, abstol = 1e-8)
    return (sol[degavg] / d0 - 1) * 100
end

# Sensitivity analysis of the solution at `t`=300 sec against parameters.
@time sensitivities = Dict(k => _calc_sens(k) for k in params)

#---
println("Relative sensitivity of average degree of mitochondrial node to parameters:")

ks = keys(sensitivities) |> collect
vs = values(sensitivities) |> collect

t = Tables.table([ks vs]; header=["Parameter", "Relative Sensitivity"])
