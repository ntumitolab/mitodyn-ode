# # Sensitivity analysis
# Influence of parameters on the average degree of mitochondrial node
using OrdinaryDiffEq
using ModelingToolkit
using Startup
using Tables
using MarkdownTables

# ## Setup the ODE system
@time "Build system" @named sys = make_model()
# Integrate the physical ODE to steady state. Do not use `SteadyStateProblem(sys, [])`:
# ModelingToolkit lowers it to an `SCCNonlinearProblem`, and `DynamicSS` then integrates
# the block residuals. The fission/fusion block is unstable under that pseudo-transient
# flow, which produces `retcode = Failure` and unphysical (negative) concentrations.
@time "Build problem" prob = ODEProblem(sys, [], (0.0, 1e5))

# Steady states
sts = unknowns(sys)
sol0 = solve(prob, FBDF(); reltol=1e-8, abstol=1e-8, maxiters=10^7)
@assert successful_retcode(sol0)

for st in sts
    println(st, " = ", sol0[st][end])
end

## Average degree of mitochondrial node
@unpack degavg = sys
d0 = sol0[degavg][end]

println("The default average degree of mitochondrial node is: ", d0)

function _calc_sens(k)
    original_value = prob.ps[k]
    _prob = remake(prob, p=[k => original_value * 1.01]) ## Increase 1% of the parameter value
    sol = solve(_prob, FBDF(); reltol = 1e-8, abstol = 1e-8, maxiters=10^7)
    successful_retcode(sol) || return missing
    return (sol[degavg][end] / d0 - 1) * 100
end

# Sensitivity analysis of the steady-state average degree against parameters.
params = parameters(sys)
@time sensitivities = Dict(k => _calc_sens(k) for k in params)

#---
println("Relative sensitivity of average degree of mitochondrial node to parameters:")

ks = keys(sensitivities) |> collect
vs = values(sensitivities) |> collect

t = Tables.table([ks vs]; header=["Parameter", "Relative Sensitivity"])
MarkdownTables.markdown_table(t)
