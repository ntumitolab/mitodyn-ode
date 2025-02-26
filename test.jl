using ComponentArrays
using OrdinaryDiffEq
using SimpleUnPack: @unpack

tspan = (0.0, 20.0)

function lorenz!(D, u, p, t; f=0.0)
    @unpack σ, ρ, β = p
    @unpack x, y, z = u

    D.x = σ*(y - x)
    D.y = x*(ρ - z) - y - f
    D.z = x*y - β*z
    return nothing
end

lorenz_p = ComponentVector(σ=10.0, ρ=28.0, β=8/3)
lorenz_ic = ComponentVector(x=1.0, y=0.0, z=0.0)
lorenz_prob = ODEProblem(lorenz!, lorenz_ic, tspan, lorenz_p)

lorenz_sol = solve(lorenz_prob)

idxmap = Dict(keys(lorenz_ic) .=> 1:length(lorenz_ic))
lorenz_sol(range(tspan[begin], tspan[end], 101), idxs=idxmap[:y])

show(err)
