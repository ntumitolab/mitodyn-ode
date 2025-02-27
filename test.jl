using LabelledArrays
using OrdinaryDiffEq
using SimpleUnPack: @unpack
import PythonPlot as plt

function lorenz!(D, u, p, t; f=0.0)
    @unpack σ, ρ, β = p
    @unpack x, y, z = u

    D.x = σ*(y - x)
    D.y = x*(ρ - z) - y - f
    D.z = x*y - β*z
    return nothing
end

tspan = (0.0, 100.0)
lorenz_p = LVector(σ=10.0, ρ=28.0, β=8/3)
lorenz_ic = LVector(x=1.0, y=0.0, z=0.0)
keys(lorenz_ic)
lorenz_prob = ODEProblem(ODEFunction(lorenz!; syms=keys(lorenz_ic)), lorenz_ic, tspan, lorenz_p)

lorenz_sol = solve(lorenz_prob)

ts = 0.0:0.01:100.0
xs = lorenz_sol(ts, idxs=:x).u
ys = lorenz_sol(ts, idxs=:y).u
zs = lorenz_sol(ts, idxs=:z).u
plt.figure()
plt.plot3D(xs, ys, zs)
plt.gcf()
