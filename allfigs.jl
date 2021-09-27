import Pkg

Pkg.activate(@__DIR__)

Pkg.instantiate()

scripts = ["fig2.jl",
           "fig3.jl",
           "fig5.jl",
           "fig6.jl",
           "fig7.jl",
           "figs1.jl",
           "figs2.jl",
           "figs3.jl"]

for src in scripts
    include(joinpath(@__DIR__, src))
end
