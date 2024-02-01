using PrettyTables, Literate, Pkg
ENV["GKSwstype"] = "100"
ENV["JULIA_DEBUG"] = "Literate"
Pkg.activate(Base.current_project())
basedir = "docs"
nbs = String[]

# Collect the list of Literate notebooks (ends with .jl)
for (root, dirs, files) in walkdir(basedir)
    for file in files
        if (endswith(file, ".jl"))
            push!(nbs, joinpath(root, file))
        end
    end
end

# Execute the notebooks in worker process(es)
ts = map(nbs) do nb
    try
        @elapsed Literate.notebook(nb, dirname(nb); mdstrings=true)
    catch e
        println(e)
        NaN
    end
end

pretty_table([nbs ts], header=["Notebook", "Elapsed (s)"])
any(isnan, ts) && error("Please check errors.")
