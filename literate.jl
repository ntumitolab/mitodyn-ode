using Pkg
using Literate
using PrettyTables
using SHA

ENV["GKSwstype"] = "100"
ENV["JULIA_DEBUG"] = "Literate"
Pkg.activate(Base.current_project())

basedir = "docs"
nbs = String[]

# Collect the list of Literate notebooks (ends with .jl)
for (root, dirs, files) in walkdir(basedir)
    for file in files
        if endswith(file, ".jl")
            nb = joinpath(root, file)
            shaval = read(nb, String) |> sha256 |> bytes2hex
            @info "$(nb): SHA256=$(shaval)"
            shafilename = splitext(nb)[1] * ".sha"
            # Cache hit
            if isfile(shafilename) && read(shafilename, String) == shaval
                @info "Notebook $(nb) cache hits and will not be executed."
            # Cache miss
            else
                write(shafilename, shaval)
                push!(nbs, nb)
            end
        # Remove notebook and sha file if the respective literate notebook does not exist
        elseif endswith(file, ".ipynb") || endswith(file, ".sha")
            filename = joinpath(root, file)
            nb = splitext(filename)[1] * ".jl"
            if !isfile(nb)
                rm(filename)
            end
        end
    end
end

# Execute notebooks
ts = map(nbs) do nb
    try
        @elapsed Literate.notebook(nb, dirname(nb); mdstrings=true)
    catch err
        println(err)
        NaN
    end
end

pretty_table([nbs ts], header=["Notebook", "Elapsed (s)"])
any(isnan, ts) && error("Please check errors.")
