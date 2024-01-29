using Literate
ENV["JULIA_DEBUG"] = "Literate"
nb = ARGS[1]
Literate.notebook(nb, dirname(nb); mdstrings=true)
