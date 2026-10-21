# Packages

Direct dependencies of the root `Project.toml` (Julia 1.13.2, resolved in `Manifest.toml`).

| Package | Version | Role in this project |
|---|---|---|
| DiffEqCallbacks | 4.19.5 | Callbacks for ODE solves (events, interruptions) |
| JSON | 1.11.0 | Processing Jupyter notebooks |
| Literate | 2.21.0 | `docs/*.jl` Literate sources → notebooks |
| MarkdownTables | 1.1.0 | Parameter tables in the docs |
| ModelingToolkit | 11.45.5 | ODE system modeling framework |
| OrdinaryDiffEq | 7.8.1 | ODE solvers |
| OrdinaryDiffEqSDIRK | 2.9.8 | SDIRK solvers for stiff/DAE problems |
| PythonPlot | 1.0.8 | Plotting via matplotlib (PythonCall) |
| SHA | 1.0.0 | Hashing (used in CI) |
| Startup | 0.1.0 | Local workspace package: `make_model()`, utils |
| SteadyStateDiffEq | 2.18.1 | Steady-state solving |
| Tables | 1.14.0 | Tables interface for results |
