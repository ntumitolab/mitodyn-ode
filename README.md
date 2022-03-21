# Mitochondrial dynamics model

Website: <https://ntumitolab.github.io/mitodyn-ode/>

## How to run the code in the cloud

Click `run in binder` on the top of web pages.

## How to run the code locally

### Install Git

Git can be downloaded and installed [here](https://git-scm.com/download/).

### Install Julia

Download and install Julia (*version 1.7+*) from [the official website](https://julialang.org/downloads/). Choose "add Julia to system PATH" during installation.

### Download the code and launch Julia

Open Git Bash and enter the following commands:

```bash
git clone https://github.com/NTUMitoLab/MitochondrialDynamics.git
cd MitochondrialDynamics
julia --project=@.
```

If that does not work, download this repositroy via the big green `Code` button => `download zip`. You will see the project folder `MitochondrialDynamics` in the extracted zip archive.

### Running the code

In the Julia terminal, type the following commands to run the code to generate figures and data. Make sure the working directory (output of `pwd()`) is the project folder `MitochondrialDynamics`.

```jl
using Pkg
Pkg.add("IJulia")
Pkg.activate(".")
Pkg.instantiate()

using IJulia
IJulia.notebook(dir=pwd())
```
