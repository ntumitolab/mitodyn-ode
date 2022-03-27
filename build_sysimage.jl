import Pkg

Pkg.add(["PackageCompiler", "IJulia"])

using PackageCompiler
PackageCompiler.create_sysimage(
    ["DifferentialEquations"];
    project=".", sysimage_path="$(pwd())/sysimage.so")

using IJulia
IJulia.installkernel("Julia SysImage", "--sysimage=$(pwd())/sysimage.so"; specname="julia-1.7")

Pkg.rm("PackageCompiler")
Pkg.gc()
