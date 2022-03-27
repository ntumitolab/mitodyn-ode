import Pkg

# Do not install default kernel
ENV["IJULIA_NODEFAULTKERNEL"] = 1

Pkg.add(["PackageCompiler", "IJulia"])

using PackageCompiler
PackageCompiler.create_sysimage(
    ["DifferentialEquations"];
    project=".", sysimage_path="$(pwd())/sysimage.so")

using IJulia
IJulia.installkernel("Julia SysImage", "--sysimage=$(pwd())/sysimage.so"; specname="julia-1.7")

Pkg.rm("PackageCompiler")
Pkg.gc()
