#===
# Supplimentary figures
## Figure S3

Changes in response to both glucose stimulation and chemical agents.

Using the Glucose-Oligomycin-FCCP protocol.
===#

using DifferentialEquations
using ModelingToolkit
using MitochondrialDynamics
using MitochondrialDynamics: second, μM, mV, mM, Hz, minute

import PyPlot as plt
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["font.size"] = 14
## rcParams["font.sans-serif"] = "Arial"
## rcParams["font.family"] = "sans-serif"

@named sys = make_model()

@unpack GlcConst, VmaxPDH, pHleak, VmaxF1, VmaxETC = sys

idxGlc = findfirst(isequal(GlcConst), parameters(sys))
idxVmaxPDH = findfirst(isequal(VmaxPDH), parameters(sys))
idxpHleak = findfirst(isequal(pHleak), parameters(sys))
idxVmaxF1 =  findfirst(isequal(VmaxF1), parameters(sys))
idxVmaxETC =  findfirst(isequal(VmaxETC), parameters(sys))

tend = 80minute
ts = range(0, tend, 401)

prob = ODEProblem(sys, [], ts[end])
probs5 = ODEProblem(sys, [], ts[end])

function remake_dm(prob; rPDH=0.5, rETC=0.75, rHL=1.4, rF1=0.5)
    p = copy(prob.p)
    p[idxVmaxETC] *= rETC
    p[idxVmaxF1] *= rF1
    p[idxpHleak] *= rHL
    p[idxVmaxPDH] *= rPDH
    return remake(prob, p=p)
end

prob_dm = remake_dm(prob)
prob_dmS5 = remake_dm(probs5)

# Define events
function add_glucose!(i)
    i.p[idxGlc] += 15mM
    set_proposed_dt!(i, 0.1)
end

add_glucose_cb = PresetTimeCallback(20minute, add_glucose!)

function add_oligomycin!(i)
    i.p[idxVmaxF1] *= 0.05
    set_proposed_dt!(i, 0.1)
end

add_oligomycin_cb = PresetTimeCallback(40minute, add_oligomycin!)

function add_rotenone!(i)
    i.p[idxVmaxETC] *= 0.05
    set_proposed_dt!(i, 0.1)
end

add_rotenone_cb = PresetTimeCallback(60minute, add_rotenone!)

function add_fccp!(i)
    i.p[idxpHleak] *= 10
    set_proposed_dt!(i, 0.1)
end

add_fccp_cb = PresetTimeCallback(60minute, add_fccp!)

cbs = CallbackSet(add_glucose_cb, add_oligomycin_cb, add_fccp_cb)
sols3 = solve(prob; callback=cbs, saveat=ts)
solDMs3 = solve(prob_dm; callback=cbs, saveat=ts)

figs3 = plot_figs2(sols3, solDMs3)
figs3

# TIFF file
figs3.savefig("FigS3-Glucose-Oligomycin-FCCP.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))

# ## Figure S5
# Baseline vs. Diabetic models using Glucose-Oligomycin-Rotenone protocol.
cbs = CallbackSet(add_glucose_cb, add_oligomycin_cb, add_rotenone_cb)
sols5 = solve(probs5; callback=cbs, saveat=ts)
sols5DM = solve(prob_dmS5; callback=cbs, saveat=ts)
figs5 = plot_figs2(sols5, sols5DM)
figs5

# TIFF file
figs5.savefig("FigS5-Glucose-Oligomycin-Rotenone.tif", dpi=300, format="tiff", pil_kwargs=Dict("compression" => "tiff_lzw"))
