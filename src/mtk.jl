using ModelingToolkit

@variables t

#=
Adenylate kinase:

2ADP <=> ATP + AMP
=#
@variables J_ADK(t) ATP_c(t) ADP_c(t) AMP_c(t)
@parameters kfAK = 1kHz / mM
@parameters kEqAK = 0.931

#=
Glucokinase (GK) Reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP (+ ADP) =>> 2G3P (+ ADP)
=#
@variables Glc(t) J_GK(t)
@parameters VmaxGK = 0.011mM * Hz
@parameters KatpGK = 0.5mM
@parameters KglcGK = 7mM
@parameters NGK = 1.7
@parameters GlcConst = 5mM

#=
Glyceraldehydes 3-phosphate dehydrogenase (GPD), lumped with pyruvate kinase (PK) to represent lower half of glycolysis

G3P + NAD + 2ADP => Pyr + NADH + 2ATP
=#
@variables J_GPD(t) G3P(t) NAD_c(t) NADH_c(t)
@parameters VmaxGPD = 0.5mM * Hz Kg3pGPD = 0.2mM KnadGPD = 0.09 KadpGPD = 20μM

#=
Lactate dehydrogenase (LDH)
The rate of lactate output is approximately 5% of the glycolytic flux at 8 mM glucose.

Pyr + NADH => Lactate + NAD
=#
@variables J_LDH(t) Pyr(t)
@parameters VmaxLDH = 1.2mM * Hz KpyrLDH = 47.5μM KnadhLDH = 1.0

# Cytosolic calcium
@variables Ca_c(t)
@parameters RestingCa = 0.09μM ActivatedCa = 0.25μM NCac = 4 KatpCac = 25.0


# Activity of AMPK
@variables AMPKactivity(t)
@parameters kAMPK = 0.07

#=
Pyruvaye dehydrogenase (PDH) as well as electron transport chain (ETC)

Pyr + 4.6NAD => CO2 + 4.6NADH
=#
@variables J_PDH(t) J_DH(t) J_CAC(t) J_FFA(t) NAD_m(t) NADH_m(t) Ca_m(t) rPDH(t)
@parameters VmaxPDH = 300μM * Hz KpyrPDH = 47.5μM KnadPDH = 81.0 U1PDH = 1.5 U2PDH = 1.1 KcaPDH = 0.05μM


#=
Electron trasnport chain (ETC)

NADH + O2 + 10H(mito) => NAD + 10H(cyto)
=#
@variables J_HR(t) J_O2(t) ΔΨm(t) rETC(t)
@parameters VmaxETC = 22mM * Hz KnadhETC = 3mM KaETC = -4.92E-3 / mV KbETC = -4.43E-3 / mV

# Proton leak
@variables J_HL(t) rHL(t)
@parameters pHleak = 2.4μM * Hz kvHleak = 0.0305 / mV

#=
F1Fo ATPase (ATP synthase) plus ANT

3ADP + 8H(cyto) => 3ATP + 8H(mito)
=#
@variables J_HF(t) J_ANT(t) rF1(t)
@parameters VmaxF1 = 8mM * Hz KadpF1 = 20μM KvF1 = 131.4mV KcaF1 = 0.165μM FmgadpF1 = 0.055

#=
Mitochondrial calcium uniporter (MCU)

Ca(cyto) <=> Ca(mito)
=#
@variables J_MCU(t)
@parameters PcaMCU = 4Hz


#=
Mitochondrial sodium calcium exchanger (NCLX)

Ca(mito) + 2Na(cyto) <=> Ca(cyto) + 2Na(mito)
=#
@variables J_NCLX(t)
@parameters Na_c = 10mM Na_m = 5mM VmaxNCLX = 75μM * Hz KnaNCLX = 8.2mM KcaNCLX = 8μM

#=
NADH shuttle

NADH(cyto) + NAD(mito) => NADH(mito) + NAD(cyto)
=#
@variables J_NADHT(t)
@parameters VmaxNADHT = 50μM * Hz Ktn_c = 0.002 Ktn_m = 16.78

# Fission-fusion rates
@variables (x(t))[1:3] x13r(t) degavg(t) (v(t))[1:2]
@parameters Kfiss1 = inv(10minute) Kfuse1 = Kfiss1 Kfiss2 = 1.5Kfiss1 Kfuse2 = 0.5Kfuse1

# Baseline consumption rates
@parameters kNADHc = 0.1Hz kNADHm = 0.1Hz kATP = 0.04Hz kATPCa = 90Hz / mM kG3P = 0.01Hz kPyr = 0.01Hz

# Conservation relationships
@parameters ΣAc = 4.5mM Σn_c = 2mM Σn_m = 2.2mM


