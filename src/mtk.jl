using ModelingToolkit

import .Utils: hill, mM, μM, ms, minute, Hz

@variables t

#=
Glucokinase (GK) Reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP (+ ADP) =>> 2G3P (+ ADP)
=#
@variables glc(t) ATP_c(t) J_GK(t)
@parameters VmaxGK = 0.011mM * Hz KatpGK = 0.5mM KglcGK = 7mM NGK = 1.7

#=
Glyceraldehydes 3-phosphate dehydrogenase (GPD), lumped with pyruvate kinase (PK)

G3P + NAD + 2ADP => Pyr + NADH + 2ATP
=#
@variables J_GPD(t) G3P(t) NAD_c(t) NADH_c(t) ADP_c(t)
@parameters VmaxGPD = 0.5mM * Hz Kg3pGPD = 0.2mM KnadGPD = 0.09 KadpGPD = 0.020mM

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

function oscillating_cac(t; ca_r=0.09μM, ka_ca=0.25μM, period=2minute, assym=5, steepness=4)
    tau = t / period
    x = assym * (tau - floor(tau))
    ca_r + ka_ca * (x * exp(1 - x))^steepness
end

@register_symbolic oscillating_cac(t)

#=
Adenylate kinase:

2ADP <=> ATP + AMP
=#
@variables J_AK(t) AMP_c(t)
@parameters kfAK = 1000Hz / mM kEqAK = 0.931

# Activity of AMPK
@variables AMPK(t)
@parameters kAMPK = 0.07

#=
Pyruvaye dehydrogenase (PDH) as well as electron transport chain (ETC)

Pyr + 4.6NAD => CO2 + 4.6NADH
=#
@variables J_PDH(t) J_DH(t) J_ETC(t) NAD_m(t) NADH_m(t) Ca_m(t)
@parameters VmaxPDH = 300μM * Hz KpyrPDH = 47.5μM KnadPDH = 81.0 U1PDH = 1.5 U2PDH = 1.1 KcaPDH = 0.05μM

function j_pdh(pyr, nad_m, nadh_m, ca_m, VMAX, K_PYR, K_NAD, U1, U2, K_CA)
    c = (hillr(ca_m, K_CA))^2
    fpCa = hillr(U2 * (1 + U1 * c))
    fNAD = hill(nad_m, nadh_m * K_NAD)
    fPyr = hill(pyr, K_PYR)
    jPDH = VMAX * fPyr * fNAD * fpCa
    return jPDH
end


#=
Electron trasnport chain (ETC)

NADH + O2 + 10H(mito) => NAD + 10H(cyto)
=#
@variables J_HR(t) J_O2(t) ΔΨm(t)
@parameters VmaxETC = 22mM * Hz KnadhETC = 3mM KaETC = -4.92E-3 / mV KbETC = -4.43E-3 / mV

# Proton leak
@variables J_HL(t)
@parameters pHleak = 2.4μM * Hz kvHleak = 0.0305 / mV

#=
F1Fo ATPase (ATP synthase) plus ANT

Reaction: 3ADP + 8H(cyto) => 3ATP + 8H(mito)
=#
@variables J_HF(t) J_ANT(t)
@parameters VmaxF1 = 8mM * Hz KadpF1 = 20μM KvF1 = 131.4mV KcaF1 = 0.165μM FmgadpF1 = 0.055

function make_model(;
    name,
    calciumEq=Ca_c ~ RestingCa + ActivatedCa * hill(ATP_c / KatpCac, ADP_c, NCac))
    eqs = [
        J_GK ~ VmaxGK * hill(ATP_c, KatpGK) * hill(glc, KglcGK),
        J_GPD ~ VmaxGPD * hill(ADP_c, KadpGPD) * hill(NAD_c / KnadGPD, NADH_c) * hill(G3P, Kg3pGPD),
        J_LDH ~ VmaxLDH * hill(Pyr, KpyrLDH) * hill(NADH_c / KnadhLDH, NAD_c),
        J_AK ~ kfAK * (ADP_c * ADP_c - ATP_c * AMP_c / kEqAK),
        AMPK ~ hill(AMP_c / ATP_c, kAMPK),
        J_PDH ~ j_pdh(Pyr, NAD_m, NADH_m, Ca_m, VmaxPDH, KpyrPDH, KnadPDH, U1PDH, U2PDH, KcaPDH),
        J_ETC ~ J_PDH,
        J_DH ~ 4.6 * J_ETC,
        J_HR ~ VmaxETC * hill(NADH_m, KnadhETC) * (1 + KaETC * ΔΨm) / (1 + KbETC * ΔΨm),
        J_O2 ~ 0.1 * J_HR,
        J_HL ~ pHleak * exp(kvHleak * ΔΨm),
        J_HF ~ VmaxF1 * hill(FmgadpF1 * ADP_c, KadpF1, 2) * hill(ΔΨm, KvF1, 8) * (1 - exp(-Ca_m / KcaF1)),
        J_ANT ~ J_HF / 3
    ]
    push!(eqs, calciumEq)
end
