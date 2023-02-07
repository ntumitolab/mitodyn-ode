using ModelingToolkit

@variables t

"""
Rate of Adenylate kinase (AdK):

2ADP <=> ATP + AMP
"""
function j_adk(atp_c, adp_c, amp_c, kf, keq)
    return kf * (adp_c * adp_c - atp_c * amp_c * keq)
end

@variables J_ADK(t) ATP_c(t) ADP_c(t) AMP_c(t)
@parameters kfAK = 1kHz / mM
@parameters kEqAK = 0.931

"""
Rate of Glucokinase (GK) and upper glycolysis reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP (+ ADP) =>> 2G3P (+ ADP)
"""
function j_gk(atp_c, glc, vmax, katp, kglc, n)
    return vmax * hill(atp_c, katp) * hill(glc, kglc, n)
end

@variables Glc(t) J_GK(t)
@parameters VmaxGK = 0.011mM * Hz
@parameters KatpGK = 0.5mM
@parameters KglcGK = 7mM
@parameters NGK = 1.7
@parameters GlcConst = 5mM

"""
Rate of Glyceraldehydes 3-phosphate dehydrogenase (GPD), lumped with pyruvate kinase (PK) to represent lower half of glycolysis

G3P + NAD + 2ADP => Pyr + NADH + 2ATP
"""
function j_gpd(g3p, nad_c, nadh_c, adp_c, vmax, kg3p, knad, kadp)
    return vmax * hill(adp_c, kadp) * hill(nad_c, nadh_c * knad) * hill(g3p, kg3p)
end

@variables J_GPD(t) G3P(t) NAD_c(t) NADH_c(t)
@parameters VmaxGPD = 0.5mM * Hz
@parameters Kg3pGPD = 0.2mM
@parameters KnadGPD = 0.09
@parameters KadpGPD = 20μM

"""
Rate of Lactate dehydrogenase (LDH)
The rate of lactate output is approximately 5% of the glycolytic flux @ glucose = 8mM.

Pyr + NADH => Lactate + NAD
"""
function j_ldh(pyr, nad_c, nadh_c, vmax, kpyr, knadh)
    return vmax * hill(pyr, kpyr) * hill(nadh_c, nad_c * knadh)
end

@variables J_LDH(t) Pyr(t)
@parameters VmaxLDH = 1.2mM * Hz
@parameters KpyrLDH = 47.5μM
@parameters KnadhLDH = 1.0

"""
Steady-state cytosolic calcium based on ATP:ADP ratio
"""
function ca_c_atp(atp_c, adp_c, resting, amplitude, n, katp)
    return resting + amplitude * hill(atp_c, adp_c * katp, n)
end

@variables Ca_c(t)
@parameters RestingCa = 0.09μM
@parameters ActivatedCa = 0.25μM
@parameters NCac = 4
@parameters KatpCac = 25.0

"""
Pyruvaye dehydrogenase (PDH) as well as electron transport chain (ETC)

Pyr + 4.6NAD => CO2 + 4.6NADH
"""
function j_pdh(pyr, nad_m, nadh_m, ca_m, VMAX, K_PYR, K_NAD, U1, U2, K_CA)
    c = (hillr(ca_m, K_CA))^2
    fpCa = hillr(U2 * (1 + U1 * c))
    fNAD = hill(nad_m, nadh_m * K_NAD)
    fPyr = hill(pyr, K_PYR)
    return VMAX * fPyr * fNAD * fpCa
end

@variables J_PDH(t) J_DH(t) J_CAC(t) J_FFA(t) NAD_m(t) NADH_m(t) Ca_m(t) rPDH(t)
@parameters VmaxPDH = 300μM * Hz
@parameters KpyrPDH = 47.5μM
@parameters KnadPDH = 81.0
@parameters U1PDH = 1.5
@parameters U2PDH = 1.1
@parameters KcaPDH = 0.05μM

"""
Rate of Electron trasnport chain (ETC)

NADH + O2 + 10H(mito) => NAD + 10H(cyto)
"""
function j_hr(nadh_m, dpsi, vmax, knadh, ka, kb)
    return vmax * hill(nadh_m, knadh) * (1 + ka * dpsi) / (1 + kb * dpsi)
end

@variables J_HR(t) J_O2(t) ΔΨm(t) rETC(t)
@parameters VmaxETC = 22mM * Hz
@parameters KnadhETC = 3mM
@parameters KaETC = -4.92E-3 / mV
@parameters KbETC = -4.43E-3 / mV

"""
Proton leak
"""
function j_hl(dpsi, phl, kv)
    return phl * exp(kv * dpsi)
end

@variables J_HL(t) rHL(t)
@parameters pHleak = 2.4μM * Hz
@parameters kvHleak = 0.0305 / mV

"""
F1Fo ATPase (ATP synthase) plumped with ANT

3ADP(cyto) + 9H(cyto) => 3ATP(cyto) + 9H(mito)
"""
function j_hf(adp_c, ca_m, dpsi, vmax, kadp, kv, kca, fmgadp, nadp=2, nv=8)
    fADP = hill(fmgadp * adp_c, kadp, nadp)
    fV = hill(dpsi, kv, nv)
    fCa = 1 - exp(-ca_m/kca)
    return vmax * fADP * fV * fCa
end

@variables J_HF(t) J_ANT(t) rF1(t)
@parameters VmaxF1 = 8mM * Hz
@parameters KadpF1 = 20μM
@parameters KvF1 = 131.4mV
@parameters KcaF1 = 0.165μM
@parameters FmgadpF1 = 0.055

"""
Mitochondrial calcium uniporter (MCU)

Ca(cyto) <=> Ca(mito)
"""
function j_uni(ca_m, ca_c, ΔΨ, P_CA)
    zvfrt = 2 * iVT * ΔΨ
    em1 = expm1(zvfrt)
    return P_CA * exprel(zvfrt, em1) * (0.341 * ca_c * (1 + em1) - 0.200 * ca_m)
end

@variables J_MCU(t)
@parameters PcaMCU = 4Hz

"""
Mitochondrial sodium calcium exchanger (NCLX)

Ca(mito) + 2Na(cyto) <=> Ca(cyto) + 2Na(mito)
"""
function j_nclx(ca_m, ca_c, na_m, na_c, VMAX, K_CA, K_NA)
    A = (na_c / K_NA)^2
    P = (na_m / K_NA)^2
    B = ca_m / K_CA
    Q = ca_c / K_CA
    AB = A * B
    PQ = P * Q
    return VMAX * (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
end

@variables J_NCLX(t)
@parameters Na_c = 10mM
@parameters Na_m = 5mM
@parameters VmaxNCLX = 75μM * Hz
@parameters KnaNCLX = 8.2mM
@parameters KcaNCLX = 8μM

"""
NADH shuttle

NADH(cyto) + NAD(mito) => NADH(mito) + NAD(cyto)
"""
function j_nadht(nad_c, nadh_c, nad_m, nadh_m, vmax, ktn_c, ktn_m)
    return vmax * hill(nadh_c, nad_c * ktn_c) * hill(nad_m, nadh_m * ktn_m)
end

@variables J_NADHT(t)
@parameters VmaxNADHT = 50μM * Hz
@parameters Ktn_c = 0.002
@parameters Ktn_m = 16.78

"""
AMPK activity
"""
function a_ampk(atp_c, amp_c, kampk)
    return hill(amp_c / atp_c, kampk)
end

@variables AMPKactivity(t)
@parameters kAMPK = 0.07

# Fission-fusion rates
@variables (x(t))[1:3] x13r(t) degavg(t) (v(t))[1:2]
@parameters Kfiss1 = inv(10minute)
@parameters Kfuse1 = Kfiss1
@parameters Kfiss2 = 1.5Kfiss1
@parameters Kfuse2 = 0.5Kfuse1

# Baseline consumption rates
@parameters kNADHc = 0.1Hz
@parameters kNADHm = 0.1Hz
@parameters kATP = 0.04Hz
@parameters kATPCa = 90Hz / mM
@parameters kG3P = 0.01Hz
@parameters kPyr = 0.01Hz

# Conservation relationships
@parameters ΣAc = 4.5mM
@parameters Σn_c = 2mM
@parameters Σn_m = 2.2mM
