import .Utils: hill, hillr, exprel, mM, μM, ms, minute, Hz, mV, iVT, iVmtx, iCmt, iVi, iVimtx, F_M

@variables t

#=
Adenylate kinase:

2ADP <=> ATP + AMP
=#
@variables J_ADK(t) ATP_c(t) ADP_c(t) AMP_c(t)
@parameters kfAK = 1000Hz / mM kEqAK = 0.931

#=
Glucokinase (GK) Reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP (+ ADP) =>> 2G3P (+ ADP)
=#
@variables Glc(t) J_GK(t)
@parameters VmaxGK = 0.011mM * Hz KatpGK = 0.5mM KglcGK = 7mM NGK = 1.7 GlcConst = 5mM

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
function j_uni(ca_m, ca_c, ΔΨ, P_CA)
    zvfrt = 2 * iVT * ΔΨ
    em1 = expm1(zvfrt)
    return P_CA * exprel(zvfrt, em1) * (0.341 * ca_c * (1 + em1) - 0.200 * ca_m)
end

#=
Mitochondrial sodium calcium exchanger (NCLX)

Ca(mito) + 2Na(cyto) <=> Ca(cyto) + 2Na(mito)
=#
@variables J_NCLX(t)
@parameters Na_c = 10mM Na_m = 5mM VmaxNCLX = 75μM * Hz KnaNCLX = 8.2mM KcaNCLX = 8μM
function j_nclx(ca_m, ca_c, na_m, na_c, VMAX, K_CA, K_NA)
    A = (na_c / K_NA)^2
    P = (na_m / K_NA)^2
    B = ca_m / K_CA
    Q = ca_c / K_CA
    AB = A * B
    PQ = P * Q
    return VMAX * (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
end

#=
NADH shuttle

NADH(cyto) + NAD(mito) => NADH(mito) + NAD(cyto)
=#
@variables J_NADHT(t)
@parameters VmaxNADHT = 50μM * Hz Ktn_c = 0.002 Ktn_m = 16.78

# Fission-fusion rates
@variables x[1:3](t) x13r(t) degavg(t) v[1:2](t)
@parameters Kfiss1 = inv(10minute) Kfuse1 = Kfiss1 Kfiss2 = 1.5Kfiss1 Kfuse2 = 0.5Kfuse1

# Baseline consumption rates
@parameters kNADHc = 0.1Hz kNADHm = 0.1Hz kATP = 0.04Hz kATPCa = 90Hz / mM kG3P = 0.01Hz kPyr = 0.01Hz

# Conservation relationships
@parameters ΣAc = 4.5mM Σn_c = 2mM Σn_m = 2.2mM

function make_model(;
    name,
    simplify=true,
    cacrhs=RestingCa + ActivatedCa * hill(ATP_c / KatpCac, ADP_c, NCac),
    glcrhs=GlcConst,
    rpdh=1,
    retc=1,
    rf1=1,
    rhleak=1,
    ffarhs=0
)
    D = Differential(t)
    eqs = [
        # Reactions
        J_GK ~ VmaxGK * hill(ATP_c, KatpGK) * hill(Glc, KglcGK, NGK),
        J_GPD ~ VmaxGPD * hill(ADP_c, KadpGPD) * hill(NAD_c / KnadGPD, NADH_c) * hill(G3P, Kg3pGPD),
        J_LDH ~ VmaxLDH * hill(Pyr, KpyrLDH) * hill(NADH_c / KnadhLDH, NAD_c),
        J_ADK ~ kfAK * (ADP_c * ADP_c - ATP_c * AMP_c / kEqAK),
        J_PDH ~ rPDH * j_pdh(Pyr, NAD_m, NADH_m, Ca_m, VmaxPDH, KpyrPDH, KnadPDH, U1PDH, U2PDH, KcaPDH),
        J_CAC ~ J_PDH + J_FFA,
        J_DH ~ 4.6 * J_CAC,
        J_HR ~ VmaxETC * rETC * hill(NADH_m, KnadhETC) * (1 + KaETC * ΔΨm) / (1 + KbETC * ΔΨm),
        J_O2 ~ 0.1 * J_HR,
        J_HL ~ pHleak * rHL * exp(kvHleak * ΔΨm),
        J_HF ~ VmaxF1 * rF1 * hill(FmgadpF1 * ADP_c, KadpF1, 2) * hill(ΔΨm, KvF1, 8) * (1 - exp(-Ca_m / KcaF1)),
        J_ANT ~ J_HF / 3,
        J_MCU ~ j_uni(Ca_m, Ca_c, ΔΨm, PcaMCU),
        J_NCLX ~ j_nclx(Ca_m, Ca_c, Na_m, Na_c, VmaxNCLX, KcaNCLX, KnaNCLX),
        J_NADHT ~ VmaxNADHT * hill(NADH_c / Ktn_c, NAD_c) * hill(NAD_m / Ktn_m, NADH_m),
        v[1] ~ Kfuse1 * J_ANT / J_HL * x[1] * x[1] - Kfiss1 * x[2],
        v[2] ~ Kfuse2 * J_ANT / J_HL * x[1] * x[2] - Kfiss2 * x[3],
        Glc ~ glcrhs,
        Ca_c ~ cacrhs,
        rPDH ~ rpdh,
        rETC ~ retc,
        rF1 ~ rf1,
        rHL ~ rhleak,
        J_FFA ~ ffarhs,
        # Conservation relationships
        ΣAc ~ ATP_c + ADP_c + AMP_c,
        Σn_c ~ NADH_c + NAD_c,
        Σn_m ~ NADH_m + NAD_m,
        1 ~ x[1] + 2x[2] + 3x[3],
        # Observables
        x13r ~ x[1] / x[3],
        degavg ~ (x[1] + 2x[2] + 3x[3]) / (x[1] + x[2] + x[3]),
        AMPKactivity ~ hill(AMP_c / ATP_c, kAMPK),
        # State variables
        D(NADH_m) ~ iVmtx * (J_DH + J_NADHT - J_O2) - kNADHm * NADH_m,
        # D(NAD_m) ~ # Conserved
        D(NADH_c) ~ iVi * (J_GPD - J_NADHT - J_LDH) - kNADHc * NADH_c,
        # D(NAD_c) ~ # Conserved
        D(Ca_m) ~ iVmtx * F_M * (J_MCU - J_NCLX),
        D(ΔΨm) ~ iCmt * (J_HR - J_HF - J_HL - J_ANT - 2J_MCU),
        D(G3P) ~ iVi * (2J_GK - J_GPD) - kG3P * G3P,
        D(Pyr) ~ iVimtx * (J_GPD - J_PDH - J_LDH) - kPyr * Pyr,
        D(ATP_c) ~ iVi * (-2J_GK + 2J_GPD + J_ANT + J_ADK) - ATP_c * (kATP + kATPCa * Ca_c),
        D(AMP_c) ~ iVi * J_ADK,
        # D(ADP_c) ~ # Conserved
        # D(x[1]) ~ # Conserved
        D(x[2]) ~ v[1] - v[2],
        D(x[3]) ~ v[2],
    ]
    sys = ODESystem(eqs, t; name,
        defaults=Dict(
            G3P => 2.8μM,
            Pyr => 8.5μM,
            NADH_c => 1μM,
            NADH_m => 60μM,
            ATP_c => 4000μM,
            AMP_c => 0μM,
            Ca_m => 0.250μM,
            ΔΨm => 100mV,
            x[2] => 0.25,
            x[3] => 0.05
        ))

    if simplify
        sys = structural_simplify(sys)
    end

    return sys
end
