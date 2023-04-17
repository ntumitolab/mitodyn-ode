module MitochondrialDynamics

using ModelingToolkit

export make_model

include("utils.jl")
include("rates.jl")



hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

"Steady-state cytosolic calcium based on ATP:ADP ratio"
function cac_eq_atp()
    @variables t Ca_c(t) ATP_c(t) ADP_c(t)
    @parameters RestingCa=0.09μM ActivatedCa=0.25μM NCac=4 KatpCac=25
    caceq = Ca_c(t) ~ RestingCa + ActivatedCa * hil(ATP_c, KatpCac * ADP_c, NCac)
    return caceq
end

function glc_eq()
    @variables t Glc(t)
    @parameters GlcConst=5mM
    glceq = Glc(t) ~ GlcConst
    return glceq
end

function make_model(;
    name,
    simplify=true,
    caceq=cac_eq_atp(),
    glceq=glc_eq(),
    rpdh=1,
    retc=1,
    rf1=1,
    rhleak=1,
    j_ffa=0,
    gk_atp_stoich::Int=2
)
    @variables t

    # Adenylate kinase (AdK)
    @variables J_ADK(t) ATP_c(t) ADP_c(t) AMP_c(t)
    @parameters kfAK=1000Hz/mM kEqAK=0.931
    adkeq = J_ADK(t) ~ kfAK * (ADP_c * ADP_c - ATP_c * AMP_c * kEqAK)

    # Glucokinase (GK)
    @variables Glc(t) J_GK(t)
    @parameters VmaxGK=0.011mM*Hz KatpGK=0.5mM KglcGK=7mM nGK=1.7
    gkeq = J_GK(t) ~ VmaxGK * hil(ATP_c, KatpGK) * hil(Glc, KglcGK, nGK)

    # Glyceraldehydes 3-phosphate dehydrogenase (GPD)
    @variables J_GPD(t) G3P(t) NAD_c(t) NADH_c(t)
    @parameters VmaxGPD=0.5mM*Hz Kg3pGPD=0.2mM KnadGPD=0.09 KadpGPD=2μM
    gpdeq = J_GPD(t) ~ VmaxGPD * hil(ADP_c, KadpGPD) * hil(NAD_c, NADH_c * KnadGPD) * hil(G3P, Kg3pGPD)

    # Lactate dehydrogenase (LDH)
    @variables J_LDH(t) Pyr(t)
    @parameters VmaxLDH=1.2mM*Hz KpyrLDH=47.5μM KnadhLDH=1
    ldheq = J_LDH(t) ~ VmaxLDH * hil(Pyr, KpyrLDH) * hil(NADH_c, NAD_c * KnadhLDH)

    # Pyruvaye dehydrogenase (PDH)
    @variables J_PDH(t) J_DH(t) J_CAC(t) J_FFA(t) NAD_m(t) NADH_m(t) Ca_m(t)
    @parameters VmaxPDH=300μM*Hz KpyrPDH=47.5μM KnadPDH=81 U1PDH=1.5 U2PDH=1.1 KcaPDH=0.05μM
    pdheq = let
        c = (hil(KcaPDH, Ca_m))^2
        fpCa = hil(1, U2PDH * (1 + U1PDH * c))
        fNAD = hil(NAD_m, NADH_m * KnadPDH)
        fPyr = hil(Pyr, KpyrPDH)
        J_PDH ~ VmaxPDH * fpCa * fNAD * fPyr
    end

    # Electron trasnport chain (ETC)
    @variables J_HR(t) J_O2(t) ΔΨm(t)
    @parameters VmaxETC=22mM*Hz KnadhETC=3mM KaETC=-4.92E-3/mV KbETC=-4.43E-3/mV
    hreq = J_HR ~ VmaxETC * hill(NADH_m, KnadhETC) * (1 + KaETC * ΔΨm) / (1 + KbETC * ΔΨm)

    # Proton leak
    @variables J_HL(t)
    @parameters pHleak=2.4μM*Hz kvHleak=0.0305/mV
    hkeq = J_HL ~ pHleak * exp(kvHleak * ΔΨm)

    # F1Fo ATPase (ATP synthase) lumped with ANT
    @variables J_HF(t) J_ANT(t)
    @parameters VmaxF1=8mM*Hz KadpF1=20μM KvF1=131.4mV KcaF1=0.165μM FmgadpF1=0.055 nadpF1=2 nvF1=8
    hfeq = let
        fADP = hil(FmgadpF1 * ADP_c, KadpF1, nadpF1)
        fV = hil(ΔΨm, KvF1, nvF1)
        fCa = 1 - exp(-Ca_m / KcaF1)
        J_HF ~ VmaxF1 * fADP * fV * fCa
    end

    # Mitochondrial calcium uniporter (MCU)
    @variables J_MCU(t)
    @parameters PcaMCU=4Hz
    mcueq = let
        zvfrt = 2 * iVT * ΔΨm
        em1 = expm1(zvfrt)
        J_MCU ~ PcaMCU * zvfrt / em1 * (0.341 * Ca_c * (1 + em1) - 0.200 * Ca_m)
    end

    # Mitochondrial sodium calcium exchanger (NCLX)
    @variables J_NCLX(t)
    @parameters Na_c=10mM Na_m=5mM VmaxNCLX=75μM*Hz KnaNCLX=8.2mM KcaNCLX=8μM
    nclxeq = let
        A = (Na_c / KnaNCLX)^2
        P = (Na_m / KnaNCLX)^2
        B = Ca_m / KcaNCLX
        Q = Ca_c / KcaNCLX
        AB = A * B
        PQ = P * Q
        J_NCLX ~ (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
    end

    # NADH shuttle
    @variables J_NADHT(t)
    @parameters VmaxNADHT=50μM*Hz =0.002 Ktn_m=16.78
    nadhteq = J_NADHT ~ VmaxNADHT * hil(NADH_c, NAD_c * Ktn_c) * hil(NAD_m, NADH_m, Ktn_m)

    # Baseline consumption rates
    @parameters kNADHc=0.1Hz kNADHm=0.1Hz kATP=0.04Hz kATPCa=90Hz/mM kG3P=0.01Hz kPyr=0.01Hz
    # Conservation relationships
    @parameters ΣAc=4.5mM Σn_c=2mM Σn_m=2.2mM

    # Fission-fusion rates
    @variables (x(t))[1:3] degavg(t)
    @parameters kfiss1=inv(10minute) kfuse1=kfiss1 kfiss2=1.5kfiss1 kfuse2=0.5kfuse1

    D = Differential(t)
    eqs = [
        # Reactions
        J_GK ~ j_gk(ATP_c, Glc, VmaxGK, KatpGK, KglcGK, nGK),
        J_GPD ~ j_gpd(G3P, NAD_c, NADH_c, ADP_c, VmaxGPD, Kg3pGPD, KnadGPD, KadpGPD),
        J_LDH ~ j_ldh(Pyr, NAD_c, NADH_c, VmaxLDH, KpyrLDH, KnadhLDH),
        J_ADK ~ j_adk(ATP_c, ADP_c, AMP_c, kfAK, kEqAK),
        J_PDH ~ rPDH * j_pdh(Pyr, NAD_m, NADH_m, Ca_m, VmaxPDH, KpyrPDH, KnadPDH, U1PDH, U2PDH, KcaPDH),
        J_CAC ~ J_PDH + j_ffa,
        J_DH ~ 4.6 * J_CAC,
        J_HR ~ rETC * j_hr(NADH_m, ΔΨm, VmaxETC, KnadhETC, KaETC, KbETC),
        J_O2 ~ 0.1 * J_HR,
        J_HL ~  rHL * j_hl(ΔΨm, pHleak, kvHleak),
        J_HF ~ rF1 * j_hf(ADP_c, Ca_m, ΔΨm, VmaxF1, KadpF1, KvF1, KcaF1, FmgadpF1),
        J_ANT ~ J_HF / 3,
        J_MCU ~ j_uni(Ca_m, Ca_c, ΔΨm, PcaMCU),
        J_NCLX ~ j_nclx(Ca_m, Ca_c, Na_m, Na_c, VmaxNCLX, KcaNCLX, KnaNCLX),
        J_NADHT ~ j_nadht(NAD_c, NADH_c, NAD_m, NADH_m, VmaxNADHT, Ktn_c, Ktn_m),
        v[1] ~ kfuse1 * J_ANT / J_HL * x[1] * x[1] - kfiss1 * x[2],
        v[2] ~ kfuse2 * J_ANT / J_HL * x[1] * x[2] - kfiss2 * x[3],
        Glc ~ glcrhs,
        Ca_c ~ cacrhs,
        rPDH ~ rpdh,
        rETC ~ retc,
        rF1 ~ rf1,
        rHL ~ rhleak,
        # Conservation relationships
        ΣAc ~ ATP_c + ADP_c + AMP_c,
        Σn_c ~ NADH_c + NAD_c,
        Σn_m ~ NADH_m + NAD_m,
        1 ~ x[1] + 2x[2] + 3x[3],
        # Observables
        x13r ~ x[1] / x[3],
        degavg ~ (x[1] + 2x[2] + 3x[3]) / (x[1] + x[2] + x[3]),
        # AMPKactivity ~ hill(AMP_c / ATP_c, kAMPK),
        # State variables
        D(NADH_m) ~ iVmtx * (J_DH + J_NADHT - J_O2) - kNADHm * NADH_m,
        # D(NAD_m) ~ # Conserved
        D(NADH_c) ~ iVi * (J_GPD - J_NADHT - J_LDH) - kNADHc * NADH_c,
        # D(NAD_c) ~ # Conserved
        D(Ca_m) ~ iVmtx * F_M * (J_MCU - J_NCLX),
        D(ΔΨm) ~ iCmt * (J_HR - J_HF - J_HL - J_ANT - 2J_MCU),
        D(G3P) ~ iVi * (2J_GK - J_GPD) - kG3P * G3P,
        D(Pyr) ~ iVimtx * (J_GPD - J_PDH - J_LDH) - kPyr * Pyr,
        D(ATP_c) ~ iVi * (-gk_atp_stoich * J_GK + 2J_GPD + J_ANT + J_ADK) - ATP_c * (kATP + kATPCa * Ca_c),
        D(AMP_c) ~ iVi * J_ADK,
        # D(ADP_c) ~ # Conserved
        # D(x[1]) ~ # Conserved
        D(x[2]) ~ v[1] - v[2],
        D(x[3]) ~ v[2],
    ]
    sys = ODESystem(eqs, t; name,
        defaults=[
            G3P => 2.9μM,
            Pyr => 8.7μM,
            NADH_c => 1μM,
            NADH_m => 57μM,
            ATP_c => 3595μM,
            AMP_c => 148μM,
            Ca_m => 0.200μM,
            ΔΨm => 92mV,
            x[2] => 0.24,
            x[3] => 0.06
        ]
    )

    if simplify
        sys = structural_simplify(sys)
    end

    return sys
end


end
