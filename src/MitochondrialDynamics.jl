module MitochondrialDynamics

using ModelingToolkit

export make_model

include("utils.jl")
include("rates.jl")


# Convineince function
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
    gk_atp_stoich::Int=2
)
    @constants C_MIT=1.812μM/mV # Mitochondrial membrane capacitance
    @constants F_M = 3E-4       # Frasction of free Ca in mitochondria
    @constants F_I = 0.01       # Fraction of free Ca in cytosol
    @constants V_I = 0.53       # Relative cytoplasmic volume
    @constants V_MT = 0.06      # Relative mitochondrial volume
    @constants V_MTX = 0.0144   # Relative mitochondrial matrix volume

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
    @variables J_PDH(t) J_DH(t) J_CAC(t) NAD_m(t) NADH_m(t) Ca_m(t)
    @parameters VmaxPDH=300μM*Hz KpyrPDH=47.5μM KnadPDH=81 U1PDH=1.5 U2PDH=1.1 KcaPDH=0.05μM J_FFA=0μM*Hz
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
    @variables (x(t))[1:3] degavg(t) v1(t) v2(t)
    @parameters kfiss1=inv(10minute) kfuse1=kfiss1 kfiss2=1.5*kfiss1 kfuse2=0.5*kfuse1

    D = Differential(t)
    # Reaction equations
    eqs = [
        caceq,
        glceq,
        adkeq,
        gkeq,
        gpdeq,
        ldheq,
        pdheq,
        hreq,
        hkeq,
        hfeq,
        mcueq,
        nclxeq,
        nadhteq,
        J_CAC ~ J_PDH + J_FFA,
        J_DH ~ 4.6 * J_CAC,
        J_O2 ~ J_HR / 10,
        J_ANT ~ J_HF / 3,
        v1 ~ kfuse1 * J_ANT / J_HL * x[1] * x[1] - kfiss1 * x[2],
        v2 ~ kfuse2 * J_ANT / J_HL * x[1] * x[2] - kfiss2 * x[3],
        # Conservation relationships
        ΣAc ~ ATP_c + ADP_c + AMP_c,
        Σn_c ~ NADH_c + NAD_c,
        Σn_m ~ NADH_m + NAD_m,
        1 ~ x[1] + 2x[2] + 3x[3],
        # Observables
        degavg ~ (x[1] + 2x[2] + 3x[3]) / (x[1] + x[2] + x[3]),
        # State variables
        V_MTX * D(NADH_m) ~ J_DH + J_NADHT - J_O2 - V_MTX * kNADHm * NADH_m,
        # D(NAD_m) ~ # Conserved
        V_I * D(NADH_c) ~ J_GPD - J_NADHT - J_LDH - V_I * kNADHc * NADH_c,
        # D(NAD_c) ~ # Conserved
        V_MTX * D(Ca_m) ~ F_M * (J_MCU - J_NCLX),
        C_MIT * D(ΔΨm) ~ J_HR - J_HF - J_HL - J_ANT - 2 * J_MCU,
        V_I * D(G3P) ~ 2J_GK - J_GPD - V_I * kG3P * G3P,
        (V_I + V_MTX) * D(Pyr) ~ (J_GPD - J_PDH - J_LDH) - (V_I + V_MTX) * kPyr * Pyr,
        V_I * D(ATP_c) ~ (-gk_atp_stoich * J_GK + 2 * J_GPD + J_ANT + J_ADK) - V_I * ATP_c * (kATP + kATPCa * Ca_c),
        V_I * D(AMP_c) ~ J_ADK,
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
