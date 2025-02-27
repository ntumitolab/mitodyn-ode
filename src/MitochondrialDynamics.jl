module MitochondrialDynamics

using LabelledArrays
using SimpleUnPack
using PythonCall  # For export TIF
import NaNMath as nm

export make_model, build_prob, indexof, extract, exportTIF

include("utils.jl")

"Average cytosolic calcium level based on the ATP:ADP ratio"
function cac_atp(atp, adp, ca_base=0.09μM, ca_act=0.25μM, katp=25)
    return ca_base + ca_act * hil(atp, adp * katp, 4)
end

function cac_atp(u::LArray, p, t)
    @unpack RestingCa, ActivatedCa, KatpCac, kEqAK = p
    @unpack aec, ΣAc = u
    atp = aec2atp(aec, kEqAK)
    adp = aec2adp(aec, kEqAK)
    return cac_atp(atp, adp, RestingCa, ActivatedCa, KatpCac)
end


"Glucokinase (GK) rate"
jGK(atp, glc, vmax, katp, kglc, nglc) = vmax * hil(atp, katp) * hil(glc, kglc, nglc)

"Glyceraldehydes 3-phosphate dehydrogenase (GPD)"
jGPD(g3p, adp, nad, nadh, vmax, kg3p, kadp, knad) = vmax * hil(adp, kadp) * hil(nad, nadh * knad) * hil(g3p, kg3p)

"Lactate dehydrogenase (LDH)"
jLDH(pyr, nad, nadh, vmax, kpyr, knadh) = vmax * hil(pyr, kpyr) * hil(nadh, nad * knadh)

"Pyruvate dehydrogenase (PDH)"
function jPDH(pyr, nad_m, nadh_m, ca_m, vmax, kpyr, knad, kca, u1, u2, rPDH=1)
    c = (hil(kca, ca_m))^2
    fpCa = hil(1, u2 * (1 + u1 * c))
    fNAD = hil(nad_m, nadh_m * knad)
    fPyr = hil(pyr, kpyr)
    return vmax * rPDH * fpCa * fNAD * fPyr
end

"Electron transport chain (ETC) and O2 consumption"
function jETC(nadh_m, dpsi, vmax, knadh, ka, kb, rETC=1, rO2=0.1)
    jhr = vmax * rETC * hil(nadh_m, knadh) * (1 + ka * dpsi) / (1 + kb * dpsi)
    return (jhr, jhr * rO2)
end

"Proton leak"
jHL(dpsi, p, kv, rHL=1) = p * rHL * exp(kv * dpsi)

"F1Fo ATPase (ATP synthase) lumped with ANT"
function jHF(adp, dpsi, ca_m, vmax, fmgadp, kadp, nadp, kv, nv, kca, rF1=1)
    jhf = vmax * rF1 * hil(fmgadp * adp, kadp, nadp) * hil(dpsi, kv, nv) * (1 - exp(-ca_m / kca))
    jant = jhf / 3
    return (jhf, jant)
end

"Mitochondrial calcium uniporter (MCU) and Mitochondrial sodium-calcium exchanger (NCLX)"
function jMitoCa(dpsi, ca_c, ca_m, na_c, na_m, PcaMCU, VmaxNCLX, KnaNCLX, KcaNCLX)
    zvfrt = 2 * iVT * dpsi
    em1 = expm1(zvfrt)
    jmcu = PcaMCU * (zvfrt / em1) * (0.341 * ca_c * (em1 + 1) - 0.2 * ca_m)
    A = (na_c / KnaNCLX)^2
    P = (na_m / KnaNCLX)^2
    B = ca_m / KcaNCLX
    Q = ca_c / KcaNCLX
    AB = A * B
    PQ = P * Q
    jnclx = VmaxNCLX * (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
    return (jmcu, jnclx)
end

"NADH shuttle"
jNADHT(nad_c, nadh_c, nad_m, nadh_m, vmax, Ktn_c, Ktn_m) = vmax * hil(nadh_c, nad_c * Ktn_c) * hil(nad_m, nadh_m * Ktn_m)

"Constant glucose equation"
function const_glc(glc=5mM)
    @variables Glc(t)
    @parameters GlcConst = glc
    return Glc ~ GlcConst
end

"Model ODE"
function mitodyn_rhs!(D, u, p, t; cacfunc = cac_atp)
    @unpack RestingCa, ActivatedCa, KatpCac, kEqAK, ΣAc, Σn_c, Σn_m = p
    @unpack aec, nad_c, nad_m, ca_m = u
    atp = ΣAc * aec2atp(aec, kEqAK)
    adp = ΣAc * aec2adp(aec, kEqAK)
    amp = ΣAc - atp - adp
    nadh_c = Σn_c - nad_c
    nadh_m = Σn_m - nad_m
    ca_c = cacfunc(u, p, t)
    @unpack Glc, VmaxGK, KatpGK, KglcGK, nGK = p
    j_gk = jGK(atp, Glc, VmaxGK, KatpGK, KglcGK, nGK)
    j_gpd =
    return (; j_gk, j_gpd, atp, adp, amp, nadh_c, nadh_m)
end


"Build parameters"
function build_params()
    return LVector(
        C_MIT=1.812μM / mV, # Mitochondrial membrane capacitance
        F_M=3E-4,       # Fraction of free Ca in mitochondria
        F_I=0.01,       # Fraction of free Ca in cytosol
        V_I=0.53,       # Relative cytoplasmic volume
        V_MT=0.06,      # Relative mitochondrial volume
        V_MTX=0.0144,   # Relative mitochondrial matrix volume
        rETC=1.0,
        rHL=1.0,
        rF1=1.0,
        rPDH=1.0,
        Glc=5mM,
        VmaxGK=0.011mM * Hz,
        KatpGK=0.5mM,
        KglcGK=7mM,
        nGK=1.7,
        VmaxGPD=0.5mM * Hz,
        Kg3pGPD=0.2mM,
        KnadGPD=0.09,
        KadpGPD=2μM,
        VmaxLDH=1.2mM * Hz,
        KpyrLDH=47.5μM,
        KnadhLDH=1.0,
        VmaxPDH=300μM * Hz,
        KpyrPDH=47.5μM,
        KnadPDH=81,
        U1PDH=1.5,
        U2PDH=1.1,
        KcaPDH=0.05μM,
        kFFA=0Hz,
        VmaxETC=22mM * Hz,
        KnadhETC=3mM,
        KaETC=-4.92E-3 / mV,
        KbETC=-4.43E-3 / mV,
        pHleak=2.4μM * Hz,
        kvHleak=0.0305 / mV,
        VmaxF1=8mM * Hz,
        KadpF1=20μM,
        KvF1=131.4mV,
        KcaF1=0.165μM,
        FmgadpF1=0.055,
        nadpF1=2,
        nvF1=8,
        PcaMCU = 4Hz,
        Na_c=10mM,
        Na_m=5mM,
        VmaxNCLX=75μM * Hz,
        KnaNCLX=8.2mM,
        KcaNCLX=8μM,
        VmaxNADHT=50μM * Hz,
        Ktn_c=0.002,
        Ktn_m=16.78,
        RestingCa=0.09μM,
        ActivatedCa=0.25μM,
        KatpCac=25,
        kEqAK = 0.931,
        ΣAc=4.56mM,
        Σn_c=2mM,
        Σn_m=2.2mM,
        kNADHc=0.1Hz,
        kNADHm=0.1Hz,
        kATP=0.04Hz,
        kATPCa=90Hz / mM,
        kG3P=0.01Hz,
        kPyr=0.01Hz
    )
end

function make_model(;
    name,
    caceq=cac_atp(),
    glceq=const_glc(5mM),
    simplify=true,
    kwargs...
)
    @constants begin
        C_MIT = 1.812μM / mV # Mitochondrial membrane capacitance
        F_M = 3E-4       # Fraction of free Ca in mitochondria
        F_I = 0.01       # Fraction of free Ca in cytosol
        V_I = 0.53       # Relative cytoplasmic volume
        V_MT = 0.06      # Relative mitochondrial volume
        V_MTX = 0.0144   # Relative mitochondrial matrix volume
        ATPstiochGK = 2
    end

    # Adjustable parameters
    @parameters (rETC=1, rHL=1, rF1=1, rPDH=1)

    # Adenylate kinase (AdK)
    @variables AEC(t) ATP_c(t) ADP_c(t) AMP_c(t)
    @parameters (kfAK=1000Hz / mM, kEqAK=0.931)

    # Glucokinase (GK)
    @variables Glc(t) J_GK(t)
    @parameters (VmaxGK=0.011mM * Hz, KatpGK=0.5mM, KglcGK=7mM, nGK=1.7, ATPstiochGK=2)

    # Glyceraldehydes 3-phosphate dehydrogenase (GPD)
    @variables J_GPD(t) G3P(t) NAD_c(t) NADH_c(t)
    @parameters (VmaxGPD=0.5mM * Hz, Kg3pGPD=0.2mM, KnadGPD=0.09, KadpGPD=2μM)

    # Lactate dehydrogenase (LDH)
    @variables J_LDH(t) Pyr(t)
    @parameters (VmaxLDH=1.2mM * Hz, KpyrLDH=47.5μM, KnadhLDH=1)

    # Pyruvate dehydrogenase (PDH)
    @variables J_PDH(t) J_DH(t) NAD_m(t) NADH_m(t) Ca_m(t) J_FFA(t)
    @parameters (VmaxPDH=300μM * Hz, KpyrPDH=47.5μM, KnadPDH=81, U1PDH=1.5, U2PDH=1.1, KcaPDH=0.05μM, kFFA=0Hz)
    pdheq = let
        c = (hil(KcaPDH, Ca_m))^2
        fpCa = hil(1, U2PDH * (1 + U1PDH * c))
        fNAD = hil(NAD_m, NADH_m * KnadPDH)
        fPyr = hil(Pyr, KpyrPDH)
        J_PDH ~ VmaxPDH * rPDH * fpCa * fNAD * fPyr
    end

    # Electron transport chain (ETC)
    @variables J_HR(t) J_O2(t) ΔΨm(t)
    @parameters (VmaxETC=22mM * Hz, KnadhETC=3mM, KaETC=-4.92E-3 / mV, KbETC=-4.43E-3 / mV)

    # Proton leak
    @variables J_HL(t)
    @parameters (pHleak=2.4μM * Hz, kvHleak=0.0305 / mV)

    # F1Fo ATPase (ATP synthase) lumped with ANT
    @variables J_HF(t) J_ANT(t)
    @parameters (VmaxF1=8mM * Hz, KadpF1=20μM, KvF1=131.4mV, KcaF1=0.165μM, FmgadpF1=0.055, nadpF1=2, nvF1=8)

    # Mitochondrial calcium uniporter (MCU)
    @variables J_MCU(t) Ca_c(t)
    @parameters PcaMCU = 4Hz
    mcueq = let
        zvfrt = 2 * iVT * ΔΨm
        em1 = expm1(zvfrt)
        J_MCU ~ PcaMCU * (zvfrt / em1) * (0.341 * Ca_c * (em1 + 1) - 0.2 * Ca_m)
    end

    # Mitochondrial sodium-calcium exchanger (NCLX)
    @variables J_NCLX(t)
    @parameters (Na_c=10mM, Na_m=5mM, VmaxNCLX=75μM * Hz, KnaNCLX=8.2mM, KcaNCLX=8μM)
    nclxeq = let
        A = (Na_c / KnaNCLX)^2
        P = (Na_m / KnaNCLX)^2
        B = Ca_m / KcaNCLX
        Q = Ca_c / KcaNCLX
        AB = A * B
        PQ = P * Q
        J_NCLX ~ VmaxNCLX * (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
    end

    # NADH shuttle
    @variables J_NADHT(t)
    @parameters (VmaxNADHT=50μM * Hz, Ktn_c=0.002, Ktn_m=16.78)

    # Baseline consumption rates
    @parameters (kNADHc=0.1Hz, kNADHm=0.1Hz, kATP=0.04Hz, kATPCa=90Hz / mM, kG3P=0.01Hz, kPyr=0.01Hz)
    # Conservation relationships
    @parameters (ΣAc=4.56mM, Σn_c=2mM, Σn_m=2.2mM)

    # Fission-fusion rates
    @variables x1(t) x2(t) x3(t) degavg(t) tiptip(t) tipside(t)
    @parameters (kfiss1=inv(10minute), kfuse1=kfiss1, kfiss2=1.5 * kfiss1, kfuse2=0.5 * kfuse1)

    eqs = [
        caceq,
        glceq,
        J_GK ~ VmaxGK * hil(ATP_c, KatpGK) * hil(Glc, KglcGK, nGK),
        J_GPD ~ VmaxGPD * hil(ADP_c, KadpGPD) * hil(NAD_c, NADH_c * KnadGPD) * hil(G3P, Kg3pGPD),
        J_LDH ~ VmaxLDH * hil(Pyr, KpyrLDH) * hil(NADH_c, NAD_c * KnadhLDH),
        pdheq,
        J_FFA ~ kFFA * NAD_m,
        J_DH ~ 4.6 * J_PDH + J_FFA,
        J_HR ~ VmaxETC * rETC * hil(NADH_m, KnadhETC) * (1 + KaETC * ΔΨm) / (1 + KbETC * ΔΨm),
        J_O2 ~ J_HR / 10,
        J_HL ~ pHleak * rHL * exp(kvHleak * ΔΨm),
        J_HF ~ VmaxF1 * rF1 * hil(FmgadpF1 * ADP_c, KadpF1, nadpF1) * hil(ΔΨm, KvF1, nvF1) * (1 - exp(-Ca_m / KcaF1)),
        J_ANT ~ J_HF / 3,
        mcueq,
        nclxeq,
        J_NADHT ~ VmaxNADHT * hil(NADH_c, NAD_c * Ktn_c) * hil(NAD_m, NADH_m * Ktn_m),
        tiptip ~ kfuse1 * J_ANT / J_HL * x1 * x1 - kfiss1 * x2,
        tipside ~ kfuse2 * J_ANT / J_HL * x1 * x2 - kfiss2 * x3,
        degavg ~ (x1 + 2x2 + 3x3) / (x1 + x2 + x3),
        # Conservation relationships
        ATP_c ~ ΣAc * aec2atp(AEC, kEqAK),
        ADP_c ~ ΣAc * aec2adp(AEC, kEqAK),
        AMP_c ~ ΣAc - ATP_c - ADP_c,
        Σn_c ~ NADH_c + NAD_c,
        Σn_m ~ NADH_m + NAD_m,
        1 ~ x1 + 2x2 + 3x3,
        # State variables
        D(NADH_m) ~ inv(V_MTX) * (J_DH + J_NADHT - J_O2) - kNADHm * NADH_m,
        # D(NAD_m) ~ # Conserved
        D(NADH_c) ~ inv(V_I) * (J_GPD - J_NADHT - J_LDH) - kNADHc * NADH_c,
        # D(NAD_c) ~ # Conserved
        D(Ca_m) ~ inv(V_MTX) * F_M * (J_MCU - J_NCLX),
        D(ΔΨm) ~ inv(C_MIT) * (J_HR - J_HF - J_HL - J_ANT - 2 * J_MCU),
        D(G3P) ~ inv(V_I) * (2J_GK - J_GPD) - kG3P * G3P,
        D(Pyr) ~ inv(V_I + V_MTX) * (J_GPD - J_PDH - J_LDH) - kPyr * Pyr,
        D(AEC) ~ (inv(V_I) * (-ATPstiochGK * J_GK + 2 * J_GPD + J_ANT) - ATP_c * (kATP + kATPCa * Ca_c)) / ΣAc,
        # D(x[1]) ~ # Conserved
        D(x2) ~ tiptip - tipside,
        D(x3) ~ tipside,
    ]

    sys = ODESystem(eqs, t; name,
        defaults=[
            G3P => 2.9μM,
            Pyr => 8.7μM,
            NADH_c => 1μM,
            NADH_m => 57μM,
            AEC => 0.9,
            Ca_m => 0.200μM,
            ΔΨm => 92mV,
            x2 => 0.24,
            x3 => 0.06
        ],
        kwargs...
    )

    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end

end # Module
