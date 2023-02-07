module MitochondrialDynamics

using ModelingToolkit

export make_model

include("utils.jl")
include("rates.jl")

function make_model(;
    name,
    simplify=true,
    cacrhs=RestingCa + ActivatedCa * hill(ATP_c, ADP_c * KatpCac, NCac),
    glcrhs=GlcConst,
    rpdh=1,
    retc=1,
    rf1=1,
    rhleak=1,
    j_ffa=0,
    gk_atp_stoich::Int=2
)
    D = Differential(t)
    eqs = [
        # Reactions
        J_GK ~ j_gk(ATP_c, Glc, VmaxGK, KatpGK, KglcGK, NGK),
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
        v[1] ~ Kfuse1 * J_ANT / J_HL * x[1] * x[1] - Kfiss1 * x[2],
        v[2] ~ Kfuse2 * J_ANT / J_HL * x[1] * x[2] - Kfiss2 * x[3],
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

end # module MitochondrialDynamics
