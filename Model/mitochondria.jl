# Mitochondria compartment descriptions
using Parameters
using Setfield
using FromFile
@from "utils.jl" import exprel, hill, hillr, iVT, Hz, mM, μM, mV

"Pyruvaye dehydrogenase (PDH)"
@with_kw struct PDH{R}
    VMAX::R = 300μM*Hz    # Max rate for PDH
    K_PYR::R = 47.5μM     # Apparent Michaelis constant for pyruvate
    K_NAD::R = 81.0       # Apparent activation constant for NAD/NADH ratio
	U1::R = 1.5           # Factor for calcium activation
	U2::R = 1.1           # Factor for calcium activation
	K_CA::R = 0.05μM      # Apparent dissociation constant for calcium
end

"Rate of Pyruvaye dehydrogenase (PDH)"
function j_pdh(pyr, nad_m, nadh_m, ca_m, p::PDH)
    @unpack VMAX, K_PYR, K_NAD, U1, U2, K_CA = p
    c = (hillr(ca_m, K_CA))^2
    fpCa = hillr(U2 * (1 + U1 * c))
    fNAD = hill(nad_m, nadh_m * K_NAD)
    fPyr = hill(pyr, K_PYR)
    jPDH = VMAX * fPyr * fNAD * fpCa
    return jPDH
end

(p::PDH)(pyr, nad_m, nadh_m, ca_m) = j_pdh(pyr, nad_m, nadh_m, ca_m, p)

"NADH production rate"
j_dh(jPDH, NADH_PYR_RATIO=46//10) = NADH_PYR_RATIO * jPDH

"""
Electron trasnport chain (ETC)

Reaction: NADH + O2 + 10H(mito) => NAD + 10H(cyto)
"""
@with_kw struct ETC{R}
    VMAX::R   = 22mM*Hz       # Maximal respiratory rate
    K_NADH::R = 3mM           # Apparent Michaelis constant for NADH
    K_A::R    = -4.92E-3/mV   # Voltage factor
    K_B::R    = -4.43E-3/mV   # Voltage factor
end

"""
Electron trasnport chain proton flux, controlled by
 - mitochodnria NADH `nadh_m`
 - mitochodnria mebrane potential `Ψ`
"""
function j_hr(nadh_m, dpsi, p::ETC)
     @unpack VMAX, K_NADH, K_A, K_B = p
    fNADH = hill(nadh_m, K_NADH)
    fΨ =  (1 + K_A * dpsi) / (1 + K_B * dpsi)
    jHR = VMAX * fNADH * fΨ
    return jHR
end

(p::ETC)(nadh_m, dpsi) = j_hr(nadh_m, dpsi, p)

"""
    j_o2(jHRes)
Rate of oxygen consumption related to respiratory rate
"""
j_o2(jHRes, O2_H_RATIO=1//10) = O2_H_RATIO * jHRes

"""
Exponential proton leak

Reaction H(cyto) => H(mito)
"""
@with_kw struct HLeakExp{R}
    P_H::R = 2.4μM*Hz      # Basal proton leak rate
    K_V::R = 0.0305/mV     # Voltage factor
end

"Rate of proton leak"
j_hl(dpsi, p::HLeakExp) = p.P_H * exp(p.K_V * dpsi)

"Rate of proton leak"
(p::HLeakExp)(dpsi) = j_hl(dpsi, p)

"Proton leak with linear profile"
@with_kw struct HLeakLin{R}
    g::R = 0.500μM*Hz/mV     # Basal proton leak rate
end

"Rate of proton leak"
j_hl(dpsi, p::HLeakLin) = p.g * dpsi

"Rate of proton leak"
(p::HLeakLin)(dpsi) = j_hl(dpsi, p)

"""
F1Fo ATPase (ATP synthase) plus ANT

Reaction: 3ADP + 8H(cyto) => 3ATP + 8H(mito)
"""
@with_kw struct F1Fo{R}
    VMAX::R = 8mM*Hz        # Max rate of F1Fo ATPase (Adjustable)
    K_ADP::R = 20μM         # Michaelis constant for ADP
    K_V::R = 131.4mV        # Half activity Voltage (fitted for mammalian mitochondria data)
    K_CA::R = 0.165μM       # Activation constant for calcium
    F_MGADP::R = 0.055      # Ratio of MgADP: Total ADP in the cytosol
end

"Proton flux through ATP synthase"
function j_hf(adp, dpsi, ca_m, p::F1Fo)
    @unpack VMAX, K_V, K_CA, K_ADP, F_MGADP = p
    mgadp = adp * F_MGADP
    fADP = hill(mgadp, K_ADP, 2)
    fΨ = hill(dpsi, K_V, 8)
    fCa = 1 - exp(-ca_m / K_CA)
    jHF = VMAX * fADP * fΨ * fCa
    return jHF
end

(pF1Fo::F1Fo)(adp, dpsi, ca_m) = j_hf(adp, dpsi, ca_m, pF1Fo)

"ATP production rate by ATP synthase"
j_ant(jHf, ATP_H_RATIO = 1//3) = jHf * ATP_H_RATIO

"Mitochondrial calcium uniporter (MCU)"
@with_kw struct MCU{R}
    αm::R = 0.200          # Activity for mitochondrial calcium
    αc::R = 0.341          # Activity for cytosolic calcium
    P_CA::R = 4Hz          # Rate of MCU (adjustable, 2159Hz in Nguyen, 2007)
end

"Calcium flux via MCU"
function j_uni(ca_m, ca_c, ΔΨ, p::MCU)
    @unpack αm, αc, P_CA = p
    zvfrt = 2 * iVT * ΔΨ
    em1 = expm1(zvfrt)
    jUni = P_CA * exprel(zvfrt, em1) * (αc * ca_c * (1 + em1) - αm * ca_m)
    return jUni
end

"Calcium flux via MCU"
(p::MCU)(ca_m, ca_c, ΔΨ) = j_uni(ca_m, ca_c, ΔΨ, p)

"Inner function for NCLX rate"
function _nclx(ca_m, ca_c, A, P, VMAX, K_CA)
    B = ca_m / K_CA
    Q = ca_c / K_CA
    AB = A * B
    PQ = P * Q
    return VMAX * (AB - PQ) / (1 + A + B + P + Q + AB + PQ)
end

"Mitochondrial sodium calcium exchanger (NCLX), electroneutral"
@with_kw struct NCLX{R}
    VMAX::R = 75μM*Hz        # Max rate of NCLX (adjustable, 18.63mM/s in Nguyen, 2007)
    K_NA::R = 8.2mM          # Dissociation constant for Na (Dash and Beard , 2009)
    K_CA::R = 8μM            # Dissociation constant for Ca (Dash and Beard , 2009)
    na_c::R = 10mM           # Default concentration of cytosolic sodium
    na_m::R = 5mM            # Default concentration of mitochondrial sodium
    A::R = (na_c / K_NA)^2   # Precalculated cytosolic sodium term
    P::R = (na_m / K_NA)^2   # Precalculated mitochondrial sodium term
end

"Rate of sodium calcium exchanger (NCLX)"
j_nclx(ca_m, ca_c, p::NCLX) = _nclx(ca_m, ca_c, p.A, p.P, p.VMAX, p.K_CA)

(p::NCLX)(ca_m, ca_c, na_m, na_c) = _nclx(ca_m, ca_c, (na_c / p.K_NA)^2, (na_m / p.K_NA)^2, p.VMAX, p.K_CA)
(p::NCLX)(ca_m, ca_c) = _nclx(ca_m, ca_c, p.A, p.P, p.VMAX, p.K_CA)

"NADH shuttle"
@with_kw struct NADHT{R}
    VMAX::R = 50μM*Hz    # Max rate for NADH shuttle
    KTN_M::R = 16.78     # Factor for mitochondrial NAD/NADH ratio
    KTN_C::R = 0.002     # Factor for cytosolic NADH/NAD ratio
end

"Rate of NADH shuttle"
function j_nadht(nad_c, nadh_c, nad_m, nadh_m, p::NADHT)
    @unpack VMAX, KTN_C, KTN_M = p
    fc = hill(nadh_c, KTN_C * nad_c)
    fm = hill(nad_m, KTN_M * nadh_m)
    jTNADH = VMAX * fc * fm
    return jTNADH
end

(p::NADHT)(nad_c, nadh_c, nad_m, nadh_m) = j_nadht(nad_c, nadh_c, nad_m, nadh_m, p)
