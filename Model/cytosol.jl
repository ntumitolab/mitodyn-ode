# Cytosolic compartment: glycolysis, cytosolic calcium relationaship, and NADH transport.
# Adapted from Fridlyand, 2010
using Parameters
using Setfield
using FromFile
@from "utils.jl" import hill, mM, μM, ms, minute, Hz

"""
Glucokinase (GK)

Reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP => 2G3P
"""
@with_kw struct GK{R}
    VMAX::R = 0.011mM * Hz  # Max rate for glucokinase
    K_ATP::R = 0.5mM        # Apparent Michaelis constant for ATP
    K_GLC::R = 7mM          # Apparent Michaelis constant for glucose
    N::R = 1.7              # Apparent Hill cosefficient for glucose
    VMAX_GLC::R = VMAX      # Max rate for glucokinase, adjusted for glucose
end

"GK flux under customized glucose concentration"
function GK(glc; VMAX=0.011mM * Hz, K_ATP = 0.5mM, K_GLC = 7mM, N = 1.7)
    GK(VMAX, K_ATP, K_GLC, N, VMAX * hill(glc, K_GLC, N))
end

"Reaction rate of Glucokinase"
j_gk(atp, glc, p::GK) = p.VMAX * hill(atp, p.K_ATP) * hill(glc, p.K_GLC, p.N)
j_gk(atp, p::GK) = p.VMAX_GLC * hill(atp, p.K_ATP)
(p::GK)(atp) = j_gk(atp, p)
(p::GK)(atp, glc) = j_gk(atp, glc, p)

"Set glucose level for Glucokinase"
setglc(glc, p::GK) = GK(glc, VMAX=p.VMAX, K_ATP=p.K_ATP, K_GLC=p.K_GLC, N=p.N)
setglc(p::GK, glc) = setglc(glc, p)

"""
Glyceraldehydes 3-phosphate dehydrogenase (GPD)

Reactions lumped with pyruvate kinase (PK)

G3P + NAD + 2ADP => Pyr + NADH + 2ATP
"""
@with_kw struct GPD{R}
    VMAX::R = 0.5mM * Hz  # Max rate for GPD
    K_G3P::R = 0.2mM      # Apparent Michaelis constant for G3P
    K_NAD::R = 0.09       # Apparent Michaelis constant for NAD/NADH ratio
    K_ADP::R = 0.020mM    # Apparent Michaelis constant for ADP
end

"Rate of Glyceraldehydes 3-phosphate dehydrogenase (GPD)"
function j_gpd(g3p, nad_c, nadh_c, adp_c, p::GPD)
    @unpack VMAX, K_G3P, K_NAD, K_ADP = p
    fADP = hill(adp_c, K_ADP)
    fNAD = hill(nad_c, K_NAD * nadh_c)
    fG3P = hill(g3p, K_G3P)
    jGPD = VMAX * fG3P * fNAD * fADP
    return jGPD
end

(p::GPD)(g3p, nad_c, nadh_c, adp_c) = j_gpd(g3p, nad_c, nadh_c, adp_c, p)


"""
Lactate dehydrogenase (LDH)

The rate of lactate output is approximately 5% of the glycolytic flux at 8 mM glucose.

Reaction:

Pyr + NADH => Lactate + NAD
"""
@with_kw struct LDH{R}
    VMAX::R = 1.2mM * Hz  # Max rate for LDH (Ad.)
    K_PYR::R = 47.5μM     # Apparent Michaelis constant for Pyruvate
    K_NADH::R = 1.0       # Apparent Michaelis constant for NADH/NAD ratio
end

"Rate of lactate dehydrogenase (LDH)"
j_ldh(pyr, nad_c, nadh_c, p::LDH) =  p.VMAX * hill(pyr, p.K_PYR) * hill(nadh_c, p.K_NADH * nad_c)

"Rate of lactate dehydrogenase (LDH)"
(p::LDH)(pyr, nad_c, nadh_c) = j_ldh(pyr, nad_c, nadh_c, p)

"Parameters for empirical cytsolic calcium dependent on ATP"
@with_kw struct CaATP{R, N}
    ca_r::R = 0.09μM        # Resting cytosolic Ca concentration
    ka_ca::R = 0.25μM       # ATP-activated calcium amplitude
    n::N = 4                # Hill coefficient for ATP/ADP ratio
    k_atp::R = 25.0         # ATP/ADP ratio @ Half maximal calcium
end

"Cytosolic calcium controlled by ATP/ADP ratio"
cacyto(adp_c, atp_c, p::CaATP) = p.ca_r + p.ka_ca * hill(atp_c, adp_c * p.k_atp, p.n)
cacyto(adp_c, atp_c, p::CaATP, t) = cacyto(adp_c, atp_c, p)
(p::CaATP)(adp_c, atp_c) = cacyto(adp_c, atp_c, p)
(p::CaATP)(adp_c, atp_c, t) = cacyto(adp_c, atp_c, p, t)

"Cytosolic calcium sinusoid oscillator"
@with_kw struct CaSine{R}
    ca_r::R   = 0.09μM           # Resting cytosolic Ca concentration
    ka_ca::R  = 0.25μM           # Calcium amplitude
    period::R = 2minute          # Calcium oscillation period
    f2::R     = 2 * inv(period)  # frequency
    a::R      = 0.5ka_ca         # Rescaled amplitude
end

"Cytosolic calcium with sinusoid waves"
cacyto(p::CaSine, t) = p.ca_r + p.a * (sinpi(p.f2 * t) + 1)
cacyto(adp_c, atp_c, p::CaSine, t) = cacyto(p, t)
(p::CaSine)(t) = cacyto(t, p)
(p::CaSine)(adp_c, atp_c, t) = p(t)

"Cytosolic calcium oscillator with smooth curves"
@with_kw struct CaOsciSmooth{R, T, S}
    ca_r::R    = 0.09μM      # Resting cytosolic Ca concentration
    ka_ca::R   = 0.25μM      # Calcium amplitude
    period::R  = 2minute     # Calcium oscillation period
    assym::T   = 5           # Assymetric factor
    steep::S   = 4           # Steepess factor
end

"Cytosolic calcium oscillation with asymmetric curves"
function cacyto(p::CaOsciSmooth, t)
    @unpack ca_r, period, ka_ca, assym, steep = p
    tau = t / period
    x = assym * (tau - floor(tau))
    return ca_r + ka_ca * (x * exp(1 - x))^steep
end

"Cytosolic calcium oscillation with asymmetric smooth curves"
cacyto(adp_c, atp_c, p::CaOsciSmooth, t) = cacyto(p, t)
(p::CaOsciSmooth)(t) = cacyto(p, t)
(p::CaOsciSmooth)(adp_c, atp_c, t) = p(t)

"""
Adenylate kinase

Reaction:

2ADP = ATP + AMP
"""
@with_kw struct AdK{R}
    KF::R  = 1000Hz/mM    # Max forward rate constant for AdK. It should be much faster than other cytosolic reactions.
    KEQ::R = 0.931        # Equilibrium constant (AMP-forming). (Goldings, 1995)
    KB::R  = KF / KEQ     # Max backward rate
end

"Rate of adenylate kinase"
j_adk(amp_c, adp_c, atp_c, p::AdK) = p.KF * adp_c * adp_c - p.KB * atp_c * amp_c
(p::AdK)(amp_c, adp_c, atp_c) = j_adk(amp_c, adp_c, atp_c, p)

@with_kw struct AMPK{R}
    K_AMP = 0.07        # Michaelis constant for AMP / ATP ratio
end

"Activity of AMPK"
k_ampk(amp, atp, p::AMPK) = hill(amp / atp, p.K_AMP)
(p::AMPK)(amp, atp) = k_ampk(amp, atp, p)
