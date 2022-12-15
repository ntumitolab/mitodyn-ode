##################################
### Units and physical constants
##################################
const second = float(1)    # second
const minute = 60second    # minute
const ms = 1e-3second      # millisecond
const Hz = inv(second)     # Herz
const kHz = 1000Hz         # kilohertz
const metre = float(1)     # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1000mL        # liter
const μL = 1E-6Liter
const pL = 1E-12Liter
const mM = float(1)
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = 1E-3mM          # micromolar
const nM = 1E-6mM          # nanomolar
const Volt = float(1)      # volt
const mV = 1E-3Volt        # millivolt
const T0 = 310.0           # Default temperature
const F = 96485.0          # Faraday constant (c / mol)
const R = 8.314            # Ideal gas constant (K/mol)
const VT = R * T0 / F      # Default thermal voltage (Volts)
const iVT = inv(VT)        # Reciprocal of thermal voltage

# Model Constants
const C_MIT = 1.812μM / mV        # Mitochondrial membrane capacitance
const F_M = 3E-4                # Frasction of free Ca in mitochondria
const F_I = 0.01                # Fraction of free Ca in cytosol
const V_I = 0.53                # Relative cytoplasmic volume
const V_MT = 0.06               # Relative mitochondrial volume
const V_MTX = 0.0144            # Relative mitochondrial matrix volume
const iVi = inv(V_I)
const iVmtx = inv(V_MTX)
const iVimtx = inv(V_MTX + V_I)
const iCmt = inv(C_MIT)

##################################
### Commonly-used functions
##################################

"""
Regular Hill function
"""
hill(x, k=one(x)) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

"""
Repressive Hill function
"""
hillr(x, k=one(x)) = hill(k, x)
hillr(x, k, n) = hill(k, x, n)

"""
    expit(x)

Logistic sigmoid function. `expit(x) = 1 / (1 + exp(-x))`
"""
expit(x) = hillr(exp(-x))

"""
    exprel(x, em1 = expm1(x))

Returns `x / (exp(x)-1)`
"""
exprel(x, em1=expm1(x)) = x / em1

##################################
### Reaction rates
##################################

"""
Rate of Adenylate kinase (AdK):

2ADP <=> ATP + AMP
"""
function j_adk(atp_c, adp_c, amp_c, kf, keq)
    return kf * (adp_c * adp_c - atp_c * amp_c * keq)
end

"""
Rate of Glucokinase (GK) and upper glycolysis reactions:

- Glc + ATP => G6P (+ ADP)
- G6P + ATP => FBP (+ ADP) =>> 2G3P (+ ADP)
"""
function j_gk(atp_c, glc, vmax, katp, kglc, n)
    return vmax * hill(atp_c, katp) * hill(glc, kglc, n)
end

"""
Rate of Glyceraldehydes 3-phosphate dehydrogenase (GPD), lumped with pyruvate kinase (PK) to represent lower half of glycolysis

G3P + NAD + 2ADP => Pyr + NADH + 2ATP
"""
function j_gpd(g3p, nad_c, nadh_c, adp_c, vmax, kg3p, knad, kadp)
    return vmax * hill(adp_c, kadp) * hill(nad_c, nadh_c * knad) * hill(g3p, kg3p)
end

"""
Rate of Lactate dehydrogenase (LDH)
The rate of lactate output is approximately 5% of the glycolytic flux @ glucose = 8mM.

Pyr + NADH => Lactate + NAD
"""
function j_ldh(pyr, nad_c, nadh_c, vmax, kpyr, knadh)
    return vmax * hill(pyr, kpyr) * hill(nadh_c, nad_c * knadh)
end

"""
Steady-state cytosolic calcium based on ATP:ADP ratio
"""
function ca_c_atp(atp_c, adp_c, resting, amplitude, n, katp)
    return resting + amplitude * hill(atp_c, adp_c * katp, n)
end

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

"""
Rate of Electron trasnport chain (ETC)

NADH + O2 + 10H(mito) => NAD + 10H(cyto)
"""
function j_hr(nadh_m, dpsi, vmax, knadh, ka, kb)
    return vmax * hill(nadh_m, knadh) * (1 + ka * dpsi) / (1 + kb * dpsi)
end

"""
Proton leak
"""
function j_hl(dpsi, phl, kv)
    return phl * exp(kv * dpsi)
end

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

"""
Mitochondrial calcium uniporter (MCU)

Ca(cyto) <=> Ca(mito)
"""
function j_uni(ca_m, ca_c, ΔΨ, P_CA)
    zvfrt = 2 * iVT * ΔΨ
    em1 = expm1(zvfrt)
    return P_CA * exprel(zvfrt, em1) * (0.341 * ca_c * (1 + em1) - 0.200 * ca_m)
end

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

"""
NADH shuttle

NADH(cyto) + NAD(mito) => NADH(mito) + NAD(cyto)
"""
function j_nadht(nad_c, nadh_c, nad_m, nadh_m, vmax, ktn_c, ktn_m)
    return vmax * hill(nadh_c, nad_c * ktn_c) * hill(nad_m, nadh_m * ktn_m)
end

"""
AMPK activity
"""
function a_ampk(atp_c, amp_c, kampk)
    return hill(amp_c / atp_c, kampk)
end
