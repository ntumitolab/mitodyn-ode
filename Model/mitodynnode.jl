using Parameters
using Setfield
using LabelledArrays
using FromFile
@from "utils.jl" import mM, μM, minute, Hz, mV, hill

const C_MIT = 1.812μM/mV        # Mitochondrial membrane capacitance
const F_M = 3E-4                # Frasction of free Ca in mitochondria
const F_I = 0.01                # Fraction of free Ca in cytosol
const V_I = 0.53                # Relative cytoplasmic volume
const V_MT = 0.06               # Relative mitochondrial volume
const V_MTX = 0.0144            # Relative mitochondrial matrix volume
const iVi = inv(V_I)
const iVmtx = inv(V_MTX)
const iVimtx = inv(V_MTX + V_I)
const iCmt = inv(C_MIT)

"""
Parameters for Mitochodrial dynamics (Node model)
"""
@with_kw struct MitoDynNode{R}
    # Fission-fusion rates
    K_FISS::R = inv(10minute)    # Basal fission rate
    K_FUSE::R = K_FISS           # Basal fusion rate
    R_FUSE23::R = 0.5            # Ratio of tip-to-side to tip-to-tip fusion rate
    R_FISS23::R = 1.5            # Ratio of X3 to X2 fission rate
    # Cytosol Parameters
    k_nadhc::R = 0.1Hz           # Baseline cytoplasm NADH consumption rate
    k_atp::R = 0.04Hz            # Baseline cytoplasm ATP hydrolysis rate
    k_atp_ca::R = 90Hz/mM        # Cytoplasm ATP hydrolysis rate activated by calcium
    k_g3p::R = 0.01Hz            # Baseline G3P consumption rate
    k_pyr::R = 0.01Hz            # Baseline Pyr consumption rate
    Σa_c::R = 4.5mM              # Free adenylate nucleotides pool (cytosolic), adjustable
    Σn_c::R = 2mM                # Free pyridine nucleotides pool (cytosolic),
    gk = GK(5mM)
    gpd = GPD()
    ldh = LDH()
    cai = CaATP()
    adk = AdK()

    # Mitochondrial Parameters
    k_nadhm::R = 0.1Hz            # Baseline mitochondrial NADH consumption rate
    Σn_m::R = 2.2mM               # Free pyridine mitochondrial nucleotides pool
	f1fo = F1Fo()
    etc = ETC()
    hleak = HLeakExp()
	pdh = PDH()
    nadht = NADHT()
    mcu = MCU()
    nclx = NCLX()
end

"Get cytosolic AMP"
ampcyto(adp_c, atp_c, p::MitoDynNode) = p.Σa_c - atp_c - adp_c
ampcyto(u, p::MitoDynNode) = ampcyto(u.adp_c, u.atp_c, p.cyto)

"Set glucose concentration"
setglc(p::MitoDynNode, glc) = MitoDynNode(p; gk = setglc(p.gk, glc))

"Get cytosolic calcium"
cacyto(adp_c, atp_c, p::MitoDynNode, t) = cacyto(adp_c, atp_c, p.cai, t)
cacyto(u, p::MitoDynNode, t) = cacyto(u.adp_c, u.atp_c, p::MitoDynNode, t)

"Degree one node by conservation"
getx1(x2, x3) = 1 - 2x2 - 3x3
getx1(u) = getx1(u.x2, u.x3)

"Average degree"
avgdeg(x2, x3, x1=getx1(x2, x3)) = (x1 + 2x2 + 3x3) / (x1 + x2 + x3)
avgdeg(u) = avgdeg(u.x2, u.x3)

"""
    model!(du, u, p::MitoDynNode, t)

ODE systems of the beta cell model.

1  Adenylate kinase reaction for 2ADP = AMP + ATP.
2. H:ATP ratio adjusted to 3 since jHf and jANT were separated
3. NCLX math description returned to Nguyen's
4. Mitochondrial dynamics with nodes
"""
function model!(du, u, p::MitoDynNode, t)
    @unpack g3p, pyr, nadh_m, nadh_c, atp_c, adp_c, ca_m, dpsi, x2, x3 = u
    # Conservative relationships
    @unpack Σn_m, Σa_c, Σn_c = p
    amp_c = Σa_c - adp_c - atp_c
    nad_c = Σn_c - nadh_c
    nad_m = Σn_m - nadh_m
    x1 = getx1(x2, x3)

    # Cytoslic rates
    @unpack cai, gk, gpd, ldh, adk = p
    ca_c = cai(adp_c, atp_c, t)
    jGK = gk(atp_c)
    jGPD = gpd(g3p, nad_c, nadh_c, adp_c)
    jLDH = ldh(pyr, nad_c, nadh_c)
    jADK = adk(amp_c, adp_c, atp_c)

    # Mitochondrial rates
    @unpack f1fo, pdh, nadht, etc, hleak, mcu, nclx = p
    jHf = f1fo(adp_c, dpsi, ca_m)
    jPDH = pdh(pyr, nad_m, nadh_m, ca_m)
    jTNADH = nadht(nad_c, nadh_c, nad_m, nadh_m)
    jHr = etc(nadh_m, dpsi)
    jHL = hleak(dpsi)
    jUni = mcu(ca_m, ca_c, dpsi)
    jNCLX = nclx(ca_m, ca_c)
    jANT = j_ant(jHf)
    jO2 = j_o2(jHr)

    # Fission-fusion rates
    @unpack K_FUSE, K_FISS, R_FUSE23, R_FISS23 = p
    kFuse1 = K_FUSE * jANT / jHL
    kFuse2 = R_FUSE23 * kFuse1
    kFiss1 = K_FISS
    kFiss2 = R_FISS23 * kFiss1

    v1 = kFuse1 * x1 * x1 - kFiss1 * x2  # Rate: X1 -> X2
    v2 = kFuse2 * x1 * x2 - kFiss2 * x3  # Rate: X2 -> X3

    # ODEs
    @unpack k_g3p, k_pyr, k_nadhc, k_nadhm, k_atp, k_atp_ca = p
    du.nadh_m = iVmtx * (46//10 * jPDH + jTNADH - jO2) - k_nadhm * nadh_m
    du.dpsi = iCmt * (jHr - jHf - jHL - jANT - 2jUni)
    du.ca_m = iVmtx * F_M * (jUni - jNCLX)
    du.g3p = iVi * (2jGK - jGPD) - k_g3p * g3p
    du.pyr = iVimtx * (jGPD - jPDH - jLDH) - k_pyr * pyr
    du.nadh_c = iVi * (jGPD - jLDH - jTNADH) - k_nadhc * nadh_c
    du.atp_c = iVi * (-2jGK + 2jGPD + jANT + jADK) - (k_atp + k_atp_ca * ca_c) * atp_c
    du.adp_c = -du.atp_c - iVi * 2jADK
    du.x2 = v1 - v2
    du.x3 = v2

    return du
end
