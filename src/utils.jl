# Units
const ms = 1                    # millisecond
const μM = 1                    # micromolar = 1
const mV = 1                    # millivolt
const Kelvin = 1                # temperature unit Kelvin
const mmol = 1                  # millimole = 1
const second = 1000ms           # second is the SI unit
const minute = 60second         # minute
const Hz = inv(second)          # Herz
const kHz = inv(ms)             # kilohertz
const m³ = mmol / μM            # cubic meter
const metre = cbrt(m³)          # meter
const cm = 0.01metre            # centimeter
const cm² = cm^2                # square centimeter
const μm = 1e-6metre            # micrometer
const mL = cm^3                 # milliliter = cubic centimeter
const Liter = 1000mL            # liter
const μL = μm^3                 # microliter
const pL = 1e-12Liter           # picoliter
const mol = 1000mmol            # mole
const mM = 1000μM               # mM is the SI unit
const Molar = 1000mM            # molarity is used in equilibrium constants
const nM = 0.001μM              # nanomolar
const μFcm⁻² = 1                # area capacitance (μF/cm²)
const μF = μFcm⁻² * cm²         # microfarad
const Farad = 1e6μF             # Farad
const μAμF = mV * inv(ms)       # common current density
const μA = μAμF * μF            # micropampere
const μAcm⁻² = μAμF * μFcm⁻²    # real current density
const Ampere = 1e6μA            # electric current unit Ampere
const Columb = Ampere * second  # electric charge unit Columb
const Volt = 1000mV             # electric potential unit Volt
const Joule = Columb * Volt     # energy unit Joule
const Seimens = Ampere / Volt   # conductance unit
const milliseimens = 0.001Seimens # milliseimens
const mScm⁻² = milliseimens / cm²
const mSμF = μAμF / mV          # conductance density
const Faraday = 96485Columb / mol # Faraday constant (columb / mol)
const T₀ = 310Kelvin            # Default temp (37C)
const RGAS = 8.314Joule / Kelvin / mol # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ / Faraday  # Thermal voltage (@37C), 26.7 mV
const iVT = inv(VT)             # Reciprocal of thermal voltage (0.037 per mV)

"""Hill and Michaelis-Menten functions"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(nm.pow(x, n), nm.pow(k, n))

"""Get the index of a symbol from an array of symbols"""
indexof(sym, syms) = findfirst(isequal(sym), syms)

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s->s[k], sim)

"""Export publication-ready TIFF file from a figure"""
exportTIF(fig, name; dpi=300) = fig.savefig(name, dpi=dpi, pil_kwargs=pydict("compression" => "tiff_lzw"))

"""
The deviation of adynylate pool steady state from energy charge
"""
aecdev(a, keq=one(a)) = (-1 + nm.sqrt(1 + 4 * (4keq - 1) * a * (1 - a))) / (8 * keq - 2)

"""
ATP proportion from energy charge
"""
aec2atp(a, keq=one(a)) = a - aecdev(a, keq)

"""
ADP proportion from energy charge
"""
aec2adp(a, keq=one(a)) = 2 * aecdev(a, keq)
