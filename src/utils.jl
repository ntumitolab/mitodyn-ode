using PythonCall

## Units and physical constants
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

"""Hill and MM functions"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(nm.pow(x, n), nm.pow(k, n))

"""Get the index of a symbol from an array of symbols"""
indexof(sym, syms) = findfirst(isequal(sym), syms)

"""Extract values from ensemble simulations by a symbol"""
extract(sims, k) = getindex.(sims.u, k)

"""Change parameter(s) from an ODESystem and returns a parameter vector."""
function change_params(sys, ps...)
    ModelingToolkit.varmap_to_vars(Dict(ps), parameters(sys); defaults = sys.defaults)
end

"""Export publication-ready TIFF file from a figure"""
exportTIF(fig, name; dpi=300) = fig.savefig(name, dpi=dpi, pil_kwargs=pydict(Dict("compression" => "tiff_lzw")))
