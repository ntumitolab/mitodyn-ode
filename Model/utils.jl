# Units
const second = float(1)    # second
const minute = 60second    # minute
const ms = 1e-3second      # millisecond
const Hz = inv(second)     # Herz
const kHz = 1e3Hz          # kilohertz
const metre = float(1)     # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1e3mL        # liter
const μL = 1E-6Liter
const pL = 1E-12Liter
const mM = float(1)
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = 1E-3mM          # micromolar
const nM = 1E-6mM          # nanomolar
const Amp = float(1)       # ampere
const mA = 1E-3Amp         # milliampere
const μA = 1E-6Amp         # micrpampere
const Volt = float(1)      # volt
const mV = 1E-3Volt        # millivolt
const T0 = 310.            # Default temp (37C)
const F = 96485.           # Faraday constant (c / mol)
const R = 8.314            # Ideal gas constant (K/mol)
const VT = R * T0 / F      # Thermal voltage (@37C) (Volts)
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)

##################################
### Commonly-used functions
##################################

"Return `1+x` with the same type as x"
p_one(x) = one(x) + x

"Return `1-x` with the same type as x"
one_m(x) = one(x) - x

"""
Regular Hill function
"""
hill(x, k = one(x)) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

"""
Repressive Hill function
"""
hillr(x, k = one(x)) = hill(k , x)
hillr(x, k, n) = hill(k, x, n)

"""
Logistic sigmoid function.
See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expit.html
"""
expit(x) = hillr(exp(-x))

"""
    exprel(x, em1 = expm1(x))

Returns `x / (exp(x)-1)` accurately when x is near zero.
See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exprel.html
Note the fraction is the inverse of `scipy.exprel()`
"""
function exprel(x, em1 = expm1(x))
    res = x / em1
    return ifelse(x ≈ zero(x), one(res), res)
end

"""
Signed sqare root
"""
sqrt_s(x) = flipsign( sqrt(abs(x)), x)

"""
Signed power
"""
pow_s(x, n) = flipsign( abs(x)^n, x)

"""
Signed Hill function
"""
hill_s(x, k, n) = hill(pow_s(x, n), pow_s(k, n))
