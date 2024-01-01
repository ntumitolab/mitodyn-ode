# Adenylate kinase (AdK) equilibrium w.r.t adenylate energy charge (AEC)

import NaNMath as nm
using PythonPlot

function aec2axp(a, keq=1)
    x = inv(8 * keq - 2) * (-1 + nm.sqrt(1 + 4 * (4keq - 1) * a * (1 - a)))
    return (atp=a - x, adp=2x, amp=1 - a - x)
end

aa = 0.0:0.01:1.0

axp = aec2axp.(aa, 0.91)

atp = getindex.(axp, :atp)
adp = getindex.(axp, :adp)
amp = getindex.(axp, :amp)

figure()
plot(aa, atp, label="ATP")
plot(aa, adp, label="ADP")
plot(aa, amp, label="AMP")
xlabel("Energy charge")
ylabel("Proportion")
legend()
grid()
gcf()
