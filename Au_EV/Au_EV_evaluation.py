import DebyeTemperature as dt
import numpy as np

# Option 1:
# evaluate via fitting to Es_cohesive vs Wigner-Seitz radius
# datafile = "Au_EV_#Job8510.txt"
# xs, ys = np.loadtxt(datafile, usecols=(0, 1), comments='$', delimiter=' ', unpack=True)
#
# WSRs, Es_total, Es_Cohesive = [], [], []
# for x, y in zip(xs, ys):
#     if x != 10:
#         WSRs.append(dt.A2WSr4Cubic(x))
#         Es_total.append(y)
#     else:
#         E_free = y
# for E in Es_total:
#     Es_Cohesive.append(E - E_free)
#
# dt.Evaluation('Au', WSRs, Es_Cohesive)

# Option 2:
# quick calculation using data from Web
PoissonRatio = 0.42
BulkModulus = 220            # GPa
qt = dt.Structure('Au')
k = dt.EvaluateCoeffK(PoissonRatio)
r0 = 1.62                    # Ang
BulkModulus_Pa = BulkModulus * 1e9     # Pa
DebyeT = dt.Calc_DebyeT(k, r0, BulkModulus_Pa, qt.MolarWeight)

print("calculated Debye temperature: ", DebyeT)
print("Debye temperature from web: ", 162)