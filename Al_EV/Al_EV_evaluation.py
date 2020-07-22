import DebyeTemperature as dt
import numpy as np

datafile = "Al_EV_#Job8502.txt"
xs, ys = np.loadtxt(datafile, usecols=(0, 1), comments='$', delimiter=' ', unpack=True)

WSRs, Es_total, Es_Cohesive = [], [], []
for x, y in zip(xs, ys):
    if x != 10:
        WSRs.append(dt.A2WSr4Cubic(x))
        Es_total.append(y)
    else:
        E_free = y
for E in Es_total:
    Es_Cohesive.append(E - E_free)

dt.Evaluation('Al', WSRs, Es_Cohesive)
