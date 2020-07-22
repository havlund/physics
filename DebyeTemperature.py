import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, fabs, sqrt
from scipy.misc import derivative
from lmfit import Minimizer, report_fit, Parameters
import MolarWeight
import json
import os.path

"""
--------------------------------------------------------------------------------------------------
# April 18th, 2019     "DebyeTemperature.py"
--------------------------------------------------------------------------------------------------
# Part 1:
# evaluate from "cohesive" energy or "binding" energy:
    * Bulk modulus
    * Debye temperature

    # evaluation procedure
    *1) prepare data of Wigner-Seitz Length (r) vs Cohesive energy (E)
    *2) fitting to arbitrary polynominal expression
    *3) determine the equilibrium state, i.e. r0 and E0
    *4) with r0 & E0, fitting to empirical model, determine Length Scale (ls)
    *5) evaluate Bulk Modulus (B) with ls & E0
    *6) evaluate Debye Temperature (DT) with B and guessed coefficient k
    Optionally,
    * step 2)-3): one can start with 3) using r0 & E0 from DFT minimizations
    * step 2)-4): one can start with 3) fitting r0, E0, and ls at one time
--------------------------------------------------------------------------------------------------
# Part 2:
# evluate from specific elastic properties, e.g. C11, C12...
    * Bulk modulus, better than that in Part 1
    * Poisson ratio
    * coefficient k, to be used for step 6) of Part 1

    # This is basically what MedeA uses to predict Bulk modulus, Debye temperature...
--------------------------------------------------------------------------------------------------
"""

k_Boltzm = 1.38064852e-23  # J/K
h_planck = 6.62607004e-34  # J.s
N_Avogadro = 6.022e23
current_folder = os.path.dirname(os.path.realpath(__file__))

# --------------------------------------------------------------------------------------------------
# Part 1. Class
#      1.1 Definition
class Structure:
    def __init__(self, formula, StruType = 'A3'):
        self._type = StruType
        self.MolarWeight = MolarWeight.Lookup_MolarWeight(formula)['AtWt'] # g/mole

#      1.2 experimental data
    def import_exp_data(self, expdata):     # data of 2 lists, E(J) and r(m)
        self._rexp = expdata[0]
        self._Eexp = expdata[1]
    def get_exp_data(self):
        return self._rexp, self._Eexp

#      1.3 empirical function and fitting function
    def E_total(self, r, ls):                         # function describing E_V
        a = (r - self.r0) / ls
        E0 = fabs(self.E0)
        return - E0 * (1 + a + 0.05 * a ** 3) * exp(-a)

    def E_total_fitting_Er(self):                 # fitting DFT data of E_V with the above function
        def fcn2min(prms, rs, Es):
            ls = prms['ls']
            E_calc = []
            for ri in rs:
                E_calc.append(self.E_total(ri, ls))
            return np.array(E_calc) - Es

        params = Parameters()
        _last_fitting = "/_last_fitting.json"
        filename = current_folder + _last_fitting
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                dct = json.load(f)
            if 'ls' in dct:
                params.add('ls', dct['ls'], vary=True)
        else:
            params.add('ls', 0.5, vary=True)
        minner = Minimizer(fcn2min, params=params, fcn_args=(self._rexp, self._Eexp))
        result = minner.minimize()
        final = self._Eexp + result.residual
        _dct = {'ls': result.params['ls'].value}
        with open(filename, 'w') as f:
            json.dump(_dct, f)
        report_fit(result)
        self.LengthScale = result.params['ls'].value
        return final

#      1.4 polynominal and its fitting function
    def Poly_Er(self, r):
        A = self.__A
        B = self.__B
        C = self.__C
        D = self.__D
        return A + B * r + C * r **2 + D * r **3

    def Poly_fitting_Er(self):
        def fcn2min(prms, rs, Es):
            A = prms['A']
            B = prms['B']
            C = prms['C']
            D = prms['D']
            E_calc = []
            E_poly = lambda r: A + B * r + C * r **2 + D * r **3
            for ri in rs:
                E_calc.append(E_poly(ri))
            return np.array(E_calc) - Es
        params = Parameters()
        params.add('A', 0, vary=True)
        params.add('B', 0, vary=True)
        params.add('C', 0, vary=True)
        params.add('D', 0, vary=True)
        minner = Minimizer(fcn2min, params=params, fcn_args=(self._rexp, self._Eexp))
        result = minner.minimize()
        final = self._Eexp + result.residual
        report_fit(result)
        self.__A = result.params['A'].value
        self.__B = result.params['B'].value
        self.__C = result.params['C'].value
        self.__D = result.params['D'].value
        return final

#   1.5. Determination of E0 and r0

# in case where data are available, they are set directly
    def set_equilibrium(self, E0, r0):
        self.E0 = E0     # equilibrium total energy
        self.r0 = r0     # equilibrium lattice parameter

# otherwise, they can be determined as the minimum of the fitting curve, here polymoninal description is used.
    def Eval_Equilibrium_Poly(self, start=1.5, accu = 1e-3):            # 1e-10 m, with an accuracy of 0.0001 nm
        done = False
        f = lambda r: derivative(self.Poly_Er, r, dx = 0.01)
        r = start
        step = accu * 1e2                                          # 1e-10 m, 0.01 nm @ accuracy = 0.0001 nm
        if f(start) > 0:
            while not done:
                r -= step
                if f(r) < 0:
                    done = True
                    intv = (r, r+step)
                    ffabmin = fabs(f(r))
                    r0 = r
        else:
            while not done:
                r += step
                if f(r) > 0:
                    done = True
                    intv = (r-step, r)
                    ffabmin = fabs(f(r))
                    r0 = r
        grids = np.linspace(intv[0], intv[1], 100)            # 100 grids in the interval of 0.01 nm
        for r in grids:
            fi = f(r)
            fifab = fabs(fi)
            if fifab < ffabmin:
                ffabmin = fifab
                r0 = r
        E0 = self.Poly_Er(r0)
        self.E0 = E0                                  # E0 in J
        self.r0 = r0                                  # r0 in m
        print("\nEvaluated equilibrium: \nr0 = {:.4f} \na0 = {:.4f} \nE0 = {:.4f}".format(r0, WSr2a4Cubic(r0), E0))

# --------------------------------------------------------------------------------------------------
# Part 2. evaluate from elastic properties
#   * Bulk modulus
#   * Poisson ratio
#   * coefficient k

def VoigtReussHill(type, C11, C12, C13, C33, C44):

    def Voigt_bound(type, C11, C12, C13, C33, C44):
        if type == 'A1':
            Bv = (C11 + 2 * C12) / 3
            Sv = (C11 - C12 + 3 * C44) / 5
            Vv = (3 * Bv - 2 * Sv) / (6 * Bv + 2 * Sv)
        elif type == 'A3':
            Bv = (2 * (C11 + C12) + C33 + 4 * C13) / 9
            Sv = (6 * (C11 - C12) + 12 * C44 + (C11 + C12 + 2 * C33 - 4 * C13)) / 30
            Vv = (3 * Bv - 2 * Sv) / (6 * Bv + 2 * Sv)
        return Bv, Sv, Vv

    def Reuss_bound(type, Bv, C11, C12, C13, C33, C44):
        if type == 'A1':
            Br = (C11 + 2 * C12) / 3
            Sr = 5 * (C11 - C12) * C44 / (4 * C44 + 3 * (C11 - C12))
            Vr = (3 * Br - 2 * Sr) / (6 * Br + 2 * Sr)
        elif type == 'A3':
            Csqr = ((C11 + C12) * C33 - 2 * C13 ** 2)
            Br = Csqr / (C11 + C12 + 2 * C33 - 4 * C13)
            Sr = 5 / 2 * (Csqr * C44) / (3 * Bv * C44 + Csqr * (2 * C44 / (C11 - C12) + 1))
            Vr = (3 * Br - 2 * Sr) / (6 * Br + 2 * Sr)
        return Br, Sr, Vr

    Bv, Sv, Vv = Voigt_bound(type, C11, C12, C13, C33, C44)
    Br, Sr, Vr = Reuss_bound(type, Bv, C11, C12, C13, C33, C44)
    BulkModulus_vrh = (Bv + Br) / 2
    PoissonRatio_vrh = (Vv + Vr) / 2
    return BulkModulus_vrh, PoissonRatio_vrh

def EvaluateCoeffK(v):
    k = ((((1+v)/(3*(1-v)))**(3/2) + 2*(2*(1+v)/(3*(1-2*v)))**(3/2)) / 3) ** (-1/3)
    return k


# ------------------------------------------------------------------------------------------------------------------
# Part 3. Calculate Debye temperature and Bulk modulus
def Calc_DebyeT(k, r0, B, Mw):
    fct = k * (h_planck/(2*pi)) / k_Boltzm * (48 * pi**5)**(1/6)
    r0 = r0 * 1e-10                  # m
    term = sqrt(r0 * B * N_Avogadro / (Mw * 1e-3))
    DebyeT = round(fct * term, 4)    # K
    return DebyeT

def Calc_BulkModulus(r0, E0, ls, unit=None):    # Emperical model, in which model coefficients depend on the units of parameters.
    E0 = fabs(E0) * 1e-19            # J
    ls = ls * 1e-10                  # m
    r0 = r0 * 1e-10                  # m
    B = E0 / (12 * pi * r0 * ls**2)   # in Pa
    B = round(B, -6)
    if unit == 'Pa' or unit == None:
        return B
    elif unit == 'GPa':
        return B * 1e-9
    elif unit == 'Bar':
        return B * 1e-5
    elif unit == 'MBar':
        return B * 1e-11

# ------------------------------------------------------------------------------------------------------------------
# Part 4. data process and fitting procedure

def A2WSr4Cubic(a):
    # Lattice constant >> Wigner_Seitz cell length
    # 4 * 4/3 * pi * r**3 = a**3
    return a / (16*pi/3) ** (1/3)

def WSr2a4Cubic(r):
    # Wigner_Seitz cell length >> Lattice constant
    # 4 * 4/3 * pi * r**3 = a**3
    return r * (16*pi/3) ** (1/3)

def WSRs4recalculation(WSRs, f_lower, f_upper):    # generate an extended list of rs for recalculation
    min = WSRs[0]
    max = WSRs[0]
    for r in WSRs:
        if r < min:
            min = r
        if r > max:
            max = r
    ex_min = f_lower * min              # widen the calculation range to better show the trend
    ex_max = f_upper * max
    return np.linspace(ex_min, ex_max, 50)

def Evaluation(element, WSRs_Ang, Es_cohesive_eV, k=0.70):    # importa data
    Es_cohesive = []
    for E in Es_cohesive_eV:
        Es_cohesive.append(E * 1.602)             # e-19 J
    Strc = Structure(element)                     # specify the element name and create an obj
    Strc.import_exp_data((WSRs_Ang, Es_cohesive)) # import exp.data for processing
    E_fit_polynomial = Strc.Poly_fitting_Er()     # Step 1: fitting with polynominal
    Strc.Eval_Equilibrium_Poly(start=1.6)         # Step 2: determine r0 and E0
    E_fit_empirical = Strc.E_total_fitting_Er()   # Step 3: fitting with empirical model, determine LengthScale
    BulkModulus = Calc_BulkModulus(Strc.r0, Strc.E0, Strc.LengthScale, 'Pa')
    DebyeT = Calc_DebyeT(k, Strc.r0, BulkModulus, Strc.MolarWeight)

    plt.plot(WSRs_Ang, E_fit_empirical, '+', label = 'fitting, empirical')
    plt.plot(WSRs_Ang, Es_cohesive, 'o', label = 'Experimental data')
    plt.plot(WSRs_Ang, E_fit_polynomial, '+', label='fitting, polynomial')

    lst_r_poly = WSRs4recalculation(WSRs_Ang, 0.9, 1.2)  # automatically generate an extended range of rs
    lst_E_poly, lst_E_empr = [], []
    for r in lst_r_poly:
        lst_E_poly.append(Strc.Poly_Er(r))
    plt.plot(lst_r_poly, lst_E_poly, '-', label='recalculation, polynomial')

    lst_r_empr = WSRs4recalculation(WSRs_Ang, 0.9, 1.5)
    for r in lst_r_empr:
        lst_E_empr.append(Strc.E_total(r, ls = Strc.LengthScale))
    plt.plot(lst_r_empr, lst_E_empr, '.-', label = 'recalculation,emperical')
    plt.xlabel('Wigner-Seitz radius, 1e-10 m')
    plt.ylabel('Cohesive energy, 1e-19 J')
    plt.legend(loc='best')

    print("Bulk Modulus: {} GPa".format(BulkModulus * 1e-9))
    print("Debye Temperature: {} K".format(DebyeT))
    plt.savefig(element + '_cohesive_energy')
    plt.show()