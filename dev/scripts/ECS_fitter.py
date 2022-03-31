from __future__ import print_function

from math import sqrt, exp
import CoolProp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import log


def viscosity_dilute(fluid, T, e_k, sigma):
    """
    T in [K], e_K in [K], sigma in [nm]
    viscosity returned is in [Pa-s]
    """
    Tstar = T / e_k
    molemass = CoolProp.CoolProp.PropsSI(fluid, 'molemass') * 1000

    # From Neufeld, 1972, Journal of Chemical Physics - checked coefficients
    OMEGA_2_2 = 1.16145 * pow(Tstar, -0.14874) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar)
    # Using the leading constant from McLinden, 2000 since the leading term from Huber 2003 gives crazy values
    eta_star = 26.692e-3 * sqrt(molemass * T) / (pow(sigma, 2) * OMEGA_2_2) / 1e6
    return eta_star


def get_psi(fluid, ref_fluid, eta, T, rhomolar, e_k, sigma_nm):

    THIS = CoolProp.AbstractState('HEOS', fluid)
    REF = CoolProp.AbstractState('HEOS', ref_fluid)

    THIS.update(CoolProp.DmolarT_INPUTS, rhomolar, T)

    def residual_for_psi(psi, REF):

        # Calculate the conformal state
        conformal_state = THIS.conformal_state(ref_fluid, -1, -1)
        # Calculate ESRR (which are based on the CONFORMAL state values)
        f = THIS.T() / conformal_state['T'];
        h = conformal_state['rhomolar'] / THIS.rhomolar();  # Must be the ratio of MOLAR densities!!

        # The F factor
        F_eta = sqrt(f) * pow(h, -2.0 / 3.0) * sqrt(THIS.molar_mass() / REF.molar_mass());

        # Dilute viscosity of fluid of interest
        eta_dilute = viscosity_dilute(fluid, T, e_k, sigma_nm)

        # Required background contribution from reference fluid
        viscosity_background_required = (eta - eta_dilute) / F_eta

        REF.update(CoolProp.DmolarT_INPUTS, conformal_state['rhomolar'] * psi, conformal_state['T'])
        visc_ref = REF.viscosity_contributions()
        residual = visc_ref['initial_density'] + visc_ref['residual']
        return residual - viscosity_background_required

    psi = scipy.optimize.newton(residual_for_psi, 1.0, args=(REF,))
    return psi


def arrayize(*args):
    return [np.array(a) for a in args]


def HFO():
    # Data from Zhao et al. dx.doi.org/10.1021/je5001457 | J. Chem. Eng. Data 2014, 59, 1366-1371
    data_R1234yf = """293.15 1109.9 32.8 12.04 0.1442 6.82
    303.09 1073.5 43.6 12.53 0.1319 5.71
    313.20 1033.6 57.8 13.16 0.1223 4.60
    323.19 990.2 76.0 13.88 0.1126 3.55
    333.14 941.4 99.7 14.82 0.1016 2.55
    343.11 883.5 132.2 16.12 0.0899 1.64
    353.08 809.6 179.9 18.17 0.0820 0.81
    358.05 761.5 214.8 19.78 0.0770 0.46
    363.05 695.7 267.7 22.44 0.0700 0.15
    365.05 657.4 301.0 24.26 0.0624 0.05"""

    # Data from Zhao et al. dx.doi.org/10.1021/je5001457 | J. Chem. Eng. Data 2014, 59, 1366-1371
    data_R1234zeE = """295.23 1172.5 24.1 12.11 0.1776 8.88
    303.19 1146.1 30.6 12.46 0.1607 7.91
    313.21 1111.1 40.8 12.93 0.1429 6.66
    323.19 1073.6 53.6 13.46 0.1319 5.48
    333.00 1033.3 69.8 14.06 0.1193 4.36
    343.05 986.7 91.3 14.82 0.1132 3.30
    353.00 924.0 119.7 15.80 0.1051 2.26
    363.12 866.8 160.4 17.28 0.0924 1.35
    373.14 776.9 225.2 19.89 0.0817 0.54"""

    for fluid, data, e_k, sigma_nm in zip(['R1234yf', 'R1234ze(E)'], [data_R1234yf, data_R1234zeE], [281.14, 292.11], [0.5328, 0.5017]):
        xx, yy, RHO, ETA, ETACP, ETARP = [], [], [], [], [], []
        for line in data.split('\n'):
            T, rhoL, rhoV, etaV, nuL, sigma = line.strip().split(' ')
            rhoL = float(rhoL)
            T = float(T)
            nuL = float(nuL)
            rhomolar = rhoL / CoolProp.CoolProp.PropsSI(fluid, 'molemass')
            eta = nuL / 1000**2 * rhoL
            psi = get_psi(fluid, 'Propane', eta, T, rhomolar, e_k, sigma_nm)
            xx.append(T)
            yy.append(psi)
            RHO.append(rhomolar)
            ETA.append(eta)
            ETACP.append(CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 0, fluid))
            ETARP.append(CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 0, 'REFPROP::' + CoolProp.CoolProp.get_fluid_param_string(fluid, 'REFPROP_name')))

        RHO, xx, ETA, ETACP, ETARP = arrayize(RHO, xx, ETA, ETACP, ETARP)
        rhor = RHO / CoolProp.CoolProp.PropsSI(fluid, 'rhomolar_critical')

        plt.plot(rhor, yy, 'o-', label='from experimental data')
        p = np.polyfit(rhor, yy, 2)
        print(p[::-1])
        plt.plot(rhor, np.polyval(p, rhor), 'o-', label='from correlation')
        plt.xlabel(r'$\rho_r$')
        plt.ylabel('$\psi$')
        plt.legend(loc='best')
        plt.show()

        plt.title(fluid)
        plt.plot(xx, (ETACP / ETA - 1) * 100, '^', label='CoolProp')
        plt.plot(xx, (ETARP / ETA - 1) * 100, 'o', label='REFPROP')
        plt.xlabel('Temperature (K)')
        plt.ylabel('$100\\times(\eta_{calc}/\eta_{exp}-1)$ (%)')
        plt.legend(loc='best')
        plt.savefig(fluid + '_deviation.pdf')
        plt.show()


def pentanes():
    # from doi 10.1021/je0202174 | J. Chem. Eng. Data 2003, 48, 1418-1421
    # T (K), rhoL (kg/m^3), rhoV (kg/m^3), eta (mPa-s)
    data_cyclopentane = """253.15 258.15 263.15 268.15 273.15 278.15 283.15 288.15 293.15 298.15 303.15 308.15 313.15 318.15 323.15 328.15 333.15 338.15 343.15 348.15 353.15
    784.64 779.53 774.59 769.77 765.12 760.20 755.32 750.27 745.02 738.63 731.97 725.15 718.32 711.59 705.11 699.08 693.40 688.44 684.25 680.96 678.71
    0.0881 0.1127 0.1443 0.1848 0.2368 0.3036 0.3894 0.4999 0.6421 0.8255 1.062 1.368 1.764 2.279 2.950 3.827 4.980 6.509 8.554 11.33 15.20
    0.7268 0.6786 0.6347 0.5930 0.5567 0.5224 0.4922 0.4646 0.4382 0.4148 0.3923 0.3714 0.3521 0.3350 0.3190 0.3048 0.2912 0.2793 0.2690 0.2590 0.2502"""

    # from doi 10.1021/je0202174 | J. Chem. Eng. Data 2003, 48, 1418-1421
    # T (K), rhoL (kg/m^3), rhoV (kg/m^3), eta (mPa-s)
    data_isopentane = """253.15 258.15 263.15 268.15 273.15 278.15 283.15 288.15 293.15 298.15 303.15 308.15 313.15 318.15 323.15 328.15 333.15 338.15 343.15 348.15 353.15
    658.32 653.55 648.73 643.87 639.01 634.15 629.35 624.63 620.05 615.69 610.87 605.63 600.05 594.23 588.24 582.18 576.13 570.18 564.41 558.92 553.79
    0.4655 0.5889 0.7372 0.9137 1.122 1.366 1.650 1.979 2.356 2.788 3.278 3.833 4.459 5.162 5.949 6.827 7.803 8.886 10.09 11.41 12.87
    0.3893 0.3661 0.3439 0.3201 0.3023 0.2859 0.2703 0.2547 0.2399 0.2289 0.2144 0.2023 0.1910 0.1813 0.1724 0.1611 0.1543 0.1480 0.1411 0.1332 0.1287"""

    fluid = ''

    def undelimit(args, delim=''):
        return [np.array([float(_) for _ in a.strip().split(delim)]) for a in args]

    from CoolProp.CoolProp import PropsSI
    for fluid, e_k, sigma_nm in zip(['CycloPentane', 'Isopentane'], [406.33, 341.06], [0.518, 0.56232]):
        xx, yy, RHO, ETA, ETACP, ETARP = [], [], [], [], [], []
        for _T, _rhoLmass, _rhoVmass, _eta_mPas in zip(*undelimit(data_cyclopentane.split('\n'), delim=' ')):
            MM = PropsSI('molemass', fluid)
            rhomolar = _rhoLmass / MM
            eta = _eta_mPas / 1000
            psi = get_psi(fluid, 'Propane', eta, _T, rhomolar, e_k, sigma_nm)
            xx.append(_T)
            yy.append(psi)
            RHO.append(rhomolar)
            try:
                ETACP.append(CoolProp.CoolProp.PropsSI('V', 'T', _T, 'Q', 0, fluid))
            except:
                ETACP.append(np.nan)
            ETA.append(eta)
            ETARP.append(CoolProp.CoolProp.PropsSI('V', 'T', _T, 'Q', 0, 'REFPROP::' + CoolProp.CoolProp.get_fluid_param_string(fluid, 'REFPROP_name')))
        xx, yy, ETACP, ETARP, RHO, = arrayize(xx, yy, ETACP, ETARP, RHO)
        rhored = CoolProp.CoolProp.PropsSI(fluid, 'rhomolar_critical')
        print('rhored', rhored)
        rhor = np.array(RHO) / rhored

        plt.title(fluid)
        plt.plot(rhor, yy, 'o-', label='from experimental data')
        p = np.polyfit(rhor, yy, 2)
        print(p[::-1])
        plt.plot(rhor, np.polyval(p, rhor), 'o-', label='from correlation')
        plt.xlabel(r'$\rho_r$')
        plt.ylabel('$\psi$')
        plt.legend(loc='best')
        plt.show()

        plt.title(fluid)
        plt.plot(xx, (ETACP / ETA - 1) * 100, '^', label='CoolProp')
        plt.plot(xx, (ETARP / ETA - 1) * 100, 'o', label='REFPROP')
        plt.xlabel('Temperature (K)')
        plt.ylabel('$100\\times(\eta_{calc}/\eta_{exp}-1)$ (%)')
        plt.legend(loc='best')
        plt.savefig(fluid + '_deviation.pdf')
        plt.show()


HFO()
pentanes()
