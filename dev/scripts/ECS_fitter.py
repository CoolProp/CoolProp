from math import sqrt,exp
import CoolProp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import log

def viscosity_dilute(fluid,T,e_k,sigma):
    """
    T in [K], e_K in [K], sigma in [nm]
    viscosity returned is in [Pa-s]
    """
    Tstar = T/e_k
    molemass = CoolProp.CoolProp.PropsSI(fluid,'molemass')*1000
    
    # From Neufeld, 1972, Journal of Chemical Physics - checked coefficients
    OMEGA_2_2 = 1.16145*pow(Tstar,-0.14874)+ 0.52487*exp(-0.77320*Tstar)+2.16178*exp(-2.43787*Tstar)
    # Using the leading constant from McLinden, 2000 since the leading term from Huber 2003 gives crazy values
    eta_star = 26.692e-3*sqrt(molemass*T)/(pow(sigma,2)*OMEGA_2_2)/1e6
    return eta_star

def get_psi(fluid, ref_fluid, eta, T, rhomolar, e_k, sigma_nm):

    THIS = CoolProp.AbstractState('HEOS', fluid)
    REF = CoolProp.AbstractState('HEOS', ref_fluid)
    RP = CoolProp.AbstractState('REFPROP', fluid)

    THIS.update(CoolProp.DmolarT_INPUTS, rhomolar, T)
    RP.update(CoolProp.DmolarT_INPUTS, rhomolar, T)

    def residual_for_psi(psi, REF):

        # Calculate the conformal state
        conformal_state = THIS.conformal_state(ref_fluid, -1, -1)
        # Calculate ESRR (which are based on the CONFORMAL state values)
        f = THIS.T()/conformal_state['T'];
        h = conformal_state['rhomolar']/THIS.rhomolar(); ## Must be the ratio of MOLAR densities!!
        
        # The F factor
        F_eta = sqrt(f)*pow(h, -2.0/3.0)*sqrt(THIS.molar_mass()/REF.molar_mass());

        # Dilute viscosity of fluid of interest
        eta_dilute = viscosity_dilute(fluid, T, e_k, sigma_nm)
        
        # Required background contribution from reference fluid
        viscosity_background_required = (eta - eta_dilute)/F_eta

        REF.update(CoolProp.DmolarT_INPUTS, conformal_state['rhomolar']*psi, conformal_state['T'])
        visc_ref = REF.viscosity_contributions()
        residual = visc_ref['initial_density'] + visc_ref['residual']
        return residual - viscosity_background_required

    psi = scipy.optimize.newton(residual_for_psi, 1.0, args = (REF,))
    return psi

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

def arrayize(*args):
    return [np.array(a) for a in args]

for fluid,data,e_k, sigma_nm in zip(['R1234yf', 'R1234ze(E)'],[data_R1234yf, data_R1234zeE],[281.14, 292.11], [0.5328, 0.5017]):
    xx, yy, RHO, ETA, ETACP, ETARP = [], [], [], [], [], []
    for line in data.split('\n'):
        T, rhoL, rhoV, etaV, nuL, sigma = line.strip().split(' ')
        rhoL = float(rhoL)
        T = float(T)
        nuL = float(nuL)
        rhomolar = rhoL/CoolProp.CoolProp.PropsSI(fluid,'molemass')
        eta = nuL/1000**2*rhoL
        psi = get_psi('R1234yf', 'Propane', eta, T, rhomolar, e_k, sigma_nm)
        xx.append(T)
        yy.append(psi)
        RHO.append(rhomolar)
        ETA.append(eta)
        ETACP.append(CoolProp.CoolProp.PropsSI('V','T',T,'Q',0,fluid))
        ETARP.append(CoolProp.CoolProp.PropsSI('V','T',T,'Q',0,'REFPROP::' + CoolProp.CoolProp.get_fluid_param_string(fluid, 'REFPROP_name')))

    RHO, xx,ETA,ETACP, ETARP = arrayize(RHO, xx, ETA, ETACP, ETARP)
    rhor = RHO/CoolProp.CoolProp.PropsSI(fluid, 'rhomolar_critical')

    plt.plot(rhor, yy, 'o-')
    p = np.polyfit(rhor, yy, 2)
    print p[::-1]
    plt.plot(rhor, np.polyval(p, rhor), 'o-')
    plt.show()

    plt.title(fluid)
    plt.plot(xx, (ETACP/ETA-1)*100,'^', label = 'CoolProp')
    plt.plot(xx, (ETARP/ETA-1)*100,'o', label = 'REFPROP')
    plt.xlabel('Temperature (K)')
    plt.ylabel('$100\\times(\eta_{calc}/\eta_{exp}-1)$ (%)')
    plt.legend(loc='best')
    plt.savefig(fluid + '_deviation.pdf')
    plt.show()