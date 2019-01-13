from math import sqrt, exp
from CoolProp.CoolProp import Props
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
from math import log
E_K = {'REFPROP-Ammonia': 386,
       'REFPROP-Argon': 143.2
       }
SIGMA = {'REFPROP-Ammonia': 0.2957,
         'REFPROP-Argon': 0.335
         }

E_K['REFPROP-Propane'] = 263.88
SIGMA['REFPROP-Propane'] = 0.49748
E_K['REFPROP-R32'] = 289.65
SIGMA['REFPROP-R32'] = 0.4098
E_K['REFPROP-R245fa'] = 329.72
SIGMA['REFPROP-R245fa'] = 0.5529


def viscosity_dilute(fluid, T, e_k, sigma):
    """
    T in [K], e_k in [K], sigma in [nm]
    viscosity returned is in [Pa-s]
    """

    Tstar = T / e_k
    molemass = Props(fluid, 'molemass')
    if fluid == 'Propane' or fluid == 'REFPROP-Propane':
        a = [0.25104574, -0.47271238, 0, 0.060836515, 0]
        theta_star = exp(a[0] * pow(log(Tstar), 0) + a[1] * pow(log(Tstar), 1) + a[3] * pow(log(Tstar), 3));
        eta_star = 0.021357 * sqrt(molemass * T) / (pow(sigma, 2) * theta_star) / 1e6;
        return eta_star

    # From Neufeld, 1972, Journal of Chemical Physics - checked coefficients
    OMEGA_2_2 = 1.16145 * pow(Tstar, -0.14874) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar)
    # Using the leading constant from McLinden, 2000 since the leading term from Huber 2003 gives crazy values
    eta_star = 26.692e-3 * sqrt(molemass * T) / (pow(sigma, 2) * OMEGA_2_2) / 1e6
    return eta_star


def viscosity_linear(fluid, T, rho, e_k, sigma):
    """
    Implements the method of Vogel 1998 (Propane) for the linear part
    """
    N_A = 6.02214129e23
    molemass = Props(fluid, 'molemass')
    Tstar = T / e_k
    b = [-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158]
    s = sum([b[i] * pow(Tstar, -0.25 * i) for i in range(7)])

    B_eta_star = s + b[7] * pow(Tstar, -2.5) + b[8] * pow(Tstar, -5.5)  # //[no units]
    B_eta = N_A * pow(sigma / 1e9, 3) * B_eta_star  # [m3/mol]
    return viscosity_dilute(fluid, T, e_k, sigma) * B_eta * rho / molemass * 1000


from PDSim.misc.datatypes import Collector
RHO = Collector()
TT = Collector()
DELTA = Collector()
TAU = Collector()
VV = Collector()
VV0 = Collector()
VV1 = Collector()
VVH = Collector()

fluid = 'REFPROP-R32'
Tc = Props(fluid, 'Tcrit')
rhoc = Props(fluid, 'rhocrit')
for T in np.linspace(290, Props(fluid, 'Tcrit') - 0.1, 100):
    rhoV = Props('D', 'T', T, 'Q', 1, fluid)
    rhoL = Props('D', 'T', T, 'Q', 0, fluid)
    rhomax = Props('D', 'T', Props(fluid, 'Tmin'), 'Q', 0, fluid)
    for rho in list(np.linspace(rhoL, rhomax, 100)):  # +list(np.linspace(rhoV,0.0001,100)):
    # for rho in list(np.linspace(rhoV,0.0001,100)):
        mu_0 = viscosity_dilute(fluid, T, E_K[fluid], SIGMA[fluid])
        mu_1 = viscosity_linear(fluid, T, rho, E_K[fluid], SIGMA[fluid])
        mu = Props('V', 'T', T, 'D', rho, fluid)
        VV << mu
        VV0 << mu_0
        VV1 << mu_1
        VVH << mu - mu_0 - mu_1
        TT << T
        RHO << rho
        DELTA << rho / rhoc
        TAU << Tc / T


def f_RHS(E, DELTA_TAU, VV):
    k = 0
    sum = 0
    DELTA = DELTA_TAU[0, :]
    TAU = DELTA_TAU[1, :]
    for i in range(2, 5):
        for j in range(3):
            sum += E[k] * DELTA**i / TAU**j
            k += 1
#    f1,f2,f3,g1,g2 = E[k],E[k+1],E[k+2],E[k+3],E[k+4]
#    DELTA0 = g1*(1+g2*np.sqrt(TAU))
#    sum += (f1+f2/TAU+f3/TAU/TAU)*(DELTA/(DELTA0-DELTA)-DELTA/DELTA0)
    print('%s %%' % np.mean(np.abs(((sum / VV - 1) * 100))))
    return sum


log_muH = np.log(VVH.v().T)

x = np.c_[DELTA.v().T, TAU.v().T].T
y = VVH.v()

linear = Model(f_RHS, extra_args=(y,))
mydata = Data(x, y)
myodr = ODR(mydata, linear, beta0=np.array([0.1] * 17),)
myoutput = myodr.run()
E = myoutput.beta
print(E)

#plt.plot(TT.vec, y,'b.',TT.vec, f_RHS(E, x, y),'r.')
# plt.show()
# plt.plot()
plt.plot(y.T, f_RHS(E, x, y))
plt.show()
