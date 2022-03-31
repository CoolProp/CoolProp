from CoolProp import CoolProp as CP
from PDSim.misc.datatypes import Collector
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.odr import *
import textwrap

fluid_REF = 'Propane'
Tcrit_REF = CP.Props(fluid_REF, 'Tcrit')
omega_REF = CP.Props(fluid_REF, "accentric")
molemass_REF = CP.Props(fluid_REF, 'molemass')
rhocrit_REF = CP.Props(fluid_REF, 'rhocrit')
Zcrit_REF = CP.DerivTerms('Z', Tcrit_REF, rhocrit_REF, fluid_REF)

fluid = 'DimethylEther'
molemass = CP.Props(fluid, 'molemass')
Ttriple = CP.Props(fluid, 'Ttriple')
Tcrit = CP.Props(fluid, 'Tcrit')
omega = CP.Props(fluid, "accentric")
rhocrit = CP.Props(fluid, 'rhocrit')
pcrit = CP.Props(fluid, 'pcrit')
Zcrit = CP.DerivTerms('Z', Tcrit, rhocrit, fluid)

N = 12

RHO, TTT, RHO0, TTT0 = Collector(), Collector(), Collector(), Collector()

rhomax = CP.Props('D', 'T', Ttriple, 'Q', 0, fluid)
# Build a database of "experimental" data
for T in np.linspace(Ttriple, Tcrit + 50, 80):
    for rho in np.linspace(1e-10, rhomax, 80):
        T0, rho0 = CP.conformal_Trho(fluid, fluid_REF, T, rho)

        p = CP.Props('P', 'T', T, 'D', rho, fluid)

        ar = CP.DerivTerms("phir", T, rho, fluid)
        ar_REF = CP.DerivTerms("phir", T0, rho0, fluid_REF)
        Z = CP.DerivTerms("Z", T, rho, fluid)
        Z_REF = CP.DerivTerms("Z", T0, rho0, fluid_REF)

        #goodstate = ((T > Tcrit and p > pcrit) or (T<Tcrit and rho > CP.rhosatL_anc(fluid,T) ))
        goodstate = (T > Tcrit or rho > CP.rhosatL_anc(fluid, T) or rho < CP.rhosatV_anc(fluid, T))
        #goodstate = True

        # Want positive value, and single-phase
        if ((T0 / T) > 0.1 and T / T0 * Tcrit_REF / Tcrit < 3 and T0 / T < 1e6 and goodstate):
            if abs((ar - ar_REF) * 2 + (Z - Z_REF)**2) > 1e-5:
                print("%s %s" % (ar - ar_REF, Z - Z_REF))
            TTT << T
            RHO << rho
            TTT0 << T0
            RHO0 << rho0

tau = Tcrit / np.array(TTT.vec)
delta = np.array(RHO.vec) / rhocrit
THETA = np.array(TTT.vec) / np.array(TTT0.vec) * Tcrit_REF / Tcrit
PHI = np.array(RHO0.vec) / np.array(RHO.vec) * rhocrit / rhocrit_REF  # Ratio of MOLAR densities - here the molar masses cancel out to make phi non-dimensional

from CoolProp.Plots.Plots import Trho
Trho(fluid)
# plt.plot(RHO.vec,TTT.vec,'.')
# plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array(RHO.vec), np.array(TTT.vec), THETA)
plt.close('all')

print('rhomin = %s' % np.min(RHO.vec))

# Define the objective function


def OBJECTIVE_theta(c, x):
    tau = x[0, :]
    delta = x[1, :]

    A1 = c[0] - c[1] * np.log(tau)
    A2 = c[2] - c[3] * np.log(tau)
    A3 = c[4] - c[5] * np.log(tau)
    A4 = c[6] - c[7] * np.log(tau)**2
    DELTA = (delta - 1)**2 + (tau - 1)**2
    PSI_theta = c[8] * delta * np.exp(-c[9] * DELTA**2)
    return 1 + (omega - omega_REF) * (A1 + A2 * np.exp(-delta**2) + A3 * np.exp(-delta**c[10]) + A4 * np.exp(-delta**c[11]) + PSI_theta)

# Define the objective function


def OBJECTIVE_phi(c, x):
    tau = x[0, :]
    delta = x[1, :]

    A1 = c[0] - c[1] * np.log(tau)
    A2 = c[2] - c[3] * np.log(tau)
    A3 = c[4] - c[5] * np.log(tau)
    A4 = c[6] - c[7] * np.log(tau)**2
    DELTA = (delta - 1)**2 + (tau - 1)**2
    PSI_theta = c[8] * delta * np.exp(-c[9] * DELTA**2)
    return Zcrit_REF / Zcrit * (1 + (omega - omega_REF) * (A1 + A2 * np.exp(-delta**2) + A3 * np.exp(-delta**c[10]) + A4 * np.exp(-delta**c[11]) + PSI_theta))


print('starting fit for theta')
XXX = np.r_[np.array(tau, ndmin=2), np.array(delta, ndmin=2)]


def fit_theta():
    mod = Model(OBJECTIVE_theta)
    mydata = Data(XXX.copy(), THETA)
    beta0 = [100 for _ in range(N)]
    myodr = ODR(mydata, mod, beta0=beta0)
    myoutput = myodr.run()
    myoutput.pprint()
    print(myoutput.sum_square)
    YFIT = OBJECTIVE_theta(myoutput.beta, XXX)
    plt.plot(THETA, YFIT, 'o', mfc='none')
    plt.show()
    ERR = YFIT - THETA
    MAE = np.mean(np.abs(YFIT / THETA - 1)) * 100
    from CoolProp.Plots.Plots import Trho
    Trho(fluid)
    plt.plot(np.array(RHO.vec)[np.abs(ERR) < 5e-2], np.array(TTT.vec)[np.abs(ERR) < 5e-2], '.')
    plt.show()

    return myoutput.beta, MAE


def fit_phi():
    mod = Model(OBJECTIVE_phi)
    mydata = Data(XXX.copy(), PHI)
    beta0 = [100 for _ in range(N)]
    myodr = ODR(mydata, mod, beta0=beta0)
    myoutput = myodr.run()
    myoutput.pprint()
    print(myoutput.sum_square)
    YFIT = OBJECTIVE_theta(myoutput.beta, XXX)
    plt.plot(PHI, YFIT, 'o', mfc='none')
    plt.show()
    ERR = YFIT - PHI
    from CoolProp.Plots.Plots import Trho
    Trho(fluid)
    plt.plot(np.array(RHO.vec)[np.abs(ERR) < 5e-2], np.array(TTT.vec)[np.abs(ERR) < 5e-2], '.')
    MAE = np.mean(np.abs(YFIT / PHI - 1)) * 100
    plt.show()

    return myoutput.beta, MAE


c, theta_MAE = fit_theta()
d, phi_MAE = fit_phi()


def write_output(c, d, theta_MAE, phi_MAE):
    import time
    from datetime import date
    cdata = ', '.join(['{val:0.16g}'.format(val=v) for v in c])
    ddata = ', '.join(['{val:0.16g}'.format(val=v) for v in d])
    name = fluid
    rhomin = np.min(RHO.vec)
    timestamp = date.fromtimestamp(time.time()).strftime("%A, %d. %B %Y")
    template = textwrap.dedent(
    """
    double {name:s}Class::viscosity_Trho(double T, double rho)
    {{
        /*
        Fitting of shape factor curves to R134a data. This method is employed because solving
        for the shape factors is computationally very expensive and not very nice
        convergence behavior is experienced.  Thus we can use the ECS method,
        but with about the execution time of a conventional viscosity correlation.

        This function code was automatically generated by the fit_shape_factor.py
        script in the dev/ folder on {timestamp:s}

        Mean absolute errors of shape factor prediction:
        theta = {theta_MAE:g} %
        phi = {phi_MAE:g} %
        */

        double e_k, sigma, tau, delta, A1, A2, A3, A4, theta, Tc, Tc0, T0, rho0;
        double DELTA, PSI_theta, psi, f, h, F_eta, M, M0, delta_omega, rho0bar;
        double B1, B2, B3, B4, PSI_phi, Zc, Zc0, rhoc0, rhoc, log_tau, phi, rhobar;

        double c[] = {{{cdata:s}}};
        double d[] = {{{ddata:s}}};

        tau = reduce.T/T;
        delta = rho/reduce.rho;

        R134aClass R134a = R134aClass();
        R134a.post_load();
        delta_omega = params.accentricfactor-R134a.params.accentricfactor;

        Zc = reduce.p/(reduce.rho*R()*reduce.T);
        Zc0 = R134a.reduce.p/(R134a.reduce.rho*R134a.R()*R134a.reduce.T);
        Tc = reduce.T;
        Tc0 = R134a.reduce.T;
        rhoc = reduce.rho;
        rhoc0 = R134a.reduce.rho;
        M = params.molemass;
        M0 = R134a.params.molemass;

        rhobar = rho/M;

        if (rho > {rhomin:g})
        {{
            DELTA = pow(delta-1,2)+pow(tau-1,2);
            log_tau = log(tau);

            A1 = c[0]-c[1]*log_tau;
            A2 = c[2]-c[3]*log_tau;
            A3 = c[4]-c[5]*log_tau;
            A4 = c[6]-c[7]*pow(log_tau,2);
            PSI_theta = c[8]*delta*exp(-c[9]*pow(DELTA,2));
            theta = 1+(delta_omega)*(A1+A2*exp(-pow(delta,2))+A3*exp(-pow(delta,c[10]))+A4*exp(-pow(delta,c[11]))+PSI_theta);

            B1 = d[0]-d[1]*log_tau;
            B2 = d[2]-d[3]*log_tau;
            B3 = d[4]-d[5]*log_tau;
            B4 = d[6]-d[7]*pow(log_tau,2);
            PSI_phi = d[8]*delta*exp(-d[9]*pow(DELTA,2));
            phi = Zc0/Zc*(1+(delta_omega)*(B1+B2*exp(-pow(delta,2))+B3*exp(-pow(delta,d[10]))+B4*exp(-pow(delta,d[11]))+PSI_phi));
        }}
        else
        {{
            // Assume unity shape factors at low density
            theta = 1.0; phi = 1.0;
        }}
        T0 = T*Tc0/theta/Tc;
        h = M/M0*rhoc0/rhoc*phi;
        rho0bar = rhobar*h;
        rho0 = M0*rho0bar;

        psi = ECS_psi_viscosity(delta);
        f = T/T0;
        F_eta = sqrt(f)*pow(h,-2.0/3.0)*sqrt(M/M0);
        ECSParams(&e_k,&sigma);
        return viscosity_dilute(T,e_k,sigma) + R134a.viscosity_background(T0,rho0*psi)*F_eta;
    }}
    """
    )

    print(template.format(**locals()))


write_output(c, d, theta_MAE, phi_MAE)
