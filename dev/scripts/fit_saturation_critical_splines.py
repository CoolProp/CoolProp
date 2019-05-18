from __future__ import print_function
import CoolProp as CoolProp
import matplotlib.pyplot as plt
import numpy as np, os, json

# Turn off the critical splines using json
jj = json.loads(CoolProp.CoolProp.get_config_as_json_string())  # Get the json values that are set already
jj['CRITICAL_SPLINES_ENABLED'] = False
CoolProp.CoolProp.set_config_as_json_string(json.dumps(jj))  # Set the values using json
# Double check that it was set properly
jj = json.loads(CoolProp.CoolProp.get_config_as_json_string()); print(jj)


class LinearFitter(object):
    def __init__(self):
        self.Nconstraints = 0
        self.A = np.zeros((2, 2))
        self.B = np.zeros((2, 1))

    def build(self):
        if (self.Nconstraints == 2):
            self.a = np.linalg.solve(self.A, self.B).squeeze();
        else:
            raise ValueError("Number of constraints[%d] is not equal to 2" % Nconstraints);

        self.c = list(np.r_[0, 0, self.a])

    def add_value_constraint(self, x, y):
        i = self.Nconstraints;
        if (i == 2):
            return false;
        self.A[i, 0] = x;
        self.A[i, 1] = 1;
        self.B[i] = y;
        self.Nconstraints += 1;

    def add_derivative_constraint(self, x, dydx):
        i = self.Nconstraints;
        if (i == 2):
            return false;
        self.A[i, 0] = 1;
        self.A[i, 1] = 0;
        self.B[i] = dydx;
        self.Nconstraints += 1;

    def evaluate(self, x):
        return np.polyval(self.a, x)


class QuadraticFitter(object):
    def __init__(self):
        self.Nconstraints = 0
        self.A = np.zeros((3, 3))
        self.B = np.zeros((3, 1))

    def build(self):
        if (self.Nconstraints == 3):
            self.a = np.linalg.solve(self.A, self.B).squeeze();
        else:
            raise ValueError("Number of constraints[%d] is not equal to 3" % Nconstraints);

        self.c = list(np.r_[0, self.a])

    def add_value_constraint(self, x, y):
        i = self.Nconstraints;
        if (i == 4):
            return false;
        self.A[i, 0] = x * x;
        self.A[i, 1] = x;
        self.A[i, 2] = 1;
        self.B[i] = y;
        self.Nconstraints += 1;

    def add_derivative_constraint(self, x, dydx):
        i = self.Nconstraints;
        if (i == 4):
            return false;
        self.A[i, 0] = 2 * x;
        self.A[i, 1] = 1;
        self.A[i, 2] = 0;
        self.B[i] = dydx;
        self.Nconstraints += 1;

    def add_second_derivative_constraint(self, x, d2ydx2):
        i = self.Nconstraints;
        if (i == 3):
            return false;
        self.A[i, 0] = 2;
        self.A[i, 1] = 0;
        self.A[i, 2] = 0;
        self.B[i] = d2ydx2;
        self.Nconstraints += 1;

    def evaluate(self, x):
        return np.polyval(self.a, x)


class CubicFitter(object):
    def __init__(self):
        self.Nconstraints = 0
        self.A = np.zeros((4, 4))
        self.B = np.zeros((4, 1))

    def build(self):
        if (self.Nconstraints == 4):
            self.a = np.linalg.solve(self.A, self.B).squeeze();
        else:
            raise ValueError("Number of constraints[%d] is not equal to 4" % Nconstraints);

        self.c = list(self.a)

    def add_value_constraint(self, x, y):
        i = self.Nconstraints;
        if (i == 4):
            return false;
        self.A[i, 0] = x * x * x;
        self.A[i, 1] = x * x;
        self.A[i, 2] = x;
        self.A[i, 3] = 1;
        self.B[i] = y;
        self.Nconstraints += 1;

    def add_derivative_constraint(self, x, dydx):
        i = self.Nconstraints;
        if (i == 4):
            return false;
        self.A[i, 0] = 3 * x * x;
        self.A[i, 1] = 2 * x;
        self.A[i, 2] = 1;
        self.A[i, 3] = 0;
        self.B[i] = dydx;
        self.Nconstraints += 1;

    def add_second_derivative_constraint(self, x, d2ydx2):
        i = self.Nconstraints;
        if (i == 4):
            return false;
        self.A[i, 0] = 6 * x;
        self.A[i, 1] = 2;
        self.A[i, 2] = 0;
        self.A[i, 3] = 0;
        self.B[i] = d2ydx2;
        self.Nconstraints += 1;

    def evaluate(self, x):
        return np.polyval(self.a, x)


class QuarticFitter(object):
    def __init__(self):
        self.Nconstraints = 0
        self.A = np.zeros((5, 5))
        self.B = np.zeros((5, 1))

    def build(self):
        if (self.Nconstraints == 5):
            self.a = np.linalg.solve(self.A, self.B);
        else:
            raise ValueError("Number of constraints[%d] is not equal to 4" % Nconstraints);

    def add_value_constraint(self, x, y):
        i = self.Nconstraints;
        if (i == 5):
            return false;
        self.A[i, 0] = x * x * x * x;
        self.A[i, 1] = x * x * x;
        self.A[i, 2] = x * x;
        self.A[i, 3] = x;
        self.A[i, 4] = 1;
        self.B[i] = y;
        self.Nconstraints += 1;

    def add_derivative_constraint(self, x, dydx):
        i = self.Nconstraints;
        if (i == 5):
            return false;
        self.A[i, 0] = 4 * x * x * x;
        self.A[i, 1] = 3 * x * x;
        self.A[i, 2] = 2 * x;
        self.A[i, 3] = 1;
        self.A[i, 4] = 0;
        self.B[i] = dydx;
        self.Nconstraints += 1;

    def add_second_derivative_constraint(self, x, d2ydx2):
        i = self.Nconstraints;
        if (i == 5):
            return false;
        self.A[i, 0] = 12 * x * x;
        self.A[i, 1] = 6 * x;
        self.A[i, 2] = x;
        self.A[i, 3] = 0;
        self.A[i, 4] = 0;
        self.B[i] = d2ydx2;
        self.Nconstraints += 1;

    def evaluate(self, x):
        return np.polyval(self.a, x)


CoolProp.CoolProp.set_debug_level(0)

if not os.path.exists('sat_spline_json.json'):

    from matplotlib.backends.backend_pdf import PdfPages
    sat_crit_spline = PdfPages('sat_crit_spline.pdf')

    V = ''
    e = {}
    e['Tc'] = []
    e['Tt'] = []
    e['dT'] = []

    JSON = {}

    for fluid in CoolProp.__fluids__:
        perfect = False
        cL = None,
        cV = None,
        if fluid in ['R410A', 'R404A', 'R507A', 'R407C', 'SES36', 'Air', 'R407F']: continue
        dT = 3
        plt.close('all')
        fig, (ax1, ax2) = plt.subplots(1, 2)
        pc = CoolProp.CoolProp.PropsSI('pcrit', fluid)
        Tc = CoolProp.CoolProp.PropsSI('Tcrit', fluid)
        Tt = CoolProp.CoolProp.PropsSI('Ttriple', fluid)
        print(fluid, Tc, Tt)
        try:

            good_T = Tc - dT
            rhomolar_crit = CoolProp.CoolProp.PropsSI('rhomolar_critical', fluid)
            ok = False  # Start assuming not ok
            for i in np.linspace(np.log10(2), 6, 5000):
                dT = 10**(-i)
                bad_T = Tc - dT

                p = CoolProp.CoolProp.PropsSI('P', 'T', Tc - dT, 'Q', 0, fluid)

                # If this run worked, good_T and good_dT are set to the working values
                good_dT = 10**(-(i - 1))
                good_T = Tc - good_dT

                # Set the flag saying this temperature was OK
                ok = True

            good_dT = 10**(-17)
            good_T = Tc - good_dT

        except ValueError as V:
            pass

        if ok:

            rhomolar_endL = CoolProp.CoolProp.PropsSI('Dmolar', 'T', good_T, 'Q', 0, fluid)
            rhomolar_endV = CoolProp.CoolProp.PropsSI('Dmolar', 'T', good_T, 'Q', 1, fluid)

            if good_dT > 1e-3:
                # Cubic or quadratic fit
                DELTAT = 1e-5

                dT_drhoL = DELTAT / (rhomolar_endL - CoolProp.CoolProp.PropsSI('Dmolar', 'T', good_T - DELTAT, 'Q', 0, fluid))
                dT_drhoV = DELTAT / (rhomolar_endV - CoolProp.CoolProp.PropsSI('Dmolar', 'T', good_T - DELTAT, 'Q', 1, fluid))

                CFL = CubicFitter()
                CFL.add_value_constraint(rhomolar_crit, Tc)
                CFL.add_derivative_constraint(rhomolar_crit, 0)
                CFL.add_value_constraint(rhomolar_endL, good_T)
                CFL.add_derivative_constraint(rhomolar_endL, dT_drhoL)
                CFL.build()

                zeros = np.roots(np.polyder(CFL.a))
                monotonic_liq = not np.any(np.logical_and(zeros > 1.000000001 * rhomolar_crit, zeros < 0.999999999 * rhomolar_endL))
                print('monotonic_liq', monotonic_liq)

                if not monotonic_liq:

                    CFL = QuadraticFitter()
                    CFL.add_value_constraint(rhomolar_crit, Tc)
                    CFL.add_derivative_constraint(rhomolar_crit, 0)
                    CFL.add_value_constraint(rhomolar_endL, good_T)
                    ##CFL.add_second_derivative_constraint(rhomolar_endL, dT_drhoL)
                    CFL.build()

                CFV = CubicFitter()
                CFV.add_value_constraint(rhomolar_crit, Tc)
                CFV.add_derivative_constraint(rhomolar_crit, 0)
                CFV.add_value_constraint(rhomolar_endV, good_T)
                CFV.add_derivative_constraint(rhomolar_endV, dT_drhoV)
                CFV.build()

                zeros = np.roots(np.polyder(CFV.a))
                monotonic_vap = not np.any(np.logical_and(zeros < 0.999999999 * rhomolar_crit, zeros > 1.000000001 * rhomolar_endV))
                print('monotonic_vap (cubic)', monotonic_vap)

                if not monotonic_vap:

                    CFV = QuadraticFitter()
                    CFV.add_value_constraint(rhomolar_crit, Tc)
                    CFV.add_derivative_constraint(rhomolar_crit, 0)
                    CFV.add_value_constraint(rhomolar_endV, good_T)
                    ##CFV.add_second_derivative_constraint(rhomolar_endV, dT_drhoV)
                    CFV.build()
                    zeros = np.roots(np.polyder(CFV.a))
                    monotonic_vap = not np.any(np.logical_and(zeros < 0.999999999 * rhomolar_crit, zeros > 1.000000001 * rhomolar_endV))
                    print('monotonic_vap', monotonic_vap)

                # Plot the cubic part on small axis
                ax1.plot(np.linspace(rhomolar_endL, rhomolar_crit), CFL.evaluate(np.linspace(rhomolar_endL, rhomolar_crit)), 'r')
                ax1.plot(np.linspace(rhomolar_endV, rhomolar_crit), CFV.evaluate(np.linspace(rhomolar_endV, rhomolar_crit)), 'r')
                # Plot the cubic part on big ais
                ax2.plot(np.linspace(rhomolar_endL, rhomolar_crit), CFL.evaluate(np.linspace(rhomolar_endL, rhomolar_crit)), 'r')
                ax2.plot(np.linspace(rhomolar_endV, rhomolar_crit), CFV.evaluate(np.linspace(rhomolar_endV, rhomolar_crit)), 'r')

                # Evaluated from EOS
                TsL = np.linspace(good_T, good_T - good_dT, 100)
                ax1.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [0] * len(TsL), fluid), TsL, 'g')
                ax1.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [1] * len(TsL), fluid), TsL, 'g')

                ax1.plot(rhomolar_endL, good_T, 'o')
                ax1.plot(rhomolar_endV, good_T, 'o')
                ax1.text(rhomolar_crit, good_T, '{T:0.2f} K'.format(T=good_T), va='top', ha='center')

                print(fluid, ':::', dT, good_T, bad_T, CoolProp.CoolProp.PropsSI('P', 'T', good_T, 'Q', 0, fluid))

                cL = CFL.c
                cV = CFV.c

            elif dT < 1e-3 and dT > 1e-6:
                # Linear fit
                LFV = LinearFitter()
                LFV.add_value_constraint(rhomolar_crit, Tc)
                LFV.add_value_constraint(rhomolar_endV, good_T)
                LFV.build()

                LFL = LinearFitter()
                LFL.add_value_constraint(rhomolar_crit, Tc)
                LFL.add_value_constraint(rhomolar_endL, good_T)
                LFL.build()

                ax1.plot(np.linspace(rhomolar_endL, rhomolar_crit), LFL.evaluate(np.linspace(rhomolar_endL, rhomolar_crit)), 'r')
                ax1.plot(np.linspace(rhomolar_endV, rhomolar_crit), LFV.evaluate(np.linspace(rhomolar_endV, rhomolar_crit)), 'r')

                e['Tc'].append(Tc)
                e['Tt'].append(Tt)
                e['dT'].append(dT)

                cL = LFL.c
                cV = LFV.c
            else:
                # Evaluated from EOS
                TsL = np.linspace(0.999 * Tc, Tc - 1e-15, 100)
                ax1.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [0] * len(TsL), fluid), TsL, 'g')
                ax1.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [1] * len(TsL), fluid), TsL, 'g')
                ax1.text(rhomolar_crit, Tc, 'PERFECT', ha='center', va='bottom')
                perfect = True

            ax1.plot(rhomolar_crit, Tc, 'o')
            ax2.plot(rhomolar_crit, Tc, 'o', ms=2)

            ax2.axhline(Tt + 1, lw=2)
            TsL = np.linspace(Tt + 1, good_T, 10000)
            ax2.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [0] * len(TsL), fluid), TsL, 'g')
            ax2.plot(CoolProp.CoolProp.PropsSI('Dmolar', 'T', TsL, 'Q', [1] * len(TsL), fluid), TsL, 'g')

            plt.title(fluid)
            plt.xlabel(r'$\rho$ [mol/m$^3$]')
            plt.ylabel(r'T [K]')
            sat_crit_spline.savefig()
            plt.close()

            if not perfect:
                JSON[fluid] = dict(cL=cL,
                                   cV=cV,
                                   _note="Coefficients for the critical cubic spline.  T = c[0]*rho^3 + c[1]*rho^2 + c[2]*rho + c[3] with rho in mol/m^3 and T in K",
                                   T_min=good_T,
                                   T_max=Tc,
                                   rhomolar_min=rhomolar_endV,
                                   rhomolar_max=rhomolar_endL
                               )
                print('\tgood_dT', good_dT)
            else:
                print('\tperfect')
        else:
            print(fluid, 'FAIL', V)

    sat_crit_spline.close()

    plt.semilogy(np.array(e['Tc']) / np.array(e['Tt']), e['dT'], 'o', mfc='none')
    plt.close()

    fp = open('sat_spline_json.json', 'w')
    json.dump(JSON, fp)
    fp.close()

else:
    # Load the generated JSON file
    jj = json.load(open('sat_spline_json.json', 'r'))

    # Iterate over the list of fluids and inject into the fluid JSON files
    for fluid, v in jj.iteritems():
        fp = open(os.path.join('..', 'fluids', fluid + '.json'), 'r')
        fluid_json = json.load(fp)
        if 'critical_region_spline' in fluid_json['EOS'][0]:
            del fluid_json['EOS'][0]['critical_region_spline']
        fluid_json['EOS'][0]['critical_region_splines'] = v
        fp.close()
        fp = open(os.path.join('..', 'fluids', fluid + '.json'), 'w')
        json.dump(fluid_json, fp)
        fp.close()
