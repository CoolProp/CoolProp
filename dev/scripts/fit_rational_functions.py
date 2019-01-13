from __future__ import division, print_function

import json
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import CoolProp
import numpy as np
import scipy.optimize
import xalglib
import os, sys


def fit_rational_polynomial(x, y, xfine, n, d):

    def obj(x, *AB):
        """
        The objective function to be minimized
        """
        A = AB[0:n + 1]
        B = list(AB[n + 1::]) + [0]
        yfit = np.polyval(A, x) / (1 + np.polyval(B, x))
        return yfit

    if d != 1:
        # n+d+1 coefficients to be solved for, select to distribute randomly
        indices = np.array(np.linspace(0, len(x) - 1, n + d + 1), dtype=int)
        xlin = x[indices]
        ylin = y[indices]

        # Solve the linear problem where Lx=R where R is the vector y, x are the unknowns A and B joined together, and L is the set of column vectors
        # L = [x^n, ... , x^1, 1, x^n,  -x^d*y,  ..., -x*y]
        R = ylin[:]
        L = np.ones((n + d + 1, n + d + 1))
        for i in range(0, n + 1):
            L[:, i] = xlin**(n - i)
        for j in range(0, d):
            L[:, n + 1 + j] = -xlin**(d - j) * ylin

        ABlin = np.linalg.solve(L, R)
        A = ABlin[0:n + 1]
        B = list(ABlin[n + 1::]) + [0]
        yfitlin = np.polyval(A, x) / (1 + np.polyval(B, x))

        AB = scipy.optimize.curve_fit(obj, x, y, p0=ABlin)[0]

        poles = np.roots(list(AB[n + 1::]) + [1])
        poles = poles[np.isreal(poles)]
        poles = poles[poles > min(x)]
        poles = poles[poles < max(x)]

    else:
        # Find a pole that is not in the range of x
        def obj2(Tpole, x, y, AB):
            B = -1 / Tpole
            A = np.polyfit(x, y * (x * B + 1), n)
            yfit = np.polyval(A, x) / (x * B + 1)
            AB[:] = list(A) + [B]  # Set so that it uses the AB passed in rather than making local variable
            rms = np.sqrt(np.sum(np.power(yfit - y, 2)))
            return rms

        AB = []
        scipy.optimize.fminbound(obj2, Tc + 0.1, 1.5 * Tc, args=(x, y, AB))

    return dict(max_abs_error=np.max(np.abs(obj(x, *AB) - y)),
                yfitnonlin=obj(xfine, *AB),
                A=AB[0:n + 1],
                B=list(AB[n + 1::]) + [1]
                )


class SplineFitter(object):
    def __init__(self):
        self.Nconstraints = 0
        self.A = np.zeros((4, 4))
        self.B = np.zeros((4, 1))

    def build(self):
        if (self.Nconstraints == 4):
            abcd = np.linalg.solve(self.A, self.B);
            self.a = abcd[0];
            self.b = abcd[1];
            self.c = abcd[2];
            self.d = abcd[3];
        else:
            raise ValueError("Number of constraints[%d] is not equal to 4" % Nconstraints);

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

    def evaluate(self, x):
        return self.a * x**3 + self.b * x**2 + self.c * x + self.d;

# See http://stackoverflow.com/a/4983359


def strictly_increasing(L):
    return all(x < y for x, y in zip(L, L[1:]))


def strictly_decreasing(L):
    return all(x > y for x, y in zip(L, L[1:]))


if not os.path.exists('hsancillaries.json'):
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('multipage.pdf')

    jj = {}
    for i, fluid in enumerate(sorted(CoolProp.__fluids__)):

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

        plt.suptitle(fluid)
        print(i, fluid)

        N = 10000
        Tc = CP.PropsSI(fluid, 'Tcrit')
        rhoc = CP.PropsSI(fluid, 'rhocrit')
        try:
            T = np.r_[np.linspace(Tc - 0.1, CP.PropsSI(fluid, 'Tmin'), N)]  # , np.linspace(Tc-0.1, Tc-1e-8, N)]
            hfg = (np.array(CP.PropsSI('Hmolar', 'T', T, 'Q', 1, fluid)) - np.array(CP.PropsSI('Hmolar', 'T', T, 'Q', 0, fluid)))
            sfg = (np.array(CP.PropsSI('Smolar', 'T', T, 'Q', 1, fluid)) - np.array(CP.PropsSI('Smolar', 'T', T, 'Q', 0, fluid)))
            if np.any(np.isnan(hfg)):
                print('nan values in hfg')
                good_values = np.isfinite(hfg)
                T = T[good_values]
                hfg = hfg[good_values]
                sfg = sfg[good_values]

        except BaseException as E:
            print(E)
            continue

        try:
            p = CP.PropsSI('P', 'T', T, 'Q', 0, fluid)
            rho = CP.PropsSI('D', 'T', T, 'Q', 0, fluid)
            Tanchor = 1.1 * Tc
            rhoanchor = 0.9 * rhoc
            hanchor_molar = CP.PropsSI('Hmolar', 'T', Tanchor, 'D', rhoanchor, fluid)
            sanchor_molar = CP.PropsSI('Smolar', 'T', Tanchor, 'D', rhoanchor, fluid)
            sL = np.array(CP.PropsSI('Smolar', 'T', T, 'Q', 0, fluid)) - sanchor_molar
            hL = np.array(CP.PropsSI('Hmolar', 'T', T, 'Q', 0, fluid)) - hanchor_molar

            x = T
            xfine = np.linspace(np.min(x), np.max(x), 5000)

            n = 7
            d = 1

            commons = dict(type="rational_polynomial", Tmin=np.min(T), Tmax=np.max(T))

            rp = fit_rational_polynomial(x, hL, xfine, n, d)
            ax1.plot(x, hL)
            ax1.plot(xfine, rp['yfitnonlin'], 'r')
            ax1.plot(xfine, rp['yfitnonlin'] + rp['max_abs_error'], 'k--')
            ax1.plot(xfine, rp['yfitnonlin'] - rp['max_abs_error'], 'k--')
            hLdict = dict(A=rp['A'][::-1], B=rp['B'][::-1], max_abs_error=rp['max_abs_error'], _note="coefficients are in increasing order; input in K, output in J/mol; value is enthalpy minus hs_anchor enthalpy", max_abs_error_units='J/mol', **commons)
            if (np.any(np.isnan(hLdict['A']))):
                print('bad A for hL')
                continue

            rp = fit_rational_polynomial(x, hfg, xfine, n, d)
            ax2.plot(x, hfg)
            ax2.plot(xfine, rp['yfitnonlin'], 'r')
            ax2.plot(xfine, rp['yfitnonlin'] + rp['max_abs_error'], 'k--')
            ax2.plot(xfine, rp['yfitnonlin'] - rp['max_abs_error'], 'k--')
            hLVdict = dict(A=rp['A'][::-1], B=rp['B'][::-1], max_abs_error=rp['max_abs_error'], _note="coefficients are in increasing order; input in K, output in J/mol; value is enthalpy minus hs_anchor enthalpy", max_abs_error_units='J/mol', **commons)

            rp = fit_rational_polynomial(x, sL, xfine, n, d)
            ax3.plot(x, sL)
            ax3.plot(xfine, rp['yfitnonlin'], 'r')
            ax3.plot(xfine, rp['yfitnonlin'] + rp['max_abs_error'], 'k--')
            ax3.plot(xfine, rp['yfitnonlin'] - rp['max_abs_error'], 'k--')
            sLdict = dict(A=rp['A'][::-1], B=rp['B'][::-1], max_abs_error=rp['max_abs_error'], _note="coefficients are in increasing order; input in K, output in J/mol/K; value is entropy minus hs_anchor entropy", max_abs_error_units='J/mol/K', **commons)
            if (np.any(np.isnan(sLdict['A']))):
                print('bad A for sL')
                continue

            rp = fit_rational_polynomial(x, sfg, xfine, n, d)
            ax4.plot(x, sfg)
            ax4.plot(xfine, rp['yfitnonlin'], 'r')
            ax4.plot(xfine, rp['yfitnonlin'] + rp['max_abs_error'], 'k--')
            ax4.plot(xfine, rp['yfitnonlin'] - rp['max_abs_error'], 'k--')
            sLVdict = dict(A=rp['A'][::-1], B=rp['B'][::-1], max_abs_error=rp['max_abs_error'], _note="coefficients are in increasing order; input in K, output in J/mol/K; value is entropy minus hs_anchor entropy", max_abs_error_units='J/mol/K', **commons)

            jj[fluid] = dict(hL=hLdict, hLV=hLVdict, sL=sLdict, sLV=sLVdict)

        except BaseException as E:
            continue
            print(E)

        pp.savefig()
        plt.close()

    pp.close()

    fp = open('hsancillaries.json', 'w')
    fp.write(json.dumps(jj, **{'indent': 2, 'sort_keys': True}))
    fp.close()
else:
    # Inject
    fp = open('hsancillaries.json', 'r')
    ancillaries = json.load(fp)
    fp.close()

    for fluid in ancillaries:

        fluid_path = '../fluids/' + fluid + '.json'

        # Open the fluid JSON file
        fp = open(fluid_path, 'r')
        j = json.load(fp)
        fp.close()

        j['ANCILLARIES']['sL'] = ancillaries[fluid]['sL']
        j['ANCILLARIES']['hL'] = ancillaries[fluid]['hL']
        j['ANCILLARIES']['sLV'] = ancillaries[fluid]['sLV']
        j['ANCILLARIES']['hLV'] = ancillaries[fluid]['hLV']

        sys.path.append('..')
        from package_json import json_options
        fp = open(fluid_path, 'w')
        fp.write(json.dumps(j, **json_options))
        fp.close()

        print('writing ' + fluid)
