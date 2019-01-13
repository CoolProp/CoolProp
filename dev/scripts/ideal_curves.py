import numpy as np
import matplotlib.pyplot as plt
import CoolProp, scipy.optimize


class CurveTracer(object):

    def __init__(self, backend, fluid, p0, T0):
        """
        p0 : Initial pressure [Pa]
        T0 : Initial temperatrure [K]
        """
        self.P = [p0]
        self.T = []
        self.AS = CoolProp.AbstractState(backend, fluid)

        # Solve for Temperature for first point
        T = scipy.optimize.newton(self.objective_T, T0, args=(p0, -1))

        self.T.append(T)

    def objective_T(self, T, p, rho_guess):
        """ Base class function """
        if rho_guess < 0:
            self.AS.update(CoolProp.PT_INPUTS, p, T)
        else:
            guesses = CoolProp.CoolProp.PyGuessesStructure()
            guesses.rhomolar = rho_guess
            self.AS.update_with_guesses(CoolProp.PT_INPUTS, p, T, guesses)
        return self.objective()

    def TPcoords(self, t, lnT, lnp, rlnT=0.1, rlnp=0.1):
        return np.exp(lnT + rlnT * np.cos(t)), np.exp(lnp + rlnp * np.sin(t))

    def obj_circle(self, t, lnT, lnp):
        T2, P2 = self.TPcoords(t, lnT, lnp)
        self.AS.update(CoolProp.PT_INPUTS, P2, T2)
        r = self.objective()
        return r

    def trace(self):
        t = self.starting_direction()
        for i in range(1000):
            try:
                lnT = np.log(self.T[-1])
                lnp = np.log(self.P[-1])
                t = scipy.optimize.brentq(self.obj_circle, t - np.pi / 2, t + np.pi / 2, args=(lnT, lnp))
                T2, P2 = self.TPcoords(t, lnT, lnp)
                self.T.append(T2)
                self.P.append(P2)
                if self.T[-1] < self.AS.keyed_output(CoolProp.iT_triple) or self.P[-1] > 1000 * self.AS.keyed_output(CoolProp.iP_critical):
                    break
            except ValueError as VE:
                print(VE)
                break

        return self.T, self.P


class IdealCurveTracer(CurveTracer):
    def __init__(self, *args, **kwargs):
        CurveTracer.__init__(self, *args, **kwargs)

    def objective(self):
        """ Z = 1 """
        return self.AS.keyed_output(CoolProp.iZ) - 1

    def starting_direction(self):
        """ Start searching directly up ( or calculate as orthogonal to gradient ) """
        return np.pi / 2.0


class BoyleCurveTracer(CurveTracer):
    def __init__(self, *args, **kwargs):
        CurveTracer.__init__(self, *args, **kwargs)

    def objective(self):
        """ dZ/dv|T = 0 """
        r = (self.AS.p() - self.AS.rhomolar() * self.AS.first_partial_deriv(CoolProp.iP, CoolProp.iDmolar, CoolProp.iT)) / (self.AS.gas_constant() * self.AS.T())
        # print self.AS.T(), self.AS.p(), r
        return r

    def starting_direction(self):
        """ Start searching directly up """
        return np.pi / 2.0


class JouleInversionCurveTracer(CurveTracer):
    def __init__(self, *args, **kwargs):
        CurveTracer.__init__(self, *args, **kwargs)

    def objective(self):
        """ dZ/dT|v = 0 """
        r = (self.AS.gas_constant() * self.AS.T() * 1 / self.AS.rhomolar() * self.AS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmolar) - self.AS.p() * self.AS.gas_constant() / self.AS.rhomolar()) / (self.AS.gas_constant() * self.AS.T())**2
        # print self.AS.T(), self.AS.p(), r
        return r

    def starting_direction(self):
        """ Start searching directly up """
        return np.pi / 2.0


class JouleThomsonCurveTracer(CurveTracer):
    def __init__(self, *args, **kwargs):
        CurveTracer.__init__(self, *args, **kwargs)

    def objective(self):
        """ dZ/dT|p = 0 """
        dvdT__constp = -self.AS.first_partial_deriv(CoolProp.iDmolar, CoolProp.iT, CoolProp.iP) / self.AS.rhomolar()**2
        r = self.AS.p() / (self.AS.gas_constant() * self.AS.T()**2) * (self.AS.T() * dvdT__constp - 1 / self.AS.rhomolar())
        # print self.AS.T(), self.AS.p(), r
        return r

    def starting_direction(self):
        """ Start searching directly up """
        return np.pi / 2.0


backend = 'HEOS'
fluid = 'R125'

kwargs = dict(lw=2)
print('Ideal')
ICT = IdealCurveTracer(backend, fluid, p0=1e5, T0=900)
T, p = ICT.trace()
plt.plot(T, p, '-', label='Ideal Curve', **kwargs)

print('Boyle')
BCT = BoyleCurveTracer(backend, fluid, p0=1e5, T0=800)
T, p = BCT.trace()
plt.plot(T, p, '-', label='Boyle Curve', **kwargs)

print('Joule Inversion')
JIT = JouleInversionCurveTracer(backend, fluid, p0=1e5, T0=1800)
T, p = JIT.trace()
plt.plot(T, p, '-', label='Joule Inversion Curve', **kwargs)

print('Joule-Thomson')
JTCT = JouleThomsonCurveTracer(backend, fluid, p0=1e5, T0=1800)
T, p = JTCT.trace()
plt.plot(T, p, '-', label='Joule-Thomson Curve', **kwargs)

print('Saturation Curve')
Tt = ICT.AS.keyed_output(CoolProp.iT_triple)
Tc = ICT.AS.keyed_output(CoolProp.iT_critical)
Ts = np.linspace(Tt, Tc - 1.e-6)
ps = CoolProp.CoolProp.PropsSI('P', 'T', Ts, 'Q', 0, backend + '::' + fluid)
plt.plot(Ts, ps, '-', label='Saturation Curve', **kwargs)

plt.yscale('log')
plt.xscale('log')
plt.xlabel('T (K)')
plt.ylabel('p (Pa)')
plt.legend(loc='best')
plt.savefig('IdealCurves.png')
plt.show()
