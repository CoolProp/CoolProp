.. _Fluid-Properties:

Pure and Pseudo-Pure fluid properties
=====================================

.. contents:: :depth: 2

Introduction
------------

Nearly all the fluids modeling in CoolProp are based on Helmholtz energy formulations.  This is a convenient construction of the equation of state because all the thermodynamic properties of interest can be obtained directly from partial derivatives of the Helmholtz energy.

It should be noted that the EOS are typically valid over the entire range of the fluid, from subcooled liquid to superheated vapor, to supercritical fluid.  

Annoyingly, different authors have selected different sets of nomenclature for the Helmholtz energy.  For consistency, the nomenclature of Lemmon will be used here.  Also, some authors present results on a mole-basis or mass-basis, further complicating comparisons.

Thermodynamic properties of Fluid
---------------------------------
In general, the EOS are based on non-dimensional terms :math:`\delta` and :math:`\tau`, where these terms are defined by

.. math::

    \delta=\rho/\rho_c
    
    \tau=T_c/T
    
where :math:`\rho_c` and :math:`T_c` are the critical density of the fluid if it is a pure fluid.  For pseudo-pure mixtures, the critical point is typically not used as the reducing state point, and often the maximum condensing temperature on the saturation curve is used instead.

The non-dimensional Helmholtz energy of the fluid is given by

.. math::

    \alpha=\alpha^0+\alpha^r
    
where :math:`\alpha^0` is the ideal-gas contribution to the Helmholtz energy, and :math:`\alpha^r` is the residual Helmholtz energy contribution which accounts for non-ideal behavior.  For a given set of :math:`\delta` and :math:`\tau`, each of the terms :math:`\alpha^0` and :math:`\alpha^r` are known.  The exact form of the Helmholtz energy terms is fluid dependent, but a relatively simple example is that of Nitrogen, which has the ideal-gas Helmholtz energy of

.. math::

    \alpha^0=\ln\delta+a_1\ln\tau+a_2+a_3\tau+a_4\tau^{-1}+a_5\tau^{-2}+a_6\tau^{-3}+a_7\ln[1-\exp(-a_8\tau)]
    
and the non-dimensional residual Helmholtz energy of

.. math::

    \alpha^r=\sum_{k=1}^{6}{N_k\delta^{i_k}\tau^{j_k}}+\sum_{k=7}^{32}{N_k\delta^{i_k}\tau^{j_k}\exp(-\delta^{l_k})}+\sum_{k=33}^{36}{N_k\delta^{i_k}\tau^{j_k}\exp(-\phi_k(\delta-1)^2-\beta_k(\tau-\gamma_k)^2)}
    
and all the terms other than :math:`\delta` and :math:`\tau` are fluid-dependent correlation parameters.

The other thermodynamic parameters can then be obtained through analytic derivatives of the Helmholtz energy terms.  For instance, the pressure is given by

.. math::

    p=\rho RT\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau} \right]
    
and the specific internal energy by

.. math::

    \frac{u}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]

and the specific enthalpy by

.. math::

    \frac{h}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] +\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+1

which can also be written as

.. math::

    \frac{h}{RT}=\frac{u}{RT}+\frac{p}{\rho RT}
    
The specific entropy is given by

.. math::

    \frac{s}{R}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]-\alpha^0-\alpha^r
    
and the specific heats at constant volume and constant pressure respectively are given by

.. math::

    \frac{c_v}{R}=-\tau^2 \left[\left(\frac{\partial^2 \alpha^0}{\partial \tau^2}\right)_{\delta}+ \left(\frac{\partial^2 \alpha^r}{\partial \tau^2}\right)_{\delta} \right]
    
    \frac{c_p}{R}=\frac{c_v}{R}+\dfrac{\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}-\delta\tau\left(\frac{\partial^2 \alpha^r}{\partial \delta\partial\tau}\right)\right]^2}{\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]}
    
The EOS is set up with temperature and density as the two independent properties, but often other inputs are known, most often temperature and pressure because they can be directly measured.  As a result, if the density is desired for a known temperature and pressure, it can be obtained iteratively.  The following algorithm is used to obtain a reasonable guess for the initial value for the iterative solver:

#. If the fluid is superheated, use a guess of ideal gas (:math:`\rho=p/(RT)`)
#. If the fluid is subcooled, use a guess of saturated liquid density
#. If the fluid is supercritical, use a guess of ideal gas (:math:`\rho=p/(RT)`)
#. No solution for density as a function of temperature and pressure if the fluid is two-phase
    
The documentation of the :mod:`CoolProp.CoolProp` module, or the :mod:`CoolProp.State` module are also available.

.. _list_of_fluids:

List of Fluids
--------------

.. note::
   You can click on the fluid name to get more information about the fluid, or click on a bracketed reference to be taken to the reference for the fluid

.. csv-table:: All the fluids included in CoolProp
   :header-rows: 1
   :widths: 40, 15, 15, 15, 15, 15, 15
   :file: PurePseudoPure.csv
   :delim: ;
   
.. include:: fluidstoc.rst.in

Ideal Curves
------------

The so-called ideal curves demonstrate the good (or not!) extrapolative behavior of the equation of state.  A few common curves are defined.  They should ideally look similar to the figure shown below.  Lemmon :cite:`Lemmon-JPCRD-2005` provides a solid coverage of the ideal curves.

Ideal Curve
^^^^^^^^^^^

.. math::

    Z = 1
    
where :math:`Z` can be given by

.. math::

    Z = \frac{pv}{RT}
    
Boyle Curve
^^^^^^^^^^^

Defined by

.. math::

    \left.\frac{\partial Z}{\partial v}\right|_{T} = 0
    
which can be expanded as (p can be expressed as a function of :math:`T` and :math:`v`)

.. math::

    \left.\frac{\partial Z}{\partial v}\right|_{T} = \dfrac{p - \rho\left.\frac{\partial p}{\partial \rho}\right|_{T}}{RT}
    
Joule Inversion Curve
^^^^^^^^^^^^^^^^^^^^^

Defined by

.. math::

    \left.\frac{\partial Z}{\partial T}\right|_{v} = 0
    
which can be expanded as (p can be expressed as a function of :math:`T` and :math:`v`)

.. math::

    \left.\frac{\partial Z}{\partial T}\right|_{v} = \frac{v}{R}\dfrac{T\left.\frac{\partial p}{\partial T}\right|_{v} - p}{T^2}
    
Joule-Thomson Curve
^^^^^^^^^^^^^^^^^^^

Defined by 

.. math::

    \left.\frac{\partial Z}{\partial T}\right|_{p} = 0
    
which can be expanded as (v can be expressed as a function of :math:`T` and :math:`p`)

.. math::

    \left.\frac{\partial Z}{\partial T}\right|_{p} = \frac{p}{R}\dfrac{T\left.\frac{\partial v}{\partial T}\right|_{p} - v}{T^2}

where

.. math::

    \left.\frac{\partial v}{\partial T}\right|_{p} = -\frac{\left.\frac{\partial \rho}{\partial T}\right|_{p}}{\rho^2}
    
    
Ideal Curves for Refrigerant R125
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import CoolProp, scipy.optimize

    class CurveTracer(object):
        
        def __init__(self, backend, fluid, p0, T0):
            """
            p0 : Initial pressure [Pa]
            T0 : Initial temperature [K]
            """
            self.P = [p0]
            self.T = []
            self.AS = CoolProp.AbstractState(backend, fluid)
            
            # Solve for Temperature for first point
            T = scipy.optimize.newton(self.objective_T, T0, args = (p0, -1))
            
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
            
        def TPcoords(self, t, lnT, lnp, rlnT = 0.1, rlnp = 0.1):
            return np.exp(lnT + rlnT*np.cos(t)), np.exp(lnp + rlnp*np.sin(t))
            
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
                    t = scipy.optimize.brentq(self.obj_circle, t-np.pi/2, t+np.pi/2, args = (lnT, lnp))
                    T2, P2 = self.TPcoords(t, lnT, lnp)
                    self.T.append(T2)
                    self.P.append(P2)
                    if self.T[-1] < self.AS.keyed_output(CoolProp.iT_triple) or self.P[-1] > 1000*self.AS.keyed_output(CoolProp.iP_critical):
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
            return np.pi/2.0
            
    class BoyleCurveTracer(CurveTracer):
        def __init__(self, *args, **kwargs):
            CurveTracer.__init__(self, *args, **kwargs)

        def objective(self):
            """ dZ/dv|T = 0 """
            r = (self.AS.p() - self.AS.rhomolar()*self.AS.first_partial_deriv(CoolProp.iP, CoolProp.iDmolar, CoolProp.iT))/(self.AS.gas_constant()*self.AS.T())
            #print self.AS.T(), self.AS.p(), r
            return r
            
        def starting_direction(self):
            """ Start searching directly up """
            return np.pi/2.0
            
    class JouleInversionCurveTracer(CurveTracer):
        def __init__(self, *args, **kwargs):
            CurveTracer.__init__(self, *args, **kwargs)

        def objective(self):
            """ dZ/dT|v = 0 """
            r = (self.AS.gas_constant()*self.AS.T()*1/self.AS.rhomolar()*self.AS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmolar)-self.AS.p()*self.AS.gas_constant()/self.AS.rhomolar())/(self.AS.gas_constant()*self.AS.T())**2
            #print self.AS.T(), self.AS.p(), r
            return r
            
        def starting_direction(self):
            """ Start searching directly up """
            return np.pi/2.0
            
    class JouleThomsonCurveTracer(CurveTracer):
        def __init__(self, *args, **kwargs):
            CurveTracer.__init__(self, *args, **kwargs)

        def objective(self):
            """ dZ/dT|p = 0 """
            dvdT__constp = -self.AS.first_partial_deriv(CoolProp.iDmolar, CoolProp.iT, CoolProp.iP)/self.AS.rhomolar()**2
            r = self.AS.p()/(self.AS.gas_constant()*self.AS.T()**2)*(self.AS.T()*dvdT__constp - 1/self.AS.rhomolar())
            #print self.AS.T(), self.AS.p(), r
            return r
            
        def starting_direction(self):
            """ Start searching directly up """
            return np.pi/2.0

    backend = 'HEOS'
    fluid = 'R125'

    kwargs = dict(lw = 2)
    print('Ideal')
    ICT = IdealCurveTracer(backend, fluid, p0 = 1e5, T0 = 900)
    T, p = ICT.trace()
    plt.plot(T, p, '-', label = 'Ideal Curve', **kwargs)

    print('Boyle')
    BCT = BoyleCurveTracer(backend, fluid, p0 = 1e5, T0 = 800)
    T, p = BCT.trace()
    plt.plot(T, p, '-', label = 'Boyle Curve', **kwargs)

    print('Joule Inversion')
    JIT = JouleInversionCurveTracer(backend, fluid, p0 = 1e5, T0 = 1800)
    T, p = JIT.trace()
    plt.plot(T, p, '-', label = 'Joule Inversion Curve', **kwargs)

    print('Joule-Thomson')
    JTCT = JouleThomsonCurveTracer(backend, fluid, p0 = 1e5, T0 = 1800)
    T, p = JTCT.trace()
    plt.plot(T, p, '-', label = 'Joule-Thomson Curve', **kwargs)

    print('Saturation Curve')
    Tt = ICT.AS.keyed_output(CoolProp.iT_triple)
    Tc = ICT.AS.keyed_output(CoolProp.iT_critical)
    Ts = np.linspace(Tt, Tc - 1.e-6)
    ps = CoolProp.CoolProp.PropsSI('P','T',Ts,'Q',0,backend + '::' + fluid)
    plt.plot(Ts, ps, '-', label = 'Saturation Curve', **kwargs)

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('T (K)')
    plt.ylabel('p (Pa)')
    plt.legend(loc = 'best')