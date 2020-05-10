.. _Humid-Air:

Humid Air Properties
********************

.. contents:: :depth: 2

If you are feeling impatient, jump to :ref:`HAProps_Sample`, or to go to the code documentation :mod:`CoolProp.HumidAirProp`, otherwise, hang in there.

The equations implemented in CoolProp are based on a publication by Hermann et al. :cite:`Herrmann2009`, which describes the outcome of the ASHRAE research project ASHREA-RP1485. 
The same source has been used in the ASHRAE Handbook 2009 to generate reference saturation property tables. The code implemented here passes all tests and reproduces the original 
data with a very high accuracy. It is applicable for pressure from 0.01 kPa up to 10 MPa, in a temperature range from -143.15 °C up to 350 °C with a humidity ratio from 0 kg of water 
up to 10 kg of water per kg of dry air.

Humid air can be modeled as a mixture of air and water vapor.  In the simplest analysis, water and air are treated as ideal gases but in principle there is interaction between the air and water molecules that must be included through the use of interaction parameters.

Because humid air is a mixture of dry air (treated as a pseudo-pure gas) and water vapor (treated as a real gas), three variables are required to fix the state by the state postulate.

In the analysis that follows, the three parameters that are ultimately needed to calculate everything else are the dry bulb temperature :math:`T`, the total pressure :math:`p`, and the molar fraction of water :math:`\psi_w`.  The molar fraction of air is simply :math:`\psi_a=1-\psi_w`.

Of course, it is not so straightforward to measure the mole fraction of water vapor molecules, so other measures are used.  There are three different variables that can be used to obtain the mole fraction of water vapor without resorting to iterative methods.

1. Humidity ratio

The humidity ratio :math:`W` is the ratio of the mass of water vapor to the mass of air in the mixture.  Thus the mole fraction of water can be obtained from

.. math::

    \psi_w=\frac{n_w}{n}=\frac{n_w}{n_a+n_w}=\frac{m_w/M_w}{m_a/M_a+m_w/M_w}=\frac{m_w}{(M_w/M_a)m_a+m_w}=\frac{1}{(M_w/M_a)/W+1}=\frac{W}{(M_w/M_a)+W}
    
or

.. math::

    \psi_w=\frac{W}{\varepsilon+W}

where the ratio of mole masses :math:`\varepsilon` is given by :math:`\varepsilon=M_w/M_a`

2. Relative Humidity

The relative humidity :math:`\varphi` is defined as the ratio of the mole fraction of water in the humid air to the saturation mole fraction of water.  Because of the presence of air with the water, the pure water saturated vapor pressure :math:`p_{w,s}` must be multiplied by an enhancement factor :math:`f` that is very close to one near atmospheric conditions.

Mathematically, the result is

.. math::

    \varphi=\frac{\psi_w}{\psi_{w,s}}

where 

.. math::

    \psi_{w,s}=\frac{fp_{w,s}}{p}
    
The product :math:`p_s` is defined by :math:`p_s=fp_{w,s}`, and :math:`p_{w,s}` is the saturation pressure of pure water (or ice) at temperature :math:`T`. This yields the result for :math:`\psi_w` of

.. math::

    \varphi=\frac{\psi_w}{p_s/p}
    
.. math::

    \psi_w=\frac{\varphi p_s}{p}

3. Dewpoint temperature

The dewpoint temperature is defined as the temperature at which the actual vapor pressure of water is equal to the saturation vapor pressure.  At the given dewpoint, the vapor pressure of water is given by

.. math::

    p_w=f(p,T_{dp})p_{w,s}(T_{dp})

and the mole fraction of water vapor is obtained from

.. math::

    \psi_w=\frac{p_w}{p}
    
Once the state has been fixed by a set of :math:`T,p,\psi_w`, any parameter of interest can be calculated

Molar Volume
------------
.. math::
    :label: eq1
    
    p=\frac{\bar R T}{\bar v}\left( 1+\frac{B_m}{\bar v}+\frac{C_m}{\bar v^2}\right)
    
The bracketed term on the right hand side is the compressibility Z factor, equal to 1 for ideal gas, and is a measure of non-ideality of the air.  The virial terms are given by
    
.. math::
    
    B_m=(1-\psi_w)^2B_{aa}+2(1-\psi_w)\psi_wB_{aw}+\psi_w^2B_{ww}
    
    C_m=(1-\psi_w)^3C_{aaa}+3(1-\psi_w)^2\psi_wC_{aaw}+3(1-\psi_w)\psi_w^2C_{aww}+\psi_w^3C_{www}
    
where the virial coefficients are described in ASRAE RP-1485 and their values are provided in :ref:`HA-Validation`.  All virial terms are functions only of temperature.

Usually the temperature is known, the water mole fraction is calculated, and :math:`\bar v` is found using iterative methods, in HAProps, using a secant solver and the first guess that the compressibility factor is 1.0.
    
Molar Enthalpy
--------------

The molar enthalpy of humid air is obtained from

.. math::

    \bar h=(1-\psi_w)\bar h_a^o+\psi_w\bar h_w^o+\bar R T \left[(B_m-T\frac{dB_m}{dT})\frac{1}{\bar v}+\left(C_m-\frac{T}{2}\frac{dC_m}{dT}\right) \frac{1}{\bar v^2}\right]

with :math:`\bar h` in kJ/kmol.  For both air and water, the full EOS is used to evaluate the enthalpy

.. math::

    \bar h_a^o=\bar h_0+\bar RT\left[ 1+\tau\left( \frac{\partial \alpha^o}{\partial \tau}\right)_{\delta}\right]

which is in kJ/kmol, using the mixture :math:`\bar v` to define the parameter :math:`\delta=1/(\bar v \bar \rho_c)` for each fluid, and using the critical molar density for the fluid obtained from :math:`\bar \rho_c=1000\rho_c/M` to give units of mol/m\ :sup:`3`\ .  The offset enthalpies for air and water are given by

.. math::
    
    \bar h_{0,a}=-7,914.149298\mbox{ kJ/kmol}
    
    \bar h_{0,w}=-0.01102303806\mbox{ kJ/kmol}
    
respectively.  The enthalpy per kg of dry air is given by

.. math::

    h=\bar h\frac{1+W}{M_{ha}}

Enhancement factor
------------------

The enhancement factor is a parameter that includes the impact of the air on the saturation pressure of water vapor.  It is only a function of temperature and pressure, but it must be iteratively obtained due to the nature of the expression for the enhancement factor.

:math:`\psi_{w,s}` is given by :math:`\psi_{w,s}=fp_{w,s}/p`, where :math:`f` can be obtained from 

.. math::

    \ln(f)=\left[ \begin{array}{l}\left [ \dfrac{(1+k_Tp_{w,s})(p-p_{w,s})-k_T\dfrac{(p^2-p_{w,s}^2)}{2}}{\overline {R} T}\right] \bar v_{w,s}+\ln[1-\beta_H(1-\psi_{w,s})p]\\
    +\left[\dfrac{(1-\psi_{w,s})^2p}{\bar R T}\right] B_{aa}-2\left[\dfrac{(1-\psi_{w,s})^2p}{\bar R T}\right]B_{aw}-\left[\dfrac{(p-p_{w,s}-(1-\psi_{w,s})^2p)}{\bar R T}\right]B_{ww} \\
    +\left[\dfrac{(1-\psi_{w,s})^3 p^2}{(\bar R T)^2}\right] C_{aaa}+\left[\dfrac{3(1-\psi_{w,s})^2[1-2(1-\psi_{w,s})]p^2}{2(\bar R T)^2}\right]C_{aaw}\\
    -\left[\dfrac{3(1-\psi_{w,s})^2\psi_{w,s}p^2}{(\bar R T)^2}\right]C_{aww}-\left[\dfrac{(3-2\psi_{w,s})\psi_{w,s}^2p^2-p_{w,s}^2}{2(\bar R T)^2}\right]C_{www}\\
    -\left[\dfrac{(1-\psi_{w,s})^2(-2+3\psi_{w,s})\psi_{w,s}p^2}{(\bar R T)^2}\right]B_{aa}B_{ww}\\
    -\left[\dfrac{2(1-\psi_{w,s})^3(-1+3\psi_{w,s})p^2}{(\bar R T)^2}\right]B_{aa}B_{aw}\\
    +\left[\dfrac{6(1-\psi_{w,s})^2\psi_{w,s}^2p^2}{(\bar R T)^2}\right]B_{ww}B_{aw}-\left[\dfrac{3(1-\psi_{w,s})^4p^2}{2(\bar R T)^2}\right]B_{aa}^2\\
    -\left[\dfrac{2(1-\psi_{w,s})^2\psi_{w,s}(-2+3\psi_{w,s})p^2}{(\bar R T)^2}\right]B_{aw}^2-\left[\dfrac{p_{w,s}^2-(4-3\psi_{w,s})(\psi_{w,s})^3p^2}{2(\bar R T)^2}\right]B_{ww}^2
    \end{array}\right]


Isothermal Compressibility
--------------------------

For water, the isothermal compressibility [in 1/Pa] is evaluated from

.. math::

    k_T=\frac{1}{\rho\frac{\partial p}{\partial \rho}}\frac{1\mbox{ kPa}}{1000\mbox{ Pa}}
    
with

.. math::

    \frac{\partial p}{\partial \rho}=RT\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]
    
in kPa/(kg/m\ :sup:`3`\ ). And for ice,

.. math::

    k_T=\left( \frac{\partial^2 g}{\partial p^2}\right) \left( \frac{\partial g}{\partial p}\right)_T^{-1}\frac{1\mbox{ kPa}}{1000\mbox{ Pa}}

.. _HAProps_Sample:

Sample HAPropsSI Code
---------------------
To use the HAPropsSI function, import it and do some calls, do something like this

.. ipython::

    #import the things you need 
    In [1]: from CoolProp.HumidAirProp import HAPropsSI
    
    #Enthalpy (J per kg dry air) as a function of temperature, pressure, 
    #    and relative humidity at dry bulb temperature T of 25C, pressure 
    #    P of one atmosphere, relative humidity R of 50%
    In [2]: h = HAPropsSI('H','T',298.15,'P',101325,'R',0.5); print(h)
    
    #Temperature of saturated air at the previous enthalpy
    In [2]: T = HAPropsSI('T','P',101325,'H',h,'R',1.0); print(T)
    
    #Temperature of saturated air - order of inputs doesn't matter
    In [2]: T = HAPropsSI('T','H',h,'R',1.0,'P',101325); print(T)

Table of Inputs/Outputs to HAPropsSI
------------------------------------

.. csv-table:: Input/Output parameters
    :header: "Parameter", "Units", "Input/Output", "Description"
    :widths: 25, 25, 25, 25
    :delim: ;
   
    ``B``, ``Twb``, ``T_wb``, ``WetBulb``; K; Input/Output; Wet-Bulb Temperature
    ``C``, ``cp``; J/kg dry air/K; Output; Mixture specific heat per unit dry air
    ``Cha``, ``cp_ha``; J/kg humid air/K; Output; Mixture specific heat per unit humid air
    ``CV``; J/kg dry air/K; Output; Mixture specific heat at constant volume per unit dry air
    ``CVha``, ``cv_ha``; J/kg humid air/K; Output; Mixture specific heat at constant volume per unit humid air
    ``D``, ``Tdp``, ``DewPoint``, ``T_dp``; K; Input/Output;	Dew-Point Temperature
    ``H``, ``Hda``, ``Enthalpy``;	J/kg dry air; Input/Output; Mixture enthalpy per dry air
    ``Hha``; J/kg humid air; Input/Output; Mixture enthalpy per humid air 
    ``K``, ``k``, ``Conductivity``; W/m/K; Output; Mixture thermal conductivity
    ``M``, ``Visc``, ``mu``;Pa-s;Output;Mixture viscosity
    ``psi_w``, ``Y``; mol water/mol humid air; Input/Output; Water mole fraction
    ``P``;Pa;Input;Pressure
    ``P_w``;Pa;Input;Partial pressure of water vapor
    ``R``, ``RH``, ``RelHum``; ; Input/Output; Relative humidity in [0, 1]
    ``S``, ``Sda``, ``Entropy``; J/kg dry air/K; Input/Output; Mixture entropy per unit dry air
    ``Sha``; J/kg humid air/K; Input/Output; Mixture entropy per unit humid air 
    ``T``, ``Tdb``, ``T_db``; K; Input/Output; Dry-Bulb Temperature
    ``V``, ``Vda``; m :math:`^3` /kg dry air; Input/Output;  Mixture volume per unit dry air
    ``Vha``; m :math:`^3` /kg humid air; Input/Output;  Mixture volume per unit humid air
    ``W``, ``Omega``, ``HumRat``; kg water/kg dry air; Input/Output; Humidity Ratio
    ``Z``; ; Output; Compressibility factor (:math:`Z = pv/(RT)`)
    
Psychrometric Chart
-------------------

.. plot::

    import numpy as np
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1,1,figsize=(10, 8))
    Tdbvec = np.linspace(-30, 55)+273.15

    # Lines of constant relative humidity
    for RH in np.arange(0.1, 1, 0.1):
        W = CP.HAPropsSI("W","R",RH,"P",101325,"T",Tdbvec)
        plt.plot(Tdbvec-273.15, W, color='k', lw = 0.5)

    # Saturation curve
    W = CP.HAPropsSI("W","R",1,"P",101325,"T",Tdbvec)
    plt.plot(Tdbvec-273.15, W, color='k', lw=1.5)

    # Lines of constant Vda
    for Vda in np.arange(0.69, 0.961, 0.01):
        R = np.linspace(0,1)
        W = CP.HAPropsSI("W","R",R,"P",101325,"Vda",Vda)
        Tdb = CP.HAPropsSI("Tdb","R",R,"P",101325,"Vda",Vda)
        plt.plot(Tdb-273.15, W, color='b', lw=1.5 if abs(Vda % 0.05) < 0.001 else 0.5)

    # Lines of constant wetbulb
    for Twb_C in np.arange(-16, 33, 2):
        if Twb_C == 0:
            continue
        R = np.linspace(0.0, 1)
        print(Twb_C)
        Tdb = CP.HAPropsSI("Tdb","R",R,"P",101325,"Twb",Twb_C+273.15)
        W = CP.HAPropsSI("W","R",R,"P",101325,"Tdb",Tdb)
        plt.plot(Tdb-273.15, W, color='r', lw=1.5 if abs(Twb_C % 10) < 0.001 else 0.5)

    plt.xlabel(r'Dry bulb temperature $T_{\rm db}$ ($^{\circ}$ C)')
    plt.ylabel(r'Humidity Ratio $W$ (kg/kg)')
    plt.ylim(0, 0.030)
    plt.xlim(-30, 55)
    # plt.show()

.. _HA-Validation:

Humid Air Validation
--------------------
Values here are obtained at documentation build-time using the Humid Air Properties module
 
.. ipython::

    In [1]: %run 'fluid_properties/Validation/HAValidation.py'
    

Verification Script
-------------------
This script, written in Python, should yield no failures::

    import CoolProp.CoolProp as CP
    import numpy as np
    import itertools
    from multiprocessing import Pool
    CP.set_config_bool(CP.DONT_CHECK_PROPERTY_LIMITS, True)

    def generate_values(TR,P=101325):
        """ Starting with T,R as inputs, generate all other values """
        T,R = TR
        psi_w = CP.HAPropsSI('psi_w','T',T,'R',R,'P',P)
        other_output_keys = ['T_wb','T_dp','Hda','Sda','Vda','Omega']
        outputs = {'psi_w':psi_w,'T':T,'P':P,'R':R}
        for k in other_output_keys:
            outputs[k] = CP.HAPropsSI(k,'T',T,'R',R,'P',P)
        return outputs

    def get_supported_input_pairs():
        """ Determine which input pairs are supported """
        good_ones = []
        inputs = generate_values((300, 0.5))
        for k1, k2 in itertools.product(inputs.keys(), inputs.keys()):
            if 'P' in [k1,k2] or k1==k2:
                continue
            args = ('psi_w', k1, inputs[k1], k2, inputs[k2], 'P', inputs['P'])
            try:
                psi_w_new = CP.HAPropsSI(*args)
                if not np.isfinite(psi_w_new):
                    raise ValueError('Returned NaN; not ok')
                good_ones.append((k1,k2))
            except BaseException as BE:
                pass
                if 'currently at least one of' in str(BE) or 'cannot provide two inputs' in str(BE):
                    pass
                else:
                    print(BE)
                    good_ones.append((k1,k2))
        return good_ones
    supported_pairs = get_supported_input_pairs()

    def calculate(inputs):
        """ For a given input, try all possible input pairs """
        errors = []
        for k1, k2 in supported_pairs:
            psi_w_input = inputs['psi_w']
            args = 'psi_w',k1,inputs[k1],k2,inputs[k2],'P',inputs['P']
            try:
                psi_w_new = CP.HAPropsSI(*args)
                if not np.isfinite(psi_w_new):
                    raise ValueError('Returned NaN; not ok')
            except BaseException as BE:
                errors.append((str(BE),args, inputs))
        return errors

    if __name__ == '__main__':
        import CoolProp
        print(CoolProp.__version__)
        TR = itertools.product(np.linspace(240, 345, 11), np.linspace(0, 1, 11))
        with Pool(processes=2) as pool:
            input_values = pool.map(generate_values, TR)
            errors = pool.map(calculate, input_values)
            for err in itertools.chain.from_iterable(errors):
                print(err)
    
..  Appendices


    Derivation for specific heat of humid air
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. math::

        p = \frac{{\bar RT}}{{\bar v}}\left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\\

    .. math::

        0 = \frac{d}{{dT}}\left[ {\frac{{\bar RT}}{{\bar v}}\left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)} \right]\\

    .. math::

        0 = \frac{{\bar RT}}{{\bar v}}\left( {0 + \frac{d}{{dT}}\left[ {\frac{{{B_m}}}{{\bar v}}} \right] + \frac{d}{{dT}}\left[ {\frac{{{C_m}}}{{{{\bar v}^2}}}} \right]} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\bar R\left( {\frac{{\bar v - T\frac{{d\bar v}}{{dT}}}}{{{{\bar v}^2}}}} \right)\\

    .. math::

        0 = \frac{{\bar RT}}{{\bar v}}\left( {0 + \frac{{\bar v\frac{{d{B_m}}}{{dT}} - {B_m}\frac{{d\bar v}}{{dT}}}}{{{{\bar v}^2}}} + \frac{{{{\bar v}^2}\frac{{d{C_m}}}{{dT}} - 2{C_m}\bar v\frac{{d\bar v}}{{dT}}}}{{{{\bar v}^4}}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\bar R\left( {\frac{{\bar v - T\frac{{d\bar v}}{{dT}}}}{{{{\bar v}^2}}}} \right)\\

    .. math::

        0 = \frac{{\bar RT}}{{\bar v}}\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} - \frac{{{B_m}}}{{{{\bar v}^2}}}\frac{{d\bar v}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}} - \frac{{2{C_m}}}{{{{\bar v}^3}}}\frac{{d\bar v}}{{dT}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\bar R\left( {\frac{1}{{\bar v}} - \frac{T}{{{{\bar v}^2}}}\frac{{d\bar v}}{{dT}}} \right)\\

    .. math::

        0 = \frac{T}{{\bar v}}\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} - \frac{{{B_m}}}{{{{\bar v}^2}}}\frac{{d\bar v}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}} - \frac{{2{C_m}}}{{{{\bar v}^3}}}\frac{{d\bar v}}{{dT}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\left( {\frac{1}{{\bar v}} - \frac{T}{{{{\bar v}^2}}}\frac{{d\bar v}}{{dT}}} \right)\\

    .. math::

        \frac{{d\bar v}}{{dT}}\left( {\frac{{{B_m}}}{{{{\bar v}^2}}}\frac{T}{{\bar v}} + \frac{{2{C_m}}}{{{{\bar v}^3}}}\frac{T}{{\bar v}} + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\frac{T}{{{{\bar v}^2}}}} \right) = \frac{T}{{\bar v}}\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\left( {\frac{1}{{\bar v}}} \right)\\

    .. math::

        \frac{{d\bar v}}{{dT}}\left( {\frac{{{B_m}}}{{{{\bar v}^2}}}T + \frac{{2T{C_m}}}{{{{\bar v}^3}}} + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\frac{T}{{\bar v}}} \right) = T\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\\

    .. math::

        \frac{{d\bar v}}{{dT}} = \frac{{T\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}}} \right) + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)}}{{\left( {\frac{{{B_m}}}{{{{\bar v}^2}}}T + \frac{{2T{C_m}}}{{{{\bar v}^3}}} + \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)\frac{T}{{\bar v}}} \right)}}\\

    .. math::

        \frac{{d\bar v}}{{dT}} = \frac{{T\left( {\frac{1}{{\bar v}}\frac{{d{B_m}}}{{dT}} + \frac{1}{{{{\bar v}^2}}}\frac{{d{C_m}}}{{dT}}} \right) + Z}}{{\left( {\frac{{{B_m}}}{{{{\bar v}^2}}}T + \frac{{2T{C_m}}}{{{{\bar v}^3}}} + Z\frac{T}{{\bar v}}} \right)}}\\

    .. math::

        \frac{{d\bar v}}{{dT}} = \frac{{\left( {\frac{{d{B_m}}}{{dT}} + \frac{1}{{\bar v}}\frac{{d{C_m}}}{{dT}}} \right) + Z\frac{{\bar v}}{T}}}{{\left( {\frac{{{B_m}}}{{\bar v}} + \frac{{2{C_m}}}{{{{\bar v}^2}}} + Z} \right)}}\\

    where

    .. math::

        Z = \left( {1 + \frac{{{B_m}}}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right)

    .. math::

        \bar h = {{\bar h}_0} + \left( {1 - {\psi _w}} \right)\bar h_a^0 + {\psi _w}\bar h_w^0 + \bar RT\left[ {\left( {{B_m} - T\frac{{d{B_m}}}{{dT}}} \right)\frac{1}{{\bar v}} + \left( {{C_m} - \frac{T}{2}\frac{{d{C_m}}}{{dT}}} \right)\frac{1}{{{{\bar v}^2}}}} \right]\\

    .. math::

        {{\bar c}_p} = \frac{{d\bar h}}{{dT}} = \frac{{\delta \bar h}}{{\delta T}} + \frac{{\delta \bar h}}{{\delta \bar v}}\frac{{\delta \bar v}}{{\delta T}}\\

    .. math::

        \frac{{\delta \bar h}}{{\delta \bar v}} = \left( {1 - {\psi _w}} \right)\frac{{d\bar h_a^0}}{{d\delta }}\frac{{d\delta }}{{d\bar v}} + {\psi _w}\frac{{d\bar h_w^0}}{{d\delta }}\frac{{d\delta }}{{d\bar v}} + \bar RT\left[ {\left( {{B_m} - T\frac{{d{B_m}}}{{dT}}} \right)\frac{{ - 1}}{{{{\bar v}^2}}} + \left( {{C_m} - \frac{T}{2}\frac{{d{C_m}}}{{dT}}} \right)\frac{{ - 2}}{{{{\bar v}^3}}}} \right]\\

    .. math::

        \frac{{\delta \bar h}}{{\delta T}} = \left( {1 - {\psi _w}} \right)\frac{{d\bar h_a^0}}{{d\tau }}\frac{{d\tau }}{{dT}} + {\psi _w}\frac{{d\bar h_w^0}}{{d\tau }}\frac{{d\tau }}{{dT}} + \bar R\left[ {\left( {{B_m} - T\frac{{d{B_m}}}{{dT}}} \right)\frac{1}{{\bar v}} + \left( {{C_m} - \frac{T}{2}\frac{{d{C_m}}}{{dT}}} \right)\frac{1}{{{{\bar v}^2}}}} \right] + \bar RT\left[ {\left( {\frac{{d{B_m}}}{{dT}} - \frac{{d{B_m}}}{{dT}} - T\frac{{{d^2}{B_m}}}{{d{T^2}}}} \right)\frac{1}{{\bar v}} + \left( {\frac{{d{C_m}}}{{dT}} - \frac{1}{2}\frac{{d{C_m}}}{{dT}} - \frac{T}{2}\frac{{{d^2}{C_m}}}{{d{T^2}}}} \right)\frac{1}{{{{\bar v}^2}}}} \right]\\

    .. math::

        \frac{{\delta \bar h}}{{\delta T}} = \left( {1 - {\psi _w}} \right)\frac{{d\bar h_a^0}}{{d\tau }}\frac{{d\tau }}{{dT}} + {\psi _w}\frac{{d\bar h_w^0}}{{d\tau }}\frac{{d\tau }}{{dT}} + \bar R\left[ {\left( {{B_m} - T\frac{{d{B_m}}}{{dT}}} \right)\frac{1}{{\bar v}} + \frac{{{C_m}}}{{{{\bar v}^2}}}} \right] + \bar R{T^2}\left[ {\left( { - \frac{{{d^2}{B_m}}}{{d{T^2}}}} \right)\frac{1}{{\bar v}} + \left( { - \frac{1}{2}\frac{{{d^2}{C_m}}}{{d{T^2}}}} \right)\frac{1}{{{{\bar v}^2}}}} \right]
