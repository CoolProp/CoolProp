

.. _Incompressibles:

Incompressible Fluids
=====================


General Introduction
--------------------

In CoolProp, the incompressible fluids are divided into three major groups.

* :ref:`Pure fluids <Pure>`.
* :ref:`Mass-based binary mixtures <MassMix>`.
* :ref:`Volume-based binary mixtures <VoluMix>`.

.. * :ref:`Mole-based binary mixtures <MoleMix>`.

The pure fluids and mass-based binary mixtures are by far the most common fluids
in this library. While the pure fluids contain data for many different kinds of
incompressible liquids, almost all of the binary mixtures are aqueous solutions.
For these liquids, the concentration always refers to the added component ranging
from 0.0 for pure water to 1.0 for no water at all. Please refer to the tables
below to find the allowed minimum and maximum concentrations. Those are likely
to be above 0.0 and below 1.0, respectively.

The first entry in the tables below is the fluid ID that can be used to call the
fluid from the high-level interface. A single PDF page showing the fit quality is
linked to that ID in case you would like to see a few more details about any
specific fluid. To get an overview over all the fits, there are also combined
documents with all the
:download:`pure fluids and all the aqueous solutions</_static/fluid_properties/incompressible/report/all_incompressibles.pdf>`.
You can read more about these reports in a dedicated
:ref:`section<FittingReports>` called :ref:`Fitting Reports<FittingReports>` below.

All incompressible fluids have an arbitrary reference state for enthalpy and entropy.
During initialisation, the reference state is defined as a temperature of 20 °C
and a pressure of 1 atm according to the U.S. National Institute of Standards and
Technology ([NIST](http://www.nist.gov)).

.. math::
   T_\text{ref} &=  293.15\,\text{K}  &= 68\,\text{°F} \\
   p_\text{ref} &=  101325\,\text{Pa} &= 14.696\,\text{psi} \\
   h_\text{ref} &=  0 & \\
   s_\text{ref} &=  0 & \\

If you use a mixture, the reference state gets updated each time you change the
composition.


Pure Fluid Examples
-------------------

Incompressible fluids only allow  for a limited subset of input variables. The
following input pairs are supported: :math:`f(p,T)`, :math:`f(p,h)`, :math:`f(p,\rho)`,
:math:`f(p,u)` and :math:`f(p,s)`. Some fluids also provide saturation state
information as :math:`f(Q,T)` with :math:`Q=0`. All functions iterate on :math:`f(p,T)` calls
internally, which makes this combination by far the fastest. However, also the
other inputs should be fast compared to the full Helmholtz-based EOS implemented
for then compressible fluids.

A call to the top-level function ``PropsSI`` can provide : density, heat capacity,
internal energy, enthalpy, entropy, viscosity and thermal conductivity. Hence,
the available output keys are: ``D``, ``C``, ``U``, ``H``, ``S``, ``V``, ``L``,
``Tmin`` and ``Tmax``.

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI

    #Density of Downtherm Q at 500 K and 1 atm.
    In [1]: PropsSI('D','T',500,'P',101325,'INCOMP::DowQ')

    #Specific heat capacity of Downtherm Q at 500 K and 1 atm
    In [1]: PropsSI('C','T',500,'P',101325,'INCOMP::DowQ')

    In [1]: PropsSI('C','D',809.0659,'P',101325,'INCOMP::DowQ')

    #Saturation pressure of Downtherm Q at 500 K
    In [1]: PropsSI('P','T',500,'Q',0,'INCOMP::DowQ')

    #Minimum temperature for Downtherm Q
    In [1]: PropsSI('Tmin','T',0,'P',0,'INCOMP::DowQ')

    #Maximum temperature for Downtherm Q
    In [1]: PropsSI('Tmax','T',0,'P',0,'INCOMP::DowQ')



Mixture Examples
----------------

Almost the same syntax can be used for mixtures. Please note that the mixture
interface developed for CoolProp 5 has not been ported to the incompressible
fluids, yet. For now, you have to use the ``PropsSI`` function with a special
composition notation. Depending on your fluid, you have to supply either the
:ref:`mass fraction<MassMix>` or the :ref:`volume fraction<VoluMix>` as additional
parameter. This is done via the fluid name by appending a dash and the
fraction of the substance other than water. The fraction notation can be in the
form of percent, ``LiBr-23%``, or as a fraction like in ``LiBr-0.23`` or
``LiBr[0.23]``, which corresponds to the new mixture syntax in CoolProp5.

..  In addition to the properties available for the pure fluids (``D``, ``C``,
  ``U``, ``H``, ``S``, ``V``, ``L``,``Tmin`` and ``Tmax``, some mixtures also
  provide the freezing temperature ``Tfreeze`` as a function of composition.


.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI

    #Density of a lithium bromide solution at 300 K and 1 atm.
    In [1]: PropsSI('D','T',300,'P',101325,'INCOMP::LiBr[0.23]')

    #Specific heat capacity of a lithium bromide solution at 300 K and 1 atm
    In [1]: PropsSI('C','T',300,'P',101325,'INCOMP::LiBr-0.23%')









.. _FittingReports:

Fitting Reports
---------------------------------------

A file with all fitting reports for the incompressible fluids can be obtained
from :download:`here </_static/fluid_properties/incompressible/report/all_incompressibles.pdf>`. These reports help you to
get an overview over the different incompressible fluids
included in CoolProp. The reports start with some basic information about
the fluid. There is the name by which it can be accessed through the
interface in the title "Fitting Report for *FluidName*" as well as a description
of what the fluid actually is, this could also be a trade name or a commonly
used non-scientific name. The next item tells you where we got the data from. This
would typically be a data sheet from a manufacturers homepage, some other software
database, a scientific publication or experimental data.

.. figure:: /_static/fluid_properties/incompressible/report/report2up.jpg
    :align: center
    :alt: Fitting reports for pure fluid and solution

    The figure above shows two examples for fitting reports generated for a pure
    fluid and a binary mixture. You can also have a look at the
    :download:`PDF version</_static/fluid_properties/incompressible/report/report2up.pdf>` of the reports side by side.

If all data is available, there is a graphs for each of the basic quantities
density :math:`\rho`, specific heat capacity :math:`c`, thermal conductivity
:math:`\lambda`, dynamic viscosity :math:`\mu`, saturation pressure
:math:`p_{sat}`, and freezing temperature :math:`T_{freeze}`. These graphs show
data points in dark blue, the fitted function from CoolProp as a red line and the
relative error in light blue dots. Note that the relative error uses the ordinate
on the right hand side while the other two data series refer to the axis on the
left hand side. In case of a solution, these graphs refer to a given concentration
that typically lies in the middle of the allowed range. Dashed red lines indicate
the limits in terms of concentration as well as the freezing temperature.



Equations
----------------------

Exp or log for visc, other poly or log poly

.. math::

    x(T) &= \sum_{i=0}^n C[i] \cdot T^i \\
    x(T) &= \exp\left( \frac{C[0]}{T+C[1]} - C[2] \right) \\
    x(T) &= \exp\left( \log  \left( \left(T+C[0]\right)^{-1} + \left( T+C[0] \right)^{-2} \right) *C[1]+C[2] \right) \\


All fluids are implemented with polynomials for density and heat capacity with typically 4 coefficients
and hence a third order polynomial. Thermal conductivity is a second order polynomial and viscosity and
vapour pressure are exponential functions.

.. math::

    \rho    &= \sum_{i=0}^n C_{\rho}[i] \cdot T^i \\
    c        &= \sum_{i=0}^n C_{c}[i] \cdot T^i \\
    u        &= \int_{0}^{1} c\left( T \right) dT
              = \sum_{i=0}^n \frac{1}{i+1} \cdot C_{c}[i]
                \cdot \left( T_1^{i+1} - T_0^{i+1} \right) \\
    s        &= \int_{0}^{1} \frac{c\left( T \right)}{T} dT
              = C_{c}[0] \cdot \ln\left(\frac{T_1}{T_0}\right)
                + \sum_{i=0}^{n-1} \frac{1}{i+1} \cdot C_{c}[i+1]
                \cdot \left( T_1^{i+1} - T_0^{i+1} \right) \\
    \lambda &= \sum_{i=0}^n C_{\lambda}[i] \cdot T^i \\
    \mu     &= \exp\left( \frac{C_{\mu}[0]}{T+C_{\mu}[1]} - C_{\mu}[2] \right) \\
    p_{sat}  &= \exp\left( \frac{C_{sat}[0]}{T+C_{sat}[1]} - C_{sat}[2] \right) \\

In some cases, the fit quality for the


Brines and Solutions
--------------------
All the brines and solutions can be accessed through the Props function. To use them, the fluid name
is something like ``"MEG-20%"`` which is a 20% by mass ethylene glycol solution. Note that these fluids
have an arbitrary reference state: Be careful with enthalpy and entropy calculations. Again, only
temperature and pressure inputs are supported directly to calculate the same subset of thermophysical
properties as above , namely: density, heat capacity, internal energy, enthalpy, entropy, viscosity
and thermal conductivity. Hence, the available output keys for the ``Props`` function are: "D", "C",
"U", "H", "S", "V", "L", "Tmin", Tmax" and "Tfreeze". An internal iteration allows us to use enthalpy
and pressure as inputs, but be aware of the reduced computational efficiency.

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI

    #Specific heat 20% mass ethylene glycol solution at 300 K and 1 atm.
    In [1]: PropsSI('C','T',300,'P',101.325,'INCOMP::MEG-20%')

For Lithium-Bromide, the publication by Patek and Klomfar from 2005 was implemented based on the
source code provided by the authors. The `paper <http://dx.doi.org/10.1016/j.ijrefrig.2005.10.007>`_
covering the equations can be found in the
`International Journal of Refrigeration <http://dx.doi.org/10.1016/j.ijrefrig.2005.10.007>`_. Data is
available for temperatures from 0 C to 225 C and for the full composition range. Use ``LiBr`` to acccess
the functions.

A number of aqueous solutions are implemented using the coefficients from Aake Melinder "Properties of
Secondary Working Fluids for Indirect Systems" published in 2010 by IIR.  According to the book, 2D
polynomials are given in a form that satisfies :math:`0 \leq i \leq 5`, :math:`0 \leq j \leq 3`
and :math:`i + j \leq 5` yielding a triangular matrix of coefficients. It is only for the freezing
temperature calculation that the implemented procedures differ from what is presented in Melinder's
book the dependency on the current temperature is removed. In CoolProp, :math:`T_{freeze}` only depends
on concentration.

==========================   ===================================================   =================   =================
Melinder Fluids              Description                                           max. T              max. x
==========================   ===================================================   =================   =================
``MEG``                      Ethylene Glycol (C2H6O2)                              +100 C              60 %
``MPG``                      Propylene Glycol (C3H8O2)                             +100 C              60 %
``MEA``                      Ethyl Alcohol, Ethanol (C2H6O)                        +40 C               60 %
``MMA``                      Methyl Alcohol, Methanol (CH4O)                       +40 C               60 %
``MGL``                      Glycerol (C3H8O3)                                     +40 C               60 %
``MAM``                      Ammonia (NH3)                                         +30 C               30 %
``MKC``                      Potassium Carbonate (K2CO3)                           +40 C               40 %
``MCA``                      Calcium Chloride (CaCl2)                              +40 C               30 %
``MMG``                      Magnesium Chloride (MgCl2)                            +40 C               30 %
``MNA``                      Sodium Chloride (NaCl)                                +40 C               23 %
``MKA``                      Potassium Acetate (CH3CO2K)                           +40 C               45 %
``MKF``                      Potassium Formate (CHKO2)                             +40 C               48 %
``MLI``                      Lithium Chloride (LiCl)                               +40 C               24 %
==========================   ===================================================   =================   =================

Furthermore, there is a number of other secondary fluids that can be accessed in the same way. Most
information is based on the data compiled by Morten Juel Skovrup in his `SecCool software <http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx>`_
provided by his employer `IPU <http://en.ipu.dk>`_. The coefficient matrix of the SecCool-based fluids
has the same structure as mentioned above.

For slurry ice, the concentration :math:`x` refers to the solid content and the heat capacity includes the heat of fusion.
It might be necessary to adjust the solid content during heat transfer. The implementation is based on the data available
in SecCool, which was originally recorded at the `Danish Technological Institute (DTI) <http://www.dti.dk/>`_.

==========================   ===================================================   =================   =================
SecCool Fluids               Description                                           max. T              max. x
==========================   ===================================================   =================   =================
``ZiAC``                     ZitrecAC (corrosion inhibitor)                        +100 C              50 %
``IceEA``                    Ethanol-water mixture with slurry ice                 -10 C               35 %
``IcePG``                    Propylene glycol-water mixture with slurry ice        -10 C               35 %
``IceNA``                    Sodium chloride-water mixture with slurry ice         -5 C                35 %
``PK2000``                   Pekasol 2000 (Potassium acetate and formate)          +100 C              100 %
==========================   ===================================================   =================   =================



In both of the above cases, :math:`i` is the exponent for the concentration :math:`x` and :math:`j`
is used with the temperature :math:`T`. Properties are modelled with the following polynomials:

.. math::

    \rho      &= \sum_{i=0}^n x^i  \cdot \sum_{j=0}^m C_{\rho}[i,j] \cdot T^j \\
    c          &= \sum_{i=0}^n x^i  \cdot \sum_{j=0}^m C_{c}[i,j] \cdot T^j \\
    u          &= \int_{0}^{1} c\left( x,T \right) dT
                = \sum_{i=0}^n x^i \cdot \sum_{j=0}^m \frac{1}{j+1} \cdot C_{c}[i,j]
                  \cdot \left( T_1^{j+1} - T_0^{j+1} \right) \\
    s          &= \int_{0}^{1} \frac{c\left( x,T \right)}{T} dT
                = \sum_{i=0}^n x^i \cdot \left(
                  C_{c}[i,0] \cdot \ln\left(\frac{T_1}{T_0}\right)
                  + \sum_{j=0}^{m-1} \frac{1}{j+1} \cdot C_{c}[i,j+1] \cdot \left( T_1^{j+1} - T_0^{j+1} \right)
                  \right) \\
    \lambda   &= \sum_{i=0}^n x^i  \cdot \sum_{j=0}^m C_{\lambda}[i,j] \cdot T^j \\
    \mu       &= \exp \left( \sum_{i=0}^n x^i  \cdot \sum_{j=0}^m C_{\mu}[i,j] \cdot T^j \right) \\
    T_{freeze} &= \sum_{i=0}^n C_{freeze}[i] \cdot x^i \\

Using a centered approach for the independent variables,
the fit quality can be enhanced. Therefore, all solutions have a reference temperature and concentration
in the original work by Melinder and Skovrup as well as in CoolProp: :math:`x = x_{real} - x_{ref}`
and :math:`T = T_{real} - T_{ref}`, this technique does not affect the calculation
of the derived quantity internal energy since the formula contains temperature differences.
However, integrating :math:`c(x,T)T^{-1}dT` for the entropy requires some changes due to
the logarithm. To structure the problem, we introduce a variable :math:`d(j,T_{real})`, which will be expressed by a third sum.
As a first step for simplification, one has to expand the the binomial :math:`(T_{real}-T_{ref})^n` to a series.
Only containing :math:`j` and :math:`T_{real}`, :math:`d` is independent from :math:`x` and can be
computed outside the loop for enhanced computational efficiency. An integration of the expanded binomial
then yields the final factor :math:`D` to be multiplied with the other coefficients and the concentration.

.. math::

    s          &= \int_{0}^{1} \frac{c\left( x,T \right)}{T} dT = \sum_{i=0}^n x^i \cdot \sum_{j=0}^m C_{c}[i,j] \cdot D(j,T_0,T_1) \\
    D          &= (-1)^j \cdot \ln \left( \frac{T_1}{T_0} \right) \cdot T_{ref}^j + \sum_{k=0}^{j-1} \binom{j}{k} \cdot \frac{(-1)^k}{j-k} \cdot \left( T_1^{j-k} - T_0^{j-k} \right) \cdot T_{ref}^k







The Different Fluids
--------------------

The fluids implemented in CoolProp cover a wide range of industrial heat
transfer media. This database has initially been developed with refrigeration
systems in mind. That is why the majority of fluids are secondary refrigerants
with application temperatures close to the freezing point of water. Besides those,
there is also incompressible water, high temperature heat transfer oils and a
molten salt mixture for extreme temperatures.

Besides the different technical data sheets and calculation tools provided by
manufactures, two specific publications provided a lot of data used for the
incompressible fluids: Åke Melinder's book *Properties of Secondary Working
Fluids for Indirect Systems* :cite:`Melinder2010` has inspired both, the work on
pure fluids and aqueous solutions. The second major source of inspiration is the
`SecCool software <http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx>`_
:cite:`Skovrup2013` software, which contains data compiled by Morten Juel
Skovrup. It is provided free of charge by his employer `IPU <http://en.ipu.dk>`_.


.. _Pure:

.. csv-table:: All incompressible pure fluids included in CoolProp
   :widths: 10, 35, 15, 20, 20
   :header-rows: 1
   :file: ../_static/fluid_properties/incompressible/table/pure-fluids.csv


There are also a number of water-based mixtures implemented in CoolProp. Most of them
are secondary heat transfer fluids, but there are also aqueous solutions of
ammonia :cite:`Melinder2010`, :download:`MAM<../_static/fluid_properties/incompressible/report/MAM_fitreport.pdf>`,
and lithium bromide :cite:`Patek2006`, :download:`LiBr<../_static/fluid_properties/incompressible/report/LiBr_fitreport.pdf>`.


.. _MassMix:

.. csv-table:: All incompressible mass-based binary mixtures included in CoolProp
   :widths: 10, 30, 12, 12, 12, 12, 12
   :header-rows: 1
   :file: ../_static/fluid_properties/incompressible/table/mass-based-fluids.csv

.. .. _MoleMix:

.. .. csv-table:: All incompressible mole-based binary mixtures included in CoolProp
   :widths: 10, 30, 12, 12, 12, 12, 12
   :header-rows: 1
   :file: ../_static/fluid_properties/incompressible/table/mole-based-fluids.csv

.. _VoluMix:

.. csv-table:: All incompressible volume-based binary mixtures included in CoolProp
   :widths: 10, 30, 12, 12, 12, 12, 12
   :header-rows: 1
   :file: ../_static/fluid_properties/incompressible/table/volume-based-fluids.csv


For slurry ice, the concentration :math:`x` refers to the solid content and the
heat capacity includes the heat of fusion. It might be necessary to adjust the
solid content during heat transfer. The implementation is based on the data
available in `SecCool<http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx>`_,
which was originally recorded at the `Danish Technological Institute (DTI) <http://www.dti.dk/>`_.


References
----------

.. bibliography:: Incompressibles.bib
   :filter: docname in docnames
   :style: unsrt
