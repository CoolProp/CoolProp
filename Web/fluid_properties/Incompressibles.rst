

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
:download:`pure fluids and all the aqueous solutions</_static/fluid_properties/Incompressibles_reports/all_incompressibles.pdf>`.
You can read more about these reports in a dedicated
:ref:`section<FittingReports>` called :ref:`Fitting Reports<FittingReports>`.

All incompressible fluids have an arbitrary reference state for enthalpy and entropy.
During initialisation, the reference state is defined as a temperature of 20 °C
and a pressure of 1 atm according to the U.S. National Institute of Standards and
Technology (`NIST <http://www.nist.gov>`_).

.. math::
   T_\text{ref} &=  293.15\:\text{K}  &=     68\:\text{°F} \\
   p_\text{ref} &=  101325\:\text{Pa} &= 14.696\:\text{psi} \\
   h_\text{ref} &=  0\:\text{J}\,\text{kg}^{-1} & \\
   s_\text{ref} &=  0\:\text{J}\,\text{kg}^{-1}\,\text{K}^{-1} & \\

.. note::
  If you use a mixture, the reference state gets updated each time you change
  the composition. Furthermore, not all temperatures can be used as reference
  temperature since the fraction :math:`T_\text{in,1} / T_\text{in,0}` occurs in the integral used to
  calculate entropy. The centered fits have a base temperature and setting
  :math:`T_\text{ref}` equal to :math:`T_\text{base}` yields :math:`T_\text{in,0}=0\:\text{K}`,
  which obviously is a problem. For non-centred fits, the base temperature is
  equal to 0 K. Read on :ref:`below<BaseValue>` for more details.


Pure Fluid Examples
-------------------

Incompressible fluids only allow  for a limited subset of input variables. The
following input pairs are supported: :math:`f(p,T)`, :math:`f(p,h)`, :math:`f(p,\rho)` 
and :math:`f(p,s)`. Some fluids also provide saturation state
information as :math:`f(Q,T)` with :math:`Q=0`. All functions iterate on :math:`f(p,T)` calls
internally, which makes this combination by far the fastest. However, also the
other inputs should be fast compared to the full Helmholtz-based EOS implemented
for then compressible fluids.

A call to the top-level function ``PropsSI`` can provide: temperature, pressure,
density, heat capacity, internal energy, enthalpy, entropy, viscosity and
thermal conductivity. Hence, the available output keys are: ``T``, ``P``, ``D``,
``C``, ``U``, ``H``, ``S``, ``V``, ``L``, ``Tmin`` and ``Tmax``.

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Specific heat capacity of Downtherm Q at 500 K and 1 atm
    In [1]: PropsSI('C','T',500,'P',101325,'INCOMP::DowQ')

    #Density of Downtherm Q at 500 K and 1 atm.
    In [1]: PropsSI('D','T',500,'P',101325,'INCOMP::DowQ')
    
    #Round trip in thermodynamic properties
    In [1]: T_init = 500.0
    
    In [2]: P_init = 101325
    
    In [3]: D_init = PropsSI('D','T',T_init,'P',P_init,'INCOMP::DowQ')
    
    In [4]: S_init = PropsSI('S','D',D_init,'P',P_init,'INCOMP::DowQ')
    
    In [5]: H_init = PropsSI('H','S',S_init,'P',P_init,'INCOMP::DowQ')
    
    In [6]: T_init = PropsSI('T','H',H_init,'P',P_init,'INCOMP::DowQ')
    
    In [7]: T_init

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
form of percent, ``LiBr-23%``, or as a fraction between 0 and 1, ``LiBr[0.23]``, which
corresponds to the new mixture syntax in CoolProp v5.

..  In addition to the properties available for the pure fluids (``D``, ``C``,
  ``U``, ``H``, ``S``, ``V``, ``L``,``Tmin`` and ``Tmax``, some mixtures also
  provide the freezing temperature ``Tfreeze`` as a function of composition.


.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI

    #Density of a lithium bromide solution at 300 K and 1 atm.
    In [1]: PropsSI('D','T',300,'P',101325,'INCOMP::LiBr[0.23]')

    #Density of a lithium bromide solution at 300 K and 1 atm.
    In [1]: PropsSI('D','T',300,'P',101325,'INCOMP::LiBr-23%')

    #Specific heat capacity of a lithium bromide solution at 300 K and 1 atm
    In [1]: PropsSI('C','T',300,'P',101325,'INCOMP::LiBr-23%')

    #Specific enthalpy of a lithium bromide solution at 300 K and 1 atm
    In [1]: PropsSI('H','T',300,'P',101325,'INCOMP::LiBr-23%')

    In [1]: PropsSI('T','H',28627,'P',101325,'INCOMP::LiBr-23%')


.. warning::
  Some mixture function have a non-monotonic behaviour, this can lead to misleading
  results when using other inputs than :math:`f(p,T)`. Keep that in mind and
  implement a way to validate the results you get from these functions. At the same
  time, mixture solvers are likely to produce errors due to the same reason...


  
Partial Derivatives
-------------------

A limited subset of partial derivatives is available for the incompressible fluids. Currently, 
the following inputs are supported by the ``PropsSI`` function: 
:math:`\left( \partial \rho / \partial p \right)_{T,x}=0`,
:math:`\left( \partial \rho / \partial h \right)_{p,x}`,
:math:`\left( \partial \rho / \partial s \right)_{p,x}`,
:math:`\left( \partial \rho / \partial T \right)_{p,x}`, 
:math:`\left( \partial h    / \partial p \right)_{T,x}`,
:math:`\left( \partial h    / \partial s \right)_{T,x}`,
:math:`\left( \partial h    / \partial T \right)_{p,x}`,
:math:`\left( \partial s    / \partial p \right)_{T,x}`,
:math:`\left( \partial s    / \partial T \right)_{p,x}`, 
and their inverse functions. 

Note that all partial derivatives require a constant concentration, which is denoted by the 
:math:`x`, but this :math:`x` is not included in the derivative string notation for ``PropsSI``: 
:math:`\left( \partial \rho / \partial T \right)_{p,x}` translates to ``d(Dmass)/d(T)|P``.

.. note::
   You can calculate other properties from the partial derivatives available. At this point, not
   all derived properties have been implemented even though some of them can be computed like the 
   isobaric expansion coefficient, which would be :math:`-\left( \partial \rho / \partial T \right)_{p,x}/\rho`.

For more general information on the partial derivatives, please have a look at the 
:ref:`documentation<partial_derivatives_high_level>` for the high level interface.


.. _FittingReports:

Fitting Reports
---------------------------------------

A file with all fitting reports for the incompressible fluids can be obtained
from :download:`here</_static/fluid_properties/Incompressibles_reports/all_incompressibles.pdf>`. These reports help you to
get an overview over the different incompressible fluids
included in CoolProp. The reports start with some basic information about
the fluid. The fluid name used in CoolProp is in the title "Fitting Report for *FluidName*"
and there is also a description of what the fluid actually consists of. The latter
could also be a trade name or a commonly used non-scientific name. The next item
tells you where we got the data from. This
would typically be a data sheet from a manufacturer's homepage, some other software
database, a scientific publication or experimental data.

.. figure:: /_static/fluid_properties/Incompressibles_reports/report2up.jpg
    :align: center
    :alt: Fitting reports for pure fluid and solution

    The figure above shows two examples for fitting reports generated for a pure
    fluid and a binary mixture. You can also have a look at the
    :download:`PDF version</_static/fluid_properties/Incompressibles_reports/report2up.pdf>` of the reports side by side.

If all data are available, there is a graph for each of the basic quantities:
density :math:`\rho`, specific heat capacity :math:`c`, thermal conductivity
:math:`\lambda`, dynamic viscosity :math:`\mu`, saturation pressure
:math:`p_\text{sat}`, and freezing temperature :math:`T_\text{freeze}`. These graphs show
data points in dark blue, the fitted function from CoolProp as a red line and the
relative error in light blue dots. Note that the relative error uses the ordinate
on the right hand side while the other two data series refer to the axis on the
left hand side. In case of a solution, these graphs refer to a given concentration
that typically lies in the middle of the allowed range. Dashed red lines indicate
the limits in terms of concentration as well as the freezing temperature.


.. _IncompThermo:

Thermodynamics of Incompressible Fluids
---------------------------------------

For an incompressible fluid, the specific at constant volume and at constant
pressure are the same allowing us to drop the subscripts, :math:`c_p=c_v=c`. Using
temperature :math:`T` and pressure :math:`p` as state variables, we can simplify
the normal thermodynamic relation as described below. working with brines and
mixtures, the concentration :math:`x` has to be considered as well. Following
the same approach as for the compressible fluids, we regard mixtures with different
compositions as independent fluids. This should be kept in mind when comparing
properties for different compositions. Setting the reference state for one
composition will always affect all fluids consisting of the same components.

.. The approach described in textbooks like Cengel and Boles :cite:`Cengel2007`
.. is that the internal energy :math:`u` only depends on temperature and does not
.. change with pressure.
.. 
.. .. Alternatively, use cancel package with \cancelto{0}{x-d} command
.. 
.. .. math::
.. 
..     du &= \overbrace{ \left( \frac{\partial u}{\partial T} \right)_p}^{=c_p=c_v=c} dT &+ \overbrace{\left( \frac{\partial u}{\partial p} \right)_T}^{\stackrel{!}{=}0} dp \\
.. 
.. By using the fourth Maxwell relation, we can extend the simplifications to the
.. entropy formulation
.. 
.. .. math::
.. 
..     ds &= \left( \frac{\partial s}{\partial T} \right)_p dT &+              \left( \frac{\partial s}{\partial p} \right)_T  dp \\
..        &= \underbrace{ \left( \frac{\partial h}{\partial T} \right)_p}_{=c_p=c_v=c} T^{-1} dT &-\underbrace{\left( \frac{\partial v}{\partial T} \right)_p}_{\stackrel{!}{=} 0} dp \\
.. 
.. As indicated by the braces above, the fluids implemented in CoolProp do also follow
.. the second common assumption of a constant specific volume :math:`v` that does
.. change neither with temperature nor with pressure. It should be highlighted, that
.. this simplification violates the integrity of the implemented equations since there
.. are changes in density as a function of temperature for all incompressible fluids.
.. 
.. Employing :math:`h=u+pv`, we can derive the impact on enthalpy as well by
.. rewriting the equation in terms of our state variables :math:`p` and :math:`T`
.. as shown by Skovrup :cite:`Skovrup1999`.
.. 
.. .. dh &= \overbrace{ \left( \frac{\partial h}{\partial T} \right)_p}^{=c_p=c_v=c} dT + \left( \frac{\partial h}{\partial p} \right)_T dp \\
.. 
.. .. math::
..     dh &= \overbrace{ \left( \frac{\partial h}{\partial T} \right)_p}^{=c_p=c_v=c} dT +              \left( \frac{\partial h}{\partial p} \right)_T         dp \\
..        &=             \left( \frac{\partial u}{\partial T} \right)_v dT               + \left( v - T \left( \frac{\partial v}{\partial T} \right)_p \right) dp \\
..        &= du + \underbrace{p dv}_{\stackrel{!}{=} 0} + v dp \quad \text{ with $v\stackrel{!}{=}v_0=$ const }  \\
.. 
.. The two assumptions used above :math:`\left( \partial v / \partial T \right)_p \stackrel{!}{=} 0`
.. and :math:`\left( \partial u / \partial T \right)_p \stackrel{!}{=} \left( \partial u / \partial T \right)_v`
.. imply that :math:`v` is constant under all circumstances. Hence, we have to use
.. the specific volume at reference conditions to calculate enthalpy from the
.. integration in :math:`T` and :math:`p`. Future work could provide a more accurate
.. formulation of entropy and enthalpy by implementing the term
.. :math:`\left( \partial v / \partial T \right)_p \neq 0`.
.. 
.. Using only polynomials for the heat capacity functions, we can derive internal
.. energy and entropy by integrating the specific heat capacity in temperature.

.. _BaseValue:

.. note::
   The internal routines for the incompressibles were updated 2015-02-10, the documentation is not fully updated. 
   We are going to add the new equation as soon as possible, probably mid-March 2015. Please be patient.

.. .. math::
.. 
..     c          &= \sum_{i=0}^n x^i \cdot \sum_{j=0}^m C_{c}[i,j] \cdot T^j \text{ yielding } \\
..     u          &= \int_{0}^{1} c\left( x,T \right) dT
..                 = \sum_{i=0}^n x^i \cdot \sum_{j=0}^m \frac{1}{j+1} \cdot C_{c}[i,j]
..                   \cdot \left( T_{1}^{j+1} - T_{0}^{j+1} \right) \text{ and } \\
..     s          &= \int_{0}^{1} \frac{c\left( x,T \right)}{T} dT
..                 = \sum_{i=0}^n x^i \cdot \left(
..                   C_{c}[i,0] \cdot \ln\left(\frac{T_{1}}{T_{0}}\right)
..                   + \sum_{j=0}^{m-1} \frac{1}{j+1} \cdot C_{c}[i,j+1] \cdot \left( T_{1}^{j+1} - T_{0}^{j+1} \right)
..                   \right) \\
..     h          &= u + v_{0} \cdot \left( p_{1} - p_{0} \right)
.. 

According to Melinder :cite:`Melinder2010` and Skovrup :cite:`Skovrup2013`,
using a centred approach for the independent variables enhances the fit quality.
Therefore, all solutions have a base temperature and concentration in the original
works as well as in CoolProp: :math:`x_\text{in} = x - x_\text{base}`
and :math:`T_\text{in} = T - T_\text{base}`, this technique does not affect the calculation
of the derived quantity internal energy since the formula contains temperature differences.
However, integrating :math:`c(x_\text{in},T_\text{in})T_\text{in}^{-1}dT_\text{in}` for the entropy requires some changes due to
the logarithm.

.. warning::
   You must **not** use the base temperature :math:`T_\text{base}`
   as reference temperature for your thermodynamic states. This will lead to an
   error caused by a division by zero during the integration carried out to
   obtain the entropy.

To structure the problem, we introduce a variable :math:`f(j,T)`,
which will be expressed by a third sum. As a first step for simplification, one
has to expand the the binomial :math:`(T-T_{base})^n` to a series. Only
containing :math:`j` and :math:`T`, :math:`f` is independent from :math:`x_\text{in}` and
can be computed outside the loop for enhanced computational efficiency. An
integration of the expanded binomial then yields the final factor :math:`F` to
be multiplied with the other coefficients and the concentration.

.. math::

    \int_{0}^{1} \left( \frac{\partial s}{\partial T} \right)_p dT          &= \int_{0}^{1} \frac{c\left( x_\text{in},T_\text{in} \right)}{T_\text{in}} dT_\text{in} = \sum_{i=0}^n x_\text{in}^i \cdot \sum_{j=0}^m C_{c}[i,j] \cdot F(j,T_\text{in,0},T_\text{in,1}) \\
    F          &= (-1)^j \cdot \ln \left( \frac{T_\text{in,1}}{T_\text{in,0}} \right) \cdot T_{base}^j + \sum_{k=0}^{j-1} \binom{j}{k} \cdot \frac{(-1)^k}{j-k} \cdot \left( T_\text{in,1}^{j-k} - T_\text{in,0}^{j-k} \right) \cdot T_{base}^k

.. _Equations:

Equations
---------

There are only four different equations used to calculate the thermophysical
properties of incompressible fluids in CoolProp:

.. math::

    f(T)  &= \exp \left( \frac{C[0]}{T+C[1]} - C[2] \right) \text{, } \\
    f(T)  &= \exp \left( \log  \left( \sum_{i=0}^l \left( T+C[0] \right)^{-i-1} \right) \cdot C[1] + C[2] \right) \text{, } \\
    f(T,x)&=             \sum_{i=0}^n x^i \cdot \sum_{j=0}^m C[i,j] \cdot T^j \text{ and } \\
    f(T,x)&= \exp \left( \sum_{i=0}^n x^i \cdot \sum_{j=0}^m C[i,j] \cdot T^j \right) \text{. } \\

Only the last two are suitable for mixtures with the input parameter :math:`x`
denoting the fraction of component other than water. Following the works of
Melinder :cite:`Melinder2010` and Skovrup :cite:`Skovrup2013`, the exponents
for the polynomials are arranged in a triangular matrix to avoid overfitting.
These conditions satisfy :math:`0 \leq i \leq n`, :math:`0 \leq j \leq m`
and :math:`i + j \leq \max(n,m)`. It is only for the freezing temperature calculation
that the implemented procedures differ from what is presented in Melinder's
book :cite:`Melinder2010`. Freezing temperature is only a function of concentration
and the dependency on the fluid temperature has been removed. For mixtures,
:math:`m=5` and :math:`n=3` are assigned as default values.
Omitting the composition term with :math:`n=0` yields the pure fluid formulations
for which we selected :math:`l=1` and :math:`m=4`.

The standard polynomials are used for the density, heat capacity and thermal
conductivity functions, while viscosity, vapour pressure and freezing temperature
are exponential functions. For exponential functions of only one variable
(:math:`\mu(T)`, :math:`p_\text{sat}(T)`, :math:`T_\text{freeze}(x)`), we start by fitting the
first equation. If the fit quality is poor, we try the second exponential function.
The exponential polynomial is used as a fall-back function for single variable
fits and it is the only function used for multivariate fits, e.g. :math:`\mu(T,x)`.

If you would like to know more about the fitting procedures, you can have a look
at this `Python notebook <http://nbviewer.ipython.org/github/CoolProp/CoolProp/blob/master/dev/incompressible_liquids/LinearAlgebra.ipynb>`_,
which describes the basics of the multivariate polynomial fits employed in this
software. Non-polynomial functions are fitted using the minimisation routines
accessible through SciPy :cite:`Jones2001`. For the extremely curious, the
Python module `CPIncomp <https://github.com/CoolProp/CoolProp/tree/master/dev/incompressible_liquids/CPIncomp>`_
contains the source code for the fits used in CoolProp as well as the code to
generate the fitting reports. Feel free to browse the code.


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
`SecCool <http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx>`_
:cite:`Skovrup2013` software, which contains data compiled by Morten Juel
Skovrup. It is provided free of charge by his employer `IPU <http://en.ipu.dk>`_.


.. _Pure:

.. csv-table:: All incompressible pure fluids included in CoolProp
   :widths: 10, 35, 13, 14, 14, 14
   :header-rows: 1
   :file: Incompressibles_pure-fluids.csv


There are also a number of water-based mixtures implemented in CoolProp. Most of them
are secondary heat transfer fluids, but there are also aqueous solutions of
ammonia :cite:`Melinder2010`, :download:`MAM<../_static/fluid_properties/Incompressibles_reports/MAM_fitreport.pdf>`,
and lithium bromide :cite:`Patek2006`, :download:`LiBr<../_static/fluid_properties/Incompressibles_reports/LiBr_fitreport.pdf>`,
which can be used to model absorption chillers.


.. _MassMix:

.. csv-table:: All incompressible mass-based binary mixtures included in CoolProp
   :widths: 10, 30, 11, 11, 11, 11, 8, 8
   :header-rows: 1
   :file: Incompressibles_mass-based-fluids.csv

.. .. _MoleMix:

.. .. csv-table:: All incompressible mole-based binary mixtures included in CoolProp
   :widths: 10, 30, 11, 11, 11, 11, 8, 8
   :header-rows: 1
   :file: Incompressibles_mole-based-fluids.csv

.. _VoluMix:

.. csv-table:: All incompressible volume-based binary mixtures included in CoolProp
   :widths: 10, 30, 11, 11, 11, 11, 8, 8
   :header-rows: 1
   :file: Incompressibles_volume-based-fluids.csv


For slurry ice, the concentration :math:`x` refers to the solid content and the
heat capacity includes the heat of fusion. It might be necessary to adjust the
solid content during heat transfer. The implementation is based on the data
available in `SecCool <https://www.ipu.dk/products/seccool>`_,
which was originally recorded at the Danish Technological Institute `(DTI) <http://www.dti.dk/>`_.


.. figure:: /_static/fluid_properties/incompressibles_consistency.jpg
    :align: center
    :alt: Validation figure showing Prandtl numbers for all incompressible fluids

    The figure above shows plots of the Prandtl numbers and associated fluid 
    properties for all incompressible fluids covering the whole temperature range
    for each fluid. The original 
    :download:`PDF version</_static/fluid_properties/incompressibles_consistency.pdf>` 
    is also available for download.


Adding New Fluids
-----------------

To add a fluid to the backend for incompressible fluids, you have to have the tabulated
property data available. Pure fluids are added to the ``PureFluids.py`` and binary 
mixtures, like aqueous mixtures, have to be added to the ``SolutionFLuids.py``. 

The basic state variable for incompressible fluids is temperature, which should be provided
as a one-dimensional numpy array with length :math:`N`. For pure fluids, all properties should match
this temperature array in size since there is a 1-to-1 relation between the temperature points
and the other quantities. 

For binary mixtures, you also need a composition vector to define the data points properly.
This composition vector is also a one-dimensional numpy array of the lenghth :math:`M` and forms the 
second axis for your property space with :math:`N \times M` data points.

Note that all properties have to have the same grid in terms of temperature and composition. 

Once you haver added your fluid data, you can regenrate the JSON files with the fitted 
parameters by running the script located at ``dev/incompressible_liquids/all_incompressibles.py``. 
Your new fluid is now part of the codebase and should be available to all CoolProp functions as
soon as you recompile the sources.


