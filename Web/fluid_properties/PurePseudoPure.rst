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
   
References
----------

.. bibliography:: ../../CoolPropBibTeXLibrary.bib
   :filter: docname in docnames
   :style: unsrt
   

.. include:: fluidstoc.rst.in
