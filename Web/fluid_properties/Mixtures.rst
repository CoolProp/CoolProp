.. _mixtures:

********
Mixtures
********

.. contents:: :depth: 2

Theoretical description
-----------------------
The mixture modeling used in CoolProp is based on the work of Kunz et al. :cite:`Kunz-BOOK-2007,Kunz-JCED-2012` and Lemmon :cite:`Lemmon-JPCRD-2000,Lemmon-JPCRD-2004,Lemmon-IJT-1999`

A mixture is composed of a number of components, and for each pair of components, it is necessary to have information for the excess Helmholtz energy term as well as the reducing function.  See below for what binary pairs are included in CoolProp.

The numerical methods required for mixtures are far more complicated than those for pure fluids, so the number of flash routines that are currently available are relatively small compared to pure fluids.

The only types of inputs that are allowed for mixtures right now are

- Pressure/quality
- Temperature/quality
- Temperature/pressure

.. Used in Python script later on
.. role:: raw-html(raw)
   :format: html

Estimating binary interaction parameters
----------------------------------------

If you have a mixture that you would like to include in your analysis, but for which no interaction parameters are currently available, there are two estimation schemes available.  These estimation schemes should be used with extreme caution.  If you use these schemes, you should definitely check that you are getting reasonable output values.  In general, these estimation schemes should be used for fluids with similar properties.

The two schemes available are

* ``linear`` - :math:`T_r` and :math:`v_r` are a linear function of molar composition between the two pure fluid values
* ``Lorentz-Berthelot`` - all interaction parameters are 1.0

Note that if  ``apply_simple_mixing_rule`` is called first in a programme before any property calculation function, the interaction parameter library is initialized with just the binary mixture that was specified, and other binary mixtures will raise an error. To augment the built-in database, call a property calculation function before calling ``apply_simple_mixing_rule``.

Here is a sample of using this in python:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    # Create a simple linear mixing rule for the Helium and Xenon pair
    In [1]: CP.apply_simple_mixing_rule('Helium', 'Xenon', 'linear')

    # Now evaluate the 50/50 mixture density at 300K and atmospheric pressure
    In [1]: CP.PropsSI('Dmass','T',300,'P',101325,'Helium[0.5]&Xenon[0.5]')
    
.. warning::

    Use with caution!! For other mixtures this can give you entirely(!) wrong predictions
    

Using your own interaction parameters
-------------------------------------

If you have your own interaction parameters that you would like to use, you can set them using the ``set_mixture_binary_pair_data`` function.  (You can also retrieve them using the ``get_mixture_binary_pair_data`` function).  You must do this before you call any other functions that use them.  Also, you must set the boolean configuration variable ``OVERWRITE_BINARY_INTERACTION`` to ``true``.  This variable is ``false`` by default and will dis-allow overwrites of parameters that already exist.  Some sample code is below.

.. ipython::
    :okexcept:

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CAS_He = CP.get_fluid_param_string('Helium','CAS')

    In [1]: CAS_Xe = CP.get_fluid_param_string('Xenon','CAS')

    @suppress
    In [1]: CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)

    # This adds a dummy entry in the library of interaction parameters if the mixture is not already there
    In [1]: CP.apply_simple_mixing_rule(CAS_He, CAS_Xe, 'linear')

    # If the binary interaction pair is already there, set the configuration flag 
    # to allow binary interaction parameters to be over-written
    In [1]: CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)
    
    # This is before setting the binary interaction parameters
    In [1]: CP.PropsSI('Dmass','T',300,'P',101325,'Helium[0.5]&Xenon[0.5]')

    In [1]: CP.set_mixture_binary_pair_data(CAS_He, CAS_Xe, 'betaT', 1.0)
    
    In [1]: CP.set_mixture_binary_pair_data(CAS_He, CAS_Xe, 'gammaT', 1.5)
    
    In [1]: CP.set_mixture_binary_pair_data(CAS_He, CAS_Xe, 'betaV', 1.0)
    
    In [1]: CP.set_mixture_binary_pair_data(CAS_He, CAS_Xe, 'gammaV', 1.5)

    # This is after setting the interaction parameters
    In [1]: CP.PropsSI('Dmass','T',300,'P',101325,'Helium[0.5]&Xenon[0.5]')
    
Once you have constructed an instance of an AbstractState using the low-level interface, you can set the interaction parameters *for only that instance* by calling the ``set_binary_interaction_double`` and ``get_binary_interaction_double`` functions.  This will have no effect on other instances, including the high-level calls, as shown below.

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    # This adds a dummy entry in the library of interaction parameters if the mixture is not already there
    In [1]: CP.apply_simple_mixing_rule(CAS_He, CAS_Xe, 'linear')
    
    # First use the high-level call to PropsSI to get the mass density of a 50/50 mixture:
    In [1]: CP.PropsSI('Dmass','T',300,'P',101325,'Helium[0.5]&Xenon[0.5]')
    
    # Now, use the low-level interface to modify the binary interaction parameters
    In [1]: AS = CP.AbstractState("HEOS","Helium&Xenon")
    
    In [1]: AS.set_binary_interaction_double(0, 1, 'betaT', 0.987)
    
    In [1]: AS.get_binary_interaction_double(0, 1, 'betaT')
    
    In [1]: AS.set_binary_interaction_double(0, 1, 'gammaT', 1.5)
    
    In [1]: AS.set_binary_interaction_double(0, 1, 'betaV', 1.0)
    
    In [1]: AS.set_binary_interaction_double(0, 1, 'gammaV', 1.5)

    In [1]: AS.set_mole_fractions([0.5,0.5])
    
    In [1]: AS.update(CP.PT_INPUTS, 101325, 300)
    
    In [1]: AS.rhomass()
    
    # Here you can see that this call to the high-level interface is untouched; 
    # giving the same ligh-level result from PropsSI() as above.
    In [1]: CP.PropsSI('Dmass','T',300,'P',101325,'Helium[0.5]&Xenon[0.5]')

And now, reset the OVERWRITE_BINARY_INTERACTION configuration variable to dis-allow overwrites.

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, False)

Phase Envelope
--------------

You can download the script that generated the following figure here: :download:`(link to script)<methane-ethane.py>`, right-click the link and then save as... or the equivalent in your browser.  You can also download this figure :download:`as a PDF<methane-ethane.pdf>`. 

.. image:: methane-ethane.png

Reducing Parameters
-------------------

From Lemmon :cite:`Lemmon-JPCRD-2000` for the properties of Dry Air, and also from Lemmon :cite:`Lemmon-JPCRD-2004` for the properties of R404A, R410A, etc.

.. math::

    \rho_r(\bar x) = \left[ \sum_{i=1}^m\frac{x_i}{\rho_{c_i}}+\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}x_ix_j\zeta_{ij}\right]^{-1}

.. math::

    T_r(\bar x) = \sum_{i=1}^mx_iT_{c_i}+\sum_{i=1}^{m-1}\sum_{j=i+1}^mx_ix_j\xi_{ij}

From the GERG 2008 formulation :cite:`Kunz-JCED-2012`

.. math::

    T_r(\bar x) = \sum_{i=1}^{N}x_i^2T_{c,i} + \sum_{i=1}^{N-1}\sum_{j=i+1}^{N}2x_ix_j\beta_{T,ij}\gamma_{T,ij}\frac{x_i+x_j}{\beta_{T,ij}^2x_i+x_j}(T_{c,i}T_{c,j})^{0.5}
    
.. math::

    \frac{1}{\rho_r(\bar x)}=v_r(\bar x) = \sum_{i=1}^{N}x_i^2\frac{1}{\rho_{c,i}} + \sum_{i=1}^{N-1}\sum_{j=i+1}^N2x_ix_j\beta_{v,ij}\gamma_{v,ij}\frac{x_i+x_j}{\beta^2_{v,ij}x_i+x_j}\frac{1}{8}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}
    
Excess Helmholtz Energy Terms
-----------------------------
From Lemmon :cite:`Lemmon-JPCRD-2004` for the properties of R404A, R410A, etc.

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \sum_{i=1}^{m-1} \sum_{j=i+1}^{m} \left [ x_ix_jF_{ij} \sum_{k}N_k\delta^{d_k}\tau^{t_k}\exp(-\delta^{l_k})\right]
    
where the terms :math:`N_k,d_k,t_k,l_k` correspond to the pair given by the indices :math:`i,j`

From Lemmon :cite:`Lemmon-JPCRD-2000` for the properties of Dry Air

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \left \lbrace \sum_{i=1}^{2} \sum_{j=i+1}^{3} x_ix_jF_{ij}\right\rbrace \left[-0.00195245\delta^2\tau^{-1.4}+0.00871334\delta^2\tau^{1.5} \right]


From Kunz and Wagner :cite:`Kunz-JCED-2012` for GERG 2008 formulation

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \sum_{i=1}^{N-1} \sum_{j=i+1}^{N} x_ix_jF_{ij}\alpha_{ij}^r(\delta,\tau)
    
where

.. math::

    \alpha_{ij}^r(\delta,\tau) = \sum_{k=1}^{K_{pol,ij}}\eta_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}+\sum_{k=K_{pol,ij}+1}^{K_{pol,ij}+K_{Exp,ij}}\eta_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}\exp[-\eta_{ij,k}(\delta-\varepsilon_{ij,k})^2-\beta_{ij,k}(\delta-\gamma_{ij,k})]
    
and is for the particular binary pair given by the indices :math:`i,j`.  This term is similar in form to other Helmholtz energy terms for pure fluids though the derivatives are slightly special.

Appendix
--------
To convert from the form from Lemmon for HFC and Air to that of GERG 2008, the following steps are required:

.. math::

    x_0T_{c0}+(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = x_0^2T_{c0}+(1-x_0)^2T_{c1} + 2x_0(1-x_0)\beta\gamma_T\frac{x_0+(1-x_0)}{\beta x_0 + (1-x_0)}\sqrt{T_{c0}T_{c1}}
    
set :math:`\beta=1`, solve for :math:`\gamma`.  Equate the terms

.. math::

    x_0T_{c0}+(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = x_0^2T_{c0}+(1-x_0)^2T_{c1} + 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Move to LHS

.. math::

    [x_0-x_0^2]T_{c0}+[(1-x_0)-(1-x_0)^2]T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}

Factor

.. math::

    x_0(1-x_0)T_{c0}+(1-x_0)[1-(1-x_0)]T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Expand

.. math::

    x_0(1-x_0)T_{c0}+x_0(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Cancel factors of :math:`x_0(1-x_0)`

.. math::

    T_{c0}+T_{c1}+\xi_{01} = 2\gamma_T\sqrt{T_{c0}T_{c1}}
    
Answer:

.. math::

    \boxed{\gamma_T = \dfrac{T_{c0}+T_{c1}+\xi_{01}}{2\sqrt{T_{c0}T_{c1}}}}
    
Same idea for the volume

.. math::

    \boxed{\gamma_v = \dfrac{v_{c0}+v_{c1}+\zeta_{01}}{\frac{1}{4}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}}}


Binary pairs
------------

.. note::
   Please hover the mouse pointer over the coefficients to get the full accuracy
   for the listed coefficients. You can also get more information on references
   that are not in bibliography.

.. csv-table:: All binary pairs included in CoolProp
   :header-rows: 1
   :file: Mixtures.csv
