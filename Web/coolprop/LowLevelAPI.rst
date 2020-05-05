.. _low_level_api:

*******************
Low Level Interface
*******************

.. contents:: :depth: 2

For more advanced use, it can be useful to have access to lower-level internals of the CoolProp code.  For simpler use, you can use the :ref:`high-level interface <high_level_api>`.  The primary reason why this low-level interface is useful is because it is much faster, and actually the high-level interface internally calls the low-level interface.  Furthermore, the low-level-interface exclusively operates using enumerated values (integers) and floating point numbers, and uses no strings.  String comparison and parsing is computationally expensive and the low-level interface allows for a very efficient execution.

At the C++ level, the code is based on the use of an :cpapi:`AbstractState` `abstract base class  <http://en.wikipedia.org/wiki/Abstract_type>`_ which defines a protocol that :ref:`the property backends <backends>` must implement.  In this way, it is very easy to extend CoolProp to connect with another completely unrelated property library, as was done for REFPROP.  As long as the interface to the library can be coerced to work within the AbstractState structure, CoolProp can interface seamlessly with the library.

In order to make most effective use of the low-level interface, you should instantiate one instance of the backend for each fluid (or mixture), and then call methods within the instance.  There is a certain amount of computational overhead in calling the constructor for the backend instance, so in order to minimize it, only call the constructor once, and pass around your class instance.

.. warning::

    While the warning about the computational overhead when generating AbstractState instances is more a recommendation, it is *required* that you allocate as few AbstractState instances as possible when using the :ref:`tabular backends (TTSE & Bicubic) <tabular_interpolation>`.

.. warning::

    In C++, the :cpapi:`AbstractState::factory` function returns a bare pointer to an AbstractState instance, you must be very careful with this instance to make sure it is appropriately destroyed.  It is HIGHLY recommended to wrap the generated instance in a shared pointer as shown in the example.  The shared pointer will take care of automatically calling the destructor of the AbstractState when needed.
    
Similar methodology is used in the other wrappers of the low-level interface to (mostly) generate 1-to-1 wrappers of the low-level functions to the target language.  Refer to the examples for each language to see how to call the low-level interface, generate an AbstractState instance, etc.

Introduction
------------

To begin with, an illustrative example of using the low-level interface in C++ is shown here:

.. literalinclude:: snippets/AbstractState1.cxx
   :language: c++

which yields the output:

.. literalinclude:: snippets/AbstractState1.cxx.output

This example demonstrates the most common application of the low-level interface.  This example could also be carried out using the high-level interface with a call like::

    PropsSI("T", "P", 101325, "Q", 1, "Water")
    
Generating Input Pairs
----------------------

A listing of the input pairs for the :cpapi:`AbstractState::update` function can be found in the source documentation at :cpapi:`CoolProp::input_pairs`.  If you know the two input variables and their values, but not their order, you can use the function :cpapi:`CoolProp::generate_update_pair` to generate the input pair.  

.. warning::

    The syntax for this function is slightly different in python since python can do multiple return arguments and C++ cannot.  
    
A simple example of this would be 

.. ipython::

    In [0]: import CoolProp
    
    In [0]: CoolProp.CoolProp.generate_update_pair(CoolProp.iT, 300, CoolProp.iDmolar, 1e-6)
    
    In [0]: CoolProp.DmolarT_INPUTS

Keyed output versus accessor functions
---------------------------------------

The simple output functions like :cpapi:`AbstractState::rhomolar` that are mapped to keys in :cpapi:`CoolProp::parameters` can be either obtained using the accessor function or by calling :cpapi:`AbstractState::keyed_output`.  The advantage of the ``keyed_output`` function is that you could in principle iterate over several keys, rather than having to hard-code calls to several accessor functions.  For instance:

.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState("HEOS", "Water")
    
    In [0]: HEOS.update(CoolProp.DmolarT_INPUTS, 1e-6, 300)
    
    In [0]: HEOS.p()
    
    In [0]: [HEOS.keyed_output(k) for k in [CoolProp.iP, CoolProp.iHmass, CoolProp.iHmolar]]
    
Things only in the low-level interface
--------------------------------------
    
You might reasonably ask at this point why we would want to use the low-level interface as opposed to the "simple" high-level interface.  In the first example, if you wanted to calculate all these output parameters using the high-level interface, it would require several calls to the pressure-quality flash routine, which is extremely slow as it requires a complex iteration to find the phases that are in equilibrium.  Furthermore, there is a lot of functionality that is only accessible through the low-level interface.  Here are a few examples of things that can be done in the low-level interface that cannot be done in the high-level interface:
   
.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState("HEOS", "Water")
    
    # Do a flash call that is a very low density state point, definitely vapor
    In [0]: %timeit HEOS.update(CoolProp.DmolarT_INPUTS, 1e-6, 300)
    
    # Specify the phase - for some inputs (especially density-temperature), this will result in a 
    # more direct evaluation of the equation of state without checking the saturation boundary
    In [0]: HEOS.specify_phase(CoolProp.iphase_gas)
    
    # We try it again - a bit faster
    In [0]: %timeit HEOS.update(CoolProp.DmolarT_INPUTS, 1e-6, 300)
    
    # Reset the specification of phase
    In [0]: HEOS.specify_phase(CoolProp.iphase_not_imposed)
    
    # A mixture of methane and ethane
    In [0]: HEOS = CoolProp.AbstractState("HEOS", "Methane&Ethane")
    
    # Set the mole fractions of the mixture
    In [0]: HEOS.set_mole_fractions([0.2,0.8])
    
    # Do the dewpoint calculation
    In [0]: HEOS.update(CoolProp.PQ_INPUTS, 101325, 1)
    
    # Liquid phase molar density
    In [0]: HEOS.saturated_liquid_keyed_output(CoolProp.iDmolar)
    
    # Vapor phase molar density
    In [0]: HEOS.saturated_vapor_keyed_output(CoolProp.iDmolar)
    
    # Liquid phase mole fractions
    In [0]: HEOS.mole_fractions_liquid()
    
    # Vapor phase mole fractions - 
    # Should be the bulk composition back since we are doing a dewpoint calculation
    In [0]: HEOS.mole_fractions_vapor()

Imposing the Phase (Optional)
-----------------------------

Each call to ``update()`` requires the phase to be determined based on the provided input pair, and may require a non-trivial flash calculation to determine if the state point is in the single-phase or two-phase region. For computational efficiency, CoolProp provide a facility to manually impose the phase before updating the state point. As demonstrated above, specifying the phase of a state point can result in a speed-up of ~10%, depending on the input pair provided and the phase region specified. If unspecified, CoolProp will determine the phase.  

Depending on the input pair, there may or may not be a speed benefit to imposing a phase, like the one demonstrated above for **DmolarT_INPUTS**.  If, for a given input pair, there is little or no speed benefit, or if the flash calculation has to be performed anyway to determine saturation values, the imposed phase may simply be ignored internally.

To specify the phase to be used, the AbstractState provides the following function:

.. ipython::

    # Imposing the single-phase gas region
    In [0]: HEOS.specify_phase(CoolProp.iphase_gas)

The specified phase parameter can be one of the CoolProp constants in the table below. 

+---------------------------------+----------------------------------------------------+
| Phase Constant                  | Phase Region                                       |
+=================================+====================================================+
| iphase_liquid                   | p < pcrit & T < Tcrit ; above saturation           |
+---------------------------------+----------------------------------------------------+
| iphase_gas                      | p < pcrit & T < Tcrit ; below saturation           |
+---------------------------------+----------------------------------------------------+
| iphase_twophase                 | p < pcrit & T < Tcrit ; mixed liquid/gas           |
+---------------------------------+----------------------------------------------------+
| iphase_supercritical_liquid     | p > pcrit & T < Tcrit                              |
+---------------------------------+----------------------------------------------------+
| iphase_supercritical_gas        | p < pcrit & T > Tcrit                              |
+---------------------------------+----------------------------------------------------+
| iphase_supercritical            | p > pcrit & T > Tcrit                              |
+---------------------------------+----------------------------------------------------+
| iphase_not_imposed              | (Default) CoolProp to determine phase              |
+---------------------------------+----------------------------------------------------+

To remove an imposed phase and allow CoolProp to determine the phase again, one of the following calls can be made.
    
.. ipython::

    # Reset the specification of phase
    In [0]: HEOS.specify_phase(CoolProp.iphase_not_imposed)

    # Or, more simply
    In [0]: HEOS.unspecify_phase()

Phase specification in the low-level interface should be performed right before the call to ``update()``.

.. warning::

   When specifying an imposed phase, it is absolutely **critical** that the input pair actually lie within the imposed phase region.  If an incorrect phase is imposed for the given input pair, ``update()`` may throw unexpected errors or incorrect results may possibly be returned from the property functions.  If the state point phase is not absolutely known, it is best to let CoolProp determine the phase at the cost of some computational efficiency.

.. _partial_derivatives_low_level:

Partial Derivatives
-------------------

It is possible to get the partial derivatives in a very computationally efficient manner using the low-level interface, using something like (python here):

For more information, see the docs: :cpapi:`CoolProp::AbstractState::first_partial_deriv` and :cpapi:`CoolProp::AbstractState::second_partial_deriv`

.. ipython::

    In [1]: import CoolProp

    In [2]: HEOS = CoolProp.AbstractState("HEOS", "Water")

    In [3]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300)

    In [4]: HEOS.cpmass()

    In [4]: HEOS.first_partial_deriv(CoolProp.iHmass, CoolProp.iT, CoolProp.iP)
    
    In [4]: %timeit HEOS.first_partial_deriv(CoolProp.iHmass, CoolProp.iT, CoolProp.iP)
    
    # See how much faster this is?
    In [4]: %timeit CoolProp.CoolProp.PropsSI('d(Hmass)/d(T)|P', 'P', 101325, 'T', 300, 'Water')
    
    In [4]: HEOS.first_partial_deriv(CoolProp.iSmass, CoolProp.iT, CoolProp.iDmass)

    # In the same way you can do second partial derivatives
    # This is the second mixed partial derivative of entropy with respect to density and temperature
    In [4]: HEOS.second_partial_deriv(CoolProp.iSmass, CoolProp.iT, CoolProp.iDmass, CoolProp.iDmass, CoolProp.iT)
    # This is the second partial derivative of entropy with respect to density at constant temperature
    In [4]: HEOS.second_partial_deriv(CoolProp.iSmass,CoolProp.iDmass,CoolProp.iT,CoolProp.iDmass,CoolProp.iT)
    
Two-Phase and Saturation Derivatives
------------------------------------

The two-phase derivatives of Thorade :cite:`Thorade-EES-2013` are implemented in the :cpapi:`CoolProp::AbstractState::first_two_phase_deriv` function, and derivatives along the saturation curve in the functions :cpapi:`CoolProp::AbstractState::first_saturation_deriv` and :cpapi:`CoolProp::AbstractState::second_saturation_deriv`.  Here are some examples of using these functions:
    
.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState("HEOS", "Water")
    
    In [0]: HEOS.update(CoolProp.QT_INPUTS, 0, 300)

    # First saturation derivative calculated analytically
    In [0]: HEOS.first_saturation_deriv(CoolProp.iP, CoolProp.iT)
    
    In [0]: HEOS.update(CoolProp.QT_INPUTS, 0, 300 + 0.001); p2 = HEOS.p()
    
    In [0]: HEOS.update(CoolProp.QT_INPUTS, 0, 300 - 0.001); p1 = HEOS.p()
    
    # First saturation derivative calculated numerically
    In [0]: (p2-p1)/(2*0.001)
     
    In [0]: HEOS.update(CoolProp.QT_INPUTS, 0.1, 300)
    
    # The d(Dmass)/d(Hmass)|P two-phase derivative
    In [0]: HEOS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP)
    
    # The d(Dmass)/d(Hmass)|P two-phase derivative using splines
    In [0]: HEOS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3)
    
An example of plotting these derivatives is here:

.. plot::

    import numpy as np
    import CoolProp
    import matplotlib.pyplot as plt

    AS = CoolProp.AbstractState('HEOS','Water')

    # Saturated liquid
    AS.update(CoolProp.PQ_INPUTS, 101325, 0)
    Ts = AS.T()
    h_fg = AS.saturated_vapor_keyed_output(CoolProp.iHmass) - AS.saturated_liquid_keyed_output(CoolProp.iHmass)
    cl = AS.cpmass()

    # Subcooled liquid
    x, y = [], []
    for T in np.linspace(Ts - 30, Ts - 0.1, 1000):
        AS.update(CoolProp.PT_INPUTS, 101325, T)
        x.append(-cl*(Ts-T)/h_fg)
        y.append(AS.first_partial_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
    plt.plot(x, y, label = 'Subcooled')

    # Two-phase derivatives (normal and splined)
    x, y1 = [], []
    for Q in np.linspace(0, 0.3, 1000):
        AS.update(CoolProp.PQ_INPUTS, 101325, Q)
        x.append(AS.Q())
        y1.append(AS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3))
    plt.plot(x, y1, label = 'Two-phase (splined)')

    # Two-phase derivatives (normal and splined)
    x, y1 = [], []
    for Q in np.linspace(0.0, 0.6, 1000):
        AS.update(CoolProp.PQ_INPUTS, 101325, Q)
        x.append(AS.Q())
        y1.append(AS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
    plt.plot(x, y1, label = 'Two-phase')

    plt.title(r'$d\rho/dh|p$')
    plt.xlabel('vapor quality (-)')
    plt.ylabel(r'$d\rho/dh|p$')
    plt.ylim(-0.005, 0.005)
    plt.legend(loc='best')
    
Reference States
----------------

To begin with, you should read :ref:`the high-level docs about the reference state <high_level_set_reference_state>`.  Those docs are also applicable to the low-level interface.  

.. warning::

    As with the high-level interface, calling ``set_reference_stateS`` (or ``set_reference_state`` in python) should be called right at the beginning of your code, and not changed again later on.
    
    Importantly, once an ``AbstractState``-derived instance has been generated from the factory function, it **DOES NOT** pick up the change in the reference state.  This is intentional, but you should watch out for this behavior.
    
Here is an example showing how to change the reference state and demonstrating the potential issues

.. ipython::

    In [1]: import CoolProp as CP
    
    # This one doesn't see the change in reference state
    In [1]: AS1 = CoolProp.AbstractState('HEOS','n-Propane'); 
    
    In [1]: AS1.update(CoolProp.QT_INPUTS, 0, 233.15); 
    
    In [1]: CoolProp.CoolProp.set_reference_state('n-Propane','ASHRAE')
    
    # This one gets the update in the reference state
    In [1]: AS2 = CoolProp.AbstractState('HEOS','n-Propane'); 
    
    In [1]: AS2.update(CoolProp.QT_INPUTS, 0, 233.15); 
    
    # Note how the AS1 has its default value (change in reference state is not seen)
    # and AS2 does see the new reference state
    In [1]: print('{0}, {1}'.format(AS1.hmass(), AS2.hmass()))
    
    # Back to the original value
    In [1]: CoolProp.CoolProp.set_reference_state('n-Propane','DEF')
    
Low-level interface using REFPROP
---------------------------------

If you have the `REFPROP library <http://www.nist.gov/srd/nist23.cfm>`_ installed, you can call REFPROP in the same way that you call CoolProp, but with ``REFPROP`` as the backend instead of ``HEOS``. For instance, as in python:

.. ipython::

    In [0]: import CoolProp

    In [0]: REFPROP = CoolProp.AbstractState("REFPROP", "WATER")
    
    In [0]: REFPROP.update(CoolProp.DmolarT_INPUTS, 1e-6, 300)
    
    In [0]: REFPROP.p(), REFPROP.hmass(), REFPROP.hmolar()
    
    In [0]: [REFPROP.keyed_output(k) for k in [CoolProp.iP, CoolProp.iHmass, CoolProp.iHmolar]]

.. _low_level_high_level:
    
Access from High-Level Interface
--------------------------------

For languages that cannot directly instantiate C++ classes or their wrappers but can call a DLL (Excel, FORTRAN, Julia, etc. ) an interface layer has been developed.  These functions can be found in :cpapi:`CoolPropLib.h`, and all these function names start with ``AbstractState``.  A somewhat limited subset of the functionality has been implemented, if more functionality is desired, please open an issue at https://github.com/CoolProp/CoolProp/issues.  Essentially, this interface generates :cpapi:`AbstractState` pointers that are managed internally, and the interface allows you to call the methods of the low-level instances.

Here is an example of the shared library usage with Julia wrapper::

    julia> import CoolProp

    julia> PT_INPUTS = CoolProp.get_input_pair_index("PT_INPUTS")
    7

    julia> cpmass = CoolProp.get_param_index("C")
    34

    julia> handle = CoolProp.AbstractState_factory("HEOS", "Water")
    0

    julia> CoolProp.AbstractState_update(handle,PT_INPUTS,101325, 300)

    julia> CoolProp.AbstractState_keyed_output(handle,cpmass)
    4180.635776569655

    julia> CoolProp.AbstractState_free(handle)

    julia> handle = CoolProp.AbstractState_factory("HEOS", "Water&Ethanol")
    1

    julia> PQ_INPUTS = CoolProp.get_input_pair_index("PQ_INPUTS")
    2

    julia> T = CoolProp.get_param_index("T")
    18

    julia> CoolProp.AbstractState_set_fractions(handle, [0.4, 0.6])

    julia> CoolProp.AbstractState_update(handle,PQ_INPUTS,101325, 0)

    julia> CoolProp.AbstractState_keyed_output(handle,T)
    352.3522142890429

    julia> CoolProp.AbstractState_free(handle)
    
The call to ``AbstractState_free`` is not strictly needed as all managed AbstractState instances will auto-deallocate.  But if you are generating thousands or millions of AbstractState instances in this fashion, you might want to tidy up periodically. 

Here is a further example in C++:

.. literalinclude:: snippets/HighLevelLowLevel.cxx
   :language: c++

which yields the output:

.. literalinclude:: snippets/HighLevelLowLevel.cxx.output

Here is a further example in C++ that shows how to obtain many output variables at the same time, either 5 common outputs, or up to 5 user-selected vectors:

.. literalinclude:: snippets/HighLevelLowLevelMulti.cxx
   :language: c++

which yields the output:

.. literalinclude:: snippets/HighLevelLowLevelMulti.cxx.output
