.. _low_level_api:

*******************
Low Level Interface
*******************

.. contents:: :depth: 2

For more advanced use, it can be useful to have access to lower-level internals of the CoolProp code.  For simpler use, you can use the :ref:`high-level interface <high_level_api>`.  The primary reason why this low-level interface is useful is because it is much faster, and actually the high-level interface internally calls the low-level interface.  Furthermore, the low-level-interface exclusively operates using enumerated values (integers) and floating point numbers, and uses no strings.  String comparison and parsing is computationally expensive and low-level interface allows for a very efficient execution.

At the C++ level, the code is based on the use of an :cpapi:`AbstractState` `abstract base class  <http://en.wikipedia.org/wiki/Abstract_type>`_ which defines a protocol that :ref:`the property backends <backends>` must implement.  In this way, it is very easy to extend CoolProp to connect with another completely unrelated property library, as was done for REFPROP.  As long as the interface to the library can be coerced to work within the AbstractState structure, CoolProp can interface seamlessly with the library.

In order to make most effective use of the low-level interface, you should instantiate one instance of the backend for each fluid (or mixture), and then call methods within the instance.  There is a certain amount of computational overhead in calling the constructor for the backend instance, so in order to minimize it, only call the constructor once, and pass around your class instance.

.. warning::

    In C++, the :cpapi:`AbstractState::factory` function returns a bare pointer to an AbstractState instance, you must be very careful with this instance to make sure it is appropriately destroyed.  It is HIGHLY recommended to wrap the generated instance in a shared pointer as shown in the example.  The shared pointer will take care of automatically calling the destructor of the AbstractState when needed.

To begin with, an illustrative example of using the low-level interface in C++ is shown here:

.. literalinclude:: snippets/AbstractState1.cxx
   :language: c++

which yields the output:

.. literalinclude:: snippets/AbstractState1.cxx.output

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

    julia>

Alternatively, the :cpapi:`AbstractState::keyed_output` function can be called with the appropriate key from :cpapi:`CoolProp::parameters`.  There should be essentially no difference in speed between these two methods.

A list of possible input pairs can be found directly in the source documentation at :cpapi:`CoolProp::input_pairs`.

Similar methodology is used in the other wrappers of the low-level interface to (mostly) generate 1-to-1 wrappers of the low-level functions to the target language.  Refer to the examples for each language to see how to call the low-level interface, generate an AbstractState instance, etc.

.. _partial_derivatives_low_level:

Partial Derivatives
-------------------

It is possible to get the partial derivatives in a very computationally efficient manner using the low-level interface, using something like (python here):

.. ipython::

    In [1]: import CoolProp

    In [2]: HEOS = CoolProp.AbstractState("HEOS", "Water")

    In [3]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300)

    In [4]: HEOS.cpmass()

    In [4]: HEOS.first_partial_deriv(CoolProp.iHmass, CoolProp.iT, CoolProp.iP)
    
Reference States
----------------

To begin with, you should read the high-level docs about the reference state: :ref:`high_level_set_reference_state`.  Those docs are also applicable to the low-level interface.  

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
    
    # Note how the first class has its default value (change in reference state is not seen)
    # and AS2 does see the new reference state
    In [1]: print AS1.hmass(), AS2.hmass()
    
    # Back to the original value
    In [1]: CoolProp.CoolProp.set_reference_state('n-Propane','DEF')
    