.. _low_level_api:

*******************
Low Level Interface
*******************

For more advanced use, it can be useful to have access to lower-level internals of the CoolProp code.  Also, this interface is strongly recommended for mixtures, where the high-level api has some serious limitations.  For simpler use, you can use the :ref:`high-level interface <high_level_api>`.

At the C++ level, the code is based on the use of an :cpapi:`AbstractState` `abstract base class  <http://en.wikipedia.org/wiki/Abstract_type>`_ which defines a protocol that :ref:`the property backends <backends>` must implement.  In this way, it is very easy to extend CoolProp to connect with another completely unrelated property library, as was done for REFPROP.

To begin with, an illustrative example of using the low-level interface is shown here:

.. _partial_derivatives_low_level:
    
Partial Derivatives
-------------------

TO BE CONTINUED


.. literalinclude:: snippets/AbstractState1.cxx
   :language: c++

which yields the output:

.. literalinclude:: snippets/AbstractState1.cxx.output

Documentation on calling the low-level interface from the wrappers will follow.
