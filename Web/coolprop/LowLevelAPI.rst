*************
Low Level API
*************

For more advanced use, it can be useful to have access to lower-level internals of the CoolProp code.  Also, this interface is strongly recommended for mixtures, where the high-level api has some serious limitations.  For simpler use, you can use the :ref:`high-level interface <high_level_api>`.

At the C++ level, the code is based on the use of an :cpapi:`AbstractState` `abstract base class  <http://en.wikipedia.org/wiki/Abstract_type>`_ which defines a protocol that :ref:`the property backends <backends>` must implement.  In this way, it is very easy to extend CoolProp to connect with another completely unrelated property library, as was done for REFPROP.

Code
----
.. literalinclude:: snippets/AbstractState1.cxx
   :language: c++

Output
------
.. literalinclude:: snippets/AbstractState1.cxx.output