.. _mixtures:

.. contents:: :depth: 2

********
Mixtures
********

Theoretical description
-----------------------
The mixture modeling used in CoolProp is based on the work of Kunz et al. :cite:`Kunz-BOOK-2007,Kunz-JCED-2012` and Lemmon :cite:`Lemmon-JPCRD-2000,Lemmon-JPCRD-2004,Lemmon-IJT-1999`

A mixture is composed of a number of components, and for each pair of components, it is necessary to have information for the excess Helmholtz energy term as well as the reducing function.  See below for what binary pairs are included in CoolProp.

The numerical methods required for mixtures are far more complicated than those for pure fluids, so the number of flash routines that are currently available are relatively small compared to pure fluids.

The only types of inputs that are allowed for mixtures right now are

- Pressure/quality
- Temperature/quality
- Temperature/pressure

Binary pairs
------------
.. csv-table:: All binary pairs included in CoolProp
   :header-rows: 1
   :file: mixture_binary_pairs.csv 

Phase Envelope
--------------
.. plot::

    import CoolProp
    import matplotlib.pyplot as plt

    HEOS = CoolProp.AbstractState('HEOS','Methane&Ethane')
    for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
        HEOS.set_mole_fractions([x0, 1 - x0])
        try:
            HEOS.build_phase_envelope("dummy")
        except ValueError as VE:
            print(VE)
        PE = HEOS.get_phase_envelope_data()
        plt.plot(PE.T, PE.p, 'o-')

    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    plt.tight_layout()
    
References
----------

.. bibliography:: ../../CoolPropBibTeXLibrary.bib
   :filter: docname in docnames
   :style: unsrt