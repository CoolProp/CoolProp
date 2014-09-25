.. _mixtures:

********
Mixtures
********

Mixtures docs will go here when they are written.

The treatment of mixtures in CoolProp as of v5 is quite rudimentary, though it will be improved in the very near future.

The only types of inputs that are allowed for mixtures right now are

- Pressure/quality
- Temperature/quality
- Temperature/pressure

.. csv-table:: All binary pairs included in CoolProp
   :header-rows: 1
   :file: mixture_binary_pairs.csv 

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

.. bibliography:: CoolPropBibTeXLibrary.bib
   :filter: docname in docnames
   :style: unsrt