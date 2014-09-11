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

.. plot::

    import CoolProp
    import matplotlib.pyplot as plt

    HEOS = CoolProp.AbstractState('HEOS','R32&R134a')
    for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
        HEOS.set_mole_fractions([x0, 1 - x0])
        try:
            HEOS.build_phase_envelope("dummy")
        except ValueError as VE:
            print VE
        PE = HEOS.get_phase_envelope_data()
        plt.plot(PE.rhomolar_vap, PE.rhomolar_liq, 'o-')

    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    plt.tight_layout()