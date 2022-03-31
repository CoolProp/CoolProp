import CoolProp
CPmod = CoolProp.CoolProp


def get_offset_NBP(name):
    # CPmod.set_debug_level(10)
    CPmod.set_reference_state(name, "RESET")
    HEOS = CoolProp.AbstractState('HEOS', name)
    HEOS.update(CoolProp.PQ_INPUTS, 101325, 0)

    gas_constant = HEOS.gas_constant() / HEOS.molar_mass()

    delta_a1 = HEOS.smass() / (gas_constant)
    delta_a2 = -HEOS.hmass() / (gas_constant * HEOS.keyed_output(CoolProp.iT_reducing))
    return delta_a1, delta_a2


if __name__ == '__main__':
    name = 'PENTANE'
    import json
    a1, a2 = get_offset_NBP(name)
    print(json.dumps({
      "a1": a1,
      "a2": a2,
      "type": "IdealGasHelmholtzLead"
    }, indent=2))
