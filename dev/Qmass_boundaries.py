import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState


def test_phase_envelope_R410A():
    print(
        "Testing mass quality calculations with R410A.MIX at 278.15 K and 0.5 quality"
    )
    print("=" * 49)

    try:
        # Create abstract state for R410A.MIX
        AS: AbstractState = AbstractState("REFPROP", "R410A.MIX")

        AS.update(CP.QT_INPUTS, 0.5, 278.15)
        print(f"  Viscosity: {AS.viscosity():.12f} Pa.s")
        print(f"  Pressure: {AS.p():.12f} Pa")

        AS.update(CP.QT_INPUTS, 0.0, 278.15)
        print(f"  Viscosity from PQ (9e5, 0.0): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.QT_INPUTS, 1.0, 278.15)
        print(f"  Viscosity from PQ (9e5, 1.0): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.QmassT_INPUTS, 0.0, 278.15)
        print(f"  Viscosity from Qmass (0.0, 278.15): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.QmassT_INPUTS, 1.0, 278.15)
        print(f"  Viscosity from Qmass (1.0, 278.15): {AS.viscosity():.12f} Pa.s")

        # -----------------------------------------

        AS.update(CP.PQ_INPUTS, 9e5, 0.5)
        print(f"  Viscosity: {AS.viscosity():.12f} Pa.s")

        AS.update(CP.PQ_INPUTS, 9e5, 0.0)
        print(f"  Viscosity from PQ (9e5, 0.0): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.PQ_INPUTS, 9e5, 1.0)
        print(f"  Viscosity from PQ (9e5, 1.0): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.PQmass_INPUTS, 9e5, 0.0)
        print(f"  Viscosity from PQmass (9e5, 0.0): {AS.viscosity():.12f} Pa.s")

        AS.update(CP.PQmass_INPUTS, 9e5, 1.0)
        print(f"  Viscosity from PQmass (9e5, 1.0): {AS.viscosity():.12f} Pa.s")

        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    test_phase_envelope_R410A()
