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

        AS.update(CP.QT_INPUTS, 0.5, 278.15)  # Saturated vapor (q=0.5 in molar basis)
        print(f"  Saturation pressure (molar basis): {AS.p():.12f} Pa")
        print(f"  Mole density: {AS.rhomolar() / 1e3:.12f} kmol/m3")
        print(f"  Mass density: {AS.rhomass():.12f} kg/m3")
        print(f"  Molar quality: {AS.Q():.12f} -")
        print(f"  Mass quality: {AS.Qmass():.12f} -")
        mf = AS.get_mole_fractions()
        print(f"  Mole fractions: {' '.join(f'{x:.12f}' for x in mf)}")
        mf_mass = AS.get_mass_fractions()
        print(f"  Mass fractions: {' '.join(f'{x:.12f}' for x in mf_mass)}")

        AS.update(
            CP.QmassT_INPUTS, 0.5, 278.15
        )  # Saturated vapor (q=0.5 in mass basis)
        print(f"  Saturation pressure (mass basis): {AS.p():.12f} Pa")
        print(f"  Mole density: {AS.rhomolar() / 1e3:.12f} kmol/m3")
        print(f"  Mass density: {AS.rhomass():.12f} kg/m3")
        print(f"  Molar quality: {AS.Q():.12f} -")
        print(f"  Mass quality: {AS.Qmass():.12f} -")

        AS.update(CP.PQ_INPUTS, 1e5, 0.5)  # Saturated vapor (q=0.5 in molar basis)
        print(f"  Saturation temperature (molar basis): {AS.T():.12f} K")
        print(f"  Mole density: {AS.rhomolar() / 1e3:.12f} kmol/m3")
        print(f"  Mass density: {AS.rhomass():.12f} kg/m3")
        print(f"  Molar quality: {AS.Q():.12f} -")
        print(f"  Mass quality: {AS.Qmass():.12f} -")

        AS.update(CP.PQmass_INPUTS, 1e5, 0.5)  # Saturated vapor (q=0.5 in mass basis)
        print(f"  Saturation temperature (mass basis): {AS.T():.12f} K")
        print(f"  Mole density: {AS.rhomolar() / 1e3:.12f} kmol/m3")
        print(f"  Mass density: {AS.rhomass():.12f} kg/m3")
        print(f"  Molar quality: {AS.Q():.12f} -")
        print(f"  Mass quality: {AS.Qmass():.12f} -")

        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_phase_envelope_R410A()
    exit(0 if success else 1)
