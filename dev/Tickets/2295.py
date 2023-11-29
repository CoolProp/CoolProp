import CoolProp as CP


p = 1e6
for fluid in ["PHE", "Water", "TVP1", "DowQ"]:
    CP.CoolProp.set_debug_level(0)

    state = CP.AbstractState("INCOMP", fluid)

    Tmax = state.trivial_keyed_output(CP.iT_max)
    Tmin = state.trivial_keyed_output(CP.iT_min)
    Tbase = (Tmax + Tmin) *  0.5

    state.update(CP.PT_INPUTS, p, Tbase + 1e-6)
    h_hi = state.hmass()
    s_hi = state.smass()
    state.update(CP.PT_INPUTS, p, Tbase - 1e-6)
    h_lo = state.hmass()
    s_lo = state.smass()

    state.update(CP.PT_INPUTS, p, Tbase)
    try:
        print(state.hmass(), "J/kg   vs", (h_hi + h_lo) * 0.5, "J/kg")
        print(state.smass(), "J/kg/K vs", (s_hi + s_lo) * 0.5, "J/kg/K")
    except ValueError as e:
        # msg = (
        #     "Exactly the middle between maximum and minimum temperature does "
        #     f"not work for {fluid}."
        # )
        # print(msg)
        print(str(e))
        print((h_hi + h_lo) * 0.5, "J/kg and", (s_hi + s_lo) * 0.5, "J/kg/K")
