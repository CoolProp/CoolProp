"""
Full-surface PARITY check of the State-capsule shim against the legacy
``CoolProp.State`` unit conventions, for single-phase points of standard fluids.

    Cython State shim -> CoolProp._capi capsule -> nanobind core -> C++ AbstractState

WHY THIS SHAPE
--------------
The legacy ``State`` getters are NOT uniformly SI -- each carries its own unit
convention: kPa/kJ for the calorific quantities, **g/mol** for ``get_MM``,
**kW/m/K** for ``get_cond``, but plain SI for ``get_rho``/``get_T``/
``get_visc``/``get_speed_sound``.

An earlier version of this test checked ``get_MM``/``get_visc``/``get_cond``
only for *finiteness*.  That let two 1000x unit errors pass silently --
``get_MM`` returned kg/mol instead of g/mol and ``get_cond`` returned W/m/K
instead of kW/m/K -- which only surfaced far downstream as a complex-sqrt
blow-up inside PDSim's isentropic-nozzle ``(k*R*T)**0.5`` (negative ``k`` from
a 1000x-too-large gas constant).

So this version:

  * asserts EVERY getter's value against the core under its documented
    convention, to rtol 1e-10 (shim and PropsSI share the backend, so a correct
    shim matches to ~machine precision); and
  * adds an independent *structural* guard on the two unit-bearing getters that
    fails if the shim returns the raw SI value -- a guard that holds regardless
    of state, so a scale error cannot hide symmetrically in both the shim and
    its PropsSI-derived expectation.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NB_BUILD = os.path.abspath(os.path.join(HERE, "..", "..", "wrappers", "Python", "nanobind", "build"))

# Single-phase points: liquid water and supercritical CO2 (Tc(CO2)=304 K < 320 K,
# so any density is single-phase) -- both unambiguous under (T, Dmass) inputs.
SINGLE_PHASE = [("Water", 320.0, 1000.0), ("CO2", 320.0, 600.0)]


def _check_fluid(CP, State, fluid, T, RHO, fails):
    """Assert the full getter surface of the shim against the core's values."""
    def SI(name):
        return CP.PropsSI(name, "T", T, "Dmass", RHO, fluid)

    def expect(label, got, want, rel=1e-10):
        if not math.isclose(got, want, rel_tol=rel):
            fails.append(f"[{fluid}] {label}: got {got:.12g}, want {want:.12g} (rel_tol {rel:g})")

    S = State.State(fluid, {"T": T, "D": RHO})

    # SI getters (no scaling)
    expect("get_T [K, SI]", S.get_T(), T, rel=1e-12)
    expect("get_rho [kg/m^3, SI]", S.get_rho(), RHO, rel=1e-12)
    expect("get_speed_sound [m/s, SI]", S.get_speed_sound(), SI("speed_of_sound"))
    expect("get_visc [Pa.s, SI]", S.get_visc(), SI("viscosity"))

    # kPa / kJ getters (SI / 1000)
    expect("get_p [kPa]", S.get_p(), SI("P") / 1000.0)
    expect("get_h [kJ/kg]", S.get_h(), SI("Hmass") / 1000.0)
    expect("get_s [kJ/kg/K]", S.get_s(), SI("Smass") / 1000.0)
    expect("get_u [kJ/kg]", S.get_u(), SI("Umass") / 1000.0)
    expect("get_cp [kJ/kg/K]", S.get_cp(), SI("Cpmass") / 1000.0)
    expect("get_cp0 [kJ/kg/K]", S.get_cp0(), SI("Cp0mass") / 1000.0)
    expect("get_cv [kJ/kg/K]", S.get_cv(), SI("Cvmass") / 1000.0)

    # Unit-bearing getters that previously hid 1000x bugs:
    si_mm = SI("molar_mass")          # kg/mol
    si_cond = SI("conductivity")      # W/m/K
    expect("get_MM [g/mol]", S.get_MM(), si_mm * 1000.0)
    expect("get_cond [kW/m/K]", S.get_cond(), si_cond / 1000.0)
    # Independent structural guards (state-agnostic): legacy g/mol is ~1000x the
    # kg/mol SI value, and legacy kW/m/K is ~1000x smaller than the W/m/K SI value.
    # A raw-SI passthrough (the historical bug) makes these equal -> guard fails.
    if not S.get_MM() > si_mm * 100.0:
        fails.append(f"[{fluid}] get_MM returned SI kg/mol ({S.get_MM():.6g}), not legacy g/mol")
    if not S.get_cond() < si_cond / 100.0:
        fails.append(f"[{fluid}] get_cond returned SI W/m/K ({S.get_cond():.6g}), not legacy kW/m/K")

    # dpdT: legacy convention is first_partial_deriv(iP,iT,iDmolar)/1000 [kPa/K]
    ref = CP.AbstractState("HEOS", fluid)
    ref.update(int(CP.DmassT_INPUTS), RHO, T)
    fpd = ref.first_partial_deriv(int(CP.iP), int(CP.iT), int(CP.iDmolar))
    expect("get_dpdT [kPa/K]", S.get_dpdT(), fpd / 1000.0)

    # cached attrs + Python property accessors mirror the getters
    expect("cached T_", S.T_, T, rel=1e-12)
    expect("cached rho_", S.rho_, RHO, rel=1e-12)
    expect("cached p_ [kPa]", S.p_, S.get_p(), rel=1e-12)
    expect("property .p", S.p, S.get_p(), rel=1e-12)
    expect("property .MM", S.MM, S.get_MM(), rel=1e-12)
    expect("property .k (cond)", S.k, S.get_cond(), rel=1e-12)

    # Props() / .pAS are SI (no scaling)
    expect("Props(iHmass) [SI]", S.Props(int(CP.iHmass)), SI("Hmass"))
    expect("pAS.keyed_output(iHmass) [SI]", S.pAS.keyed_output(int(CP.iHmass)), SI("Hmass"))
    expect("pAS.rhomass [SI]", S.pAS.rhomass(), RHO, rel=1e-12)
    expect("pAS.cpmass [SI]", S.pAS.cpmass(), SI("Cpmass"))
    expect("pAS.cvmass [SI]", S.pAS.cvmass(), SI("Cvmass"))

    # copy() preserves state
    expect("copy().get_h [kJ/kg]", S.copy().get_h(), SI("Hmass") / 1000.0)

    names = S.pAS.fluid_names()
    if not (names and fluid.lower() in names[0].lower()):
        fails.append(f"[{fluid}] fluid_names mismatch: {names}")


def main():
    r = subprocess.run([sys.executable, "setup.py", "build_ext", "--inplace"],
                       cwd=HERE, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write("State shim FAILED TO BUILD:\n" + r.stdout[-4000:] + r.stderr[-4000:])
        return 1

    sys.path.insert(0, NB_BUILD)   # nanobind core (exports _capi)
    sys.path.insert(0, HERE)       # State.*.so
    import CoolProp as CP          # noqa: E402
    import State                   # noqa: E402

    fails = []
    for fluid, T, RHO in SINGLE_PHASE:
        _check_fluid(CP, State, fluid, T, RHO, fails)

    # --- structural exercises (input pairs, two-phase, hot loop, cimport) ---
    def SI_water(name):
        return CP.PropsSI(name, "T", 320.0, "Dmass", 1000.0, "Water")

    # ph round-trip (kPa/kJ in) returns to the same state
    S = State.State("Water", {"T": 320.0, "D": 1000.0})
    S.update_ph(S.get_p(), S.get_h())
    if not math.isclose(S.get_rho(), 1000.0, rel_tol=1e-6):
        fails.append(f"update_ph round-trip rho {S.get_rho()}")

    # two-phase construction + quality
    S2 = State.State("Water", {"T": 373.0, "Q": 0.4})
    if not math.isclose(S2.get_Q(), 0.4, rel_tol=1e-9):
        fails.append(f"two-phase get_Q {S2.get_Q()}")
    if not math.isclose(S2.get_p(), CP.PropsSI("P", "T", 373.0, "Q", 0.4, "Water") / 1000.0, rel_tol=1e-10):
        fails.append(f"two-phase get_p {S2.get_p()}")

    # hot loop (PDSim inner pattern): update_Trho + getters stay finite/consistent
    S3 = State.State("Water", {"T": 320.0, "D": 1000.0})
    acc = 0.0
    for i in range(1000):
        S3.update_Trho(320.0 + 0.001 * i, 1000.0)
        acc += S3.get_p() + S3.get_h() + S3.pAS.keyed_output(int(CP.iSmass))
    if not (math.isfinite(acc) and acc != 0.0):
        fails.append(f"hot_loop acc {acc}")

    # --- mixtures (GH #3151) -------------------------------------------------
    # HEOS has no (T, D) flash for mixtures (DHSU_T_flash does not support
    # mixtures yet), so a bare (T, D) mixture construction must raise that REAL
    # error -- not the bogus "'str' object has no attribute 'decode'" that the
    # broken _raise_if_error used to mask it with.
    MIX = "HEOS::R32[0.7]&R125[0.3]"
    try:
        State.State(MIX, {"T": 300.0, "D": 30.0})
        fails.append("mixture (T,D) without phase unexpectedly succeeded")
    except ValueError as e:
        msg = str(e)
        if "decode" in msg or "AttributeError" in msg:
            fails.append(f"mixture error still masked by decode bug: {msg}")
        elif "mixtures" not in msg:
            fails.append(f"mixture (T,D) raised an unexpected error: {msg}")
    except Exception as e:  # not even a ValueError -> the mask is back
        fails.append(f"mixture (T,D) raised {type(e).__name__}, not a clear ValueError: {e}")

    # With a phase hint the (T, D) update is a direct evaluation and works; the
    # composition must be honoured (fractions forwarded to the backend).
    Smix = State.State(MIX, {"T": 300.0, "D": 30.0}, phase="gas")
    if not math.isclose(Smix.get_T(), 300.0, rel_tol=1e-12):
        fails.append(f"mixture get_T {Smix.get_T()}")
    if not math.isclose(Smix.get_rho(), 30.0, rel_tol=1e-12):
        fails.append(f"mixture get_rho {Smix.get_rho()}")
    if sorted(n.lower() for n in Smix.pAS.fluid_names()) != ["r125", "r32"]:
        fails.append(f"mixture fluid_names {Smix.pAS.fluid_names()}")

    # copy() must forward the composition AND reproduce the phase (PDSim's
    # guess_outlet_temp path) -- bit-for-bit the same state as the source.
    Cmix = Smix.copy()
    for name, a, b in (("T", Smix.get_T(), Cmix.get_T()),
                       ("p", Smix.get_p(), Cmix.get_p()),
                       ("rho", Smix.get_rho(), Cmix.get_rho()),
                       ("h", Smix.get_h(), Cmix.get_h()),
                       ("s", Smix.get_s(), Cmix.get_s())):
        if not math.isclose(a, b, rel_tol=1e-12):
            fails.append(f"mixture copy() {name}: source {a:.12g} != copy {b:.12g}")
    # the copy must have lifted the imposed phase: a later update on a different
    # (supported) input pair auto-detects, i.e. does not crash with a stale phase.
    try:
        Cmix.update({"P": Cmix.get_p(), "T": Cmix.get_T() + 20.0})
        _ = Cmix.get_h()  # lazy calc under the (now auto-detected) phase
    except Exception as e:
        fails.append(f"mixture copy() left an imposed phase (later update failed): {e}")

    # cimport-level: PDSim's `cdef State` + cdef method/.pAS calls compile & run
    import cimport_check   # noqa: E402
    cc = cimport_check.check("Water", 320.0, 1000.0)
    if not math.isclose(cc["p"], SI_water("P") / 1000.0, rel_tol=1e-10):
        fails.append(f"cimport cdef get_p {cc['p']}")
    if not math.isclose(cc["rhomass"], 1000.0, rel_tol=1e-9):
        fails.append(f"cimport cdef pas.rhomass {cc['rhomass']}")

    if fails:
        sys.stderr.write("\n*** SHIM PARITY MISMATCHES ***\n  " + "\n  ".join(fails) + "\n")
        return 1
    print("STATE CAPSULE: PASS — full getter surface matches the legacy State "
          "conventions to rtol 1e-10 for single-phase Water & CO2")
    print("  (incl. get_MM in g/mol and get_cond in kW/m/K — the 1000x conventions "
          "that the old finiteness-only check missed)")
    print("  + cimport-level: `cdef State` + cdef method/.pAS calls compile & run (PDSim pattern)")
    return 0


def test_capsule():
    """pytest entry point."""
    assert main() == 0


if __name__ == "__main__":
    raise SystemExit(main())
