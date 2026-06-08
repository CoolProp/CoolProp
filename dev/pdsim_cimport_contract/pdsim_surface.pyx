# distutils: language = c++
# cython: language_level=3
"""
PDSim cimport contract — a frozen miniature of how ibell/pdsim consumes
CoolProp's *Cython* surface at the C level.

Why this file exists
--------------------
PDSim does not use the high-level Python API; it ``cimport``s CoolProp's
``State`` / ``AbstractState`` cdef classes and the ``constants_header`` enum
module, then calls their methods at C speed inside its ODE loops.  That makes
those declarations a binary/source *contract*.  The v8 reorg (dropping the
hand-written Cython interface in favour of a nanobind core + a thin frozen
``State`` shim) MUST preserve this surface or PDSim breaks.

This module reproduces *exactly* the members PDSim touches (extracted from
ibell/pdsim @ main).  Two guard rails:

  * COMPILE-TIME: if a v8 change drifts the cimportable ``State`` /
    ``AbstractState`` / ``constants_header`` surface, this file fails to
    compile -> the test errors out.  That is the "PDSim need not change"
    guarantee, proven mechanically.
  * BEHAVIOURAL: test_pdsim_contract.py runs the functions below and asserts
    the values match the high-level ``PropsSI`` reference, proving any future
    capsule-forwarding implementation is numerically correct.

Keep this list in sync with SURFACE.md.  Anything added here should be a
deliberate, reviewed widening of the contract.
"""
from CoolProp.State cimport State                   # the cdef class (CoolProp.State path)
cimport CoolProp.constants_header as constants      # enum module (constants.iXxx)
from CoolProp.constants_header cimport parameters   # enum used as a C type


cpdef dict exercise(object Fluid, double T, double rho):
    """Touch every single-phase State + .pAS member PDSim uses; return values."""
    cdef State S = State(Fluid, {'T': T, 'D': rho})     # dict construction
    cdef State S2

    # --- the three update variants PDSim calls ---
    S.update_Trho(T, rho)
    cdef double p_from_Trho = S.get_p()
    S.update({'T': T, 'D': rho})                        # dict update
    S.update_ph(S.get_p(), S.get_h())                   # ph update

    # --- direct cdef attribute reads (cached doubles + py attrs) ---
    cdef double Tc = S.T_
    cdef double pc = S.p_
    cdef double rc = S.rho_
    cdef object fluid_name = S.Fluid
    cdef object phase = S.phase

    # --- the State getter surface ---
    cdef dict getters = {
        'T': S.get_T(), 'p': S.get_p(), 'h': S.get_h(), 'rho': S.get_rho(),
        's': S.get_s(), 'u': S.get_u(), 'cp': S.get_cp(), 'cp0': S.get_cp0(),
        'cv': S.get_cv(), 'MM': S.get_MM(), 'dpdT': S.get_dpdT(),
        'visc': S.get_visc(), 'cond': S.get_cond(),
        'speed_sound': S.get_speed_sound(),
    }

    # --- Props(enum) ---
    cdef double h_via_Props = S.Props(constants.iHmass)

    # --- reach THROUGH .pAS to the AbstractState (PDSim does this 24x) ---
    cdef double rho_pAS = S.pAS.rhomass()
    cdef double cp_pAS = S.pAS.cpmass()
    cdef double cv_pAS = S.pAS.cvmass()
    cdef double T_pAS = S.pAS.T()
    cdef double p_pAS = S.pAS.p()
    cdef double ko_h = S.pAS.keyed_output(constants.iHmass)
    cdef object names = S.pAS.fluid_names()
    S.pAS.update(constants.PT_INPUTS, pc * 1000.0, Tc)  # AbstractState.update (SI: Pa)
    # mirror State.get_dpdT exactly: first_partial_deriv(iP,iT,iDmolar) in SI
    cdef double dpdT_pAS = S.pAS.first_partial_deriv(
        constants.iP, constants.iT, constants.iDmolar)

    # --- copy() returns a fresh State ---
    S2 = S.copy()
    cdef double h_copy = S2.get_h()

    return {
        'p_from_Trho': p_from_Trho, 'T_': Tc, 'p_': pc, 'rho_': rc,
        'Fluid': fluid_name, 'phase': phase, 'getters': getters,
        'h_via_Props': h_via_Props, 'rho_pAS': rho_pAS, 'cp_pAS': cp_pAS,
        'cv_pAS': cv_pAS, 'T_pAS': T_pAS, 'p_pAS': p_pAS, 'ko_h': ko_h,
        'names': list(names), 'dpdT_pAS': dpdT_pAS, 'h_copy': h_copy,
    }


cpdef dict exercise_twophase(object Fluid, double T, double Q):
    """get_Q + two-phase construction (the one getter undefined single-phase)."""
    cdef State S = State(Fluid, {'T': T, 'Q': Q})
    return {'Q': S.get_Q(), 'T': S.get_T(), 'p': S.get_p()}


cpdef dict exercise_refrigeration(object Fluid, double P_kPa, double T):
    """The refrigeration-cycle surface PDSim cimports, widened for v8
    (bd CoolProp-r9sq.26): bare ``State(Fluid, None)`` construction, the
    cimport-level ``set_Fluid``, the integer ``Phase()``, and the saturation
    quantities ``Tsat`` / ``subcooling`` / ``superheat`` (P in legacy kPa)."""
    cdef State S = State(Fluid, None)               # bare construction (no state update)
    S.set_Fluid(Fluid, b'HEOS')                     # cimport-level set_Fluid
    S.update({'P': P_kPa, 'T': T})
    cdef long ph = S.Phase()                        # integer phase flag
    return {
        'phase': ph,
        'Tsat': S.get_Tsat(1.0),
        'subcooling': S.get_subcooling(),
        'superheat': S.get_superheat(),
    }


cpdef double hot_loop(object Fluid, double T, double rho, int n):
    """Mirror PDSim's inner pattern: repeated update_Trho + getters + .pAS."""
    cdef State S = State(Fluid, {'T': T, 'D': rho})
    cdef double acc = 0.0
    cdef int i
    for i in range(n):
        S.update_Trho(T + 0.001 * i, rho)
        acc += S.get_p() + S.get_h() + S.pAS.keyed_output(constants.iSmass)
    return acc
