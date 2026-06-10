#ifndef COOLPROP_DETAIL_STATE_CAPI_H
#define COOLPROP_DETAIL_STATE_CAPI_H

// C-ABI bridge for the frozen Cython `State` compat shim (the PDSim cimport
// surface).  A non-Cython core (here, the nanobind module) exports a pointer to
// this table as the `CoolProp._capi` PyCapsule; the thin `State` Cython shim
// grabs it at import and forwards every call through it -- so PDSim keeps
// cimporting `State` while the Cython AbstractState wrapper goes away.
//
// extern "C": only long/double/void* cross the boundary -- no C++ name
// mangling or std layout leaks, so the table is a stable ABI.  The opaque
// handle carries a `std::shared_ptr<AbstractState>` directly (no long-handle
// registry lookup).  See dev/state_capsule/ for the shim + contract harness.
//
// Why this indirection (measured, not assumed): the consuming shim is ~200 KB
// and contains NO CoolProp C++ -- it forwards to the single copy compiled into
// the core, so the wheel is not doubled.  A `State` that wrapped the C++
// AbstractState directly would be ~5 ns/read faster (~1-2% on PDSim's
// read-heavy ODE loop) but would embed its own ~7 MB of CoolProp, ~doubling the
// wheel.  The capsule trades that ~2% for no duplication + link-free
// decoupling; a consumer that can link CoolProp may instead wrap AbstractState
// directly against this same C-ABI.

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct  // NOLINT(modernize-use-using) -- C-ABI table; typedef kept C-compatible
    {
        void* (*make)(const char* backend, const char* fluids);
        void (*destroy)(void* handle);
        void (*update)(void* handle, long input_pair, double value1, double value2);
        double (*keyed_output)(void* handle, long key);
        double (*first_partial_deriv)(void* handle, long Of, long Wrt, long Constant);
        // CoolProp throws C++ exceptions; they cannot cross this extern "C" boundary.
        // Each call above catches them and stashes the message; last_error() returns
        // it (or NULL) so the Cython shim can re-raise a Python exception -- matching
        // the legacy State's `except *` behaviour.  Cleared on the next successful call.
        const char* (*last_error)();
        // Set the composition (mole fractions) of a mixture handle, so the shim's
        // set_Fluid can honour bracketed strings like "R32[0.5]&R134a[0.5]".
        // Appended after last_error to keep the existing field offsets stable.
        void (*set_mole_fractions)(void* handle, const double* fractions, long n);
        // Impose a phase (a `phases` enum value) on the handle, so the shim's
        // ``State(..., phase=...)`` can force the gas/liquid root like legacy.
        void (*specify_phase)(void* handle, long phase);
        // Lift a previously-imposed phase so the backend resumes auto-detecting it.
        // (specify_phase(iphase_not_imposed) is NOT a substitute: it also clobbers
        // the cached _phase, breaking lazy property evaluation.)  copy() uses this
        // to reproduce a mixture state under an imposed phase, then restore
        // auto-detect on the returned handle.  Appended last to keep field offsets
        // stable.
        void (*unspecify_phase)(void* handle);
    } CoolProp_StateCAPI;

#ifdef __cplusplus
}
#endif

#endif  // COOLPROP_DETAIL_STATE_CAPI_H
