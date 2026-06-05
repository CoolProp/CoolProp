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
    } CoolProp_StateCAPI;

#ifdef __cplusplus
}
#endif

#endif  // COOLPROP_DETAIL_STATE_CAPI_H
