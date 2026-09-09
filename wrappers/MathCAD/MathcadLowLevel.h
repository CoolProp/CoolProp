// MathcadLowLevel.h : Low-Level (AbstractState) API functions for the Mathcad wrapper.
//
// This is an included implementation fragment, not a standalone header -- it
// is meant to be #include'd from exactly one place, partway through
// CoolPropMathcad.cpp, after that file has already set up:
//   - mcadincl.h (and the MC_STRING/STRING substitution trick above it)
//   - CoolProp/CoolProp.h (for CoolProp::set_error_string, strsplit) and
//     CoolProp/Configuration.h (for get_config_string/LIST_STRING_DELIMITER)
//   - CoolProp/DataStructures.h (for get_parameter_information,
//     get_input_pair_short_desc, get_phase_index -- used below for
//     proactive index/name validation rather than substring-sniffing
//     whatever exception AbstractState's own dispatch happens to throw)
//   - enum EC and CPErrorMessageTable (for BAD_HANDLE, BAD_PARAMETER,
//     BAD_INPUT_PAIR, BAD_FLUID, BAD_BACKEND, BAD_PHASE, NOT_MIXTURE,
//     BAD_FRACTION_SUM, UNEQUAL_LENGTH, TOO_MANY_OUTPUTS, INV_PARAMETER_IDX,
//     INV_INPUT_PAIR_STR, INV_INPUT_PAIR_IDX, MAKELRESULT)
//   - the general Mathcad wrapper helpers: CheckRealOrError,
//     CheckRealArrayOrError, get_nan, AllocateToMathcadArray
// It is kept separate from CoolPropMathcad.cpp purely to keep that file from
// growing unbounded as more Low-Level functions are added -- see that file's
// own "Low-Level (AbstractState) API Functions" include point for the usage
// banner (authoring patterns, handle-as-scalar rationale, etc.), which is
// not repeated here.
//
// These wrap CoolProp's handle-based low-level C API
// (CoolProp/CoolPropLib.h), together with MathcadStateGuard (a small,
// Mathcad-SDK-independent registry giving AS_factory() "replace, don't
// leak" semantics across worksheet recalculations).

#ifndef MATHCAD_LOWLEVEL_H
#define MATHCAD_LOWLEVEL_H

#include "CoolProp/CoolPropLib.h"
#include "MathcadStateGuard.h"

// Fixed size of the local errcode/message_buffer used by every Low-Level
// (AbstractState) wrapper function below to call into CoolPropLib.h.
constexpr long AS_ERR_BUFFER_LEN = 500;

// Helper: round a Mathcad complex scalar's real part to the nearest integer
// and return it as a long.  Used both for AbstractState handles and for the
// input-pair/parameter index arguments to the Low-Level API functions below
// -- Mathcad has no integer type, so all of these are carried through a
// worksheet as an ordinary real scalar (imag == 0, checked separately via
// CheckRealOrError).
static inline long RoundToLong(LPCCOMPLEXSCALAR val) {
    return static_cast<long>(std::llround(val->real));
}

// Helper: validate that a Mathcad complex scalar's real part is finite, then
// round it to a long -- combines the NaN/Inf guard with RoundToLong() so
// every Handle/InputPairIdx/ParamIdx conversion below gets both checks
// together, rather than risking one being added without the other.  Mathcad
// can legitimately carry NaN (this file's own AS_props_multi uses get_nan()
// to mark a failed point in its output array), so a non-finite value
// reaching one of these arguments is a real, reachable failure mode, not a
// hypothetical one: std::llround() on NaN/Inf, or the subsequent narrowing
// static_cast<long> of an out-of-range magnitude, is undefined behavior
// that in practice aliases to handle 0 (or another live handle/index)
// rather than erroring -- silently operating on the WRONG object is worse
// than a crash.  `code` is the caller's EC value for this specific
// argument (BAD_HANDLE/INV_INPUT_PAIR_IDX/INV_PARAMETER_IDX), so the Custom
// Error stays specific to what was actually being converted.
static LRESULT ToLongOrError(LPCCOMPLEXSCALAR val, EC code, unsigned int position, long* out) {
    if (!std::isfinite(val->real)) return MAKELRESULT(code, position);
    *out = RoundToLong(val);
    return 0;
}

// Helper: translate a CoolPropLib.h AbstractState_* errcode/message_buffer
// pair (errcode != 0) into a Mathcad EC error code, mirroring the
// sentinel-based error handling CP_Props1SI/CP_PropsSI/etc. use in
// CoolPropMathcad.cpp.
//
// BAD_HANDLE covers the two CoolProp::HandleError messages
// ("could not get handle" / "could not free handle") AbstractStateLibrary
// (src/CoolPropLib.cpp) raises for a stale/invalid handle -- both are
// prefixed "HandleError:" by HandleException() there.  Everything else
// becomes the generic LOWLEVEL_ERROR, with the low-level API's own message
// routed through set_error_string() so get_global_param_string("errstring")
// still retrieves full detail.  Unlike HandleProps1SIError/HandlePropsSIError
// in CoolPropMathcad.cpp, this does NOT attempt substring-matching for a more
// specific code: the low-level API's message text does not follow the same
// shape as the ValueError text those helpers were written against.
static LRESULT TranslateASError(const char* message, unsigned int position) {
    CoolProp::set_error_string(message);
    if (std::strncmp(message, "HandleError:", 12) == 0) {
        return MAKELRESULT(BAD_HANDLE, position);
    }
    return MAKELRESULT(LOWLEVEL_ERROR, position);
}

// Helper: translate an AS_factory() creation failure (from
// MathcadStateGuard::get_or_create(), i.e. AbstractState_factory()/
// AbstractState::factory()) into a Mathcad EC error code, distinguishing an
// invalid Backend string (argument 1) from an invalid Fluids string
// (argument 2) by substring-matching AbstractState::factory()'s own wording
// -- the same best-effort substring-matching approach
// HandleProps1SIError/HandlePropsSIError in CoolPropMathcad.cpp already use
// for this class of error.  Falls back to TranslateASError's generic
// handling (BAD_HANDLE/LOWLEVEL_ERROR, pointed at argument 1) for anything
// not recognized.
static LRESULT TranslateFactoryError(const char* message) {
    if (std::strstr(message, "Invalid backend name") != nullptr) {
        CoolProp::set_error_string(message);
        return MAKELRESULT(BAD_BACKEND, 1);
    }
    if (std::strstr(message, "not found") != nullptr) {
        CoolProp::set_error_string(message);
        return MAKELRESULT(BAD_FLUID, 2);
    }
    return TranslateASError(message, 1);
}

// Helper: true if idx is a registered CoolProp::parameters index (i.e. a
// value AS_param_index could have returned).  Proactive validation via the
// same lookup get_param_index()/keyed_output() ultimately rely on, rather
// than substring-matching whatever exception AbstractState::keyed_output()'s
// internal dispatch happens to throw for an unregistered key -- which,
// depending on exactly how unregistered, can throw one of two differently
// worded messages (see AbstractState::keyed_output()'s default case and
// get_parameter_information() in src/DataStructures.cpp).
static inline bool IsValidParamIndex(long idx) {
    try {
        // get_parameter_information() is [[nodiscard]]; only its throw-or-not
        // behavior matters here, so explicitly discard the returned string.
        (void)CoolProp::get_parameter_information(static_cast<int>(idx), "short");
        return true;
    } catch (...) {
        return false;
    }
}

// Helper: true if idx is a registered CoolProp::input_pairs index (i.e. a
// value AS_input_pair_index could have returned).  Same rationale as
// IsValidParamIndex above.
static inline bool IsValidInputPairIndex(long idx) {
    try {
        CoolProp::get_input_pair_short_desc(static_cast<CoolProp::input_pairs>(idx));
        return true;
    } catch (...) {
        return false;
    }
}

// Process-wide registry giving AS_factory() "get-or-create" (memoized)
// semantics across worksheet recalculations: recalculating the same
// AS_factory() cell with the same (Backend, Fluids) returns the SAME live
// handle rather than rebuilding the backend, so it doesn't pay construction
// cost -- up to 80-140 ms for tabular backends -- on every recalculation.
// See MathcadStateGuard.h.
static MathcadStateGuard as_state_guard;

// This code executes the user function CP_AS_factory, which is a wrapper for
// AbstractState_factory(), used to get (or, the first time, create) a
// persistent low-level fluid/mixture state and return an integer handle (as
// a real scalar) for use by the other AS_* functions below.  Recalculating
// this cell with the same Backend/Fluids returns the SAME handle rather than
// rebuilding the backend -- see MathcadStateGuard.h.
static LRESULT CP_AS_factory(LPCOMPLEXSCALAR Handle,  // output: handle for use by the other AS_* functions
                             LPCMCSTRING Backend,     // backend to use, e.g. "HEOS", "REFPROP", "BICUBIC&HEOS"
                             LPCMCSTRING Fluids)       // '&' delimited list of fluids
{
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];

    long handle = as_state_guard.get_or_create(Backend->str, Fluids->str, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateFactoryError(msg);

    Handle->real = static_cast<double>(handle);
    Handle->imag = 0;

    // normal return
    return 0;
}

// Fixed tolerance for the AS_set_fractions sum-to-1.0 check below.  Loose
// enough to tolerate ordinary floating-point representation of fractions a
// user typed to a handful of decimal places, tight enough to catch a
// genuinely missing/mistyped component.
constexpr double AS_FRACTION_SUM_TOLERANCE = 1e-6;

// This code executes the user function CP_AS_set_fractions, which is a wrapper for
// AbstractState_set_fractions(), used to set the mole/mass/volume fractions for a
// mixture handle created by AS_factory.  Returns Handle unchanged so downstream
// cells that use this call's return value depend on it.
static LRESULT CP_AS_set_fractions(LPCOMPLEXSCALAR HandleOut,   // output: Handle, unchanged
                                   LPCCOMPLEXSCALAR Handle,     // AbstractState handle from AS_factory
                                   LPCCOMPLEXARRAY Fractions)   // mole/mass/volume fractions
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealArrayOrError(Fractions, 2);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    // Look up how many fluids this handle actually has, so the two most
    // common misuses -- calling AS_set_fractions on a pure fluid at all, or
    // passing the wrong number of fractions for the mixture -- get a clear,
    // specific error instead of whatever message
    // AbstractState::set_mole_fractions()/set_mass_fractions() happens to
    // throw, routed through the generic LOWLEVEL_ERROR fallback.
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    char namesBuf[AS_ERR_BUFFER_LEN];
    AbstractState_fluid_names(handle, namesBuf, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    const std::string delimiter = CoolProp::get_config_string(LIST_STRING_DELIMITER);
    const std::string namesStr(namesBuf);
    if (namesStr.find(delimiter) == std::string::npos) {
        // A single name with no delimiter -- a pure (or pseudo-pure) fluid.
        return MAKELRESULT(NOT_MIXTURE, 1);
    }
    const long nFluids = static_cast<long>(strsplit(namesStr, delimiter[0]).size());
    if (static_cast<long>(Fractions->rows) != nFluids) {
        return MAKELRESULT(UNEQUAL_LENGTH, 2);
    }

    std::vector<double> fracVec(Fractions->hReal[0], Fractions->hReal[0] + Fractions->rows);

    double sum = 0.0;
    for (double f : fracVec) {
        sum += f;
    }
    if (std::fabs(sum - 1.0) > AS_FRACTION_SUM_TOLERANCE) {
        return MAKELRESULT(BAD_FRACTION_SUM, 2);
    }

    AbstractState_set_fractions(handle, fracVec.data(), static_cast<long>(fracVec.size()), &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_specify_phase, which is a wrapper for
// AbstractState_specify_phase(), used to impose a phase on a handle created by
// AS_factory for all subsequent AS_update/AS_props/AS_props_multi calls -- call
// this before any of those, once per handle.  Returns Handle unchanged so
// downstream cells that use this call's return value depend on it.
static LRESULT CP_AS_specify_phase(LPCOMPLEXSCALAR HandleOut,  // output: Handle, unchanged
                                   LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                   LPCMCSTRING Phase)          // phase name: "phase_liquid", "phase_gas", "phase_twophase", "phase_supercritical",
                                                               // "phase_supercritical_gas", "phase_supercritical_liquid", "phase_critical_point",
                                                               // "phase_unknown", or "phase_not_imposed" (CoolProp::phases, DataStructures.h)
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    // Validate the phase name proactively (rather than substring-matching
    // whatever AbstractState_specify_phase()'s internal
    // CoolProp::get_phase_index() call throws) so a bad Phase string is
    // unambiguously reported as BAD_PHASE, pointed at argument 2.
    try {
        CoolProp::get_phase_index(Phase->str);
    } catch (const CoolProp::ValueError& e) {
        CoolProp::set_error_string(e.what());
        return MAKELRESULT(BAD_PHASE, 2);
    }

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    AbstractState_specify_phase(handle, Phase->str, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_unspecify_phase, which is a wrapper for
// AbstractState_unspecify_phase(), used to remove a phase imposed by
// AS_specify_phase from a handle created by AS_factory.  Returns Handle
// unchanged so downstream cells that use this call's return value depend on it.
static LRESULT CP_AS_unspecify_phase(LPCOMPLEXSCALAR HandleOut,   // output: Handle, unchanged
                                     LPCCOMPLEXSCALAR Handle)     // AbstractState handle from AS_factory
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    AbstractState_unspecify_phase(handle, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_free, which is a wrapper for
// AbstractState_free(), used to explicitly release a handle created by
// AS_factory.  Safe as the last statement of a Mathcad program block;
// calling it from an independent worksheet cell is discouraged since nothing
// guarantees it runs after every reader of the same handle -- rely on
// AS_factory's registry guard to bound leakage there instead (see README.md).
static LRESULT CP_AS_free(LPCOMPLEXSCALAR Dummy,     // output (dummy value, 0 on success)
                          LPCCOMPLEXSCALAR Handle)   // AbstractState handle to release
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    AbstractState_free(handle, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Dummy->real = 0;
    Dummy->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_param_index, which is a wrapper for
// get_param_index(), used to resolve an output parameter name (e.g. "T", "Dmolar")
// to the integer index AS_get/AS_props/AS_props_multi expect -- resolve once
// anywhere in the worksheet, reuse many times, so a Mathcad array formula
// evaluating many points never marshals a string per point.
static LRESULT CP_AS_param_index(LPCOMPLEXSCALAR Index,  // output: parameter index
                                 LPCMCSTRING Name)        // parameter name, e.g. "T", "Dmolar", "Hmass"
{
    long idx = get_param_index(Name->str);
    if (idx < 0) return MAKELRESULT(BAD_PARAMETER, 1);

    Index->real = static_cast<double>(idx);
    Index->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_input_pair_index, which is a wrapper for
// get_input_pair_index(), used to resolve an input pair name (e.g. "PT_INPUTS") to
// the integer index AS_update/AS_props/AS_props_multi expect.
static LRESULT CP_AS_input_pair_index(LPCOMPLEXSCALAR Index,  // output: input pair index
                                      LPCMCSTRING Name)        // input pair name, e.g. "PT_INPUTS", "HmassP_INPUTS"
{
    long idx = get_input_pair_index(Name->str);
    if (idx < 0) return MAKELRESULT(INV_INPUT_PAIR_STR, 1);

    Index->real = static_cast<double>(idx);
    Index->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_update, which is a wrapper for
// AbstractState_update(), used to move a handle created by AS_factory to a
// new input point without reading any output yet.  Returns Handle unchanged
// (same convention as AS_set_fractions) so subsequent AS_get() calls that use
// this call's return value as their own Handle argument are guaranteed by
// Mathcad's dependency tracking to see the updated state.  An alternative to
// AS_props/AS_props_multi when several outputs are wanted from the same
// point: update once here, then call AS_get as many times as needed without
// re-running the flash for each one.
static LRESULT CP_AS_update(LPCOMPLEXSCALAR HandleOut,      // output: Handle, unchanged
                            LPCCOMPLEXSCALAR Handle,        // AbstractState handle from AS_factory
                            LPCCOMPLEXSCALAR InputPairIdx,  // input pair index, from AS_input_pair_index
                            LPCCOMPLEXSCALAR Value1,        // first input value
                            LPCCOMPLEXSCALAR Value2)        // second input value
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealOrError(InputPairIdx, 2);
    if (r) return r;
    r = CheckRealOrError(Value1, 3);
    if (r) return r;
    r = CheckRealOrError(Value2, 4);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;
    long inputPair;
    r = ToLongOrError(InputPairIdx, INV_INPUT_PAIR_IDX, 2, &inputPair);
    if (r) return r;
    if (!IsValidInputPairIndex(inputPair)) return MAKELRESULT(INV_INPUT_PAIR_IDX, 2);

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    AbstractState_update(handle, inputPair, Value1->real, Value2->real, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_get, which is a wrapper for
// AbstractState_keyed_output(), used to read one output parameter from a
// handle's CURRENT state -- pair with AS_update to move the state once and
// then read as many outputs as needed with separate AS_get calls, as a
// leaner alternative to AS_props/AS_props_multi.
static LRESULT CP_AS_get(LPCOMPLEXSCALAR Prop,       // output: the requested value
                         LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory/AS_update
                         LPCCOMPLEXSCALAR ParamIdx)  // output parameter index, from AS_param_index
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealOrError(ParamIdx, 2);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;
    long paramIdx;
    r = ToLongOrError(ParamIdx, INV_PARAMETER_IDX, 2, &paramIdx);
    if (r) return r;
    if (!IsValidParamIndex(paramIdx)) return MAKELRESULT(INV_PARAMETER_IDX, 2);

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    double value = AbstractState_keyed_output(handle, paramIdx, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Prop->real = value;
    Prop->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_props, which fuses
// AbstractState_update() and AbstractState_keyed_output() into one call for
// one-off/interactive use against a handle created by AS_factory -- so a
// single Mathcad cell is atomic regardless of which authoring pattern is in
// use.  All arguments after Handle are pre-resolved integer indices (from
// AS_input_pair_index/AS_param_index): the hot path here is pure numeric, no
// string marshalling per call.  See AS_update/AS_get above for a leaner
// alternative when several outputs are wanted from the same point.
static LRESULT CP_AS_props(LPCOMPLEXSCALAR Prop,           // output: computed value
                           LPCCOMPLEXSCALAR Handle,        // AbstractState handle from AS_factory
                           LPCCOMPLEXSCALAR InputPairIdx,  // input pair index, from AS_input_pair_index
                           LPCCOMPLEXSCALAR Value1,        // first input value
                           LPCCOMPLEXSCALAR Value2,        // second input value
                           LPCCOMPLEXSCALAR ParamIdx)      // output parameter index, from AS_param_index
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealOrError(InputPairIdx, 2);
    if (r) return r;
    r = CheckRealOrError(Value1, 3);
    if (r) return r;
    r = CheckRealOrError(Value2, 4);
    if (r) return r;
    r = CheckRealOrError(ParamIdx, 5);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;
    long inputPair;
    r = ToLongOrError(InputPairIdx, INV_INPUT_PAIR_IDX, 2, &inputPair);
    if (r) return r;
    long paramIdx;
    r = ToLongOrError(ParamIdx, INV_PARAMETER_IDX, 5, &paramIdx);
    if (r) return r;
    if (!IsValidInputPairIndex(inputPair)) return MAKELRESULT(INV_INPUT_PAIR_IDX, 2);
    if (!IsValidParamIndex(paramIdx)) return MAKELRESULT(INV_PARAMETER_IDX, 5);

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];

    AbstractState_update(handle, inputPair, Value1->real, Value2->real, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    double value = AbstractState_keyed_output(handle, paramIdx, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Prop->real = value;
    Prop->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_props_multi, which fuses a
// batched AbstractState_update_and_5_out() call into one Mathcad call for
// array/vectorized use against a handle created by AS_factory -- the main
// "minimum overhead" workhorse: one native flash loop evaluates every point,
// with no per-point Mathcad call and no per-point string marshalling.
static LRESULT CP_AS_props_multi(LPCOMPLEXARRAY Prop,            // output: matrix, one row per input point, one column per requested output
                                 LPCCOMPLEXSCALAR Handle,        // AbstractState handle from AS_factory
                                 LPCCOMPLEXSCALAR InputPairIdx,  // input pair index, from AS_input_pair_index
                                 LPCCOMPLEXARRAY Value1Array,    // first input values (Array)
                                 LPCCOMPLEXARRAY Value2Array,    // second input values (Array)
                                 LPCCOMPLEXARRAY ParamIdxArray)  // 1-5 output parameter indices, from AS_param_index (Array)
{
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealOrError(InputPairIdx, 2);
    if (r) return r;
    r = CheckRealArrayOrError(Value1Array, 3);
    if (r) return r;
    r = CheckRealArrayOrError(Value2Array, 4);
    if (r) return r;
    r = CheckRealArrayOrError(ParamIdxArray, 5);
    if (r) return r;

    // make sure input arrays are the same length (mirrors CP_PropsSImulti in CoolPropMathcad.cpp)
    if (Value1Array->rows != Value2Array->rows) {
        return MAKELRESULT(UNEQUAL_LENGTH, 4);
    }

    const long nOut = static_cast<long>(ParamIdxArray->rows);
    if (nOut < 1) return MAKELRESULT(BAD_PARAMETER, 5);
    if (nOut > 5) return MAKELRESULT(TOO_MANY_OUTPUTS, 5);

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;
    long inputPair;
    r = ToLongOrError(InputPairIdx, INV_INPUT_PAIR_IDX, 2, &inputPair);
    if (r) return r;
    if (!IsValidInputPairIndex(inputPair)) return MAKELRESULT(INV_INPUT_PAIR_IDX, 2);

    const long N = static_cast<long>(Value1Array->rows);
    std::vector<double> Value1Vec(Value1Array->hReal[0], Value1Array->hReal[0] + N);
    std::vector<double> Value2Vec(Value2Array->hReal[0], Value2Array->hReal[0] + N);

    // AbstractState_update_and_5_out always evaluates exactly 5 output keys per
    // point (src/CoolPropLib.cpp) -- pad any unused slots with a repeat of the
    // first requested (valid) key rather than a sentinel, so a bad/unused
    // padding index can never throw mid-point and blank out a later,
    // genuinely requested output for that same point: all 5 keyed_output()
    // calls for one point share a single try/catch, in order.  Validate each
    // genuinely requested (non-padding) index up front -- both that it's
    // finite (a raw array read, not a COMPLEXSCALAR, so not covered by
    // ToLongOrError above; NaN/Inf here hits the same std::llround() UB
    // ToLongOrError's own doc comment describes) and that it's a registered
    // parameter -- so a bad one is a clean BAD_PARAMETER Custom Error
    // instead of a LOWLEVEL_ERROR surfacing from deep inside the batch call.
    long outputs[5];
    for (long i = 0; i < nOut; ++i) {
        const double entry = ParamIdxArray->hReal[0][i];
        if (!std::isfinite(entry)) return MAKELRESULT(INV_PARAMETER_IDX, 5);
        outputs[i] = static_cast<long>(std::llround(entry));
        if (!IsValidParamIndex(outputs[i])) return MAKELRESULT(INV_PARAMETER_IDX, 5);
    }
    for (long i = nOut; i < 5; ++i) {
        outputs[i] = outputs[0];
    }

    // Pre-fill with NaN (not left default-zero) so a point the batch call
    // silently leaves untouched -- e.g. because the flash failed for that
    // point -- reads as an obvious sentinel in the returned matrix, matching
    // CP_PropsSImulti's _HUGE -> NaN convention in CoolPropMathcad.cpp.
    const double NaN = get_nan();
    std::vector<double> out1(N, NaN), out2(N, NaN), out3(N, NaN), out4(N, NaN), out5(N, NaN);

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN];
    AbstractState_update_and_5_out(handle, inputPair, Value1Vec.data(), Value2Vec.data(), N, outputs, out1.data(), out2.data(), out3.data(),
                                    out4.data(), out5.data(), &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    // Transpose into row-per-point order: AllocateToMathcadArray/hReal are
    // column-row indexed (dest->hReal[col][row]), and Mathcad users expect a
    // results TABLE -- N rows (one per input pair), up to 5 columns (one per
    // requested output) -- not one row per output with the points running
    // across columns.
    const std::vector<double>* const outArrays[5] = {&out1, &out2, &out3, &out4, &out5};
    std::vector<std::vector<double>> IO(static_cast<size_t>(N), std::vector<double>(static_cast<size_t>(nOut)));
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < nOut; ++j) {
            IO[i][j] = (*outArrays[j])[i];
        }
    }

    // Copy the results into the output complex array Prop (see AllocateToMathcadArray in CoolPropMathcad.cpp).
    LRESULT rc = AllocateToMathcadArray(Prop, IO);
    if (rc) return rc;

    // normal return
    return 0;
}

// Return every currently-live Low-Level state's handle, as a column vector,
// and AS_list_states() (below)'s Backend|Fluids string, in the same order.
// See MathcadStateGuard::snapshot()'s comment for the ordering guarantee and
// its one caveat (a factory/free call landing between this call and
// AS_list_states() on the same worksheet -- rely on Recalculate Worksheet
// for a guaranteed-consistent pair, same advice already given for the
// worksheet-level-handle authoring pattern).
//
// `Trigger` is intentionally unused -- Mathcad Prime's Custom Function
// syntax requires at least one argument, and there is no *natural* one to
// key a whole-registry snapshot on, so this exists purely to satisfy that
// requirement (any real scalar works, e.g. a literal 0). It does double
// duty, though: wiring it to a handle already on the sheet (rather than a
// bare literal) gives Mathcad a real dependency edge, so this call re-runs
// whenever THAT handle's defining cell does -- a narrower, more reliable
// trigger than waiting for a full Recalculate Worksheet.
static LRESULT CP_AS_list_handles(LPCOMPLEXARRAY Handles, LPCCOMPLEXSCALAR Trigger) {
    (void)Trigger;
    auto snap = as_state_guard.snapshot();
    if (snap.empty()) return MAKELRESULT(NO_ACTIVE_STATES, 1);

    std::vector<std::vector<double>> Vec;
    Vec.reserve(snap.size());
    for (const auto& kv : snap) {
        Vec.push_back(std::vector<double>(1, static_cast<double>(kv.second)));
    }
    return AllocateToMathcadArray(Handles, Vec);
}

// Companion to AS_list_handles(): the same live states' "Backend|Fluids" key
// (the same string AS_factory()'s two arguments were joined into), ";"-
// delimited, in the same order AS_list_handles() returns their handles.
// `Trigger` is unused -- see CP_AS_list_handles()'s comment above.
static LRESULT CP_AS_list_states(LPMCSTRING States, LPCCOMPLEXSCALAR Trigger) {
    (void)Trigger;
    auto snap = as_state_guard.snapshot();
    if (snap.empty()) return MAKELRESULT(NO_ACTIVE_STATES, 1);

    std::vector<std::string> keys;
    keys.reserve(snap.size());
    for (const auto& kv : snap) {
        keys.push_back(kv.first);
    }
    States->str = AllocMathcadString(strjoin(keys, ";"));
    return 0;
}

// ********************************************************************************************************
// Fill out FUNCTIONINFO structures for the Low-Level (AbstractState) API functions above
// ********************************************************************************************************

FUNCTIONINFO ASFactory = {
  const_cast<char*>("AS_factory"),                                                                  // Name by which Mathcad will recognize the function
  const_cast<char*>("Backend, Fluids"),                                                             // Description of input parameters
  const_cast<char*>("Creates a persistent Low-Level fluid/mixture state and returns a handle"),      // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_factory,                                                                       // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar (the handle)
  2,                                                                                                 // Number of arguments
  {MC_STRING, MC_STRING}                                                                            // Argument types
};

FUNCTIONINFO ASSetFractions = {
  const_cast<char*>("AS_set_fractions"),                                                            // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Fractions"),                                                            // Description of input parameters
  const_cast<char*>("Sets the mole/mass/volume fractions for a mixture Handle; returns Handle"),     // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_set_fractions,                                                                  // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar (Handle, unchanged)
  2,                                                                                                 // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_ARRAY}                                                                    // Argument types
};

FUNCTIONINFO ASSpecifyPhase = {
  const_cast<char*>("AS_specify_phase"),                                                             // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Phase"),                                                                // Description of input parameters
  const_cast<char*>("Imposes a phase on a Low-Level state Handle for subsequent updates; returns Handle"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_specify_phase,                                                                  // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar (Handle, unchanged)
  2,                                                                                                 // Number of arguments
  {COMPLEX_SCALAR, MC_STRING}                                                                        // Argument types
};

FUNCTIONINFO ASUnspecifyPhase = {
  const_cast<char*>("AS_unspecify_phase"),                                                           // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle"),                                                                       // Description of input parameters
  const_cast<char*>("Removes a phase imposed by AS_specify_phase from a Low-Level state Handle; returns Handle"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_unspecify_phase,                                                                // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar (Handle, unchanged)
  1,                                                                                                 // Number of arguments
  {COMPLEX_SCALAR}                                                                                   // Argument types
};

FUNCTIONINFO ASFree = {
  const_cast<char*>("AS_free"),                                                                      // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle"),                                                                       // Description of input parameters
  const_cast<char*>("Releases a Low-Level state Handle created by AS_factory"),                     // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_free,                                                                           // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar (dummy value)
  1,                                                                                                 // Number of arguments
  {COMPLEX_SCALAR}                                                                                   // Argument types
};

FUNCTIONINFO ASParamIndex = {
  const_cast<char*>("AS_param_index"),                                                               // Name by which Mathcad will recognize the function
  const_cast<char*>("Name"),                                                                         // Description of input parameters
  const_cast<char*>("Returns the integer index for an output parameter name, e.g. \"T\", \"Dmolar\""),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_param_index,                                                                    // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar
  1,                                                                                                 // Number of arguments
  {MC_STRING}                                                                                        // Argument types
};

FUNCTIONINFO ASInputPairIndex = {
  const_cast<char*>("AS_input_pair_index"),                                                          // Name by which Mathcad will recognize the function
  const_cast<char*>("Name"),                                                                         // Description of input parameters
  const_cast<char*>("Returns the integer index for an input pair name, e.g. \"PT_INPUTS\""),         // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_input_pair_index,                                                               // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                   // Returns a Mathcad complex scalar
  1,                                                                                                 // Number of arguments
  {MC_STRING}                                                                                        // Argument types
};

FUNCTIONINFO ASUpdate = {
  const_cast<char*>("AS_update"),                                                                                  // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Input Pair Index, Input Property 1, Input Property 2"),                               // Description of input parameters
  const_cast<char*>("Updates a Low-Level state Handle to a new input point (no output read); returns Handle"),     // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_update,                                                                                       // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                                  // Returns a Mathcad complex scalar (Handle, unchanged)
  4,                                                                                                                // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_SCALAR}                                                 // Argument types
};

FUNCTIONINFO ASGet = {
  const_cast<char*>("AS_get"),                                                                                     // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Output Parameter Index"),                                                             // Description of input parameters
  const_cast<char*>("Returns one output parameter from a Low-Level state Handle's current point"),                // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get,                                                                                          // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                                  // Returns a Mathcad complex scalar
  2,                                                                                                                // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                                 // Argument types
};

FUNCTIONINFO ASProps = {
  const_cast<char*>("AS_props"),                                                                                          // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Input Pair Index, Input Property 1, Input Property 2, Output Parameter Index"),              // Description of input parameters
  const_cast<char*>("Updates a Low-Level state Handle and returns one output parameter"),                                 // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_props,                                                                                                // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                                          // Returns a Mathcad complex scalar
  5,                                                                                                                       // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_SCALAR}                                        // Argument types
};

FUNCTIONINFO ASPropsMulti = {
  const_cast<char*>("AS_props_multi"),                                                                                                    // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Input Pair Index, Input Property 1 (Array), Input Property 2 (Array), Output Parameter Indices (Array)"),     // Description of input parameters
  const_cast<char*>("Updates a Low-Level state Handle for a range of inputs and returns up to 5 output parameters as a table (row per input point, column per output)"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_props_multi,                                                                                                          // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                                                           // Returns a Mathcad complex array
  5,                                                                                                                                        // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR, COMPLEX_ARRAY, COMPLEX_ARRAY, COMPLEX_ARRAY}                                                            // Argument types
};

FUNCTIONINFO ASListHandles = {
  const_cast<char*>("AS_list_handles"),                                                              // Name by which Mathcad will recognize the function
  const_cast<char*>("Trigger"),                                                                       // Description of input parameters (unused -- see function comment)
  const_cast<char*>("Returns a column vector of all currently-live Low-Level state Handles"),         // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_list_handles,                                                                    // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                       // Returns a Mathcad complex array
  1,                                                                                                   // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR}                                                                                     // Argument types
};

FUNCTIONINFO ASListStates = {
  const_cast<char*>("AS_list_states"),                                                                                // Name by which Mathcad will recognize the function
  const_cast<char*>("Trigger"),                                                                                      // Description of input parameters (unused -- see function comment)
  const_cast<char*>("Returns \";\"-delimited \"Backend|Fluids\" for all currently-live states, matching AS_list_handles() order"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_list_states,                                                                                     // Pointer to the function code.
  MC_STRING,                                                                                                          // Returns a Mathcad string
  1,                                                                                                                  // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR}                                                                                                    // Argument types
};

#endif  // MATHCAD_LOWLEVEL_H
