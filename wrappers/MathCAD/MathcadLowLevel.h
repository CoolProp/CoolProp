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
//     BAD_FRACTION_SUM, ZERO_FRACTION_SUM, UNEQUAL_LENGTH, TOO_MANY_OUTPUTS, INV_PARAMETER_IDX,
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

#include <algorithm>  // std::max_element, used by CP_AS_pe_tmax/CP_AS_pe_pmax
#include <limits>     // std::numeric_limits<long>, used by TryRoundToLong()
#include <mutex>      // std::mutex/std::scoped_lock, used by g_as_mutex below

// Fixed size of the local errcode/message_buffer used by every Low-Level
// (AbstractState) wrapper function below to call into CoolPropLib.h.
constexpr long AS_ERR_BUFFER_LEN = 500;

// Helper: round a plain double to the nearest integer and return it as a
// long, validating BOTH that it's finite AND that the rounded value is
// actually representable in a long -- the single checked-conversion path
// shared by ToLongOrError() below (for COMPLEXSCALAR Handle/InputPairIdx/
// ParamIdx arguments) and CP_AS_props_multi()'s own per-entry ParamIdxArray
// loop (a raw array read, not a COMPLEXSCALAR, so it can't go through
// ToLongOrError directly). Mathcad has no integer type, so all of these
// arrive as an ordinary real scalar -- and Mathcad can legitimately carry
// NaN (this file's own AS_props_multi uses get_nan() to mark a failed point
// in its output array) or an arbitrarily large finite magnitude, so both
// are real, reachable failure modes, not hypothetical ones: std::llround()
// on a non-finite value, or narrowing an in-range-for-long-long-but-
// out-of-range-for-`long` result via static_cast<long>, is undefined
// behavior that in practice aliases to handle 0 (or another live handle/
// index) rather than erroring -- silently operating on the WRONG object is
// worse than a crash. Returns false (leaving *out untouched) for either
// failure; the caller maps that to whichever Custom Error fits its
// argument.
static inline bool TryRoundToLong(double real, long* out) {
    if (!std::isfinite(real)) return false;
    // Reject anything outside `long`'s range (with 0.5 of slack either side
    // for correct rounding at the boundary) BEFORE calling std::llround() --
    // not just after narrowing its result.  A finite double can be far
    // outside even long long's representable range (roughly +-9.2e18,
    // against double's own range out to ~1.8e308), and the standard leaves
    // std::llround()'s return value unspecified in that case: checking the
    // post-round long long against long's range (as this function used to)
    // relies on llround() saturating in practice rather than being
    // guaranteed to by the standard.  Pre-filtering to (roughly) long's own
    // range keeps the value llround() actually sees always representable in
    // long long, so the call itself is well-defined, not just its result's
    // subsequent range check.
    constexpr double kLongMin = static_cast<double>((std::numeric_limits<long>::min)());
    constexpr double kLongMax = static_cast<double>((std::numeric_limits<long>::max)());
    if (real < kLongMin - 0.5 || real > kLongMax + 0.5) return false;
    const long long rounded = std::llround(real);
    // NOT dead code: the pre-check above only bounds `real`, not `rounded`.
    // At the exact half-integer boundary (real == kLongMax + 0.5, say),
    // the pre-check's strict `>` lets it through, but round-half-away-from-
    // zero then produces kLongMax + 1 -- one past what fits in `long`. This
    // check is what actually rejects that case; removing it as "redundant"
    // with the pre-check reintroduces the narrowing bug this function
    // exists to close.
    if (rounded < static_cast<long long>((std::numeric_limits<long>::min)())
        || rounded > static_cast<long long>((std::numeric_limits<long>::max)())) {
        return false;
    }
    *out = static_cast<long>(rounded);
    return true;
}

// Helper: validate/convert a Mathcad complex scalar's real part into a long
// via TryRoundToLong() above, for every Handle/InputPairIdx/ParamIdx
// argument below. `code` is the caller's EC value for this specific
// argument (BAD_HANDLE/INV_INPUT_PAIR_IDX/INV_PARAMETER_IDX), so the Custom
// Error stays specific to what was actually being converted.
static LRESULT ToLongOrError(LPCCOMPLEXSCALAR val, EC code, unsigned int position, long* out) {
    if (!TryRoundToLong(val->real, out)) return MAKELRESULT(code, position);
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
//
// Every `msg`/`namesBuf` buffer passed in here (and to every AbstractState_*
// call throughout this file) is zero-initialized at its declaration --
// HandleException() (src/CoolPropLib.cpp) only memcpy's into message_buffer
// when the formatted text fits; on its "didn't fit" path (errcode==2, not
// reachable in practice at AS_ERR_BUFFER_LEN=500 but not provably
// unreachable either) the buffer is left untouched, and the strncmp() below
// must not read uninitialized stack memory in that case.
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

// Process-wide registry giving AS_factory() "get-or-create" (cached)
// semantics across worksheet recalculations: recalculating the same
// AS_factory() call with the same (Backend, Fluids) returns the SAME live
// handle rather than rebuilding the backend, so it doesn't pay construction
// cost -- up to 80-140 ms for tabular backends -- on every recalculation.
// See MathcadStateGuard.h.
static MathcadStateGuard as_state_guard;

// Process-wide lock serializing every CP_AS_* call below that touches a
// handle or the as_state_guard registry (every one of them except the three
// pure lookups that never reference a handle: CP_AS_param_index,
// CP_AS_input_pair_index, CP_AS_generate_update_pair).
//
// Why this is needed: MathcadStateGuard's own mutex (see MathcadStateGuard.h)
// protects only ITS registry map -- it is released before AbstractState_*
// functions ever touch the underlying AbstractState object.
// AbstractStateLibrary::get() (src/CoolPropLib.cpp) is the same story: its
// mutex protects the handle table lookup, then releases before the caller
// reads or mutates the object the returned pointer refers to. Mathcad Prime
// can call custom DLL functions from more than one worksheet-recalculation
// thread at once (nothing in the Mathcad SDK registration used here declares
// these functions as needing single-threaded dispatch), so two callbacks
// that happen to reference the SAME handle -- e.g. one cell's AS_props and
// another's AS_get, both reading/writing the shared AbstractState h refers
// to -- could otherwise interleave: an update from one call landing between
// another call's update and its output read, silently returning a value for
// the WRONG input point rather than erroring. AbstractState itself documents
// no thread-safety guarantee for concurrent calls on one instance, so this
// is a real, not hypothetical, race.
//
// Fix: a single global lock, held for the ENTIRE body of every handle-
// touching CP_AS_* function (acquired as each function's first statement,
// released automatically on every return path via std::scoped_lock's
// destructor) -- not just around the individual AbstractState_* calls
// inside it, and not two narrower locks taken separately by, say,
// CP_AS_update and CP_AS_get, which would still let a third call interleave
// between them. This fully serializes the Low-Level API surface: correct
// regardless of whether Mathcad Prime's multi-threaded worksheet
// calculation is enabled, at the cost of no longer benefiting from it for
// these calls specifically (they queue rather than overlap). Given
// AS_props_multi already batches an entire array's worth of points into one
// call specifically to avoid needing many parallel Low-Level calls, that
// cost is small next to the correctness it buys.
//
// Note this is a single PROCESS-wide lock, not one per handle: while one
// call is in progress, every other handle-touching call blocks too, even
// against a completely unrelated handle/fluid -- e.g. CP_AS_build_phase_envelope
// tracing one handle's envelope holds this lock for that whole trace, which
// briefly serializes an unrelated AS_get on a different handle behind it.
// A deliberate simplicity-over-throughput trade-off: per-handle locking
// would avoid that unrelated-handle stall, but a single mutex is trivially
// easy to reason about (see the deadlock argument above) and each call here
// is a thin wrapper, not a Mathcad-facing hot loop -- AS_props_multi is
// still the answer for genuinely high-throughput array evaluation.
static std::mutex g_as_mutex;

// This code executes the user function CP_AS_factory, which is a wrapper for
// AbstractState_factory(), used to get (or, the first time, create) a
// persistent low-level fluid/mixture state and return an integer handle (as
// a real scalar) for use by the other AS_* functions below.  Recalculating
// this call with the same Backend/Fluids returns the SAME handle rather than
// rebuilding the backend -- see MathcadStateGuard.h.
static LRESULT CP_AS_factory(LPCOMPLEXSCALAR Handle,  // output: handle for use by the other AS_* functions
                             LPCMCSTRING Backend,     // backend to use, e.g. "HEOS", "REFPROP", "BICUBIC&HEOS"
                             LPCMCSTRING Fluids)       // '&' delimited list of fluids
{
    std::scoped_lock lock(g_as_mutex);
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};

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
// equations that use this call's return value depend on it.
static LRESULT CP_AS_set_fractions(LPCOMPLEXSCALAR HandleOut,   // output: Handle, unchanged
                                   LPCCOMPLEXSCALAR Handle,     // AbstractState handle from AS_factory
                                   LPCCOMPLEXARRAY Fractions)   // mole/mass/volume fractions
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealArrayOrError(Fractions, 2);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    // Look up how many fluids this handle actually has, so the two most
    // common misuses -- calling AS_set_fractions on a pure fluid, or
    // passing the wrong number of fractions for the mixture -- get a clear,
    // specific error instead of whatever message
    // AbstractState::set_mole_fractions()/set_mass_fractions() happens to
    // throw, routed through the generic LOWLEVEL_ERROR fallback.
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    char namesBuf[AS_ERR_BUFFER_LEN] = {};
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

// Helper: molar mass (kg/mol) of each fluid in `handle`'s mixture, in the
// same order AbstractState_fluid_names() lists them -- shared by
// CP_AS_mole_to_mass_fractions()/CP_AS_mass_to_mole_fractions() below.
//
// Deliberately avoids adding any new CoolPropLib.h export: per-component
// molar mass isn't exposed through the handle itself, but each component's
// NAME is (AbstractState_fluid_names(), already used by CP_AS_set_fractions
// above), and CoolProp::Props1SI() -- a plain, handle-independent,
// name-based lookup already used by CP_Props1SI() elsewhere in this file --
// resolves "molar_mass" for any of them directly. No new shared C API
// surface needed for what is, underneath, the same computation
// AbstractState::calc_mass_fractions() already does in the C++ API (mass_i
// = mm_i * x_i / sum(mm_j * x_j)) -- just re-derived here from
// already-exposed building blocks instead of wrapping that C++-only method.
static LRESULT GetComponentMolarMasses(long handle, std::vector<double>* molarMasses, unsigned int position) {
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    char namesBuf[AS_ERR_BUFFER_LEN] = {};
    AbstractState_fluid_names(handle, namesBuf, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, position);

    const std::string delimiter = CoolProp::get_config_string(LIST_STRING_DELIMITER);
    const std::string namesStr(namesBuf);
    const std::vector<std::string> names =
      (namesStr.find(delimiter) == std::string::npos) ? std::vector<std::string>{namesStr} : strsplit(namesStr, delimiter[0]);

    molarMasses->clear();
    molarMasses->reserve(names.size());
    for (const auto& name : names) {
        double mm = CoolProp::Props1SI(name, "molar_mass");
        if (!ValidNumber(mm)) {
            std::string emsg = CoolProp::get_global_param_string("errstring");
            CoolProp::set_error_string(emsg);
            return MAKELRESULT(LOWLEVEL_ERROR, position);
        }
        molarMasses->push_back(mm);
    }
    return 0;
}

// This code executes the user function CP_AS_mole_to_mass_fractions, which
// converts an arbitrary mole-fraction composition to the equivalent mass
// fractions for `handle`'s mixture (component identities and molar masses
// come from the handle; the fractions to convert are a separate argument,
// not whatever happens to already be set on the handle -- so this is usable
// as a preprocessing step before AS_set_fractions, not just as a read-back).
// Self-normalizing: divides by the actual weighted sum rather than assuming
// MoleFractions already sums to 1, so a not-quite-normalized input still
// produces a correctly-normalized result.
static LRESULT CP_AS_mole_to_mass_fractions(LPCOMPLEXARRAY MassFractions,  // output: column vector of mass fractions
                                            LPCCOMPLEXSCALAR Handle,      // AbstractState handle from AS_factory (for component identities)
                                            LPCCOMPLEXARRAY MoleFractions)  // column vector of mole fractions to convert
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealArrayOrError(MoleFractions, 2);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    std::vector<double> molarMasses;
    r = GetComponentMolarMasses(handle, &molarMasses, 1);
    if (r) return r;

    if (static_cast<size_t>(MoleFractions->rows) != molarMasses.size()) {
        return MAKELRESULT(UNEQUAL_LENGTH, 2);
    }

    double denom = 0.0;
    for (size_t i = 0; i < molarMasses.size(); ++i) {
        denom += molarMasses[i] * MoleFractions->hReal[0][i];
    }
    // Guard the normalization divide -- an all-zero (or exactly canceling)
    // MoleFractions input, trivially reachable even for a pure fluid via
    // MoleFractions=[0], would otherwise silently produce NaN in every
    // output element instead of a diagnosable Custom Error.
    if (denom == 0.0) {
        return MAKELRESULT(ZERO_FRACTION_SUM, 2);
    }

    std::vector<std::vector<double>> Vec(molarMasses.size());
    for (size_t i = 0; i < molarMasses.size(); ++i) {
        Vec[i] = {molarMasses[i] * MoleFractions->hReal[0][i] / denom};
    }
    return AllocateToMathcadArray(MassFractions, Vec);
}

// This code executes the user function CP_AS_mass_to_mole_fractions -- the
// inverse of CP_AS_mole_to_mass_fractions() above (mole_i = (w_i / mm_i) /
// sum(w_j / mm_j)). See that function's comment for the shared rationale.
static LRESULT CP_AS_mass_to_mole_fractions(LPCOMPLEXARRAY MoleFractions,  // output: column vector of mole fractions
                                            LPCCOMPLEXSCALAR Handle,      // AbstractState handle from AS_factory (for component identities)
                                            LPCCOMPLEXARRAY MassFractions)  // column vector of mass fractions to convert
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;
    r = CheckRealArrayOrError(MassFractions, 2);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    std::vector<double> molarMasses;
    r = GetComponentMolarMasses(handle, &molarMasses, 1);
    if (r) return r;

    if (static_cast<size_t>(MassFractions->rows) != molarMasses.size()) {
        return MAKELRESULT(UNEQUAL_LENGTH, 2);
    }

    double denom = 0.0;
    for (size_t i = 0; i < molarMasses.size(); ++i) {
        denom += MassFractions->hReal[0][i] / molarMasses[i];
    }
    // See the matching guard in CP_AS_mole_to_mass_fractions() above.
    if (denom == 0.0) {
        return MAKELRESULT(ZERO_FRACTION_SUM, 2);
    }

    std::vector<std::vector<double>> Vec(molarMasses.size());
    for (size_t i = 0; i < molarMasses.size(); ++i) {
        Vec[i] = {(MassFractions->hReal[0][i] / molarMasses[i]) / denom};
    }
    return AllocateToMathcadArray(MoleFractions, Vec);
}

// This code executes the user function CP_AS_specify_phase, which is a wrapper for
// AbstractState_specify_phase(), used to impose a phase on a handle created by
// AS_factory for all subsequent AS_update/AS_props/AS_props_multi calls -- call
// this before any of those, once per handle.  Returns Handle unchanged so
// downstream equations that use this call's return value depend on it.
static LRESULT CP_AS_specify_phase(LPCOMPLEXSCALAR HandleOut,  // output: Handle, unchanged
                                   LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                   LPCMCSTRING Phase)          // phase name: "phase_liquid", "phase_gas", "phase_twophase", "phase_supercritical",
                                                               // "phase_supercritical_gas", "phase_supercritical_liquid", "phase_critical_point",
                                                               // "phase_unknown", or "phase_not_imposed" (CoolProp::phases, DataStructures.h)
{
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
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
// unchanged so downstream equations that use this call's return value depend on it.
static LRESULT CP_AS_unspecify_phase(LPCOMPLEXSCALAR HandleOut,   // output: Handle, unchanged
                                     LPCCOMPLEXSCALAR Handle)     // AbstractState handle from AS_factory
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    AbstractState_unspecify_phase(handle, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_get_phase, which is a wrapper
// for AbstractState_phase() -- the read-only complement to
// AS_specify_phase()/AS_unspecify_phase() (which impose/remove a phase
// constraint): returns the phase the CURRENT point actually is in right
// now, as a string, e.g. "phase_liquid" -- the exact same string
// AS_specify_phase()'s Phase argument accepts, so the two round-trip.
// Useful for worksheet branching -- e.g. checking the state is actually
// two-phase before calling AS_get_sat_liquid/AS_mole_fractions_liquid --
// without relying on those raising a LOWLEVEL_ERROR to find out.
static LRESULT CP_AS_get_phase(LPMCSTRING PhaseStr,        // output: phase name, e.g. "phase_liquid"
                               LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                               LPCCOMPLEXSCALAR Trigger)   // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    int phase = AbstractState_phase(handle, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    PhaseStr->str = AllocMathcadString(CoolProp::get_phase_short_desc(static_cast<CoolProp::phases>(phase)));

    // normal return
    return 0;
}

// This code executes the user function CP_AS_free, which is a wrapper for
// AbstractState_free(), used to explicitly release a handle created by
// AS_factory.  Safe as the last statement of a Mathcad program block;
// calling it from an independent worksheet equation is discouraged since nothing
// guarantees it runs after every reader of the same handle -- rely on
// AS_factory's registry guard to bound leakage there instead (see README.md).
static LRESULT CP_AS_free(LPCOMPLEXSCALAR Dummy,     // output (dummy value, 0 on success)
                          LPCCOMPLEXSCALAR Handle)   // AbstractState handle to release
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
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
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
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
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
    double value = AbstractState_keyed_output(handle, paramIdx, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Prop->real = value;
    Prop->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_get_sat_liquid, which is a
// wrapper for AbstractState_saturated_liquid_keyed_output(), used to read
// one output parameter from the SATURATED LIQUID side of a handle's current
// two-phase state -- e.g. after an AS_update/AS_props call with a Q (quality)
// input. Distinct from AS_get, which reads the bulk/overall state.
static LRESULT CP_AS_get_sat_liquid(LPCOMPLEXSCALAR Prop,       // output: the requested value
                                    LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory/AS_update
                                    LPCCOMPLEXSCALAR ParamIdx)  // output parameter index, from AS_param_index
{
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
    double value = AbstractState_saturated_liquid_keyed_output(handle, paramIdx, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Prop->real = value;
    Prop->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_AS_get_sat_vapor, which is a
// wrapper for AbstractState_saturated_vapor_keyed_output() -- see
// CP_AS_get_sat_liquid()'s comment above; identical except it reads the
// SATURATED VAPOR side of the current two-phase state.
static LRESULT CP_AS_get_sat_vapor(LPCOMPLEXSCALAR Prop,       // output: the requested value
                                   LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory/AS_update
                                   LPCCOMPLEXSCALAR ParamIdx)  // output parameter index, from AS_param_index
{
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
    double value = AbstractState_saturated_vapor_keyed_output(handle, paramIdx, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    Prop->real = value;
    Prop->imag = 0;

    // normal return
    return 0;
}

// Helper: fetch a component-indexed vector -- mole fractions, or anything
// else sharing the same (values, maxN, N, errcode, message_buffer,
// buffer_length) "checked memory" contract most AbstractState_get_*()
// getters use -- into `out`, sized from the mixture's ACTUAL component
// count instead of a fixed guess. Shared by
// CP_AS_get_mole_fractions()/CP_AS_mole_fractions_liquid()/
// CP_AS_mole_fractions_vapor() below, so this probe/fetch dance (mirroring
// FetchPhaseEnvelope()'s, for the same reason -- see that function's
// comment for the full ordering argument) only needs writing once.
//
// `call` is a callable taking (double* values, long maxN, long* N, long*
// errcode, char* message_buffer) that forwards straight into one of those
// AbstractState_get_mole_fractions*() C functions with `handle` (and, for
// the satState variant, its "liquid"/"gas" argument) already bound --
// buffer_length is always AS_ERR_BUFFER_LEN here, so it's baked in by the
// caller's lambda rather than threaded through this helper's own signature.
//
// Two calls into `call`: the first probes with maxN=0 and a null buffer.
// Every getter this is used with writes *N* before comparing it to maxN and
// throwing if too small (see e.g. AbstractState_get_mole_fractions() in
// src/CoolPropLib.cpp) -- so maxN=0 reliably throws (a mixture always has
// >=1 component) but only AFTER *N already holds the true count, and before
// the write loop that would touch the output buffer ever runs -- making a
// null buffer safe for this probe, the same way FetchPhaseEnvelope()'s
// length/component probes rely on. The second call then allocates `out` to
// that real size and fetches for real, instead of capping at a fixed bound
// that would otherwise reject any mixture with more components than that
// guess (a 64-component bound is generous for realistic worksheets, but not
// a hard CoolProp limit -- a large predefined mixture or a many-component
// custom blend can exceed it).
template <typename Fn>
static LRESULT FetchComponentVector(Fn&& call, std::vector<double>* out) {
    long N = 0;
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    call(nullptr, 0, &N, &errcode, msg);
    if (N <= 0) {
        // Should be unreachable in practice -- see FetchPhaseEnvelope()'s
        // matching comment for the same reasoning -- but report whatever
        // this probe actually raised rather than proceeding with a bogus
        // zero count.
        if (errcode) return TranslateASError(msg, 1);
        return MAKELRESULT(LOWLEVEL_ERROR, 1);
    }

    out->resize(static_cast<size_t>(N));
    long finalN = 0;
    errcode = 0;
    char finalMsg[AS_ERR_BUFFER_LEN] = {};
    call(out->data(), N, &finalN, &errcode, finalMsg);
    if (errcode) return TranslateASError(finalMsg, 1);
    out->resize(static_cast<size_t>(finalN));
    return 0;
}

// This code executes the user function CP_AS_get_mole_fractions, which is a
// wrapper for AbstractState_get_mole_fractions() -- the handle's current
// BULK mole fractions (whatever AS_set_fractions last set, or the trivial
// [1] for a pure fluid). Distinct from AS_mole_fractions_liquid/vapor below,
// which read the saturated liquid/vapor side of a two-phase point, not the
// overall composition. Useful to read back what AS_set_fractions actually
// applied, or the composition of a handle built from a predefined-mixture
// string.
static LRESULT CP_AS_get_mole_fractions(LPCOMPLEXARRAY Fractions,   // output: column vector of mole fractions
                                        LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                        LPCCOMPLEXSCALAR Trigger)   // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    std::vector<double> fracVec;
    r = FetchComponentVector(
      [handle](double* buf, long maxN, long* N, long* errcode, char* msg) {
          AbstractState_get_mole_fractions(handle, buf, maxN, N, errcode, msg, AS_ERR_BUFFER_LEN);
      },
      &fracVec);
    if (r) return r;

    std::vector<std::vector<double>> Vec(fracVec.size());
    for (size_t i = 0; i < fracVec.size(); ++i) {
        Vec[i] = {fracVec[i]};
    }
    return AllocateToMathcadArray(Fractions, Vec);
}

// This code executes the user function CP_AS_mole_fractions_liquid, which is
// a wrapper for AbstractState_get_mole_fractions_satState() with
// saturated_state="liquid" -- the SATURATED LIQUID side's mole fractions at
// a handle's current two-phase state, as a column vector. Requires the
// state to actually be in the two-phase region (0 <= quality <= 1);
// AbstractState_get_mole_fractions_satState() enforces that itself and its
// message surfaces via LOWLEVEL_ERROR if not.
//
// Like AS_get_mole_fractions() above, sizes its output from the mixture's
// actual component count via FetchComponentVector() rather than a fixed
// guess -- see that helper's comment.
//
// `Trigger` is unused by this function's body, but NOT because Mathcad
// requires an argument -- Handle already satisfies that on its own. The
// real reason: Handle's own value never changes when the AbstractState it
// names is mutated in place -- AS_update/AS_props/AS_specify_phase all echo
// Handle back unchanged, by design (see AS_update's comment) -- so an equation
// whose only input is Handle gives Mathcad's dependency graph nothing to
// key a recalculation on when the underlying point moves. Wire Trigger to
// whatever value actually drives the state you want reflected here -- e.g.
// the quality value fed into the AS_update/AS_props call that put the state
// in the two-phase region this function reads -- so this equation re-evaluates
// whenever that does, instead of needing a full Recalculate Worksheet. If
// this equation already references the freshly-reassigned Handle from that same
// update (the normal chaining idiom), that alone may already provide the
// dependency edge; Trigger is the explicit fallback for call shapes where
// it doesn't.
static LRESULT CP_AS_mole_fractions_liquid(LPCOMPLEXARRAY Fractions,   // output: column vector of mole fractions
                                           LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                           LPCCOMPLEXSCALAR Trigger)   // unused -- see comment above for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    std::vector<double> fracVec;
    r = FetchComponentVector(
      [handle](double* buf, long maxN, long* N, long* errcode, char* msg) {
          AbstractState_get_mole_fractions_satState(handle, "liquid", buf, maxN, N, errcode, msg, AS_ERR_BUFFER_LEN);
      },
      &fracVec);
    if (r) return r;

    std::vector<std::vector<double>> Vec(fracVec.size());
    for (size_t i = 0; i < fracVec.size(); ++i) {
        Vec[i] = {fracVec[i]};
    }
    return AllocateToMathcadArray(Fractions, Vec);
}

// This code executes the user function CP_AS_mole_fractions_vapor -- see
// CP_AS_mole_fractions_liquid()'s comment above; identical except it wraps
// AbstractState_get_mole_fractions_satState() with saturated_state="gas"
// (the SATURATED VAPOR side).
static LRESULT CP_AS_mole_fractions_vapor(LPCOMPLEXARRAY Fractions,   // output: column vector of mole fractions
                                          LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                          LPCCOMPLEXSCALAR Trigger)   // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    std::vector<double> fracVec;
    r = FetchComponentVector(
      [handle](double* buf, long maxN, long* N, long* errcode, char* msg) {
          AbstractState_get_mole_fractions_satState(handle, "gas", buf, maxN, N, errcode, msg, AS_ERR_BUFFER_LEN);
      },
      &fracVec);
    if (r) return r;

    std::vector<std::vector<double>> Vec(fracVec.size());
    for (size_t i = 0; i < fracVec.size(); ++i) {
        Vec[i] = {fracVec[i]};
    }
    return AllocateToMathcadArray(Fractions, Vec);
}

// This code executes the user function CP_AS_generate_update_pair, which is
// a wrapper for CoolProp::generate_update_pair() (DataStructures.h) -- given
// two output-parameter indices in EITHER order, resolves which
// CoolProp::input_pairs the combination corresponds to (if any) and returns
// its name as a string, e.g. "PT_INPUTS", for further use with
// AS_input_pair_index/AS_update/AS_props/AS_props_multi. Unlike those,
// generate_update_pair() takes no Handle -- it is a pure lookup over the two
// parameter keys, not tied to any particular fluid/mixture state, so this
// function needs no Trigger argument either: its two real ParamIdx
// arguments already give Mathcad everything it needs to know when to re-run
// this equation.
//
// No Value1/Value2 arguments: generate_update_pair()'s own implementation
// picks the pair purely from key1/key2 (a long chain of
// match_pair(key1, key2, ...) checks against the two keys, nothing else) --
// its value1/value2 parameters exist only to be copied into out1/out2 in
// the resolved pair's order, which this function doesn't surface anyway
// (see below). Passing them would be dead weight, so this only takes the
// two indices. The dummy 0.0s below stand in for the unused value
// arguments generate_update_pair()'s signature still requires.
static LRESULT CP_AS_generate_update_pair(LPMCSTRING PairName,        // output: resolved input pair name, e.g. "PT_INPUTS"
                                          LPCCOMPLEXSCALAR ParamIdx1,  // first output parameter index, from AS_param_index
                                          LPCCOMPLEXSCALAR ParamIdx2)  // second output parameter index, from AS_param_index
{
    LRESULT r = CheckRealOrError(ParamIdx1, 1);
    if (r) return r;
    r = CheckRealOrError(ParamIdx2, 2);
    if (r) return r;

    long idx1;
    r = ToLongOrError(ParamIdx1, INV_PARAMETER_IDX, 1, &idx1);
    if (r) return r;
    if (!IsValidParamIndex(idx1)) return MAKELRESULT(INV_PARAMETER_IDX, 1);

    long idx2;
    r = ToLongOrError(ParamIdx2, INV_PARAMETER_IDX, 2, &idx2);
    if (r) return r;
    if (!IsValidParamIndex(idx2)) return MAKELRESULT(INV_PARAMETER_IDX, 2);

    double out1, out2;
    CoolProp::input_pairs pair =
      CoolProp::generate_update_pair(static_cast<CoolProp::parameters>(idx1), 0.0, static_cast<CoolProp::parameters>(idx2), 0.0, out1, out2);
    if (pair == CoolProp::INPUT_PAIR_INVALID) {
        return MAKELRESULT(NO_SUCH_INPUT_PAIR, 1);
    }

    PairName->str = AllocMathcadString(CoolProp::get_input_pair_short_desc(pair));

    // normal return
    return 0;
}

// This code executes the user function CP_AS_props, which fuses
// AbstractState_update() and AbstractState_keyed_output() into one call for
// one-off/interactive use against a handle created by AS_factory -- so a
// single Mathcad equation is atomic regardless of which authoring pattern is in
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
    std::scoped_lock lock(g_as_mutex);
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
    char msg[AS_ERR_BUFFER_LEN] = {};

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
    std::scoped_lock lock(g_as_mutex);
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
    // genuinely requested (non-padding) index up front via TryRoundToLong()
    // -- a raw array read, not a COMPLEXSCALAR, so it can't go through
    // ToLongOrError() above directly, but needs the exact same finite-AND-
    // in-range check that helper does -- and that it's a registered
    // parameter, so a bad one is a clean BAD_PARAMETER Custom Error instead
    // of a LOWLEVEL_ERROR surfacing from deep inside the batch call.
    long outputs[5];
    for (long i = 0; i < nOut; ++i) {
        const double entry = ParamIdxArray->hReal[0][i];
        if (!TryRoundToLong(entry, &outputs[i])) return MAKELRESULT(INV_PARAMETER_IDX, 5);
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
    char msg[AS_ERR_BUFFER_LEN] = {};
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
// whenever THAT handle's defining equation does -- a narrower, more reliable
// trigger than waiting for a full Recalculate Worksheet.
static LRESULT CP_AS_list_handles(LPCOMPLEXARRAY Handles, LPCCOMPLEXSCALAR Trigger) {
    std::scoped_lock lock(g_as_mutex);
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
    std::scoped_lock lock(g_as_mutex);
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

// This code executes the user function CP_AS_backend_name, which returns
// the short backend string (e.g. "HEOS", "REFPROP", "BICUBIC&HEOS") a
// single handle is using -- the same short form AS_factory's Backend
// argument takes. AS_list_states() already reports this for every
// currently-open handle at once (as the "Backend|Fluids" half of its key),
// so this is mainly a convenience when you only have one specific Handle in
// scope and don't want to fetch/parse the whole registry listing.
//
// Deliberately does NOT just return AbstractState_backend_name()'s value
// as-is: that call returns CoolProp's internal C++ implementation class
// name for the backend (get_backend_string() in src/DataStructures.cpp
// literally maps, e.g., HEOS_BACKEND_MIX to "HelmholtzEOSMixtureBackend"),
// which is correct but reads as an implementation detail to a Mathcad user
// expecting the same short string they typed into AS_factory. Recovering
// that short string requires CoolProp's backend-family lookup tables,
// which are private to DataStructures.cpp (no public header declares them,
// unlike get_phase_short_desc()/get_input_pair_short_desc() for phases and
// input pairs) -- rather than reaching into that internal machinery,
// MathcadStateGuard already remembers the short string verbatim: it's the
// "Backend" half of the "Backend|Fluids" key this handle was registered
// under. AbstractState_backend_name() is still called first, purely so
// this function validates/errors on a dead handle exactly like every other
// AS_* function does; its result becomes the fallback if, for any reason,
// the handle isn't found in the registry snapshot (confirmed live above,
// so this shouldn't normally happen, but a long name beats no answer).
static LRESULT CP_AS_backend_name(LPMCSTRING BackendStr,     // output: backend name, e.g. "HEOS"
                                  LPCCOMPLEXSCALAR Handle,   // AbstractState handle from AS_factory
                                  LPCCOMPLEXSCALAR Trigger)  // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    char backendBuf[AS_ERR_BUFFER_LEN] = {};
    AbstractState_backend_name(handle, backendBuf, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    for (const auto& kv : as_state_guard.snapshot()) {
        if (kv.second == handle) {
            const std::string& key = kv.first;
            auto sep = key.find('|');
            BackendStr->str = AllocMathcadString(sep == std::string::npos ? key : key.substr(0, sep));
            return 0;
        }
    }
    BackendStr->str = AllocMathcadString(std::string(backendBuf));  // fallback -- see comment above

    // normal return
    return 0;
}

// This code executes the user function CP_AS_build_phase_envelope, which is a
// wrapper for AbstractState_build_phase_envelope(), used to trace the phase
// envelope (dew/bubble curve) for a handle created by AS_factory before any
// call to AS_get_phase_envelope_data() on it.  Returns Handle unchanged so
// downstream equations that use this call's return value depend on it.
static LRESULT CP_AS_build_phase_envelope(LPCOMPLEXSCALAR HandleOut,  // output: Handle, unchanged
                                          LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                                          LPCMCSTRING Level)          // refinement level -- CoolProp recommends "none" (skip refining)
{
    std::scoped_lock lock(g_as_mutex);
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    AbstractState_build_phase_envelope(handle, Level->str, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    HandleOut->real = Handle->real;
    HandleOut->imag = 0;

    // normal return
    return 0;
}

// Holds one handle's traced phase envelope, fetched via FetchPhaseEnvelope()
// below. Per-component compositions (x/y) are deliberately not kept here --
// see FetchPhaseEnvelope()'s comment.
struct PhaseEnvelopeTPRho
{
    std::vector<double> T, P, rhomolar_vap, rhomolar_liq;
};

// Helper: fetch the phase envelope traced by a prior AS_build_phase_envelope
// call for `handle` -- shared by CP_AS_get_phase_envelope_data() and the
// cricondentherm/cricondenbar max-point functions below, so the probe/fetch
// dance only needs writing once. On success, returns 0 and `out` is sized to
// the actual point count; on failure, returns the LRESULT to propagate
// (BAD_HANDLE / LOWLEVEL_ERROR / PHASE_ENVELOPE_NOT_BUILT) and `out` is left
// however std::vector::resize leaves it (unspecified contents, not read by
// any caller that checks the return value first).
//
// Mathcad (or, for the two max-point functions, this function itself) must
// allocate its own output storage before it can be filled, but neither the
// point count nor the component count is known until AFTER the envelope has
// been built -- so this makes THREE calls into the same underlying C API
// function, each learning one more thing than the last.
//
// Call 1 (length probe): length=0, maxComponents=0, no output buffers.
// Looking at AbstractState_get_phase_envelope_data_checkedMemory()'s own
// implementation (src/CoolPropLib.cpp), it always writes *actual_length
// before it can either succeed or throw on the length check, so this probe
// safely reports the true point count even though it also reports an error
// (0 is essentially never a big enough buffer). *actual_components is NOT
// learned here -- that check, and the write to it, only happens AFTER the
// length check already passed, so a length=0 probe never reaches it.
//
// Call 2 (component probe): length=<the real point count from call 1>,
// maxComponents=0, still no output buffers. With the length check now
// passing, execution reaches the components check -- which writes
// *actual_components BEFORE comparing it to maxComponents=0 and throwing
// (a mixture always has >=1 component, so this throw is as reliable as call
// 1's). Crucially, the x/y write loop is reached only AFTER that comparison
// passes, so passing null x/y here is exactly as safe as passing null
// T/P/rhomolar_vap/rhomolar_liq was in call 1. This is what lets this
// function size its composition buffers from the mixture's ACTUAL component
// count instead of guessing/capping at a fixed bound -- a mixture with more
// components than that fixed bound would otherwise make call 3 below throw
// a LOWLEVEL_ERROR instead of returning the envelope,
// even though CoolProp itself has no trouble tracing it.
//
// Call 3 (the real fetch): length and maxComponents both sized from what
// calls 1-2 actually learned, output buffers sized to match. The
// per-component x/y buffers this call still requires are allocated and
// immediately discarded -- no caller of this helper needs them: an N x
// Ncomp matrix per phase is a meaningfully different, more complex shape
// than what any of them return.
//
// All three calls happen while the caller (every CP_AS_* function that
// reaches this helper) already holds g_as_mutex, so they're atomic with
// respect to any other Low-Level call on this or any other handle -- no
// other thread's AS_build_phase_envelope()/AS_free() etc. can run between
// them and change what these three calls see out from under this function.
static LRESULT FetchPhaseEnvelope(long handle, PhaseEnvelopeTPRho* out) {
    long probe_length = 0, probe_components = 0;
    long errcode = 0;
    char msg[AS_ERR_BUFFER_LEN] = {};
    AbstractState_get_phase_envelope_data_checkedMemory(handle, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &probe_length,
                                                        &probe_components, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (probe_length <= 0) {
        // errcode==0 here means AS->get_phase_envelope_data() genuinely
        // returned zero points (build_phase_envelope() was never called for
        // this handle) -- not an error the underlying call raised itself.
        // A non-zero errcode with probe_length still 0 means the exception
        // happened before the length was even written (e.g. a dead handle),
        // so TranslateASError's usual BAD_HANDLE/LOWLEVEL_ERROR mapping applies.
        if (errcode) return TranslateASError(msg, 1);
        return MAKELRESULT(PHASE_ENVELOPE_NOT_BUILT, 1);
    }

    long actual_components = 0;
    long component_probe_length = 0;
    long component_errcode = 0;
    char component_msg[AS_ERR_BUFFER_LEN] = {};
    AbstractState_get_phase_envelope_data_checkedMemory(handle, probe_length, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                        &component_probe_length, &actual_components, &component_errcode, component_msg,
                                                        AS_ERR_BUFFER_LEN);
    if (actual_components <= 0) {
        // Should be unreachable in practice -- get_phase_envelope_data()
        // always returns one composition entry per fluid in the mixture,
        // and a handle always has at least one -- but if this probe somehow
        // didn't throw (a future CoolProp version tolerating
        // maxComponents==0?) or threw before *actual_components was even
        // written (e.g. the handle went dead between call 1 and here),
        // report whatever this probe actually raised rather than dividing
        // by a bogus zero component count below.
        if (component_errcode) return TranslateASError(component_msg, 1);
        return MAKELRESULT(LOWLEVEL_ERROR, 1);
    }

    out->T.resize(static_cast<size_t>(probe_length));
    out->P.resize(static_cast<size_t>(probe_length));
    out->rhomolar_vap.resize(static_cast<size_t>(probe_length));
    out->rhomolar_liq.resize(static_cast<size_t>(probe_length));
    std::vector<double> x(static_cast<size_t>(probe_length) * static_cast<size_t>(actual_components));
    std::vector<double> y(x.size());

    long final_length = 0, final_components = 0;
    errcode = 0;
    AbstractState_get_phase_envelope_data_checkedMemory(handle, probe_length, actual_components, out->T.data(), out->P.data(),
                                                        out->rhomolar_vap.data(), out->rhomolar_liq.data(), x.data(), y.data(), &final_length,
                                                        &final_components, &errcode, msg, AS_ERR_BUFFER_LEN);
    if (errcode) return TranslateASError(msg, 1);

    out->T.resize(static_cast<size_t>(final_length));
    out->P.resize(static_cast<size_t>(final_length));
    out->rhomolar_vap.resize(static_cast<size_t>(final_length));
    out->rhomolar_liq.resize(static_cast<size_t>(final_length));
    return 0;
}

// This code executes the user function CP_AS_get_phase_envelope_data, which
// returns the traced phase envelope (from a prior AS_build_phase_envelope
// call) as a table: one row per envelope point, columns T, P, rhomolar_vap,
// rhomolar_liq. Per-component compositions (x/y) are not surfaced by this
// function -- see FetchPhaseEnvelope()'s comment.
static LRESULT CP_AS_get_phase_envelope_data(LPCOMPLEXARRAY Data,       // output: N rows x {T, P, rhomolar_vap, rhomolar_liq}
                                             LPCCOMPLEXSCALAR Handle,   // AbstractState handle from AS_factory
                                             LPCCOMPLEXSCALAR Trigger)  // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    PhaseEnvelopeTPRho pe;
    r = FetchPhaseEnvelope(handle, &pe);
    if (r) return r;

    std::vector<std::vector<double>> Vec(pe.T.size(), std::vector<double>(4));
    for (size_t i = 0; i < pe.T.size(); ++i) {
        Vec[i][0] = pe.T[i];
        Vec[i][1] = pe.P[i];
        Vec[i][2] = pe.rhomolar_vap[i];
        Vec[i][3] = pe.rhomolar_liq[i];
    }
    return AllocateToMathcadArray(Data, Vec);
}

// This code executes the user function CP_AS_pe_tmax, the cricondentherm --
// the point on the phase envelope traced by a prior AS_build_phase_envelope
// call with the highest temperature. Returns a 2-element column vector
// [T; P]. Implemented as a max-scan over FetchPhaseEnvelope()'s T column,
// entirely on this side of the C API -- CoolProp itself tracks this same
// point internally (PhaseEnvelopeData::iTsat_max, set in
// PhaseEnvelopeRoutines::finalize(), src/Backends/Helmholtz/
// PhaseEnvelopeRoutines.cpp) but does not expose it through the public
// Low-Level C API this wrapper is built on, and extending that shared
// surface is out of scope here -- see this function's entry in
// MathcadWrappers.rst for what that means for exactness of the result.
static LRESULT CP_AS_pe_tmax(LPCOMPLEXARRAY Point,       // output: 2-element column vector [T; P]
                             LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                             LPCCOMPLEXSCALAR Trigger)   // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    PhaseEnvelopeTPRho pe;
    r = FetchPhaseEnvelope(handle, &pe);
    if (r) return r;

    size_t imax = static_cast<size_t>(std::max_element(pe.T.begin(), pe.T.end()) - pe.T.begin());
    std::vector<std::vector<double>> Vec = {{pe.T[imax]}, {pe.P[imax]}};
    return AllocateToMathcadArray(Point, Vec);
}

// This code executes the user function CP_AS_pe_pmax, the cricondenbar --
// the point on the phase envelope traced by a prior AS_build_phase_envelope
// call with the highest pressure. Returns a 2-element column vector [T; P].
// See CP_AS_pe_tmax()'s comment above (mirrors it exactly, scanning P
// instead of T; CoolProp's internal counterpart is
// PhaseEnvelopeData::ipsat_max).
static LRESULT CP_AS_pe_pmax(LPCOMPLEXARRAY Point,       // output: 2-element column vector [T; P]
                             LPCCOMPLEXSCALAR Handle,    // AbstractState handle from AS_factory
                             LPCCOMPLEXSCALAR Trigger)   // unused -- see CP_AS_mole_fractions_liquid()'s comment for why this argument exists
{
    std::scoped_lock lock(g_as_mutex);
    (void)Trigger;
    LRESULT r = CheckRealOrError(Handle, 1);
    if (r) return r;

    long handle;
    r = ToLongOrError(Handle, BAD_HANDLE, 1, &handle);
    if (r) return r;

    PhaseEnvelopeTPRho pe;
    r = FetchPhaseEnvelope(handle, &pe);
    if (r) return r;

    size_t imax = static_cast<size_t>(std::max_element(pe.P.begin(), pe.P.end()) - pe.P.begin());
    std::vector<std::vector<double>> Vec = {{pe.T[imax]}, {pe.P[imax]}};
    return AllocateToMathcadArray(Point, Vec);
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

FUNCTIONINFO ASBuildPhaseEnvelope = {
  const_cast<char*>("AS_build_phase_envelope"),                                                       // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Level"),                                                                 // Description of input parameters
  const_cast<char*>("Traces the phase envelope for a Low-Level state Handle; returns Handle"),        // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_build_phase_envelope,                                                            // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                      // Returns a Mathcad complex scalar (Handle, unchanged)
  2,                                                                                                    // Number of arguments
  {COMPLEX_SCALAR, MC_STRING}                                                                          // Argument types
};

FUNCTIONINFO ASGetPhaseEnvelopeData = {
  const_cast<char*>("AS_get_phase_envelope_data"),                                                                    // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                              // Description of input parameters
  const_cast<char*>("Returns the traced phase envelope as a table: T, P, rhomolar_vap, rhomolar_liq (one row per point)"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get_phase_envelope_data,                                                                        // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                                     // Returns a Mathcad complex array
  2,                                                                                                                  // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                                    // Argument types
};

FUNCTIONINFO ASPeTmax = {
  const_cast<char*>("AS_pe_tmax"),                                                                    // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                // Description of input parameters
  const_cast<char*>("Cricondentherm: [T; P] at the phase envelope's highest-temperature point"),      // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_pe_tmax,                                                                          // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                       // Returns a Mathcad complex array (2-element column vector)
  2,                                                                                                    // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                     // Argument types
};

FUNCTIONINFO ASPePmax = {
  const_cast<char*>("AS_pe_pmax"),                                                                    // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                // Description of input parameters
  const_cast<char*>("Cricondenbar: [T; P] at the phase envelope's highest-pressure point"),           // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_pe_pmax,                                                                          // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                       // Returns a Mathcad complex array (2-element column vector)
  2,                                                                                                    // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                     // Argument types
};

FUNCTIONINFO ASGetSatLiquid = {
  const_cast<char*>("AS_get_sat_liquid"),                                                                          // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Output Parameter Index"),                                                             // Description of input parameters
  const_cast<char*>("Returns one output parameter from the saturated LIQUID side of a Low-Level state Handle's current point"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get_sat_liquid,                                                                               // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                                  // Returns a Mathcad complex scalar
  2,                                                                                                                // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                                 // Argument types
};

FUNCTIONINFO ASGetSatVapor = {
  const_cast<char*>("AS_get_sat_vapor"),                                                                           // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Output Parameter Index"),                                                             // Description of input parameters
  const_cast<char*>("Returns one output parameter from the saturated VAPOR side of a Low-Level state Handle's current point"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get_sat_vapor,                                                                                // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                                  // Returns a Mathcad complex scalar
  2,                                                                                                                // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                                 // Argument types
};

FUNCTIONINFO ASMoleFractionsLiquid = {
  const_cast<char*>("AS_mole_fractions_liquid"),                                                        // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                  // Description of input parameters
  const_cast<char*>("Returns the saturated LIQUID side's mole fractions at a Low-Level state Handle's current point"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_mole_fractions_liquid,                                                              // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                         // Returns a Mathcad complex array (column vector)
  2,                                                                                                      // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                       // Argument types
};

FUNCTIONINFO ASMoleFractionsVapor = {
  const_cast<char*>("AS_mole_fractions_vapor"),                                                         // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                  // Description of input parameters
  const_cast<char*>("Returns the saturated VAPOR side's mole fractions at a Low-Level state Handle's current point"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_mole_fractions_vapor,                                                               // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                         // Returns a Mathcad complex array (column vector)
  2,                                                                                                      // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                       // Argument types
};

FUNCTIONINFO ASGenerateUpdatePair = {
  const_cast<char*>("AS_generate_update_pair"),                                                                                // Name by which Mathcad will recognize the function
  const_cast<char*>("Output Parameter Index 1, Output Parameter Index 2"),                                                     // Description of input parameters
  const_cast<char*>("Resolves two output parameters to the CoolProp input pair name they form, e.g. \"PT_INPUTS\""),           // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_generate_update_pair,                                                                                     // Pointer to the function code.
  MC_STRING,                                                                                                                   // Returns a Mathcad string
  2,                                                                                                                            // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                                             // Argument types
};

FUNCTIONINFO ASMoleToMassFractions = {
  const_cast<char*>("AS_mole_to_mass_fractions"),                                                       // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, MoleFractions"),                                                            // Description of input parameters
  const_cast<char*>("Converts a mole-fraction composition to the equivalent mass fractions for Handle's mixture"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_mole_to_mass_fractions,                                                             // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                         // Returns a Mathcad complex array (column vector)
  2,                                                                                                      // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_ARRAY}                                                                        // Argument types
};

FUNCTIONINFO ASMassToMoleFractions = {
  const_cast<char*>("AS_mass_to_mole_fractions"),                                                       // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, MassFractions"),                                                            // Description of input parameters
  const_cast<char*>("Converts a mass-fraction composition to the equivalent mole fractions for Handle's mixture"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_mass_to_mole_fractions,                                                             // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                         // Returns a Mathcad complex array (column vector)
  2,                                                                                                      // Number of arguments
  {COMPLEX_SCALAR, COMPLEX_ARRAY}                                                                        // Argument types
};

FUNCTIONINFO ASGetPhase = {
  const_cast<char*>("AS_get_phase"),                                                                    // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                  // Description of input parameters
  const_cast<char*>("Returns the phase name of a Low-Level state Handle's current point, e.g. \"phase_liquid\""),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get_phase,                                                                          // Pointer to the function code.
  MC_STRING,                                                                                              // Returns a Mathcad string
  2,                                                                                                      // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                       // Argument types
};

FUNCTIONINFO ASGetMoleFractions = {
  const_cast<char*>("AS_get_mole_fractions"),                                                            // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                   // Description of input parameters
  const_cast<char*>("Returns a Low-Level state Handle's current bulk mole fractions"),                    // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_get_mole_fractions,                                                                  // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                          // Returns a Mathcad complex array (column vector)
  2,                                                                                                       // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                        // Argument types
};

FUNCTIONINFO ASBackendName = {
  const_cast<char*>("AS_backend_name"),                                                                  // Name by which Mathcad will recognize the function
  const_cast<char*>("Handle, Trigger"),                                                                   // Description of input parameters
  const_cast<char*>("Returns the backend name (e.g. \"HEOS\") a Low-Level state Handle is using"),        // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_AS_backend_name,                                                                        // Pointer to the function code.
  MC_STRING,                                                                                               // Returns a Mathcad string
  2,                                                                                                       // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR, COMPLEX_SCALAR}                                                                        // Argument types
};

#endif  // MATHCAD_LOWLEVEL_H
