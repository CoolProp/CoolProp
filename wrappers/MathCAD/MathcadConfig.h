// MathcadConfig.h : CoolProp global Configuration get/set functions for the
// Mathcad wrapper.
//
// This is an included implementation fragment, not a standalone header -- it
// is meant to be #include'd from exactly one place, partway through
// CoolPropMathcad.cpp, after that file has already set up:
//   - mcadincl.h (and the MC_STRING/STRING substitution trick above it)
//   - CoolProp/CoolProp.h and CoolProp/Configuration.h (for
//     configuration_keys, config_string_to_key,
//     get_config_bool/int/double/string, set_config_bool/int/double/string)
//   - enum EC and CPErrorMessageTable (for BAD_CONFIG_KEY,
//     RESTRICTED_CONFIG_KEY, BAD_CONFIG_TYPE, BAD_CONFIG_BOOL_VALUE,
//     BAD_CONFIG_INT_VALUE, BAD_CONFIG_DOUBLE_VALUE,
//     BAD_CONFIG_RESTRICTED_VALUE, MUST_BE_REAL, MAKELRESULT)
//   - the general Mathcad wrapper helpers: CheckRealOrError, AllocMathcadString
// It is kept separate from CoolPropMathcad.cpp (and from MathcadLowLevel.h)
// purely to keep each file from growing unbounded as more functions are
// added -- see MathcadLowLevel.h's own top comment for the same rationale.
// Deliberately independent of MathcadLowLevel.h: configuration get/set
// applies equally to the high-level (PropsSI/HAPropsSI/...) and Low-Level
// (AS_*) surface -- it is not itself part of either -- so these functions
// are NOT prefixed with "AS_" and this file has no dependency on
// MathcadStateGuard, g_as_mutex, or anything else scoped to AbstractState
// handles.
//
// These wrap CoolProp's global Configuration API
// (https://coolprop.org/coolprop/Configuration.html), which stores one
// process-wide value per configuration_keys entry, typed at registration
// time as exactly one of bool/int/double/string (see
// include/CoolProp/detail/configuration_keys.h's CONFIGURATION_KEYS_ENUM
// X-macro -- the single source of truth for every key's name, type, and
// default). Mathcad has no tagged/variant type, so -- unlike a single
// generic "config_get"/"config_set" -- this is 4 type-specific getter/
// setter pairs, each Mathcad-facing function committing to exactly one
// COMPLEX_SCALAR-or-MC_STRING return type as Mathcad's FUNCTIONINFO
// registration requires:
//
//   config_get_bool(Key)               -> 1 or 0 (not Mathcad "true"/"false")
//   config_set_bool(Key, Value)        -> "Set" on success; Value must be
//                                          exactly 0 or 1
//   config_get_int(Key)                -> integer-as-double
//   config_set_int(Key, Value)         -> "Set" on success; Value is
//                                          rounded to the nearest int, same
//                                          finite+range checked conversion
//                                          MathcadLowLevel.h's
//                                          TryRoundToLong() uses for `long`
//                                          -- and, for the handful of int
//                                          keys with a documented small
//                                          legal set (e.g.
//                                          MIXTURE_STABILITY_ALGORITHM: 0
//                                          or 1 only), must exactly match
//                                          one of them -- see
//                                          IsLegalRestrictedIntValue()
//   config_get_double(Key)             -> double
//   config_set_double(Key, Value)      -> "Set" on success; Value must be
//                                          finite (not NaN/Infinity)
//   config_get_string(Key)             -> string
//   config_set_string(Key, Value)      -> "Set" on success
//
// Every config_set_*() returns the Mathcad string "Set" on success (never
// a dummy number) -- each function's own `resultText` local is initialized
// to "Fail" up front and reassigned to "Set" only immediately before the
// true success return, so a future code path that somehow reached the
// final AllocMathcadString() call without actually succeeding would surface
// as an honest "Fail" string rather than a stale/uninitialized value; every
// actual failure path returns a Custom Error before reaching that point at
// all, so "Fail" is never really seen in practice.
//
// Calling the getter/setter for the WRONG type on a given key (e.g.
// config_get_bool("TABULAR_NX"), an int-valued key) is a Custom Error
// (BAD_CONFIG_TYPE), not a silent misread -- CoolProp's own
// ConfigurationItem::check_data_type() already refuses this with a
// same-shaped "type does not match" exception; this just gives it a
// specific Mathcad error code instead of falling through to a generic one.
//
// Two keys are refused by every config_set_*() function (but remain
// readable via the matching config_get_*()): FLOAT_PUNCTUATION and
// LIST_STRING_DELIMITER. Both are relied on by this wrapper's OWN string
// parsing -- FLOAT_PUNCTUATION controls the decimal separator CoolProp uses
// when formatting/parsing numbers in strings; LIST_STRING_DELIMITER is the
// separator GetComponentMolarMasses() (MathcadLowLevel.h) already assumes
// when splitting a handle's "&"-joined fluid-name list via
// CoolProp::get_config_string(LIST_STRING_DELIMITER). Changing either at
// runtime would silently corrupt string parsing elsewhere in this same DLL,
// not just whatever the caller intended -- so these two are read-only from
// Mathcad, surfaced as RESTRICTED_CONFIG_KEY rather than allowed to quietly
// break unrelated functions.
//
// None of the eight take a Trigger argument, unlike several AS_* functions
// in MathcadLowLevel.h that read handle state Mathcad's dependency graph
// can't otherwise see changing (see CP_AS_mole_fractions_liquid()'s
// comment there for that rationale). Deliberately different here: these
// eight functions read/write a single process-wide value with NO handle to
// scope it to, so relying on Mathcad's dependency graph to sequence
// config_get_*() correctly after some earlier config_set_*() is fragile by
// construction -- Mathcad recalculates by region/dependency order, not
// top-to-bottom source order, so which one "wins" for a same-key get/set
// pair with no explicit dependency between them is exactly the kind of
// out-of-order surprise a Trigger could paper over in one specific case
// without fixing the general problem. Rather than lean on that, these
// functions are meant to be called sequentially -- e.g. every
// config_set_*() call for a worksheet's desired configuration placed
// together near the top, in one Mathcad program block, or as ordinary
// worksheet regions followed by an explicit **Recalculate Worksheet**
// (Ctrl-F5 or Ctrl-F9) before anything reads the result -- the same
// "known, deliberate order" authoring discipline AS_factory()'s own
// worksheet-level pattern documents in MathcadLowLevel.h, just without a
// handle to chain through. A Trigger argument would only invite relying on
// automatic recalculation for something these functions are specifically
// not safe to use that way for.
//
// A ninth function, get_config_as_json_string(Trigger), exists alongside
// the eight above for viewing the WHOLE configuration at once -- see the
// comment further down, just above CP_get_config_as_json_string(), for its
// own design. It keeps a Trigger argument purely because it has no other
// argument to satisfy Mathcad's one-argument-minimum with (config_get_*()
// above still has Key for that) -- not because it participates in the
// dependency graph any more safely than the eight above do; the same
// sequential-call, Recalculate-Worksheet discipline applies to it too.

#include <limits>  // std::numeric_limits<int>, used by TryRoundToInt()

// Helper: round a plain double to the nearest integer and return it as an
// int, validating BOTH that it's finite AND that the rounded value is
// representable in an int. Mirrors MathcadLowLevel.h's TryRoundToLong() --
// see that function's own comment for the full undefined-behavior
// rationale (std::llround() on a non-finite value, or narrowing an
// in-range-for-long-long-but-out-of-range-for-`int` result, is UB that in
// practice could alias a different int value rather than erroring) --
// reimplemented here (targeting `int` instead of `long`) rather than
// shared across the two files, since this file is deliberately independent
// of MathcadLowLevel.h (see this file's own top comment).
static inline bool TryRoundToInt(double real, int* out) {
    if (!std::isfinite(real)) return false;
    // Reject anything outside `int`'s range (with 0.5 of slack either side
    // for correct rounding at the boundary) BEFORE calling std::llround() --
    // see TryRoundToLong()'s matching comment in MathcadLowLevel.h for why
    // this has to happen before the call, not just checked on its result.
    constexpr double kIntMin = static_cast<double>((std::numeric_limits<int>::min)());
    constexpr double kIntMax = static_cast<double>((std::numeric_limits<int>::max)());
    if (real < kIntMin - 0.5 || real > kIntMax + 0.5) return false;
    const long long rounded = std::llround(real);
    // NOT dead code -- see TryRoundToLong()'s matching comment: at the
    // exact half-integer boundary the pre-check above can still let
    // round-half-away-from-zero produce kIntMax+1/kIntMin-1, which only
    // this check catches.
    if (rounded < static_cast<long long>((std::numeric_limits<int>::min)()) || rounded > static_cast<long long>((std::numeric_limits<int>::max)())) {
        return false;
    }
    *out = static_cast<int>(rounded);
    return true;
}

// Helper: resolve a Mathcad-supplied configuration key NAME to CoolProp's
// enum, shared by all eight config_get_*/config_set_* functions below.
// CoolProp::config_string_to_key() throws a message-less
// CoolProp::ValueError() for an unrecognized name (verified in
// src/Configuration.cpp: every recognized name is matched by the
// X-macro-generated if-chain, with a bare throw as the final fallthrough)
// -- proactively caught here and translated to BAD_CONFIG_KEY rather than
// left to surface as an unhelpful UNKNOWN.
static LRESULT ResolveConfigKey(const char* name, unsigned int position, configuration_keys* out) {
    try {
        *out = CoolProp::config_string_to_key(name);
        return 0;
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_KEY, position);
    }
}

// Helper: reject the two configuration keys that would break THIS
// wrapper's own assumptions if changed from their defaults -- shared by
// all four config_set_*() functions below. config_get_*() never needs
// this: reading either key back is harmless, only mutating one is the
// problem -- see this file's own top comment for the full rationale.
static LRESULT CheckConfigKeyNotRestricted(configuration_keys key, unsigned int position) {
    if (key == FLOAT_PUNCTUATION || key == LIST_STRING_DELIMITER) {
        return MAKELRESULT(RESTRICTED_CONFIG_KEY, position);
    }
    return 0;
}

// This code executes the user function CP_config_get_bool, which is a
// wrapper for CoolProp::get_config_bool() -- see this file's top comment
// for the overall config_get_*/config_set_* design.
static LRESULT CP_config_get_bool(LPCOMPLEXSCALAR Result,  // output: 1 if true, 0 if false
                                  LPCMCSTRING Key)         // configuration key name, e.g. "NORMALIZE_GAS_CONSTANTS"
{
    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;

    bool value;
    try {
        value = CoolProp::get_config_bool(key);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    Result->real = value ? 1.0 : 0.0;
    Result->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_bool, which is a
// wrapper for CoolProp::set_config_bool() -- see this file's top comment.
static LRESULT CP_config_set_bool(LPMCSTRING Dummy,        // output: "Set" on success
                                  LPCMCSTRING Key,         // configuration key name
                                  LPCCOMPLEXSCALAR Value)  // 1 (true) or 0 (false) -- no other value is accepted
{
    std::string resultText = "Fail";  // overwritten to "Set" only at the true success point below

    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;
    r = CheckConfigKeyNotRestricted(key, 1);
    if (r) return r;

    r = CheckRealOrError(Value, 2);
    if (r) return r;
    if (Value->real != 0.0 && Value->real != 1.0) {
        return MAKELRESULT(BAD_CONFIG_BOOL_VALUE, 2);
    }

    try {
        CoolProp::set_config_bool(key, Value->real != 0.0);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    resultText = "Set";
    Dummy->str = AllocMathcadString(resultText);

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_int, which is a
// wrapper for CoolProp::get_config_int() -- see this file's top comment.
static LRESULT CP_config_get_int(LPCOMPLEXSCALAR Result,  // output: integer value, as a real scalar
                                 LPCMCSTRING Key)         // configuration key name, e.g. "TABULAR_NX"
{
    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;

    int value;
    try {
        value = CoolProp::get_config_int(key);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    Result->real = static_cast<double>(value);
    Result->imag = 0;

    // normal return
    return 0;
}

// Helper: does `key` have a small, documented set of legal integer values
// rather than accepting any finite in-range integer? Currently only
// MIXTURE_STABILITY_ALGORITHM ("0: legacy, 1: Michelsen (default)", per its
// own description in configuration_keys.h) -- a single named key handled
// directly here rather than a general per-key legal-value registry, since
// that's the simplest thing that solves the actual problem for exactly one
// key today. If more restricted int keys are added later, this function
// (and IsLegalRestrictedIntValue() below) is the place to extend, and
// BAD_CONFIG_RESTRICTED_VALUE's single fixed message (CPErrorMessageTable
// in CoolPropMathcad.cpp) would need to become key-specific at that point
// too -- deliberately not built out now for a set of one.
//
// REFPROP_ERROR_THRESHOLD was deliberately considered and left OUT of this
// restricted set: unlike MIXTURE_STABILITY_ALGORITHM, it isn't a small
// enumerated choice -- it's a threshold compared against REFPROP's own
// `ierr` output across many different internal Fortran subroutines
// (src/Backends/REFPROP/REFPROPMixtureBackend.cpp), where the SIGN carries
// the primary meaning (negative = warning-only, positive = hard error) and
// specific magnitudes are assigned per-subroutine by REFPROP itself, not
// enumerated anywhere in CoolProp's own source or public documentation.
static bool IsRestrictedIntKey(configuration_keys key) {
    return key == MIXTURE_STABILITY_ALGORITHM;
}

// Helper: is `intValue` one of `key`'s legal values? Only meaningful when
// IsRestrictedIntKey(key) is true; CP_config_set_int() below only calls
// this after that check has already passed.
static bool IsLegalRestrictedIntValue(configuration_keys key, int intValue) {
    if (key == MIXTURE_STABILITY_ALGORITHM) {
        return intValue == 0 || intValue == 1;
    }
    return true;
}

// This code executes the user function CP_config_set_int, which is a
// wrapper for CoolProp::set_config_int() -- see this file's top comment.
static LRESULT CP_config_set_int(LPMCSTRING Dummy,        // output: "Set" on success
                                 LPCMCSTRING Key,         // configuration key name
                                 LPCCOMPLEXSCALAR Value)  // integer value (rounded to the nearest int)
{
    std::string resultText = "Fail";  // overwritten to "Set" only at the true success point below

    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;
    r = CheckConfigKeyNotRestricted(key, 1);
    if (r) return r;

    r = CheckRealOrError(Value, 2);
    if (r) return r;
    int intValue;
    if (!TryRoundToInt(Value->real, &intValue)) {
        return MAKELRESULT(BAD_CONFIG_INT_VALUE, 2);
    }

    if (IsRestrictedIntKey(key)) {
        // A restricted key's Value must be an EXACT whole number, not just
        // close enough to round -- e.g. 0.5 would otherwise round to 1 (a
        // legal MIXTURE_STABILITY_ALGORITHM choice) and silently succeed
        // without the caller having actually typed a whole number. Checked
        // BEFORE the legal-value-set check below, so a non-integer Value
        // always surfaces the more fundamental "not a whole number" error
        // (BAD_CONFIG_INT_VALUE) rather than "not a legal choice"
        // (BAD_CONFIG_RESTRICTED_VALUE) -- unrestricted int keys are
        // unaffected either way, since this whole branch only runs for a
        // key IsRestrictedIntKey() already said yes to.
        if (Value->real != static_cast<double>(intValue)) {
            return MAKELRESULT(BAD_CONFIG_INT_VALUE, 2);
        }
        if (!IsLegalRestrictedIntValue(key, intValue)) {
            return MAKELRESULT(BAD_CONFIG_RESTRICTED_VALUE, 2);
        }
    }

    try {
        CoolProp::set_config_int(key, intValue);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    resultText = "Set";
    Dummy->str = AllocMathcadString(resultText);

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_double, which is a
// wrapper for CoolProp::get_config_double() -- see this file's top comment.
static LRESULT CP_config_get_double(LPCOMPLEXSCALAR Result,  // output: the configuration value
                                    LPCMCSTRING Key)         // configuration key name, e.g. "R_U_CODATA"
{
    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;

    double value;
    try {
        value = CoolProp::get_config_double(key);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    Result->real = value;
    Result->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_double, which is a
// wrapper for CoolProp::set_config_double() -- see this file's top comment.
static LRESULT CP_config_set_double(LPMCSTRING Dummy,        // output: "Set" on success
                                    LPCMCSTRING Key,         // configuration key name
                                    LPCCOMPLEXSCALAR Value)  // the value to set -- must be finite
{
    std::string resultText = "Fail";  // overwritten to "Set" only at the true success point below

    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;
    r = CheckConfigKeyNotRestricted(key, 1);
    if (r) return r;

    r = CheckRealOrError(Value, 2);
    if (r) return r;
    if (!std::isfinite(Value->real)) {
        return MAKELRESULT(BAD_CONFIG_DOUBLE_VALUE, 2);
    }

    try {
        CoolProp::set_config_double(key, Value->real);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    resultText = "Set";
    Dummy->str = AllocMathcadString(resultText);

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_string, which is a
// wrapper for CoolProp::get_config_string() -- see this file's top comment.
static LRESULT CP_config_get_string(LPMCSTRING Result,  // output: the configuration value
                                    LPCMCSTRING Key)    // configuration key name, e.g. "ALTERNATIVE_REFPROP_PATH"
{
    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;

    std::string value;
    try {
        value = CoolProp::get_config_string(key);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    Result->str = AllocMathcadString(value);

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_string, which is a
// wrapper for CoolProp::set_config_string() -- see this file's top comment.
// Setting ALTERNATIVE_REFPROP_PATH/_HMX_BNC_PATH/_LIBRARY_PATH additionally
// forces REFPROP to unload (CoolProp::force_unload_REFPROP(), inside
// set_config_string() itself, src/Configuration.cpp) so the next REFPROP
// call re-loads from the new path -- not something this wrapper needs to
// special-case, it's already handled underneath.
static LRESULT CP_config_set_string(LPMCSTRING Dummy,   // output: "Set" on success
                                    LPCMCSTRING Key,    // configuration key name
                                    LPCMCSTRING Value)  // the string to set
{
    std::string resultText = "Fail";  // overwritten to "Set" only at the true success point below

    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;
    r = CheckConfigKeyNotRestricted(key, 1);
    if (r) return r;

    try {
        CoolProp::set_config_string(key, std::string(Value->str));
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    resultText = "Set";
    Dummy->str = AllocMathcadString(resultText);

    // normal return
    return 0;
}

// This code executes the user function get_config_as_json_string, which is
// a verbatim wrapper for CoolProp::get_config_as_json_string() -- the
// entire CoolProp Configuration (every key in
// include/CoolProp/detail/configuration_keys.h) as one JSON object string,
// e.g. {"NORMALIZE_GAS_CONSTANTS":true,"TABULAR_NX":200,...}. Named to
// match that C++/Python function directly (the same name
// CoolProp.CoolProp.get_config_as_json_string() already uses in Python)
// rather than this file's own "config_..." convention -- there's no
// per-key selection here to distinguish it from, so reusing the name
// already familiar from the other bindings is the least surprising choice.
//
// Returns the raw JSON with no reformatting -- an earlier version of this
// function tried inserting a Tab after each key, but the values didn't end
// up aligned anyway (key length varies too much for a single Tab stop to
// line up more than a few of them), so that was dropped as not actually
// worth the extra code. For a readable, one-pair-per-line view, use a
// Mathcad program block to split this string on ',' and rejoin with a
// line-break character built via vec2str() -- see the worked example
// linked in this function's entry in MathcadWrappers.rst.
static LRESULT CP_get_config_as_json_string(LPMCSTRING Result, LPCCOMPLEXSCALAR Trigger) {
    (void)Trigger;

    std::string json;
    try {
        json = CoolProp::get_config_as_json_string();
    } catch (...) {
        return MAKELRESULT(UNKNOWN, 1);
    }

    Result->str = AllocMathcadString(json);

    // normal return
    return 0;
}

// ********************************************************************************************************
// FUNCTIONINFO structs for the nine functions above.
// ********************************************************************************************************

FUNCTIONINFO ConfigGetBool = {
  const_cast<char*>("config_get_bool"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key"),              // Description of input parameters
  const_cast<char*>(
    "Returns a boolean CoolProp configuration value (1 or 0) for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_bool,                                                      // Pointer to the function code.
  COMPLEX_SCALAR,                                                                       // Returns a Mathcad complex scalar
  1,                                                                                    // Number of arguments
  {MC_STRING}                                                                           // Argument types
};

FUNCTIONINFO ConfigSetBool = {
  const_cast<char*>("config_set_bool"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),       // Description of input parameters
  const_cast<char*>(
    "Sets a boolean CoolProp configuration value (Value must be 1 or 0); returns \"Set\" on success"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_bool,                                                                      // Pointer to the function code.
  MC_STRING,                                                                                            // Returns a Mathcad string ("Set")
  2,                                                                                                    // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                           // Argument types
};

FUNCTIONINFO ConfigGetInt = {
  const_cast<char*>("config_get_int"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key"),             // Description of input parameters
  const_cast<char*>(
    "Returns an integer CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_int,                                               // Pointer to the function code.
  COMPLEX_SCALAR,                                                               // Returns a Mathcad complex scalar
  1,                                                                            // Number of arguments
  {MC_STRING}                                                                   // Argument types
};

FUNCTIONINFO ConfigSetInt = {
  const_cast<char*>("config_set_int"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),      // Description of input parameters
  const_cast<char*>(
    "Sets an integer CoolProp configuration value (rounded to the nearest int); returns \"Set\" on success"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_int,                                                                              // Pointer to the function code.
  MC_STRING,                                                                                                   // Returns a Mathcad string ("Set")
  2,                                                                                                           // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                                  // Argument types
};

FUNCTIONINFO ConfigGetDouble = {
  const_cast<char*>("config_get_double"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key"),                // Description of input parameters
  const_cast<char*>(
    "Returns a double-valued CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_double,                                                 // Pointer to the function code.
  COMPLEX_SCALAR,                                                                    // Returns a Mathcad complex scalar
  1,                                                                                 // Number of arguments
  {MC_STRING}                                                                        // Argument types
};

FUNCTIONINFO ConfigSetDouble = {
  const_cast<char*>("config_set_double"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),         // Description of input parameters
  const_cast<char*>(
    "Sets a double-valued CoolProp configuration value (must be finite); returns \"Set\" on success"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_double,                                                                    // Pointer to the function code.
  MC_STRING,                                                                                            // Returns a Mathcad string ("Set")
  2,                                                                                                    // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                           // Argument types
};

FUNCTIONINFO ConfigGetString = {
  const_cast<char*>("config_get_string"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key"),                // Description of input parameters
  const_cast<char*>(
    "Returns a string-valued CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_string,                                                 // Pointer to the function code.
  MC_STRING,                                                                         // Returns a Mathcad string
  1,                                                                                 // Number of arguments
  {MC_STRING}                                                                        // Argument types
};

FUNCTIONINFO ConfigSetString = {
  const_cast<char*>("config_set_string"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),         // Description of input parameters
  const_cast<char*>(
    "Sets a string-valued CoolProp configuration value; returns \"Set\" on success"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_string,                                                   // Pointer to the function code.
  MC_STRING,                                                                           // Returns a Mathcad string ("Set")
  2,                                                                                   // Number of arguments
  {MC_STRING, MC_STRING}                                                               // Argument types
};

FUNCTIONINFO GetConfigAsJsonString = {
  const_cast<char*>("get_config_as_json_string"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Trigger"),                    // Description of input parameters (unused -- see CP_get_config_as_json_string()'s comment)
  const_cast<char*>(
    "Returns the entire CoolProp configuration as one raw JSON string"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_config_as_json_string,                              // Pointer to the function code.
  MC_STRING,                                                              // Returns a Mathcad string
  1,                                                                      // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {COMPLEX_SCALAR}                                                        // Argument types
};
