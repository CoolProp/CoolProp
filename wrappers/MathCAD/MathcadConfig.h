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
//     BAD_CONFIG_INT_VALUE, BAD_CONFIG_DOUBLE_VALUE, MUST_BE_REAL,
//     MAKELRESULT)
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
//   config_get_bool(Key, Trigger)     -> 1 or 0 (not Mathcad "true"/"false")
//   config_set_bool(Key, Value)       -> dummy 0 on success; Value must be
//                                         exactly 0 or 1
//   config_get_int(Key, Trigger)      -> integer-as-double
//   config_set_int(Key, Value)        -> dummy 0 on success; Value is
//                                         rounded to the nearest int, same
//                                         finite+range checked conversion
//                                         MathcadLowLevel.h's
//                                         TryRoundToLong() uses for `long`
//   config_get_double(Key, Trigger)   -> double
//   config_set_double(Key, Value)     -> dummy 0 on success; Value must be
//                                         finite (not NaN/Infinity)
//   config_get_string(Key, Trigger)   -> string
//   config_set_string(Key, Value)     -> dummy 0 on success
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
// `Trigger` on the four getters exists for the same reason it does on
// several AS_* functions in MathcadLowLevel.h (see
// CP_AS_mole_fractions_liquid()'s comment there for the full rationale):
// Key's own value never changes between recalculations, so a getter whose
// only argument is Key gives Mathcad's dependency graph nothing to key a
// recalculation on when a DIFFERENT cell's config_set_*() call mutates the
// same process-wide value out from under it. Wire Trigger to something that
// actually changes when you need the getter to re-run -- e.g. the setter
// call's own dummy return value -- or use Recalculate Worksheet. The four
// setters don't need a Trigger: each already takes Key/Value as real
// arguments, and the normal case (editing Key or Value) already gives
// Mathcad a natural recalculation edge; sequencing a setter relative to
// unrelated reads still follows the same two authoring patterns documented
// for AS_factory() (a Mathcad program block, or Recalculate Worksheet).

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
static LRESULT CP_config_get_bool(LPCOMPLEXSCALAR Result,    // output: 1 if true, 0 if false
                                  LPCMCSTRING Key,           // configuration key name, e.g. "NORMALIZE_GAS_CONSTANTS"
                                  LPCCOMPLEXSCALAR Trigger)  // unused -- see this file's top comment for why this argument exists
{
    (void)Trigger;
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
static LRESULT CP_config_set_bool(LPCOMPLEXSCALAR Dummy,   // output: dummy value (0) on success
                                  LPCMCSTRING Key,         // configuration key name
                                  LPCCOMPLEXSCALAR Value)  // 1 (true) or 0 (false) -- no other value is accepted
{
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

    Dummy->real = 0;
    Dummy->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_int, which is a
// wrapper for CoolProp::get_config_int() -- see this file's top comment.
static LRESULT CP_config_get_int(LPCOMPLEXSCALAR Result,    // output: integer value, as a real scalar
                                 LPCMCSTRING Key,           // configuration key name, e.g. "TABULAR_NX"
                                 LPCCOMPLEXSCALAR Trigger)  // unused -- see this file's top comment for why this argument exists
{
    (void)Trigger;
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

// This code executes the user function CP_config_set_int, which is a
// wrapper for CoolProp::set_config_int() -- see this file's top comment.
static LRESULT CP_config_set_int(LPCOMPLEXSCALAR Dummy,   // output: dummy value (0) on success
                                 LPCMCSTRING Key,         // configuration key name
                                 LPCCOMPLEXSCALAR Value)  // integer value (rounded to the nearest int)
{
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

    try {
        CoolProp::set_config_int(key, intValue);
    } catch (...) {
        return MAKELRESULT(BAD_CONFIG_TYPE, 1);
    }

    Dummy->real = 0;
    Dummy->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_double, which is a
// wrapper for CoolProp::get_config_double() -- see this file's top comment.
static LRESULT CP_config_get_double(LPCOMPLEXSCALAR Result,    // output: the configuration value
                                    LPCMCSTRING Key,           // configuration key name, e.g. "R_U_CODATA"
                                    LPCCOMPLEXSCALAR Trigger)  // unused -- see this file's top comment for why this argument exists
{
    (void)Trigger;
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
static LRESULT CP_config_set_double(LPCOMPLEXSCALAR Dummy,   // output: dummy value (0) on success
                                    LPCMCSTRING Key,         // configuration key name
                                    LPCCOMPLEXSCALAR Value)  // the value to set -- must be finite
{
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

    Dummy->real = 0;
    Dummy->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_get_string, which is a
// wrapper for CoolProp::get_config_string() -- see this file's top comment.
static LRESULT CP_config_get_string(LPMCSTRING Result,         // output: the configuration value
                                    LPCMCSTRING Key,           // configuration key name, e.g. "ALTERNATIVE_REFPROP_PATH"
                                    LPCCOMPLEXSCALAR Trigger)  // unused -- see this file's top comment for why this argument exists
{
    (void)Trigger;
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
static LRESULT CP_config_set_string(LPCOMPLEXSCALAR Dummy,  // output: dummy value (0) on success
                                    LPCMCSTRING Key,        // configuration key name
                                    LPCMCSTRING Value)      // the string to set
{
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

    Dummy->real = 0;
    Dummy->imag = 0;

    // normal return
    return 0;
}

// ********************************************************************************************************
// FUNCTIONINFO structs for the eight functions above.
// ********************************************************************************************************

FUNCTIONINFO ConfigGetBool = {
  const_cast<char*>("config_get_bool"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Trigger"),     // Description of input parameters
  const_cast<char*>(
    "Returns a boolean CoolProp configuration value (1 or 0) for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_bool,                                                      // Pointer to the function code.
  COMPLEX_SCALAR,                                                                       // Returns a Mathcad complex scalar
  2,                           // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {MC_STRING, COMPLEX_SCALAR}  // Argument types
};

FUNCTIONINFO ConfigSetBool = {
  const_cast<char*>("config_set_bool"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),       // Description of input parameters
  const_cast<char*>(
    "Sets a boolean CoolProp configuration value (Value must be 1 or 0); returns a dummy 0"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_bool,                                                             // Pointer to the function code.
  COMPLEX_SCALAR,                                                                              // Returns a Mathcad complex scalar (dummy)
  2,                                                                                           // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                  // Argument types
};

FUNCTIONINFO ConfigGetInt = {
  const_cast<char*>("config_get_int"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Trigger"),    // Description of input parameters
  const_cast<char*>(
    "Returns an integer CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_int,                                               // Pointer to the function code.
  COMPLEX_SCALAR,                                                               // Returns a Mathcad complex scalar
  2,                                                                            // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {MC_STRING, COMPLEX_SCALAR}                                                   // Argument types
};

FUNCTIONINFO ConfigSetInt = {
  const_cast<char*>("config_set_int"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),      // Description of input parameters
  const_cast<char*>(
    "Sets an integer CoolProp configuration value (rounded to the nearest int); returns a dummy 0"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_int,                                                                     // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                     // Returns a Mathcad complex scalar (dummy)
  2,                                                                                                  // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                         // Argument types
};

FUNCTIONINFO ConfigGetDouble = {
  const_cast<char*>("config_get_double"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Trigger"),       // Description of input parameters
  const_cast<char*>(
    "Returns a double-valued CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_double,                                                 // Pointer to the function code.
  COMPLEX_SCALAR,                                                                    // Returns a Mathcad complex scalar
  2,                                                                                 // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {MC_STRING, COMPLEX_SCALAR}                                                        // Argument types
};

FUNCTIONINFO ConfigSetDouble = {
  const_cast<char*>("config_set_double"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),         // Description of input parameters
  const_cast<char*>(
    "Sets a double-valued CoolProp configuration value (must be finite); returns a dummy 0"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_double,                                                           // Pointer to the function code.
  COMPLEX_SCALAR,                                                                              // Returns a Mathcad complex scalar (dummy)
  2,                                                                                           // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                  // Argument types
};

FUNCTIONINFO ConfigGetString = {
  const_cast<char*>("config_get_string"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Trigger"),       // Description of input parameters
  const_cast<char*>(
    "Returns a string-valued CoolProp configuration value for the given key name"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_get_string,                                                 // Pointer to the function code.
  MC_STRING,                                                                         // Returns a Mathcad string
  2,                                                                                 // Number of arguments (Mathcad requires >= 1; Trigger is unused)
  {MC_STRING, COMPLEX_SCALAR}                                                        // Argument types
};

FUNCTIONINFO ConfigSetString = {
  const_cast<char*>("config_set_string"),  // Name by which Mathcad will recognize the function
  const_cast<char*>("Key, Value"),         // Description of input parameters
  const_cast<char*>(
    "Sets a string-valued CoolProp configuration value; returns a dummy 0"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_string,                                          // Pointer to the function code.
  COMPLEX_SCALAR,                                                             // Returns a Mathcad complex scalar (dummy)
  2,                                                                          // Number of arguments
  {MC_STRING, MC_STRING}                                                      // Argument types
};
