// MathcadConfig.h : CoolProp global Configuration functions for the Mathcad
// wrapper.
//
// This is an included implementation fragment, not a standalone header -- it
// is meant to be #include'd from exactly one place, partway through
// CoolPropMathcad.cpp, after that file has already set up:
//   - mcadincl.h (and the MC_STRING/STRING substitution trick above it)
//   - CoolProp/CoolProp.h and CoolProp/Configuration.h (for
//     configuration_keys, config_string_to_key,
//     get_config_bool/int/double/string, set_config_bool/int/double/string,
//     CoolProp::set_warning_string())
//   - enum EC and CPErrorMessageTable (for BAD_CONFIG_KEY,
//     RESTRICTED_CONFIG_KEY, BAD_CONFIG_TYPE, BAD_CONFIG_BOOL_VALUE,
//     BAD_CONFIG_INT_VALUE, BAD_CONFIG_DOUBLE_VALUE,
//     BAD_CONFIG_RESTRICTED_VALUE, MUST_BE_REAL, MAKELRESULT)
//   - the general Mathcad wrapper helpers: CheckRealOrError, AllocMathcadString
// It is kept separate from CoolPropMathcad.cpp (and from MathcadLowLevel.h)
// purely to keep each file from growing unbounded as more functions are
// added -- see MathcadLowLevel.h's own top comment for the same rationale.
// Deliberately independent of MathcadLowLevel.h: configuration access applies
// equally to the high-level (PropsSI/HAPropsSI/...) and Low-Level (AS_*)
// surface -- it is not itself part of either -- so these functions are NOT
// prefixed with "AS_" and this file has no dependency on MathcadStateGuard,
// g_as_mutex, or anything else scoped to AbstractState handles (this file
// has its own dedicated mutex, g_config_set_mutex, below).
//
// These wrap CoolProp's global Configuration API
// (https://coolprop.org/coolprop/Configuration.html), which stores one
// process-wide value per configuration_keys entry, typed at registration
// time as exactly one of bool/int/double/string (see
// include/CoolProp/detail/configuration_keys.h's CONFIGURATION_KEYS_ENUM
// X-macro -- the single source of truth for every key's name, type, and
// default). Mathcad has no tagged/variant type, so -- unlike a single
// generic "config_get"/"config_set" -- this is 4 type-specific getter/setter
// pairs, each Mathcad-facing function committing to exactly one
// COMPLEX_SCALAR-or-MC_STRING return type as Mathcad's FUNCTIONINFO
// registration requires:
//
//   config_get_bool(Key)               -> 1 or 0 (not Mathcad "true"/"false")
//   config_set_bool(Key, Value)        -> "Set" or "Not Set" -- see below
//   config_get_int(Key)                -> integer-as-double
//   config_set_int(Key, Value)         -> "Set" or "Not Set" -- see below
//   config_get_double(Key)             -> double
//   config_set_double(Key, Value)      -> "Set" or "Not Set" -- see below
//   config_get_string(Key)             -> string
//   config_set_string(Key, Value)      -> "Set" or "Not Set" -- see below
//
// Calling the getter/setter for the WRONG type on a given key (e.g.
// config_get_bool("TABULAR_NX"), an int-valued key) is a Custom Error
// (BAD_CONFIG_TYPE), not a silent misread -- CoolProp's own
// ConfigurationItem::check_data_type() already refuses this with a
// same-shaped "type does not match" exception; this just gives it a
// specific Mathcad error code instead of falling through to a generic one.
//
// EACH KEY CAN ONLY BE SET ONCE PER MATHCAD PRIME SESSION -- BY DESIGN, NOT
// A LIMITATION TO WORK AROUND.
// An earlier revision of this file let config_set_*() mutate a key on every
// call. Manual testing in Mathcad Prime surfaced a correctness hazard
// specific to this wrapper's use of CoolProp's global Configuration API:
// every configuration key is ONE process-wide value, shared by every
// worksheet open in that Mathcad Prime instance -- and by every worksheet
// subsequently opened in it -- until that Mathcad process exits and the DLL
// is unloaded (separate Mathcad Prime instances on the same machine are
// unaffected; each loads its own copy of this DLL with its own independent
// CoolProp Configuration singleton). Combined with Mathcad's dependency-
// graph recalculation order (region/dependency order, not top-to-bottom
// source order, and not guaranteed stable across recalculations or across
// which worksheets happen to be open), an unconstrained setter call in one
// worksheet could silently change results in a completely unrelated open
// worksheet, in an order that isn't reliably reproducible from the
// worksheet's own contents.
//
// The fix here doesn't try to make an unconstrained write safe -- there is
// still no handle to scope a write to (unlike AS_factory's handles), so
// nothing in this wrapper could make arbitrary repeated mutation safe
// without fighting Mathcad's own evaluation model. Instead, each key gets
// AT MOST ONE effective write per session: the FIRST successful
// config_set_*() call for a given key (from any worksheet, whichever one
// Mathcad happens to evaluate first) wins and actually applies. That
// decision is latched process-wide via g_config_locked_keys below. Every
// config_get_*()/config_set_*() call in this file holds g_config_set_mutex
// for its ENTIRE body, not just the g_config_locked_keys lookup -- CoolProp's
// Configuration has no locking of its own, so a lock scoped only around the
// claim registry would leave the actual CoolProp::get_config_*()/
// set_config_*() calls unsynchronized against each other, and a losing
// call's compare-read could then run concurrently with the winning call's
// write to that same key: a real data race on Configuration's internal
// storage, not just a stale-read risk. See g_config_set_mutex's own comment
// below for the rest of this rationale. Every later call for that same key:
//   - If it passes the SAME value the winning call used: treated as a
//     harmless re-affirmation, not an error -- this is what makes
//     recalculating the same config_set_*() program-block statement safe
//     (Mathcad re-executes a program block's statements on every
//     recalculation, so a worksheet's own call to a key it already won
//     must keep succeeding). Returns "Set", same as the original call.
//   - If it passes a DIFFERENT value: the configuration is NOT changed --
//     the winning value stays in effect. Returns the Mathcad string
//     "Not Set" (a normal return value, not a Custom Error, so it doesn't
//     interrupt worksheet evaluation with a red error box) AND calls
//     CoolProp::set_warning_string() with a message identifying the key,
//     the value that was ignored, and the value actually in effect --
//     readable via get_global_param_string("warnstring") (CoolProp.h;
//     read-then-cleared, thread-safe, the same shared "outbox" PropsSI's
//     own error path uses -- see the comment on error_string/warning_string
//     in src/CoolProp.cpp).
//
// This does NOT make the underlying ordering deterministic: which
// worksheet's call happens to run first still depends on Mathcad's
// recalculation order, not on anything visible in either worksheet's own
// contents. What it changes is that a conflicting later call is now
// detectable (a different return string, plus an inspectable warning)
// instead of silently overwriting a value another worksheet was relying on.
// A worksheet that discards config_set_*()'s return value without checking
// it will still not SEE that it lost -- checking the return, or checking
// get_global_param_string("warnstring"), is what surfaces a conflict.
//
// Two keys are refused by every config_set_*() function unconditionally,
// with a hard Custom Error, regardless of the once-per-session logic above
// (but remain readable via the matching config_get_*()): FLOAT_PUNCTUATION
// and LIST_STRING_DELIMITER. Both are relied on by this wrapper's OWN
// string parsing -- FLOAT_PUNCTUATION controls the decimal separator
// CoolProp uses when formatting/parsing numbers in strings;
// LIST_STRING_DELIMITER is the separator GetComponentMolarMasses()
// (MathcadLowLevel.h) already assumes when splitting a handle's "&"-joined
// fluid-name list via CoolProp::get_config_string(LIST_STRING_DELIMITER).
// Changing either at runtime -- even just once -- would silently corrupt
// string parsing elsewhere in this same DLL, not just whatever the caller
// intended, so these two stay read-only from Mathcad no matter what.
//
// One int key only accepts a small, documented set of legal values rather
// than any finite in-range integer: MIXTURE_STABILITY_ALGORITHM ("0:
// legacy, 1: Michelsen (default)", per its own description in
// configuration_keys.h). See IsRestrictedIntKey()/IsLegalRestrictedIntValue()
// below -- this validation runs before the once-per-session claim logic, so
// an invalid Value is always rejected with a Custom Error regardless of
// whether the key has already been claimed this session.
//
// Also supported, and unaffected by all of the above: setting a
// configuration key via a Windows environment variable named "COOLPROP_"
// followed by the key name (e.g. COOLPROP_MIXTURE_STABILITY_ALGORITHM=1)
// before launching Mathcad Prime -- read exactly once by CoolProp's own
// Configuration::possibly_set_from_env() when the singleton is first
// constructed (include/CoolProp/Configuration.h), before any worksheet's
// calculations run. Doing it this way means every worksheet in the session
// starts from that value already applied, with nothing left for
// config_set_*() to contend over for that key -- the FIRST config_set_*()
// call to run against an env-set key still only "claims" it if that call's
// Value matches, or diverges from, whatever is already active, following
// the exact same first-call/compare rules described above.
//
// None of the four getters (nor the four setters) take a Trigger argument,
// unlike several AS_* functions in MathcadLowLevel.h that read handle state
// Mathcad's dependency graph can't otherwise see changing (see
// CP_AS_mole_fractions_liquid()'s comment there for that rationale). A
// config_get_*() call for a key that's already been claimed this session
// only ever sees that one fixed value for the rest of the session, so
// there's no analogous "changed after this Handle was created" case for a
// Trigger to cover here.
//
// A ninth function, get_config_as_json_string(Trigger), exists alongside
// the eight above for viewing the WHOLE configuration at once -- see the
// comment further down, just above CP_get_config_as_json_string(), for its
// own design. It keeps a Trigger argument purely because it has no other
// argument to satisfy Mathcad's one-argument-minimum with (config_get_*()
// above still has Key for that).

#include <limits>  // std::numeric_limits<int>, used by TryRoundToInt()
#include <mutex>   // std::mutex/std::scoped_lock, used by g_config_set_mutex below -- this
                   // file is declared independent of MathcadLowLevel.h above, so it can't
                   // rely on that file's own #include <mutex> having already run first
#include <set>     // g_config_locked_keys, below

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

// The once-per-session claim registry described in this file's top
// comment: which configuration_keys have already had a successful
// config_set_*() call in this Mathcad Prime session.
//
// g_config_set_mutex guards EVERY CoolProp::Configuration access this file
// makes: config_get_*()'s read, and config_set_*()'s whole claim-decide-
// then-write-or-compare sequence, each held for the call's full duration --
// not just the g_config_locked_keys lookup/mutation. Deliberately separate
// from MathcadLowLevel.h's g_as_mutex (that one is scoped to AbstractState
// handle operations; this one has nothing to do with handles).
//
// Why the whole call, not just the registry: CoolProp's Configuration/
// ConfigurationItem (include/CoolProp/Configuration.h) has no locking of
// its own. A lock scoped only around g_config_locked_keys's insert/erase
// (as an earlier revision of this file did, via separate TryClaimConfigKey/
// ReleaseConfigKeyClaim helpers, each taking and releasing the lock on
// their own) leaves a real window: thread A wins the claim and starts
// writing via CoolProp::set_config_*(); thread B, in the same moment,
// loses the claim and immediately reads back via CoolProp::get_config_*()
// to decide "same value" vs "different value" -- both calls running fully
// unsynchronized against each other, on the SAME key's storage. That's not
// just a stale-read risk: for a string-valued key, it's a concurrent
// unsynchronized std::string assignment racing a read of that same
// std::string, which is undefined behavior, not merely a logic bug.
// Holding one lock across each call's entire body closes this for every
// caller in this file; a std::set::insert()'s return value (whether the
// insertion actually happened) still doubles as the compare-and-set
// primitive deciding winner vs. loser, now just evaluated inside that
// wider critical section instead of its own narrow one.
//
// Trade-off (matching g_as_mutex's own documented one in
// MathcadLowLevel.h): this briefly serializes calls even against UNRELATED
// keys, and even plain config_get_*() reads against any in-flight
// config_set_*() call. Simplicity/correctness over throughput is the right
// call here -- these are quick, infrequent calls, not a hot loop like
// AS_props_multi.
//
// A claim that turns out not to succeed (CoolProp::set_config_*() itself
// throws, e.g. BAD_CONFIG_TYPE) is released by erasing the key from
// g_config_locked_keys before returning, still under the same lock, so a
// call that never actually took effect doesn't permanently lock the key
// out for a later, valid call. This assumes CoolProp::set_config_*() only
// throws BEFORE mutating any state, never partway through -- true today
// for all four (ConfigurationItem::set_bool/set_integer/set_double/
// set_string all call check_data_type(), which throws on a mismatch,
// before touching the stored value; set_config_string()'s post-write
// force_unload_REFPROP() step, src/Configuration.cpp, swallows its own
// errors rather than throwing). If a future change ever made one of these
// throw AFTER a partial write, releasing the claim here would let a later
// call re-attempt against an already-mutated key -- worth re-checking if
// that assumption ever changes.
static std::set<configuration_keys> g_config_locked_keys;
static std::mutex g_config_set_mutex;

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
    {
        // Locked for the same reason every config_set_*() call is -- see
        // g_config_set_mutex's own comment above.
        std::scoped_lock lock(g_config_set_mutex);
        try {
            value = CoolProp::get_config_bool(key);
        } catch (...) {
            return MAKELRESULT(BAD_CONFIG_TYPE, 1);
        }
    }

    Result->real = value ? 1.0 : 0.0;
    Result->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_bool, a wrapper for
// CoolProp::set_config_bool() -- see this file's top comment for the
// once-per-session claim/compare design this and the other three setters
// share.
static LRESULT CP_config_set_bool(LPMCSTRING Dummy,        // output: "Set" or "Not Set"
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
    const bool newValue = (Value->real != 0.0);

    std::string resultText;
    {
        // Held for the whole claim-decide-then-write-or-compare sequence --
        // see g_config_set_mutex's own comment above for why this can't be
        // narrowed to just the g_config_locked_keys lookup.
        std::scoped_lock lock(g_config_set_mutex);
        if (g_config_locked_keys.insert(key).second) {
            try {
                CoolProp::set_config_bool(key, newValue);
            } catch (...) {
                g_config_locked_keys.erase(key);
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            resultText = "Set";
        } else {
            bool current;
            try {
                current = CoolProp::get_config_bool(key);
            } catch (...) {
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            if (current == newValue) {
                resultText = "Set";
            } else {
                CoolProp::set_warning_string(format("config_set_bool(\"%s\", %s) ignored -- already set to %s earlier in this Mathcad Prime "
                                                    "session; configuration keys can only be set once per session",
                                                    Key->str, newValue ? "1" : "0", current ? "1" : "0"));
                resultText = "Not Set";
            }
        }
    }

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
    {
        std::scoped_lock lock(g_config_set_mutex);
        try {
            value = CoolProp::get_config_int(key);
        } catch (...) {
            return MAKELRESULT(BAD_CONFIG_TYPE, 1);
        }
    }

    Result->real = static_cast<double>(value);
    Result->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_int, a wrapper for
// CoolProp::set_config_int() -- see CP_config_set_bool()'s comment above
// for the shared once-per-session claim/compare design.
static LRESULT CP_config_set_int(LPMCSTRING Dummy,        // output: "Set" or "Not Set"
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
    int newValue;
    if (!TryRoundToInt(Value->real, &newValue)) {
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
        if (Value->real != static_cast<double>(newValue)) {
            return MAKELRESULT(BAD_CONFIG_INT_VALUE, 2);
        }
        if (!IsLegalRestrictedIntValue(key, newValue)) {
            return MAKELRESULT(BAD_CONFIG_RESTRICTED_VALUE, 2);
        }
    }

    std::string resultText;
    {
        std::scoped_lock lock(g_config_set_mutex);
        if (g_config_locked_keys.insert(key).second) {
            try {
                CoolProp::set_config_int(key, newValue);
            } catch (...) {
                g_config_locked_keys.erase(key);
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            resultText = "Set";
        } else {
            int current;
            try {
                current = CoolProp::get_config_int(key);
            } catch (...) {
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            if (current == newValue) {
                resultText = "Set";
            } else {
                CoolProp::set_warning_string(format("config_set_int(\"%s\", %d) ignored -- already set to %d earlier in this Mathcad Prime "
                                                    "session; configuration keys can only be set once per session",
                                                    Key->str, newValue, current));
                resultText = "Not Set";
            }
        }
    }

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
    {
        std::scoped_lock lock(g_config_set_mutex);
        try {
            value = CoolProp::get_config_double(key);
        } catch (...) {
            return MAKELRESULT(BAD_CONFIG_TYPE, 1);
        }
    }

    Result->real = value;
    Result->imag = 0;

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_double, a wrapper for
// CoolProp::set_config_double() -- see CP_config_set_bool()'s comment above
// for the shared once-per-session claim/compare design.
static LRESULT CP_config_set_double(LPMCSTRING Dummy,        // output: "Set" or "Not Set"
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
    const double newValue = Value->real;

    std::string resultText;
    {
        std::scoped_lock lock(g_config_set_mutex);
        if (g_config_locked_keys.insert(key).second) {
            try {
                CoolProp::set_config_double(key, newValue);
            } catch (...) {
                g_config_locked_keys.erase(key);
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            resultText = "Set";
        } else {
            double current;
            try {
                current = CoolProp::get_config_double(key);
            } catch (...) {
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            // Exact equality, not a tolerance compare: safe for the
            // documented usage (a literal Value in the config_set_double
            // call parses to a bit-identical double on every recalculation,
            // so the intended "same call, run again" case always compares
            // equal). The gap this doesn't cover is a computed Value
            // expression whose upstream calculation isn't guaranteed
            // bit-reproducible across recalculations -- a 1-ULP drift there
            // could spuriously report "Not Set" for what the caller
            // considers a no-op re-run. Not addressed here since it
            // requires picking an arbitrary tolerance; worth revisiting if
            // it turns out to matter in practice.
            if (current == newValue) {
                resultText = "Set";
            } else {
                // %.17g, not %g: %g's default 6-significant-digit rounding
                // can print two genuinely different doubles identically
                // (e.g. 1000000 vs 1000001 both show as "1e+06"), which
                // would make this warning useless for telling them apart.
                // %.17g is enough digits to round-trip any double exactly.
                CoolProp::set_warning_string(format("config_set_double(\"%s\", %.17g) ignored -- already set to %.17g earlier in this Mathcad "
                                                    "Prime session; configuration keys can only be set once per session",
                                                    Key->str, newValue, current));
                resultText = "Not Set";
            }
        }
    }

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
    {
        std::scoped_lock lock(g_config_set_mutex);
        try {
            value = CoolProp::get_config_string(key);
        } catch (...) {
            return MAKELRESULT(BAD_CONFIG_TYPE, 1);
        }
    }

    Result->str = AllocMathcadString(value);

    // normal return
    return 0;
}

// This code executes the user function CP_config_set_string, a wrapper for
// CoolProp::set_config_string() -- see CP_config_set_bool()'s comment above
// for the shared once-per-session claim/compare design.
// Setting ALTERNATIVE_REFPROP_PATH/_HMX_BNC_PATH/_LIBRARY_PATH additionally
// forces REFPROP to unload (CoolProp::force_unload_REFPROP(), inside
// set_config_string() itself, src/Configuration.cpp) so the next REFPROP
// call re-loads from the new path -- not something this wrapper needs to
// special-case, it's already handled underneath, and only happens on the
// one claiming call, not on every re-affirming recalculation.
static LRESULT CP_config_set_string(LPMCSTRING Dummy,   // output: "Set" or "Not Set"
                                    LPCMCSTRING Key,    // configuration key name
                                    LPCMCSTRING Value)  // the string to set
{
    configuration_keys key;
    LRESULT r = ResolveConfigKey(Key->str, 1, &key);
    if (r) return r;
    r = CheckConfigKeyNotRestricted(key, 1);
    if (r) return r;

    const std::string newValue(Value->str);

    std::string resultText;
    {
        std::scoped_lock lock(g_config_set_mutex);
        if (g_config_locked_keys.insert(key).second) {
            try {
                CoolProp::set_config_string(key, newValue);
            } catch (...) {
                g_config_locked_keys.erase(key);
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            resultText = "Set";
        } else {
            std::string current;
            try {
                current = CoolProp::get_config_string(key);
            } catch (...) {
                return MAKELRESULT(BAD_CONFIG_TYPE, 1);
            }
            if (current == newValue) {
                resultText = "Set";
            } else {
                CoolProp::set_warning_string(format("config_set_string(\"%s\", \"%s\") ignored -- already set to \"%s\" earlier in this Mathcad "
                                                    "Prime session; configuration keys can only be set once per session",
                                                    Key->str, newValue.c_str(), current.c_str()));
                resultText = "Not Set";
            }
        }
    }

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
// Returns the raw JSON with no reformatting. Attempts to manually insert
// line breaks into the string failed because:
// 1. Mathcad doesn't handle CR & LF characters well, and
// 2. The interface strips Unicode line breaks, Char(133) from the return string.
// Any parsing/formatting of this string for readability has to be performed
// in Mathcad after the string is returned.
static LRESULT CP_get_config_as_json_string(LPMCSTRING Result, LPCCOMPLEXSCALAR Trigger) {
    (void)Trigger;

    std::string json;
    {
        // Locked for the same reason every other function in this file is --
        // see g_config_set_mutex's own comment above.
        std::scoped_lock lock(g_config_set_mutex);
        try {
            json = CoolProp::get_config_as_json_string();
        } catch (...) {
            return MAKELRESULT(UNKNOWN, 1);
        }
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
    "Sets a boolean CoolProp configuration value once per Mathcad session (Value must be 1 or 0); returns \"Set\", "
    "or \"Not Set\" if Key was already set to a different value earlier this session"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_bool,                                                       // Pointer to the function code.
  MC_STRING,                                                                             // Returns a Mathcad string ("Set" or "Not Set")
  2,                                                                                     // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                            // Argument types
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
    "Sets an integer CoolProp configuration value once per Mathcad session (rounded to the nearest int); returns "
    "\"Set\", or \"Not Set\" if Key was already set to a different value earlier this session"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_int,                                                                 // Pointer to the function code.
  MC_STRING,                                                                                      // Returns a Mathcad string ("Set" or "Not Set")
  2,                                                                                              // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                                     // Argument types
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
    "Sets a double-valued CoolProp configuration value once per Mathcad session (must be finite); returns \"Set\", "
    "or \"Not Set\" if Key was already set to a different value earlier this session"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_double,                                                     // Pointer to the function code.
  MC_STRING,                                                                             // Returns a Mathcad string ("Set" or "Not Set")
  2,                                                                                     // Number of arguments
  {MC_STRING, COMPLEX_SCALAR}                                                            // Argument types
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
    "Sets a string-valued CoolProp configuration value once per Mathcad session; returns \"Set\", or \"Not Set\" if "
    "Key was already set to a different value earlier this session"),  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_config_set_string,                                   // Pointer to the function code.
  MC_STRING,                                                           // Returns a Mathcad string ("Set" or "Not Set")
  2,                                                                   // Number of arguments
  {MC_STRING, MC_STRING}                                               // Argument types
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
