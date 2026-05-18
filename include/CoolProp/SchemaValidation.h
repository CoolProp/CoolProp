#ifndef COOLPROP_SCHEMA_VALIDATION_H
#define COOLPROP_SCHEMA_VALIDATION_H

#include <string>

#if !defined(SWIG)
#    include "rapidjson_include.h"
#endif

namespace CoolProp {

#if !defined(SWIG)

// Strict-mode JSON Schema validation.
//
// Validates `instance` against `schema` using RapidJSON's built-in
// validator.  On a validation failure, throws CoolProp::ValueError with
// a message that quotes both the schema-pointer and the instance-pointer
// of the failing assertion (e.g. "/critical_patch/tolerance does not
// satisfy schema /properties/critical_patch/properties/tolerance: type
// expected number, got string").
//
// The schema is consumed as a parsed Document so multiple validations
// against the same schema can share parse cost.  Caller owns both
// documents.
void validate_against_schema(const rapidjson::Value& instance, const rapidjson::Document& schema);

// Convenience overload — parse `schema_json` once and validate.  Useful
// for the per-backend entry points where the schema is a constant
// string literal compiled into the binary.  Throws ValueError if
// schema_json itself is invalid JSON.
void validate_against_schema(const rapidjson::Value& instance, const std::string& schema_json);

// Produce a deterministic string representation of `value` for hashing
// and reproducibility.  Object keys are sorted recursively; arrays
// preserve order; numbers use RapidJSON's default formatting (already
// deterministic within a single build).  Strings are UTF-8 verbatim
// (no NFC normalisation performed — see note).
//
// This is *not* strict RFC 8785 (JCS): number normalisation skips
// the ECMAScript JSON.stringify rules, and string content isn't
// NFC-normalised.  It's "canonical within a CoolProp process" — the
// same logical-options-value produces the same bytes inside one
// binary, which is what cache hashing needs.  The form is not part
// of any external on-wire format, so this scope is sufficient.
//
// Returns the canonical string (UTF-8, no trailing newline).
std::string to_canonical_json(const rapidjson::Value& value);

// Convenience overload — parse `json_str` and canonicalise.  Throws
// ValueError if `json_str` is invalid JSON.
std::string to_canonical_json(const std::string& json_str);

#endif  // !SWIG

// Schema-validation helpers that work on JSON strings — exposed even
// to SWIG since they take only std::string.  Internally they parse,
// validate, and (for to_canonical_json_str) re-serialise.
//
// Mirrors the typed API above for callers that don't want to touch
// rapidjson directly (e.g. wrappers, tests).
void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json);

std::string to_canonical_json_str(const std::string& json_str);

}  // namespace CoolProp

#endif  // COOLPROP_SCHEMA_VALIDATION_H
