#ifndef COOLPROP_FACTORY_OPTIONS_H
#define COOLPROP_FACTORY_OPTIONS_H

#include <string>

namespace CoolProp {

// Result of splitting a factory string of the form
//
//     <backend-and-fluid>[?<options>]
//
// into its two halves.  `options_json` is the *raw* options payload as a
// JSON string, populated by either an inline `?{...}` blob or a
// `?@<path>` indirection that reads the file at <path>.  No JSON parsing
// or schema validation happens here — that lives one layer up in the
// schema validator.  When the factory string has no `?` suffix or the
// suffix is empty, `options_json` is the empty string and downstream
// callers treat that as "no options".
struct FactoryOptionsParseResult
{
    std::string clean_string;  ///< Factory string with any `?<suffix>` stripped.
    std::string options_json;  ///< Raw options JSON (empty if none supplied).
};

// Split a factory string on its first `?` and resolve the suffix.
//
//   - No `?`                         → options_json = ""
//   - `?` then empty / whitespace    → options_json = ""
//   - `?@<path>`                     → options_json = contents of <path>
//                                      (FileIOError on read failure)
//   - anything else                  → options_json = the verbatim tail
//
// The parser deliberately does *not* validate JSON, look at the schema,
// or interpret the suffix in any other way.  Any `?` characters after
// the first one (e.g. inside a JSON string value, inside a URL, inside
// the `@<path>` filename itself) are kept verbatim — the contract is
// `find_first_of('?')` exactly once.
FactoryOptionsParseResult parse_factory_options(const std::string& factory_string);

}  // namespace CoolProp

#endif  // COOLPROP_FACTORY_OPTIONS_H
