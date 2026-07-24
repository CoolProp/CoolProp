#ifndef COOLPROP_SCHEMAS_TABULAR_OPTIONS_H
#define COOLPROP_SCHEMAS_TABULAR_OPTIONS_H

// JSON Schema for the tabular backends' (BICUBIC / TTSE) factory-string
// options blob.
//
// Compiled into the binary as a string literal so the validator has no
// runtime file dependency.  Matches
// docs/superpowers/specs/2026-05-16-backend-options-string-design.md.
//
// Strict-mode: `additionalProperties: false` at every level rejects
// unknown keys at factory time so typos surface immediately rather
// than silently defaulting.

namespace CoolProp {

inline constexpr const char kTabularOptionsSchemaJson[] = R"JSON({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "title": "Tabular (BICUBIC/TTSE) backend options",
  "type": "object",
  "additionalProperties": false,
  "properties": {
    "schema": {
      "type": "integer",
      "const": 1,
      "description": "Schema version. Bump when the layout below changes."
    },
    "grid": {
      "type": "object",
      "additionalProperties": false,
      "required": ["Nx", "Ny"],
      "description": "Per-instance single-phase grid sizing. Overrides the TABULAR_NX / TABULAR_NY configuration globals for this instance only. Both axes must be given together so the instance's resolution never depends on the globals' later values.",
      "properties": {
        "Nx": {"type": "integer", "minimum": 2, "maximum": 100000,
                "description": "Number of grid points along the x axis (enthalpy for the PH table, temperature for the PT table)."},
        "Ny": {"type": "integer", "minimum": 2, "maximum": 100000,
                "description": "Number of grid points along the y axis (log-spaced pressure for both tables)."}
      }
    }
  }
})JSON";

}  // namespace CoolProp

#endif  // COOLPROP_SCHEMAS_TABULAR_OPTIONS_H
