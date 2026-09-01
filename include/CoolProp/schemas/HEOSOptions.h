#ifndef COOLPROP_SCHEMAS_HEOS_OPTIONS_H
#define COOLPROP_SCHEMAS_HEOS_OPTIONS_H

// JSON Schema for the HEOS backend's factory-string options blob.
//
// Compiled into the binary as a string literal so the validator has no runtime file dependency.
// See Web/coolprop/BackendOptions.rst for the grammar and the wider mechanism.
//
// Strict-mode: `additionalProperties: false` at every level rejects unknown keys at factory time,
// so a typo surfaces immediately rather than silently defaulting.

namespace CoolProp {

inline constexpr const char kHEOSOptionsSchemaJson[] = R"JSON({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "title": "HEOS backend options",
  "type": "object",
  "additionalProperties": false,
  "properties": {
    "schema": {
      "type": "integer",
      "const": 1,
      "description": "Schema version. Bump when the layout below changes."
    },
    "RES": {
      "type": "object",
      "additionalProperties": false,
      "description": "Residual entropy scaling transport models. Opt-in per property: HEOS already has reference correlations for most fluids, and RES must never silently displace one.",
      "properties": {
        "viscosity": {
          "type": "boolean",
          "description": "Use the RES viscosity model (Martinek 2025) instead of the fluid's reference correlation. Throws at evaluation time if the fluid carries no RES viscosity parameters."
        },
        "conductivity": {
          "type": "boolean",
          "description": "Use the RES thermal-conductivity model (Li 2024) instead of the fluid's reference correlation. Throws at evaluation time if the fluid carries no RES conductivity parameters."
        },
        "mixture_critical_enhancement": {
          "type": "boolean",
          "description": "Apply the Olchowy-Sengers critical enhancement to MIXTURES as well as pure fluids. Off by default: it needs the mixture critical point, which HEOS has to solve for -- slowly, and not always successfully -- and the physical case for a critical enhancement in mixtures is not well established. Pure fluids always get the enhancement and are unaffected by this."
        }
      }
    }
  }
})JSON";

}  // namespace CoolProp

#endif  // COOLPROP_SCHEMAS_HEOS_OPTIONS_H
