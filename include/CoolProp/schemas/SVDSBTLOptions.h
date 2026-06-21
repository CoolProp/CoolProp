#ifndef COOLPROP_SCHEMAS_SVDSBTL_OPTIONS_H
#define COOLPROP_SCHEMAS_SVDSBTL_OPTIONS_H

// JSON Schema for the SVDSBTL backend's factory-string options blob.
//
// Compiled into the binary as a string literal so the validator has no
// runtime file dependency.  Matches
// docs/superpowers/specs/2026-05-16-backend-options-string-design.md.
//
// Strict-mode: `additionalProperties: false` at every level rejects
// unknown keys at factory time so typos surface immediately rather
// than silently defaulting.

namespace CoolProp {

inline constexpr const char kSVDSBTLOptionsSchemaJson[] = R"JSON({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "title": "SVDSBTL backend options",
  "type": "object",
  "additionalProperties": false,
  "properties": {
    "schema": {
      "type": "integer",
      "const": 1,
      "description": "Schema version. Bump when the layout below changes."
    },
    "prebuild": {
      "type": "boolean",
      "description": "When true, eagerly build every supported input-pair surface (PT, HmassP, DmassT, PSmass) at construction instead of lazy-loading the secondary pairs on first query. Opt-in complement to the default lazy loading: materializes the whole fluid up front (docs / benchmarking / warm-cache pre-fill) and turns a build/env failure into a loud construction-time error instead of a silent blank later."
    },
    "pmin": {
      "type": "number",
      "exclusiveMinimum": 0,
      "description": "Lower pressure bound (absolute Pa) for the subcritical PT / HmassP / PSmass surfaces. Defaults to the fluid's triple-point pressure. Must be >= p_triple: the subcritical regions are bounded by the liquid-vapour saturation curve, which does not exist below the triple line. The DmassT surface is temperature-indexed and unaffected."
    },
    "grid": {
      "type": "object",
      "additionalProperties": false,
      "description": "Per-region SVD grid sizing.",
      "properties": {
        "NT": {"type": "integer", "minimum": 2, "maximum": 100000,
                "description": "Number of points along the secondary (non-log) axis."},
        "NR": {"type": "integer", "minimum": 2, "maximum": 100000,
                "description": "Number of points along the primary (log-p) axis."},
        "rank": {"type": "integer", "minimum": 1, "maximum": 10000,
                  "description": "SVD truncation rank."}
      }
    },
    "properties": {
      "type": "object",
      "additionalProperties": false,
      "description": "Per-output-property knobs.",
      "properties": {
        "transport": {
          "type": "string",
          "enum": ["auto", "on", "off"],
          "description": "Whether to tabulate viscosity / conductivity. 'auto' probes the source backend for transport support."
        }
      }
    },
    "critical_patch": {
      "type": "object",
      "additionalProperties": false,
      "description": "HEOS-fallback patch covering the rank-truncation gap near the critical point.",
      "properties": {
        "mode": {
          "type": "string",
          "enum": ["auto", "off", "fixed"],
          "description": "How the patch bbox is chosen: 'auto' = build-time calibration; 'off' = no patch; 'fixed' = use the explicit bbox below."
        },
        "source": {
          "type": ["string", "null"],
          "enum": [null, "HEOS", "REFPROP", "IF97"],
          "description": "Backend used inside the patch. null = same as the SVDSBTL truth source."
        },
        "tolerance": {
          "type": "number",
          "exclusiveMinimum": 0,
          "description": "Target residual for auto-calibration."
        },
        "metric": {
          "type": "string",
          "enum": ["D", "A", "H"],
          "description": "Property whose residual drives auto-calibration."
        },
        "bbox": {
          "type": ["array", "null"],
          "minItems": 4,
          "maxItems": 4,
          "items": {"type": "number"},
          "description": "[Tlo, Thi, plo, phi] as multipliers of (Tcrit, pcrit). Only honoured when mode == 'fixed'."
        }
      }
    }
  }
})JSON";

}  // namespace CoolProp

#endif  // COOLPROP_SCHEMAS_SVDSBTL_OPTIONS_H
