import { useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { fitFluid, parseTable } from "../lib/incompfit";

/**
 * Dialog for defining a new incompressible (pure) fluid from tabulated data.
 *
 * The user pastes a table (from a spreadsheet or as CSV) with a temperature
 * column and any of density / heat capacity / conductivity / viscosity, all
 * in SI units. Fitting happens live in the frontend (src/lib/incompfit.ts,
 * mirroring the offline Python pipeline); saving registers the fluid with
 * the CoolProp core via add_fluids_as_JSON("INCOMP", ...) and persists the
 * definition so it is re-registered on the next start.
 */

export const USER_FLUIDS_STORAGE_KEY = "coolprop.userIncompressibleFluids";

export function loadUserFluidDefinitions(): Record<string, string> {
  try {
    return JSON.parse(localStorage.getItem(USER_FLUIDS_STORAGE_KEY) ?? "{}");
  } catch {
    return {};
  }
}

/** Re-register every persisted user fluid with the core (called on startup).
 * Returns the names that were registered successfully. */
export async function registerUserFluids(): Promise<string[]> {
  const names: string[] = [];
  const defs = loadUserFluidDefinitions();
  for (const [name, json] of Object.entries(defs)) {
    try {
      await invoke("add_incompressible_fluid", { fluidJson: json });
      names.push(name);
    } catch (e) {
      console.error(`re-registering user fluid ${name} failed:`, e);
    }
  }
  return names;
}

const PLACEHOLDER = [
  "T\trho\tcp\tk\tmu",
  "290\t1001\t4181\t0.59\t0.0011",
  "300\t998\t4184\t0.61\t0.00089",
  "310\t994\t4187\t0.62\t0.00069",
  "…\t(SI units: K, kg/m³, J/kg/K, W/m/K, Pa·s)",
].join("\n");

const PROPERTY_LABELS: Record<string, string> = {
  density: "Density",
  specific_heat: "Heat capacity",
  conductivity: "Conductivity",
  viscosity: "Viscosity",
};

interface Props {
  onSaved: (name: string) => void;
  onCancel: () => void;
}

export default function AddFluidDialog({ onSaved, onCancel }: Props) {
  const [name, setName] = useState("");
  const [description, setDescription] = useState("");
  const [reference, setReference] = useState("");
  const [tableText, setTableText] = useState("");
  const [saveError, setSaveError] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);

  const preview = useMemo(() => {
    if (!tableText.trim() || !name.trim()) return null;
    try {
      const table = parseTable(tableText);
      return { ok: true as const, ...fitFluid({ name: name.trim(), description, reference, table }) };
    } catch (e) {
      return { ok: false as const, message: String(e instanceof Error ? e.message : e) };
    }
  }, [tableText, name, description, reference]);

  const handleSave = async () => {
    if (!preview || !preview.ok) return;
    setSaving(true);
    setSaveError(null);
    try {
      const json = JSON.stringify(preview.fluidJson);
      await invoke("add_incompressible_fluid", { fluidJson: json });
      const defs = loadUserFluidDefinitions();
      defs[name.trim()] = json;
      localStorage.setItem(USER_FLUIDS_STORAGE_KEY, JSON.stringify(defs));
      onSaved(name.trim());
    } catch (e) {
      setSaveError(String(e));
    } finally {
      setSaving(false);
    }
  };

  return (
    <div className="modal-overlay" onClick={(e) => { if (e.target === e.currentTarget) onCancel(); }}>
      <div className="modal-card" style={{ maxWidth: 560 }}>
        <div className="modal-title">Add Incompressible Fluid</div>

        <div className="modal-body">
          <div className="field-row">
            <div className="field-group" style={{ flex: 1 }}>
              <label>Name</label>
              <input
                type="text"
                value={name}
                placeholder="MyCoolant"
                onChange={(e) => setName(e.target.value)}
              />
            </div>
            <div className="field-group" style={{ flex: 2 }}>
              <label>Description (optional)</label>
              <input type="text" value={description} onChange={(e) => setDescription(e.target.value)} />
            </div>
          </div>

          <div className="field-group">
            <label>Data source / reference (optional)</label>
            <input
              type="text"
              value={reference}
              placeholder="manufacturer datasheet, publication, …"
              onChange={(e) => setReference(e.target.value)}
            />
          </div>

          <div className="field-group">
            <label>
              Property table — paste from a spreadsheet. Columns: T plus any of
              rho, cp, k, mu (SI units); a header row is recognised, blank cells
              mean “no data”. Density and cp are required.
            </label>
            <textarea
              value={tableText}
              placeholder={PLACEHOLDER}
              rows={9}
              spellCheck={false}
              style={{ fontFamily: "monospace", fontSize: 12, whiteSpace: "pre" }}
              onChange={(e) => setTableText(e.target.value)}
            />
          </div>

          {preview && !preview.ok && <div className="error-msg">{preview.message}</div>}

          {preview && preview.ok && (
            <div className="field-group">
              <label>Fit preview</label>
              <table className="results-table">
                <thead>
                  <tr><th>Property</th><th>Points</th><th>Degree</th><th>RMS (rel.)</th></tr>
                </thead>
                <tbody>
                  {Object.entries(preview.report).map(([prop, r]) => (
                    <tr key={prop}>
                      <td>{PROPERTY_LABELS[prop] ?? prop}</td>
                      <td className="num">{r.points}</td>
                      <td className="num">{r.degree}</td>
                      <td className="num">{(r.nrms * 100).toPrecision(2)}%</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {preview.warnings.length > 0 && (
                <div style={{ fontSize: 11, color: "#a60" }}>
                  {preview.warnings.map((w, i) => <div key={i}>⚠ {w}</div>)}
                </div>
              )}
              <div style={{ fontSize: 11, color: "#666" }}>
                Enthalpy and entropy are derived exactly from the density and
                heat-capacity fits. Saturation pressure is unavailable (no
                vapour-pressure data), so the fluid is treated as liquid at
                any pressure.
              </div>
            </div>
          )}

          {saveError && <div className="error-msg">{saveError}</div>}
        </div>

        <div className="modal-footer">
          <button className="btn-secondary" onClick={onCancel}>Cancel</button>
          <button className="primary" onClick={handleSave} disabled={!preview || !preview.ok || saving}>
            {saving ? "Saving…" : "Save fluid"}
          </button>
        </div>
      </div>
    </div>
  );
}
