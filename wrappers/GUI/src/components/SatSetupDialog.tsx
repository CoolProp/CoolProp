import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { Basis } from "../App";

export interface SatConfig {
  fluid: string;
  backend: string;
  byTemp: boolean;
  minVal: number;
  maxVal: number;
  nPoints: number;
}

interface FluidLimits {
  t_min: number;
  t_max: number;
  p_min: number;
  p_max: number;
}

interface Props {
  fluids: string[];
  basis: Basis;
  initial?: Partial<SatConfig>;
  onConfirm: (config: SatConfig) => void;
  onCancel?: () => void;
}

type StepMode = "n_points" | "step";

export default function SatSetupDialog({ fluids, basis, initial, onConfirm, onCancel }: Props) {
  const [fluid,   setFluid]   = useState(initial?.fluid   ?? "Water");
  const [byTemp,  setByTemp]  = useState(initial?.byTemp  ?? true);
  const [minVal,  setMinVal]  = useState(String(initial?.minVal ?? ""));
  const [maxVal,  setMaxVal]  = useState(String(initial?.maxVal ?? ""));
  const [stepMode, setStepMode] = useState<StepMode>("n_points");
  const [nPoints, setNPoints] = useState(String(initial?.nPoints ?? 20));
  const [stepSize, setStepSize] = useState("");
  const [limitsError, setLimitsError] = useState<string | null>(null);

  // Auto-populate min/max from EOS limits when fluid or axis changes
  useEffect(() => {
    if (!fluid) return;
    setLimitsError(null);
    invoke<FluidLimits>("get_fluid_limits", { fluid })
      .then((lim) => {
        if (byTemp) {
          setMinVal(lim.t_min.toFixed(2));
          setMaxVal(lim.t_max.toFixed(2));
          setStepSize("5");
        } else {
          setMinVal(lim.p_min.toExponential(3));
          setMaxVal(lim.p_max.toExponential(3));
          setStepSize(((lim.p_max - lim.p_min) / 50).toExponential(2));
        }
      })
      .catch((e) => setLimitsError(String(e)));
  }, [fluid, byTemp]);

  const handleConfirm = () => {
    const min = parseFloat(minVal);
    const max = parseFloat(maxVal);
    if (isNaN(min) || isNaN(max) || max <= min) return;

    let n: number;
    if (stepMode === "n_points") {
      n = parseInt(nPoints, 10);
    } else {
      const step = parseFloat(stepSize);
      if (isNaN(step) || step <= 0) return;
      // n cells of width `step` → n+1 points; bias toward "include endpoints"
      n = Math.max(2, Math.min(500, Math.round((max - min) / step) + 1));
    }
    if (isNaN(n) || n < 2) return;

    onConfirm({ fluid, backend: "HEOS", byTemp, minVal: min, maxVal: max, nPoints: n });
  };

  const label = byTemp ? "T (K)" : "P (Pa)";
  const stepUnit = byTemp ? "K" : "Pa";

  const previewN = (() => {
    if (stepMode === "n_points") return null;
    const min = parseFloat(minVal), max = parseFloat(maxVal), step = parseFloat(stepSize);
    if (!isFinite(min) || !isFinite(max) || !isFinite(step) || step <= 0 || max <= min) return null;
    return Math.max(2, Math.min(500, Math.round((max - min) / step) + 1));
  })();

  return (
    <div className="modal-overlay" onClick={(e) => { if (e.target === e.currentTarget) onCancel?.(); }}>
      <div className="modal-card">
        <div className="modal-title">Saturation Table Setup</div>

        <div className="modal-body">
          <div className="field-group">
            <label>Fluid</label>
            <select value={fluid} onChange={(e) => setFluid(e.target.value)}>
              {(fluids.length > 0 ? fluids : ["Water"]).map((f) => (
                <option key={f} value={f}>{f}</option>
              ))}
            </select>
          </div>

          <div className="field-group">
            <label>Basis</label>
            <div className="readonly-value">{basis.charAt(0).toUpperCase() + basis.slice(1)}</div>
          </div>

          <div className="field-group">
            <label>Sweep by</label>
            <select value={byTemp ? "T" : "P"} onChange={(e) => setByTemp(e.target.value === "T")}>
              <option value="T">Temperature (K)</option>
              <option value="P">Pressure (Pa)</option>
            </select>
          </div>

          {limitsError && <div className="error-msg">{limitsError}</div>}

          <div className="field-row">
            <div className="field-group" style={{ flex: 1 }}>
              <label>Min {label}</label>
              <input type="number" value={minVal} onChange={(e) => setMinVal(e.target.value)} />
            </div>
            <div className="field-group" style={{ flex: 1 }}>
              <label>Max {label}</label>
              <input type="number" value={maxVal} onChange={(e) => setMaxVal(e.target.value)} />
            </div>
          </div>

          <div className="field-row">
            <div className="field-group" style={{ flex: 1 }}>
              <label>Step mode</label>
              <select value={stepMode} onChange={(e) => setStepMode(e.target.value as StepMode)}>
                <option value="n_points">Number of points</option>
                <option value="step">Fixed step</option>
              </select>
            </div>
            {stepMode === "n_points" ? (
              <div className="field-group" style={{ width: 100 }}>
                <label>Points</label>
                <input type="number" value={nPoints} min={2} max={500}
                  onChange={(e) => setNPoints(e.target.value)} />
              </div>
            ) : (
              <div className="field-group" style={{ width: 140 }}>
                <label>Step ({stepUnit})</label>
                <input type="number" value={stepSize}
                  onChange={(e) => setStepSize(e.target.value)} />
              </div>
            )}
          </div>

          {stepMode === "step" && previewN !== null && (
            <div style={{ fontSize: 11, color: "#666", marginTop: -4 }}>
              → {previewN} points (capped at 500)
            </div>
          )}
        </div>

        <div className="modal-footer">
          {onCancel && (
            <button className="btn-secondary" onClick={onCancel}>Cancel</button>
          )}
          <button className="primary" onClick={handleConfirm}>Generate</button>
        </div>
      </div>
    </div>
  );
}
