import { useState, useEffect, useRef } from "react";
import type { ReactNode } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { Basis } from "../App";
import type { SatConfig } from "./SatSetupDialog";

interface SatPoint {
  independent: number;
  liquid: Record<string, number>;
  vapor: Record<string, number>;
}

interface ColDef {
  key: string;
  headerNode: ReactNode;
  headerText: string; // plain-text fallback for clipboard
  get: (row: SatPoint) => number | undefined;
  fmt: (v: number) => string;
}

const f6 = (v: number) => v.toPrecision(6);
const f5 = (v: number) => v.toPrecision(5);
const f4 = (v: number) => v.toPrecision(4);

/** Render `sym<sub>sub</sub> (units)` and return both the JSX and a plain-text twin. */
function H(sym: string, sub: string | null, units: string): { node: ReactNode; text: string } {
  const node = (
    <>
      {sym}
      {sub !== null && <sub>{sub}</sub>}
      {" "}({units})
    </>
  );
  const text = sub !== null ? `${sym}_${sub} (${units})` : `${sym} (${units})`;
  return { node, text };
}

function getColumns(basis: Basis): ColDef[] {
  const mass = basis === "mass";
  const Trow  = H("T", null, "K");
  const Prow  = H("P", null, "kPa");
  const rhoL  = H("ρ", "L", mass ? "kg/m³" : "mol/m³");
  const rhoV  = H("ρ", "V", mass ? "kg/m³" : "mol/m³");
  const hL    = H("h", "L", mass ? "kJ/kg" : "J/mol");
  const hV    = H("h", "V", mass ? "kJ/kg" : "J/mol");
  const sL    = H("s", "L", mass ? "kJ/kg·K" : "J/mol·K");
  const sV    = H("s", "V", mass ? "kJ/kg·K" : "J/mol·K");
  const cpL   = H("cp", "L", mass ? "kJ/kg·K" : "J/mol·K");
  const cpV   = H("cp", "V", mass ? "kJ/kg·K" : "J/mol·K");
  const muL   = H("μ", "L", "μPa·s");
  const muV   = H("μ", "V", "μPa·s");
  const lamL  = H("λ", "L", "mW/m·K");
  const lamV  = H("λ", "V", "mW/m·K");

  return [
    { key: "T",   headerNode: Trow.node, headerText: Trow.text,
      get: (r) => r.liquid.T, fmt: f6 },
    { key: "P",   headerNode: Prow.node, headerText: Prow.text,
      get: (r) => r.liquid.P !== undefined ? r.liquid.P / 1000 : undefined, fmt: f6 },
    { key: "rL",  headerNode: rhoL.node, headerText: rhoL.text,
      get: (r) => mass ? r.liquid.Dmass : r.liquid.Dmolar, fmt: f6 },
    { key: "rV",  headerNode: rhoV.node, headerText: rhoV.text,
      get: (r) => mass ? r.vapor.Dmass  : r.vapor.Dmolar, fmt: f6 },
    { key: "hL",  headerNode: hL.node, headerText: hL.text,
      get: (r) => { const v = mass ? r.liquid.Hmass : r.liquid.Hmolar; return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f6 },
    { key: "hV",  headerNode: hV.node, headerText: hV.text,
      get: (r) => { const v = mass ? r.vapor.Hmass  : r.vapor.Hmolar;  return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f6 },
    { key: "sL",  headerNode: sL.node, headerText: sL.text,
      get: (r) => { const v = mass ? r.liquid.Smass : r.liquid.Smolar; return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f5 },
    { key: "sV",  headerNode: sV.node, headerText: sV.text,
      get: (r) => { const v = mass ? r.vapor.Smass  : r.vapor.Smolar;  return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f5 },
    { key: "cpL", headerNode: cpL.node, headerText: cpL.text,
      get: (r) => { const v = mass ? r.liquid.Cpmass : r.liquid.Cpmolar; return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f5 },
    { key: "cpV", headerNode: cpV.node, headerText: cpV.text,
      get: (r) => { const v = mass ? r.vapor.Cpmass  : r.vapor.Cpmolar;  return v !== undefined ? (mass ? v/1000 : v) : undefined; }, fmt: f5 },
    { key: "muL", headerNode: muL.node, headerText: muL.text,
      get: (r) => r.liquid.viscosity !== undefined ? r.liquid.viscosity * 1e6 : undefined, fmt: f4 },
    { key: "muV", headerNode: muV.node, headerText: muV.text,
      get: (r) => r.vapor.viscosity  !== undefined ? r.vapor.viscosity  * 1e6 : undefined, fmt: f4 },
    { key: "lL",  headerNode: lamL.node, headerText: lamL.text,
      get: (r) => r.liquid.conductivity !== undefined ? r.liquid.conductivity * 1000 : undefined, fmt: f4 },
    { key: "lV",  headerNode: lamV.node, headerText: lamV.text,
      get: (r) => r.vapor.conductivity  !== undefined ? r.vapor.conductivity  * 1000 : undefined, fmt: f4 },
  ];
}

function getParams(basis: Basis): string[] {
  const mass = basis === "mass";
  return [
    mass ? "Dmass" : "Dmolar",
    mass ? "Hmass" : "Hmolar",
    mass ? "Smass" : "Smolar",
    mass ? "Cpmass" : "Cpmolar",
    "viscosity",
    "conductivity",
  ];
}

interface Props {
  stateId: number;
  config: SatConfig;
  basis: Basis;
}

export default function SaturationTable({ stateId, config, basis }: Props) {
  const [rows, setRows] = useState<SatPoint[]>([]);
  const [computing, setComputing] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [copied, setCopied] = useState(false);
  const copyTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Compute whenever stateId, config, or basis changes
  useEffect(() => {
    setComputing(true);
    setError(null);
    setRows([]);
    invoke<SatPoint[]>("compute_saturation_table", {
      id: stateId,
      byTemperature: config.byTemp,
      minVal: config.minVal,
      maxVal: config.maxVal,
      nPoints: config.nPoints,
      params: getParams(basis),
    })
      .then(setRows)
      .catch((e) => setError(String(e)))
      .finally(() => setComputing(false));
  }, [stateId, config, basis]);

  const copyToClipboard = () => {
    const cols = getColumns(basis);
    const header = cols.map((c) => c.headerText).join("\t");
    const body = rows
      .map((row) => cols.map((c) => { const v = c.get(row); return v !== undefined && isFinite(v) ? c.fmt(v) : ""; }).join("\t"))
      .join("\n");
    navigator.clipboard.writeText(header + "\n" + body).then(() => {
      setCopied(true);
      if (copyTimerRef.current) clearTimeout(copyTimerRef.current);
      copyTimerRef.current = setTimeout(() => setCopied(false), 2000);
    });
  };

  const columns = getColumns(basis);

  return (
    <div className="sat-layout">
      <div className="sat-controls">
        <div className="sat-info">
          <span className="sat-info-fluid">{config.fluid}</span>
          <span className="sat-info-detail">
            {config.backend} · {config.byTemp ? "T" : "P"} sweep ·{" "}
            {config.byTemp
              ? `${config.minVal.toFixed(1)}–${config.maxVal.toFixed(1)} K`
              : `${config.minVal.toExponential(2)}–${config.maxVal.toExponential(2)} Pa`}
            {" · "}{config.nPoints} pts
          </span>
        </div>

        {rows.length > 0 && (
          <button className="primary" onClick={copyToClipboard}>
            {copied ? "Copied!" : "Copy TSV"}
          </button>
        )}

        {computing && <span className="sat-status">Computing…</span>}
        {error && <div className="error-msg">{error}</div>}
      </div>

      {rows.length > 0 && (
        <div className="sat-table-wrap">
          <table className="sat-table">
            <thead>
              <tr>{columns.map((c) => <th key={c.key}>{c.headerNode}</th>)}</tr>
            </thead>
            <tbody>
              {rows.map((row, i) => (
                <tr key={i}>
                  {columns.map((c) => {
                    const v = c.get(row);
                    return <td key={c.key}>{v !== undefined && isFinite(v) ? c.fmt(v) : "—"}</td>;
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
