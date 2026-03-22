import { useState, useCallback, useEffect } from "react";
import { useAbstractState } from "../hooks/useAbstractState";
import type { Basis } from "../App";

const BACKENDS = ["HEOS", "INCOMP", "REFPROP"];

interface InputPairDef {
  label: string;
  pair: string;
  v1Label: string;
  v2Label: string;
  v1Default: string;
  v2Default: string;
}

function getInputPairs(basis: Basis): InputPairDef[] {
  const mass = basis === "mass";
  return [
    {
      label: "P, T",
      pair: "PT_INPUTS",
      v1Label: "P (Pa)", v2Label: "T (K)",
      v1Default: "101325", v2Default: "300",
    },
    {
      label: "P, h",
      pair: mass ? "HmassP_INPUTS" : "HmolarP_INPUTS",
      v1Label: mass ? "h (J/kg)" : "h (J/mol)", v2Label: "P (Pa)",
      v1Default: mass ? "100000" : "8000", v2Default: "101325",
    },
    {
      label: "P, s",
      pair: mass ? "PSmass_INPUTS" : "PSmolar_INPUTS",
      v1Label: "P (Pa)", v2Label: mass ? "s (J/kg·K)" : "s (J/mol·K)",
      v1Default: "101325", v2Default: mass ? "1000" : "100",
    },
    {
      label: "ρ, T",
      pair: mass ? "DmassT_INPUTS" : "DmolarT_INPUTS",
      v1Label: mass ? "ρ (kg/m³)" : "ρ (mol/m³)", v2Label: "T (K)",
      v1Default: mass ? "1000" : "55000", v2Default: "300",
    },
    {
      label: "P, Q",
      pair: "PQ_INPUTS",
      v1Label: "P (Pa)", v2Label: "Q (–)",
      v1Default: "101325", v2Default: "0.5",
    },
    {
      label: "Q, T",
      pair: "QT_INPUTS",
      v1Label: "Q (–)", v2Label: "T (K)",
      v1Default: "0.5", v2Default: "373.15",
    },
  ];
}

interface OutputDef {
  label: string;
  param: string;
  unit: string;
}

function getOutputParams(basis: Basis): OutputDef[] {
  const mass = basis === "mass";
  return [
    { label: "T",              param: "T",                          unit: "K"        },
    { label: "P",              param: "P",                          unit: "Pa"       },
    { label: "ρ",              param: mass ? "Dmass"   : "Dmolar",  unit: mass ? "kg/m³"   : "mol/m³"  },
    { label: "h",              param: mass ? "Hmass"   : "Hmolar",  unit: mass ? "J/kg"    : "J/mol"   },
    { label: "s",              param: mass ? "Smass"   : "Smolar",  unit: mass ? "J/kg·K"  : "J/mol·K" },
    { label: "u",              param: mass ? "Umass"   : "Umolar",  unit: mass ? "J/kg"    : "J/mol"   },
    { label: "Cv",             param: mass ? "Cvmass"  : "Cvmolar", unit: mass ? "J/kg·K"  : "J/mol·K" },
    { label: "Cp",             param: mass ? "Cpmass"  : "Cpmolar", unit: mass ? "J/kg·K"  : "J/mol·K" },
    { label: "Speed of sound", param: "speed_of_sound",             unit: "m/s"      },
    { label: "Viscosity",      param: "viscosity",                  unit: "Pa·s"     },
    { label: "Conductivity",   param: "conductivity",               unit: "W/m·K"    },
    { label: "Q",              param: "Q",                          unit: "–"        },
    { label: "Phase",          param: "Phase",                      unit: "–"        },
  ];
}

const PHASE_NAMES: Record<number, string> = {
  0: "liquid", 1: "supercritical", 2: "supercritical_gas",
  3: "supercritical_liquid", 4: "critical_point", 5: "gas",
  6: "two_phase", 7: "unknown", 8: "not_imposed",
};

interface Props {
  fluids: string[];
  basis: Basis;
}

export default function PropertyCalculator({ fluids, basis }: Props) {
  const [backend, setBackend] = useState("HEOS");
  const [fluid, setFluid] = useState("Water");
  const [pairIdx, setPairIdx] = useState(0);
  const [v1, setV1] = useState(getInputPairs(basis)[0].v1Default);
  const [v2, setV2] = useState(getInputPairs(basis)[0].v2Default);
  const [results, setResults] = useState<Record<string, number>>({});
  const [calcError, setCalcError] = useState<string | null>(null);
  const [computing, setComputing] = useState(false);

  const state = useAbstractState(backend, fluid);

  // Reset input values and clear results when basis changes
  useEffect(() => {
    const pairs = getInputPairs(basis);
    setV1(pairs[pairIdx].v1Default);
    setV2(pairs[pairIdx].v2Default);
    setResults({});
    setCalcError(null);
  }, [basis]); // intentionally omit pairIdx — only trigger on basis flip

  const handlePairChange = (idx: number) => {
    const pairs = getInputPairs(basis);
    setPairIdx(idx);
    setV1(pairs[idx].v1Default);
    setV2(pairs[idx].v2Default);
  };

  const compute = useCallback(async () => {
    if (!state.ready) return;
    setComputing(true);
    setCalcError(null);
    try {
      const ip = getInputPairs(basis)[pairIdx];
      await state.update(ip.pair, parseFloat(v1), parseFloat(v2));
      const out: Record<string, number> = {};
      for (const { param } of getOutputParams(basis)) {
        try { out[param] = await state.get(param); } catch { /* param unavailable */ }
      }
      setResults(out);
    } catch (e) {
      setCalcError(String(e));
    } finally {
      setComputing(false);
    }
  }, [state, basis, pairIdx, v1, v2]);

  const inputPairs = getInputPairs(basis);
  const outputParams = getOutputParams(basis);
  const ip = inputPairs[pairIdx];

  return (
    <div className="calc-layout">
      <div className="calc-controls">
        <div className="field-group">
          <label>Backend</label>
          <select value={backend} onChange={(e) => setBackend(e.target.value)}>
            {BACKENDS.map((b) => <option key={b} value={b}>{b}</option>)}
          </select>
        </div>

        <div className="field-group">
          <label>Fluid</label>
          <select value={fluid} onChange={(e) => setFluid(e.target.value)}>
            {(fluids.length > 0 ? fluids : ["Water"]).map((f) => (
              <option key={f} value={f}>{f}</option>
            ))}
          </select>
        </div>

        {state.error && <div className="error-msg">{state.error}</div>}

        <div className="field-group">
          <label>Input pair</label>
          <select
            value={pairIdx}
            onChange={(e) => handlePairChange(Number(e.target.value))}
          >
            {inputPairs.map((p, i) => (
              <option key={p.pair} value={i}>{p.label}</option>
            ))}
          </select>
        </div>

        <div className="field-group">
          <label>{ip.v1Label}</label>
          <input type="number" value={v1} onChange={(e) => setV1(e.target.value)} />
        </div>

        <div className="field-group">
          <label>{ip.v2Label}</label>
          <input type="number" value={v2} onChange={(e) => setV2(e.target.value)} />
        </div>

        <button
          className="primary"
          onClick={compute}
          disabled={!state.ready || computing}
        >
          {computing ? "Computing…" : "Calculate"}
        </button>

        {calcError && <div className="error-msg">{calcError}</div>}
      </div>

      <div className="calc-results">
        <table className="results-table">
          <thead>
            <tr>
              <th>Property</th>
              <th style={{ textAlign: "right" }}>Value</th>
              <th>Unit</th>
            </tr>
          </thead>
          <tbody>
            {outputParams.map(({ label, param, unit }) => {
              const val = results[param];
              let display = "—";
              if (val !== undefined) {
                display = param === "Phase"
                  ? (PHASE_NAMES[Math.round(val)] ?? String(Math.round(val)))
                  : val.toPrecision(7);
              }
              return (
                <tr key={param}>
                  <td>{label}</td>
                  <td className="num">{display}</td>
                  <td className="unit">{unit}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
