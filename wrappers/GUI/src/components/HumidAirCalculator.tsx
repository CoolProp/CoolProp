import { useState, useCallback, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import Plot from "react-plotly.js";
import { useDraggableSplit } from "../hooks/useDraggableSplit";

// ── Input parameter definitions ───────────────────────────────────────────────

interface HAInputDef {
  name: string;   // HAPropsSI key
  label: string;  // display label
  unit: string;
  defaultVal: string;
}

const HA_INPUTS: HAInputDef[] = [
  { name: "T", label: "Dry-bulb temp",  unit: "K",      defaultVal: "293.15" },
  { name: "R", label: "Rel. humidity",  unit: "0–1",    defaultVal: "0.5"    },
  { name: "W", label: "Humidity ratio", unit: "kg/kg_da", defaultVal: "0.0073" },
  { name: "D", label: "Dew-point",      unit: "K",      defaultVal: "282"    },
  { name: "B", label: "Wet-bulb",       unit: "K",      defaultVal: "287"    },
];

// ── Output definitions ────────────────────────────────────────────────────────

interface HAOutputDef {
  key: string;
  label: string;
  unit: string;
}

const HA_OUTPUTS: HAOutputDef[] = [
  { key: "T",   label: "Dry-bulb temp",   unit: "K"          },
  { key: "B",   label: "Wet-bulb",        unit: "K"          },
  { key: "D",   label: "Dew-point",       unit: "K"          },
  { key: "R",   label: "Rel. humidity",   unit: "–"          },
  { key: "W",   label: "Humidity ratio",  unit: "kg/kg_da"   },
  { key: "H",   label: "Enthalpy",        unit: "J/kg_da"    },
  { key: "S",   label: "Entropy",         unit: "J/kg_da·K"  },
  { key: "V",   label: "Specific volume", unit: "m³/kg_da"   },
  { key: "C",   label: "Cp",              unit: "J/kg_da·K"  },
  { key: "M",   label: "Viscosity",       unit: "Pa·s"       },
  { key: "K",   label: "Conductivity",    unit: "W/m·K"      },
];

const OUTPUT_KEYS = HA_OUTPUTS.map((o) => o.key);

// ── Psychrometric chart helpers ───────────────────────────────────────────────

const T_RANGE = { min: 233.15, max: 333.15, n: 60 };

/** Compute W (humidity ratio) from T and R at pressure P (default 101325 Pa). */
async function computeW(T: number, R: number, P = 101325): Promise<number | null> {
  try {
    const res = await invoke<Record<string, number>>("compute_humid_air", {
      n1: "P", v1: P,
      n2: "T", v2: T,
      n3: "R", v3: R,
      outputs: ["W"],
    });
    return res["W"] ?? null;
  } catch {
    return null;
  }
}

// ── Component ─────────────────────────────────────────────────────────────────

export default function HumidAirCalculator() {
  const [pressure, setPressure] = useState("101325");
  const [input2Idx, setInput2Idx] = useState(0); // index into HA_INPUTS (T by default)
  const [input3Idx, setInput3Idx] = useState(1); // index into HA_INPUTS (R by default)
  const [v2, setV2] = useState(HA_INPUTS[0].defaultVal);
  const [v3, setV3] = useState(HA_INPUTS[1].defaultVal);

  const [results, setResults] = useState<Record<string, number> | null>(null);
  const [calcError, setCalcError] = useState<string | null>(null);
  const [computing, setComputing] = useState(false);

  // Chart traces: iso-RH lines + saturation curve, computed once on mount
  const [isoTraces, setIsoTraces] = useState<Plotly.Data[]>([]);
  const [chartError, setChartError] = useState<string | null>(null);
  const chartReady = useRef(false);

  // Iso-property lines through the most recently computed state point
  const [pointIsolines, setPointIsolines] = useState<Plotly.Data[]>([]);
  const [showIsolines, setShowIsolines] = useState(true);

  const { width: leftWidth, startDrag } = useDraggableSplit(260, 200, 600);

  // Live hover readout: properties at the mouse position over the chart.
  const [hoverProps, setHoverProps] = useState<
    { T: number; W: number; B: number; D: number; R: number; H: number; V: number } | null
  >(null);
  const graphDivRef = useRef<HTMLDivElement | null>(null);

  // Build chart background once
  useEffect(() => {
    let cancelled = false;
    (async () => {
      const RH_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
      const Ts: number[] = [];
      for (let i = 0; i < T_RANGE.n; i++) {
        Ts.push(T_RANGE.min + (i / (T_RANGE.n - 1)) * (T_RANGE.max - T_RANGE.min));
      }

      const traces: Plotly.Data[] = [];

      for (const R of RH_LEVELS) {
        const xs: number[] = [];
        const ys: number[] = [];
        for (const T of Ts) {
          const W = await computeW(T, R);
          if (W !== null && W >= 0) {
            xs.push(T);
            ys.push(W);
          }
        }
        if (cancelled) return;
        if (xs.length === 0) continue;

        const isSat = R === 1.0;
        const color = isSat ? "#1c2b3a" : "#bbb";
        const width = isSat ? 1.5 : 1;

        // Line trace
        traces.push({
          x: xs,
          y: ys,
          type: "scatter",
          mode: "lines",
          line: { color, width },
          showlegend: false,
          hoverinfo: "skip",
        } as Plotly.Data);

        // Label at rightmost point
        const labelText = isSat ? "Sat." : `${Math.round(R * 100)}%`;
        traces.push({
          x: [xs[xs.length - 1]],
          y: [ys[ys.length - 1]],
          type: "scatter",
          mode: "text",
          text: [labelText],
          textposition: "middle right",
          textfont: { size: 10, color },
          showlegend: false,
          hoverinfo: "skip",
        } as Plotly.Data);
      }

      if (!cancelled) {
        if (traces.length === 0) {
          setChartError("Failed to compute psychrometric chart background.");
        } else {
          setIsoTraces(traces);
          chartReady.current = true;
        }
      }
    })().catch((e: unknown) => {
      if (!cancelled) setChartError(String(e));
    });
    return () => { cancelled = true; };
  }, []);

  const handleInput2Change = (idx: number) => {
    setInput2Idx(idx);
    setV2(HA_INPUTS[idx].defaultVal);
  };

  const handleInput3Change = (idx: number) => {
    setInput3Idx(idx);
    setV3(HA_INPUTS[idx].defaultVal);
  };

  const compute = useCallback(async () => {
    const in2 = HA_INPUTS[input2Idx];
    const in3 = HA_INPUTS[input3Idx];
    if (in2.name === in3.name) {
      setCalcError(`Input 2 and Input 3 are both "${in2.label}" — they must be different.`);
      return;
    }
    setComputing(true);
    setCalcError(null);
    try {
      const res = await invoke<Record<string, number>>("compute_humid_air", {
        n1: "P",     v1: parseFloat(pressure),
        n2: in2.name, v2: parseFloat(v2),
        n3: in3.name, v3: parseFloat(v3),
        outputs: OUTPUT_KEYS,
      });
      setResults(res);
    } catch (e) {
      setCalcError(String(e));
    } finally {
      setComputing(false);
    }
  }, [pressure, input2Idx, input3Idx, v2, v3]);

  // Recompute isolines whenever the state point, pressure, or toggle changes.
  // Lines are clipped at the saturation curve by filtering points where R > 1.
  useEffect(() => {
    if (!showIsolines || !results) {
      setPointIsolines([]);
      return;
    }
    let cancelled = false;
    (async () => {
      const P = parseFloat(pressure);
      const Bp = results["B"], Dp = results["D"], Vp = results["V"], Wp = results["W"];

      const Ts: number[] = [];
      for (let i = 0; i < T_RANGE.n; i++) {
        Ts.push(T_RANGE.min + (i / (T_RANGE.n - 1)) * (T_RANGE.max - T_RANGE.min));
      }

      const isoFor = async (held: "B" | "V", heldVal: number) => {
        const xs: number[] = [];
        const ys: number[] = [];
        for (const T of Ts) {
          if (cancelled) return { xs: [], ys: [] };
          try {
            const r = await invoke<Record<string, number>>("compute_humid_air", {
              n1: "P", v1: P,
              n2: "T", v2: T,
              n3: held, v3: heldVal,
              outputs: ["W", "R"],
            });
            const W = r["W"], R = r["R"];
            // Drop points outside the saturation envelope (R > 1).
            if (W !== undefined && R !== undefined && W >= 0 && R <= 1.0 && isFinite(W)) {
              xs.push(T);
              ys.push(W);
            }
          } catch { /* skip */ }
        }
        return { xs, ys };
      };

      const traces: Plotly.Data[] = [];
      const addLine = (xs: number[], ys: number[], color: string, label: string) => {
        if (xs.length < 2) return;
        traces.push({
          x: xs, y: ys,
          type: "scatter", mode: "lines",
          line: { color, width: 1.5, dash: "dash" },
          showlegend: false, hoverinfo: "skip",
        } as Plotly.Data);
        traces.push({
          x: [xs[0]], y: [ys[0]],
          type: "scatter", mode: "text",
          text: [label], textposition: "middle left",
          textfont: { size: 10, color },
          showlegend: false, hoverinfo: "skip",
        } as Plotly.Data);
      };

      if (Bp !== undefined && isFinite(Bp)) {
        const { xs, ys } = await isoFor("B", Bp);
        if (!cancelled) addLine(xs, ys, "#2980b9", `B=${Bp.toFixed(2)} K`);
      }
      if (Vp !== undefined && isFinite(Vp)) {
        const { xs, ys } = await isoFor("V", Vp);
        if (!cancelled) addLine(xs, ys, "#27ae60", `v=${Vp.toFixed(3)} m³/kg`);
      }
      // Iso-dewpoint == iso-W (at fixed P): horizontal segment from T=Dp rightward.
      if (Dp !== undefined && Wp !== undefined && isFinite(Wp)) {
        addLine(
          [Dp, T_RANGE.max],
          [Wp, Wp],
          "#c0392b",
          `D=${Dp.toFixed(2)} K`,
        );
      }

      if (!cancelled) setPointIsolines(traces);
    })();
    return () => { cancelled = true; };
  }, [results, pressure, showIsolines]);

  // Live hover: converts cursor position to (T, W) using Plotly's axis
  // calibrators, calls HAPropsSI for the surrounding properties, and only
  // displays them when the point is below the saturation curve. Throttled to
  // ~80 ms so we don't flood the FFI on rapid mousemoves.
  useEffect(() => {
    const gd = graphDivRef.current;
    if (!gd) return;
    let lastT = 0;
    let inflight = false;
    let pendingTimer: ReturnType<typeof setTimeout> | null = null;

    const tryCompute = (T: number, W: number) => {
      if (inflight) return;
      inflight = true;
      invoke<Record<string, number>>("compute_humid_air", {
        n1: "P", v1: parseFloat(pressure),
        n2: "T", v2: T,
        n3: "W", v3: W,
        outputs: ["B", "D", "R", "H", "V"],
      })
        .then((res) => {
          // Only show if we're at or below saturation (R ≤ 1).
          if (
            res["R"] !== undefined && res["R"] <= 1.0 && isFinite(res["R"]) &&
            isFinite(res["B"] ?? NaN) && isFinite(res["D"] ?? NaN)
          ) {
            setHoverProps({
              T, W,
              B: res["B"], D: res["D"], R: res["R"], H: res["H"], V: res["V"],
            });
          } else {
            setHoverProps(null);
          }
        })
        .catch(() => setHoverProps(null))
        .finally(() => { inflight = false; });
    };

    const handleMove = (ev: MouseEvent) => {
      const layout = (gd as unknown as { _fullLayout?: any })._fullLayout;
      if (!layout || !layout.xaxis || !layout.yaxis) return;
      const rect = gd.getBoundingClientRect();
      const margin = layout._size;
      const xPx = ev.clientX - rect.left - margin.l;
      const yPx = ev.clientY - rect.top - margin.t;
      // Outside the plotting area
      if (xPx < 0 || yPx < 0 || xPx > margin.w || yPx > margin.h) {
        setHoverProps(null);
        return;
      }
      const T = layout.xaxis.p2c(xPx);
      const W = layout.yaxis.p2c(yPx);
      // The saturation envelope reaches W ≈ 0.3 at the top of the T range,
      // and the chart's y-axis can scroll well above that. Only filter the
      // physically impossible (W < 0); supersaturation is filtered by the
      // R ≤ 1 check on the result.
      if (T < T_RANGE.min || T > T_RANGE.max + 8 || W < 0) {
        setHoverProps(null);
        return;
      }
      const now = Date.now();
      if (now - lastT < 80) {
        if (pendingTimer) clearTimeout(pendingTimer);
        pendingTimer = setTimeout(() => tryCompute(T, W), 80);
        return;
      }
      lastT = now;
      tryCompute(T, W);
    };
    const handleLeave = () => setHoverProps(null);

    gd.addEventListener("mousemove", handleMove);
    gd.addEventListener("mouseleave", handleLeave);
    return () => {
      gd.removeEventListener("mousemove", handleMove);
      gd.removeEventListener("mouseleave", handleLeave);
      if (pendingTimer) clearTimeout(pendingTimer);
    };
  // graphDivRef is set after first paint; we re-run on pressure change so the
  // closure picks up the new P. Re-attaching listeners is cheap.
  }, [pressure, isoTraces.length]);

  // State point trace (shown when results exist)
  const stateTrace: Plotly.Data | null =
    results && results["T"] !== undefined && results["W"] !== undefined
      ? {
          x: [results["T"]],
          y: [results["W"]],
          type: "scatter",
          mode: "markers",
          marker: { color: "#e74c3c", size: 10, symbol: "circle" },
          showlegend: false,
          hovertemplate: "T=%{x:.2f} K<br>W=%{y:.5f} kg/kg<extra></extra>",
        }
      : null;

  const plotData: Plotly.Data[] = stateTrace
    ? [...isoTraces, ...pointIsolines, stateTrace]
    : isoTraces;

  const in2 = HA_INPUTS[input2Idx];
  const in3 = HA_INPUTS[input3Idx];

  return (
    <div
      className="ha-layout"
      style={{ gridTemplateColumns: `${leftWidth}px 6px 1fr` }}
    >
      {/* ── Left panel: inputs + results ── */}
      <div className="ha-left">
        <div className="calc-controls">
          <div className="field-group">
            <label>P (Pa)</label>
            <input
              type="number"
              value={pressure}
              onChange={(e) => setPressure(e.target.value)}
            />
          </div>

          <div className="field-group">
            <label>Input 2</label>
            <select
              value={input2Idx}
              onChange={(e) => handleInput2Change(Number(e.target.value))}
            >
              {HA_INPUTS.map((inp, i) => (
                <option key={inp.name} value={i}>
                  {inp.label} ({inp.name})
                </option>
              ))}
            </select>
          </div>

          <div className="field-group">
            <label>
              {in2.label} ({in2.unit})
            </label>
            <input
              type="number"
              value={v2}
              onChange={(e) => setV2(e.target.value)}
            />
          </div>

          <div className="field-group">
            <label>Input 3</label>
            <select
              value={input3Idx}
              onChange={(e) => handleInput3Change(Number(e.target.value))}
            >
              {HA_INPUTS.map((inp, i) => (
                <option key={inp.name} value={i}>
                  {inp.label} ({inp.name})
                </option>
              ))}
            </select>
          </div>

          <div className="field-group">
            <label>
              {in3.label} ({in3.unit})
            </label>
            <input
              type="number"
              value={v3}
              onChange={(e) => setV3(e.target.value)}
            />
          </div>

          <button
            className="primary"
            onClick={compute}
            disabled={computing}
          >
            {computing ? "Computing…" : "Calculate"}
          </button>

          <label className="checkbox-row">
            <input
              type="checkbox"
              checked={showIsolines}
              onChange={(e) => setShowIsolines(e.target.checked)}
            />
            Show isolines
          </label>

          {calcError && <div className="error-msg">{calcError}</div>}
        </div>

        {results && (
          <div className="calc-results" style={{ marginTop: 12 }}>
            <table className="results-table">
              <thead>
                <tr>
                  <th>Property</th>
                  <th style={{ textAlign: "right" }}>Value</th>
                  <th>Unit</th>
                </tr>
              </thead>
              <tbody>
                {HA_OUTPUTS.map(({ key, label, unit }) => {
                  const val = results[key];
                  const display =
                    val !== undefined ? val.toPrecision(7) : "—";
                  return (
                    <tr key={key}>
                      <td>{label}</td>
                      <td className="num">{display}</td>
                      <td className="unit">{unit}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        )}
      </div>

      <div
        className="calc-divider"
        onMouseDown={startDrag}
        role="separator"
        aria-orientation="vertical"
      />

      {/* ── Right panel: psychrometric chart ── */}
      <div className="ha-chart">
        {chartError && (
          <div className="error-msg" style={{ padding: 16 }}>{chartError}</div>
        )}
        <Plot
          data={plotData}
          layout={{
            margin: { l: 60, r: 70, t: 30, b: 50 },
            xaxis: {
              title: { text: "Dry-bulb temperature (K)" },
              range: [T_RANGE.min, T_RANGE.max + 8],
            },
            yaxis: {
              title: { text: "Humidity ratio W (kg/kg_da)" },
              rangemode: "tozero",
            },
            plot_bgcolor: "#fafafa",
            paper_bgcolor: "#fff",
            font: { size: 12 },
          }}
          config={{ responsive: true, displayModeBar: false }}
          style={{ width: "100%", height: "100%" }}
          useResizeHandler
          onInitialized={(_, gd) => { graphDivRef.current = gd as HTMLDivElement; }}
          onUpdate={(_, gd) => { graphDivRef.current = gd as HTMLDivElement; }}
        />
        {hoverProps && (
          <div className="ha-hover-readout">
            <div><span>T</span><span>{hoverProps.T.toFixed(2)} K</span></div>
            <div><span>W</span><span>{hoverProps.W.toFixed(5)} kg/kg</span></div>
            <div><span>RH</span><span>{(hoverProps.R * 100).toFixed(1)} %</span></div>
            <div><span>B</span><span>{hoverProps.B.toFixed(2)} K</span></div>
            <div><span>D</span><span>{hoverProps.D.toFixed(2)} K</span></div>
            <div><span>h</span><span>{(hoverProps.H / 1000).toFixed(2)} kJ/kg</span></div>
            <div><span>v</span><span>{hoverProps.V.toFixed(4)} m³/kg</span></div>
          </div>
        )}
      </div>
    </div>
  );
}
