import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import PropertyCalculator from "./components/PropertyCalculator";
import SaturationView from "./components/SaturationView";
import HumidAirCalculator from "./components/HumidAirCalculator";

type Tab = "calculator" | "saturation" | "humidair" | "diagram";
export type Basis = "mass" | "molar";

export default function App() {
  const [tab, setTab] = useState<Tab>("calculator");
  const [fluids, setFluids] = useState<string[]>([]);
  const [basis, setBasis] = useState<Basis>("mass");

  useEffect(() => {
    invoke<string[]>("get_fluids_list").then(setFluids).catch(console.error);
  }, []);

  return (
    <div className="app">
      <header className="app-header">
        <span className="app-title">CoolProp</span>
        <nav className="tab-bar">
          {(["calculator", "saturation", "humidair", "diagram"] as Tab[]).map((t) => (
            <button
              key={t}
              className={"tab-btn" + (tab === t ? " active" : "")}
              onClick={() => setTab(t)}
            >
              {t === "humidair" ? "Humid Air" : t.charAt(0).toUpperCase() + t.slice(1)}
            </button>
          ))}
        </nav>
        {tab !== "humidair" && (
          <div className="seg-ctrl">
            {(["mass", "molar"] as Basis[]).map((b) => (
              <button
                key={b}
                className={"seg-btn" + (basis === b ? " active" : "")}
                onClick={() => setBasis(b)}
              >
                {b.charAt(0).toUpperCase() + b.slice(1)}
              </button>
            ))}
          </div>
        )}
      </header>
      <main className="app-main">
        {/* All tabs stay mounted so state (results, isolines, sat-table panels) persists. */}
        <div className="tab-pane" hidden={tab !== "calculator"}>
          <PropertyCalculator fluids={fluids} basis={basis} />
        </div>
        <div className="tab-pane" hidden={tab !== "saturation"}>
          <SaturationView fluids={fluids} basis={basis} />
        </div>
        <div className="tab-pane" hidden={tab !== "humidair"}>
          <HumidAirCalculator />
        </div>
        <div className="tab-pane" hidden={tab !== "diagram"}>
          <div className="placeholder">Diagram — coming soon</div>
        </div>
      </main>
    </div>
  );
}
