import { useState, useEffect, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import PropertyCalculator from "./components/PropertyCalculator";
import SaturationView from "./components/SaturationView";
import HumidAirCalculator from "./components/HumidAirCalculator";
import AboutModal from "./components/AboutModal";
import AddFluidDialog, { registerUserFluids } from "./components/AddFluidDialog";
import UpdateChecker from "./components/UpdateChecker";
import SponsorSplash from "./components/SponsorSplash";
import { SPONSOR_URL } from "./constants";

type Tab = "calculator" | "saturation" | "humidair" | "diagram";
export type Basis = "mass" | "molar";

export default function App() {
  const [tab, setTab] = useState<Tab>("calculator");
  const [fluids, setFluids] = useState<string[]>([]);
  const [incompFluids, setIncompFluids] = useState<string[]>([]);
  const [basis, setBasis] = useState<Basis>("mass");
  const [aboutOpen, setAboutOpen] = useState(false);
  const [addFluidOpen, setAddFluidOpen] = useState(false);

  const refreshIncompFluids = useCallback(() => {
    invoke<string[]>("get_incompressible_fluids_list").then(setIncompFluids).catch(console.error);
  }, []);

  useEffect(() => {
    invoke<string[]>("get_fluids_list").then(setFluids).catch(console.error);
    // Re-register any user-defined incompressible fluids persisted from
    // earlier sessions, then load the (now complete) INCOMP fluid list.
    registerUserFluids().finally(refreshIncompFluids);
  }, [refreshIncompFluids]);

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
        <div className="header-right">
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
          <a
            className="tab-btn sponsor-btn"
            href={SPONSOR_URL}
            target="_blank"
            rel="noreferrer"
            title="Support CoolProp on GitHub Sponsors"
          >
            💚 Sponsor
          </a>
          <button
            className="tab-btn about-btn"
            onClick={() => setAboutOpen(true)}
            title="About / third-party notices"
          >
            About
          </button>
        </div>
      </header>
      <main className="app-main">
        {/* All tabs stay mounted so state (results, isolines, sat-table panels) persists. */}
        <div className="tab-pane" hidden={tab !== "calculator"}>
          <PropertyCalculator
            fluids={fluids}
            incompFluids={incompFluids}
            basis={basis}
            onAddFluid={() => setAddFluidOpen(true)}
          />
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
      {aboutOpen && <AboutModal onClose={() => setAboutOpen(false)} />}
      {addFluidOpen && (
        <AddFluidDialog
          onSaved={() => {
            setAddFluidOpen(false);
            refreshIncompFluids();
          }}
          onCancel={() => setAddFluidOpen(false)}
        />
      )}
      <UpdateChecker />
      <SponsorSplash />
    </div>
  );
}
