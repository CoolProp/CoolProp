import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { Basis } from "../App";
import SatSetupDialog, { type SatConfig } from "./SatSetupDialog";
import SaturationTable from "./SaturationTable";

interface Panel {
  id: string;
  stateId: number;
  config: SatConfig;
}

interface Props {
  fluids: string[];
  basis: Basis;
}

const newId = () =>
  (typeof crypto !== "undefined" && "randomUUID" in crypto)
    ? crypto.randomUUID()
    : `${Date.now()}-${Math.random().toString(16).slice(2)}`;

function panelLabel(c: SatConfig): string {
  return `${c.fluid} (${c.byTemp ? "T" : "P"})`;
}

export default function SaturationView({ fluids, basis }: Props) {
  const [panels, setPanels] = useState<Panel[]>([]);
  const [activeId, setActiveId] = useState<string | null>(null);
  const [dialogOpen, setDialogOpen] = useState(false);

  // Free all AbstractStates when the view unmounts (e.g. app close).
  const panelsRef = useRef<Panel[]>([]);
  panelsRef.current = panels;
  useEffect(() => {
    return () => {
      for (const p of panelsRef.current) {
        invoke("free_state", { id: p.stateId }).catch(() => {});
      }
    };
  }, []);

  const addPanel = async (config: SatConfig) => {
    try {
      const stateId = await invoke<number>("create_state", {
        backend: config.backend,
        fluid: config.fluid,
      });
      const id = newId();
      setPanels((prev) => [...prev, { id, stateId, config }]);
      setActiveId(id);
      setDialogOpen(false);
    } catch (e) {
      console.error("create_state failed:", e);
    }
  };

  const closePanel = async (panelId: string) => {
    const p = panels.find((pp) => pp.id === panelId);
    if (!p) return;
    invoke("free_state", { id: p.stateId }).catch(console.error);
    setPanels((prev) => {
      const next = prev.filter((pp) => pp.id !== panelId);
      if (activeId === panelId) {
        setActiveId(next.length > 0 ? next[next.length - 1].id : null);
      }
      return next;
    });
  };

  return (
    <div className="sat-view">
      <div className="sat-tabs">
        {panels.map((p) => (
          <div
            key={p.id}
            className={"sat-tab" + (p.id === activeId ? " active" : "")}
            onClick={() => setActiveId(p.id)}
          >
            <span className="sat-tab-label">{panelLabel(p.config)}</span>
            <button
              className="sat-tab-close"
              onClick={(e) => { e.stopPropagation(); closePanel(p.id); }}
              aria-label="Close panel"
              title="Close"
            >
              ×
            </button>
          </div>
        ))}
        <button className="sat-tab-new" onClick={() => setDialogOpen(true)}>
          + New
        </button>
      </div>

      <div className="sat-panel-stack">
        {panels.map((p) => (
          <div
            key={p.id}
            className="sat-panel"
            hidden={p.id !== activeId}
          >
            <SaturationTable
              stateId={p.stateId}
              config={p.config}
              basis={basis}
            />
          </div>
        ))}
        {panels.length === 0 && (
          <div className="placeholder">
            No saturation tables open. Click <strong>+ New</strong> to create one.
          </div>
        )}
      </div>

      {dialogOpen && (
        <SatSetupDialog
          fluids={fluids}
          basis={basis}
          onConfirm={addPanel}
          onCancel={() => setDialogOpen(false)}
        />
      )}
    </div>
  );
}
