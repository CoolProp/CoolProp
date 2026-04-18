import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { Basis } from "../App";
import SatSetupDialog, { type SatConfig } from "./SatSetupDialog";
import SaturationTable from "./SaturationTable";

interface ActiveState {
  id: number;
  config: SatConfig;
}

interface Props {
  fluids: string[];
  basis: Basis;
}

export default function SaturationView({ fluids, basis }: Props) {
  const [active, setActive] = useState<ActiveState | null>(null);
  const [dialogOpen, setDialogOpen] = useState(true);
  const stateIdRef = useRef<number | null>(null);

  // Free AbstractState on unmount
  useEffect(() => {
    return () => {
      if (stateIdRef.current !== null) {
        invoke("free_state", { id: stateIdRef.current }).catch(console.error);
      }
    };
  }, []);

  const handleConfirm = async (config: SatConfig) => {
    // Free the previous state if any
    if (stateIdRef.current !== null) {
      await invoke("free_state", { id: stateIdRef.current }).catch(console.error);
      stateIdRef.current = null;
    }
    try {
      const id = await invoke<number>("create_state", {
        backend: config.backend,
        fluid: config.fluid,
      });
      stateIdRef.current = id;
      setActive({ id, config });
      setDialogOpen(false);
    } catch (e) {
      console.error("create_state failed:", e);
    }
  };

  return (
    <>
      {/* Table is always rendered once created; dialog overlays it */}
      {active ? (
        <SaturationTable
          stateId={active.id}
          config={active.config}
          basis={basis}
          onNewTable={() => setDialogOpen(true)}
        />
      ) : (
        <div className="placeholder">Configure a saturation table to begin.</div>
      )}

      {dialogOpen && (
        <SatSetupDialog
          fluids={fluids}
          basis={basis}
          initial={active?.config}
          onConfirm={handleConfirm}
          onCancel={() => setDialogOpen(false)}
        />
      )}
    </>
  );
}
