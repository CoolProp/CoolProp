import { useState, useEffect, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";

export interface AbstractStateHandle {
  id: number | null;
  ready: boolean;
  error: string | null;
  update: (inputPair: string, v1: number, v2: number) => Promise<void>;
  get: (param: string) => Promise<number>;
}

export function useAbstractState(
  backend: string,
  fluid: string
): AbstractStateHandle {
  const [id, setId] = useState<number | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!fluid) return;

    let stateId: number | null = null;
    let cancelled = false;

    invoke<number>("create_state", { backend, fluid })
      .then((newId) => {
        if (cancelled) {
          invoke("free_state", { id: newId }).catch(console.error);
        } else {
          stateId = newId;
          setId(newId);
          setError(null);
        }
      })
      .catch((err) => {
        if (!cancelled) setError(String(err));
      });

    return () => {
      cancelled = true;
      if (stateId !== null) {
        invoke("free_state", { id: stateId }).catch(console.error);
        setId(null);
      }
    };
  }, [backend, fluid]);

  const update = useCallback(
    async (inputPair: string, v1: number, v2: number): Promise<void> => {
      if (id === null) throw new Error("AbstractState not ready");
      await invoke("update_state", { id, inputPair, v1, v2 });
    },
    [id]
  );

  const get = useCallback(
    async (param: string): Promise<number> => {
      if (id === null) throw new Error("AbstractState not ready");
      return invoke<number>("get_property", { id, param });
    },
    [id]
  );

  return { id, ready: id !== null, error, update, get };
}
