import { useEffect, useState } from "react";
import { check, type Update } from "@tauri-apps/plugin-updater";
import { relaunch } from "@tauri-apps/plugin-process";

type Status = "idle" | "downloading" | "ready" | "error";

export default function UpdateChecker() {
  const [update, setUpdate] = useState<Update | null>(null);
  const [status, setStatus] = useState<Status>("idle");
  const [progress, setProgress] = useState<{ done: number; total: number } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [dismissed, setDismissed] = useState(false);

  useEffect(() => {
    let cancelled = false;
    const id = setTimeout(() => {
      // Run a short delay after launch so the splash settles before any
      // dialog. Failures are silent (network down, dev build, etc.).
      check()
        .then((u) => {
          if (!cancelled && u && u.available) setUpdate(u);
        })
        .catch(() => { /* swallow */ });
    }, 5000);
    return () => {
      cancelled = true;
      clearTimeout(id);
    };
  }, []);

  if (!update || dismissed) return null;

  const install = async () => {
    setStatus("downloading");
    setError(null);
    try {
      await update.downloadAndInstall((event) => {
        switch (event.event) {
          case "Started":
            setProgress({ done: 0, total: event.data.contentLength ?? 0 });
            break;
          case "Progress":
            setProgress((p) =>
              p ? { ...p, done: p.done + event.data.chunkLength } : p
            );
            break;
          case "Finished":
            setStatus("ready");
            break;
        }
      });
      await relaunch();
    } catch (e) {
      setStatus("error");
      setError(String(e));
    }
  };

  const pct =
    progress && progress.total > 0
      ? Math.round((progress.done / progress.total) * 100)
      : null;

  return (
    <div className="update-banner" role="status">
      <div className="update-banner-text">
        {status === "idle" && (
          <>Update available: <strong>{update.version}</strong></>
        )}
        {status === "downloading" && (
          <>Downloading {update.version}{pct !== null ? ` — ${pct}%` : "…"}</>
        )}
        {status === "ready" && <>Restarting to apply {update.version}…</>}
        {status === "error" && (
          <>Update failed: <span className="update-error">{error}</span></>
        )}
      </div>
      <div className="update-banner-actions">
        {status === "idle" && (
          <>
            <button className="primary" onClick={install}>Install</button>
            <button className="btn-secondary" onClick={() => setDismissed(true)}>
              Later
            </button>
          </>
        )}
        {status === "error" && (
          <button className="btn-secondary" onClick={() => setDismissed(true)}>
            Dismiss
          </button>
        )}
      </div>
    </div>
  );
}
