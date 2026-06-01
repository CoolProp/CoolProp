import { useState } from "react";
import { SPONSOR_URL } from "../constants";
import { APP_VERSION } from "../generated/notices";

interface Props {
  /** Semver string; the integer before the first "." is the major version.
   *  Defaults to the GUI bundle version generated at build time. */
  version?: string;
}

function majorOf(version: string): string {
  const major = version.split(".")[0];
  return /^\d+$/.test(major) ? major : "0";
}

export default function SponsorSplash({ version = APP_VERSION }: Props) {
  const storageKey = `coolprop.sponsorSplash.seen.major.${majorOf(version)}`;

  const [open, setOpen] = useState<boolean>(() => {
    try {
      return localStorage.getItem(storageKey) !== "1";
    } catch {
      return true;
    }
  });

  function markSeen() {
    try {
      localStorage.setItem(storageKey, "1");
    } catch {
      /* ignore unavailable storage */
    }
  }

  function dismiss() {
    markSeen();
    setOpen(false);
  }

  if (!open) return null;

  return (
    <div
      className="modal-overlay"
      onClick={(e) => {
        if (e.target === e.currentTarget) dismiss();
      }}
    >
      <div className="modal-card sponsor-splash-card">
        <div className="sponsor-splash-heart">💚</div>
        <div className="modal-title">Enjoying CoolProp?</div>
        <div className="modal-body">
          <p>
            CoolProp is free and open-source, maintained on volunteer time. If it
            helps your work, please consider sponsoring its development.
          </p>
        </div>
        <div className="modal-footer">
          <button onClick={dismiss}>Maybe later</button>
          <a
            className="primary sponsor-splash-cta"
            href={SPONSOR_URL}
            target="_blank"
            rel="noreferrer"
            onClick={markSeen}
          >
            Sponsor on GitHub
          </a>
        </div>
      </div>
    </div>
  );
}
