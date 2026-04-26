import { COOLPROP_VERSION, COOLPROP_GIT_HASH, NOTICES } from "../generated/notices";

interface Props {
  onClose: () => void;
}

export default function AboutModal({ onClose }: Props) {
  return (
    <div
      className="modal-overlay"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div className="modal-card about-card">
        <div className="modal-title">About CoolProp Desktop</div>

        <div className="modal-body">
          <p>
            <strong>CoolProp:</strong>{" "}
            <span className="about-version">{COOLPROP_VERSION}</span>{" "}
            <span className="about-hash">(build {COOLPROP_GIT_HASH})</span>
          </p>
          <p>
            Built with Tauri, React, and{" "}
            <a
              href="https://github.com/CoolProp/CoolProp"
              target="_blank"
              rel="noreferrer"
            >
              CoolProp
            </a>. Released under the MIT license.
          </p>

          <details className="about-notices">
            <summary>Third-party notices</summary>
            <pre className="about-notices-pre">{NOTICES}</pre>
          </details>
        </div>

        <div className="modal-footer">
          <button className="primary" onClick={onClose}>Close</button>
        </div>
      </div>
    </div>
  );
}
