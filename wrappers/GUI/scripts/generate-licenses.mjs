#!/usr/bin/env node
// Generates `wrappers/GUI/src/generated/notices.ts` with:
//   - third-party license notices (cargo metadata + license-checker)
//   - the CoolProp version pulled from CMakeLists.txt
//   - the current git short hash
//   - APP_VERSION (kept in lockstep with CoolProp's MAJOR.MINOR.PATCH)
//
// Also writes the synced version into wrappers/GUI/package.json,
// wrappers/GUI/src-tauri/Cargo.toml, and wrappers/GUI/src-tauri/tauri.conf.json
// so the bundle metadata (.dmg/.msi/.deb version strings) matches what the
// About modal displays.
//
// Run via `npm run build:licenses` (chained from prebuild/predev/pretest).

import { execSync } from "node:child_process";
import { writeFileSync, mkdirSync, readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const guiDir = resolve(here, "..");
const tauriDir = resolve(guiDir, "src-tauri");
const repoRoot = resolve(guiDir, "../..");
const outFile = resolve(guiDir, "src/generated/notices.ts");

// ── CoolProp version + git hash ─────────────────────────────────────────────
const cmakeText = readFileSync(resolve(repoRoot, "CMakeLists.txt"), "utf8");
const grab = (key) => cmakeText.match(new RegExp(`set\\(${key}\\s+([^\\)\\s]+)\\)`))?.[1];
const major = grab("COOLPROP_VERSION_MAJOR");
const minor = grab("COOLPROP_VERSION_MINOR");
const patch = grab("COOLPROP_VERSION_PATCH");
const revision = grab("COOLPROP_VERSION_REVISION") ?? "";
if (!major || !minor || !patch) {
  throw new Error("Failed to parse COOLPROP_VERSION_{MAJOR,MINOR,PATCH} from CMakeLists.txt");
}
const semver = `${major}.${minor}.${patch}`;
const fullCoolPropVersion = revision && revision !== "" ? `${semver}-${revision}` : semver;

let gitHash = "unknown";
try {
  gitHash = execSync("git rev-parse --short HEAD", { cwd: repoRoot, encoding: "utf8" }).trim();
} catch { /* not a git checkout — leave as 'unknown' */ }

console.log(`CoolProp version: ${fullCoolPropVersion} (${gitHash})`);

// ── Sync GUI bundle version into package.json / Cargo.toml / tauri.conf.json
function syncJsonVersion(path, newVersion) {
  const json = JSON.parse(readFileSync(path, "utf8"));
  if (json.version === newVersion) return false;
  json.version = newVersion;
  writeFileSync(path, JSON.stringify(json, null, 2) + "\n");
  return true;
}
function syncCargoTomlVersion(path, newVersion) {
  const orig = readFileSync(path, "utf8");
  // Replace only the [package] section's version line; assume it's the first
  // `version = "x.y.z"` after `[package]`.
  const updated = orig.replace(
    /(\[package\][\s\S]*?version\s*=\s*)"[^"]*"/,
    `$1"${newVersion}"`
  );
  if (updated === orig) return false;
  writeFileSync(path, updated);
  return true;
}

const synced = [
  syncJsonVersion(resolve(guiDir, "package.json"), semver),
  syncJsonVersion(resolve(tauriDir, "tauri.conf.json"), semver),
  syncCargoTomlVersion(resolve(tauriDir, "Cargo.toml"), semver),
];
if (synced.some(Boolean)) {
  console.log(`Bundle version synced to ${semver} in ${synced.filter(Boolean).length} file(s)`);
}

const pkg = JSON.parse(readFileSync(resolve(guiDir, "package.json"), "utf8"));

// ── Collect Rust + npm dependency metadata ──────────────────────────────────
console.log("Collecting Rust crate metadata via `cargo metadata`…");
const cargoMeta = JSON.parse(
  execSync("cargo metadata --format-version 1 --frozen --offline 2>/dev/null || cargo metadata --format-version 1", {
    cwd: tauriDir,
    encoding: "utf8",
    maxBuffer: 64 * 1024 * 1024,
  })
);

const rustDeps = cargoMeta.packages
  .filter((p) => p.name !== "coolprop-gui")
  .map((p) => ({
    name: p.name,
    version: p.version,
    license: p.license || p.license_file || "(unspecified)",
    repository: p.repository || "",
  }))
  .sort((a, b) => a.name.localeCompare(b.name));

console.log(`  → ${rustDeps.length} Rust crates`);

console.log("Collecting npm package metadata via license-checker…");
let npmDeps = [];
try {
  const out = execSync("npx --no-install license-checker --production --json", {
    cwd: guiDir,
    encoding: "utf8",
    maxBuffer: 64 * 1024 * 1024,
  });
  npmDeps = Object.entries(JSON.parse(out))
    .filter(([k]) => !k.startsWith(`${pkg.name}@`))
    .map(([k, v]) => {
      const at = k.lastIndexOf("@");
      return {
        name: at > 0 ? k.slice(0, at) : k,
        version: at > 0 ? k.slice(at + 1) : "",
        license: Array.isArray(v.licenses) ? v.licenses.join(", ") : v.licenses || "(unspecified)",
        repository: v.repository || "",
      };
    })
    .sort((a, b) => a.name.localeCompare(b.name));
  console.log(`  → ${npmDeps.length} npm packages`);
} catch (e) {
  console.warn("  ! license-checker not available — npm section will be empty.");
  console.warn(`    (${e.message})`);
}

const fmt = (d) => {
  const repo = d.repository ? ` ([source](${d.repository}))` : "";
  return `- **${d.name}** ${d.version} — ${d.license}${repo}`;
};

const md = `# Third-party notices

CoolProp Desktop builds on the open-source components listed below. Full
license texts are available in each project's repository (linked where
metadata exposes the URL).

## Rust crates (${rustDeps.length})

${rustDeps.map(fmt).join("\n")}

## npm packages (${npmDeps.length})

${npmDeps.length ? npmDeps.map(fmt).join("\n") : "_(license-checker unavailable at build time)_"}
`;

mkdirSync(dirname(outFile), { recursive: true });
const content = `// AUTO-GENERATED by scripts/generate-licenses.mjs — do not edit.
export const APP_VERSION = ${JSON.stringify(semver)};
export const APP_NAME = ${JSON.stringify(pkg.name)};
export const COOLPROP_VERSION = ${JSON.stringify(fullCoolPropVersion)};
export const COOLPROP_GIT_HASH = ${JSON.stringify(gitHash)};
export const NOTICES = ${JSON.stringify(md)};
`;
writeFileSync(outFile, content);
console.log(`Wrote ${outFile} (${content.length} bytes)`);
