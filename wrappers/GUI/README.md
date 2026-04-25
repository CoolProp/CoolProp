# CoolProp Desktop GUI

A cross-platform desktop application for interactive CoolProp calculations,
built with [Tauri 2](https://tauri.app/) (Rust backend + native WebView) and
React/TypeScript on the frontend. CoolProp itself is compiled as a static
library and linked into the app bundle.

Tabs:

- **Calculator** — `AbstractState`-based property lookup at a single point.
- **Saturation** — multi-panel saturation tables (T or P sweep) per fluid.
- **Humid Air** — `HAPropsSI` calculator with an interactive psychrometric
  chart and optional isolines (wet-bulb / volume / dew-point) through the
  computed state point.

## Prerequisites

You need three toolchains installed:

1. **Node.js** ≥ 18 (for the Vite-based frontend and the Tauri CLI).
2. **Rust** stable (≥ 1.78). Install via [rustup](https://rustup.rs/).
3. **CMake** ≥ 3.16 and a C++17 compiler (used by `build.rs` to compile
   CoolProp as a static library).

Platform-specific WebView dependencies follow the [Tauri prerequisites
guide](https://tauri.app/start/prerequisites/):

- **macOS** — Xcode Command Line Tools (`xcode-select --install`). The build
  pins `CMAKE_OSX_DEPLOYMENT_TARGET=10.15` because CoolProp's
  `CPfilepaths.cpp` uses `std::filesystem`.
- **Windows** — MSVC (Visual Studio Build Tools 2019+) and the WebView2
  runtime (preinstalled on Windows 11; install separately on Windows 10).
- **Linux** — `webkit2gtk-4.1`, `libssl-dev`, `libgtk-3-dev`, `libayatana-appindicator3-dev`,
  and `librsvg2-dev`. See the Tauri Linux prerequisites for the exact package
  names on your distro.

## Install

From this directory (`wrappers/GUI/`):

```sh
npm install
```

This pulls down the Tauri CLI and the React/Plotly frontend dependencies. No
Rust crates are fetched yet — that happens on first `tauri` invocation.

## Run in dev mode

```sh
npm run tauri dev
```

What this does, in order:

1. Starts the Vite dev server on `http://localhost:1420` for the frontend.
2. Runs `cargo run` against `src-tauri/`, which triggers `build.rs` to
   configure and build CoolProp as a static library
   (`-DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_LIB`) under
   `src-tauri/target/debug/build/coolprop-gui-*/out/build/`.
3. Links the resulting `libCoolProp.a` (or `CoolProp.lib`) into the Rust
   crate and launches the Tauri app pointed at the Vite dev URL.

The first build takes several minutes because CoolProp's translation units
are not pre-built; subsequent rebuilds are incremental and fast (typically
~10 s of relink time after a frontend or Rust-side change). Frontend edits
are picked up by Vite HMR without a Rust rebuild.

## Build a release bundle

```sh
npm run tauri build
```

Produces a bundle for the host platform under
`src-tauri/target/release/bundle/`:

- macOS: `.app` and `.dmg`
- Windows: `.msi` and `.exe`
- Linux: `.deb`, `.rpm`, and `.AppImage`

CI (`.github/workflows/gui_builder.yml`) builds these bundles on every push
to GUI-affecting paths and uploads them as artifacts.

### Running an unsigned macOS bundle (CI artifact)

CI artifacts are not currently signed or notarized. macOS Gatekeeper will
block them on first launch ("App is damaged" or "from an unidentified
developer"). To launch:

```sh
xattr -dr com.apple.quarantine /path/to/CoolProp.app
open /path/to/CoolProp.app
```

The same applies if you run a `.app` extracted from a CI-built `.dmg` —
strip quarantine on the installed copy after dragging it to Applications.

### Code signing (production)

Both macOS and Windows signing are wired into the workflow but stay inactive
until the relevant repo secrets are defined. With nothing configured, the
build still produces functional unsigned bundles.

#### macOS — Developer ID + notarization

One-time setup:

1. **Get an Apple Developer ID.** Sign up at
   [developer.apple.com](https://developer.apple.com/) ($99/yr at the time
   of writing). The membership doesn't have to be in the project's name —
   any maintainer's account works.
2. **Create a Developer ID Application certificate.** In Xcode → Settings →
   Accounts → "Manage Certificates…" → `+` → *Developer ID Application*. (Or
   use the developer portal directly: Certificates, Identifiers & Profiles
   → Certificates → `+` → *Developer ID Application*.)
3. **Export the cert as `.p12`.** In Keychain Access, right-click the cert
   → *Export*, choose Personal Information Exchange (.p12), set a password.
4. **Base64-encode it for GitHub.**
   `base64 -i Cert.p12 | pbcopy` — paste into the GitHub secret value.
5. **Generate a notarytool API key (recommended over app-passwords).**
   App Store Connect → Users and Access → Integrations → App Store Connect
   API → `+` → grant "Developer" role. Download the `.p8` once (you can't
   get it back) and note the Issuer ID and Key ID. Base64 the `.p8`:
   `base64 -i AuthKey_XXXXXXXXXX.p8 | pbcopy`.
6. **Add the GitHub secrets** (Repo Settings → Secrets and variables →
   Actions → New repository secret):

| Secret                       | Value                                                              |
|------------------------------|--------------------------------------------------------------------|
| `APPLE_CERTIFICATE`          | base64 of the exported `.p12`                                      |
| `APPLE_CERTIFICATE_PASSWORD` | password used during `.p12` export                                 |
| `APPLE_SIGNING_IDENTITY`     | full identity string, e.g. `Developer ID Application: J. Doe (ABCDE12345)` |
| `APPLE_TEAM_ID`              | Team ID, e.g. `ABCDE12345`                                         |

Then **one** of the two notarization flows:

- **API key flow (recommended):** `APPLE_API_KEY` (base64 of the `.p8`),
  `APPLE_API_ISSUER` (Issuer ID UUID), `APPLE_API_KEY_PATH` (path the
  workflow should write the decoded `.p8` to, e.g. `~/AuthKey.p8`).
- **App-password flow:** `APPLE_ID` (Apple-ID email), `APPLE_PASSWORD`
  (app-specific password from appleid.apple.com → Sign-In and Security).

Once those exist, the workflow's *Configure macOS signing* step forwards
them to `tauri-action`, which signs the `.app`, staples the notarization
ticket, and produces a `.dmg` that opens cleanly without the `xattr`
workaround.

#### Windows — Authenticode signing

Two paths supported:

- **Direct certificate (any Authenticode cert).** Set:

| Secret                          | Value                                       |
|---------------------------------|---------------------------------------------|
| `WINDOWS_CERTIFICATE`           | base64 of a `.pfx` containing the cert + key |
| `WINDOWS_CERTIFICATE_PASSWORD`  | password used during the `.pfx` export       |

  Encode with PowerShell: `[Convert]::ToBase64String([IO.File]::ReadAllBytes('cert.pfx')) | Set-Clipboard`.
  The workflow's *Sign Windows .exe and .msi* step decodes the cert into a
  temp `.pfx`, locates `signtool.exe` from the bundled Windows Kits, and
  signs the `target/release/*.exe`, `bundle/msi/*.msi`, and any
  `bundle/nsis/*.exe` produced by `tauri build`. Timestamp server is
  DigiCert's free endpoint.

- **SignPath OSS (recommended for OSS projects without a paid cert).**
  Free for qualifying OSS via [signpath.io](https://signpath.io/products/foundation).
  Replace the *Sign Windows .exe and .msi* step with the
  `signpath/github-action-submit-signing-request` action — see SignPath's
  [GitHub action docs](https://about.signpath.io/documentation/build-system-integration#github-actions)
  for the exact YAML; you'll set `SIGNPATH_*` secrets instead of
  `WINDOWS_CERTIFICATE`.

- **Azure Trusted Signing** is also a fit (pay-per-signature, no HSM
  hardware) — analogous post-build step.

#### Linux

GPG-sign release artifacts out-of-band; Flathub builds sign their own.

### Releases

Tag the GUI separately from the CoolProp library. The CI workflow watches
for tags matching `gui-v*` and creates a draft GitHub Release with all
platform bundles attached:

```sh
# Bump versions in three places (must match):
#   wrappers/GUI/package.json                 .version
#   wrappers/GUI/src-tauri/Cargo.toml         [package].version
#   wrappers/GUI/src-tauri/tauri.conf.json    .version

git tag gui-v0.1.0
git push origin gui-v0.1.0
```

After the workflow finishes, find the draft release on GitHub, edit the
release notes, and publish. macOS bundles in the release are signed and
notarized only if the macOS signing secrets are configured (see above);
Windows bundles only if the Windows signing secrets are.

## Regenerating the app icon

The app icon is the canonical CoolProp water phase-diagram logo. To
regenerate it from scratch (requires the `CoolProp` Python package,
`matplotlib`, and `scipy`):

```sh
python scripts/regenerate-icons.py
```

This renders the logo at 1024×1024 and runs `npx tauri icon` to produce the
full per-platform icon set under `src-tauri/icons/`.

## Troubleshooting

- **`failed to run custom build command for coolprop-gui`** — typically
  means CMake or a C++ toolchain is missing. Verify with `cmake --version`
  and that a working C++ compiler is on PATH.
- **macOS: `error: target ... requires a deployment target of 10.15+`** —
  the build sets `CMAKE_OSX_DEPLOYMENT_TARGET=10.15` automatically; if you
  see this in your own scripts, set the same env var.
- **Linux WebKit errors at runtime** — install the platform WebKit packages
  listed above; Tauri will not bundle them.
- **Window doesn't appear / process exits immediately on macOS** — launch
  from a regular terminal (not a backgrounded shell). Tauri/NSApp needs an
  interactive session to attach the window to.
