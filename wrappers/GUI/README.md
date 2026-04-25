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

Signing is wired into the workflow but inactive by default. It activates
when the GitHub repo defines the relevant secrets.

**macOS — Developer ID + notarization.** Add these repo secrets:

| Secret                       | Source                                                       |
|------------------------------|--------------------------------------------------------------|
| `APPLE_CERTIFICATE`          | base64 of the exported `.p12` Developer ID Application cert  |
| `APPLE_CERTIFICATE_PASSWORD` | password used during the `.p12` export                       |
| `APPLE_SIGNING_IDENTITY`     | e.g. `Developer ID Application: Your Name (TEAMID)`          |

Plus one of the two notarization flows:

- App-password: `APPLE_ID`, `APPLE_PASSWORD` (an app-specific password
  generated at appleid.apple.com), `APPLE_TEAM_ID`.
- API key (preferred for CI): `APPLE_API_ISSUER`, `APPLE_API_KEY`,
  `APPLE_API_KEY_PATH` (path to the `.p8` written by the workflow).

With those set, `tauri-action` will codesign the `.app`, staple the
notarization ticket, and produce a `.dmg` that opens cleanly without the
`xattr` workaround.

**Windows — code signing.** Not yet wired in. The intended path for OSS
projects is either [SignPath.io](https://signpath.io/) (free for qualifying
OSS) or [Azure Trusted Signing](https://learn.microsoft.com/azure/trusted-signing/)
(pay-per-signature, no HSM). Both can be invoked as a post-build step
against the `.msi` and `.exe` produced by `tauri build`.

**Linux** — sign the release artifacts with GPG out-of-band; Flathub builds
sign their own.

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
