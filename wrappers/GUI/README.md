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

Produces a signed-ready bundle for the host platform under
`src-tauri/target/release/bundle/`:

- macOS: `.app` and `.dmg`
- Windows: `.msi` and `.exe`
- Linux: `.deb`, `.rpm`, and `.AppImage`

Code signing and notarization are not configured here — see the Tauri docs
for `tauri.conf.json > bundle > macOS` (entitlements, signing identity) and
the Windows code-signing options.

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
