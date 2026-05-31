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

### Downloading a build

For end users the intended channel is a **GitHub Release** (created by
pushing a `gui-v*` tag): each installer is attached as a *direct* download —
the macOS `.dmg`, the Windows `.msi` / `-setup.exe`, the Linux `.deb` /
`.rpm` / `.AppImage`.

The per-build **CI artifacts** (`coolprop-gui-<os>`) are for testing. GitHub
always wraps an artifact in a `.zip`, so you download e.g.
`coolprop-gui-macos-latest.zip` and unzip it. The artifact is curated to the
installers only:

- **macOS** — open `dmg/CoolProp_*.dmg`, then drag the app to Applications.
  Do **not** try to launch a `.app` pulled straight out of a zip: a macOS app
  bundle loses the bits that make it launchable through a plain zip, so the
  `.dmg` is the correct, self-contained format. (`*.app.tar.gz` is the
  auto-updater payload, not a download.)
- **Windows** — run `msi/CoolProp_*.msi` (or `nsis/CoolProp_*-setup.exe`).
- **Linux** — install `deb/*.deb` or `rpm/*.rpm`, or run `appimage/*.AppImage`.

### Gatekeeper / SmartScreen on unsigned builds

A build produced without active signing is **unsigned**, and macOS Gatekeeper
blocks an unsigned `.dmg`/`.app` on first launch ("App is damaged" / "from an
unidentified developer"). Clear quarantine to run it:

```sh
xattr -dr com.apple.quarantine /Applications/CoolProp.app
open /Applications/CoolProp.app
```

Signed + notarized builds (release tags) open with no workaround. See **Code
signing (production)** below for what enables signing.

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
5. **Generate an App Store Connect API key (preferred over app-passwords —
   it isn't tied to a person and won't break on 2FA / password changes).**
   App Store Connect → Users and Access → Integrations → App Store Connect
   API → `+` → grant the **Developer** role. Download the `.p8` once (you
   can't get it back), and note the **Issuer ID** (the UUID above the keys
   table) and the **Key ID** (the short identifier in the Key ID column).
   Base64 the `.p8`: `base64 -i AuthKey_XXXXXXXXXX.p8 | pbcopy`.
6. **Add the GitHub secrets** (Repo Settings → Secrets and variables →
   Actions → New repository secret):

| Secret                       | Value                                                              |
|------------------------------|--------------------------------------------------------------------|
| `APPLE_CERTIFICATE`          | base64 of the exported `.p12`                                      |
| `APPLE_CERTIFICATE_PASSWORD` | password used during `.p12` export                                 |
| `APPLE_SIGNING_IDENTITY`     | full identity string, e.g. `Developer ID Application: J. Doe (ABCDE12345)` |
| `APPLE_TEAM_ID`              | Team ID, e.g. `ABCDE12345`                                         |

Then **one** of the two notarization flows:

- **API-key flow (preferred):** three secrets —

  | Secret                 | Value                                                       |
  |------------------------|-------------------------------------------------------------|
  | `APPLE_API_ISSUER`     | Issuer ID (UUID above the keys table)                       |
  | `APPLE_API_KEY`        | **Key ID** — the short identifier, *not* the key body       |
  | `APPLE_API_KEY_BASE64` | base64 of the downloaded `AuthKey_*.p8`                      |

  The *Configure macOS signing* step decodes `APPLE_API_KEY_BASE64` to a
  file on the runner and sets `APPLE_API_KEY_PATH` itself — you do **not**
  add `APPLE_API_KEY_PATH` as a secret.

- **App-specific-password flow (alternative):** `APPLE_ID` (Apple-ID email)
  and `APPLE_PASSWORD` (app-specific password from appleid.apple.com →
  Sign-In and Security). Uses `APPLE_TEAM_ID` from the table above.

Once those exist, the workflow's *Configure macOS signing* step forwards
them to `tauri-action`, which signs the `.app`, staples the notarization
ticket, and produces a `.dmg` that opens cleanly without the `xattr`
workaround.

#### Windows — Authenticode signing via SignPath

The workflow signs Windows installers through **SignPath OSS** (free for
qualifying open-source projects via
[signpath.io](https://signpath.io/products/foundation)). SignPath keeps the
certificate in its HSM and signs an *already-uploaded* GitHub artifact, so the
*Sign Windows installers via SignPath* step uploads the unsigned bundle,
submits it via `signpath/github-action-submit-signing-request`, and overlays
the signed installers back onto the bundle dir before the Windows smoke test
and the release upload. It activates when the `SIGNPATH_API_TOKEN` secret is
present, and no-ops to an unsigned build when it isn't (forks, secret unset).

**Signing is opt-in** to avoid a SignPath email on every GUI build (SignPath
has no per-state notification filter). A signing request fires only when:

- the triggering commit's message contains the marker **`[gui-sign]`**, or
- the build is for a **`gui-v*` release tag** (releases are always signed).

Routine GUI commits therefore build unsigned and send no email; add
`[gui-sign]` to a commit message (or push a release tag) when you want to
exercise/validate SignPath. macOS notarization is *not* gated this way — it
runs on every build whenever the Apple secrets are present.

| Setting                | Where                | Value                                              |
|------------------------|----------------------|----------------------------------------------------|
| `SIGNPATH_API_TOKEN`   | repo **secret**      | SignPath REST API token for the CI user            |
| organization-id        | hard-coded in YAML   | CoolProp's SignPath org GUID                        |
| project-slug           | hard-coded in YAML   | `coolprop`                                          |
| `SIGNPATH_POLICY_SLUG` | repo **variable**    | signing policy slug; defaults to `test-signing`     |

The policy defaults to `test-signing`, which exercises the full pipeline using
SignPath's **untrusted test certificate** — this validates the wiring but will
*not* clear SmartScreen. SignPath Foundation issues the trusted release
certificate only after manually approving the project; once that lands, set the
repo variable `SIGNPATH_POLICY_SLUG=release-signing` (Settings → Secrets and
variables → Actions → Variables) to go live with no code change.

**Signing the installed app binary too.** The CI step uploads the whole
bundle (the `.msi` and the NSIS `-setup.exe`), so the *installed* executable
(`coolprop-gui.exe`, embedded inside those installers) can be signed in the
same request — no CI change needed. Enable it portal-side by ticking **"Sign
nested files"** on the SignPath *artifact configuration*; SignPath then signs
both the outer installer and the PE files it contains in a single pass. The
workflow's post-overlay `Get-AuthenticodeSignature` check only verifies the
outer installers; trust SignPath's nested-signing for the embedded binary.

**Alternatives** (would require swapping the signing step): a direct
Authenticode `.pfx` via `signtool` (set base64 `WINDOWS_CERTIFICATE` +
`WINDOWS_CERTIFICATE_PASSWORD` secrets), or **Azure Trusted Signing**
(pay-per-signature, no HSM hardware).

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

### Auto-updater

The app uses `@tauri-apps/plugin-updater` to check for new releases on
launch (5-second delay after window appears) and prompt the user to
install. Manifest URL: the GitHub Release latest.json asset.

Activating updates requires a one-time signing-keypair setup:

1. **Generate the keypair** on a workstation (NOT in CI):
   ```sh
   cd wrappers/GUI
   npx @tauri-apps/cli signer generate -w /path/to/coolprop-update.key
   ```
   This emits a private key file and prints the matching public key.

2. **Replace the placeholder pubkey** in
   `wrappers/GUI/src-tauri/tauri.conf.json` (`plugins.updater.pubkey`)
   with the printed public key, **and** add `"createUpdaterArtifacts":
   true` to the `bundle` block in the same file. Commit that change —
   the public key is safe to publish. The flag is left off by default so
   builds don't fail when the signing secrets aren't configured yet.

3. **Add the GitHub secrets** so CI can sign release artifacts:

| Secret                                | Value                                           |
|---------------------------------------|-------------------------------------------------|
| `TAURI_SIGNING_PRIVATE_KEY`           | full contents of the private key file           |
| `TAURI_SIGNING_PRIVATE_KEY_PASSWORD`  | password used during keypair generation         |

4. **Back up the private key** offline. Losing it means future releases
   can't be auto-installed by users running existing builds; they'd
   have to manually download the new version.

With those in place, every `gui-v*` tag pushes both the bundles and a
signed `latest.json` to the GitHub Release. Older builds will detect the
release, verify the signature against the embedded public key, download,
and self-update.

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
