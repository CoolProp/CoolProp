# Design: GitHub Sponsors promotion across CoolProp surfaces

**Date:** 2026-05-31
**Status:** Approved (visual mockups reviewed and locked in)

## Goal

The CoolProp GitHub organization is now eligible for GitHub Sponsors. Surface
that fact tastefully but persistently across every place a user encounters the
project: the repository, the documentation/website, and the desktop GUI. The
tone is *tasteful everywhere* — visible and repeated, but never a blocking
interruption.

## Single source of truth

All sponsor links point to:

```
https://github.com/sponsors/CoolProp
```

Define this once per surface (a constant in the GUI, a substitution/variable in
the docs) so it can be changed in one place.

## Surfaces

Delivered as **one PR** covering all four surfaces.

### 1. Repository — `.github/FUNDING.yml`

New file enabling the green **Sponsor** button on the repo page and sidebar:

```yaml
github: [CoolProp]
```

Zero runtime risk; GitHub validates the file automatically.

### 2 & 3. Website / docs (`Web/`, Sphinx + pydata-sphinx-theme)

The `Web/` directory is the Sphinx source for coolprop.org. Three additions:

**(a) Persistent nav link.** Add one `<li>` to
`Web/_templates/navbar-center.html` (the hand-written nav `<ul>`):

```html
<li class="nav-item">
  <a class="nav-link" href="https://github.com/sponsors/CoolProp"
     target="_blank" rel="noopener">💚 Sponsor</a>
</li>
```

**(b) Dismissable site-wide top banner.** A self-contained CSS + JS pair added
to `Web/_static/`:

- `sponsor-banner.css` — styles a slim green banner pinned at the top of the
  page body, plus a × dismiss control.
- `sponsor-banner.js` — on `DOMContentLoaded`, if the dismissal key is absent
  from `localStorage`, inject the banner at the top of the content. The ×
  click sets the key so the banner does not reappear on subsequent page loads.

Register both in `Web/conf.py`:

```python
html_css_files = ['sponsor-banner.css']           # new
html_js_files  = [('3Dmol-min.js', {'priority': 450}),
                  'sponsor-banner.js']             # append
```

Banner copy:

> 💚 **CoolProp is free & open-source**, maintained on volunteer time. If it
> helps your work, please consider sponsoring → *Sponsor*

The `localStorage` key is versioned, e.g. `coolprop.sponsorBanner.dismissed.v1`,
so the banner can be intentionally re-shown later by bumping the suffix.

**(c) Dedicated page.** New `Web/online/sponsor.rst` titled "Support CoolProp",
explaining what sponsorship funds (ongoing maintenance, REFPROP-grade
validation, build/CI infrastructure, the desktop GUI) with a prominent link to
the Sponsors page. Add it to the `Web/online/index.rst` toctree and cross-link
it from the **Help** section of `Web/index.rst`.

### 4. Desktop GUI (`wrappers/GUI/`, Tauri 2 + React 18)

External links in the GUI already use plain
`<a target="_blank" rel="noreferrer">` (see `AboutModal.tsx`), which Tauri opens
in the system browser. **No new Tauri plugin, capability, or dependency is
required.**

**(a) Header link.** Add a `💚 Sponsor` button to the `header-right` group in
`src/App.tsx`, between the Mass/Molar segmented control and the **About**
button. Styled to match the existing header controls.

**(b) About modal.** Add a "Support CoolProp 💚" link to
`src/components/AboutModal.tsx`, near the existing CoolProp/MIT-license line.

**(c) One-time startup splash.** New `src/components/SponsorSplash.tsx`, reusing
the existing `modal-overlay` / `modal-card` styles. Behavior:

- On app mount, read `localStorage` key `coolprop.sponsorSplash.seen.v1`.
- If unset, render the splash; set the key when the user dismisses it (either
  button), so it shows **once per install**.
- Buttons: **Maybe later** (ghost, closes) and **Sponsor on GitHub** (primary,
  opens the Sponsors URL via the same `<a target="_blank">` pattern, then
  closes).

Splash copy:

> 💚 **Enjoying CoolProp?**
> CoolProp is free and open-source, maintained on volunteer time. If it helps
> your work, please consider sponsoring its development.

Wire `<SponsorSplash />` into `App.tsx` alongside the existing
`<UpdateChecker />`.

Define the Sponsors URL as a shared constant (e.g. in a small
`src/constants.ts`) referenced by the header link, the About modal, and the
splash.

## Visual design (approved)

- Accent color: GitHub-Sponsors green (`#2ea043`), light fill `#e8f6ec`.
- Wording: heart emoji 💚 + the word "Sponsor".
- GUI splash and web banner are both **dismissable and remembered**; the web
  banner additionally persists a small always-visible nav link, and the GUI a
  small always-visible header chip.

## Testing

- **FUNDING.yml** — validated automatically by GitHub.
- **Web** — build the Sphinx site locally (`Web/`); confirm: banner renders,
  × dismiss persists across page loads, nav link present, "Support CoolProp"
  page builds and is reachable from the toctree and the Help section.
- **GUI** — vitest unit tests for `SponsorSplash`:
  - renders when the `seen` key is absent;
  - does not render when the key is set;
  - dismissing (either button) sets the key.
  Confirm the existing `App` / `PropertyCalculator` tests stay green.

## Out of scope / YAGNI

- No blocking/modal interrupts on the website (banner only).
- No recurring or rate-limited re-prompting; both the GUI splash and web banner
  are strictly once-until-dismissed.
- No analytics/telemetry on sponsor clicks.
- No changes to the legacy Python/wxPython GUI under `wrappers/Python`.
