# GitHub Sponsors Promotion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote `https://github.com/sponsors/CoolProp` tastefully across four surfaces — the repo Sponsor button, the docs/website (nav link + dismissable banner + dedicated page), and the Tauri desktop GUI (header link + About entry + once-per-major-version splash).

**Architecture:** Four independent surfaces in one PR. The repo button is a static YAML file. The website changes are theme-agnostic static assets (CSS+JS injected via `conf.py`) plus reStructuredText. The GUI changes are React components reusing existing modal/header styles, with a single shared URL constant and `localStorage`-gated visibility — no new Tauri plugin, capability, or npm dependency (external links use the existing `<a target="_blank" rel="noreferrer">` pattern).

**Tech Stack:** GitHub FUNDING.yml; Sphinx + pydata-sphinx-theme (Python, RST, vanilla CSS/JS); Tauri 2 + React 18 + TypeScript + Vitest/Testing-Library.

**Spec:** `docs/superpowers/specs/2026-05-31-sponsor-promotion-design.md`

**Canonical URL (used everywhere):** `https://github.com/sponsors/CoolProp`

---

## File Structure

**Repo**
- Create: `.github/FUNDING.yml` — enables the repo Sponsor button.

**Website (`Web/`)**
- Create: `Web/_static/sponsor-banner.css` — banner + dismiss-button styling.
- Create: `Web/_static/sponsor-banner.js` — inject banner unless dismissed (localStorage).
- Create: `Web/sponsor.rst` — "Support CoolProp" page.
- Modify: `Web/_templates/navbar-center.html` — persistent nav Sponsor link.
- Modify: `Web/conf.py` — register the CSS/JS via `html_css_files`/`html_js_files`.
- Modify: `Web/contents.rst` — add `sponsor.rst` to the master toctree.
- Modify: `Web/index.rst` — cross-link the page from the Help section.

**GUI (`wrappers/GUI/`)**
- Create: `src/constants.ts` — `SPONSOR_URL` shared constant.
- Create: `src/components/SponsorSplash.tsx` — once-per-major-version splash.
- Create: `src/components/SponsorSplash.test.tsx` — vitest unit tests.
- Modify: `src/App.tsx` — header Sponsor link + mount `<SponsorSplash/>`.
- Modify: `src/components/AboutModal.tsx` — "Support CoolProp" link.
- Modify: `src/App.css` — `.sponsor-btn` / splash styles.

---

## Task 1: Repo Sponsor button (`.github/FUNDING.yml`)

**Files:**
- Create: `.github/FUNDING.yml`

- [ ] **Step 1: Create the FUNDING.yml file**

Create `.github/FUNDING.yml` with exactly:

```yaml
# Enables the "Sponsor" button on the repository.
# https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/displaying-a-sponsor-button-in-your-repository
github: [CoolProp]
```

- [ ] **Step 2: Verify YAML is well-formed**

Run: `python -c "import yaml,sys; print(yaml.safe_load(open('.github/FUNDING.yml')))"`
Expected: `{'github': ['CoolProp']}`

- [ ] **Step 3: Commit**

```bash
git add .github/FUNDING.yml
git commit -m "feat(repo): add FUNDING.yml for GitHub Sponsors button"
```

---

## Task 2: Website — dedicated "Support CoolProp" page

**Files:**
- Create: `Web/sponsor.rst`
- Modify: `Web/contents.rst`
- Modify: `Web/index.rst:49-57` (Help section)

- [ ] **Step 1: Create the page**

Create `Web/sponsor.rst`:

```rst
.. _sponsor:

****************
Support CoolProp
****************

CoolProp is free, open-source software, developed and maintained on
volunteer time. If CoolProp helps your work, please consider sponsoring its
continued development.

`💚 Sponsor CoolProp on GitHub <https://github.com/sponsors/CoolProp>`_

Your sponsorship directly supports:

* Ongoing maintenance and bug fixes
* High-accuracy validation against reference data (including NIST REFPROP)
* Build and continuous-integration infrastructure across many platforms
* Development of the cross-platform desktop GUI

Thank you for helping keep CoolProp free and open for everyone.
```

- [ ] **Step 2: Add the page to the master toctree**

In `Web/contents.rst`, add `sponsor.rst` to the second (visible) `toctree`. The
block currently reads:

```rst
.. toctree::

    general_information.rst
    online/index.rst
    fluid_properties/index.rst
    coolprop/index.rst
    develop/index.rst
    zbibliography
```

Change it to:

```rst
.. toctree::

    general_information.rst
    online/index.rst
    fluid_properties/index.rst
    coolprop/index.rst
    develop/index.rst
    sponsor.rst
    zbibliography
```

- [ ] **Step 3: Cross-link from the Help section**

In `Web/index.rst`, the Help section (lines 49-57) ends with the dev-docs
bullet. After the line:

```rst
* `Docs for development version of CoolProp <http://www.coolprop.org/dev/>`_
```

add a new bullet:

```rst
* (**Support the project**) :ref:`Sponsor CoolProp <sponsor>`
```

- [ ] **Step 4: Commit**

```bash
git add Web/sponsor.rst Web/contents.rst Web/index.rst
git commit -m "docs(web): add Support CoolProp sponsor page + cross-links"
```

---

## Task 3: Website — persistent nav Sponsor link

**Files:**
- Modify: `Web/_templates/navbar-center.html`

- [ ] **Step 1: Add the nav link**

`Web/_templates/navbar-center.html` is a hand-written `<ul class="navbar-nav">`.
Add a final `<li>` before the closing `</ul>`. The file currently ends:

```html
  <li class="nav-item"><a class="nav-link" href="{{ pathto('coolprop/wrappers/index') }}">Interfaces</a></li>
</ul>
```

Change it to:

```html
  <li class="nav-item"><a class="nav-link" href="{{ pathto('coolprop/wrappers/index') }}">Interfaces</a></li>
  <li class="nav-item"><a class="nav-link" href="https://github.com/sponsors/CoolProp" target="_blank" rel="noopener">💚 Sponsor</a></li>
</ul>
```

- [ ] **Step 2: Commit**

```bash
git add Web/_templates/navbar-center.html
git commit -m "docs(web): add persistent Sponsor link to nav bar"
```

---

## Task 4: Website — dismissable site-wide banner (CSS)

**Files:**
- Create: `Web/_static/sponsor-banner.css`

- [ ] **Step 1: Create the banner stylesheet**

Create `Web/_static/sponsor-banner.css`:

```css
/* Sponsor banner injected by sponsor-banner.js on every docs page. */
#coolprop-sponsor-banner {
  display: flex;
  align-items: center;
  gap: 10px;
  background: #e8f6ec;
  color: #1b5e2a;
  border-bottom: 1px solid #cfe8d6;
  padding: 9px 16px;
  font-size: 0.9rem;
  line-height: 1.4;
}
#coolprop-sponsor-banner .cp-sponsor-text { flex: 1 1 auto; }
#coolprop-sponsor-banner .cp-sponsor-text strong { color: #14532d; }
#coolprop-sponsor-banner a.cp-sponsor-cta {
  flex: 0 0 auto;
  background: #2ea043;
  color: #fff;
  text-decoration: none;
  font-weight: 600;
  border-radius: 6px;
  padding: 4px 12px;
  white-space: nowrap;
}
#coolprop-sponsor-banner button.cp-sponsor-dismiss {
  flex: 0 0 auto;
  background: none;
  border: none;
  color: #6b8f74;
  font-size: 1.15rem;
  line-height: 1;
  cursor: pointer;
  padding: 0 4px;
}
#coolprop-sponsor-banner button.cp-sponsor-dismiss:hover { color: #14532d; }
```

- [ ] **Step 2: Commit**

```bash
git add Web/_static/sponsor-banner.css
git commit -m "docs(web): add sponsor banner stylesheet"
```

---

## Task 5: Website — dismissable site-wide banner (JS)

**Files:**
- Create: `Web/_static/sponsor-banner.js`

- [ ] **Step 1: Create the banner script**

Create `Web/_static/sponsor-banner.js`:

```javascript
/* Injects a dismissable Sponsor banner at the top of every docs page.
   Dismissal is remembered in localStorage (versioned key) so it does not
   reappear on subsequent page loads. Theme-agnostic: prepends to <body>. */
(function () {
  "use strict";
  var DISMISS_KEY = "coolprop.sponsorBanner.dismissed.v1";
  var SPONSOR_URL = "https://github.com/sponsors/CoolProp";

  function dismissed() {
    try { return window.localStorage.getItem(DISMISS_KEY) === "1"; }
    catch (e) { return false; }
  }
  function remember() {
    try { window.localStorage.setItem(DISMISS_KEY, "1"); } catch (e) { /* ignore */ }
  }

  function build() {
    if (dismissed() || document.getElementById("coolprop-sponsor-banner")) return;

    var bar = document.createElement("div");
    bar.id = "coolprop-sponsor-banner";

    var text = document.createElement("span");
    text.className = "cp-sponsor-text";
    text.innerHTML =
      "💚 <strong>CoolProp is free &amp; open-source</strong>, " +
      "maintained on volunteer time. If it helps your work, please consider sponsoring →";

    var cta = document.createElement("a");
    cta.className = "cp-sponsor-cta";
    cta.href = SPONSOR_URL;
    cta.target = "_blank";
    cta.rel = "noopener";
    cta.textContent = "Sponsor";

    var x = document.createElement("button");
    x.className = "cp-sponsor-dismiss";
    x.type = "button";
    x.setAttribute("aria-label", "Dismiss sponsor banner");
    x.innerHTML = "&times;";
    x.addEventListener("click", function () {
      remember();
      if (bar.parentNode) bar.parentNode.removeChild(bar);
    });

    bar.appendChild(text);
    bar.appendChild(cta);
    bar.appendChild(x);
    document.body.insertBefore(bar, document.body.firstChild);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", build);
  } else {
    build();
  }
})();
```

- [ ] **Step 2: Commit**

```bash
git add Web/_static/sponsor-banner.js
git commit -m "docs(web): add sponsor banner injection script"
```

---

## Task 6: Website — register banner assets in conf.py

**Files:**
- Modify: `Web/conf.py:326` (`html_js_files`) and add `html_css_files`

- [ ] **Step 1: Register the CSS and JS**

In `Web/conf.py`, the `html_js_files` line currently reads (around line 326):

```python
html_js_files = [('3Dmol-min.js', {'priority': 450})]
```

Replace it with:

```python
html_js_files = [('3Dmol-min.js', {'priority': 450}), 'sponsor-banner.js']
html_css_files = ['sponsor-banner.css']
```

- [ ] **Step 2: Build the docs and verify the banner renders**

Run (from `Web/`):
```bash
cd Web && python -m sphinx -b html -j auto . _build/html 2>&1 | tail -20; cd ..
```
Expected: build completes; `Web/_build/html/_static/sponsor-banner.css` and
`sponsor-banner.js` exist, and `Web/_build/html/index.html` contains both
`sponsor-banner.css` and `sponsor-banner.js` references.

Verify references:
```bash
grep -o "sponsor-banner\.\(css\|js\)" Web/_build/html/index.html | sort -u
```
Expected output:
```
sponsor-banner.css
sponsor-banner.js
```

> Note: if the full Sphinx build is too slow locally (it executes notebooks),
> the grep over a single built page is the key assertion. The
> `html_css_files`/`html_js_files` mechanism is standard Sphinx and copies
> `_static` files automatically.

- [ ] **Step 3: Commit**

```bash
git add Web/conf.py
git commit -m "docs(web): register sponsor banner css/js in Sphinx config"
```

---

## Task 7: GUI — shared SPONSOR_URL constant

**Files:**
- Create: `wrappers/GUI/src/constants.ts`

- [ ] **Step 1: Create the constant**

Create `wrappers/GUI/src/constants.ts`:

```typescript
/** Canonical GitHub Sponsors URL for the CoolProp organization. */
export const SPONSOR_URL = "https://github.com/sponsors/CoolProp";
```

- [ ] **Step 2: Type-check**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx tsc --noEmit && cd ../..
```
Expected: no errors.

- [ ] **Step 3: Commit**

```bash
git add wrappers/GUI/src/constants.ts
git commit -m "feat(GUI): add shared SPONSOR_URL constant"
```

---

## Task 8: GUI — SponsorSplash component (TDD)

**Files:**
- Create: `wrappers/GUI/src/components/SponsorSplash.test.tsx`
- Create: `wrappers/GUI/src/components/SponsorSplash.tsx`

The component shows a one-time splash, gated by a `localStorage` key that embeds
the **major version**. It accepts the version as a prop (defaulting to the
generated `APP_VERSION`) so it is deterministically testable. The major version
is the integer before the first `.` of the semver string.

- [ ] **Step 1: Write the failing tests**

Create `wrappers/GUI/src/components/SponsorSplash.test.tsx`:

```tsx
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import SponsorSplash from "./SponsorSplash";

const KEY_V7 = "coolprop.sponsorSplash.seen.major.7";
const KEY_V8 = "coolprop.sponsorSplash.seen.major.8";

beforeEach(() => {
  localStorage.clear();
});

describe("SponsorSplash", () => {
  it("renders when no seen-key is set for the current major version", () => {
    render(<SponsorSplash version="7.2.1" />);
    expect(screen.getByText(/Enjoying CoolProp/i)).toBeInTheDocument();
  });

  it("does not render when the current major-version key is set", () => {
    localStorage.setItem(KEY_V7, "1");
    render(<SponsorSplash version="7.2.1" />);
    expect(screen.queryByText(/Enjoying CoolProp/i)).not.toBeInTheDocument();
  });

  it("renders again when only a previous major-version key is set", () => {
    localStorage.setItem(KEY_V7, "1");
    render(<SponsorSplash version="8.0.0" />);
    expect(screen.getByText(/Enjoying CoolProp/i)).toBeInTheDocument();
    expect(localStorage.getItem(KEY_V8)).toBeNull();
  });

  it("dismissing via 'Maybe later' hides it and sets the current major key", async () => {
    const user = userEvent.setup();
    render(<SponsorSplash version="7.2.1" />);
    await user.click(screen.getByRole("button", { name: /Maybe later/i }));
    expect(screen.queryByText(/Enjoying CoolProp/i)).not.toBeInTheDocument();
    expect(localStorage.getItem(KEY_V7)).toBe("1");
  });

  it("the Sponsor link points at the Sponsors URL and sets the seen key", async () => {
    const user = userEvent.setup();
    render(<SponsorSplash version="7.2.1" />);
    const link = screen.getByRole("link", { name: /Sponsor on GitHub/i });
    expect(link).toHaveAttribute("href", "https://github.com/sponsors/CoolProp");
    await user.click(link);
    expect(localStorage.getItem(KEY_V7)).toBe("1");
  });
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx vitest run src/components/SponsorSplash.test.tsx; cd ../..
```
Expected: FAIL — `Failed to resolve import "./SponsorSplash"` (file not created yet).

- [ ] **Step 3: Implement the component**

Create `wrappers/GUI/src/components/SponsorSplash.tsx`:

```tsx
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx vitest run src/components/SponsorSplash.test.tsx; cd ../..
```
Expected: PASS — all 5 tests green.

- [ ] **Step 5: Commit**

```bash
git add wrappers/GUI/src/components/SponsorSplash.tsx wrappers/GUI/src/components/SponsorSplash.test.tsx
git commit -m "feat(GUI): once-per-major-version SponsorSplash component"
```

---

## Task 9: GUI — wire splash + header link into App.tsx

**Files:**
- Modify: `wrappers/GUI/src/App.tsx:1-7` (imports), `:51-58` (header-right), `:75-77` (mount)

- [ ] **Step 1: Add imports**

In `wrappers/GUI/src/App.tsx`, the import block (lines 1-7) ends with:

```tsx
import UpdateChecker from "./components/UpdateChecker";
```

Add after it:

```tsx
import SponsorSplash from "./components/SponsorSplash";
import { SPONSOR_URL } from "./constants";
```

- [ ] **Step 2: Add the header Sponsor link**

In the `header-right` div, the About button currently reads:

```tsx
          <button
            className="tab-btn about-btn"
            onClick={() => setAboutOpen(true)}
            title="About / third-party notices"
          >
            About
          </button>
```

Insert this `<a>` immediately **before** that `<button>`:

```tsx
          <a
            className="tab-btn sponsor-btn"
            href={SPONSOR_URL}
            target="_blank"
            rel="noreferrer"
            title="Support CoolProp on GitHub Sponsors"
          >
            💚 Sponsor
          </a>
```

- [ ] **Step 3: Mount the splash**

Near the end of the component, the JSX currently reads:

```tsx
      {aboutOpen && <AboutModal onClose={() => setAboutOpen(false)} />}
      <UpdateChecker />
    </div>
```

Change it to:

```tsx
      {aboutOpen && <AboutModal onClose={() => setAboutOpen(false)} />}
      <UpdateChecker />
      <SponsorSplash />
    </div>
```

- [ ] **Step 4: Type-check and run the full GUI test suite**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx tsc --noEmit && npx vitest run; cd ../..
```
Expected: tsc clean; all tests pass (existing PropertyCalculator/hook tests + the 5 SponsorSplash tests).

- [ ] **Step 5: Commit**

```bash
git add wrappers/GUI/src/App.tsx
git commit -m "feat(GUI): add header Sponsor link and mount SponsorSplash"
```

---

## Task 10: GUI — Support link in AboutModal

**Files:**
- Modify: `wrappers/GUI/src/components/AboutModal.tsx:1` (import), `:22-31` (body)

- [ ] **Step 1: Add the import**

`AboutModal.tsx` line 1 currently reads:

```tsx
import { COOLPROP_VERSION, COOLPROP_GIT_HASH, NOTICES } from "../generated/notices";
```

Add a second import line after it:

```tsx
import { SPONSOR_URL } from "../constants";
```

- [ ] **Step 2: Add the Support paragraph**

After the "Built with Tauri…MIT license." `<p>` (which ends at line 31 with
`</p>`), and before the `<details className="about-notices">` block, insert:

```tsx
          <p>
            Support CoolProp:{" "}
            <a href={SPONSOR_URL} target="_blank" rel="noreferrer">
              💚 Sponsor on GitHub
            </a>
          </p>
```

- [ ] **Step 3: Type-check**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx tsc --noEmit; cd ../..
```
Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add wrappers/GUI/src/components/AboutModal.tsx
git commit -m "feat(GUI): add Support CoolProp link to About modal"
```

---

## Task 11: GUI — splash & sponsor-button styles

**Files:**
- Modify: `wrappers/GUI/src/App.css` (append a new section; reference existing `.about-btn` at line 336)

- [ ] **Step 1: Append styles**

Append to the end of `wrappers/GUI/src/App.css`:

```css

/* ── Sponsor link + splash ──────────────────────────── */

.sponsor-btn {
  font-size: 12px;
  color: #2ea043;
  border: 1px solid #2ea043;
  background: #e8f6ec;
  text-decoration: none;
  display: inline-flex;
  align-items: center;
  gap: 4px;
}
.sponsor-btn:hover {
  background: #d6efdd;
}

.modal-card.sponsor-splash-card {
  width: 400px;
  text-align: center;
}
.sponsor-splash-heart {
  font-size: 34px;
  line-height: 1;
  margin-bottom: 8px;
}
.modal-card.sponsor-splash-card .modal-footer {
  justify-content: center;
}
a.primary.sponsor-splash-cta {
  text-decoration: none;
  display: inline-flex;
  align-items: center;
}
```

- [ ] **Step 2: Build the frontend to confirm CSS/TS compile together**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npm run build; cd ../..
```
Expected: `tsc && vite build` completes with no errors (a `dist/` is produced).

> Note: `npm run build` runs `prebuild` → `build:licenses`, which regenerates
> `src/generated/notices.ts` from `CMakeLists.txt`. This is expected.

- [ ] **Step 3: Commit**

```bash
git add wrappers/GUI/src/App.css
git commit -m "style(GUI): sponsor button and splash styling"
```

---

## Task 12: Final verification & PR prep

- [ ] **Step 1: Full GUI test + typecheck sweep**

Run (from `wrappers/GUI/`):
```bash
cd wrappers/GUI && npx tsc --noEmit && npx vitest run; cd ../..
```
Expected: tsc clean; all suites pass.

- [ ] **Step 2: Confirm no stray generated/build artifacts are staged**

Run:
```bash
git status --short
```
Expected: clean (no `Web/_build/`, no `wrappers/GUI/dist/`, no
`wrappers/GUI/src/generated/`). If any appear, confirm they are gitignored;
do not commit them.

- [ ] **Step 3: Pre-PR adversarial review (MANDATORY per CLAUDE.md)**

Invoke the `superpowers:code-reviewer` subagent against the diff between
`feat/sponsor-promotion` and `origin/master`. Address or justify each blocking
finding before proceeding.

- [ ] **Step 4: Push and open the PR**

```bash
git push -u origin feat/sponsor-promotion
gh pr create --title "feat: promote GitHub Sponsors across repo, docs, and GUI" \
  --body "Implements the approved sponsor-promotion spec: FUNDING.yml, docs nav link + dismissable banner + Support page, and GUI header link + About entry + once-per-major-version splash. Spec: docs/superpowers/specs/2026-05-31-sponsor-promotion-design.md"
```

---

## Self-Review Notes

- **Spec coverage:** FUNDING.yml (Task 1) ✓; web nav link (Task 3) ✓; web banner CSS/JS + registration (Tasks 4-6) ✓; web Support page + cross-links (Task 2) ✓; GUI header link (Task 9) ✓; GUI About entry (Task 10) ✓; GUI once-per-major-version splash (Task 8) ✓; shared URL constant (Task 7, used by 8/9/10) ✓; tests (Task 8) ✓.
- **Spec deviation (intentional):** the dedicated page is `Web/sponsor.rst` added to the master `contents.rst` toctree, rather than `Web/online/sponsor.rst` + `online/index.rst` toctree — because `online/index.rst` has no toctree directive. Same user-facing result; cleaner integration.
- **Type consistency:** `SPONSOR_URL` (constants.ts) imported identically in SponsorSplash/App/AboutModal; `storageKey` format `coolprop.sponsorSplash.seen.major.<N>` matches the test keys `KEY_V7`/`KEY_V8`; `version` prop defaults to `APP_VERSION` from the generated notices module (same module AboutModal already imports).
