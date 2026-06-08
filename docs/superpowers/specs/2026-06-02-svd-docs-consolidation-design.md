# SVD validation docs improvements — design

**Date:** 2026-06-02
**Issue:** CoolProp-ndf2
**Status:** approved (design). Two independently-shippable deliverables, each its own PR.

Two related improvements to the SVD documentation, shipped as **separate PRs**:

- **Part A — consolidation** (toctree structure). Independent of any open PR.
- **Part B — DT validation panels** (content). Builds on #3068 (merged) and the
  DT-indexed surface from #2984 (merged).

---

## Part A — consolidate the SVD docs under `SVDSBTL.rst`

### Goal

Group the four scattered SVD pages into one coherent section by making
`SVDSBTL.rst` a section landing page with a nested toctree, so the sidebar
shows a single **SVDSBTL** entry that expands to its sub-pages.

### Current state

`Web/coolprop/index.rst` has one flat `toctree`. The four SVD pages are split:
`SVDSBTL.rst` at line 14 (among the backends); `SVDComponents.ipynb`,
`SVDSBTLValidation.ipynb`, `SVDSBTLHeatExchangerDemo.ipynb` at lines 25–27
(among unrelated notebooks). `SVDSBTL.rst` links `SVDComponents` but not the
Validation or HX-demo pages.

### Design

Target tree:

```
coolprop/
├─ SVDSBTL                       ← landing page
│   ├─ SVDComponents             (theory / internals)
│   ├─ SVDSBTLValidation         (accuracy)
│   └─ SVDSBTLHeatExchangerDemo  (application)
```

1. `Web/coolprop/index.rst`: remove the three notebook entries (lines 25–27);
   leave `SVDSBTL.rst` at line 14.
2. `Web/coolprop/SVDSBTL.rst`: add a nested toctree at the end:
   ```rst
   .. toctree::
       :maxdepth: 1
       :caption: SVDSBTL in depth

       SVDComponents
       SVDSBTLValidation
       SVDSBTLHeatExchangerDemo
   ```
   Order: theory → accuracy → application.
3. Add a short intro near the top of `SVDSBTL.rst` linking all three sub-pages
   (it currently surfaces only `SVDComponents`).

### Link integrity

Files do not move — only the toctree nesting changes — so every existing
`:doc:`/`:ref:` (IF97.rst, Tabular.rst, LowLevelAPI.rst, HighLevelAPI.rst,
changelog.rst) stays valid. Verification greps all cross-references.

---

## Part B — add DT (ρ, T) validation panels

### Goal

Add a third validation panel per fluid to `SVDSBTLValidation.ipynb` exercising
the DT-indexed surface (`DmassT_INPUTS`, landed in #2984).

### Key design point

The page currently validates **density error** `|Δρ/ρ|` for `(P,T)` and `(H,P)`
inputs. For DT, **density is the input** (exact by construction), so the DT
panel validates **pressure**: `|Δp/p|` of `SVDSBTL&HEOS` vs `HEOS`, both via
`DmassT_INPUTS` (the surface returns P directly). This mirrors the DT PR's own
`fig1301`.

### Design

- `INPUT_PAIRS = ['PT', 'HP', 'DT']` → 24 panels (8 fluids × 3).
- New `_dt_axes` helper: density range from a robust interior sample of the
  (T, ρ) rectangle (same pattern as the fixed `_hp_axes`), ρ on a **log** axis.
- `compute_error_grid` DT branch: grid is T (linear, y) × ρ (log, x); per cell
  `ref.update(DmassT_INPUTS, ρ, T)` and `svd.update(...)`, error
  `= |p_svd − p_ref| / p_ref`. **No saturation mask** — the two-phase dome is a
  valid DT region (SVDSBTL routes it via dome lever-rule; HEOS gives p_sat(T)),
  so it shows real (near-zero) error rather than a NaN wedge.
- `plot_from_grid` DT branch: x = ρ (log), y = T/Tc (linear), colorbar
  `log10 |Δp/p|` (z-range −12…−2 as for the other panels).
- Broaden the intro framing ("ρ-error for (P,T)/(H,P); p-error for (ρ,T)") and
  add a DT bullet to the Notes.

### Cost

A third surface (the DmassT surface) builds lazily per fluid on first DT query;
amortised by the surface cache from #3069. More panels = heavier first build,
fast thereafter.

### Sequencing

Branches off master (has #2984 DT surface + #3068 reworked compute cell + #3069
cache). Separate PR from Part A.

---

## Out of scope

- No prose/figure changes to the existing four pages beyond the framing/intro
  edits above.
- No merging of notebooks into a single page.

## Verification (both parts)

- `make html` builds with no new Sphinx `toctree`/`:doc:`/`:ref:` warnings.
- Part A: SVDSBTL sidebar entry expands to the three sub-pages; they no longer
  appear at the top level of `coolprop/`; grep confirms no broken refs.
- Part B: the notebook executes end-to-end with 0 errors; DT panels render real
  `|Δp/p|` fields (≈1e-8 single-phase, near-zero in the dome) for CO2 and a
  spot-check fluid; generator byte-stable / idempotent; `nbformat` valid.
