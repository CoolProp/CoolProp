"""Generate the SVDSBTL validation notebook (8 fluids x {PT, HP} inputs, rho error).

Writes ``Web/coolprop/SVDSBTLValidation.ipynb`` as a *stripped* notebook
(no embedded outputs).  The docs build re-executes it in place on every
``make html``: ``Web/conf.py`` walks every ``.ipynb`` and runs
``jupyter nbconvert --execute`` (timeout bumped to 1 h, because SVDSBTL
surfaces lazily build the first time each (fluid, input_pair) is queried
in a given environment).  The committed copy therefore intentionally
carries no outputs — they are regenerated at build time.

This file is the single source of truth for the notebook; re-run it after
any edit and commit the regenerated ``.ipynb`` (the two must stay in sync).

Docs build environment must have:
  - CoolProp Python wrapper built from current source (``pip install .`` at repo root)
  - plotly (>= 5), joblib installed
  - Writable cache directory at ``~/.CoolProp/SVDTables-docs/`` (or whatever
    ``COOLPROP_ALTERNATIVE_SVDTABLES_DIRECTORY`` points to). Surfaces persist
    across builds; first build of a fresh (fluid, input_pair) takes ~20-80 s,
    subsequent builds reuse the serialized surface.

Local verification (renders the HTML the docs build would produce):
  jupyter nbconvert --execute --to html \\
    --ExecutePreprocessor.timeout=3600 \\
    Web/coolprop/SVDSBTLValidation.ipynb
"""

from pathlib import Path

import nbformat as nbf

FLUIDS = ['Water', 'Argon', 'Hydrogen', 'R1234yf', 'R245fa', 'D6', 'CO2', 'Helium']
INPUT_PAIRS = ['PT', 'HP']
N_X, N_Y = 300, 300        # 90k cells per (fluid, input_pair)

nb = nbf.v4.new_notebook()

# Deterministic cell IDs so re-running this generator produces a byte-stable
# .ipynb diff (nbformat's default is a random 8-hex id per call, which dirties
# every cell on every regen).
_cell_counter = [0]
def _md(text):
    _cell_counter[0] += 1
    return nbf.v4.new_markdown_cell(text, id=f'cell-{_cell_counter[0]:03d}')
def _code(src):
    _cell_counter[0] += 1
    return nbf.v4.new_code_cell(src, id=f'cell-{_cell_counter[0]:03d}')

nb.cells = [
    _md(
        "# SVDSBTL backend — validation against HEOS\n\n"
        "This page compares the **SVDSBTL** surrogate backend against the "
        "**HEOS** Helmholtz-energy reference for `rho` as a function of two "
        "different input pairs: `(P, T)` and `(H, P)`. Each fluid below shows "
        "`log10 | (rho_SVDSBTL - rho_HEOS) / rho_HEOS |` over single-phase "
        "rectangles in the two coordinate systems.\n\n"
        "The heatmaps are interactive — hover for the local error value, drag to pan, "
        "scroll to zoom, double-click to reset. The figures are static HTML (no kernel, "
        "no Binder, no JupyterLite); the interactivity is plotly.js embedded by nbsphinx.\n\n"
        "**Backends:** `HEOS` (reference) vs `SVDSBTL&HEOS` (rank-truncated SVD surrogate over HEOS).  \n"
        "**Inputs:** `PT_INPUTS` and `HmassP_INPUTS`. **Output:** `rhomass()` [kg/m^3].  \n"
        "**Fluids:** " + ", ".join(f"`{f}`" for f in FLUIDS) + ".\n"
    ),
    _code(
        "import os\n"
        "from pathlib import Path\n"
        "\n"
        "# Route SVDSBTL surface I/O to a docs-dedicated cache directory so the\n"
        "# build does not collide with developer interactive work in\n"
        "# ~/.CoolProp/SVDTables/. Setting this *before* importing CoolProp lets\n"
        "# the C++ config-init read the env var on first access, and joblib's\n"
        "# loky workers inherit the env so they hit the same cache.\n"
        "_cache = os.environ.setdefault(\n"
        "    'COOLPROP_ALTERNATIVE_SVDTABLES_DIRECTORY',\n"
        "    str(Path.home() / '.CoolProp' / 'SVDTables-docs'),\n"
        ")\n"
        "Path(_cache).mkdir(parents=True, exist_ok=True)\n"
        "\n"
        "import numpy as np\n"
        "import plotly.io as pio\n"
        "import plotly.graph_objects as go\n"
        "import CoolProp.CoolProp as CP\n"
        "from CoolProp.CoolProp import PT_INPUTS, HmassP_INPUTS, QT_INPUTS\n"
        "\n"
        "# Also set via the explicit config setter in the parent process. NOTE:\n"
        "# this only mutates the parent's Configuration singleton; joblib's loky\n"
        "# workers spawn fresh interpreters with their own singletons and rely\n"
        "# entirely on env-var inheritance (which is why os.environ.setdefault\n"
        "# above MUST stay before the CoolProp import — do not remove it).\n"
        "CP.set_config_string(CP.ALTERNATIVE_SVDTABLES_DIRECTORY, _cache)\n"
        "\n"
        "pio.renderers.default = 'notebook_connected'\n"
        "print(f'SVDSBTL cache directory: {_cache}')\n"
        "\n"
        f"FLUIDS = {FLUIDS!r}\n"
        f"INPUT_PAIRS = {INPUT_PAIRS!r}\n"
        f"N_X, N_Y = {N_X}, {N_Y}    # grid resolution per axis ({N_X*N_Y//1000}k cells per panel)\n"
        "BAND_REL = 0.02            # half-width band around saturation line to mask (PT only)\n"
    ),
    _code(
        "def _hp_axes(ref, Tmin, Tmax, pmin, pmax, N_Y):\n"
        "    \"\"\"Probe the four corners of the (T, p) rectangle to estimate an\n"
        "    enthalpy range that covers the single-phase region for the HP plot.\n"
        "    \"\"\"\n"
        "    hs = []\n"
        "    for T_q, p_q in [(Tmin, pmin), (Tmin, pmax), (Tmax, pmin), (Tmax, pmax)]:\n"
        "        try:\n"
        "            ref.update(PT_INPUTS, p_q, T_q)\n"
        "            hs.append(ref.hmass())\n"
        "        except Exception:\n"
        "            pass\n"
        "    if not hs:\n"
        "        return np.linspace(0.0, 1e6, N_Y)\n"
        "    h_min, h_max = min(hs), max(hs)\n"
        "    if h_max == h_min:\n"
        "        # All surviving probes returned the same enthalpy. Expand by a\n"
        "        # small relative delta so np.linspace produces a non-degenerate\n"
        "        # axis; the panel will be uninformative but at least renders.\n"
        "        delta = max(abs(h_min) * 1e-6, 1e-6)\n"
        "        h_min -= delta\n"
        "        h_max += delta\n"
        "    return np.linspace(h_min, h_max, N_Y)\n"
        "\n"
        "def compute_error_grid(fluid, input_pair, N_X=N_X, N_Y=N_Y):\n"
        "    \"\"\"Return (x_axis, y_axis, log10_err, Tc, pc, x_label, y_label, err_msg).\n"
        "\n"
        "    input_pair is 'PT' or 'HP'.\n"
        "      'PT' -> x = p (log, in Pa), y = T (linear, in K).\n"
        "             Mask a +/- BAND_REL band around the saturation line.\n"
        "      'HP' -> x = p (log, in Pa), y = h (linear, in J/kg).\n"
        "             No explicit saturation mask; the SVDSBTL surface is undefined\n"
        "             inside the dome and queries there are caught and set to NaN.\n"
        "\n"
        "    SVDSBTL surface build failures are caught and reported via err_msg;\n"
        "    the returned grid is all-NaN in that case.\n"
        "    \"\"\"\n"
        "    ref = CP.AbstractState('HEOS', fluid)\n"
        "    Tc = ref.T_critical()\n"
        "    pc = ref.p_critical()\n"
        "    Tmin = max(ref.Tmin(), 0.6 * Tc)\n"
        "    Tmax = min(ref.Tmax(), 1.4 * Tc)\n"
        "    pmin = max(1e3, 0.01 * pc)\n"
        "    pmax = 3.0 * pc\n"
        "\n"
        "    p_axis = np.geomspace(pmin, pmax, N_X)\n"
        "    if input_pair == 'PT':\n"
        "        y_axis = np.linspace(Tmin, Tmax, N_Y)\n"
        "        y_label = 'reduced temperature T / T_c'\n"
        "        # Saturation line p_sat(T) for masking.\n"
        "        psat = np.full_like(y_axis, np.nan)\n"
        "        for i, T in enumerate(y_axis):\n"
        "            if T < Tc:\n"
        "                try:\n"
        "                    ref.update(QT_INPUTS, 0.0, T)\n"
        "                    psat[i] = ref.p()\n"
        "                except Exception:\n"
        "                    pass\n"
        "    elif input_pair == 'HP':\n"
        "        y_axis = _hp_axes(ref, Tmin, Tmax, pmin, pmax, N_Y)\n"
        "        y_label = 'specific enthalpy h [kJ/kg]'\n"
        "        psat = None\n"
        "    else:\n"
        "        raise ValueError(f'unknown input_pair: {input_pair!r}')\n"
        "\n"
        "    x_label = 'reduced pressure p / p_c'\n"
        "    err = np.full((N_Y, N_X), np.nan)\n"
        "\n"
        "    try:\n"
        "        svd = CP.AbstractState('SVDSBTL&HEOS', fluid)\n"
        "    except Exception as e:\n"
        "        return (p_axis, y_axis, err, Tc, pc, x_label, y_label,\n"
        "                f'SVDSBTL surface build failed: {e}')\n"
        "\n"
        "    for i, y in enumerate(y_axis):\n"
        "        for j, p in enumerate(p_axis):\n"
        "            if psat is not None and np.isfinite(psat[i]) and abs(p - psat[i]) / psat[i] < BAND_REL:\n"
        "                continue\n"
        "            try:\n"
        "                if input_pair == 'PT':\n"
        "                    ref.update(PT_INPUTS, p, y)\n"
        "                    svd.update(PT_INPUTS, p, y)\n"
        "                else:  # HP\n"
        "                    ref.update(HmassP_INPUTS, y, p)\n"
        "                    svd.update(HmassP_INPUTS, y, p)\n"
        "                rho_ref = ref.rhomass()\n"
        "                rho_svd = svd.rhomass()\n"
        "                if rho_ref > 0:\n"
        "                    err[i, j] = abs(rho_svd - rho_ref) / rho_ref\n"
        "            except Exception:\n"
        "                pass\n"
        "    # Floor zero / near-zero errors at 1e-16 so exact-match cells render at\n"
        "    # the colorbar floor (best-possible accuracy) instead of going to -inf\n"
        "    # and being shown as NaN/white — which would falsely suggest a failed\n"
        "    # query. SVDSBTL inside the two-phase dome returns the same density as\n"
        "    # HEOS (both go through the same saturation lookup), so err is 0 there.\n"
        "    with np.errstate(divide='ignore', invalid='ignore'):\n"
        "        log_err = np.log10(np.maximum(err, 1e-16))\n"
        "    return p_axis, y_axis, log_err, Tc, pc, x_label, y_label, None\n"
    ),
    _code(
        "def plot_from_grid(fluid, input_pair, grid):\n"
        "    p_axis, y_axis, log_err, Tc, pc, x_label, y_label, err_msg = grid\n"
        "    pr = p_axis / pc\n"
        "    title = f'{fluid} [{input_pair}] \\u2014 Tc={Tc:.2f} K, pc={pc/1e5:.3f} bar'\n"
        "    if err_msg is not None:\n"
        "        title += '  [SVDSBTL surface unavailable]'\n"
        "    if input_pair == 'PT':\n"
        "        # x = p (log), y = T (linear). z is stored as [y_index, x_index].\n"
        "        Tr = y_axis / Tc\n"
        "        heat = go.Heatmap(\n"
        "            x=pr, y=Tr, z=log_err,\n"
        "            colorscale='Viridis', zmin=-12, zmax=-2,\n"
        "            colorbar=dict(title='log10 |\\u0394\\u03c1/\\u03c1|'),\n"
        "            hovertemplate=('p/p_c=%{x:.3f}<br>'\n"
        "                           'T/T_c=%{y:.3f}<br>'\n"
        "                           'log10 err=%{z:.2f}<extra></extra>'),\n"
        "        )\n"
        "        xaxis = dict(title=x_label, type='log')\n"
        "        yaxis = dict(title=y_label)\n"
        "    else:  # HP — Mollier convention: x = h (linear), y = p (log).\n"
        "        h_disp = y_axis / 1e3   # J/kg -> kJ/kg for display\n"
        "        heat = go.Heatmap(\n"
        "            x=h_disp, y=pr, z=log_err.T,\n"
        "            colorscale='Viridis', zmin=-12, zmax=-2,\n"
        "            colorbar=dict(title='log10 |\\u0394\\u03c1/\\u03c1|'),\n"
        "            hovertemplate=('h=%{x:.1f} kJ/kg<br>'\n"
        "                           'p/p_c=%{y:.3f}<br>'\n"
        "                           'log10 err=%{z:.2f}<extra></extra>'),\n"
        "        )\n"
        "        xaxis = dict(title=y_label)                     # 'specific enthalpy h [kJ/kg]'\n"
        "        yaxis = dict(title='reduced pressure p / p_c', type='log')\n"
        "    fig = go.Figure(heat)\n"
        "    fig.update_layout(\n"
        "        title=title,\n"
        "        xaxis=xaxis,\n"
        "        yaxis=yaxis,\n"
        "        width=720, height=480,\n"
        "        margin=dict(l=70, r=20, t=50, b=60),\n"
        "    )\n"
        "    if err_msg is not None:\n"
        "        print(err_msg)\n"
        "    return fig\n"
    ),
    _md(
        "### Parallel grid computation\n\n"
        "Each `(fluid, input_pair)` pair is an independent unit of work, so we compute "
        "all 12 of them concurrently via `joblib` "
        "(process-based workers; `n_jobs=min(len(FLUIDS)*len(INPUT_PAIRS), cpu_count)`). "
        "Wall-clock and per-panel timings are printed below.\n"
    ),
    _code(
        "import time, os\n"
        "from joblib import Parallel, delayed\n"
        "\n"
        "def compute_timed(fluid, input_pair):\n"
        "    t0 = time.perf_counter()\n"
        "    grid = compute_error_grid(fluid, input_pair)\n"
        "    return (fluid, input_pair), grid, time.perf_counter() - t0\n"
        "\n"
        "tasks = [(f, ip) for f in FLUIDS for ip in INPUT_PAIRS]\n"
        "n_jobs = min(len(tasks), os.cpu_count() or 1)\n"
        "t0 = time.perf_counter()\n"
        "results = Parallel(n_jobs=n_jobs, backend='loky')(\n"
        "    delayed(compute_timed)(f, ip) for f, ip in tasks\n"
        ")\n"
        "wall = time.perf_counter() - t0\n"
        "\n"
        "GRIDS = {key: grid for key, grid, _ in results}\n"
        "print(f'wall-clock: {wall:.2f} s   (n_jobs={n_jobs}, '\n"
        "      f'cores={os.cpu_count()}, panels={len(tasks)})')\n"
        "for key, _, dt in results:\n"
        "    print(f'  {key[0]:9s} {key[1]:>3s}  {dt:6.2f} s')\n"
        "print(f'  sum         {sum(dt for _, _, dt in results):6.2f} s   '\n"
        "      f'speedup vs serial: {sum(dt for _, _, dt in results) / wall:.2f}x')\n"
    ),
]

for fluid in FLUIDS:
    nb.cells.append(_md(f"## {fluid}"))
    for ip in INPUT_PAIRS:
        nb.cells.append(_md(f"### {ip}"))
        nb.cells.append(_code(
            f"plot_from_grid({fluid!r}, {ip!r}, GRIDS[({fluid!r}, {ip!r})])"
        ))

nb.cells.append(_md(
    "## Notes\n\n"
    "- **PT panels** mask a 2-percent band around the saturation line, since "
    "(T, p) inside that band straddles the dome at ULP scale and SVDSBTL is "
    "single-phase only.\n"
    "- **HP panels** have no explicit mask: the two-phase dome is a true 2D "
    "region in (h, p), and SVDSBTL surface queries inside it throw — caught "
    "and reported as NaN. The white wedge in each HP heatmap is the dome.\n"
    "- `R1234yf`, `R245fa`, and `D6` have narrow native `Tmin/Tmax` envelopes "
    "that may clip the plot domain below the nominal `0.6 Tc` lower bound.\n"
    "- The 12 panels are computed concurrently via `joblib` (process workers; "
    "`n_jobs=min(panels, cpu_count)`). Wall-clock and per-panel breakdown are "
    "in the parallel-compute cell above.\n"
    "- This notebook is re-executed on every `make html` build "
    "(`Web/conf.py` pre-executes all `.ipynb` via `jupyter nbconvert --execute`). "
    "Surfaces are cached at the directory printed in the first code cell "
    "(default `~/.CoolProp/SVDTables-docs/`, override via "
    "`COOLPROP_ALTERNATIVE_SVDTABLES_DIRECTORY`). First build of a fresh "
    "(fluid, input_pair) takes 20-80 s; subsequent builds reuse the serialized "
    "surface and run in ~1 s.\n"
    "- Source generator: `Web/coolprop/_gen/gen_SVDSBTLValidation.py`.\n"
))

out = str(Path(__file__).resolve().parents[1] / 'SVDSBTLValidation.ipynb')
with open(out, 'w') as f:
    nbf.write(nb, f)

n_code = sum(1 for c in nb.cells if c.cell_type == 'code')
n_md = sum(1 for c in nb.cells if c.cell_type == 'markdown')
print(f'wrote {out}  ({n_code} code cells, {n_md} markdown cells, stripped)')
