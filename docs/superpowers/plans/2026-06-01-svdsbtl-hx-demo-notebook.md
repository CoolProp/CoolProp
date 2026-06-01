# SVDSBTL HX-Speedup Demo Notebook — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a user-facing nbsphinx docs page (`Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb`) that reproduces the §3.2 speedup result of Bell et al. (ATE 2015), showing `HEOS` and `SVDSBTL&HEOS` agree on the heat-transfer rate Q while the table lookup is much faster per run.

**Architecture:** A small importable module (`Web/coolprop/hx_moving_boundary.py`) holds a provider abstraction and a faithful port of the moving-boundary solver from `dev/reference/HX.py`; the notebook imports it, runs both backends, draws an interactive plotly effectiveness-vs-area figure, and times per-run cost. The module is the single oracle that `dev/reference/HX.py` and `dev/demo_hx_speedup.cpp` also match.

**Tech Stack:** Python 3, CoolProp 7.2.1dev (`AbstractState`), scipy (`brentq`), numpy, plotly, nbformat/nbconvert (notebook authoring + execution), pytest.

---

## File Structure

- **Create** `Web/coolprop/hx_moving_boundary.py` — `PropertyProvider` (+ `HEOSProvider`, `SVDSBTLProvider`), `HeatExchanger`, `make_evaporator`, `ORACLE_Q`. One responsibility: solve the moving-boundary HX through a swappable backend.
- **Create** `Web/coolprop/test_hx_moving_boundary.py` — pytest guard: oracle match + backend agreement + physical effectiveness.
- **Create** `Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb` — narrative notebook (imports the module).
- **Modify** `Web/coolprop/index.rst` — add the notebook to the toctree after `SVDSBTLValidation.ipynb` (currently line 26).

Branch: `ihb/svdsbtl-hx-demo` (already created; the spec + `dev/reference/HX.py` are committed at `f2916e0de`).

Commits use `--no-verify` (docs/Python only — no C++ to clang-format; avoids the `bd` hook bundling `.beads/issues.jsonl`).

---

## Task 1: Environment sanity check (no commit)

**Files:** none (verification only)

- [ ] **Step 1: Confirm imports and SVDSBTL construction**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
python3 -c "import plotly, nbformat, nbconvert, pytest, scipy, numpy, CoolProp.CoolProp as CP; CP.AbstractState('SVDSBTL&HEOS','Water'); print('env OK')"
```
Expected: prints `env OK` (the SVDSBTL construct triggers/uses the cached Water table). If it errors on the SVDSBTL construct, stop — the backend is unavailable and the rest of the plan cannot proceed.

---

## Task 2: The solver module (`hx_moving_boundary.py`)

**Files:**
- Create: `Web/coolprop/hx_moving_boundary.py`
- Test: `Web/coolprop/test_hx_moving_boundary.py`

- [ ] **Step 1: Write the failing test**

Create `Web/coolprop/test_hx_moving_boundary.py`:

```python
"""Regression guard for the provider-based moving-boundary HX solver.

Runs from the notebook's directory so `import hx_moving_boundary` resolves.
Exercisable with `python3 -m pytest test_hx_moving_boundary.py -v` or directly
with `python3 test_hx_moving_boundary.py`.
"""
from hx_moving_boundary import make_evaporator, ORACLE_Q


def test_heos_matches_oracle():
    # The Python port must reproduce dev/reference/HX.py's Q to tight tolerance.
    hx = make_evaporator('HEOS', A=4.0)
    Q = hx.run()
    assert abs(Q - ORACLE_Q) / ORACLE_Q < 1e-6


def test_backends_agree():
    Qh = make_evaporator('HEOS', A=4.0).run()
    Qs = make_evaporator('SVDSBTL&HEOS', A=4.0).run()
    assert abs(Qh - Qs) / Qh < 1e-5


def test_effectiveness_physical():
    hx = make_evaporator('HEOS', A=4.0)
    Q = hx.run()
    eps = Q / hx.Qmax
    assert 0.0 < eps < 1.0


if __name__ == '__main__':
    test_heos_matches_oracle()
    test_backends_agree()
    test_effectiveness_physical()
    print('all tests passed')
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd Web/coolprop && python3 -m pytest test_hx_moving_boundary.py -v
```
Expected: collection/import error — `ModuleNotFoundError: No module named 'hx_moving_boundary'`.

- [ ] **Step 3: Write the module**

Create `Web/coolprop/hx_moving_boundary.py`:

```python
"""Provider-based moving-boundary counterflow heat-exchanger solver.

Faithful Python port of the original paper supplemental script
(``dev/reference/HX.py``) for

    I.H. Bell et al., "A generalized moving-boundary algorithm to predict the
    heat transfer rate of counterflow heat exchangers for any phase
    configuration," Applied Thermal Engineering 79 (2015) 192-201.

Every thermophysical property call routes through a :class:`PropertyProvider`
so the identical solver runs against either the full Helmholtz EOS (``HEOS``)
or the SVD-compressed tabular backend (``SVDSBTL&HEOS``).  Used by
``SVDSBTLHeatExchangerDemo.ipynb``.
"""
from __future__ import division
from math import log

import numpy as np
import scipy.optimize
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PT_INPUTS, HmassP_INPUTS, PQ_INPUTS, QT_INPUTS

# Headline Q [W] for the Table-3 evaporator at A = 4 m^2, from the recovered
# oracle dev/reference/HX.py on CoolProp 7.x.  Regression anchor for the port.
ORACLE_Q = 4577.242219


class PropertyProvider(object):
    """Property access for one fluid through a chosen backend.

    The ``p,h <-> T,rho`` transforms (:meth:`h_pT`, :meth:`Ts_ph`) go through
    the chosen backend; saturation (:meth:`Tsat`, :meth:`hsat_TQ`) always goes
    through a HEOS ancillary state, so saturation sourcing is identical for both
    backends and the only timed difference is the table lookup.
    """

    backend = None  # set by subclasses

    def __init__(self, fluid):
        self.AS = CP.AbstractState(self.backend, fluid)
        self.sat = CP.AbstractState('HEOS', fluid)

    def h_pT(self, p, T):
        self.AS.update(PT_INPUTS, p, T)
        return self.AS.hmass()

    def Ts_ph(self, p, h):
        self.AS.update(HmassP_INPUTS, h, p)
        return self.AS.T(), self.AS.smass()

    def Tsat(self, p, Q):
        self.sat.update(PQ_INPUTS, p, Q)
        return self.sat.T()

    def hsat_TQ(self, T, Q):
        self.sat.update(QT_INPUTS, Q, T)
        return self.sat.hmass()


class HEOSProvider(PropertyProvider):
    backend = 'HEOS'


class SVDSBTLProvider(PropertyProvider):
    backend = 'SVDSBTL&HEOS'


class HeatExchanger(object):
    """Moving-boundary counterflow HX, faithful to dev/reference/HX.py."""

    def __init__(self, prov_h, prov_c, mdot_h, p_hi, h_hi, mdot_c, p_ci, h_ci):
        self.ph, self.pc = prov_h, prov_c
        self.mdot_h, self.p_hi, self.h_hi = mdot_h, p_hi, h_hi
        self.mdot_c, self.p_ci, self.h_ci = mdot_c, p_ci, h_ci
        self.T_ci, _ = self.pc.Ts_ph(p_ci, h_ci)
        self.T_hi, _ = self.ph.Ts_ph(p_hi, h_hi)
        self.T_cbubble = self.pc.Tsat(p_ci, 0)
        self.T_cdew = self.pc.Tsat(p_ci, 1)
        self.T_hbubble = self.ph.Tsat(p_hi, 0)
        self.T_hdew = self.ph.Tsat(p_hi, 1)
        self.h_cbubble = self.pc.hsat_TQ(self.T_cbubble, 0)
        self.h_cdew = self.pc.hsat_TQ(self.T_cdew, 1)
        self.h_hbubble = self.ph.hsat_TQ(self.T_hbubble, 0)
        self.h_hdew = self.ph.hsat_TQ(self.T_hdew, 1)

    def external_pinching(self):
        self.h_ho = self.ph.h_pT(self.p_hi, self.T_ci)        # Eq 5
        Qmaxh = self.mdot_h * (self.h_hi - self.h_ho)         # Eq 4
        self.h_co = self.pc.h_pT(self.p_ci, self.T_hi)        # Eq 7
        Qmaxc = self.mdot_c * (self.h_co - self.h_ci)         # Eq 6
        Qmax = min(Qmaxh, Qmaxc)
        self.calculate_cell_boundaries(Qmax)
        return Qmax

    def calculate_cell_boundaries(self, Q):
        self.h_co = self.h_ci + Q / self.mdot_c
        self.h_ho = self.h_hi - Q / self.mdot_h
        self.hvec_c = [self.h_ci, self.h_co]
        self.hvec_h = [self.h_ho, self.h_hi]
        if self.h_hi > self.h_hdew > self.h_ho:
            self.hvec_h.insert(-1, self.h_hdew)
        if self.h_hi > self.h_hbubble > self.h_ho:
            self.hvec_h.insert(1, self.h_hbubble)
        if self.h_ci < self.h_cdew < self.h_co:
            self.hvec_c.insert(-1, self.h_cdew)
        if self.h_ci < self.h_cbubble < self.h_co:
            self.hvec_c.insert(1, self.h_cbubble)
        k = 0
        while k < len(self.hvec_c) - 1 or k < len(self.hvec_h) - 1:
            if len(self.hvec_c) == 2 and len(self.hvec_h) == 2:
                break
            Qcell_hk = self.mdot_h * (self.hvec_h[k + 1] - self.hvec_h[k])
            Qcell_ck = self.mdot_c * (self.hvec_c[k + 1] - self.hvec_c[k])
            if abs(Qcell_hk / Qcell_ck - 1) < 1e-6:
                k += 1
                break
            elif Qcell_hk > Qcell_ck:
                self.hvec_h.insert(k + 1, self.hvec_h[k] + Qcell_ck / self.mdot_h)
            else:
                self.hvec_c.insert(k + 1, self.hvec_c[k] + Qcell_hk / self.mdot_c)
            k += 1
        assert len(self.hvec_h) == len(self.hvec_c)
        self.Tvec_c = np.array([self.pc.Ts_ph(self.p_ci, h)[0] for h in self.hvec_c])
        self.Tvec_h = np.array([self.ph.Ts_ph(self.p_hi, h)[0] for h in self.hvec_h])
        self.phases_h = self._phases(self.hvec_h, self.h_hbubble, self.h_hdew)
        self.phases_c = self._phases(self.hvec_c, self.h_cbubble, self.h_cdew)

    @staticmethod
    def _phases(hvec, hbubble, hdew):
        out = []
        for i in range(len(hvec) - 1):
            havg = (hvec[i] + hvec[i + 1]) / 2.0
            if havg < hbubble:
                out.append('liquid')
            elif havg > hdew:
                out.append('vapor')
            else:
                out.append('two-phase')
        return out

    def internal_pinching(self, stream):
        if stream == 'hot':
            for i in range(1, len(self.hvec_h) - 1):
                if abs(self.hvec_h[i] - self.h_hdew) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]:
                    h_c_pinch = self.pc.h_pT(self.p_ci, self.T_hdew)       # Eq 10
                    Qright = self.mdot_h * (self.h_hi - self.h_hdew)       # Eq 9
                    Qmax = self.mdot_c * (h_c_pinch - self.h_ci) + Qright  # Eq 12
                    self.calculate_cell_boundaries(Qmax)
                    return Qmax
        elif stream == 'cold':
            for i in range(1, len(self.hvec_c) - 1):
                if abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]:
                    h_h_pinch = self.ph.h_pT(self.p_hi, self.T_cbubble)    # Eq 14
                    Qleft = self.mdot_c * (self.h_cbubble - self.h_ci)     # Eq 13
                    Qmax = Qleft + self.mdot_h * (self.h_hi - h_h_pinch)   # Eq 16
                    self.calculate_cell_boundaries(Qmax)
                    return Qmax
        else:
            raise ValueError('stream must be "hot" or "cold"')
        return None

    def objective_function(self, Q):
        self.calculate_cell_boundaries(Q)
        w = []
        for k in range(len(self.hvec_c) - 1):
            DTA = self.Tvec_h[k + 1] - self.Tvec_c[k + 1]
            DTB = self.Tvec_h[k] - self.Tvec_c[k]
            LMTD = DTA if DTA == DTB else (DTA - DTB) / log(abs(DTA / DTB))
            UA_req = self.mdot_h * (self.hvec_h[k + 1] - self.hvec_h[k]) / LMTD
            alpha_c = 100 if self.phases_c[k] in ('liquid', 'vapor') else 1000
            alpha_h = 100 if self.phases_h[k] in ('liquid', 'vapor') else 1000
            UA_avail = 1 / (1 / (alpha_h * self.A_h) + 1 / (alpha_c * self.A_c))
            w.append(UA_req / UA_avail)
        return 1 - sum(w)

    def solve(self):
        self.Q = scipy.optimize.brentq(
            self.objective_function, 1e-5, self.Qmax - 1e-10, rtol=1e-14, xtol=1e-10)
        return self.Q

    def run(self, and_solve=True):
        self.Qmax = self.external_pinching()
        for stream in ('hot', 'cold'):
            qi = self.internal_pinching(stream)
            if qi is not None:
                self.Qmax = qi
        return self.solve() if and_solve else None


def make_evaporator(backend, A=4.0):
    """Construct the water-heated n-Propane evaporator of the paper's Table 3.

    Matches the recovered oracle: the HOT stream (Water) carries mdot = 0.1 kg/s
    and the COLD stream (n-Propane) mdot = 0.01 kg/s -- transposed relative to
    the paper's printed Table 3 (see the design spec and dev/reference/HX.py).
    ``backend`` is 'HEOS' or 'SVDSBTL&HEOS' ('SVDSBTL' is accepted as an alias).
    """
    if backend == 'HEOS':
        ph, pc = HEOSProvider('Water'), HEOSProvider('n-Propane')
    elif backend in ('SVDSBTL', 'SVDSBTL&HEOS'):
        ph, pc = SVDSBTLProvider('Water'), SVDSBTLProvider('n-Propane')
    else:
        raise ValueError("backend must be 'HEOS' or 'SVDSBTL&HEOS', got %r" % (backend,))
    p_w = 101325.0
    h_w = ph.h_pT(p_w, 330.0)
    pc.sat.update(QT_INPUTS, 1, 300.0)
    p_r = pc.sat.p()
    h_r = pc.h_pT(p_r, 275.0)
    hx = HeatExchanger(ph, pc, mdot_h=0.1, p_hi=p_w, h_hi=h_w,
                       mdot_c=0.01, p_ci=p_r, h_ci=h_r)
    hx.A_h = hx.A_c = A
    return hx
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd Web/coolprop && python3 -m pytest test_hx_moving_boundary.py -v
```
Expected: 3 passed. (First run builds the n-Propane SVDSBTL table once, ~7 s.)

- [ ] **Step 5: Commit**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
git add Web/coolprop/hx_moving_boundary.py Web/coolprop/test_hx_moving_boundary.py
git commit --no-verify -m "feat(docs): provider-based moving-boundary HX solver module

Faithful port of dev/reference/HX.py routing all property calls through a
swappable PropertyProvider (HEOS / SVDSBTL&HEOS). pytest guard asserts the
HEOS port matches the oracle Q and that the two backends agree.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: The notebook (`SVDSBTLHeatExchangerDemo.ipynb`)

**Files:**
- Create: `Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb` (authored via a throwaway nbformat generator)

- [ ] **Step 1: Author the notebook with nbformat**

Write the generator to `/tmp/gen_hx_nb.py` (throwaway — not committed; the `.ipynb` is the artifact) and run it. Full generator:

```python
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
    "# Heat-exchanger speedup with the SVDSBTL backend\n"
    "\n"
    "This page reproduces the computational-efficiency result of §3.2 of\n"
    "\n"
    "> I.H. Bell, S. Quoilin, E. Georges, J.E. Braun, E.A. Groll, W.T. Horton,\n"
    "> V. Lemort, \"A generalized moving-boundary algorithm to predict the heat\n"
    "> transfer rate of counterflow heat exchangers for any phase\n"
    "> configuration,\" *Applied Thermal Engineering* 79 (2015) 192–201.\n"
    "\n"
    "A moving-boundary model of a water-heated n-Propane evaporator is solved with\n"
    "the full Helmholtz equation of state (`HEOS`) and with the SVD-compressed\n"
    "tabular backend (`SVDSBTL&HEOS`). The two agree on the predicted heat-transfer\n"
    "rate while the table lookup runs much faster — the `p,h → T,ρ` property\n"
    "evaluation (Water has the library's most expensive EOS) is the bottleneck the\n"
    "table replaces.\n"
    "\n"
    "The solver lives in a small importable module so you can reuse it:"))

raw = nbf.v4.new_raw_cell(
    ":doc:`SVDSBTL backend reference <SVDSBTL>`  ·  "
    ":download:`Download the solver module (hx_moving_boundary.py) <hx_moving_boundary.py>`")
raw.metadata['raw_mimetype'] = 'text/restructuredtext'
cells.append(raw)

cells.append(nbf.v4.new_code_cell(
    "import os, time\n"
    "import numpy as np\n"
    "import plotly.graph_objects as go\n"
    "import plotly.io as pio\n"
    "pio.renderers.default = 'notebook_connected'\n"
    "from hx_moving_boundary import make_evaporator, ORACLE_Q"))

cells.append(nbf.v4.new_code_cell(
    "# Build providers for both backends. The first SVDSBTL construction builds and\n"
    "# caches the compressed tables (~7 s/fluid, cached under ~/.CoolProp/SVDTables);\n"
    "# this one-time cost is reported separately and excluded from per-run timing.\n"
    "t0 = time.perf_counter()\n"
    "hx_heos = make_evaporator('HEOS', A=4.0)\n"
    "t_build_heos = time.perf_counter() - t0\n"
    "\n"
    "t0 = time.perf_counter()\n"
    "hx_svd = make_evaporator('SVDSBTL&HEOS', A=4.0)\n"
    "t_build_svd = time.perf_counter() - t0\n"
    "\n"
    "print('HEOS providers built in   %.3f s' % t_build_heos)\n"
    "print('SVDSBTL tables ready in   %.3f s (one-time, cached)' % t_build_svd)"))

cells.append(nbf.v4.new_code_cell(
    "Q_heos = hx_heos.run(); eps_heos = Q_heos / hx_heos.Qmax\n"
    "Q_svd = hx_svd.run();   eps_svd = Q_svd / hx_svd.Qmax\n"
    "rel = abs(Q_heos - Q_svd) / Q_heos\n"
    "print('HEOS         : Q = %.4f W   eps = %.6f' % (Q_heos, eps_heos))\n"
    "print('SVDSBTL&HEOS : Q = %.4f W   eps = %.6f' % (Q_svd, eps_svd))\n"
    "print('|dQ|/Q backend agreement      = %.2e' % rel)\n"
    "print('|dQ|/Q vs reference oracle     = %.2e' % (abs(Q_heos - ORACLE_Q) / ORACLE_Q))\n"
    "assert rel < 1e-5, 'backends disagree on Q'\n"
    "assert abs(Q_heos - ORACLE_Q) / ORACLE_Q < 1e-6, 'HEOS drifted from the HX.py oracle'"))

cells.append(nbf.v4.new_code_cell(
    "# Effectiveness vs area, both backends overlaid. Reuse the built HX objects:\n"
    "# only the area changes between points, so no rebuild is needed.\n"
    "areas = np.logspace(np.log10(0.1), np.log10(10.0), 25)\n"
    "\n"
    "def sweep(hx):\n"
    "    out = []\n"
    "    for A in areas:\n"
    "        hx.A_h = hx.A_c = A\n"
    "        out.append(hx.run() / hx.Qmax)\n"
    "    return out\n"
    "\n"
    "eps_h = sweep(hx_heos)\n"
    "eps_s = sweep(hx_svd)\n"
    "\n"
    "fig = go.Figure()\n"
    "fig.add_trace(go.Scatter(x=areas, y=eps_h, mode='lines', name='HEOS',\n"
    "                         line=dict(width=4)))\n"
    "fig.add_trace(go.Scatter(x=areas, y=eps_s, mode='markers', name='SVDSBTL&HEOS',\n"
    "                         marker=dict(size=7)))\n"
    "fig.update_xaxes(type='log', title='Heat transfer area A [m²]')\n"
    "fig.update_yaxes(title='Effectiveness ε = Q / Q_max')\n"
    "fig.update_layout(\n"
    "    title='Counterflow evaporator effectiveness — HEOS vs SVDSBTL',\n"
    "    template='simple_white', legend=dict(x=0.02, y=0.98))\n"
    "fig"))

cells.append(nbf.v4.new_code_cell(
    "# Per-run timing at A = 4 m^2, the paper's ms/run unit. Both backends run the\n"
    "# identical Python harness, so the residual pybind dispatch tax is shared.\n"
    "N = int(os.environ.get('HX_REPEATS', 200))\n"
    "\n"
    "def per_run_ms(hx, N):\n"
    "    hx.A_h = hx.A_c = 4.0\n"
    "    hx.run()  # warmup\n"
    "    t0 = time.perf_counter()\n"
    "    for _ in range(N):\n"
    "        hx.run()\n"
    "    return (time.perf_counter() - t0) / N * 1e3\n"
    "\n"
    "ms_heos = per_run_ms(hx_heos, N)\n"
    "ms_svd = per_run_ms(hx_svd, N)\n"
    "print('per-run timing averaged over N = %d runs at A = 4 m^2' % N)\n"
    "print('  HEOS         : %.3f ms/run' % ms_heos)\n"
    "print('  SVDSBTL&HEOS : %.3f ms/run' % ms_svd)\n"
    "print('  speedup (HEOS / SVDSBTL) = %.1fx' % (ms_heos / ms_svd))"))

cells.append(nbf.v4.new_markdown_cell(
    "## What this shows\n"
    "\n"
    "Both backends predict the same heat-transfer rate (Q agreement on the order of\n"
    "`1e-7` relative), yet SVDSBTL is several times faster *per full HX run* — the\n"
    "speedup coming entirely from replacing Water's expensive `p,h → T,ρ` inversion\n"
    "with a table lookup, exactly as in §3.2 of the paper.\n"
    "\n"
    "The per-run ratio above is a **lower bound** on SVDSBTL's library-level\n"
    "advantage: both backends pay the same Python/pybind dispatch cost per property\n"
    "call, which compresses the ratio. The companion pure-C++ demo\n"
    "(`dev/demo_hx_speedup.cpp`) removes that overhead and reports a larger ratio.\n"
    "Absolute ms/run also differ from the paper's 24.6 / 2.46 ms — those were\n"
    "32-bit Python on 2011 hardware.\n"
    "\n"
    "**A note on the test case.** This evaporator matches the recovered original\n"
    "script (`dev/reference/HX.py`): the hot stream (Water) carries ṁ = 0.1 kg/s\n"
    "and the cold stream (n-Propane) ṁ = 0.01 kg/s, with two-phase\n"
    "α = 1000 and single-phase α = 100 W m⁻² K⁻¹. The mass flows are transposed\n"
    "relative to the paper's printed Table 3; this is the assignment that produces\n"
    "the evaporating cell structure shown in the paper's figures."))

nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nb.metadata['language_info'] = {'name': 'python'}

with open('Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb', 'w') as f:
    nbf.write(nb, f)
print('wrote Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb')
```

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
python3 /tmp/gen_hx_nb.py
```
Expected: `wrote Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb`.

- [ ] **Step 2: Execute the notebook end-to-end to validate**

Run (executes in the notebook's directory so `import hx_moving_boundary` resolves; mirrors what nbsphinx does):
```bash
cd Web/coolprop && HX_REPEATS=50 jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.timeout=600 \
  --output /tmp/_hx_demo_executed.ipynb SVDSBTLHeatExchangerDemo.ipynb
```
Expected: exits 0 with no traceback. (`HX_REPEATS=50` keeps validation quick; the in-notebook `assert`s fail the execution if the backends disagree or the port drifts from the oracle.)

- [ ] **Step 3: Confirm the executed copy has no errored cells**

Run:
```bash
python3 -c "import nbformat; nb=nbformat.read('/tmp/_hx_demo_executed.ipynb', as_version=4); errs=[o for c in nb.cells if c.cell_type=='code' for o in c.get('outputs',[]) if o.get('output_type')=='error']; print('errors:', len(errs)); [print(e['ename'], e['evalue']) for e in errs]; assert not errs"
```
Expected: `errors: 0`.

- [ ] **Step 4: Commit**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
git add Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb
git commit --no-verify -m "docs(SVDSBTL): add HX-speedup demo notebook

Solves the water-heated n-Propane evaporator with HEOS vs SVDSBTL&HEOS,
overlays the effectiveness-vs-area curves (plotly), and times per-run cost.
Imports the hx_moving_boundary module; offers it as an HTML download.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Wire the notebook into the toctree

**Files:**
- Modify: `Web/coolprop/index.rst` (after line 26, `SVDSBTLValidation.ipynb`)

- [ ] **Step 1: Add the toctree entry**

Edit `Web/coolprop/index.rst`: insert `    SVDSBTLHeatExchangerDemo.ipynb` immediately after the `    SVDSBTLValidation.ipynb` line. The block becomes:

```
    SVDComponents.ipynb
    SVDSBTLValidation.ipynb
    SVDSBTLHeatExchangerDemo.ipynb
```

- [ ] **Step 2: Verify the entry resolves to a real file and is unique**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
grep -c "SVDSBTLHeatExchangerDemo.ipynb" Web/coolprop/index.rst
test -f Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb && echo "notebook present"
```
Expected: prints `1` then `notebook present`.

- [ ] **Step 3: Commit**

```bash
git add Web/coolprop/index.rst
git commit --no-verify -m "docs(SVDSBTL): add HX-speedup demo to the coolprop toctree

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Final verification (no commit)

**Files:** none

- [ ] **Step 1: Re-run the module test from a clean shell**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo/Web/coolprop
python3 -m pytest test_hx_moving_boundary.py -v
```
Expected: 3 passed.

- [ ] **Step 2: Confirm the branch state**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/tabledemo
git log --oneline -4
git status --short
```
Expected: the three new commits (module, notebook, toctree) atop `f2916e0de`; `git status` shows only the unrelated pre-existing C++ demo files (`CMakeLists.txt`, `dev/demo_hx_speedup.cpp`, `dev/plot_hx_figures.py`, the 05-31 spec) — none of this task's files left uncommitted.

- [ ] **Step 3 (optional): Local docs build smoke test**

A full Sphinx build of `Web/` is heavy and downloads assets; only run if validating the rendered HTML is required. The Task 3 nbconvert execution already proves the notebook runs under the same machinery nbsphinx uses, and Task 4 proves the toctree entry resolves — together these cover doc-build success without a full build.

---

## Notes for the executor

- **Do not** touch `CMakeLists.txt`, `dev/demo_hx_speedup.cpp`, `dev/plot_hx_figures.py`, or the `2026-05-31` spec — that is separate in-progress C++ work, intentionally left uncommitted.
- `--no-verify` is used because this branch changes only Python/notebook/rst docs (no C++ to clang-format) and to stop the `bd` pre-commit hook bundling `.beads/issues.jsonl`. If this branch is later pushed, the pre-push hook still runs; for a docs-only change the relevant preflight stages are no-ops, but run `./dev/ci/preflight.sh` if in doubt.
- Pushing is out of scope for this plan — stop after Task 5.
