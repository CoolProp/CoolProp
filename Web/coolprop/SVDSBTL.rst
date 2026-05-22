.. _SVDSBTL:

*************************************************
SVDSBTL --- SVD-Compressed Tabular Lookup Backend
*************************************************

SVDSBTL is a tabular-lookup backend that combines an atlas of
analytically-bounded regions with a low-rank SVD of each region's
property surfaces.  It produces **sub-microsecond** property
evaluations (~200-400 ns per probe in the batched ``fast_evaluate``
path) at IAPWS G13-15 conformance for water, and a single-digit
percent accuracy ceiling for the multi-fluid HEOS-backed presets,
at a disk footprint of ~10 MB per (fluid, input pair, source backend).

This page is the **user-facing** page: how to use the backend, what it
guarantees, and when to pick it.  For the underlying SVD math
(rank-:math:`r` reconstruction, cubic Hermite slopes, axis transforms),
see :doc:`SVDComponents`.

.. _SVDSBTL-architecture:

What it is
==========

A traditional bicubic table covers the full :math:`(p, h)` or
:math:`(p, T)` envelope with one coarse grid and stores the source
backend's property values at each cell.  Two known limitations of
that approach: (a) the **saturation boundary** can't be expressed on
a regular grid without padding cells on either side, so a single
table either crosses the dome (with the discontinuity smearing into
the bicubic interpolant) or excludes a stripe around it (so the
backend can't answer there); and (b) **memory grows quadratically**
with grid resolution — the bicubic table for a single fluid at
useful accuracy is hundreds of MB.  SVDSBTL replaces the single
bicubic table with three layered components that address both:

1. **Region atlas.**  The thermodynamic envelope is partitioned into
   ~3-7 disjoint regions in :math:`(p, h)` or :math:`(p, T)` — LIQUID
   (subcooled), VAPOR (superheated), SUPER (supercritical), with
   further splits for IF97 (R1/R3, R2/R3 boundary curves).  Each
   region's secondary axis is normalised to :math:`[0, 1]` via two
   analytic boundary curves — saturation dome on the subcritical
   side, low-T / high-T isotherms on the supercritical side.  Lookups
   first AABB-test against region bounding boxes (cheap), then
   evaluate the precise boundary curves on the candidate.

2. **Per-region SVD.**  Inside each normalised region, every tabulated
   property (:math:`\rho, h, s, u, w, \eta, \lambda`) is represented
   as a rank-:math:`r` SVD of a dense :math:`N_p \times N_h` (default
   :math:`200 \times 800`, :math:`r = 20`) sample grid.  Hermite
   bicubic interpolation between SVD grid points gives :math:`C^1`
   continuity across the region.  Density and transport properties
   use a **log transform**
   (:math:`\rho = \exp(\sum_k S_k U_k V_k)`) so their multi-decade
   dynamic range fits the SVD's smoothness assumption.

3. **Critical-patch HEOS fallback.**  Near the critical point the
   rank-truncated SVD cannot resolve the cusp in
   :math:`(\partial \rho / \partial p)_h` — a single dense region
   needs unbounded rank.  SVDSBTL identifies a small rectangular
   patch in :math:`(T, p)` around the critical point (auto-calibrated
   per fluid; see :ref:`Critical-patch bbox <SVDSBTL-critpatch>`
   below) and routes queries inside the patch directly through the
   source backend.  Outside the patch the SVD owns the answer; this
   keeps the SVD's per-cell rank low everywhere else.

Queries through the public ``AbstractState`` interface dispatch to
whichever of these three handles the input pair: the atlas locates the
region, the per-region SVD returns the property, and the patch
overrides when the query is in the critical neighbourhood.  Two-phase
dome queries lever-rule-blend between the source backend's saturation
endpoints (or :doc:`SuperAncillary` when the source exposes one).

Constructor & supported pairs
=============================

SVDSBTL is constructed via the factory string
``SVDSBTL&<source>``, where ``<source>`` is the truth-source backend
that the SVD samples from at table-build time:

.. code-block:: python

    import CoolProp.CoolProp as CP

    AS = CP.AbstractState("SVDSBTL&IF97", "Water")     # IAPWS G13-15
    AS = CP.AbstractState("SVDSBTL&HEOS", "Methane")   # multi-fluid HEOS
    AS = CP.AbstractState("SVDSBTL&REFPROP", "Water")  # REFPROP truth

The currently-supported source backends are ``HEOS``, ``REFPROP``, and
``IF97`` (water only).  ``SVDSBTL`` *without* a source backend
(``SVDSBTL`` alone) is intentionally not supported — there is no
default truth source, since the choice changes both the sampled
property values and the cache filename.

Two input pairs are currently tabulated per fluid:

* ``HmassP_INPUTS`` (the dominant industrial use case)
* ``PT_INPUTS`` (where the user already has :math:`(p, T)` and wants to
  skip the inverse :math:`T(p, h)` solve)

``HmolarP_INPUTS`` is served from the ``HmassP_INPUTS`` table by
multiplying by molar mass.  Two-phase ``PQ_INPUTS`` and ``QT_INPUTS``
route directly through the source's saturation line — no separate
table is required.

A planned follow-up adds ``HmassSmass`` and ``Dmass-Umass`` for the
remaining G13-15 input-pair coverage (Tables 14-17); see bd
``CoolProp-phv``.

The on-disk table cache
=======================

The first construction for a given ``(fluid, source, input_pair,
options)`` tuple samples ~160,000 :math:`(p, h)` or :math:`(p, T)`
points from the source backend, runs SVD on each region's property
matrix, and persists the result as a zlib-compressed binary file
under

.. code-block:: text

    $HOME/.CoolProp/SVDTables/<Fluid>.<Source>.<InputPair>.<OptHash>.svd.bin.z

Subsequent constructions for the same tuple load from this cache.

* **Build time:** ~10-60 s per (fluid, input pair) on a 2024-era laptop,
  HEOS-source.  REFPROP-source is several times slower (REFPROP flash
  per sample point).
* **Disk footprint:** ~5-15 MB per file depending on region count.
* **Cache invalidation:** the on-wire revision is checked on load.
  Changes to the serializer format, sampling code, or preset geometry
  bump the revision and force a rebuild.  Source-code paths that
  affect the cached output are hashed into the GHA cache key in
  ``.github/workflows/test_catch2.yml`` so CI rebuilds occur
  automatically.

The cache directory is `XDG_DATA_HOME`-aware on Linux and uses
``~/.CoolProp/SVDTables/`` by default elsewhere.  The path is
configurable via the ``ALTERNATIVE_SVDTABLES_DIRECTORY`` configuration
key (see :ref:`configuration`).

.. _SVDSBTL-regimes:

Three performance regimes
=========================

SVDSBTL exposes three distinct call patterns with very different
per-call cost.  Picking the right one matters far more than picking
the backend.

1. **AbstractState.update — reuse the instance (canonical path).**
   Construct the ``AbstractState`` once, then call ``update`` /
   property accessors in a loop.  Table load is amortised over the
   full workload:

   .. ipython::

       In [1]: import CoolProp.CoolProp as CP

       In [1]: from CoolProp import AbstractState

       In [1]: import numpy as np

       In [1]: AS = AbstractState("SVDSBTL&IF97", "Water")

       In [1]: h_arr = np.linspace(2.0e5, 3.0e6, 5)

       In [1]: p_arr = np.full_like(h_arr, 1.0e6)

       In [1]: def loop():
          ...:     for h, p in zip(h_arr, p_arr):
          ...:         AS.update(CP.HmassP_INPUTS, h, p)
          ...:         AS.rhomass(); AS.T()

       In [1]: %timeit loop()

   Native C++ per-call wall time is ~200 ns; the bulk of the
   Python-visible cost is the wrapper marshalling, and even in Python
   the per-call cost is hundreds of nanoseconds — far faster than
   HEOS's :math:`(p, h)` flash.  This is what most user code should
   look like.

2. **fast_evaluate (batched, no instance mutation).**
   ``AbstractState.fast_evaluate(input_pair, val1[], val2[],
   outputs[], out_buffer, status)`` writes N property values for M
   probes directly into a caller-supplied buffer.  No per-point
   ``AbstractState`` cache mutation, no Python wrapper marshalling
   per probe, and shared :math:`(p, h)` locate / Hermite-basis setup
   across the N requested properties:

   .. ipython::

       In [1]: import CoolProp.CoolProp as CP

       In [1]: from CoolProp import AbstractState

       In [1]: import numpy as np

       In [1]: AS = AbstractState("SVDSBTL&IF97", "Water")

       In [1]: N = 10_000

       In [1]: h = np.random.uniform(1e5, 3e6, N)

       In [1]: p = np.random.uniform(1e5, 3e7, N)

       In [1]: outputs = np.array([CP.iT, CP.iDmass, CP.iSmass], dtype=np.int32)

       In [1]: out = np.empty((N, 3)); status = np.empty(N, dtype=np.int32)

       In [1]: %timeit AS.fast_evaluate(CP.HmassP_INPUTS, h, p, outputs, out, status)

   Per-point wall time: **~340-400 ns** for a 4-output query on Apple
   Silicon (10k-probe batch).  Marginal cost per added surface output
   ~57 ns after locate / basis-weight setup is paid once per probe.
   For CFD-scale workloads (≥100k probes per timestep) this is the
   only sensible API.

3. **PropsSI — the high-level wrapper (off by default).**
   ``PropsSI`` is the high-level entry point that resolves the
   backend and constructs an ``AbstractState`` on every call.  For
   SVDSBTL that means a zlib-decompressed ``.svd.bin.z`` load (~80 ms
   per ``AbstractState`` construction).  For a one-off interactive
   query this is fine; for a throughput workload it is catastrophic
   — *every* ``PropsSI`` call pays the table-load cost again.

   SVDSBTL therefore opts *out* of ``PropsSI`` by default.  Enable
   only for interactive / scripting use where the construction cost
   is irrelevant:

   .. code-block:: python

       CP.set_config_bool(CP.ALLOW_SVDSBTL_IN_PROPSSI, True)
       CP.PropsSI("D", "T", 500.0, "P", 1e6, "SVDSBTL&IF97::Water")

   This gate is :ref:`enforced architecturally <SVDSBTL-config>` — the
   same pattern the bicubic and TTSE backends use, for the same
   reason.

The :ref:`profile figure on the IF97 page <IF97-Conformance>` shows
all three regimes side-by-side, against HEOS and IF97 baselines.

.. _SVDSBTL-critpatch:

Critical-patch bbox
===================

Near the critical point the rank-:math:`r` SVD's reconstruction error
diverges — the :math:`(\partial \rho / \partial p)_h` cusp at the
critical singularity needs unbounded rank to represent, and even a
generous rank like :math:`r = 20` leaves multi-percent error in
:math:`\rho` and ten-percent-class error in :math:`w` within a
fraction of :math:`T_c` from the critical point.  SVDSBTL covers this
slice with a small axis-aligned rectangular **patch** in
:math:`(T, p)`: every query whose state lands inside the patch
forwards to the source backend (HEOS / REFPROP / IF97) for the
calculation; everything outside the patch goes through the SVD.

**Patch shape per fluid.**  The patch is computed automatically at
the first backend construction for a given (fluid, source backend,
options) combination, via a binary-search shrink loop:

1. Start from the Water-sized default
   :math:`[0.95, 1.05] T_c \times [0.75, 1.15] p_c`.
2. For each of the four axes :math:`(T_\mathrm{lo},
   T_\mathrm{hi}, p_\mathrm{lo}, p_\mathrm{hi})`, binary-search-shrink
   the multiplier toward 1.0 (= the critical point) until the SVD's
   reconstruction at a strip of probes *just outside* the candidate
   patch boundary passes the calibration error budget against the
   source backend (1% relative in :math:`\rho, h, s`; 5% relative in
   :math:`w`).
3. The result is persisted to a sidecar file alongside the
   ``.svd.bin.z`` table caches:

   .. code-block:: text

       $HOME/.CoolProp/SVDTables/<Fluid>.<Source>.critpatch.<OptHash>.bin

   Subsequent constructions load the cached multipliers; the
   calibration probe is not re-run.

The calibrator never widens the patch beyond the Water default — so
Water keeps its tested bbox and other fluids tighten only as much as
their critical-singularity footprint allows.  Typical results
(post-CoolProp-5ni / -dxd):

.. list-table::
   :header-rows: 1
   :widths: 22 18 18 18 18

   * - Fluid (source)
     - :math:`T_\mathrm{lo} / T_c`
     - :math:`T_\mathrm{hi} / T_c`
     - :math:`p_\mathrm{lo} / p_c`
     - :math:`p_\mathrm{hi} / p_c`
   * - Water (HEOS)
     - 0.970
     - 1.050
     - 0.996
     - 1.150
   * - Water (IF97)
     - 0.984
     - 1.027
     - 0.996
     - 1.150
   * - CarbonDioxide (HEOS)
     - 0.961
     - 1.050
     - 0.863
     - 1.150
   * - R134a (HEOS)
     - 0.995
     - 1.013
     - 0.996
     - 1.080

These numbers are regenerated at backend construction; the exact
values are sensitive to the SVD grid size and rank (and so move when
the options blob changes).  ``CoolProp-3c4`` tracks the residual
accuracy ceiling at extreme thin-support cells *outside* the patch,
which the calibrator deliberately ignores (axis-aligned bboxes can't
cover those scattered failures).

For the **HmassP_INPUTS** input pair the :math:`(p, h)` bbox is
derived from the :math:`(T, p)` patch by walking the perimeter
through the source backend's :math:`h(T, p)` and taking the
axis-aligned envelope of the resulting :math:`h` values — that way
any :math:`(T, p)` in the canonical patch round-trips into a
:math:`(h, p)` that the patch test also fires on.  For Water/HEOS
the resulting :math:`(p, h)` bbox is roughly
:math:`p \in [22.0, 25.4]` MPa, :math:`h \in [1.55, 2.93]` MJ/kg.

**Overrides.**  Two escape hatches:

* ``critical_patch.bbox`` in the options JSON — set explicit
  :math:`(T_\mathrm{lo}, T_\mathrm{hi}, p_\mathrm{lo},
  p_\mathrm{hi})` multipliers (skips auto-calibration entirely).
  Useful for pinning the patch to a known-good shape across versions.
* ``critical_patch.mode = "off"`` — disables the patch.  Queries near
  the critical point will return the SVD's rank-truncated reconstruction
  with the corresponding accuracy degradation.  Useful for benchmark
  studies of the SVD itself.

* ``critical_patch.source`` — override the source backend that serves
  in-patch queries.  E.g. ``SVDSBTL&IF97`` with
  ``critical_patch.source = "HEOS"`` uses IF97 for the SVD truth
  source but HEOS for in-patch queries (which gives full IAPWS-95
  accuracy at the critical singularity rather than IF97's R3
  formulation).

The ``set_critical_bbox_multipliers(T_lo, T_hi, p_lo, p_hi)`` C++
entry point also lets callers override the patch at runtime without
rebuilding the backend.

.. _SVDSBTL-config:

Configuration keys
==================

``ALLOW_SVDSBTL_IN_PROPSSI`` (default ``false``)
    Gate for ``PropsSI`` routing.  See :ref:`Three performance regimes
    <SVDSBTL-regimes>`.  Set to ``true`` for interactive
    use; leave ``false`` for throughput code.

``ALTERNATIVE_SVDTABLES_DIRECTORY`` (default empty)
    Override the cache directory.  Useful when ``$HOME`` is read-only
    or on shared workstations where a centrally-managed cache is
    preferable.

Backend options JSON (``factory("SVDSBTL&HEOS", "Water",
'{"grid": {"NT": 200, "NR": 800, "rank": 20}}')``) lets per-instance
overrides tune the SVD grid and rank; see :doc:`BackendOptions`.  Any
non-default options change the table's ``OptHash`` and therefore its
cache filename, so multiple grid sizes for the same fluid coexist
peacefully on disk.

.. _SVDSBTL-decision-guide:

When to use SVDSBTL
===================

SVDSBTL pays off when one of these is true:

* You evaluate a fixed-formulation EOS at **>10k probes per
  AbstractState instance**.  The first construction pays a
  ~80 ms table-load cost (or, for a cold cache, a few seconds of
  table-build cost); both amortise over the subsequent queries.
* You need :math:`(p, h)` lookup (typical for power-cycle,
  refrigeration, CFD); HEOS's :math:`(T, \rho)` formulation makes
  :math:`(p, h)` inputs 5-10× slower than :math:`(p, T)`.
* You want conformance to IAPWS G13-15 for water — both the
  *backward* :math:`T(p, h)` equations and the *forward* equations
  (SVDSBTL&IF97 returns forward-consistent values via the
  TOMS748-polished SVD).
* You need a non-pre-tabulated input pair (e.g. :math:`(p, s)`-driven
  turbine stages) at IF97 conformance: SVDSBTL's atlas + per-region
  SVD generalises beyond the input pair that IF97 ships native
  backward equations for.

SVDSBTL is **not** the right pick when:

* Your workload is interactive single-shot ``PropsSI`` calls.  Use
  HEOS or IF97 directly — the ~80 ms per-call table-load overhead
  dominates.
* You need transport properties for a fluid whose source backend
  lacks them.  SVDSBTL just samples the source; if the source's HEOS
  model has no :math:`\eta` / :math:`\lambda` correlation, neither
  does the SVDSBTL table for that fluid.  (This is a source-backend
  limitation, not an SVDSBTL one — but it bites SVDSBTL users
  because they may pick SVDSBTL specifically to skip the
  per-call source flash, then discover the property they want isn't
  in the table.)

Accuracy envelope
=================

* **SVDSBTL&IF97::Water**: conformant to IAPWS G13-15 Tables 8-13 for
  :math:`T, v, s, w` in regions R1, R2, R5; R3 close to but not
  uniformly inside budget (the critical-patch fallback covers the
  near-critical R3 cells).  Transport properties :math:`\eta,
  \lambda` exceed G13-15 budgets in R3 — those budgets are
  intrinsically tight (:math:`10^{-5}` relative) and the rank-20 SVD
  has known headroom for tighter ranks.  See
  :ref:`IF97-Conformance` for the per-region fail maps and exact
  numbers.

* **SVDSBTL&HEOS::<Fluid>**: ~10^-3 relative error on
  :math:`\rho, h, s` across the full subcritical envelope.
  Peak errors at extreme region corners (Hermite thin-support cells)
  reach ~1% — investigating tighter axis bounds; see bd
  ``CoolProp-3c4``.

* **SVDSBTL&REFPROP::<Fluid>**: matches REFPROP truth on tested
  probes within REFPROP's own numerical-noise floor on the sat
  curves; same single-phase envelope as HEOS-source.  Two-phase
  endpoints currently fall through to REFPROP's PQ flash (~ms per
  probe); ``CoolProp-077`` will replace that with an on-disk
  saturation-surrogate spline.

Deviation plot
==============

The figure below shows the density deviation of ``SVDSBTL&HEOS`` from
the underlying HEOS truth source over the full :math:`(h, p)`
envelope of R245fa.  Sample density 20,000 random :math:`(h, p)`
points (matching the :ref:`tabular-interpolation accuracy figure
<tabular_interpolation>` for the BICUBIC / TTSE backends).

.. plot::

    import matplotlib
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.ticker
    import matplotlib.pyplot as plt
    import numpy as np
    import random

    import CoolProp
    import CoolProp.CoolProp as CP

    Ref = 'R245fa'
    EOS = CoolProp.AbstractState('HEOS', Ref)
    SVD = CoolProp.AbstractState('SVDSBTL&HEOS', Ref)

    HHH, PPP, EEE = [], [], []
    cNorm = colors.LogNorm(vmin=1e-12, vmax=10)
    for _ in range(20000):
        h = random.uniform(150000, 590000)
        p = 10 ** random.uniform(np.log10(100000), np.log10(7000000))
        try:
            EOS.update(CoolProp.HmassP_INPUTS, h, p)
            SVD.update(CoolProp.HmassP_INPUTS, h, p)
            err = abs(SVD.rhomolar() / EOS.rhomolar() - 1) * 100
        except Exception:
            continue
        HHH.append(h); PPP.append(p); EEE.append(err)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes((0.13, 0.13, 0.74, 0.78))
    sc = ax.scatter(HHH, PPP, s=8, c=EEE, edgecolors='none',
                    cmap='jet', norm=cNorm)
    ax.set_yscale('log')
    ax.set_xlabel('Enthalpy [J/kg]')
    ax.set_ylabel('Pressure [Pa]')
    ax.set_title('SVDSBTL&HEOS density deviation — ' + Ref)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label(r'$|\rho_{\mathrm{SVDSBTL}} / \rho_{\mathrm{HEOS}} - 1|\ \%$')

For the IF97-sourced SVDSBTL on water, the per-property fail-map
figure in :ref:`IF97-Conformance` shows the equivalent — same
:math:`(h, p)`-coloured deviation pattern but classified against the
IAPWS G13-15 conformance budgets rather than reported as raw
relative error.

For the underlying SVD compression math — why a rank-:math:`r` SVD
of a smooth function on a rectangle gives :math:`O(r^{-\infty})`
error decay, why the natural cubic spline beats finite-difference
slopes for the Hermite per-mode interpolation, and the
``f(x, y) = \exp(-x^2 - y^2)`` rank-2 worked example — see
:doc:`SVDComponents`.

Programmatic access
===================

Beyond the high-level ``AbstractState`` API, the C++ entry point
``CoolProp::SVDSBTLBackend`` exposes:

* ``fast_evaluate`` for batched cache-bypassing evaluation (also
  bound to Python — see the example above)
* ``registered_input_pairs()`` for introspection (which surfaces did
  this instance load)
* ``set_critical_bbox_multipliers(T_lo_mult, T_hi_mult, p_lo_mult,
  p_hi_mult)`` for ad-hoc critical-patch tuning

All other ``AbstractState`` methods (``rhomass``, ``hmass``, ``smass``,
``viscosity``, ``conductivity``, ``T_critical``, ``Ttriple``, etc.)
work identically to other backends.
