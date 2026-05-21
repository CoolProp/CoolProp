.. _SVDSBTL:

*************************************************
SVDSBTL --- SVD-Compressed Tabular Lookup Backend
*************************************************

SVDSBTL is a tabular-lookup backend that combines an atlas of
analytically-bounded regions with a low-rank SVD of each region's
property surfaces.  It produces single-digit-microsecond property
evaluations at IAPWS G13-15 conformance for water, and a single-digit
percent accuracy ceiling for the multi-fluid HEOS-backed presets, at a
disk footprint of ~10 MB per (fluid, input pair, source backend).

This page is the **user-facing** page: how to use the backend, what it
guarantees, and when to pick it.  For the underlying SVD math
(rank-:math:`r` reconstruction, cubic Hermite slopes, axis transforms),
see :doc:`SVDComponents`.

.. _SVDSBTL-architecture:

What it is
==========

A traditional bicubic table covers the full :math:`(p, h)` or
:math:`(p, T)` envelope with one coarse grid and stores the source
backend's property values at each cell.  SVDSBTL replaces that with
three layered components:

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
   continuity across the region.  Density / transport properties
   ride a log transform (:math:`\rho = \exp(\sum_k S_k U_k V_k)`) so
   their multi-decade dynamic range fits the SVD's smoothness
   assumption.

3. **Critical-patch HEOS fallback.**  Near the critical point the
   rank-truncated SVD cannot resolve the cusp in
   :math:`(\partial \rho / \partial p)_h` — a single dense region
   needs unbounded rank.  SVDSBTL identifies a small rectangular
   patch in :math:`(T, p)` around the critical point (default
   :math:`[0.95, 1.05] \cdot T_c \times [0.75, 1.15] \cdot p_c` for
   Water) and routes queries inside the patch directly through the
   source backend.  Outside the patch the SVD owns the answer.  This
   keeps the SVD's per-cell rank low everywhere it is responsible.

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

A follow-up adds ``HmassSmass`` and ``Dmass-Umass`` for the remaining
G13-15 input-pair coverage (Tables 14-17); see bd
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

1. **PropsSI (low-level wrapper, off by default).**
   ``PropsSI`` rebuilds an ``AbstractState`` instance per call,
   including a zlib-decompressed ``.svd.bin.z`` load (~80 ms per
   ``AbstractState`` construction).  For a one-off query this is
   cheap-enough; for a throughput workload it is catastrophic —
   *every* ``PropsSI`` call pays the table-load cost again.

   SVDSBTL therefore opts *out* of ``PropsSI`` by default.  Enable
   only for interactive / scripting use where the construction cost
   is irrelevant:

   .. code-block:: python

       CP.set_config_bool(CP.ALLOW_SVDSBTL_IN_PROPSSI, True)
       CP.PropsSI("D", "T", 500.0, "P", 1e6, "SVDSBTL&IF97::Water")

   This gate is :ref:`enforced architecturally <SVDSBTL-config>` — the
   same pattern the bicubic and TTSE backends use, for the same
   reason.

2. **AbstractState.update (reuse the instance).**
   Construct the ``AbstractState`` once, then call ``update`` /
   property accessors in a loop.  Table load is amortised over the
   full workload:

   .. code-block:: python

       AS = CP.AbstractState("SVDSBTL&IF97", "Water")
       for h, p in probes:
           AS.update(CP.HmassP_INPUTS, h, p)
           rho.append(AS.rhomass())
           T.append(AS.T())

   Per-call wall time: **~1-2 μs** on a modern CPU (includes Python
   wrapper marshalling; native C++ closer to ~200 ns).  This is
   what most user code should look like.

3. **fast_evaluate (batched, no instance mutation).**
   ``AbstractState.fast_evaluate(input_pair, val1[], val2[],
   outputs[], out_buffer, status)`` writes N property values for M
   probes directly into a caller-supplied buffer.  No per-point
   ``AbstractState`` cache mutation, no Python wrapper marshalling
   per probe, and shared :math:`(p, h)` locate / Hermite-basis setup
   across the N requested properties:

   .. code-block:: python

       import numpy as np
       AS = CP.AbstractState("SVDSBTL&IF97", "Water")
       N = 10_000
       h = np.random.uniform(1e5, 3e6, N)
       p = np.random.uniform(1e5, 3e7, N)
       outputs = np.array([CP.iT, CP.iDmass, CP.iSmass], dtype=np.int32)
       out = np.empty((N, 3))
       status = np.empty(N, dtype=np.int32)
       AS.fast_evaluate(CP.HmassP_INPUTS, h, p, outputs, out, status)

   Per-point wall time: **~340-400 ns** for a 4-output query on Apple
   Silicon (10k-probe batch).  Marginal cost per added surface output
   ~57 ns after locate / basis-weight setup is paid once per probe.
   For CFD-scale workloads (≥100k probes per timestep) this is the
   only sensible API.

The :ref:`profile figure on the IF97 page <IF97-Conformance>` shows
all three regimes side-by-side, against HEOS and IF97 baselines.

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

* You evaluate a fixed-formulation EOS at **>10k probes per process
  lifetime**.  The table-build amortises out below that point.
* You need :math:`p, h` lookup (typical for power-cycle, refrigeration,
  CFD); HEOS's :math:`(T, \rho)` formulation makes :math:`(p, h)`
  inputs 5-10× slower than :math:`(p, T)`.
* You want conformance to IAPWS G13-15 for water (not just the
  G13-15 *backward* :math:`T(p, h)` accuracy of IF97, but the
  *forward* equations as well — SVDSBTL&IF97 inherits IF97's
  :math:`±25` mK :math:`T(p, h)` floor by construction).
* You need to use IF97 forward equations from non-pre-IF97 input
  coordinates (e.g. :math:`(p, s)`-driven turbine stages) without
  paying IF97's :math:`(p, h) \to T \to (p, T)` inversion cost.

SVDSBTL is **not** the right pick when:

* Your workload is interactive single-shot ``PropsSI`` calls.  Use
  HEOS or IF97 directly.
* You need accuracy in the critical region better than HEOS's
  :math:`±0.5` mK (SVDSBTL's critical-patch HEOS fallback inherits
  HEOS's accuracy in that slice, so this is rarely a hard limit, but
  if you are running ultra-precise EOS verification, the SVD path is
  not the right tool).
* You need transport properties for a fluid whose HEOS model lacks
  them — SVDSBTL's table will not include :math:`\eta` or
  :math:`\lambda` for that fluid, and queries will return NaN.

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
