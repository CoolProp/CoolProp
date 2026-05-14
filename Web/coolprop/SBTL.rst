.. _sbtl_interpolation:

******************************************
SBTL Coordinate-Aligned Spline Interpolation
******************************************

The SBTL (Spline-Based Table Look-Up) backend provides fast thermodynamic
property evaluation for **pure fluids only** via two coordinate-aligned
cubic B-spline tables that put the saturation curve on a coordinate axis.
This eliminates the cell-straddling errors that otherwise dominate accuracy
near the dome.

SBTL supports only ``PT_INPUTS``, ``HmolarP_INPUTS``, and ``HmassP_INPUTS``
as input pairs.  Mixtures (including pseudo-pure refrigerant blends like
R407C, R410A) are rejected at construction; use the BICUBIC or TTSE
backends for tabular mixture support.

On-disk caching
---------------

The first call to ``factory("SBTL&HEOS", "<fluid>")`` builds the six
coordinate-aligned tables (three normph + three normpt) from scratch
— this is HEOS-bound at every grid corner and takes ~20–30 s per
fluid.  The tables (corner values *and* the per-cell Hermite bicubic
coefficient vectors) are then msgpack-packed and zlib-compressed to
``$HOME/.CoolProp/Tables/HelmholtzEOSBackend(<fluid>[<mole_frac>])/``
in six files named ``sbtl_norm{ph,pt}_{liquid,vapor,super}.bin.z``.
Subsequent calls deserialise from those files in ~1–2 s, skipping
the full HEOS sweep entirely.  Per-fluid cache size is ≈ 200 MB
compressed.

The cache directory is shared with BICUBIC and TTSE for the same
fluid; the file-name prefixes (``single_phase_log*`` vs ``sbtl_*``)
prevent collisions.  SBTL is subject to the same
``MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB`` warning and the same
``ALTERNATIVE_TABLES_DIRECTORY`` config override as the other tabular
backends.  Cache files keyed by the current ``TABULAR_NX``/``TABULAR_NY``
grid size are evicted on mismatch and rebuilt at the requested
resolution.

Backbone interpolation
----------------------

For each cell :math:`[x_i, x_{i+1}] \times [y_j, y_{j+1}]`, every property
:math:`z` is represented as a 16-coefficient bi-cubic polynomial in the
cell-local coordinates:

.. math::

    z(\hat{x}, \hat{y}) = \sum_{m=0}^{3} \sum_{n=0}^{3} a_{mn}\, \hat{x}^m \hat{y}^n

The per-cell coefficient set :math:`a_{mn}` is built by a **two-pass overlay**
inside :func:`SBTLBackend::build_bspline_coeffs`:

#. **Pass 1 — global C\ :sup:`2` cubic B-spline backbone** (Kunick C\ :sup:`2`).
   A single 2-D not-a-knot natural cubic B-spline is fit through the
   :math:`N_x \times N_y` matrix of property values at the grid corners,
   and each cell's :math:`a_{mn}` is extracted from the global fit.
   This pass writes *every* cell.  Continuity across cell boundaries is
   C\ :sup:`2` because the underlying B-spline is.

#. **Pass 2 — C\ :sup:`1` Hermite-bicubic overlay**
   (:func:`SBTLBackend::build_normph_hermite_alphas` /
   ``build_normpt_hermite_alphas``).  For every cell, if all four corners
   yield finite HEOS values, finite first partials (:math:`\partial \rho /
   \partial h`, :math:`\partial \rho / \partial p`, etc.) **and** a finite
   mixed second partial :math:`\partial^2 \rho / \partial h \partial p`,
   the cell's :math:`a_{mn}` from pass 1 is **replaced** by a Hermite
   bicubic that exactly interpolates the corner values + their first
   derivatives (with the cross-deriv term filling the 16-th coefficient).
   Continuity across cell boundaries is then C\ :sup:`1` (values + first
   derivatives match at shared corners by construction; second
   derivatives may jump).

So the final per-cell polynomial is:

* **C\ :sup:`1` Hermite bicubic** in the bulk of every coordinate-aligned
  subcritical region — including the cells touching the saturation
  boundary, which use the H-superancillary's analytic
  :math:`\mathrm{d}h_{\mathrm{sat}}/\mathrm{d}p` to fill the row-boundary
  derivative.
* **C\ :sup:`2` B-spline backbone**, untouched by the overlay, anywhere
  pass 2 declines to run — typically a thin ring around the critical
  point where HEOS partials return non-finite values, and any cell where
  the EOS rejected the corner input.  These cells fall through to the
  HEOS fallback box at lookup time for high accuracy; the cubic-B-spline
  fallback handles queries outside the box.

(The earlier version of the docs described only pass 1 and was technically
incomplete: most of the table the user ever queries is the C\ :sup:`1`
Hermite overlay, not the C\ :sup:`2` B-spline backbone.)

The B-spline backbone uses the same coordinate-aligned :math:`(\hat x,
\log p)` grid for both pass-1 and pass-2 polynomials, with three region
flavours:

* **PH tables** (``NormalizedPHTable``) — left/right column boundaries
  are ``h_lo(P)`` (cold-isotherm or saturated-vapor curve) and
  ``h_hi(P)`` (hot-isotherm or saturated-liquid curve), so column
  boundaries lie exactly on the saturation curve.
* **PT tables** (``NormalizedPTTable``) — same construction with
  ``T_lo(P)`` and ``T_hi(P)``.

The coordinate-aligned tables guarantee that no cell straddles the
saturation curve: ``xnorm = 1`` for the LIQUID region is exactly the
saturated-liquid curve, ``xnorm = 0`` for the VAPOR region is exactly the
saturated-vapor curve.  This is the central design idea borrowed from the
IAPWS G13-15 guideline.

``HmassP_INPUTS`` / ``HmolarP_INPUTS``
---------------------------------------

For pure fluids, ``(h, P)`` queries route through the appropriate
``NormalizedPHTable``:

1. If ``P >= P_crit``: use the SUPER table.
2. Otherwise compute ``h_sat,L(P)`` and ``h_sat,V(P)`` via the H-superancillary
   (machine precision); route to LIQUID, VAPOR, or two-phase based on
   ``h`` vs the sat bounds.
3. A small HEOS-fallback box around ``(h_crit, p_crit)`` (``±0.1 MPa × ±30 % h_crit``)
   handles the cusp where no smooth polynomial backbone can faithfully
   reproduce :math:`\rho(h, P)` at the critical point.

``PT_INPUTS``
-------------

For pure fluids, ``(T, P)`` queries route through the
``NormalizedPTTable`` family in the same pattern:

1. If ``P >= P_crit``: use the SUPER table.
2. Otherwise compute ``T_sat(P)`` via ``PQ_INPUTS`` (cached single-deep);
   route to LIQUID if ``T <= T_sat,L``, VAPOR if ``T >= T_sat,V``, or
   throw if strictly inside the dome (``PT_INPUTS`` is ambiguous there).
3. A larger HEOS-fallback box around ``(T_crit, p_crit)``
   (``[0.75·p_c, 1.75·p_c] × [0.90·T_c, 1.30·T_c]``) handles the cusp and
   the supercritical-shoulder where :math:`\rho(T, P)` drops steeply just
   above ``p_crit``.

The PT path puts the saturation curve on a coordinate axis (``xnorm = 1`` for
LIQUID, ``xnorm = 0`` for VAPOR) so cells never straddle the dome.  Random-PT
median accuracy is at the 1e-8 floor for well-conditioned states.

Phase specification on the saturation curve
---------------------------------------------

A ``PT_INPUTS`` query at ``T = T_sat(P)`` is formally ambiguous between
saturated liquid and saturated vapor.  Specify the phase before the
update to disambiguate:

.. ipython::
    :okexcept:

    In [0]: import CoolProp

    In [1]: SBTL = CoolProp.AbstractState("SBTL&HEOS", "CO2")

    In [2]: HEOS = CoolProp.AbstractState("HEOS", "CO2")

    # Find the exact saturation temperature at the chosen pressure.
    In [3]: P = 15.8e5

    In [4]: HEOS.update(CoolProp.PQ_INPUTS, P, 0.0); T_sat = HEOS.T()

    In [5]: round(T_sat, 3)

    # Default: T = T_sat at this P routes to the LIQUID table.
    In [6]: SBTL.update(CoolProp.PT_INPUTS, P, T_sat); round(SBTL.rhomass(), 1)

    # Same query, but specify gas → routes to the VAPOR table.
    In [7]: SBTL.specify_phase(CoolProp.iphase_gas)

    In [8]: SBTL.update(CoolProp.PT_INPUTS, P, T_sat); round(SBTL.rhomass(), 1)

    In [9]: SBTL.unspecify_phase()

The mapping from ``imposed_phase_index`` to the table that handles the
query:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Imposed phase
     - Routing
   * - ``iphase_not_imposed``
     - Auto: route by T vs T_sat (LIQUID at T ≤ T_sat,L; VAPOR at T ≥ T_sat,V; SUPER if p ≥ p_c)
   * - ``iphase_liquid``
     - LIQUID table at p < p_c; SUPER at p ≥ p_c
   * - ``iphase_gas``
     - VAPOR  table at p < p_c; SUPER at p ≥ p_c
   * - ``iphase_supercritical``
     - SUPER  table regardless of p
   * - ``iphase_supercritical_liquid``
     - SUPER  table (caller knows the T < T_c / p ≥ p_c quadrant)
   * - ``iphase_supercritical_gas``
     - SUPER  table (T > T_c, p ≥ p_c)
   * - ``iphase_critical_point``
     - SUPER  table; the critbox fallback catches the cusp
   * - ``iphase_twophase``
     - Rejected: PT in the dome is ambiguous; use PQ_INPUTS or specify_phase(iphase_liquid|iphase_gas)

The same contract applies to ``HmolarP_INPUTS`` / ``HmassP_INPUTS`` at
``h = h_sat,L(P)`` or ``h = h_sat,V(P)``.

Other input pairs
-----------------

SBTL currently supports only ``PT_INPUTS``, ``HmolarP_INPUTS``, and
``HmassP_INPUTS``.  Other input pairs (``DmolarT_INPUTS``,
``DmolarUmolar_INPUTS``, ``SmolarP_INPUTS``, etc.) throw cleanly; the caller
should switch to the HEOS backend or to BICUBIC / TTSE.  The Kunick-style
``(D, U)`` and ``(D, T)`` direct-flash tables documented in
``dev/sbtl_du_kunick_redesign_plan.md`` are queued for a follow-up PR.

Mixtures
--------

Mixtures (real or pseudo-pure) are rejected at SBTL construction.  The
coordinate-aligned tables require a 1-D saturation locus
(``T_sat,L(P) = T_sat,V(P)``) which only holds for single-component fluids;
mixture support would route through phase-envelope bubble/dew lines and is
planned as a follow-up.  For tabular mixture properties, use BICUBIC or TTSE.

Accuracy comparison
-------------------

Density obtained for R245fa using the equation of state, TTSE, BICUBIC,
and SBTL:

.. ipython::
    :okexcept:

    In [0]: import CoolProp

    In [1]: HEOS = CoolProp.AbstractState("HEOS",         "R245fa")

    In [2]: TTSE = CoolProp.AbstractState("TTSE&HEOS",    "R245fa")

    In [3]: BICU = CoolProp.AbstractState("BICUBIC&HEOS", "R245fa")

    In [4]: SBTL = CoolProp.AbstractState("SBTL&HEOS",    "R245fa")

    In [5]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300); TTSE.update(CoolProp.PT_INPUTS, 101325, 300); BICU.update(CoolProp.PT_INPUTS, 101325, 300); SBTL.update(CoolProp.PT_INPUTS, 101325, 300)

    In [6]: print(HEOS.rhomolar(), TTSE.rhomolar(), BICU.rhomolar(), SBTL.rhomolar())

The figures below quantify density error in the single-phase region for
five fluids spanning very different molecular size and critical-point
locations: **Argon** (cryogenic, simple), **CarbonDioxide** (industrial
reference, narrow liquid range), **Water** (high :math:`p_{\rm crit}`,
strong hydrogen bonding), **R245fa** (refrigerant), and **D6** (siloxane,
low :math:`p_{\rm crit}\approx 0.96\,{\rm MPa}`).

For each fluid we draw 3,000 random points :math:`(T, p)` uniformly in
:math:`T` and log-uniformly in :math:`p` over the table's full envelope,
compute :math:`h={\rm HEOS}(T,p)`, then look up :math:`\rho` from each
tabular backend with both ``HmassP_INPUTS`` and ``PT_INPUTS``.  Error is
:math:`|\rho_{\rm table}/\rho_{\rm EOS}-1|`.  Two-phase samples are
rejected at HEOS lookup.

``HmassP_INPUTS`` error map
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SBTL routes through ``NormalizedPHTable``; BICUBIC uses Hermite bi-cubic
on a ``logph`` grid; TTSE uses bilinear Taylor expansion.  SBTL eliminates
the cell-straddling streaks along the saturation curve that BICUBIC shows,
at the cost of larger error in a small band immediately below
:math:`p_{\rm crit}` (most visible for R245fa).  See the discussion below
the PT plot for the mechanism.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    FLUIDS = ['Argon', 'CarbonDioxide', 'Water', 'R245fa', 'D6']
    N = 6000
    fig, axes = plt.subplots(len(FLUIDS), 3, figsize=(11, 14))
    fig.subplots_adjust(left=0.08, right=0.90, top=0.97, bottom=0.04,
                        wspace=0.32, hspace=0.50)
    cNorm = colors.LogNorm(vmin=1e-10, vmax=10)
    cmap = plt.get_cmap('jet')
    for ir, Ref in enumerate(FLUIDS):
        EOS  = CoolProp.AbstractState('HEOS',         Ref)
        TTSE = CoolProp.AbstractState('TTSE&HEOS',    Ref)
        BICU = CoolProp.AbstractState('BICUBIC&HEOS', Ref)
        SBTL = CoolProp.AbstractState('SBTL&HEOS',    Ref)
        Tmin  = CP.PropsSI(Ref, 'Tmin') + 0.5
        Tmax  = CP.PropsSI(Ref, 'Tmax') - 0.5
        Tcrit = CP.PropsSI(Ref, 'Tcrit')
        pmin  = CP.PropsSI(Ref, 'ptriple') * 1.05
        pmax  = CP.PropsSI(Ref, 'pmax') * 0.99
        Tsat = np.linspace(max(Tmin, CP.PropsSI(Ref, 'Ttriple') + 0.1),
                           Tcrit - 0.01, 200)
        psat = CP.PropsSI('P', 'T', Tsat, 'Q', 0, Ref)
        hLs  = CP.PropsSI('Hmass', 'T', Tsat, 'Q', 0, Ref)
        hVs  = CP.PropsSI('Hmass', 'T', Tsat, 'Q', 1, Ref)
        random.seed(0)
        H, P, eT, eB, eS = [], [], [], [], []
        for _ in range(N):
            T = random.uniform(Tmin, Tmax)
            p = 10**random.uniform(np.log10(pmin), np.log10(pmax))
            try:
                EOS.update(CoolProp.PT_INPUTS, p, T)
                h = EOS.hmass(); rE = EOS.rhomolar()
                TTSE.update(CoolProp.HmassP_INPUTS, h, p)
                BICU.update(CoolProp.HmassP_INPUTS, h, p)
                SBTL.update(CoolProp.HmassP_INPUTS, h, p)
                H.append(h); P.append(p)
                eT.append(max(abs(TTSE.rhomolar()/rE - 1)*100, 1e-12))
                eB.append(max(abs(BICU.rhomolar()/rE - 1)*100, 1e-12))
                eS.append(max(abs(SBTL.rhomolar()/rE - 1)*100, 1e-12))
            except ValueError:
                pass
        for ax, e, lab in zip(axes[ir], (eT, eB, eS),
                              ('TTSE', 'Bicubic', 'SBTL')):
            ax.scatter(np.array(H)/1e3, P, c=e, s=3, cmap=cmap,
                       norm=cNorm, edgecolors='none')
            ax.plot(hLs/1e3, psat, 'k-', lw=1.0)
            ax.plot(hVs/1e3, psat, 'k-', lw=1.0)
            ax.set_yscale('log')
            ax.set_title(f'{Ref}: {lab}', fontsize=9)
            ax.set_xlabel('h [kJ/kg]', fontsize=8)
            ax.tick_params(labelsize=7)
            if ax is axes[ir][0]:
                ax.set_ylabel('p [Pa]', fontsize=8)
    sm = plt.cm.ScalarMappable(norm=cNorm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.92, 0.04, 0.020, 0.93])
    fig.colorbar(sm, cax=cax,
                 label=r'$|\rho_{\rm table}/\rho_{\rm EOS}-1|\times 100\,[\%]$')

``PT_INPUTS`` error map
~~~~~~~~~~~~~~~~~~~~~~~

SBTL's ``NormalizedPTTable`` places the saturation curve on a coordinate
axis (``xnorm = 1`` for LIQUID, ``xnorm = 0`` for VAPOR) so cells never
straddle the dome.  The BICUBIC panels show distinctive red streaks along
the saturation curve from cell-straddling; SBTL eliminates them.  In the
bulk single-phase region both backends are below
:math:`10^{-6}`.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    FLUIDS = ['Argon', 'CarbonDioxide', 'Water', 'R245fa', 'D6']
    N = 6000
    fig, axes = plt.subplots(len(FLUIDS), 2, figsize=(8, 14))
    fig.subplots_adjust(left=0.11, right=0.88, top=0.97, bottom=0.04,
                        wspace=0.30, hspace=0.50)
    cNorm = colors.LogNorm(vmin=1e-10, vmax=10)
    cmap = plt.get_cmap('jet')
    for ir, Ref in enumerate(FLUIDS):
        EOS  = CoolProp.AbstractState('HEOS',         Ref)
        BICU = CoolProp.AbstractState('BICUBIC&HEOS', Ref)
        SBTL = CoolProp.AbstractState('SBTL&HEOS',    Ref)
        Tmin  = CP.PropsSI(Ref, 'Tmin') + 0.5
        Tmax  = CP.PropsSI(Ref, 'Tmax') - 0.5
        Tcrit = CP.PropsSI(Ref, 'Tcrit')
        pmin  = CP.PropsSI(Ref, 'ptriple') * 1.05
        pmax  = CP.PropsSI(Ref, 'pmax') * 0.99
        Tsat = np.linspace(max(Tmin, CP.PropsSI(Ref, 'Ttriple') + 0.1),
                           Tcrit - 0.01, 200)
        psat = CP.PropsSI('P', 'T', Tsat, 'Q', 0, Ref)
        random.seed(0)
        TT, PP, eB, eS = [], [], [], []
        for _ in range(N):
            T = random.uniform(Tmin, Tmax)
            p = 10**random.uniform(np.log10(pmin), np.log10(pmax))
            try:
                EOS.update(CoolProp.PT_INPUTS, p, T); rE = EOS.rhomolar()
                BICU.update(CoolProp.PT_INPUTS, p, T)
                SBTL.update(CoolProp.PT_INPUTS, p, T)
                TT.append(T); PP.append(p)
                eB.append(max(abs(BICU.rhomolar()/rE - 1)*100, 1e-12))
                eS.append(max(abs(SBTL.rhomolar()/rE - 1)*100, 1e-12))
            except ValueError:
                pass
        for ax, e, lab in zip(axes[ir], (eB, eS), ('Bicubic', 'SBTL')):
            ax.scatter(TT, PP, c=e, s=3, cmap=cmap, norm=cNorm,
                       edgecolors='none')
            ax.plot(Tsat, psat, 'k-', lw=1.0)
            ax.set_yscale('log')
            ax.set_title(f'{Ref}: {lab}', fontsize=9)
            ax.set_xlabel('T [K]', fontsize=8)
            ax.tick_params(labelsize=7)
            if ax is axes[ir][0]:
                ax.set_ylabel('p [Pa]', fontsize=8)
    sm = plt.cm.ScalarMappable(norm=cNorm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.90, 0.04, 0.022, 0.93])
    fig.colorbar(sm, cax=cax,
                 label=r'$|\rho_{\rm table}/\rho_{\rm EOS}-1|\times 100\,[\%]$')

Near-critical PH behaviour
~~~~~~~~~~~~~~~~~~~~~~~~~~

Along the saturation curve (:math:`\eta = 0` boundary) the vapour-side
density behaves as

.. math::

    \rho_{\rm sat,V}(p) \approx \rho_{\rm crit} - C\,(p_{\rm crit} - p)^{\beta},
    \qquad \beta = 1/2

with :math:`\partial\rho_{\rm sat}/\partial p \to \infty` as
:math:`p \to p_{\rm crit}`.  The exponent :math:`\beta = 1/2` is the
classical mean-field value (CoolProp's pure-fluid Helmholtz EOS is
analytic in :math:`(\tau, \delta)` and therefore reproduces mean-field
scaling near critical, not the experimental 3D-Ising
:math:`\beta \approx 0.326`).  A cubic Hermite cannot reproduce a
square-root-style cusp; its residual at the cell midpoint is bounded by
the cell's :math:`\log p` span times derivatives of
:math:`\rho_{\rm sat}` that diverge near critical.

This is **not** a coordinate-normalization issue — the H-superancillary
delivers :math:`h_{\rm sat}(p)` to near-machine precision, so the
:math:`\eta = (h - h_{\rm sat}(p))/(h_{\rm hi}(p) - h_{\rm sat}(p))`
mapping at lookup time is exact.  The error is in the *property surface*
being interpolated, controlled by how short the cell's :math:`\log p`
span is at the cusp.

The subcritical PH grid uses a two-zone log-uniform :math:`p`-axis
layout that concentrates ~60 % of rows into
:math:`[0.5\,p_{\rm crit},\,p_{\rm crit}]`, shrinking the cell
:math:`\log p` span at the cusp by ~6× relative to a single
log-uniform grid.  Random-PH error at :math:`p = 0.93\,p_{\rm crit}`
on R245fa drops from ~0.59 % (single log-uniform grid) to ~10⁻⁵ %
(cusp-concentrated grid) at the dome, with similar improvements
elsewhere in :math:`(0.5,\,1.0)\,p_{\rm crit}`.  Residual error
remains in a narrow band right at the critical point where the
HEOS-fallback box takes over.

The supercritical region (``SUPER`` table) is unaffected — no dome, no
cusp — and stays at the bulk single-phase accuracy.  PT lookups avoid
the issue too because the PT normalization makes :math:`T` (not the
density) the dome coordinate, and :math:`T_{\rm sat}(p)` is smooth where
:math:`\rho_{\rm sat}(p)` is not.  Mitigation paths exist (finer cells
in the top 15 % of :math:`p`, or a non-polynomial correction at
:math:`\eta = 0`) but are out of scope for this PR.

Speed comparison
----------------

The benchmark measures wall time per ``AbstractState::update`` call from
Python for each backend, on the same set of randomly drawn single-phase
points outside the critical-fallback box.  All four backends see an
identical query sequence for each fluid; queries that any backend cannot
evaluate are screened out beforehand so the timing loops never throw.
N=20,000 queries per fluid, median of five replicates.

Numbers include the Cython wrapper overhead (~0.3 µs per call); the
*relative* comparison between backends is unaffected since all four pay
the same overhead.  Hardware: Apple M-series.  HEOS uses analytic
:math:`\alpha^r(\tau,\delta)` derivatives; ``HmassP_INPUTS`` requires a
Newton inversion which is responsible for the 10-50× HmassP/PT cost
ratio for HEOS.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import numpy as np
    import random
    import statistics
    import time

    FLUIDS = ['Argon', 'CarbonDioxide', 'Water', 'R245fa', 'D6']
    BACKENDS = [('HEOS', 'HEOS'),
                ('TTSE&HEOS', 'TTSE'),
                ('BICUBIC&HEOS', 'BICUBIC'),
                ('SBTL&HEOS', 'SBTL')]
    N = 20000
    N_REPS = 5
    results = {}
    for Ref in FLUIDS:
        EOS = CoolProp.AbstractState('HEOS', Ref)
        Tmin = CP.PropsSI(Ref, 'Tmin') + 1.0
        Tmax = CP.PropsSI(Ref, 'Tmax') - 1.0
        pmin = CP.PropsSI(Ref, 'ptriple') * 1.1
        pmax = CP.PropsSI(Ref, 'pmax') * 0.95
        pcrit = CP.PropsSI(Ref, 'pcrit')
        Tcrit = CP.PropsSI(Ref, 'Tcrit')
        random.seed(42)
        arrT, arrP, arrH = [], [], []
        while len(arrT) < N:
            T = random.uniform(Tmin, Tmax)
            p = 10**random.uniform(np.log10(pmin), np.log10(pmax))
            if abs(p / pcrit - 1) < 0.06 and abs(T / Tcrit - 1) < 0.06:
                continue
            try:
                EOS.update(CoolProp.PT_INPUTS, p, T)
                arrT.append(T); arrP.append(p); arrH.append(EOS.hmass())
            except ValueError:
                continue
        arrT = np.array(arrT); arrP = np.array(arrP); arrH = np.array(arrH)
        screens = [CoolProp.AbstractState(b, Ref) for b, _ in BACKENDS]
        mask = np.ones(N, dtype=bool)
        for i in range(N):
            for AS_s in screens:
                try:
                    AS_s.update(CoolProp.PT_INPUTS, arrP[i], arrT[i])
                    AS_s.update(CoolProp.HmassP_INPUTS, arrH[i], arrP[i])
                except ValueError:
                    mask[i] = False; break
        arrT = arrT[mask]; arrP = arrP[mask]; arrH = arrH[mask]
        Ng = len(arrT)
        for backend, name in BACKENDS:
            AS = CoolProp.AbstractState(backend, Ref)
            AS.update(CoolProp.PT_INPUTS, arrP[0], arrT[0])
            for input_label, pair, X, Y in [
                ('PT', CoolProp.PT_INPUTS, arrP, arrT),
                ('HmassP', CoolProp.HmassP_INPUTS, arrH, arrP),
            ]:
                trials = []
                for _ in range(N_REPS):
                    t0 = time.perf_counter()
                    for i in range(Ng):
                        AS.update(pair, X[i], Y[i])
                    trials.append((time.perf_counter() - t0) / Ng)
                results[(Ref, name, input_label)] = statistics.median(trials)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(FLUIDS))
    width = 0.20
    colors_b = {'HEOS': '#555', 'TTSE': '#1f77b4',
                'BICUBIC': '#ff7f0e', 'SBTL': '#2ca02c'}
    for ax, ipair, ymax in zip(axes, ['PT', 'HmassP'], [12, 100]):
        for k, (_, name) in enumerate(BACKENDS):
            ys = [results[(f, name, ipair)] * 1e6 for f in FLUIDS]
            ax.bar(x + (k - 1.5) * width, ys, width, color=colors_b[name],
                   label=name, edgecolor='black', linewidth=0.4)
            for xi, yi in zip(x + (k - 1.5) * width, ys):
                ax.text(xi, yi + ymax * 0.012, f'{yi:.2f}',
                        ha='center', va='bottom', fontsize=7.5)
        ax.set_xticks(x); ax.set_xticklabels(FLUIDS, rotation=15)
        ax.set_ylabel('Time per call [µs]')
        ax.set_title(f'{ipair}_INPUTS')
        ax.legend(loc='upper right' if ipair == 'PT' else 'upper center',
                  fontsize=8, ncols=4 if ipair == 'HmassP' else 1)
        ax.set_ylim(0, ymax)
        ax.grid(axis='y', alpha=0.3)
    fig.suptitle(
        f'Per-call cost (median of {N_REPS} runs of N={N//1000}k '
        f'single-phase queries; critical-fallback box excluded)',
        fontsize=10)
    fig.subplots_adjust(top=0.86, bottom=0.13, left=0.07, right=0.97,
                        wspace=0.22)

Takeaways:

* **TTSE and BICUBIC** are the fastest tabular backends at ~0.4–0.7 µs/call
  for both input pairs.  Both ultimately evaluate a 16-coefficient
  polynomial; their cost is dominated by Python-wrapper overhead.
* **SBTL** is ~0.5–1.5 µs/call.  For ``PT_INPUTS`` it is 2–3× slower than
  BICUBIC: the per-query saturation lookup for region routing plus the
  Chebyshev evaluation for the :math:`\eta`-coordinate normalization are
  the unavoidable cost of putting the dome on a coordinate axis.  For
  ``HmassP_INPUTS`` SBTL is *competitive with or faster than* BICUBIC —
  the coordinate-aligned PH table inverts the lookup directly, while
  BICUBIC pays for the dome-straddling cell-bump logic.
* **HEOS direct** is the slow path: ~2–9 µs for ``PT_INPUTS``, ~17–93 µs
  for ``HmassP_INPUTS`` (the latter is iterated; the former is closed-form).
* The headline speedup vs HEOS is **3–10×** on ``PT_INPUTS`` depending on
  fluid (HEOS's own PT cost varies 4× across the five fluids), and
  **20–130×** on ``HmassP_INPUTS`` since the tabular backends skip HEOS's
  Newton inversion.

For queries inside the critical-region HEOS-fallback box, SBTL bypasses
the spline and calls HEOS directly (~10 µs), trading speed for the
machine-precision answer that no smooth polynomial can give at the cusp.

More information
----------------

The coordinate-aligned table design is described in:

    International Association for the Properties of Water and Steam,
    *Guideline on the Fast Calculation of Steam and Water Properties with
    the Spline-Based Table Look-Up Method (SBTL)*, IAPWS G13-15 (2015).
    `<https://iapws.org/documents/release/SBTL>`_

Per-region implementation notes live alongside the SBTL backend source
in the ``dev/`` directory of the CoolProp source tree:

* ``dev/sbtl_normalized_ph_design.md`` — PH coordinate-aligned table.
* ``dev/sbtl_normalized_pt_design.md`` — PT coordinate-aligned table.
* ``dev/sbtl_pt_outstanding_work.md`` — current status and follow-up work.
* ``dev/sbtl_du_kunick_redesign_plan.md`` — the queued ``(D, U)``
  follow-up PR.

The tabular data shares the on-disk format of the BICUBIC and TTSE
backends; see :ref:`tabular_interpolation` for details on the directory
structure and serialisation format.
