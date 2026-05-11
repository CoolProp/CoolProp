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

The grid storage and on-disk caching reuse the same machinery as BICUBIC and
TTSE (``HOME/.CoolProp/Tables``), and SBTL is subject to the same
directory-size warnings.

Backbone interpolation
----------------------

For each cell :math:`[x_i, x_{i+1}] \times [y_j, y_{j+1}]`, every property
:math:`z` is represented as a 16-coefficient bi-cubic polynomial in the
cell-local coordinates:

.. math::

    z(\hat{x}, \hat{y}) = \sum_{m=0}^{3} \sum_{n=0}^{3} a_{mn}\, \hat{x}^m \hat{y}^n

The coefficients come from one of two families:

* **Coordinate-aligned PH tables** (``NormalizedPHTable``, three flavours
  LIQUID / VAPOR / SUPER) — coefficients computed from a 2-D natural cubic
  B-spline (Kunick C\ :sup:`2`) over a uniform :math:`(\hat x, \log p)` grid
  whose left/right column boundaries are the ``h_lo(P)`` and ``h_hi(P)``
  isotherms / saturation isobars.
* **Coordinate-aligned PT tables** (``NormalizedPTTable``, same three flavours)
  — same B-spline backbone with ``T_lo(P)`` and ``T_hi(P)`` replacing
  ``h_lo`` / ``h_hi``.

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

.. code-block:: c++

   auto AS = AbstractState::factory("SBTL&HEOS", "CO2");
   const double P = 15.8e5;                       // T_sat ≈ 246.2 K
   const double T_sat = 246.207;
   // Default: T = T_sat routes to LIQUID by convention
   AS->update(PT_INPUTS, P, T_sat);               // rho() ≈ 1063 (sat-L)
   AS->specify_phase(iphase_gas);
   AS->update(PT_INPUTS, P, T_sat);               // rho() ≈ 41   (sat-V)
   AS->unspecify_phase();

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

    In [0]: import CoolProp

    In [1]: HEOS = CoolProp.AbstractState("HEOS",         "R245fa")

    In [2]: TTSE = CoolProp.AbstractState("TTSE&HEOS",    "R245fa")

    In [3]: BICU = CoolProp.AbstractState("BICUBIC&HEOS", "R245fa")

    In [4]: SBTL = CoolProp.AbstractState("SBTL&HEOS",    "R245fa")

    In [5]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300); TTSE.update(CoolProp.PT_INPUTS, 101325, 300); BICU.update(CoolProp.PT_INPUTS, 101325, 300); SBTL.update(CoolProp.PT_INPUTS, 101325, 300)

    In [6]: print(HEOS.rhomolar(), TTSE.rhomolar(), BICU.rhomolar(), SBTL.rhomolar())

The figure below shows density errors for ``HmassP_INPUTS`` across the
single-phase region of R245fa.  SBTL routes through the coordinate-aligned
``NormalizedPHTable``; BICUBIC uses Hermite bi-cubic on a ``logph`` grid;
TTSE uses bilinear Taylor expansion.  SBTL and BICUBIC both substantially
outperform TTSE; SBTL additionally avoids the cell-straddling errors that
BICUBIC sees just above and below the saturation curve.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.ticker
    import numpy as np
    import random

    fig = plt.figure(figsize=(14, 5))
    ax1 = fig.add_axes((0.05, 0.12, 0.24, 0.80))
    ax2 = fig.add_axes((0.35, 0.12, 0.24, 0.80))
    ax3 = fig.add_axes((0.65, 0.12, 0.24, 0.80))

    Ref = 'R245fa'

    SBTL    = CoolProp.AbstractState('SBTL&HEOS',    Ref)
    BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS', Ref)
    TTSE    = CoolProp.AbstractState('TTSE&HEOS',    Ref)
    EOS     = CoolProp.AbstractState('HEOS',         Ref)

    T  = np.linspace(CP.PropsSI(Ref, 'Tmin') + 0.1, CP.PropsSI(Ref, 'Tcrit') - 0.01, 300)
    pV = CP.PropsSI('P',     'T', T, 'Q', 1, Ref)
    hL = CP.PropsSI('Hmass', 'T', T, 'Q', 0, Ref)
    hV = CP.PropsSI('Hmass', 'T', T, 'Q', 1, Ref)

    HHH, PPP, E_TTSE, E_BICU, E_SBTL = [], [], [], [], []

    cNorm = colors.LogNorm(vmin=1e-12, vmax=10)
    random.seed(0)

    for _ in range(40000):
        h = random.uniform(150000, 590000)
        p = 10**random.uniform(np.log10(100000), np.log10(7000000))
        try:
            EOS.update(CoolProp.HmassP_INPUTS, h, p)
            rhoEOS = EOS.rhomolar()

            TTSE.update(CoolProp.HmassP_INPUTS, h, p)
            BICUBIC.update(CoolProp.HmassP_INPUTS, h, p)
            SBTL.update(CoolProp.HmassP_INPUTS, h, p)

            HHH.append(h); PPP.append(p)
            E_TTSE.append(abs(TTSE.rhomolar()    / rhoEOS - 1) * 100)
            E_BICU.append(abs(BICUBIC.rhomolar() / rhoEOS - 1) * 100)
            E_SBTL.append(abs(SBTL.rhomolar()    / rhoEOS - 1) * 100)
        except ValueError:
            pass

    kw = dict(s=8, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
    SC = ax1.scatter(HHH, PPP, c=E_TTSE, **kw)
    ax2.scatter(HHH, PPP, c=E_BICU, **kw)
    ax3.scatter(HHH, PPP, c=E_SBTL, **kw)

    titles = ['Error from TTSE', 'Error from Bicubic', 'Error from SBTL']
    for ax, title in zip([ax1, ax2, ax3], titles):
        ax.set_title(title)
        ax.set_xlim(250000, 550000)
        ax.set_ylim(100000, 7000000)
        ax.set_yscale('log')
        ax.plot(hL, pV, 'k', lw=4)
        ax.plot(hV, pV, 'k', lw=4)
        ticks_h = [150000, 250000, 350000, 450000, 550000]
        ticks_p = [100000, 200000, 400000, 600000, 800000, 1000000, 2000000, 4000000, 6000000]
        ax.set_xticks(ticks_h); ax.set_xticklabels([str(t//1000) for t in ticks_h])
        ax.set_yticks(ticks_p); ax.set_yticklabels([str(t//1000) for t in ticks_p])
        ax.tick_params(axis='y', which='minor', left='off')
        ax.set_xlabel('Enthalpy [kJ/kg]')
        ax.set_ylabel('Pressure [kPa]')

    cbar_ax = fig.add_axes([0.91, 0.15, 0.025, 0.70])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'$(\rho/\rho_{\rm EOS}-1)\times 100$ [%]')

The next figure shows density errors for ``PT_INPUTS`` over the same
window.  SBTL's coordinate-aligned ``NormalizedPTTable`` matches BICUBIC's
Hermite bi-cubic in well-conditioned regions and additionally avoids
the cell-straddling errors that BICUBIC sees just above and below the
saturation curve.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_axes((0.07, 0.13, 0.36, 0.78))
    ax2 = fig.add_axes((0.50, 0.13, 0.36, 0.78))

    Ref = 'R245fa'

    SBTL    = CoolProp.AbstractState('SBTL&HEOS',    Ref)
    BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS', Ref)
    EOS     = CoolProp.AbstractState('HEOS',         Ref)

    Ttpl  = EOS.Ttriple()
    Tcrit = EOS.T_critical()
    Tmax  = EOS.Tmax()
    ptpl  = CP.PropsSI(Ref, 'ptriple')
    pmax  = CP.PropsSI(Ref, 'pmax')

    # Saturation dome in (T, P) coordinates
    T_sat = np.linspace(Ttpl + 0.1, Tcrit - 0.01, 300)
    p_sat = CP.PropsSI('P', 'T', T_sat, 'Q', 0, Ref)

    random.seed(0)
    TTT, PPP, E_BICU, E_SBTL = [], [], [], []
    for _ in range(20000):
        T = random.uniform(Ttpl + 1, Tmax)
        P = 10**random.uniform(np.log10(ptpl * 1.05), np.log10(pmax))
        try:
            EOS.update(CoolProp.PT_INPUTS, P, T)
            if EOS.phase() == CP.iphase_twophase:
                continue
            rho_eos = EOS.rhomolar()
            BICUBIC.update(CoolProp.PT_INPUTS, P, T)
            SBTL.update(CoolProp.PT_INPUTS, P, T)
            TTT.append(T); PPP.append(P)
            E_BICU.append(abs(BICUBIC.rhomolar() / rho_eos - 1) * 100)
            E_SBTL.append(abs(SBTL.rhomolar()    / rho_eos - 1) * 100)
        except (ValueError, RuntimeError):
            pass

    cNorm = colors.LogNorm(vmin=1e-12, vmax=10)
    kw = dict(s=4, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
    SC = ax1.scatter(TTT, PPP, c=E_BICU, **kw)
    ax2.scatter(TTT, PPP, c=E_SBTL, **kw)

    for ax, title in zip([ax1, ax2], ['Error from BICUBIC PT', 'Error from SBTL PT']):
        ax.set_title(title)
        ax.set_yscale('log')
        ax.plot(T_sat, p_sat, 'k', lw=3)
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Pressure [Pa]')

    cbar_ax = fig.add_axes([0.89, 0.15, 0.025, 0.74])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'$(\rho/\rho_{\rm EOS}-1)\times 100$ [%]')

Speed comparison
----------------

For ``HmassP_INPUTS`` and ``PT_INPUTS``, SBTL and BICUBIC have similar
per-call cost since both ultimately evaluate a 16-coefficient bi-cubic
polynomial.  The coordinate-aligned tables add a small constant overhead
(one ``T_sat(P)`` lookup for the subcritical routing decision, cached
single-deep for repeat queries at the same P).  In a typical pipeline
workload the cache hits and SBTL averages ≈ 1.3 µs per ``PT_INPUTS``
call vs ≈ 10.7 µs for HEOS direct (≈ 8× speedup).

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
