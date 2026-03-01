.. _sbtl_interpolation:

*******************************
SBTL Bi-Quadratic Interpolation
*******************************

The SBTL (Spline-Based Table Look-Up) method uses bi-quadratic (degree-2) Lagrange polynomials to interpolate thermodynamic properties from precomputed tables.  It uses the same pressure–enthalpy and pressure–temperature grids as the :ref:`BICUBIC and TTSE backends <tabular_interpolation>` (stored in ``HOME/.CoolProp/Tables``), and is subject to the same directory-size warnings.

For pure fluids, SBTL also builds dedicated **DU-space tables** and **DT-space tables** that make ``DmolarUmolar_INPUTS`` and ``DmolarT_INPUTS`` direct forward evaluations with no Newton iteration — the fastest supported input pairs.

SBTL Interpolation
------------------

For each cell :math:`[x_i, x_{i+1}] \times [y_j, y_{j+1}]`, the property :math:`z` is represented as a bi-quadratic polynomial in normalized coordinates :math:`\hat{x}, \hat{y} \in [0,1]`:

.. math::

    z(\hat{x}, \hat{y}) = \sum_{m=0}^{2} \sum_{n=0}^{2} a_{mn}\, \hat{x}^m \hat{y}^n

.. math::

    \hat{x} = \frac{x - x_i}{x_{i+1} - x_i}, \qquad \hat{y} = \frac{y - y_j}{y_{j+1} - y_j}

The 9 coefficients :math:`a_{mn}` are obtained by fitting a degree-2 Lagrange polynomial through the 3×3 stencil of node values centred at cell :math:`(i,j)`.  The polynomial agrees exactly with the tabulated values at all nine stencil nodes.

Compared to bicubic interpolation — which fits a degree-3 polynomial by matching both values and derivatives at the four cell corners — SBTL requires only node values.  Table construction is therefore faster and inversion (finding :math:`\hat{x}` given a target :math:`z` and known :math:`\hat{y}`) reduces to a 1-D quadratic equation solved analytically.  The trade-off is that the lower polynomial degree makes SBTL somewhat less accurate than BICUBIC for pressure–enthalpy and pressure–temperature lookups; both backends are substantially more accurate than TTSE.  The main advantage of SBTL over BICUBIC is the dedicated DU-space tables described below.

``DmolarUmolar_INPUTS`` — Direct Flash
---------------------------------------

For pure fluids, SBTL builds additional **DU-space tables** storing :math:`T(D,U)` and :math:`p(D,U)` as precomputed bi-quadratic polynomials in molar density *D* and molar internal energy *U*.  When ``DmolarUmolar_INPUTS`` is used the property is obtained by a single polynomial evaluation:

1. **O(1) phase detection**: *D* is compared against cached triple-point saturation densities.  States outside this density range are guaranteed single-phase.
2. **Single-phase**: The appropriate liquid (:math:`D > D_{\rm crit}`) or gas (:math:`D \le D_{\rm crit}`) DU table is evaluated directly.
3. **Two-phase check**: If the density falls in the potentially two-phase range, a 1-D Newton iteration on a precomputed saturation cache (O(1) lookup per evaluation via log-indexed interpolation) determines the saturation pressure, and the quality follows from the lever rule.

All other input pairs (e.g. ``PT_INPUTS``, ``HmassP_INPUTS``) use the standard pressure–temperature and pressure–enthalpy tables, exactly as BICUBIC and TTSE do.

``DmolarT_INPUTS`` — Direct Flash
-----------------------------------

For pure fluids, SBTL builds additional **DT-space tables** storing :math:`p(D,T)`, :math:`h(D,T)`, :math:`s(D,T)` and :math:`u(D,T)` as precomputed bi-quadratic polynomials in molar density *D* and temperature *T*.  When ``DmolarT_INPUTS`` is used properties are obtained by a single polynomial evaluation:

1. **Phase detection**: at the given temperature *T*, the density *D* is compared against saturated liquid and vapour densities from the cached saturation table.  States with :math:`D > D_{\rm sat,L}(T)` are liquid; :math:`D < D_{\rm sat,V}(T)` are gas.
2. **Single-phase**: the liquid (:math:`D > D_{\rm crit}`) or gas (:math:`D \le D_{\rm crit}`) DT table is evaluated directly.
3. **Two-phase**: if :math:`D_{\rm sat,V}(T) \le D \le D_{\rm sat,L}(T)`, the quality follows from the lever rule on the saturation table.

Because *D* and *T* are the natural coordinates of both tables, no Newton iteration or polynomial inversion is required even near the saturation boundary — in contrast to the pressure–temperature table path, which must invert :math:`\rho(T,p)` for *p*, an ill-conditioned step for nearly-incompressible liquids.

Accuracy comparison
-------------------

Here is a simple comparison of accuracy: density is obtained for R245fa using the equation of state, TTSE, Bicubic, and SBTL interpolation.

.. ipython::

    In [0]: import CoolProp

    In [1]: HEOS = CoolProp.AbstractState("HEOS",         "R245fa")

    In [2]: TTSE = CoolProp.AbstractState("TTSE&HEOS",    "R245fa")

    In [3]: BICU = CoolProp.AbstractState("BICUBIC&HEOS", "R245fa")

    In [4]: SBTL = CoolProp.AbstractState("SBTL&HEOS",    "R245fa")

    In [5]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300); TTSE.update(CoolProp.PT_INPUTS, 101325, 300); BICU.update(CoolProp.PT_INPUTS, 101325, 300); SBTL.update(CoolProp.PT_INPUTS, 101325, 300)

    In [6]: print(HEOS.rhomolar(), TTSE.rhomolar(), BICU.rhomolar(), SBTL.rhomolar())

The figure below shows density errors for ``HmassP_INPUTS`` across the single-phase region of R245fa.  SBTL is less accurate than BICUBIC here because the degree-2 polynomial captures less curvature than degree-3; both are substantially better than TTSE.  The primary use case for SBTL — ``DmolarUmolar_INPUTS`` — is shown separately below.

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
    ax1 = fig.add_axes((0.06, 0.12, 0.26, 0.80))
    ax2 = fig.add_axes((0.38, 0.12, 0.26, 0.80))
    ax3 = fig.add_axes((0.70, 0.12, 0.26, 0.80))

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

    cbar_ax = fig.add_axes([0.97, 0.15, 0.02, 0.70])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'$(\rho/\rho_{\rm EOS}-1)\times 100$ [%]')

The ``DmolarUmolar_INPUTS`` fast path can also be checked by round-tripping through HEOS:

.. ipython::

    In [7]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300)

    In [8]: D = HEOS.rhomolar(); U = HEOS.umolar()

    In [9]: SBTL.update(CoolProp.DmolarUmolar_INPUTS, D, U)

    In [10]: print('T:', SBTL.T(), '  p:', SBTL.p())

The figure below gives a comprehensive view of SBTL accuracy for ``DmolarUmolar_INPUTS`` across the single-phase liquid and gas regions of R245fa.  The saturation dome (black curve) separates liquid (upper-left, high density) from gas (lower-right, low density).

.. note::

    Relative pressure error is large for compressed liquid states at low pressure (upper-left corner), where the equation of state is nearly incompressible: a small change in density corresponds to a large pressure change, so any fixed-resolution table will show elevated :math:`\Delta p / p`.  Temperature accuracy is good throughout.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    Ref = 'R245fa'

    SBTL = CoolProp.AbstractState('SBTL&HEOS', Ref)
    EOS  = CoolProp.AbstractState('HEOS',      Ref)

    Ttpl  = EOS.Ttriple()
    Tcrit = EOS.T_critical()
    Tmax  = EOS.Tmax()
    ptpl  = CP.PropsSI(Ref, 'ptriple')
    pmax  = CP.PropsSI(Ref, 'pmax')

    # Saturation dome in (u_mass, rho_mass) space
    T_sat  = np.linspace(Ttpl + 0.1, Tcrit - 0.01, 300)
    u_satL = np.zeros_like(T_sat); rho_satL = np.zeros_like(T_sat)
    u_satV = np.zeros_like(T_sat); rho_satV = np.zeros_like(T_sat)
    for k, T in enumerate(T_sat):
        EOS.update(CP.QT_INPUTS, 0, T)
        u_satL[k] = EOS.umass() / 1000; rho_satL[k] = EOS.rhomass()
        EOS.update(CP.QT_INPUTS, 1, T)
        u_satV[k] = EOS.umass() / 1000; rho_satV[k] = EOS.rhomass()

    # Random single-phase test points
    random.seed(1)
    UUU, RHO, errT, errp = [], [], [], []
    attempts = 0
    while len(UUU) < 20000 and attempts < 300000:
        attempts += 1
        T = random.uniform(Ttpl + 1, Tmax)
        P = 10**random.uniform(np.log10(ptpl * 1.01), np.log10(pmax))
        try:
            EOS.update(CP.PT_INPUTS, P, T)
            if EOS.phase() == CP.iphase_twophase:
                continue
            D = EOS.rhomolar(); U = EOS.umolar()
            SBTL.update(CP.DmolarUmolar_INPUTS, D, U)
            errT.append(abs(SBTL.T() / T  - 1) * 100)
            errp.append(abs(SBTL.p() / P  - 1) * 100)
            UUU.append(EOS.umass() / 1000)
            RHO.append(EOS.rhomass())
        except:
            pass

    UUU = np.array(UUU); RHO = np.array(RHO)
    errT = np.array(errT); errp = np.array(errp)

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_axes((0.09, 0.12, 0.38, 0.80))
    ax2 = fig.add_axes((0.57, 0.12, 0.38, 0.80))

    cNorm = colors.LogNorm(vmin=1e-8, vmax=1e-1)
    kw = dict(s=4, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
    SC = ax1.scatter(UUU, RHO, c=errT, **kw)
    ax2.scatter(UUU, RHO, c=errp, **kw)

    for ax, title in zip([ax1, ax2],
                         [r'Error in $T(D,U)$', r'Error in $p(D,U)$']):
        ax.plot(u_satL, rho_satL, 'k', lw=3)
        ax.plot(u_satV, rho_satV, 'k', lw=3)
        ax.set_yscale('log')
        ax.set_xlabel('Specific internal energy [kJ/kg]')
        ax.set_ylabel(r'Density [kg/m$^3$]')
        ax.set_title(title)

    cbar_ax = fig.add_axes([0.97, 0.15, 0.02, 0.70])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'Relative error $\times\,100$ [%]')

The plot below shows the same temperature error data for ``DmolarUmolar_INPUTS`` in molar :math:`(D,U)` coordinates — molar density on the *y*-axis versus molar internal energy on the *x*-axis.  This is the natural coordinate system of the DU tables and makes the saturation dome appear in the orientation directly relevant to table-lookup phase detection.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    Ref = 'R245fa'

    SBTL = CoolProp.AbstractState('SBTL&HEOS', Ref)
    EOS  = CoolProp.AbstractState('HEOS',      Ref)

    Ttpl  = EOS.Ttriple()
    Tcrit = EOS.T_critical()
    Tmax  = EOS.Tmax()
    ptpl  = CP.PropsSI(Ref, 'ptriple')
    pmax  = CP.PropsSI(Ref, 'pmax')

    # Saturation dome in (U_molar [kJ/mol], D_molar [mol/m³]) space
    T_sat  = np.linspace(Ttpl + 0.1, Tcrit - 0.01, 300)
    u_satL = np.zeros_like(T_sat); D_satL = np.zeros_like(T_sat)
    u_satV = np.zeros_like(T_sat); D_satV = np.zeros_like(T_sat)
    for k, T in enumerate(T_sat):
        EOS.update(CP.QT_INPUTS, 0, T)
        u_satL[k] = EOS.umolar() / 1000; D_satL[k] = EOS.rhomolar()
        EOS.update(CP.QT_INPUTS, 1, T)
        u_satV[k] = EOS.umolar() / 1000; D_satV[k] = EOS.rhomolar()

    # Random single-phase test points (same seed as existing DU plot)
    random.seed(1)
    UUU_mol, DDD_mol, errT = [], [], []
    attempts = 0
    while len(UUU_mol) < 20000 and attempts < 300000:
        attempts += 1
        T = random.uniform(Ttpl + 1, Tmax)
        P = 10**random.uniform(np.log10(ptpl * 1.01), np.log10(pmax))
        try:
            EOS.update(CP.PT_INPUTS, P, T)
            if EOS.phase() == CP.iphase_twophase:
                continue
            D = EOS.rhomolar(); U = EOS.umolar()
            SBTL.update(CP.DmolarUmolar_INPUTS, D, U)
            errT.append(abs(SBTL.T() / T - 1) * 100)
            UUU_mol.append(U / 1000)
            DDD_mol.append(D)
        except:
            pass

    UUU_mol = np.array(UUU_mol); DDD_mol = np.array(DDD_mol)
    errT    = np.array(errT)

    fig = plt.figure(figsize=(6, 5))
    ax  = fig.add_axes((0.14, 0.12, 0.76, 0.80))

    cNorm = colors.LogNorm(vmin=1e-8, vmax=1e-1)
    kw = dict(s=4, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
    SC = ax.scatter(UUU_mol, DDD_mol, c=errT, **kw)
    ax.plot(u_satL, D_satL, 'k', lw=3)
    ax.plot(u_satV, D_satV, 'k', lw=3)
    ax.set_yscale('log')
    ax.set_xlabel(r'Molar internal energy [kJ/mol]')
    ax.set_ylabel(r'Molar density [mol/m$^3$]')
    ax.set_title(r'Error in $T(D,U)$ — molar $(D,\,U)$ coordinates')

    cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.70])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'Relative error $\times\,100$ [%]')

The figure below gives a comprehensive view of SBTL accuracy for ``DmolarT_INPUTS`` across the single-phase liquid and gas regions of R245fa, plotted directly in the :math:`(T,\rho)` input-coordinate plane.  The saturation dome (black curve) bounds the two-phase region; all test points are single-phase.

.. note::

    Relative pressure error is elevated in the compressed-liquid region (high density, moderate temperature) because the nearly-incompressible equation of state means a large pressure change corresponds to only a tiny density change.  Temperature and enthalpy accuracy are good throughout.

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import random

    Ref = 'R245fa'

    SBTL = CoolProp.AbstractState('SBTL&HEOS', Ref)
    EOS  = CoolProp.AbstractState('HEOS',      Ref)

    Ttpl  = EOS.Ttriple()
    Tmax  = EOS.Tmax()
    Tcrit = EOS.T_critical()
    ptpl  = CP.PropsSI(Ref, 'ptriple')
    pmax  = CP.PropsSI(Ref, 'pmax')

    # Saturation dome in (T [K], rho_mass [kg/m³]) space
    T_sat    = np.linspace(Ttpl + 0.1, Tcrit - 0.01, 300)
    rho_satL = np.zeros_like(T_sat)
    rho_satV = np.zeros_like(T_sat)
    for k, T in enumerate(T_sat):
        EOS.update(CP.QT_INPUTS, 0, T)
        rho_satL[k] = EOS.rhomass()
        EOS.update(CP.QT_INPUTS, 1, T)
        rho_satV[k] = EOS.rhomass()

    # Random single-phase test points sampled from (T, P) space
    random.seed(2)
    TTT, RHO, errP, errH = [], [], [], []
    attempts = 0
    while len(TTT) < 20000 and attempts < 300000:
        attempts += 1
        T = random.uniform(Ttpl + 1, Tmax)
        P = 10**random.uniform(np.log10(ptpl * 1.01), np.log10(pmax))
        try:
            EOS.update(CP.PT_INPUTS, P, T)
            if EOS.phase() == CP.iphase_twophase:
                continue
            D     = EOS.rhomolar()
            P_ref = EOS.p()
            H_ref = EOS.hmass()
            SBTL.update(CP.DmolarT_INPUTS, D, T)
            errP.append(abs(SBTL.p()     / P_ref - 1) * 100)
            errH.append(abs(SBTL.hmass() / H_ref - 1) * 100)
            TTT.append(T)
            RHO.append(EOS.rhomass())
        except:
            pass

    TTT  = np.array(TTT);  RHO  = np.array(RHO)
    errP = np.array(errP); errH = np.array(errH)

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_axes((0.09, 0.12, 0.38, 0.80))
    ax2 = fig.add_axes((0.57, 0.12, 0.38, 0.80))

    cNorm = colors.LogNorm(vmin=1e-8, vmax=1e-1)
    kw = dict(s=4, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
    SC = ax1.scatter(TTT, RHO, c=errP, **kw)
    ax2.scatter(TTT, RHO, c=errH, **kw)

    for ax, title in zip([ax1, ax2],
                         [r'Error in $p(D,T)$', r'Error in $H(D,T)$']):
        ax.plot(T_sat, rho_satL, 'k', lw=3)
        ax.plot(T_sat, rho_satV, 'k', lw=3)
        ax.set_yscale('log')
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel(r'Density [kg/m$^3$]')
        ax.set_title(title)

    cbar_ax = fig.add_axes([0.97, 0.15, 0.02, 0.70])
    CB = fig.colorbar(SC, cax=cbar_ax)
    CB.set_label(r'Relative error $\times\,100$ [%]')

Speed comparison
----------------

The primary motivation for SBTL is computational speed.  The table below gives Computing Time Ratios (CTR = t\ :sub:`reference EOS` / t\ :sub:`SBTL`) from the IAPWS G13-15 benchmark for water, where IAPWS-95 (REFPROP 10) is used as the reference and phase is pre-specified so that phase-detection overhead is excluded from both timings.

.. list-table::
   :header-rows: 1
   :widths: 40 20 20

   * - Input pair
     - CTR liquid
     - CTR gas
   * - ``DmolarUmolar_INPUTS``
     - 251×
     - 410×
   * - ``DmolarT_INPUTS``
     - comparable to DmolarUmolar
     - comparable to DmolarUmolar
   * - ``HmassP_INPUTS``
     - comparable to BICUBIC
     - comparable to BICUBIC

For ``DmolarUmolar_INPUTS`` and ``DmolarT_INPUTS``, SBTL uses dedicated DU-space and DT-space tables respectively; each lookup is a single bi-quadratic polynomial evaluation with no iteration.  For all other input pairs, SBTL and BICUBIC have similar performance because both evaluate a fixed-order polynomial on the same underlying pressure–enthalpy and pressure–temperature grids.

More Information
----------------

The SBTL method is described in the IAPWS guideline:

    Kruse, A. & Knobloch, K., *IAPWS Guideline on the SBTL Method for Fast Calculations of Water and Steam Properties*, IAPWS G13-15 (2015). `<http://www.iapws.org/relguide/SBTL-2015.pdf>`_

The tabular data used by SBTL is stored in the same format as the BICUBIC and TTSE backends; see :ref:`tabular_interpolation` for details on the directory structure and serialisation format.
