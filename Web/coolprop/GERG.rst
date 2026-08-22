.. _gerg_backend:

*********************************************
GERG-2004 and GERG-2008 Equations of State
*********************************************

.. contents:: :depth: 2

Introduction
============

CoolProp provides two *strict* implementations of the GERG wide-range equations
of state for natural gases and other mixtures: ``GERG2004``
:cite:`Kunz-BOOK-2007` and ``GERG2008`` :cite:`Kunz-JCED-2012`.  They are
separate backend families, not options on the default multi-fluid backend.

"Strict" means exactly what it says:

* only the components each model publishes are admitted,
* only the parameters published with each model are carried — the pure-fluid
  equations of state, the ideal-gas coefficients, the binary reducing
  parameters, and the departure functions all come from the model's own
  tables rather than from CoolProp's fluid and mixture libraries,
* the gas constant is GERG's :math:`R = 8.314472\ \mathrm{J\,mol^{-1}\,K^{-1}}`,
  not the CODATA value CoolProp normally uses, and
* anything the model does not cover raises an exception instead of quietly
  answering from somewhere else.

A number obtained from the ``GERG2008`` backend is a GERG-2008 number: no
parameter in it is borrowed from another correlation.  The refusals are
model-level rules enforced through the documented API, not a sandbox — see
*The strictness rules are model-level, not C++-level* below.

Usage
=====

Change the backend name, exactly as for the other backends.  In the
:ref:`high-level interface <high_level_api>`::

    PropsSI("Dmolar", "T", 300, "P", 1e6, "GERG2008::METHANE[0.9]&NITROGEN[0.1]")
    # -> 406.9418737 mol/m^3

and in the :ref:`low-level interface <low_level_api>`::

    AS = CP.AbstractState("GERG2008", "Methane&Nitrogen")
    AS.set_mole_fractions([0.9, 0.1])
    AS.update(CP.PT_INPUTS, 1e6, 300)
    AS.rhomolar()   # -> 406.9418737 mol/m^3

Pure fluids work the same way — ``GERG2008::Methane`` is the GERG-2008
pure-methane equation of state, which is *not* the Setzmann-Wagner reference
equation CoolProp's ``HEOS::Methane`` uses.

``GERG2004`` behaves identically apart from its smaller component set and the
handful of parameters GERG-2008 revised.

Components
==========

GERG-2004 defines 18 components; GERG-2008 keeps all 18 and adds three more,
for 21.  GERG-2008 also revises the pure-fluid equation of state for carbon
monoxide and isopentane, so those two give different numbers under the two
models.

Component names resolve through CoolProp's normal alias and CAS machinery, so
``CO2``, ``R744``, ``CarbonDioxide`` and ``124-38-9`` all reach the same GERG
component.  The canonical CoolProp names are listed below.

============================  =========  =========  ==================
CoolProp name                 GERG-2004  GERG-2008  CAS
============================  =========  =========  ==================
``Methane``                   yes        yes        74-82-8
``Nitrogen``                  yes        yes        7727-37-9
``CarbonDioxide``             yes        yes        124-38-9
``Ethane``                    yes        yes        74-84-0
``n-Propane``                 yes        yes        74-98-6
``n-Butane``                  yes        yes        106-97-8
``IsoButane``                 yes        yes        75-28-5
``n-Pentane``                 yes        yes        109-66-0
``Isopentane``                yes        yes [1]_   78-78-4
``n-Hexane``                  yes        yes        110-54-3
``n-Heptane``                 yes        yes        142-82-5
``n-Octane``                  yes        yes        111-65-9
``Hydrogen``                  yes        yes        1333-74-0
``Oxygen``                    yes        yes        7782-44-7
``CarbonMonoxide``            yes        yes [1]_   630-08-0
``Water``                     yes        yes        7732-18-5
``Helium``                    yes        yes        7440-59-7
``Argon``                     yes        yes        7440-37-1
``HydrogenSulfide``           no         yes        7783-06-4
``n-Nonane``                  no         yes        111-84-2
``n-Decane``                  no         yes        124-18-5
============================  =========  =========  ==================

.. [1] GERG-2008 replaces the GERG-2004 pure-fluid equation of state for this
   component.  The component itself is present in both models.

Reducing parameters (:math:`\beta`, :math:`\gamma`) exist for **every** binary
pair — 153 for GERG-2004 and 210 for GERG-2008 — so no pair of components is
unsupported.  Fifteen of the GERG-2008 pairs additionally carry a departure
function; the rest take :math:`F_{ij} = 0` and contribute no departure term.

Requesting a fluid outside the model's set raises ``ValueError``::

    AbstractState("GERG2008", "R134a")
    # ValueError: [R134a] (CAS 811-97-2) is not a component of this GERG model

    AbstractState("GERG2004", "n-Decane")
    # ValueError: [n-Decane] is a GERG-2008 component but not a GERG-2004 component

Range of validity
=================

Kunz & Wagner give two ranges for the mixture model:

===============  ==========================  ===================
Range            Temperature                 Pressure
===============  ==========================  ===================
Normal           90 K to 450 K               :math:`p \le` 35 MPa
Extended         60 K to 700 K               :math:`p \le` 70 MPa
===============  ==========================  ===================

**The backend enforces the extended range**, not the normal one.  Inside the
normal range the model's uncertainties are the ones quoted in the publications;
between the normal and the extended limits the model still evaluates but with
larger (and less well characterised) uncertainty.  Choosing the extended range
as the enforced one means the backend does not refuse to answer questions the
authors consider answerable; it does not mean every answer inside it carries
the headline accuracy.

A temperature outside the enforced range raises ``OutOfRangeError``::

    AS.update(CP.PT_INPUTS, 1e6, 800)
    # OutOfRangeError: Temperature [800 K] is outside the GERG range of validity [60, 700] K

Two details of the enforcement are worth knowing.

**The check is on temperature only.**  Pressure is not checked, deliberately.
The 70 MPa ceiling is an operating-envelope statement about the mixture model,
whereas the underlying single-phase equation of state is legitimately evaluated
at :math:`(T, \rho)` points whose pressure is far outside it — including inside
the two-phase dome, where the single-phase equation returns a very large or
negative pressure by construction.  A pressure check cannot distinguish those
two cases, so there is none.  Like every other property-limit guard in
CoolProp, the temperature check is skipped when the
``DONT_CHECK_PROPERTY_LIMITS`` configuration flag is set.

**For a mixture, the enforced limits are the intersection of the per-component
limits** — the largest of the components' ``Tmin`` values and the smallest of
their ``Tmax`` values.  Each GERG component carries ``Tmax = 700 K`` and
``Tmin = min(60 K, T_c)``; the ``Tmin`` cap exists because taking 60 K
literally would put ``Tmin`` *above* the reducing temperature for helium
(5.1953 K) and hydrogen (33.19 K), which is self-contradictory.  So in practice
a mixture's enforced range is exactly 60–700 K — the published mixture-model
range — unless *every* component has :math:`T_c` below 60 K, as in a
helium/hydrogen mixture (33.19–700 K).

Note this is **not** the same as what ``Tmin()`` and ``Tmax()`` report.  Those
are CoolProp's generic mixture accessors, which return a mole-fraction-weighted
sum of the per-component limits; the GERG range check deliberately does not use
them, because that sum is not normalised by :math:`\sum_i x_i` and would
therefore scale the enforced bounds with the magnitude of the composition
vector — a composition summing to 2 would have accepted 900 K.  The enforced
bounds are composition-independent; ``Tmin()``/``Tmax()`` are unchanged and
still report the weighted values.

**Density is not enforced either.**  Each GERG fluid carries
``EOS.limits.rhomax = 1e6`` mol/m³, because ``EquationOfState::limits`` has the
field and leaving it at its default would be worse; nothing in this backend
checks it.  Treat it as metadata, not as a guard.

The same cap has a consequence for the **pure** fluids: because helium's and
hydrogen's ``Tmin`` are 5.1953 K and 33.19 K rather than 60 K, ``GERG2008::Helium``
will happily return a density at 6 K and ``GERG2008::Hydrogen`` at 34 K, with no
error — both far below GERG's own published 60 K floor.  Nothing flags this, so
treat any pure helium or hydrogen result below 60 K as extrapolation.

What these backends deliberately do not provide
===============================================

Each of the following is a deliberate refusal, not a gap waiting to be filled
in.  In every case the alternative would be to return a plausible number
computed from something that is not GERG, labelled as GERG.

Transport properties
--------------------

Viscosity, thermal conductivity and surface tension are **not part of either
GERG model**.  CoolProp has correlations for all three, and they are perfectly
good numbers — from a different model.  Returning them through a backend named
``GERG2008`` invites misattribution, so they throw::

    AS.viscosity()
    # NotImplementedError: Transport properties are not part of the
    # GERG-2004/GERG-2008 models. Use the HEOS backend if you want
    # CoolProp's transport correlations.

    AS.surface_tension()
    # NotImplementedError: Surface tension is not part of the
    # GERG-2004/GERG-2008 models.

If you need transport properties for a natural gas, use the ``HEOS`` backend
(or REFPROP) and be explicit in your reporting that the transport numbers come
from a different source than the thermodynamic ones.

Superancillaries
----------------

GERG fluids carry no :doc:`superancillary </coolprop/SuperAncillary>` Chebyshev
expansions.  This one is load-bearing, because the obvious shortcut is a silent
correctness bug: CoolProp's saturation path returns a superancillary value as
*the answer*, not as an iteration guess.  Attaching CoolProp's existing
methane superancillary to a GERG methane fluid would therefore return
Setzmann-Wagner saturation densities labelled GERG-2008, with no warning and no
iteration to correct them.

Genuine superancillaries would have to be fitted against the 23 distinct GERG
pure equations of state.  That is a reasonable follow-on and purely additive,
but it has not been done.  Until it is, pure-fluid saturation goes through the
classical ancillary-seeded VLE solver, which converges to GERG-consistent
values — correct, just slower than the superancillary path.

Mutable binary interaction parameters
-------------------------------------

``set_binary_interaction_double`` and ``set_binary_interaction_string`` throw::

    AS.set_binary_interaction_double(0, 1, "betaT", 1.0)
    # ValueError: GERG binary interaction parameters are fixed by the published
    # model and cannot be modified. A mixture with altered beta/gamma is not GERG.

``AbstractState`` declares four public mutator routes to the same data — an
index-keyed and a CAS-keyed overload of each of those two functions.  Only the
**index-keyed** pair is overridden here, with the ``ValueError`` above; the two
**CAS-keyed** overloads are inherited from ``AbstractState``, which already
throws ``NotImplementedError`` for them unconditionally.  So the exception type
you see depends on which overload you called — both refuse, and a test pins
that the inherited pair keeps refusing.

A mixture whose :math:`\beta` or :math:`\gamma` has been altered is not GERG,
so the setters refuse rather than producing a mutant model that still answers
to the name.  If you want to adjust interaction parameters, that is what the
``HEOS`` backend is for.

Data GERG does not publish, derived here
===========================================

Two quantities CoolProp needs are not part of the published GERG models.
Borrowing them from a different equation of state is what the rest of this
page refuses to do, so where one can be *derived from the GERG equation
itself* it is derived instead.

.. _gerg_acentric_factor:

Acentric factor
----------------

GERG publishes no acentric factor, but CoolProp needs one.  Three separate
pieces of machinery read :math:`\omega` *directly* off the equation-of-state
object: the Wilson K-factor seed used by every mixture VLE flash and phase
envelope, ``T_DP_PengRobinson``, and the SRK density estimate ``solver_rho_Tp``
starts from.  None of them checks the value first, which is why leaving it at
an ``inf`` sentinel used to break mixture saturation.

Borrowing CoolProp's tabulated value for the same fluid would work numerically
and would break the property this backend exists to have — that nothing in it
traces back to a different equation of state.  So :math:`\omega` is **derived
from the GERG equation itself**, by Pitzer's defining relation

.. math::

    \omega = -1 - \log_{10}\!\left( \frac{p_{\mathrm{sat}}(0.7\,T_c)}{p_c} \right)

where :math:`p_{\mathrm{sat}}(0.7\,T_c)` comes from a *converged* saturation
solve on the GERG pure equation (not from the fitted ``pV`` ancillary, which is
only ever a seed), and :math:`p_c` is the pressure the equation produces at
GERG's tabulated reducing state — the same number ``p_critical()`` reports, so
that the Wilson correlation reads a consistent :math:`p_c`, :math:`T_c` and
:math:`\omega` off one fluid.  There are 23 values for 21 components, because
GERG-2008 moves the reducing parameters of carbon monoxide and isopentane.

``acentric_factor()`` therefore returns a number for a pure GERG fluid (it used
to throw ``NotImplementedError``), and ``PropsSI("acentric", ...)`` returns it
too.  A *mixture* still throws ``ValueError``, but that is CoolProp's rule for
every backend, not a GERG restriction.

The values agree with published acentric factors to better than 0.004 across
all 23 — reassuring, but not a validation target: GERG's reducing point is a
fitted quantity rather than the true critical point, and these come from the
shortened technical form rather than each fluid's reference equation.

**Helium and hydrogen** need one deliberate step outside the rules.  Their
enforced ``Tmin`` equals their reducing temperature (5.1953 K and 33.19 K), so
:math:`0.7\,T_c` — 3.6367 K and 23.2330 K — is below the range ``update()``
lets a caller reach.  That range is the published operating envelope of the
*mixture* model; the pure equations themselves are perfectly well behaved
there (the same VLE solver traces both down to :math:`0.30\,T_c` when fitting
their ancillaries).  Because :math:`\omega` is a definitional constant of the
equation rather than a state anyone queries, it is evaluated there once,
offline, and only the resulting scalar ships.  No run-time guard was weakened.

The table is generated by ``dev/gerg/compute_acentric.py`` into
``src/Backends/GERG/GERGAcentric.h`` and is re-derived and diffed by
``dev/gerg/verify_transcription.py``; the test suite additionally re-derives
every row from CoolProp's own saturation solver and checks it reproduces the
committed value.


Known limitations
=================

These are real and current.  They are stated here rather than discovered later.

.. _gerg_supported_inputs:

Which input pairs actually work
--------------------------------

.. warning::

   For a **pure** GERG fluid, the pressure-plus-caloric input pairs
   (``HmolarP`` / ``HmassP``, ``PSmolar`` / ``PSmass``, ``PUmolar`` /
   ``PUmass``) **do not work at all**.  Through ``PropsSI`` this is not an
   exception you can catch — it is ``inf`` plus an error string.

Measured on a 4x4 grid, :math:`p \in \{10^5, 10^6, 10^7, 5 \times 10^7\}` Pa
and :math:`T \in \{200, 300, 500, 650\}` K, round-tripping a state first
obtained from ``PT_INPUTS`` on the same backend:

======================  =============  =============  =============  =============
Backend / fluid         ``HmolarP``    ``PSmolar``    ``PUmolar``    ``DmolarP``
======================  =============  =============  =============  =============
``HEOS::Methane``       16/16          16/16          16/16          16/16
``GERG2008::Methane``   **0/16**       **0/16**       **0/16**       16/16
``GERG2008::Nitrogen``  **0/16**       **0/16**       **0/16**       16/16
``GERG2008::CO2``       **0/16**       **0/16**       **0/16**       16/16
======================  =============  =============  =============  =============

So::

    h = PropsSI("Hmolar", "P", 1e6, "T", 300, "GERG2008::Methane")
    PropsSI("T", "P", 1e6, "Hmolar", h, "GERG2008::Methane")
    # -> inf   (and an errstring; NOT an exception)

``DmolarP`` used to be in the same state (3/16, 2/16, 5/16 above).  It failed
through ``FlashRoutines::T_DP_PengRobinson``, which reads the acentric factor;
now that the backend :ref:`derives one from its own equation
<gerg_acentric_factor>` that pair round-trips everywhere on this grid.

The three caloric pairs are blocked by two *different* causes, neither of which
is the acentric factor:

1. **No triple-point pressure.**  GERG publishes no triple point, so
   ``p_triple()`` is the internal ``inf`` sentinel.  The single-phase
   pressure-plus-caloric flash brackets temperature, and its gas-phase branch
   tests ``p < p_triple()`` to decide whether a saturation temperature exists.
   With ``p_triple()`` infinite that test is always true, so the bracket's
   lower end drops to ``Tmin()`` = 60 K instead of the saturation temperature
   (149.1 K for methane at 1 MPa).  At 60 K and 1 MPa methane is a compressed
   liquid while the solver has been told the phase is gas, and the density
   iteration diverges to a negative root.
2. **The bracket's upper end is outside the model's range.**  The same solver
   takes ``1.5 * Tmax()`` = 1050 K as its upper bound, and this backend's own
   range guard rejects any temperature above 700 K — including the internal
   ``update()`` calls the bracketing solver makes.  Clearing cause 1 in a local
   experiment moves the failure straight to this one.

**What does work**, and is covered by the test suite: ``PT``, ``DmolarT`` /
``DmassT``, ``DmolarP`` / ``DmassP``, and — for pure fluids — ``QT``, ``PQ``
and ``DmolarQ``.

The caloric-pair restriction affects **pure fluids only**.  The same
``HmolarP`` flash on a GERG *mixture* succeeds, because the mixture code path
takes a different density-guess route.  Until this is fixed, obtain pure-fluid
states from ``PT`` or ``DmolarT`` inputs, or use ``HEOS`` when you need a
caloric input pair.

Mixture ``DQ`` and ``DmolarP`` flashes
---------------------------------------

Mixture saturation, phase envelopes and VLE flashes **work** — this section
used to say they did not.  A ``QT`` or ``PQ`` flash on a GERG mixture and
``build_phase_envelope()`` on a GERG mixture all converge::

    AS = CP.AbstractState("GERG2008", "Methane&Ethane")
    AS.set_mole_fractions([0.9, 0.1])
    AS.update(CP.QT_INPUTS, 0.0, 150)   # bubble point: p = 9.2408e5 Pa
    AS.update(CP.QT_INPUTS, 1.0, 150)   # dew point:    p = 9.3875e4 Pa
    AS.build_phase_envelope("")         # 222 points

They did not before.  CoolProp's mixture VLE machinery seeds itself with
Wilson K-factors and an SRK density estimate, both of which read the acentric
factor, and GERG publishes none — so the fluids carried the internal
:math:`\omega = +\infty` sentinel, the Wilson pressure estimate collapsed to
exactly zero, and the flash failed several solver frames later with
``solver_rho_Tp was unable to find a solution for T=150, p=0, with guess value
nan``.  The fix was to *derive* an acentric factor from GERG's own equation
rather than borrow one; see :ref:`Acentric factor <gerg_acentric_factor>`.

**Pure-fluid** saturation was never affected and is unchanged: ``QT``, ``PQ``
and ``DQ`` inputs on ``GERG2008::Methane`` converge normally.  All 23 GERG pure
equations of state carry a fitted ancillary; the ones reachable through the
public API are those whose saturation curve lies inside the enforced
temperature range, which is 19 of the 21 GERG-2008 components — helium and
hydrogen are the exceptions and throw ``OutOfRangeError`` (see below).

A mixture ``DQ`` or ``DmolarP`` flash still fails, but for an unrelated and
pre-existing reason — ``DQ_flash not ready for mixtures`` /
``DP_flash not ready for mixtures`` — which fails identically on ``HEOS`` and
is not a GERG limitation.  See
:ref:`Which input pairs actually work <gerg_supported_inputs>` for the
pure-fluid caloric-input restriction, which has a different root cause.

Compositions with two or more exactly-zero mole fractions
----------------------------------------------------------

The natural way to hand over a natural-gas analysis is the full 21-name
GERG-2008 component list with a mole fraction for each, most of them exactly
zero.  Until this release *every* such composition returned NaN for *every*
property, silently.  It is now fixed as far as it mathematically can be, and
the boundary is worth understanding, because it is a property of the function
rather than a choice about how much work to do.

Why there is a boundary at all
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reducing function sums pairwise terms

.. math::

    f_{Y,ij}(x_i, x_j) = \frac{x_i x_j (x_i + x_j)}{\beta_{Y,ij}^2 x_i + x_j}

which is :math:`0/0` as soon as two mole fractions are exactly zero.  The
numerator is cubic and the denominator linear, so :math:`f_{Y,ij}` is
**homogeneous of degree 2**:

.. math::

    f_{Y,ij}(t x_i, t x_j) = t^2 f_{Y,ij}(x_i, x_j), \qquad t > 0.

Differentiating that identity :math:`k` times in the compositions lowers the
degree by :math:`k`, so the :math:`k`-th composition derivative is homogeneous
of degree :math:`2 - k`.  Restricted to **non-negative** mole fractions the
denominator is bounded away from zero on the unit directions —
:math:`\beta^2 a + b \ge \min(\beta^2, 1)` for :math:`a, b \ge 0` with
:math:`a^2 + b^2 = 1` — so that direction set is compact and each derivative is
bounded on it by some :math:`M`.  Writing an approach to the corner as
:math:`t \mathbf{u}` with :math:`\mathbf{u}` a unit direction,

.. math::

    \left| \partial^k f_{Y,ij}(t \mathbf{u}) \right| \le M \, t^{\,2-k}.

That single inequality decides everything:

* :math:`k = 0` and :math:`k = 1` — the bound vanishes as :math:`t \to 0`, and
  it does so **uniformly in the direction** :math:`\mathbf{u}`.  The limit is
  therefore :math:`0` no matter how the corner is approached, and defining the
  value to be :math:`0` there is the unique *continuous extension*.  It is not
  a convention or a fudge: the extended function is genuinely :math:`C^1`.
  Physically this is just the statement that a pair of components which are
  both absent contributes nothing, and that adding a little of *both* changes
  the reducing state only at second order in that amount.
* :math:`k \ge 2` — the exponent is zero or negative.  A degree-zero function
  is *constant along each ray* but generally *different between rays*, so it
  has no limit at the corner and no value is correct.  Measured at
  :math:`\beta = 1.3`, :math:`\partial^2 f / \partial x_i^2` is
  :math:`-7.09 \times 10^{-2}` approaching along :math:`(1,1)`,
  :math:`-8.24 \times 10^{-1}` along :math:`(1,9)` and
  :math:`-3.24 \times 10^{-4}` along :math:`(9,1)` — identical at
  :math:`t = 10^{-3}`, :math:`10^{-5}` and :math:`10^{-7}`, so this is genuine
  direction dependence rather than slow convergence.

So the value and the first derivatives are repairable and are repaired; the
second and higher derivatives are not, and are deliberately left as NaN.
Substituting any number for them would replace an honest NaN with a silently
direction-dependent one.

.. note::

   The uniform bound is what fails outside the physical domain.  With
   :math:`\beta = 1` and :math:`x_i = -x_j` the denominator cancels while both
   mole fractions are non-zero, and the direction set is no longer bounded away
   from the singularity.  That corner is unreachable from a valid composition
   and is a separate pre-existing defect, not something these guards address.

What works
~~~~~~~~~~

The reducing state and everything that flows from it — :math:`T_r`,
:math:`\rho_r`, :math:`\alpha^0`, :math:`\alpha^r`, :math:`p`, :math:`\rho`,
:math:`c_v`, :math:`c_p`, speed of sound, :math:`h`, :math:`s`, molar mass —
and, since this release, the **first** composition derivatives and therefore
``fugacity()`` and ``fugacity_coefficient()``.  All of these are identical, to
the last digit, to the same composition with the zero components trimmed away::

    AS = CP.AbstractState("GERG2008", "Methane&Nitrogen&Ethane&Propane")
    AS.set_mole_fractions([0.9, 0.1, 0.0, 0.0])
    AS.update(CP.PT_INPUTS, 1e6, 300)
    AS.rhomolar()                 # 406.942     correct
    AS.fugacity_coefficient(0)    # 0.983354    correct

Five call sites contain the :math:`0/0`, and all five are now guarded: three
helper functions — ``f_Y_ij`` itself and its two first-derivative helpers
``dfYkidxi__constxk`` / ``dfYikdxi__constxk`` — plus, added in this release,
the two *inlined* expressions in the ``XN_DEPENDENT`` branch of
``dYrdxi__constxj``, which duplicate those helpers rather than calling them.  That last one
was the reason a single trailing zero was so destructive: the inlined loop runs
over every component, so with ``x[N-1] == 0`` one further zero anywhere in the
composition made the whole first derivative NaN.

What still does not work
~~~~~~~~~~~~~~~~~~~~~~~~

The second and higher composition derivatives, per the argument above.  These
are not reachable from the Python API; at the C++ level, with ``z`` the mole
fraction vector used above,

.. code-block:: cpp

    heos->Reducing->d2Trdxidxj(z, 0, 1, XN_DEPENDENT);   // nan, and correctly so

In practice this means **phase envelopes**, which need
:math:`\partial \ln \varphi_i / \partial x_j`.  Both backends throw
``Unable to calculate at least 4 points in phase envelope; quitting`` on a
composition with two exact zeros, and both build the envelope normally once the
zero-mole-fraction components are trimmed out.  ``PQ`` and ``QT`` flashes
succeed with the zeros present.

If you need phase envelopes, **trim the zero-mole-fraction components out of
the composition** rather than passing the full component list.  The trimmed
result is exact.

**This is not GERG-specific.**  The reducing function is shared, so the same
behaviour — before and after — applies to the default ``HEOS`` backend, and it
predates these backends entirely.  It is tracked as
`GitHub #1677 <https://github.com/CoolProp/CoolProp/issues/1677>`_.

Saturation states below the enforced ``Tmin``
----------------------------------------------

The published range of validity, 60–700 K, is applied **verbatim to every
component**.  GERG states it for the mixture model as a whole and publishes no
per-component lower limit, so the backend does not invent one.

Seven components have a fitted saturation ancillary whose low-temperature end
lies *below* 60 K: methane (57.17 K), nitrogen (37.86 K), oxygen (46.41 K),
carbon monoxide (39.86 K), argon (45.24 K), hydrogen (9.96 K) and helium
(1.56 K).  ``get_state("triple_liquid")`` — and, from C++, ``calc_Tmin_sat()``
/ ``calc_pmin_sat()`` on the Helmholtz backend — therefore report a state that
``update()`` will refuse to evaluate::

    AS = CP.AbstractState("GERG2008", "Methane")
    AS.update(CP.QT_INPUTS, 0.0, 58.0)
    # OutOfRangeError: Temperature [58 K] is outside the GERG range of validity [60, 700] K

The ancillary data below 60 K is real and was traced with teqp; it is simply
not reachable through the public API, because the model's own range of validity
stops first.

Helium and hydrogen are supercritical-only as pure fluids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Helium is critical at 5.1953 K and hydrogen at 33.19 K, both far below the
60 K lower limit.  As **pure** fluids in these backends they are therefore
supercritical across the entire valid range, and no saturation state can be
requested at all::

    AS = CP.AbstractState("GERG2008", "Helium")
    AS.update(CP.QT_INPUTS, 0.5, 4.5)
    # OutOfRangeError: Temperature [4.5 K] is outside the GERG range of validity [60, 700] K

That is the physics rather than a defect: at 60 K helium *is* supercritical.
Single-phase properties work normally, and both are ordinary components in a
**mixture**, where the mixture reducing temperature governs rather than the
pure-component critical temperature.

An earlier version of this backend capped ``Tmin`` at each component's own
critical temperature, which gave helium and hydrogen a degenerate range whose
lower end equalled ``Tcrit`` — no contradiction in the numbers, but a limit the
model does not publish, and a silently empty saturation dome instead of a clear
refusal.  If you need subcritical helium or hydrogen, use ``HEOS``.

Properties GERG does not define at all
---------------------------------------

GERG publishes no triple point.  ``Ttriple()`` returns 0, ``p_triple()``
returns ``inf``, and ``get_state("triple_liquid")`` returns the
low-temperature end of the fitted saturation curve under a name CoolProp
inherited — none of them is a GERG triple point, because there is no such
thing in these models.  The ``p_triple()`` sentinel is not inert: it is one of
the two reasons the pure-fluid caloric input pairs still fail (see
:ref:`Which input pairs actually work <gerg_supported_inputs>`).

``set_reference_stateS`` and ``set_reference_stateD`` are not available
-----------------------------------------------------------------------

Both throw ``NotImplementedError`` on these backends::

    CP.set_reference_state("GERG2008::Methane", "NBP")
    # NotImplementedError: set_reference_stateS is not implemented for the
    # GERG2008 backend. ...

    CP.set_reference_state("GERG2008::Methane", 300.0, 1000.0, 0.0, 0.0)
    # NotImplementedError: set_reference_stateD is not implemented for the
    # GERG2008 backend. ...

CoolProp applies a reference-state change by writing an offset into the global
fluid-library entry for the fluid, and the GERG backends do not read that
library.  Before this was made explicit the ``S`` call was a **silent no-op** —
it returned without error and without effect, and did not even validate the
reference-state string.  The ``D`` variant, which ends in the same
library write, failed in two different unhelpful ways instead: with a backend
prefix it raised an internal ``key [GERG2008::Methane] was not found in
string_to_index_map in JSONFluidLibrary``, and without one it silently adjusted
``HEOS``.  See *Reference state* below for what to do instead.

The throw covers every spelling that resolves to a GERG family, including
``GERG2008Backend::``, ``GERG2008?<options>::`` and composed strings such as
``BICUBIC&GERG2008::``.  It does **not** cover a bare fluid name with no
``::`` prefix: ``set_reference_stateS("Methane", "NBP")`` resolves to the
default backend and adjusts the ``HEOS`` fluid library, which has no effect on
a GERG state and raises no error.  That is inherent to the string API — the
name carries no backend information — so always pass the prefix.

The strictness rules are model-level, not C++-level
----------------------------------------------------

The throws listed under *What these backends deliberately do not provide*
cover the documented API.  They are not a sandbox.

Two further routes that are reachable from a plain ``AbstractState`` — with no
downcast — are also **closed**, and are not listed above only because they are
not "properties GERG does not publish":

* ``change_EOS(i, "SRK")`` throws ``ValueError``.  Swapping the equation of
  state out from under a GERG backend leaves an object that still answers to
  the name ``GERG2008`` while computing from a cubic; it was confirmed to
  install an SRK cubic silently before this refusal existed.
* ``update_with_guesses(...)`` applies the same range check as ``update()``.
  It is a separate virtual entry point and previously accepted
  :math:`T = 900` K.

What remains open needs C++ that already holds a concrete backend pointer or
reaches through a public data member — for example
``Reducing->set_binary_interaction_double(...)`` on the reducing-function
object, direct assignment into ``residual_helmholtz->Excess.F[i][j]``, or
``update_DmolarT_direct()``, which bypasses the range check by design because
it is what the backend uses to build its own fluids.  These are enumerated in
``GERGBackend.h`` and pinned by tests, so closing or reopening one shows up as
a test change.  The strictness rules exist to stop a *plausible mistake*, not a
determined one.

Tabular backends wrapping GERG
-------------------------------

``BICUBIC&GERG2008`` and ``TTSE&GERG2008`` are **not supported**.  They do not
fail loudly, which is the problem: the table build completes without error and
persists a cache under ``~/.CoolProp/Tables/GERG2008Backend(...)``, and then
every lookup inside the model's own range is rejected::

    AS = CP.AbstractState("BICUBIC&GERG2008", "Methane")
    AS.update(CP.PT_INPUTS, 1e6, 300)
    # ValueError: inputs are not in range, p=1e+06 Pa, T=300 K

Delete the cache directory if you have created one.  Use the ``GERG2008``
backend directly.

Both halves of that — the factory succeeding, and every subsequent lookup
being rejected — are pinned by a test, so if a future change makes these
wrappers work, this section starts failing rather than quietly going stale.

Helium and hydrogen have no reachable saturation state
-------------------------------------------------------

Because each component's ``Tmin`` is capped at its own reducing temperature
(see *Range of validity* above), helium and hydrogen end up with
``Tmin == T_c`` — 5.1953 K and 33.19 K respectively.  Their entire subcritical
saturation curve is therefore below the enforced lower temperature limit, and
any saturation call on them raises ``OutOfRangeError``::

    AS = CP.AbstractState("GERG2008", "Helium")
    AS.update(CP.QT_INPUTS, 0.0, 4.67577)
    # OutOfRangeError: Temperature [4.67577 K] is outside the GERG range of
    # validity [5.1953, 700] K

This is consistent with the model: GERG's published lower limit is 60 K, so
neither fluid's saturation curve is inside the model's range in the first
place.  Both components exist in GERG because they appear as dilute
constituents of natural gas, not because GERG is a helium or hydrogen
saturation model.

The one place this limit is deliberately stepped past is the
:ref:`acentric factor <gerg_acentric_factor>`, which is defined at
:math:`0.7\,T_c` and is a constant of the equation rather than a state a
caller queries.  It is evaluated offline, once, and only the resulting scalar
is compiled in; the run-time guard above is unchanged.

Reference state
---------------

.. warning::

   These backends report :math:`h = 0` and :math:`s = 0` for the **ideal gas**
   at 298.15 K and 101325 Pa, per component.

That is not any of CoolProp's usual reference states (IIR, ASHRAE, NBP), and it
is the *ideal-gas* value at that state, not the real-fluid value.  The
published integration constants are deliberately discarded and recomputed to
satisfy it, exactly as teqp does; without that recomputation :math:`h` and
:math:`s` would not match teqp even though :math:`p`, :math:`c_v` and the speed
of sound did.

The practical consequence: **enthalpy and entropy from a GERG backend differ
from the same property from ``HEOS`` by a large offset.**  For methane at
101325 Pa,

===========  ============================  ============================
:math:`T`    :math:`h_{GERG} - h_{HEOS}`   :math:`s_{GERG} - s_{HEOS}`
===========  ============================  ============================
300 K        -14614.05 J/mol               -107.1158 J/mol/K
400 K        -14614.26 J/mol               -107.1164 J/mol/K
500 K        -14613.81 J/mol               -107.1154 J/mol/K
===========  ============================  ============================

This is correct and it matches teqp.  The offset is essentially constant — it
varies by about 3 parts in :math:`10^5` across the table above, which is the
genuine difference between the two models' ideal-gas and residual parts, not a
reference-state effect.

**Any cross-check of a GERG backend against ``HEOS`` must compare *differences*
in** :math:`h` **and** :math:`s` **, never absolute values.**  Comparing
absolute enthalpies will make a correct implementation look catastrophically
wrong.  Neither ``set_reference_stateS`` nor ``set_reference_stateD`` is
available as an escape hatch here
(see above) — subtract your own offset, or take a single reference point from
each backend and compare everything relative to it.

Relationship to CoolProp's default HEOS mixture model
=====================================================

CoolProp's default multi-fluid (``HEOS``) backend descends from the same
Kunz & Wagner formulation, so it is structurally GERG-shaped.  It is
nevertheless a **different model** and gives different numbers.  Two independent
reasons:

**Pure-fluid equations of state.**  This is the load-bearing one.  ``HEOS``
ships the *reference* equation of state for each fluid — Setzmann-Wagner for
methane, Span-Wagner for carbon dioxide, and so on.  GERG uses its own
shortened technical form, 12 to 24 terms with a largely shared exponent set.
These are different equations producing different numbers.

**Binary interaction parameters.**  All 210 GERG-2008 binary pairs are present
in CoolProp's binary-pair library, and 194 of them are the GERG values.  The
other **16 are later refits** — 15 from Gernert's thesis and 1 from
Tkaczuk et al. (2020) — which is exactly what a general-purpose library should
prefer, and exactly what a backend named ``GERG2008`` must not use.  Twelve of
those 16 shift the mixture reducing temperature by more than 0.03 K at
:math:`z = (0.35, 0.65)`; four leave :math:`\beta`/:math:`\gamma` unchanged but
attach a *different* departure function.

**The gas constant.**  GERG specifies
:math:`R = 8.314472\ \mathrm{J\,mol^{-1}\,K^{-1}}`.  Under CoolProp's default
``NORMALIZE_GAS_CONSTANTS`` configuration the ``HEOS`` backend uses the CODATA
value for mixtures, which rescales :math:`p`, :math:`\alpha^{ig}`, :math:`c_v`
and :math:`w` by about :math:`1.1 \times 10^{-6}`.  The GERG backends override
this.

The scale of the difference is small but systematic.  For
90 % methane / 10 % nitrogen at 300 K and 1 MPa:

======================  ==========================
Backend                 :math:`\rho` [mol/m³]
======================  ==========================
``GERG2008``            406.9418737
``HEOS``                406.9525552
======================  ==========================

— a relative difference of :math:`2.6 \times 10^{-5}`, which is small compared
to a typical custody-transfer tolerance but is not numerical noise.  If you
need GERG numbers, ask for GERG.

Note also that this is *not* the same thing as REFPROP's ``REFPROP_USE_GERG``
flag, which substitutes GERG pure-fluid equations of state inside REFPROP's own
mixture model.  These backends implement GERG whole.

Implementation notes and references
===================================

The reference implementation is `teqp <https://github.com/usnistgov/teqp>`_,
specifically ``include/teqp/models/GERG/GERG.hpp``.  The coefficient tables in
CoolProp's GERG backend are transcribed from it, with a durable verification
script (``dev/gerg/verify_transcription.py``) that checks every table family
against teqp, and the test suite compares against teqp-generated reference
values at relative tolerances of 1e-12 for :math:`\alpha^r` and
:math:`\alpha^{ig}` and 1e-10 for :math:`p`, :math:`c_v` and the speed of
sound.  teqp's own values are in turn checked against the AGA8 reference
implementation.

The design of record for these backends — including the rationale for every
strictness rule above — is
``docs/superpowers/specs/2026-07-25-gerg-strict-backend-design.md`` in the
CoolProp repository.

The models themselves are published in:

* GERG-2004: Kunz, Klimeck, Wagner & Jaeschke, *The GERG-2004 Wide-Range
  Equation of State for Natural Gases and Other Mixtures*, GERG Technical
  Monograph 15, VDI Verlag, 2007 :cite:`Kunz-BOOK-2007`.
* GERG-2008: Kunz & Wagner, *J. Chem. Eng. Data* **57**, 3032-3091, 2012
  :cite:`Kunz-JCED-2012`.
