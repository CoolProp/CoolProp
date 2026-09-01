.. _RES:

*********************************
Residual Entropy Scaling (RES)
*********************************

.. contents:: :depth: 2

Residual entropy scaling correlates a transport property against the *residual*
entropy :math:`s^{\rm res}` of the fluid rather than against temperature and
density directly.  Because the residual entropy comes from the equation of
state, one set of per-fluid coefficients covers gas, liquid and supercritical
states with a single expression, and the same coefficients can be re-fitted for
a different equation of state.

CoolProp ships the models of two papers:

* **Viscosity** — Martinek, Bell, Herzog, Richter, Yang, *Entropy Scaling of
  Viscosity IV — Application to 124 Industrially Important Fluids*,
  J. Chem. Eng. Data **70**, 727–742 (2025),
  `10.1021/acs.jced.4c00451 <https://doi.org/10.1021/acs.jced.4c00451>`_
* **Thermal conductivity** — Li, Duan, Yang, *Linking Thermal Conductivity to
  Equations of State Using the Residual Entropy Scaling Theory*,
  Ind. Eng. Chem. Res. **63**, 18160–18175 (2024),
  `10.1021/acs.iecr.4c02946 <https://doi.org/10.1021/acs.iecr.4c02946>`_

with the coefficients for the cubic equations of state from Yang, *Viscosity and
Thermal Conductivity Models of 151 Common Fluids Based on Residual Entropy
Scaling and Cubic Equations of State*, ACS Omega (2025),
`10.1021/acsomega.4c10815 <https://doi.org/10.1021/acsomega.4c10815>`_.

Using it
========

On ``HEOS`` the model is **opt-in**, per property, through the
:ref:`backend options <backend-options>` suffix on the factory string.  It is
never enabled implicitly: most HEOS fluids already have a reference transport
correlation, and silently replacing one would change results for existing
callers::

    from CoolProp.CoolProp import PropsSI, AbstractState

    PropsSI("V", "T", 400, "P", 5e6, 'HEOS::Propane?{"RES":{"viscosity":true}}')

    AS = AbstractState("HEOS", 'Propane?{"RES":{"viscosity":true,"conductivity":true}}')
    AS.update(CoolProp.PT_INPUTS, 5e6, 400)
    AS.viscosity(), AS.conductivity()

On the **cubic** backends (``PR``, ``SRK``) there is nothing to opt into: those
equations of state carry no transport model at all, so RES is used automatically
wherever the fluid has coefficients, and ``viscosity()`` keeps raising
``NotImplementedError`` where it does not::

    PropsSI("V", "T", 400, "P", 5e6, "PR::Propane")

Options
-------

====================================== ========= ===============================
Key                                    Default   Meaning
====================================== ========= ===============================
``RES.viscosity``                      ``false`` Use the RES viscosity model.
``RES.conductivity``                   ``false`` Use the RES conductivity model.
``RES.mixture_critical_enhancement``   ``false`` Apply the critical enhancement
                                                 to mixtures as well as to pure
                                                 fluids.
====================================== ========= ===============================

Unknown keys and wrong types are rejected when the state is constructed, so a
typo surfaces immediately rather than leaving the model quietly switched off.

Which fluids are covered
========================

122 fluids carry RES parameters.  Not all of them on every backend:

* **PR / SRK** — 112 fluids.  Coefficients were regressed against those
  equations of state directly (Yang 2025), so they apply exactly.
* **HEOS** — 108 fluids.  The Helmholtz coefficients are the *REFPROP* fits from
  the two papers, used here as an approximation.  For 14 fluids they do not
  transfer well enough, and those entries are deliberately **withheld**: the
  backend raises rather than returning a plausible but wrong number.

  The withheld fluids are BENZENE, D5, HEPTANE, MD3M, MD4M, R1123, R1224yd(Z),
  R1233zd(E), R1234yf, R1243zf, R13, R161, R41 and VinylChloride.

The critical-enhancement parameters are **not** per-backend: the same
``phi0``, ``Gamma``, ``q_D``, ``t_ref`` and ``gamma`` are used on HEOS, PR and
SRK.  That follows the sources rather than departing from them — Yang 2025
applies the enhancement to the cubics too, and the table published with it is
the same as Li's to rounding for 147 of 151 fluids (the differences are
1.2415 against 1.242 and the like; in four cases Yang drops an enhancement Li
keeps, never the reverse).  Making them backend-specific would not repay the
extra table.

The exclusion list is measured — see
``HEOS_TRANSFER_EXCLUDE`` in ``dev/convert_RES_csv_to_json.py`` for the
procedure and ``dev/RES_reference/README.md`` for the reference implementations
it is measured against.  Across the 108 fluids that are shipped, the median
deviation from the same model evaluated on REFPROP is 0.00004 % (viscosity) and
0.00007 % (conductivity) at states where the residual term dominates.

Where CoolProp differs from the papers
======================================

**The dilute terms always use the fitted polynomials.** Viscosity and thermal
conductivity are both expressed as a sum of a dilute term and the residual term.  
Martinek 2025 defaults to REFPROP for the dilute term of pure fluids with a fitted
polynomial as fallback and uses the polynomials + Wilke rule for mixtures.  Li 2024
does the reverse (polynomials for pure fluids, REFPROP for mixtures).  CoolProp always 
uses the polynomials for self-consistency and because the cubic backends have no 
native viscosity model.

**The critical enhancement uses the RES viscosity.**  Li 2024 feeds REFPROP's
own viscosity model into the Olchowy–Sengers term.  CoolProp uses the RES
viscosity instead for the same reasons given above.

A third point is worth knowing but is not a divergence: opting into RES
*viscosity* alone slightly changes the reference *conductivity*, because the
reference conductivity's own critical enhancement consumes the fluid's
viscosity.

Mixtures
========

Both models extend to mixtures: the dilute term is combined with a Wilke rule
and the residual coefficients are mole-fraction averaged.  Two caveats:

* The critical enhancement is **off** for mixtures unless
  ``mixture_critical_enhancement`` is set.  It needs the mixture critical point,
  which HEOS has to solve for, and the physical case for a critical enhancement
  in mixtures is not well established in the first place.  Activating it may increase
  accuracy in the critical region for some systems, but costs a significant increase
  in evaluation time and adds a risk of failure when the mixture critical point 
  cannot be located.
* Only binary mixtures are covered by the published validation.
* Everything else about mixtures — the Wilke rule for the dilute term, the
  mole-fraction averaging of the residual coefficients — works without the
  enhancement and is validated in the same way as the pure fluids.

Re-fitting for a different alpha function
=========================================

The cubic coefficients were regressed against each backend's default alpha
function.  Changing it with ``set_cubic_alpha_C`` invalidates them, and RES then
refuses rather than using coefficients that no longer apply.  Supply a re-fit
through the existing per-fluid parameter interface::

    for k, value in enumerate(new_viscosity_coefficients):        # n0, n1, n2
        AS.set_fluid_parameter_double(0, f"RES_viscosity_n{k}", value)
    AS.set_fluid_parameter_double(0, "RES_viscosity_xita", 1.0)

    for k, value in enumerate(new_conductivity_coefficients):     # n0 .. n3
        AS.set_fluid_parameter_double(0, f"RES_conductivity_n{k}", value)
    AS.set_fluid_parameter_double(0, "RES_conductivity_xita", 1.0)

**Every** coefficient of a property has to be written before that property will
evaluate again — all three (viscosity) or four (conductivity) ``n`` values *and*
``xita``.

The two properties are tracked separately: re-fitting the viscosity does not
release the conductivity.  Note that the conductivity's critical enhancement
consumes the RES viscosity, so a conductivity-only re-fit is refused until the
viscosity is re-fitted too.

``get_fluid_parameter_double`` reads the same keys back, on pure fluids and
mixtures alike.
