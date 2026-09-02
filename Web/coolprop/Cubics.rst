.. _cubic_backend:

************************
Cubic Equations of State
************************

.. contents:: :depth: 2

Introduction
============

CoolProp (as of version 6) comes with two standard cubic equations of state: Soave-Redlich-Kwong (SRK) and Peng-Robinson (PR).  These two equations of state can be expressed in a common, generalized form:

.. math::

  p = \frac{RT}{v-b} + \frac{a}{(v+\Delta_1b)(v+\Delta_2b)}

where for pure fluids, :math:`a` and :math:`b` are not composition dependent, whereas for mixtures, they have composition dependence.  These cubic EOS can be converted to a form that is compatible with the multi-fluid model used in CoolProp according to the analysis in Bell and Jager :cite:`Bell-JRN-2016`.

The motivations for the use of cubic EOS are twofold:

* The only required information for the EOS are :math:`T_c`, :math:`p_c`, and the acentric factor of the pure fluids
* They are much more computationally efficient (see below)

Caveats
-------

.. warning:: NOT ALL PROPERTIES ARE AVAILABLE AS INPUTS/OUTPUTS

Only a limited subset of properties are available currently. You can do:

* Flash calculations with TP, PQ, DT, QT inputs
* Calculation of mixture critical point(s)
* Calculation of some mixture flashes

Pure Fluids
===========

Usage
-----

As a user, in the :ref:`high-level interface <high_level_api>`, all that you have to do to evaluate properties using the SRK or PR backends is to change the backend name.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # The multi-parameter Helmholtz backend
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"HEOS::Propane")

    # SRK
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"SRK::Propane")

    # Peng-Robinson
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"PR::Propane")

The same holds for the :ref:`low-level interface <low_level_api>`:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("SRK", "Propane"); AS.update(CP.QT_INPUTS, 0, 300); print(AS.p())

    In [0]: AS = CP.AbstractState("PR", "Propane"); AS.update(CP.QT_INPUTS, 0, 300); print(AS.p())

All the fluids available in CoolProp are also available through the cubic equations of state.  The fluids can be extended according to the analysis shown below

Speed
-----

The increase in speed for evaluating properties from a cubic EOS is one of the primary motivations.  Here we show an example of the speedup to VLE calculations:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # The multi-parameter Helmholtz backend
    In [0]: AS = CP.AbstractState("HEOS", "Propane")

    In [0]: %timeit AS.update(CP.QT_INPUTS, 0, 300)

    # The cubic SRK backend
    In [0]: AS = CP.AbstractState("SRK", "Propane")

    In [0]: %timeit AS.update(CP.QT_INPUTS, 0, 300)

And here, we run the PT flash

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # The multi-parameter Helmholtz backend
    In [0]: AS = CP.AbstractState("HEOS", "Propane")

    In [0]: AS.specify_phase(CP.iphase_gas)

    In [0]: %timeit AS.update(CP.PT_INPUTS, 101325, 300)

    # The cubic SRK backend
    In [0]: AS = CP.AbstractState("SRK", "Propane")

    In [0]: AS.specify_phase(CP.iphase_gas)

    In [0]: %timeit AS.update(CP.PT_INPUTS, 101325, 300)

As you can see, the speed difference is quite significant

Superancillaries for SRK and PR
-------------------------------

For an additional speedup on saturation calculations, piecewise Chebyshev *superancillary* expansions (Bell and Deiters :cite:`Bell-IECR-2021`, ported from teqp) are available for pure fluids with the SRK and PR backends.  They return :math:`p^{\mathrm{sat}}(T)`, :math:`\rho'(T)`, and :math:`\rho''(T)` directly from closed-form polynomials — no iterative flash — at an accuracy better than 10\ :sup:`-3` relative to the cubic EOS result over the full saturation dome up to :math:`T_c`.

Two entry points are provided on the :cpapi:`AbstractState <CoolProp::AbstractState>` (low-level interface only):

* :cpapi:`update_QT_pure_superanc(Q, T) <CoolProp::AbstractState::update_QT_pure_superanc>` — populates the full two-phase state (same fields as ``update(QT_INPUTS, Q, T)``) using the superancillary instead of a full VLE flash.
* :cpapi:`saturation_ancillary(param, Q, iT, T) <CoolProp::AbstractState::saturation_ancillary>` — returns a single saturation property (``iP`` or ``iDmolar``) without updating state.

.. warning:: Superancillaries are only available for **pure** fluids with the **SRK** and **PR** backends, and only with **T** as the independent variable (``T`` must satisfy :math:`T \le T_c`).  They are *not* used transparently by ``PropsSI`` — you must use the low-level interface explicitly.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("PR", "Propane")

    # Full VLE flash through the cubic EOS
    In [0]: %timeit AS.update(CP.QT_INPUTS, 1, 300)

    # Same state populated from the superancillary
    In [0]: %timeit AS.update_QT_pure_superanc(1, 300)

    # Single saturation property (no state update)
    In [0]: %timeit AS.saturation_ancillary(CP.iP, 0, CP.iT, 300)

On a typical modern laptop this corresponds to roughly a 5--6x speedup for ``update_QT_pure_superanc`` relative to the full cubic flash, and roughly 9--10x for ``saturation_ancillary`` on a single property; exact numbers are hardware-dependent.

Mixtures
========

Interaction Parameters
----------------------

For mixtures, cubic EOS (and their modifications) are heavily used.  Cubic EOS allow for reasonable prediction of VLE with only one adjustable parameter, at least for reasonable mixtures.  CoolProp does not come with any values for the :math:`k_{ij}` parameter, but it is straightforward to add it yourself.

.. warning:: The ability to adjust interaction parameters is ONLY available in the :ref:`low-level interface <low_level_api>`

.. warning:: When you call ``set_binary_interaction_double``, it only applies to the given instance of the AbstractState

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("SRK", "Methane&Ethane")

    In [0]: AS.set_mole_fractions([0.5, 0.5])

    In [0]: AS.update(CP.QT_INPUTS, 0, 120); print(AS.p())

    In [0]: AS.set_binary_interaction_double(0,1,"kij",-0.05)

    In [0]: AS.update(CP.QT_INPUTS, 0, 120); print(AS.p())

Critical Points
---------------

Following the approach of Bell *et al.*, it is possible to calculate critical points of mixtures from cubic EOS.  For instance, here is how to calculate all the critical points (there can be more than one) that are found for an equimolar mixture of methane and ethane:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("SRK", "Methane&Ethane")

    In [0]: AS.set_mole_fractions([0.5, 0.5])

    In [0]: pts = AS.all_critical_points()

    In [0]: [(pt.T, pt.p, pt.rhomolar, pt.stable) for pt in pts]

Cubics in multi-fluid model
---------------------------

The cubic equations of state can also be used to replace a single fluid in a multi-fluid model (the GERG-like model) (see Bell and Jager :cite:`Bell-JRN-2016`), by appending either ``-SRK`` or ``-PengRobinson`` to the fluid name.  For instance, you could calculate the NBP with both the conventional multi-fluid model (as in GERG), or translating both to cubic equations of state

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # With the normal multi-fluid model
    In [0]: CP.PropsSI('T', 'P', 101325, 'Q', 0, 'HEOS::Methane-SRK[0.4]&Ethane-SRK[0.6]')

    # With both fluids using cubic translations in the multi-fluid model
    In [0]: CP.PropsSI('T', 'P', 101325, 'Q', 0, 'HEOS::Methane[0.4]&Ethane[0.6]')

Detailed Example
----------------

Here we plot phase envelopes and critical points for an equimolar methane/ethane mixture

.. plot::

    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt

    # Increase the starting pressure a bit, behavior at very low pressure is problematic
    CP.set_config_double(CP.PHASE_ENVELOPE_STARTING_PRESSURE_PA, 1e4)

    SRK = CP.AbstractState('SRK','Methane&Ethane')
    SRK.set_mole_fractions([0.5, 1 - 0.5])
    for kij, c in zip([0.0, 0.1],['r','b']):
        
        # Set the interaction parameter
        SRK.set_binary_interaction_double(0, 1, "kij", kij)

        # Some VLE calculations
        for p in [1e5, 1e6]:
            SRK.update(CP.PQ_INPUTS, p, 0)
            plt.plot(SRK.T(), SRK.p(), '<', color = c)

            SRK.update(CP.PQ_INPUTS, p, 1)
            plt.plot(SRK.T(), SRK.p(), '>', color = c)

        # Phase envelope
        SRK.build_phase_envelope("")
        PE = SRK.get_phase_envelope_data()
        plt.plot(PE.T, PE.p, '-', label = '$k_{ij} = $' + str(kij), color = c)

        # Critical point
        pts = SRK.all_critical_points()
        for pt in pts:
          plt.plot(pt.T, pt.p, '*', color = c)

    # A phase envelope calculated with SRK transformations in a multi-fluid model
    HEOS = CP.AbstractState('HEOS','Methane-SRK&Ethane-SRK')
    HEOS.set_mole_fractions([0.5, 0.5])
    HEOS.build_phase_envelope("none")
    PE = HEOS.get_phase_envelope_data()
    plt.plot(PE.T, PE.p, '-', label = 'SRK with transformations in multi-fluid', color = 'g')

    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

Volume Translation
==================

Cubic equations of state predict saturated liquid densities poorly with average deviations of about 
16 % and 8 % for SRK and PR respectively :cite:`Jaubert-FPE-2016`.  The Péneloux correction :cite:`Peneloux-FPE-1982` 
fixes this at almost no cost by shifting the molar volume by a constant:

.. math::

  v = v^{\mathrm{EOS}} - c

Equivalently, the Helmholtz energy is evaluated at a shifted volume,
:math:`a^{vs}(T,v) = a^{\mathrm{EOS}}(T,v+c)`, which is how it is implemented in CoolProp.  For a
mixture the translation is composition dependent through the **linear** mixing rule

.. math::

  c_m(\mathbf{x}) = \sum_i x_i c_i

This rule is not a choice.  Privat et al. :cite:`Privat-FPE-2016` show that requiring the translation
to leave phase equilibria unchanged forces the arithmetic-mean combining rule for :math:`c_{ij}`,
which collapses to the linear form above -- so there is no :math:`c_{ij}` to fit or store.

What changes and what does not
------------------------------

A constant translation is *not* a fudge factor applied to one output; it changes some properties by
an exactly predictable amount and leaves others untouched.  The table below follows Jaubert et al.
:cite:`Jaubert-FPE-2016` and Privat et al. :cite:`Privat-FPE-2016` and is pinned by the test suite 
at round-off level.

=====================================================  =============================================
Quantity                                               Effect of the translation
=====================================================  =============================================
:math:`v`                                              :math:`v - c_m(\mathbf{z})`
:math:`h`, :math:`g`                                   shifted by :math:`-p\,c_m(\mathbf{z})`
:math:`\mu_i`                                          shifted by :math:`-p\,c_i`, using the
                                                       **pure-component** :math:`c_i`
:math:`\ln \varphi_i`                                  shifted by :math:`-p\,c_i/(RT)`, the same
                                                       pure-component :math:`c_i`
:math:`B` (second virial coefficient)                  :math:`B - c_m(\mathbf{z})`
:math:`s`, :math:`u`, :math:`a`, :math:`c_p`,          **unchanged**
:math:`c_v`
:math:`-(\partial v/\partial p)_T`,                    **unchanged**
:math:`(\partial v/\partial T)_p`
:math:`p^{\mathrm{sat}}`, :math:`p_{\mathrm{bub}}`,    **unchanged**
:math:`p_{\mathrm{dew}}`, :math:`K_i`, phase
compositions
:math:`w` (speed of sound)                             **changes** (:math:`w^2 \propto v^2`)
:math:`\Delta_{\mathrm{vap}} S`,                       **unchanged**
:math:`\Delta_{\mathrm{vap}} U`,
:math:`\Delta_{\mathrm{vap}} A`,
:math:`\Delta_{\mathrm{vap}} c_p`,
:math:`\Delta_{\mathrm{vap}} c_v`
:math:`\Delta_{\mathrm{vap}} V`,                       unchanged for a pure fluid; for a mixture
:math:`\Delta_{\mathrm{vap}} H`                        shifted by
                                                       :math:`-[c_m(\mathbf{y})-c_m(\mathbf{x})]`
                                                       and
                                                       :math:`-p[c_m(\mathbf{y})-c_m(\mathbf{x})]`
=====================================================  =============================================

Two consequences are worth spelling out.  Because :math:`\ln \varphi_i` shifts by the same amount in
both phases at the same :math:`(T,p)`, the shift cancels in the equifugacity condition: **VLE is
invariant for arbitrary, unequal** :math:`c_i`, not merely for equal ones.  And because the two
coexisting phases have different compositions, :math:`\Delta_{\mathrm{vap}} V` and
:math:`\Delta_{\mathrm{vap}} H` are the one place a mixture behaves differently from a pure fluid.

.. warning:: CoolProp implements a **temperature-independent** :math:`c` only.  A :math:`c(T)` would
             break most of the "unchanged" rows above -- it perturbs :math:`s`, :math:`u`,
             :math:`c_p` and :math:`c_v` as well, which is precisely why constant translations are
             the standard choice.

Setting the translation
-----------------------

:math:`c` is a per-component fluid parameter, in m\ :sup:`3`/mol, reached through
:cpapi:`set_fluid_parameter_double <CoolProp::AbstractState::set_fluid_parameter_double>` and
:cpapi:`get_fluid_parameter_double <CoolProp::AbstractState::get_fluid_parameter_double>` under the
names ``c``, ``cm`` or ``c_m``.  The spelling ``c_all`` broadcasts one value to every component.
Setting a value that would drive :math:`b_i - c_i \le 0` raises a ``ValueError``, since past that
point the repulsive pole moves through the liquid branch and the model returns finite nonsense
rather than failing loudly.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("PR", "n-Propane")

    In [0]: AS.set_fluid_parameter_double(0, "cm", -3.7e-6)

    In [0]: AS.get_fluid_parameter_double(0, "cm")


Adding your own fluids
======================

The cubic fluids in CoolProp are defined based on a JSON format, which could yield something like this for a **FAKE** (illustrative) fluid.  For instance if we had the file ``fake_fluid.json`` (download it here: :download:`fake_fluid.json`):

.. literalinclude:: fake_fluid.json

The JSON-formatted fluid information is validated against the `JSON schema <https://json-schema.org>`_, which can be obtained by calling the ``get_global_param_string`` function:

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP

    # Just the first bit of the schema (it's a bit large; see below for the complete schema)
    In [0]: CP.get_global_param_string("cubic_fluids_schema")[0:60]

A `JSON schema <https://json-schema.org>`_ defines the structure of the data, and places limits on the values, defines what values must be provided, etc..  There are a number of online tools that can be used to experiment with and validate JSON schemas, for instance https://www.jsonschemavalidator.net .

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP, json

    # Normally you would read it from a file!
    In [0]: fake_fluids = [
       ...:                   { 
       ...:                     "CAS": "000-0-00",
       ...:                     "Tc": 400.0,
       ...:                     "Tc_units": "K",
       ...:                     "acentric": 0.1,
       ...:                     "aliases": [
       ...:                     ],
       ...:                     "molemass": 0.04,
       ...:                     "molemass_units": "kg/mol", 
       ...:                     "name": "FAKEFLUID",
       ...:                     "pc": 4000000.0,
       ...:                     "pc_units": "Pa"
       ...:                   }
       ...:               ]

    # Adds fake fluid for both SRK and PR backends since they need the same parameters from the same internal library
    In [0]: CP.add_fluids_as_JSON("PR", json.dumps(fake_fluids))

    In [0]: CP.PropsSI("T","P",101325,"Q",0,"SRK::FAKEFLUID")

    # Put in a placeholder interaction parameter (all beta and gamma values are 1.0)
    In [0]: CP.apply_simple_mixing_rule("000-0-00", CP.get_fluid_param_string("Ethane","CAS"), "Lorentz-Berthelot")

    # Once a fluid has been added to the cubic library, it can be used in the multi-fluid model
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"HEOS::FAKEFLUID-SRK[0.3]&Ethane-SRK[0.7]")

Complete fluid schema
---------------------

For completeness, here is the whole `JSON schema <https://json-schema.org>`_ used for the cubic backends:

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP

    In [0]: print(CP.get_global_param_string("cubic_fluids_schema"))