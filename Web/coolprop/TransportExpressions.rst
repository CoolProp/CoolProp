.. _transport_expressions:

**************************************
Transport-Property Expression Language
**************************************

A viscosity or thermal-conductivity correlation can be written directly in a fluid's JSON
file as a formula, instead of being implemented in C++.  A block with
``"type": "expression"`` is compiled once when the fluid loads and evaluated at every
call.  Many of the reference correlations shipped with CoolProp, among them xenon,
krypton, nitrogen, ammonia and the refrigerants R-32 and R-245fa, are written this way,
and their fluid files in ``dev/fluids/`` are working examples.

Where a block can go
====================

Expression blocks are accepted in five places in a pure or pseudo-pure fluid's
``TRANSPORT`` section:

=========================================  =========================================
JSON location                              Contribution
=========================================  =========================================
``viscosity.dilute``                       dilute-gas viscosity, Pa s
``viscosity.initial_density``              initial-density term, Pa s
``viscosity.higher_order``                 residual viscosity, Pa s
``conductivity.dilute``                    dilute-gas conductivity, W/(m K)
``conductivity.residual``                  residual conductivity, W/(m K)
=========================================  =========================================

Each block returns its own contribution in base SI units, and CoolProp adds the
contributions together (for conductivity, with the critical enhancement, which is not
available as an expression).  In particular the ``initial_density`` block returns a viscosity, not the
Rainwater-Friend coefficient :math:`B_\eta`; if a correlation is written in terms of
:math:`B_\eta`, the block has to multiply by the dilute viscosity and the density
itself.  Expression blocks can sit alongside the built-in correlation types; the
xenon conductivity, for example, uses expressions for the dilute and residual terms and
the built-in ``simplified_Olchowy_Sengers`` type for the critical enhancement.

When ``viscosity`` or ``conductivity`` is a list of models, only the first is loaded, so
an expression block in a later entry is never compiled and its errors are never
reported.  Try such a block from Python (see `Trying a block from Python`_) before
relying on it.

A block also runs when the fluid is a mixture component, evaluated for that component
at the mixture's temperature and density, and when the fluid is the reference fluid of
another fluid's extended-corresponding-states model.

Block structure
===============

.. code-block:: none

    "residual": {
        "type": "expression",
        "note": "Eq. (3) and Table 3 of Velliadou et al. (2021) ...",
        "formula": "let Tr = T/Tc\nlet rhor = Dmolar/rhoc\nsum(i: (B1[i] + B2[i]*Tr)*rhor^d[i])",
        "state_variables": ["T", "Dmolar"],
        "constants": {"Tc": 289.733, "rhoc": 8400.0},
        "arrays": {
            "B1": [0.00694552, 0.00876111, -0.01199, 0.00684476, -0.00102229],
            "B2": [-7.32747e-05, -0.00268366, 0.00563598, -0.00314076, 0.000605394],
            "d": [1, 2, 3, 4, 5]
        }
    }

This is the residual thermal conductivity of xenon, as shipped.

``formula``
    The expression.  See `Formula syntax`_.
``state_variables``
    The thermodynamic quantities the formula reads from the current state.  See
    `State variables`_.
``constants``
    Named scalars, in SI units.
``arrays``
    Named coefficient vectors, which can only be read inside a ``sum``.
``note``, ``BibTeX``
    Ignored, like any other key apart from ``type`` (and ``hardcoded``, which the
    loader checks before ``type``).  The formula has no comment syntax,
    so ``note`` is where to record the source, unit conversions, and any departure
    from the equation as printed.

Formula syntax
==============

A formula is zero or more ``let`` statements followed by one final expression, whose
value is the block's result.  Statements are separated by newlines (``\n`` inside the
JSON string) or by semicolons.  A newline always ends a statement, even inside
parentheses, so one statement cannot be wrapped over several lines; use ``let`` to break
a long formula up instead.  Column numbers in error messages count from the start of
the whole formula, not from the start of the line.

.. code-block:: none

    let L = ln(T/T_ref)
    1e-3*lambda_ref*exp(sum(i: b[i]*L^p[i]))

**Arithmetic.**  ``+ - * / ^`` and parentheses.  ``^`` binds tightest and is
right-associative, so ``2^3^2`` is 512, and unary minus binds looser than ``^``, so
``-2^2`` is -4.  Numbers are decimal literals with an optional exponent (``8.4e3``, ``.5``).
The decimal separator is always ``.``, whatever the host program's locale; hexadecimal
literals and literals outside the range of a double are errors.

**Functions.**  ``exp``, ``ln``, ``log10``, ``sqrt``, ``abs``, ``sinh``, ``cosh``,
``tanh``, ``sin``, ``cos``, ``atan``, each of one argument, and ``pow(x, y)``.
They follow the C++ standard library at the edges of their domains, so ``ln(0)``
gives ``-inf`` and ``sqrt(-1)`` gives ``nan``.  Nothing is guarded, so that a
formula reproduces a hard-coded C++ correlation bit for bit.

**Sums.**  ``sum(i: body)`` sums ``body`` over the coefficient arrays.  Inside the
body, ``a[i]`` is the ``i``-th element of array ``a``, and bare ``i`` is the index
itself, counted from 0.  The number of terms is the length of the arrays subscripted
in the body, which must all be the same length.  Arrays can be subscripted only by
the index of the enclosing sum, and sums cannot be nested.  The body has to subscript at
least one non-empty array.  The index can have any name, but inside the body it hides
any ``let``, state variable or constant of the same name, so use one the formula does
not otherwise need.

**Names.**  Identifiers start with a letter or underscore and are case-sensitive.
A name is resolved first as a ``let``, then as a declared state variable, then as a
constant.  A ``let`` may not reuse the name of a declared state variable.  It may
reuse the name of a constant, which then hides that constant from the rest of the
formula, so it is best avoided.  A ``let`` can rebind an earlier ``let``
(``let x = 1; let x = x + 1``).  ``let`` and ``sum`` are reserved words; function
names are recognized only when called, so a constant may be named ``exp``.

State variables
===============

A formula can only read the quantities its block lists in ``state_variables``, using
CoolProp's own parameter names:

======================  ============================================================
Name                    Quantity
======================  ============================================================
``T``                   temperature, K
``P``                   pressure, Pa
``Dmolar``              molar density, mol/m\ :sup:`3`
``Dmass``               mass density, kg/m\ :sup:`3`
``molar_mass``          molar mass, kg/mol
``Smolar_residual``     residual molar entropy, J/(mol K)
``Bvirial``             second virial coefficient, m\ :sup:`3`/mol
``dBvirial_dT``         its temperature derivative, m\ :sup:`3`/(mol K)
======================  ============================================================

The declaration exists so that a block claims only the names it uses.  A correlation
that never needs pressure is free to call its exponent array ``p``, which is what most
viscosity papers call it.  The list is checked when the fluid loads:

* a name the formula reads without declaring it is an error, and the message says to
  add it to ``state_variables``;
* a declared name the formula never reads is an error, since each one costs a property
  evaluation at every call;
* a declared name that is also a constant or an array is an error.

Anything outside the table is refused, including CoolProp's aliases (``D``, ``S`` and
so on) and every other output.  Two refusals are deliberate.  Transport properties
would re-enter the correlation being defined.  The critical point and the EOS reducing
state depend on configuration (with superancillaries enabled, the critical point is the
numerical one), so a correlation's reducing constants belong in ``constants``, at the
values its authors regressed against.  The set is extended when a correlation needs a
new quantity.

A larger example
================

The krypton viscosity of Polychroniadou et al. (2022) is an entropy-scaling
correlation, so its residual block reads the residual entropy and the second virial
coefficient.  Because a block cannot read another block's result, it recomputes the
dilute viscosity ``eta0`` from the same coefficients, computes the total viscosity, and
subtracts ``eta0``, so that the ``dilute`` and ``higher_order`` contributions add up to
the paper's total:

.. code-block:: none

    let L = ln(T/T_ref)
    let eta0 = eta_ref*exp(sum(i: a[i]*L^np[i]))
    let NA = R/k_B
    let mass = M/NA
    let root_mkT = (mass*k_B*T)^0.5
    let Theta2 = (Bvirial + T*dBvirial_dT)/NA
    let eplus_dilute = eta0/root_mkT*Theta2^(2/3)
    let splus = -Smolar_residual/R
    let eplus_res = f*(exp(sum(i: d[i]*splus^dp[i])) - 1)
    (Dmolar*NA)^(2/3)*root_mkT*(eplus_res + eplus_dilute)/splus^(2/3) - eta0

with ``"state_variables": ["T", "Bvirial", "dBvirial_dT", "Smolar_residual",
"Dmolar"]``.  The complete block is in ``dev/fluids/Krypton.json``.

Trying a block from Python
==========================

``CoolProp.CoolProp.Expression`` compiles a block from the same JSON text that goes in
the fluid file and evaluates it at the state of any ``AbstractState``, so a correlation
can be checked against its paper without rebuilding CoolProp or editing a fluid file:

.. code-block:: python

    import json
    from CoolProp import CoolProp as CP

    block = {
        "formula": "let Tr = T/Tc\nlet rhor = Dmolar/rhoc\nsum(i: (B1[i] + B2[i]*Tr)*rhor^d[i])",
        "state_variables": ["T", "Dmolar"],
        "constants": {"Tc": 289.733, "rhoc": 8400.0},
        "arrays": {
            "B1": [0.00694552, 0.00876111, -0.01199, 0.00684476, -0.00102229],
            "B2": [-7.32747e-05, -0.00268366, 0.00563598, -0.00314076, 0.000605394],
            "d": [1, 2, 3, 4, 5],
        },
    }
    expr = CP.Expression(json.dumps(block))
    print(expr.required_inputs())       # ['T', 'Dmolar']

    AS = CP.AbstractState("HEOS", "Xenon")
    AS.update(CP.DmassT_INPUTS, 1200.0, 300.0)
    print(expr.evaluate(AS) * 1e3)      # 11.0620..., the paper's check value is 11.0621 mW/(m K)

``required_inputs()`` lists the state variables in the order the formula first reads
them.  Compilation errors raise ``ValueError`` with the same diagnostic the fluid loader
reports when it loads the fluid.  ``evaluate()`` raises if a state variable reads back as non-finite, which
usually means the ``AbstractState`` was never updated.  The C++ equivalent is
``CoolProp::expression::ExpressionBlock`` in
``include/CoolProp/expression/ExpressionBlock.h``.
