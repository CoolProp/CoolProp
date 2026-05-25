.. _backend-options:

***************
Backend Options
***************

Per-instance backend configuration is passed through the factory string
as a JSON suffix.  The mechanism is opt-in per backend, immutable after
construction, available identically from every language that wraps
CoolProp, and works through both the :ref:`high-level <high_level_api>`
and :ref:`low-level <low_level_api>` interfaces.

.. _backend-options-grammar:

Grammar
=======

Append ``?<options>`` to the factory string:

.. code-block:: text

   SVDSBTL&HEOS::Water                                  ← no options, all schema defaults
   SVDSBTL&HEOS::Water?{"critical_patch":"off"}         ← inline JSON
   SVDSBTL&HEOS::Water?@/path/to/cfg.json               ← file indirection
   HEOS::Water?{}                                       ← empty options (no-op)

The parser splits on the **first** ``?`` and treats the entire remainder
as a single JSON value (or, if it begins with ``@``, as a file path
read verbatim).  Any subsequent ``?`` characters — inside a JSON string
value, in a URL embedded in a value, in a file path after ``@`` — are
preserved as-is.

Options live in either the backend or the fluid half of the factory
string (whichever your call site finds convenient).  If both halves
carry an options suffix, factory construction throws a ``ValueError``
to flag the typo.

.. _backend-options-highlevel:

Use from the high-level interfaces
==================================

The factory string is parsed in exactly one place inside CoolProp
(``AbstractState::factory``).  Every wrapper that takes a fluid name
forwards it verbatim, so the ``?<options>`` syntax works from any
language with no per-wrapper changes:

**Python**

.. code-block:: python

   from CoolProp.CoolProp import PropsSI
   PropsSI("D", "T", 300, "P", 1e5,
           'SVDSBTL&HEOS::Water?{"critical_patch":"off"}')

**FORTRAN**

.. code-block:: fortran

   d = PropsSI('D'//c_null_char, 'T'//c_null_char, 300.0_dp, &
               'P'//c_null_char, 1.0e5_dp, &
               'SVDSBTL&HEOS::Water?@/opt/coolprop/h2o.json'//c_null_char)

**Excel / LibreOffice Calc**

.. code-block:: text

   =PropsSI("D";"T";300;"P";1E5;"SVDSBTL&HEOS::Water?@C:\Configs\h2o.json")

**MATLAB**

.. code-block:: matlab

   d = py.CoolProp.CoolProp.PropsSI('D','T',300,'P',1e5, ...
           'SVDSBTL&HEOS::Water?@~/coolprop/opts.json');

For Excel / FORTRAN / similar interfaces where embedding inline JSON
gets painful (quote escaping, cell length limits, etc.), the
``?@<path>`` form is the recommended path — write the JSON once to a
file, reference the path on every call.

.. _backend-options-schema:

Strict validation
=================

Every options-aware backend ships a JSON Schema (Draft 2020-12) that
the factory validates against on every call:

* Unknown keys throw — typos surface immediately rather than silently
  defaulting.
* Type mismatches throw — strings vs numbers vs booleans are checked
  per the schema.
* Required keys missing throw with the dotted path of the missing field.
* Enum violations throw with the offending value plus the allowed set.

Schemas live alongside the headers in
``include/CoolProp/schemas/<backend>.schema.json``.  Each opted-in
backend exposes its schema as a constant string literal compiled into
the binary — no runtime file lookups, no version-skew between docs and
runtime.

The schema is itself versioned (``"schema": 1`` at the top level) so
the layout can evolve with explicit migration support.

.. _backend-options-canonical:

Reproducibility — ``build_options_json()`` and cache keys
=========================================================

Two callers need a single canonical form per logically-equal options
value: cache filename hashing (so two callers passing the same options
in different key orders hit the same cache file) and reproducibility
logging.  The shared helper:

* Sorts object keys recursively.
* Preserves array order.
* Uses RapidJSON's default formatting for scalars (deterministic
  within a single CoolProp build).

It is **not** strict `RFC 8785 (JCS)
<https://datatracker.ietf.org/doc/html/rfc8785>`_ — no NFC string
normalisation, no ECMAScript number rounding — but it's
deterministic within a single CoolProp build, which is all the cache
key needs.

Every ``AbstractState`` exposes:

.. code-block:: cpp

   virtual std::string build_options_json() const;

Returns the canonical-JSON string the instance was built with
(defaults expanded).  Default is the empty string; backends that
accept options override to return the canonical form so callers can
copy-paste it into a fresh ``factory()`` call and reproduce the
construction byte-for-byte.

For caching backends (SVDSBTL today), the cache filename incorporates
a 16-hex-character FNV-1a 64 prefix of the canonical bytes alongside
fluid name, source backend, and symbolic input-pair name.  FNV-1a 64
is the same hash family CoolProp already uses for the superancillary
``source_eos_hash`` freshness stamp — no new dependencies, identical
determinism guarantees across compilers.

.. _backend-options-opting-in:

Opting a backend in
===================

Backends that want to consume options publish a JSON Schema and
override the options-aware ``AbstractStateGenerator`` virtual:

.. code-block:: cpp

   class MyBackendGenerator : public AbstractStateGenerator {
    public:
       AbstractState* get_AbstractState(
           const std::vector<std::string>& fluid_names,
           const std::string& options_json) override {
           rapidjson::Document opts;
           if (!options_json.empty()) opts.Parse(options_json);
           else                       opts.SetObject();
           validate_against_schema(opts, kMyBackendOptionsSchemaJson);
           // ...construct using the parsed options...
       }
       AbstractState* get_AbstractState(
           const std::vector<std::string>& fluid_names) override {
           return get_AbstractState(fluid_names, "");
       }
   };

The default implementation of the options-aware overload forwards to
the no-options overload when ``options_json`` is empty / ``"{}"`` and
throws ``NotImplementedError`` otherwise — so backends that haven't
opted in reject options loudly instead of silently dropping them.

.. _backend-options-migration:

Migration from Configuration globals
====================================

The ``Configuration`` mechanism is process-global and stringly-typed.
Several backend knobs that are currently Configuration globals are
better expressed as per-instance options.  Over time the following
mappings will land as deprecations of the corresponding globals:

.. list-table::
   :header-rows: 1
   :widths: 28 32 40

   * - Backend
     - Today (Configuration global)
     - After migration (options key)
   * - SVDSBTL
     - *(none — new backend)*
     - ``grid``, ``properties``, ``critical_patch``
   * - BICUBIC
     - ``BICUBIC_GRID_NT``, ``BICUBIC_GRID_NR``
     - ``grid.NT``, ``grid.NR``
   * - TTSE
     - ``TTSE_GRID_NT``, ``TTSE_GRID_NR``
     - ``grid.NT``, ``grid.NR``
   * - REFPROP
     - ``ALTERNATIVE_REFPROP_LIBRARY_PATH``,
       ``ALTERNATIVE_REFPROP_HMX_BNC_PATH``,
       ``ALTERNATIVE_REFPROP_PATH``
     - ``library_path``, ``hmx_path``, ``root_path``

Each follow-up PR ships its backend's schema, an opt-in for the
options-aware overload, and a deprecation warning emitted whenever
the corresponding Configuration global is set.  Existing scripts that
rely on the Configuration globals keep working through at least one
release cycle.

.. _backend-options-rationale:

Design rationale
================

A short tour of why the design lands where it does:

* **Why not mutators?**  ``set_<knob>()`` methods make caching brittle
  (cache key would have to reflect every post-construction setting)
  and pessimise reproducibility.  The factory-string form is the
  single source of truth, baked into the cache key, and round-trips
  through ``build_options_json()``.
* **Why not Configuration globals?**  Process-global state forbids
  mixing two backends with different policies in the same process,
  silently inherits across unrelated code, and rides outside the
  cache key.  The factory-string form is per-instance and explicit.
* **Why not a builder API?**  Wrapper bridges (SWIG, Cython,
  Mathematica) would each need per-method bindings.  The string-only
  factory entry-point requires zero new wrapper code.
* **Why JSON over query-string** (``?critical_patch=off&NT=200``)?
  JSON gives typed nested structures (``critical_patch.tolerance``
  next to ``critical_patch.bbox`` without flat-key collisions), and
  the schema evolution story is straightforward.  The high-level
  ``?@path.json`` indirection sidesteps the inline-quoting tax.
* **Why split on the first** ``?`` **only?**  Any other rule re-introduces
  parser ambiguity for URLs and regexes inside JSON string values.
  JSON's quoting rules are the source of truth for everything after
  the first ``?``.

See the design document
``docs/superpowers/specs/2026-05-16-backend-options-string-design.md``
for the full decision trail.
