.. _Rust:

************
Rust Wrapper
************

.. contents:: :depth: 2

rfluids (3-party wrapper)
=========================

``rfluids`` provides a safe, idiomatic Rust interface to CoolProp.

The crate is published on `crates.io <https://crates.io/crates/rfluids>`_ and
the API documentation is available on `docs.rs <https://docs.rs/rfluids>`_.

Installation
------------

Add ``rfluids`` to your ``Cargo.toml``: ::

    [dependencies]
    rfluids = "<version>"

Or use ``cargo add``: ::

    cargo add rfluids

Supported Platforms
-------------------

* Linux AArch64
* Linux x86-64
* macOS AArch64
* macOS x86-64
* Windows AArch64
* Windows x86-64

Benefits
--------

* **Idiomatic Rust interface:** ergonomic high-level abstractions, plus direct access to CoolProp's high-level and low-level APIs through the ``native`` module
* **Type-safety:** leverages Rust's type system with the `Typestate Pattern <https://en.wikipedia.org/wiki/Typestate_analysis>`_ to prevent invalid operations at compile time
* **Configuration management:** flexible control over CoolProp configuration via type-safe builder, with optional ``serde`` support for loading from configuration files
* **Batteries included:** pre-compiled CoolProp dynamic libraries for all supported platforms
* **Binding regeneration:** optional ``regen-bindings`` feature for regenerating FFI bindings when needed
* **Comprehensive documentation:** API docs include practical examples for common fluid, mixture, humid air, configuration, and native CoolProp use cases

Examples
--------

All calculations are performed in SI units.

Specific heat **[J/kg/K]** of saturated water vapor at *1 atm*: ::

    use approx::assert_relative_eq;
    use rfluids::prelude::*;

    let mut water_vapor = Fluid::from(Pure::Water)
        .in_state(FluidInput::pressure(101_325.0), FluidInput::quality(1.0))?;
    assert_relative_eq!(water_vapor.specific_heat()?, 2_079.937_085_633_241, max_relative = 1e-6);

Dynamic viscosity **[Pa·s]** of propylene glycol aqueous solution
with *60 %* mass fraction at *100 kPa* and *-20 °C*: ::

    use approx::assert_relative_eq;
    use rfluids::prelude::*;

    let mut propylene_glycol = Fluid::from(BinaryMixKind::MPG.with_fraction(0.6)?)
        .in_state(FluidInput::pressure(100e3), FluidInput::temperature(253.15))?;
    assert_relative_eq!(
        propylene_glycol.dynamic_viscosity()?,
        0.139_073_910_539_388_78,
        max_relative = 1e-6
    );

Density **[kg/m³]** of ethanol aqueous solution
(with ethanol *40 %* mass fraction) at *200 kPa* and *4 °C*: ::

    use approx::assert_relative_eq;
    use rfluids::prelude::*;

    let mut mix =
        Fluid::try_from(CustomMix::mass_based([(Pure::Water, 0.6), (Pure::Ethanol, 0.4)])?)?
            .in_state(FluidInput::pressure(200e3), FluidInput::temperature(277.15))?;
    assert_relative_eq!(mix.density()?, 883.392_277_162_775_9, max_relative = 1e-6);

Wet-bulb temperature **[K]** of humid air at *300 m* above sea level,
*30 °C* and *50 %* relative humidity: ::

    use approx::assert_relative_eq;
    use rfluids::prelude::*;

    let mut humid_air = HumidAir::new().in_state(
        HumidAirInput::altitude(300.0)?,
        HumidAirInput::temperature(303.15),
        HumidAirInput::rel_humidity(0.5),
    )?;
    assert_relative_eq!(
        humid_air.wet_bulb_temperature()?,
        295.067_569_033_474_57,
        max_relative = 1e-6
    );

For questions, more examples, and source code, see
`rfluids on GitHub <https://github.com/portyanikhin/rfluids>`_.
