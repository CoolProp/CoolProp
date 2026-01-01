import Module from './coolprop.js'
var coolprop = await Module();

// Test basic functions
console.log('32 F in K:', coolprop.F2K(32.0));
console.log('NBP of water in K:', coolprop.PropsSI('T','P',101325,'Q',0,'Water'));

var Tnbp = coolprop.PropsSI('T','P',101325,'Q',0,'Water');
if (Tnbp < 373.124 || Tnbp > 373.125){
  console.error("NBP is invalid");
  process.exit(1);
}

// Test AbstractState
console.log("Testing AbstractState...");
try {
    // Correct way to use the factory function exposed globally or on Module
    var AS = coolprop.factory("HEOS", "Water");
    console.log("Backend name:", AS.backend_name());

    // Test update
    console.log("Updating state to P=101325 Pa, Q=0");
    AS.update(coolprop.input_pairs.PQ_INPUTS, 101325, 0);

    // Test property accessors
    var T = AS.T();
    console.log("T:", T);
    if (Math.abs(T - 373.124) > 0.1) {
        console.error("AbstractState T is invalid");
        process.exit(1);
    }

    console.log("P:", AS.p());
    console.log("rhomolar:", AS.rhomolar());
    console.log("rhomass:", AS.rhomass());
    console.log("hmolar:", AS.hmolar());
    console.log("smolar:", AS.smolar());
    console.log("cpmolar:", AS.cpmolar());
    console.log("viscosity:", AS.viscosity());
    console.log("conductivity:", AS.conductivity());

    // Test parameters enum usage
    console.log("T (via keyed_output):", AS.keyed_output(coolprop.parameters.iT));
    console.log("P (via keyed_output):", AS.keyed_output(coolprop.parameters.iP));

    // Test phases enum
    console.log("Phase (enum value):", coolprop.phases.iphase_twophase.value);

    // Test specific functions
    console.log("Critical T:", AS.T_critical());
    console.log("Gas constant:", AS.gas_constant());
    console.log("Molar mass:", AS.molar_mass());

    // Clean up if necessary (JS GC handles it mostly, but for C++ objects...)
    AS.delete();

    console.log("AbstractState tests passed!");

} catch (e) {
    console.error("AbstractState test failed:", e);
    process.exit(1);
}
