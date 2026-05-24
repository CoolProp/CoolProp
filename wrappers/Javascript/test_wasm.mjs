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

// Test mixture composition via set_mole_fractions(VectorDouble).
// The set is verified by cross-checking the resulting density against
// PropsSI with an explicit-composition fluid string; a bit-exact match
// confirms the VectorDouble was correctly threaded into AbstractState.
console.log("Testing set_mole_fractions on a binary mixture...");
try {
    var z = new coolprop.VectorDouble();
    z.push_back(0.4);
    z.push_back(0.6);

    var ASmix = coolprop.factory("HEOS", "Methane&Ethane");
    if (!ASmix.using_mole_fractions()) {
        console.error("HEOS mixture should report using_mole_fractions() === true");
        process.exit(1);
    }
    ASmix.set_mole_fractions(z);
    ASmix.update(coolprop.input_pairs.PT_INPUTS, 1e6, 250);

    var rho = ASmix.rhomass();
    var rhoRef = coolprop.PropsSI('D', 'P', 1e6, 'T', 250, 'HEOS::Methane[0.4]&Ethane[0.6]');
    console.log("mixture rhomass (set_mole_fractions):", rho);
    console.log("mixture rhomass (PropsSI ref):       ", rhoRef);
    if (Math.abs(rho - rhoRef) > 1e-8 * Math.abs(rhoRef)) {
        console.error("set_mole_fractions did not produce a state matching PropsSI reference");
        process.exit(1);
    }

    z.delete();
    ASmix.delete();
    console.log("set_mole_fractions test passed!");
} catch (e) {
    console.error("set_mole_fractions test failed:", e);
    process.exit(1);
}
