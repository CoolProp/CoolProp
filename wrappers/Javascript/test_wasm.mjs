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

    // Round-trip the composition through get_mole_fractions(). The
    // VectorDouble::get(i) → std::optional<double> binding is registered
    // explicitly in src/emscripten_interface.cxx; embind surfaces it as a
    // plain JS number when the value is present.
    var got = ASmix.get_mole_fractions();
    var got0 = got.get(0);
    var got1 = got.get(1);
    console.log("mixture composition round-trip:", got0, got1);
    if (Math.abs(got0 - 0.4) > 1e-12 || Math.abs(got1 - 0.6) > 1e-12) {
        console.error("get_mole_fractions did not round-trip composition");
        process.exit(1);
    }
    got.delete();

    z.delete();
    ASmix.delete();
    console.log("set_mole_fractions test passed!");
} catch (e) {
    console.error("set_mole_fractions test failed:", e);
    process.exit(1);
}

// Sanity-check that enum entries added in the binding-coverage follow-up
// resolve (catches accidental drops if these enums are edited later).
console.log("Testing enum coverage...");
try {
    var enumProbes = [
        ["parameters.iQmass",                  coolprop.parameters.iQmass],
        ["input_pairs.QmassT_INPUTS",          coolprop.input_pairs.QmassT_INPUTS],
        ["input_pairs.PQmass_INPUTS",          coolprop.input_pairs.PQmass_INPUTS],
        ["input_pairs.DmolarQmass_INPUTS",     coolprop.input_pairs.DmolarQmass_INPUTS],
        ["backend_families.INVALID_BACKEND_FAMILY", coolprop.backend_families.INVALID_BACKEND_FAMILY],
        ["backend_families.SVDSBTL_BACKEND_FAMILY", coolprop.backend_families.SVDSBTL_BACKEND_FAMILY],
    ];
    for (var i = 0; i < enumProbes.length; i++) {
        var name = enumProbes[i][0], v = enumProbes[i][1];
        if (v === undefined || v === null || !Number.isFinite(v.value)) {
            console.error("enum entry missing or non-numeric:", name, "=", v);
            process.exit(1);
        }
        console.log("  " + name + " =", v.value);
    }
    console.log("Enum coverage test passed!");
} catch (e) {
    console.error("Enum coverage test failed:", e);
    process.exit(1);
}

// Test the build_phase_envelope -> get_phase_envelope_data round-trip.
// The value_object<PhaseEnvelopeData> binding (X-macro driven) is the
// highest-risk surface from PR #2664 and was not covered by its tests.
console.log("Testing phase envelope round-trip...");
try {
    var ASenv = coolprop.factory("HEOS", "Methane&Ethane");
    var ze = new coolprop.VectorDouble();
    ze.push_back(0.5);
    ze.push_back(0.5);
    ASenv.set_mole_fractions(ze);
    ASenv.build_phase_envelope("");
    var env = ASenv.get_phase_envelope_data();

    // env.T and env.p are VectorDouble fields populated by the envelope tracer.
    var nT = env.T.size();
    var np = env.p.size();
    console.log("envelope points: T=" + nT + ", p=" + np);
    if (nT < 5 || np !== nT) {
        console.error("envelope returned too few points or T/p size mismatch:", nT, np);
        process.exit(1);
    }
    // Sample the middle point; both must be finite and in plausible ranges
    // for Methane&Ethane. Exercises VectorDouble::get(i) -> optional<double>.
    var Tmid = env.T.get(Math.floor(nT/2));
    var pmid = env.p.get(Math.floor(np/2));
    console.log("envelope mid-point: T=" + Tmid + " K, p=" + pmid + " Pa");
    if (!(Number.isFinite(Tmid) && Tmid > 50 && Tmid < 400)) {
        console.error("envelope mid-T out of plausible range:", Tmid);
        process.exit(1);
    }
    if (!(Number.isFinite(pmid) && pmid > 1 && pmid < 1e8)) {
        console.error("envelope mid-p out of plausible range:", pmid);
        process.exit(1);
    }

    ze.delete();
    ASenv.delete();
    console.log("Phase envelope test passed!");
} catch (e) {
    console.error("Phase envelope test failed:", e);
    process.exit(1);
}

// Test derivative bindings. The point is to exercise each binding's signature
// end-to-end against real backend code; tolerances are loose.
console.log("Testing derivative bindings...");
try {
    var ASd = coolprop.factory("HEOS", "Water");

    // first_partial_deriv: dH/dT|P at single-phase point should equal Cp_molar.
    ASd.update(coolprop.input_pairs.PT_INPUTS, 1e6, 300);
    var dHdT = ASd.first_partial_deriv(coolprop.parameters.iHmolar,
                                       coolprop.parameters.iT,
                                       coolprop.parameters.iP);
    var cp = ASd.cpmolar();
    console.log("dH/dT|P:", dHdT, " cp_molar:", cp);
    if (!Number.isFinite(dHdT) || Math.abs(dHdT - cp) > 1e-6 * Math.abs(cp)) {
        console.error("first_partial_deriv(H,T,P) should equal cpmolar");
        process.exit(1);
    }

    // first_saturation_deriv: dp/dT along the saturation line, Clausius-Clapeyron > 0.
    ASd.update(coolprop.input_pairs.PQ_INPUTS, 101325, 0);
    var dpsatdT = ASd.first_saturation_deriv(coolprop.parameters.iP,
                                             coolprop.parameters.iT);
    console.log("dp/dT|sat at NBP:", dpsatdT);
    if (!(Number.isFinite(dpsatdT) && dpsatdT > 0)) {
        console.error("first_saturation_deriv(P,T) should be positive at the NBP");
        process.exit(1);
    }

    // first_two_phase_deriv_splined: smooth two-phase derivative; just check finite.
    // x_end is the upper bound of the spline region; backend requires Q < x_end.
    ASd.update(coolprop.input_pairs.PQ_INPUTS, 1e5, 0.05);
    var dDdh_spl = ASd.first_two_phase_deriv_splined(coolprop.parameters.iDmolar,
                                                    coolprop.parameters.iHmolar,
                                                    coolprop.parameters.iP,
                                                    0.1);
    console.log("d rho_molar / dh |P (splined, Q=0.05, x_end=0.1):", dDdh_spl);
    if (!Number.isFinite(dDdh_spl)) {
        console.error("first_two_phase_deriv_splined returned non-finite");
        process.exit(1);
    }

    ASd.delete();
    console.log("Derivative bindings test passed!");
} catch (e) {
    console.error("Derivative bindings test failed:", e);
    process.exit(1);
}
