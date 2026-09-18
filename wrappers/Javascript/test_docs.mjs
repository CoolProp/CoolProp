// Extended tests for the CoolProp WebAssembly module, exercising the
// examples that appear in the public documentation:
//   - HighLevelAPI.rst   (PropsSI, imposed phase, Props1SI, derivatives,
//                         mixtures, incompressible fluids, IF97)
//   - examples.rst       (HAPropsSI humid-air examples)
//   - wrappers/Javascript README (npm package usage pattern)
// Run with:  node test_docs.mjs   (next to coolprop.js / coolprop.wasm)
import Module from './coolprop.js'
var coolprop = await Module();

var failures = 0;

function check(what, got, expected, relTol) {
    relTol = relTol || 1e-6;
    var ok = Number.isFinite(got) && Math.abs(got - expected) <= relTol * Math.abs(expected);
    if (ok) {
        console.log('  OK   ' + what + ' = ' + got);
    } else {
        console.error('  FAIL ' + what + ': got ' + got + ', expected ~' + expected);
        failures++;
    }
    return ok;
}

// ---------------------------------------------------------------------------
// HighLevelAPI.rst - "PropsSI function": enthalpies and latent heat of water
// at 1 atm. H_V - H_L must be the latent heat of vaporization (~2.257 MJ/kg).
// ---------------------------------------------------------------------------
console.log('Testing documented PropsSI examples (water at 1 atm)...');
try {
    var H_V = coolprop.PropsSI('H', 'P', 101325, 'Q', 1, 'Water');
    var H_L = coolprop.PropsSI('H', 'P', 101325, 'Q', 0, 'Water');
    console.log('  H_V =', H_V, ' H_L =', H_L);
    check('latent heat of vaporization of water at 1 atm', H_V - H_L, 2.2569e6, 1e-3);
} catch (e) {
    console.error('  FAIL PropsSI examples threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// HighLevelAPI.rst - "Imposing the Phase": phase-qualified input keys.
//   PropsSI('D','T|liquid',461.1,'P',5e6,'Water')
//   PropsSI('D','T',597.9,'P|gas',5e6,'Water')
// ---------------------------------------------------------------------------
console.log('Testing imposed phase syntax (T|liquid, P|gas)...');
try {
    var Dliq = coolprop.PropsSI('D', 'T|liquid', 461.1, 'P', 5e6, 'Water');
    console.log('  rho(T=461.1 K, p=5 MPa | liquid) =', Dliq);
    if (!(Number.isFinite(Dliq) && Dliq > 500 && Dliq < 1200)) {
        console.error('  FAIL imposed-liquid density out of plausible range');
        failures++;
    }
    var Dgas = coolprop.PropsSI('D', 'T', 597.9, 'P|gas', 5e6, 'Water');
    console.log('  rho(T=597.9 K, p=5 MPa | gas) =', Dgas);
    if (!(Number.isFinite(Dgas) && Dgas > 1 && Dgas < 100)) {
        console.error('  FAIL imposed-gas density out of plausible range');
        failures++;
    }
} catch (e) {
    console.error('  FAIL imposed phase threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// HighLevelAPI.rst - "Trivial inputs": Props1SI('Tcrit','Water') and the
// PropsSI dummy-argument equivalent, both must agree with T_critical().
// ---------------------------------------------------------------------------
console.log('Testing Props1SI / trivial inputs...');
try {
    var Tcrit1 = parseFloat(coolprop.Props1SI('Water', 'Tcrit'));
    var Tcrit2 = coolprop.PropsSI('Tcrit', '', 0, '', 0, 'Water');
    var pcrit = coolprop.PropsSI('pcrit', '', 0, '', 0, 'Water');
    console.log('  Tcrit =', Tcrit1, '/', Tcrit2, ' pcrit =', pcrit);
    check('Props1SI Tcrit of water', Tcrit1, 647.096, 1e-5);
    check('PropsSI Tcrit of water (dummy args)', Tcrit2, 647.096, 1e-5);
    check('PropsSI pcrit of water', pcrit, 22.064e6, 1e-5);

    var ASt = coolprop.factory('HEOS', 'Water');
    var TcritAS = ASt.T_critical();
    ASt.delete();
    check('AbstractState.T_critical agrees with Props1SI', TcritAS, Tcrit1, 1e-9);
} catch (e) {
    console.error('  FAIL Props1SI threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// HighLevelAPI.rst - "First Partial Derivatives": the string-encoded
// derivative d(Hmass)/d(T)|P must equal c_p evaluated directly.
// ---------------------------------------------------------------------------
console.log('Testing string-encoded partial derivatives...');
try {
    var cp_direct = coolprop.PropsSI('C', 'P', 101325, 'T', 300, 'Water');
    var cp_deriv = coolprop.PropsSI('d(Hmass)/d(T)|P', 'P', 101325, 'T', 300, 'Water');
    console.log('  C =', cp_direct, ' d(Hmass)/d(T)|P =', cp_deriv);
    check('c_p of water at 300 K/1 atm', cp_direct, 4180.6, 1e-3);
    check('d(Hmass)/d(T)|P equals C', cp_deriv, cp_direct, 1e-9);

    // Second derivative d(d(Hmass)/d(T)|P)/d(Hmass)|P must be finite
    var d2 = coolprop.PropsSI('d(d(Hmass)/d(T)|P)/d(Hmass)|P', 'P', 101325, 'T', 300, 'Water');
    console.log('  d(d(Hmass)/d(T)|P)/d(Hmass)|P =', d2);
    if (!Number.isFinite(d2)) {
        console.error('  FAIL second derivative is not finite');
        failures++;
    }
} catch (e) {
    console.error('  FAIL derivative examples threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// HighLevelAPI.rst - "Pre-defined mixtures" / "Mixtures with composition":
// Air.mix, R410A and an explicit HEOS mixture with composition.
// ---------------------------------------------------------------------------
console.log('Testing mixture examples from the docs...');
try {
    check('rho of Air.mix at 300 K/1 atm', coolprop.PropsSI('D', 'P', 101325, 'T', 300, 'Air.mix'), 1.1766, 1e-3);
    check('rho of R410A at 300 K/1 atm', coolprop.PropsSI('D', 'T', 300, 'P', 101325, 'R410A'), 2.9869, 1e-3);
    check('rho of HEOS::R32[0.697615]&R125[0.302385]',
          coolprop.PropsSI('D', 'T', 300, 'P', 101325, 'HEOS::R32[0.697615]&R125[0.302385]'),
          2.9869, 1e-3);
    // IF97 industrial formulation for water/steam
    check('IF97 density of steam at 400 K (Q=1)',
          coolprop.PropsSI('D', 'T', 400, 'Q', 1, 'IF97::Water'), 1.3697, 1e-2);
    // Incompressible mixture: aqueous ethylene glycol 20%
    check('cp of INCOMP::MEG-20% at 298.15 K',
          coolprop.PropsSI('C', 'T', 298.15, 'P', 101325, 'INCOMP::MEG-20%'), 3896, 1e-2);
} catch (e) {
    console.error('  FAIL mixture examples threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// examples.rst - "Sample HAPropsSI Code": humid air at STP, and the inverse
// solve for the saturation temperature at the same enthalpy (R=1).
// ---------------------------------------------------------------------------
console.log('Testing HAPropsSI examples (humid air)...');
try {
    var h = coolprop.HAPropsSI('H', 'T', 298.15, 'P', 101325, 'R', 0.5);
    console.log('  h(humid air, 25 C, 1 atm, R=0.5) =', h);
    check('humid-air enthalpy at STP, R=0.5', h, 50423, 1e-3);
    var Tsat = coolprop.HAPropsSI('T', 'P', 101325, 'H', h, 'R', 1.0);
    console.log('  T(saturated, same h) =', Tsat);
    if (!(Number.isFinite(Tsat) && Tsat > 273 && Tsat < 298.15)) {
        console.error('  FAIL saturation temperature out of plausible range:', Tsat);
        failures++;
    }
    // Dew point of the same state must be below the dry-bulb temperature
    var Tdp = coolprop.HAPropsSI('D', 'T', 298.15, 'P', 101325, 'R', 0.5);
    console.log('  dew point =', Tdp);
    if (!(Number.isFinite(Tdp) && Tdp < 298.15)) {
        console.error('  FAIL dew point should be below dry-bulb temperature');
        failures++;
    }
} catch (e) {
    console.error('  FAIL HAPropsSI examples threw:', e);
    failures++;
}

// ---------------------------------------------------------------------------
// npm package README usage pattern: default-export factory, high-level call,
// low-level AbstractState call.
// ---------------------------------------------------------------------------
console.log('Testing README usage pattern...');
try {
    check('F2K(32)', coolprop.F2K(32), 273.15, 1e-12);
    var AS = coolprop.factory('HEOS', 'Water');
    AS.update(coolprop.input_pairs.PQ_INPUTS, 101325, 0);
    check('AS.rhomass() at NBP', AS.rhomass(), 958.37, 1e-4);
    AS.delete();
} catch (e) {
    console.error('  FAIL README pattern threw:', e);
    failures++;
}

if (failures > 0) {
    console.error('\n' + failures + ' documentation example test(s) FAILED');
    process.exit(1);
}
console.log('\nAll documentation example tests passed!');
