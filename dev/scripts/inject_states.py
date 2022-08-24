import sys, os
import json, CoolProp, sys

CP = CoolProp.CoolProp
CP.set_config_bool(CP.DONT_CHECK_PROPERTY_LIMITS, True)

def inject_hsanchor(fluid, i, json_data):

    Tanchor = 1.1 * CoolProp.CoolProp.PropsSI('Tcrit', fluid)
    rhoanchor = 0.9 * CoolProp.CoolProp.PropsSI('rhomolar_critical', fluid)
    hanchor_molar = CP.PropsSI('Hmolar', 'T', Tanchor, 'Dmolar', rhoanchor, fluid)
    sanchor_molar = CP.PropsSI('Smolar', 'T', Tanchor, 'Dmolar', rhoanchor, fluid)
    p = CP.PropsSI('P', 'T', Tanchor, 'Dmolar', rhoanchor, fluid)

    json_data['EOS'][i]['STATES']['hs_anchor'] = {
          "T": Tanchor,
          "T_units": "K",
          "hmolar": hanchor_molar,
          "hmolar_units": "J/mol",
          "p": p,
          "p_units": "Pa",
          "rhomolar": rhoanchor,
          "rhomolar_units": "mol/m^3",
          "smolar": sanchor_molar,
          "smolar_units": "J/mol/K"
    }


def inject_triples(fluid, i, json_data):

    # Ttriple = json_data['STATES']['triple_liquid']['T']
    Ttriple = json_data['EOS'][i]['STATES']['sat_min_vapor']['T']

    triple_liquid = {
          "T": Ttriple,
          "T_units": "K",
          "hmolar": CP.PropsSI('Hmolar', 'T', Ttriple, 'Q', 0, fluid),
          "hmolar_units": "J/mol",
          "p": CP.PropsSI('P', 'T', Ttriple, 'Q', 0, fluid),
          "p_units": "Pa",
          "rhomolar": CP.PropsSI('Dmolar', 'T', Ttriple, 'Q', 0, fluid),
          "rhomolar_units": "mol/m^3",
          "smolar": CP.PropsSI('Smolar', 'T', Ttriple, 'Q', 0, fluid),
          "smolar_units": "J/mol/K"
    }
    triple_vapor = {
          "T": Ttriple,
          "T_units": "K",
          "hmolar": CP.PropsSI('Hmolar', 'T', Ttriple, 'Q', 1, fluid),
          "hmolar_units": "J/mol",
          "p": CP.PropsSI('P', 'T', Ttriple, 'Q', 1, fluid),
          "p_units": "Pa",
          "rhomolar": CP.PropsSI('Dmolar', 'T', Ttriple, 'Q', 1, fluid),
          "rhomolar_units": "mol/m^3",
          "smolar": CP.PropsSI('Smolar', 'T', Ttriple, 'Q', 1, fluid),
          "smolar_units": "J/mol/K"
    }
    json_data['STATES']['triple_vapor'] = triple_vapor
    json_data['STATES']['triple_liquid'] = triple_liquid
    json_data['EOS'][i]['STATES']['sat_min_vapor'] = triple_vapor
    json_data['EOS'][i]['STATES']['sat_min_liquid'] = triple_liquid


def inject_critical(fluid, json_data):
    """
    Inject values for pressure, enthalpy, entropy for critical point
    """
    T = json_data['STATES']['critical']['T']
    D = json_data['STATES']['critical']['rhomolar']

    json_data['STATES']['critical'] = {
          "T": T,
          "T_units": "K",
          "hmolar": CP.PropsSI('Hmolar', 'T', T, 'Dmolar', D, fluid),
          "hmolar_units": "J/mol",
          "p": CP.PropsSI('P', 'T', T, 'Dmolar', D, fluid),
          "p_units": "Pa",
          "rhomolar": D,
          "rhomolar_units": "mol/m^3",
          "smolar": CP.PropsSI('Smolar', 'T', T, 'Dmolar', D, fluid),
          "smolar_units": "J/mol/K"
    }


def inject_reducing(fluid, i, json_data):
    """
    Inject values for pressure, enthalpy, entropy for reducing point
    """
    T = json_data['EOS'][i]['STATES']['reducing']['T']
    D = json_data['EOS'][i]['STATES']['reducing']['rhomolar']

    json_data['EOS'][i]['STATES']['reducing'] = {
          "T": T,
          "T_units": "K",
          "hmolar": CP.PropsSI('Hmolar', 'T', T, 'Dmolar', D, fluid),
          "hmolar_units": "J/mol",
          "p": CP.PropsSI('P', 'T', T, 'Dmolar', D, fluid),
          "p_units": "Pa",
          "rhomolar": D,
          "rhomolar_units": "mol/m^3",
          "smolar": CP.PropsSI('Smolar', 'T', T, 'Dmolar', D, fluid),
          "smolar_units": "J/mol/K"
    }


def inject_acentric(fluid, i, json_data):
    """
    Inject values for acentric factor for i-th EOS
    """
    Tc = json_data['STATES']['critical']['T']
    pc = json_data['STATES']['critical']['p']
    p = CoolProp.CoolProp.PropsSI('P', 'T', Tc * 0.7, "Q", 0, fluid)
    import math
    json_data['EOS'][i]['acentric'] = -math.log10(p / pc) - 1


def inject_states(fluid):
    fluid_path = '../fluids/' + fluid + '.json'

    # AS =

    # Open the fluid JSON file
    with open(fluid_path, 'r') as fp:
        json_data = json.load(fp)

    # Inject states
    inject_hsanchor(fluid, 0, json_data)
    inject_triples(fluid, 0, json_data)
    inject_critical(fluid, json_data)
    inject_reducing(fluid, 0, json_data)
    inject_acentric(fluid, 0, json_data)

    # Write back to file
    sys.path.append('..')
    from package_json import json_options
    with open(fluid_path, 'w') as fp:
        fp.write(json.dumps(json_data, **json_options))

if __name__ == '__main__':


    for fld in ['MD2M','MD3M','MD4M','R1243zf','R1234ze(Z)','Neon','HydrogenChloride','HeavyWater']:
      inject_states(fld)