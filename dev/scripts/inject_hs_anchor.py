import json, CoolProp, sys

CP = CoolProp.CoolProp

def inject(fluid):
    
    fluid_path = '../fluids/'+fluid+'.json'
    
    # Open the fluid JSON file
    fp = open(fluid_path, 'r')
    j = json.load(fp)
    fp.close()
    
    Tanchor = 1.1*CoolProp.CoolProp.PropsSI('Tcrit','D2O')
    rhoanchor = 0.9*CoolProp.CoolProp.PropsSI('rhomolar_critical','D2O')
    hanchor_molar = CP.PropsSI('Hmolar','T',Tanchor,'D',rhoanchor,fluid)
    sanchor_molar = CP.PropsSI('Smolar','T',Tanchor,'D',rhoanchor,fluid)
    p = CP.PropsSI('P','T',Tanchor,'D',rhoanchor,fluid)
    
    j['EOS'][0]['STATES']['hs_anchor'] = {
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
    
    sys.path.append('..')
    from package_json import json_options
    fp = open(fluid_path, 'w')
    fp.write(json.dumps(j, **json_options))
    fp.close()
    
    print('writing '+ fluid)

inject('HeavyWater')