import CoolProp, json
PropsSI = CoolProp.CoolProp.PropsSI
data = []
for fluid in CoolProp.__fluids__:
    if CoolProp.CoolProp.get_fluid_param_string(fluid,"pure") == 'true':
        data.append(dict(Tc = PropsSI('Tcrit', fluid),
                         Tc_units = "K",
                         pc = PropsSI('pcrit', fluid),
                         pc_units = "Pa",
                         acentric = PropsSI('acentric', fluid),
                         molemass = PropsSI('molemass', fluid),
                         molemass_units = "kg/mol",
                         name = fluid.upper(),
                         CAS = CoolProp.CoolProp.get_fluid_param_string(fluid,"CAS"),
                         aliases = list(set([_.upper() for _ in CoolProp.CoolProp.get_fluid_param_string(fluid,"aliases").split(', ')]))
                         ))

with open('all_cubic_fluids.json','w') as fp:
    json.dump(data, fp, indent = 2, sort_keys = True)
    