import CoolProp
masses = {'Cl': 35.45, 'H': 1.008, 'C': 12.011, 'O': 15.999, 'N': 14.007, 'Ar': 39.948, 'S': 32.06, 'Si': 28.085, 'F': 18.998, 'D': 2.01410, 'He': 4.0026, 'Kr': 83.798, 'Ne': 20.180, 'I': 126.90, 'Xe': 131.29}
for fluid in CoolProp.__fluids__:
    formula = CoolProp.CoolProp.get_fluid_param_string(fluid, "formula")
    if formula == 'N/A': continue
    bits = [_ for _ in formula.split('}') if _]
    m = 0
    for b in bits:
        el, c = b.replace('_', '').split('{')
        m += masses[el] * int(c) / 1000.0
    err = m / CoolProp.CoolProp.PropsSI(fluid, "M") - 1
    if abs(err) > 1e-3:
        print("%s %s" % (fluid, m / CoolProp.CoolProp.PropsSI(fluid, "M") - 1))
