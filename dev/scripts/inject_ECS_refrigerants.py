XXX = 'XXX'
import json

# From McLinden, IJR, 2000
McLinden_sigma_ek_f_data = {
'R11': [0.5447, 363.61, 1.4000e-3],
'R12': [0.5186, 297.24, 1.3440e-3],
'R13': [0.4909, 233.36, 1.3200e-3],
'R22': [0.4666, 284.72, 7.7817e-4, +1.2536e-6],
'R23': [0.4430, 230.83, 6.0570e-4, +1.8604e-6],
'R32': [0.4098, 289.65, 8.1980e-4, +2.2352e-7],
'R114': [0.5770, 323.26, 1.3200e-3],
'R115': [0.5476, 272.53, 1.3200e-3],
'R125': [0.5101, 261.39, 1.2565e-3, +2.2296e-6],
'R142b': [0.5320, 316.64, 1.3200e-3],
'R143a': [0.5025, 267.10, 1.0066e-3, +1.3729e-6]
}

# From McLinden, IJR, 2000
McLinden_psi_cond_data = {
'R11': [1.0724, -0.0226720],
'R12': [0.9910, 0.0029509],
'R13': [1.4078, -0.2634600, 0.037978],
'R22': [1.0750, -0.0385740],
'R23': [1.3801, -0.2797500, 0.048798],
'R32': [1.2325, -0.0883940],
'R114': [1.0961, -0.0348990],
'R115': [1.0338, -0.0020661],
'R125': [1.0369, -0.0030368],
'R142b': [1.6808, -0.8395440, 0.321957, -0.039706],
'R143a': [1.1779, -0.2054100, 0.064870, -0.006473]}

# From Klein, 1997
Klein_psi_visc_data = {
'R11': [1.171, -0.0592],
'R12': [1.087, -0.0379],
'R113': [1.212, -0.0569],
'R22': [1.106, -0.0491],
'R32': [0.898, +0.0099],
'R123': [1.155, -0.0513],
'R125': [1.110, -0.0332],
'R143a': [1.134, -0.0801],
'R152a': [0.807, +0.0496]}

template = {
"conductivity": {
      "BibTeX": "McLinden-IJR-2000",
      "f_int": {
        "T_reducing": 1.0,
        "T_reducing_units": "K",
        "a": [
          XXX
        ],
        "t": [
          XXX
        ]
      },
      "psi": {
        "a": [
          XXX
        ],
        "rhomolar_reducing": XXX,
        "rhomolar_reducing_units": "mol/m^3",
        "t": [
          XXX
        ]
      },
      "q_D": 1999999999.9999998,
      "q_D_units": "m",
      "reference_fluid": "R134a",
      "type": "ECS"
    },
    "viscosity": {
      "BibTeX": "Klein-IJR-1997",
      "epsilon_over_k": XXX,
      "epsilon_over_k_units": "K",
      "psi": {
        "a": [
          XXX
        ],
        "rhomolar_reducing": XXX,
        "rhomolar_reducing_units": "mol/m^3",
        "t": [
          XXX
        ]
      },
      "reference_fluid": "R134a",
      "sigma_eta": XXX,
      "sigma_eta_units": "m",
      "type": "ECS"
    }
}

import CoolProp

# Find all fluids with both a viscosity and conductivity curve
fluids = set(('R11,R12,R113,R22,R32,R123,R125,R143a,R152a'.split(','))).intersection(('R11,R12,R13,R22,R23,R32,R114,R115,R125,R142b,R143a'.split(',')))

for fluid in fluids:
    try:
        v = CoolProp.CoolProp.PropsSI('V', 'T', 300, 'Q', 0, fluid)
        l = CoolProp.CoolProp.PropsSI('L', 'T', 300, 'Q', 0, fluid)
        # print 'GOOD', fluid, v, l
    except ValueError:

        print('BAD %s' % fluid)

        new = template.copy()
        new['viscosity']['sigma_eta'] = McLinden_sigma_ek_f_data[fluid][0] / 1e9
        new['viscosity']['epsilon_over_k'] = McLinden_sigma_ek_f_data[fluid][1]
        a = Klein_psi_visc_data[fluid]
        new['viscosity']['psi']['a'] = a
        new['viscosity']['psi']['t'] = range(len(a))
        new['viscosity']['psi']['rhomolar_reducing'] = CoolProp.CoolProp.PropsSI(fluid, 'rhomolar_reducing')

        a = McLinden_sigma_ek_f_data[fluid][2::]
        new['conductivity']['f_int']['a'] = a
        new['conductivity']['f_int']['t'] = range(len(a))
        a = McLinden_psi_cond_data[fluid]
        new['conductivity']['psi']['a'] = a
        new['conductivity']['psi']['t'] = range(len(a))
        new['conductivity']['psi']['rhomolar_reducing'] = CoolProp.CoolProp.PropsSI(fluid, 'rhomolar_reducing')

        fname = '../fluids/' + fluid + '.json'
        with open(fname, 'r') as fp:
            jj = json.load(fp)
        jj['TRANSPORT'] = new
        with open(fname, 'w') as fp:
            fp.write(json.dumps(jj))
