data = """CO2-H2O  1.030538  0.828472  1.021392  0.895156  1
CO2-N2  0.994140013  1.107654104  1.022709642  1.047578256  1
CO2-O2  1.000000  1.031986  1.000000  1.084460  0
CO2-Ar  1.027147  0.968781  1.001378  1.029710  1
CO2-CO  0.993245  1.068392  1.030855  1.245499  0
H2O-N2  0.954149  0.805147  1.079628  0.733443  1
H2O-O2  0.798046  0.807842  0.972576  0.873460  0.6017
H2O-Ar  0.679104  0.921000  0.940398  1.050952  0
H2O-CO  1.045927  0.823984  1.063348  0.766756  0.9897
N2-O2  0.997190589  0.995157044  0.999521770  0.997082328  0
N2-Ar  0.999442  0.989311  1.006697  1.001549  0
N2-CO  1.002409  0.994100  1.000000  1.001317  0
O2-Ar  0.999039  0.988822  1.006502  1.001341  0
O2-CO  1.000000  1.000000  1.000000  1.000000  0
CO-Ar  1.000000000  0.954215746  1.000000000  1.159720623  0"""

namedict = dict(O2='Oxygen', N2='Nitrogen', CO2='CarbonDioxide', CO='CarbonMonoxide', H2O='Water', Ar='Argon')

import CoolProp
CASdict = {namedict[n]: CoolProp.CoolProp.get_fluid_param_string(namedict[n], "CAS") for n in namedict}
functiondict = {'CO2-H2O': 'CarbonDioxide-Water',
               'CO2-N2': 'CarbonDioxide-Nitrogen',
               'CO2-Ar': 'CarbonDioxide-Argon',
               'H2O-N2': 'GeneralizedAirWater',
               'H2O-O2': 'GeneralizedAirWater',
               'H2O-CO': 'GeneralizedAirWater'}
out = []
for line in data.split('\n'):
    pair, betaT, betaV, gammaT, gammaV, F = line.split('  ')
    n1, n2 = pair.split('-')

    out.append(dict(BibTeX='Gernert-Thesis-2013',
                    F=float(F),
                    betaT=float(betaT),
                    betaV=float(betaV),
                    gammaT=float(gammaT),
                    gammaV=float(gammaV),
                    Name1=namedict[n1],
                    Name2=namedict[n2],
                    CAS1=CASdict[namedict[n1]],
                    CAS2=CASdict[namedict[n2]]))

    if F != '0':
        out[-1]['function'] = functiondict[pair]

import json, sys
sys.path.append('..')
from package_json import json_options
print(json.dumps(out, **json_options))
