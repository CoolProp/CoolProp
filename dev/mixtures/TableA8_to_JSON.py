from CoolProp.CoolProp import get_fluid_param_string

# Map from chemical formula to name
name_map = {'CH4': 'Methane',
            'N2': 'Nitrogen',
            'O2': 'Oxygen',
            'CO2': 'CarbonDioxide',
            'CO': 'CarbonMonoxide',
            'H2S': 'HydrogenSulfide',
            'H2': 'Hydrogen',
            'H2O': 'Water',
            'He': 'Helium',
            'Ar': 'Argon',
            'C2H6': 'Ethane',
            'C3H8': 'Propane',
            'nC4H10': 'n-Butane',
            'iC4H10': 'IsoButane',
            'nC5H12': 'n-Pentane',
            'iC5H12': 'Isopentane',
            'nC6H14': 'n-Hexane',
            'nC7H16': 'n-Heptane',
            'nC8H18': 'n-Octane',
            'nC9H20': 'n-Nonane',
            'nC10H22': 'n-Decane',
            }

F_factors = {}
lines = open('KunzWagner2012_TableA6.txt', 'r').readlines()
for line in lines:
    names, F = line.strip().split(' ')
    name1, name2 = names.split('/')
    CAS1 = get_fluid_param_string(name1, 'CAS')
    CAS2 = get_fluid_param_string(name2, 'CAS')
    F_factors[(CAS1, CAS2)] = F
    F_factors[(CAS2, CAS1)] = F

lines = open('KunzWagner2012_TableA8.txt', 'r').readlines()

template = """{{"Name1" : "{Name1:s}",
"Name2" : "{Name2:s}",
"CAS1" : "{CAS1:s}",
"CAS2" : "{CAS2:s}",
"betaV" : {betaV:s},
"gammaV" : {gammaV:s},
"betaT" : {betaT:s},
"gammaT" : {gammaT:s},
"F" : {F:s}
}},"""

for line in lines:
    vals = line.strip().split(' ')

    if len(vals) == 6:
        names, a, betav, gammav, betaT, gammaT = vals
    else:
        names, betav, gammav, betaT, gammaT = vals

    name1, name2 = names.split('-')
    name1 = name_map[name1]
    name2 = name_map[name2]

    CAS1 = get_fluid_param_string(name1, 'CAS')
    CAS2 = get_fluid_param_string(name2, 'CAS')

    if (CAS1, CAS2) in F_factors:
        F = F_factors[(CAS1, CAS2)]
    else:
        F = '0.0'

    print(template.format(Name1=name1,
                          Name2=name2,
                          CAS1=CAS1,
                          CAS2=CAS2,
                          betaV=betav,
                          gammaV=gammav,
                          betaT=betaT,
                          gammaT=gammaT,
                          F=F))
