import json, CoolProp.CoolProp as CP

with open('Bell2016Fluids.txt', 'r') as fp:
    lines = fp.readlines()
name2CAS = {}
for line in lines:
    Name, longname, CAS, N = line.strip().split(' ')
    name2CAS[Name] = CAS

with open('Bell2016Coefficients.txt', 'r') as fp:
    lines = fp.readlines()

with open('mixture_binary_pairs.json', 'r') as fp:
    jj = json.load(fp)

CAS_pairs = []
for pair in jj:
    CAS_pairs.append((pair['CAS1'], pair['CAS2']))

for line in lines:
    Name1, Name2, betaT, gammaT, MARE, N = line.strip().split(' ')
    CAS1 = name2CAS[Name1]
    CAS2 = name2CAS[Name2]
    if ((CAS1, CAS2)) not in CAS_pairs and ((CAS2, CAS1)) not in CAS_pairs:
        entry = {
            "BibTeX": "Bell-JCED-2016",
            "CAS1": CAS1,
            "CAS2": CAS2,
            "F": 0.0,
            "Name1": Name1,
            "Name2": Name2,
            "betaT": float(betaT),
            "gammaT": float(gammaT),
            "betaV": 1.0,
            "gammaV": 1.0
            }
        jj.append(entry)
    else:
        print('Skipping %s %s' % (Name1, Name2))
with open('mixture_binary_pairs.json', 'w') as fp:
    json.dump(jj, fp, indent=2, sort_keys=True)
