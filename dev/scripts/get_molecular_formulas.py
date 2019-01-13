import urllib, sys, CoolProp, json

sys.path.append('C:\RDKit_2014_09_2')
from rdkit import Chem

from collections import Counter

for fluid in CoolProp.__fluids__:
    CAS = CoolProp.CoolProp.get_fluid_param_string(fluid, "CAS")
    print(fluid, CAS, end=' ')
    if '.ppf' in CAS or '.PPF' in CAS or 'o' in CAS or 'p' in CAS:
        print('')
        continue

    txt = urllib.urlretrieve('http://cactus.nci.nih.gov/chemical/structure/' + CAS + '/file?format=mol')

    with open(txt[0], 'r') as fp:
        contents = fp.read()

    if '<h1>Page not found (404)</h1>' in contents:
        print('MISSING FLUID')
        continue

    mol = Chem.MolFromMolBlock(contents)
    mol = Chem.AddHs(mol)

    c = dict(Counter(atom.GetSymbol() for atom in mol.GetAtoms()))
    formula = ''.join([k + '_{' + str(c[k]) + '}' for k in sorted(c.keys())])

    fname = '../fluids/' + fluid + '.json'
    with open(fname, 'r') as fp:
        jj = json.load(fp)
    jj['FORMULA'] = formula
    with open(fname, 'w') as fp:
        fp.write(json.dumps(jj))
