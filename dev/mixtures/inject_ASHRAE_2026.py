from pathlib import Path
import numpy as np
import pandas 
import CoolProp.CoolProp as CP
import json

new_entries = []

_existing = Path(__file__).parent / 'predefined_mixtures.json'
_existing = json.load(open(_existing, 'r')) if _existing.exists() else {}
existing = {entry['name']: entry for entry in _existing}

def _get_molar_mass(component):
    try:
        return CP.PropsSI('M', component)
    except Exception:
        if component == 'R143b':
            return 0.084041 # From R143a, which has the same molar mass because they are isomers
        if component == 'RE170':
            return CP.PropsSI('M', 'DME')
        if component == 'R1132a':
            return 0.064035
        if component == 'R1132(E)':
            return 0.064035
        raise ValueError(f'Unknown component: {component}')
        
for ir, row in pandas.read_csv(Path(__file__).parent / 'ASHRAE_predefined_2026.tsv', sep='\t', comment='#').iterrows():
    # print (row[0], row[1])

    # Convert mass fractions to mole fractions
    components, mass_fractions_string = row[1].removeprefix('R-').rsplit('(', 1)
    components = ['R' + component.strip() for component in components.split('/')]
    mass_fractions = [float(x) for x in mass_fractions_string.strip().removesuffix(')').split('/')]
    total_mass_fraction = sum(mass_fractions)
    mass_fractions = [x / total_mass_fraction for x in mass_fractions]
    mole_fractions = []
    for i, component in enumerate(components):
        molar_mass = _get_molar_mass(component)
        mole_fraction = mass_fractions[i] / molar_mass
        mole_fractions.append(mole_fraction)
    total_mole_fraction = sum(mole_fractions)
    mole_fractions = [x / total_mole_fraction for x in mole_fractions]

    name = f'R{row[0]}'
    if name in existing:
        print (f'{name} already exists, skipping')
        old_mole_fractions = existing[name]['mole_fractions']
        new_mole_fractions = mole_fractions
        if not np.allclose(old_mole_fractions, new_mole_fractions):
            print (f'{name} has different mole fractions, skipping')
            print (f'Old: {old_mole_fractions}')
            print (f'New: {new_mole_fractions}')
        continue

    try:
        CP.AbstractState('HEOS', name)
    except BaseException as BE:
        msg = str(BE)
        print(msg)

    new_entries.append({
        'fluids': components,
        'mole_fractions': mole_fractions,   
        'name': f'R{row[0]}'
    })
print(json.dumps(new_entries, indent=2))

contents = json.load((Path(__file__).parent / 'predefined_mixtures.json').open('r')) if (Path(__file__).parent / 'predefined_mixtures.json').exists() else []
with open(Path(__file__).parent / 'predefined_mixtures.json', 'w') as f:
    f.write(json.dumps(contents + new_entries, indent=2))