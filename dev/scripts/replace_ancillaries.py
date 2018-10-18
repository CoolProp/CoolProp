#!/Users/ian/anaconda/bin/python

import json, glob, CoolProp

for fluid in glob.glob('../fluids/*.json'):
    with open(fluid, 'r') as fp:
        jj = json.load(fp)

    pL = jj['ANCILLARIES'].pop('pL')
    pV = jj['ANCILLARIES'].pop('pV')
    # Keep the one with the lower error
    if pL['max_abserror_percentage'] < pV['max_abserror_percentage']:
        pS = pL
    else:
        pS = pV

    pseudo_pure = jj['EOS'][0]['pseudo_pure']
    if pseudo_pure:
        print('-----------------PSEUDO (SKIPPING !!!) %s' % fluid)
    else:
        print(fluid)
        jj['ANCILLARIES']['pS'] = pS
        with open(fluid, 'w') as fp:
            json.dump(jj, fp)
