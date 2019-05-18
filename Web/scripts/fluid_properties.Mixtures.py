from __future__ import print_function

from CPWeb.BibtexTools import getCitationOrAlternative, getBibtexParser
import CoolProp
import os.path

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
csvfile = os.path.join(web_dir, 'fluid_properties', 'Mixtures.csv')


def merge_args(*args):
    return " :raw-html:`<br/>` ".join(list(args))


def printCoeff(number):
    if number is None or \
       len(str(number).strip()) < 1:
        return " "
    number = float(number)
    short = "{0:.4e}".format(number)
    long = "{0:.14e}".format(number)
    return u':raw-html:`<span title="{1}">{0}</span>`'.format(short, long)


class Dossier:
    def __init__(self):
        self.data = {}

    def add(self, key, value):
        if key not in self.data:
            self.data[key] = []
        self.data[key].append(value)


d = Dossier()

pairs = CoolProp.get('mixture_binary_pairs_list')
print(len(pairs.split(',')))
for pair in pairs.split(','):
    CAS1, CAS2 = pair.split('&')
    d.add('CAS1', CAS1)
    d.add('CAS2', CAS2)
    for key in ['name1', 'name2', 'F', 'function', 'BibTeX', 'xi', 'zeta', 'betaT', 'betaV', 'gammaT', 'gammaV']:
        try:
            d.add(key, CoolProp.CoolProp.get_mixture_binary_pair_data(CAS1, CAS2, key))
        except BaseException as BE:
            d.add(key, '')

import pandas
df = pandas.DataFrame(d.data)
df = df.sort_values(by=['BibTeX', 'name1'], ascending=[0, 1])

bibtexer = getBibtexParser()  # filename = '../../../CoolPropBibTeXLibrary.bib')

with open(csvfile, 'w') as fp:
    header = 'Ref.,Name1,Name2,function,F,'
    header += merge_args("xi", "zeta,")
    header += merge_args("betaT", "betaV,")
    header += merge_args("gammaT", "gammaV")
    header += '\n'
    fp.write(header)

    for index, row in df.iterrows():
        text = ','.join([ \
          getCitationOrAlternative(bibtexer, row['BibTeX']),
          row['name1'],
          row['name2'],
          row['function'],
          row['F'],
          merge_args(printCoeff(row['xi']), printCoeff(row['zeta'])),
          merge_args(printCoeff(row['betaT']), printCoeff(row['betaV'])),
          merge_args(printCoeff(row['gammaT']), printCoeff(row['gammaV']))
        ]) + '\n'
        fp.write(text)
