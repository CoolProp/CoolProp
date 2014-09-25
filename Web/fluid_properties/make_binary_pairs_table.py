import CoolProp

class Dossier:
    def __init__(self):
        self.data = {}
    def add(self, key, value):
        if key not in self.data:
            self.data[key] = []
        self.data[key].append(value)
            
d = Dossier()

pairs = CoolProp.get('mixture_binary_pairs_list')
for pair in pairs.split(','):
    CAS1, CAS2 = pair.split('&')
    d.add('CAS1', CAS1)
    d.add('CAS2', CAS2)
    for key in ['name1','name2','F','function','BibTeX','xi','zeta','betaT','betaV','gammaT','gammaV']:
        try:
            d.add(key, CoolProp.CoolProp.get_mixture_binary_pair_data(CAS1, CAS2, key))
        except BaseException as BE:
            d.add(key, '')
            
import pandas
df = pandas.DataFrame(d.data)
df = df.sort(['BibTeX','name1'], ascending = [0, 1])

with open('mixture_binary_pairs.csv','w') as fp:
    fp.write('Ref.,Name1,Name2,function,F,xi,zeta,betaT,betaV,gammaT,gammaV\n')
    for index, row in df.iterrows():
        text = ','.join([':cite:`'+row['BibTeX']+'`', row['name1'], row['name2'], row['function'], row['F'], row['xi'], row['zeta'], row['betaT'], row['betaV'], row['gammaT'], row['gammaV']])+'\n'
        fp.write(text)