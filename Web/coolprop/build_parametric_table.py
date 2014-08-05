import CoolProp5
import pandas
grouping = dict()
grouping2 = []
# Group aliases
for parameter in CoolProp5.get('parameter_list').split(','):
    
    index = CoolProp5.CoolProp.get_parameter_index(parameter)
    units = CoolProp5.CoolProp.get_parameter_information(index, 'units').replace('-',' ')
    IO = CoolProp5.CoolProp.get_parameter_information(index, 'IO')
    long = CoolProp5.CoolProp.get_parameter_information(index, 'long')
    short = CoolProp5.CoolProp.get_parameter_information(index, 'short')
    
    RHS = (units, IO, long)
    if RHS not in grouping:
        grouping[RHS] = [parameter]
    else:
        grouping[RHS].append(parameter)
    
for k, v in grouping.iteritems():
    grouping2.append([', '.join(['``'+_+'``' for _ in v])] + list(k))
    
headers = ['Parameter','Units','Input/Output','Description']

df3 = pandas.DataFrame(grouping2, columns = headers)
df4 = df3.sort_index(by = ['Input/Output', 'Parameter'])
grouping2 = [row for row in df4.values]

N = []
for i in range(len(grouping2[0])):
    N.append(max([len(el[i]) for el in grouping2]))

for i in range(len(N)):
    if N[i] < len(headers[i]):
        N[i] = len(headers[i])
        
top_line = '='*N[0] + ' ' + '='*N[1] + ' ' + '='*N[2] + ' ' + '='*N[3]
header = ' '.join([h.ljust(n) for h,n in zip(headers,N)])

fp = open('parameter_table.rst.in', 'w')
fp.write('.. constructed with the build_parameter_table.py script in this folder \n\n')
fp.write(top_line + '\n')
fp.write(header + '\n')
fp.write(top_line + '\n')
for line in grouping2:
    fp.write(' '.join([h.ljust(n) for h,n in zip(line,N)]) + '\n')
fp.write(top_line + '\n')
fp.close()

