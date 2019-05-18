import CoolProp
import pandas
import six
grouping = dict()
grouping2 = []
# Group aliases
for parameter in CoolProp.get('parameter_list').split(','):

    index = CoolProp.CoolProp.get_parameter_index(parameter)
    units = CoolProp.CoolProp.get_parameter_information(index, 'units').replace('-', ' ')
    IO = CoolProp.CoolProp.get_parameter_information(index, 'IO')
    long = CoolProp.CoolProp.get_parameter_information(index, 'long')
    short = CoolProp.CoolProp.get_parameter_information(index, 'short')
    trivial = str(CoolProp.CoolProp.is_trivial_parameter(index))

    RHS = (units, IO, trivial, long)
    if RHS not in grouping:
        grouping[RHS] = [parameter]
    else:
        grouping[RHS].append(parameter)

for k, v in six.iteritems(grouping):
    grouping2.append([', '.join(['``' + _ + '``' for _ in v])] + list(k))

headers = ['Parameter', 'Units', 'Input/Output', 'Trivial', 'Description']

df3 = pandas.DataFrame(grouping2, columns=headers)
df4 = df3.sort_values(by=['Input/Output', 'Parameter'])
grouping2 = [row for row in df4.values]

N = []
for i in range(len(grouping2[0])):
    N.append(max([len(el[i]) for el in grouping2]))

for i in range(len(N)):
    if N[i] < len(headers[i]):
        N[i] = len(headers[i])

top_line = '=' * N[0] + ' ' + '=' * N[1] + ' ' + '=' * N[2] + ' ' + '=' * N[3] + ' ' + '=' * N[4]
header = ' '.join([h.ljust(n) for h, n in zip(headers, N)])

fp = open('../coolprop/parameter_table.rst.in', 'w')
fp.write('.. constructed with the coolprop.parameter_table.py script in the web/scripts folder \n\n')
fp.write(top_line + '\n')
fp.write(header + '\n')
fp.write(top_line + '\n')
for line in grouping2:
    fp.write(' '.join([h.ljust(n) for h, n in zip(line, N)]) + '\n')
fp.write(top_line + '\n')
fp.close()
