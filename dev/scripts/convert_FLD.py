import CoolProp, json


def extract_coeffs(lines):
    N_normal, C_normal, N_critical, C_critical, N_spare, C_spare = [int(el) for el in lines[0].split('!')[0].strip().split(' ') if el]

    terms = []

    assert(N_normal > 0)
    assert(C_normal == 4)

    chunk = []
    for line in lines[1:N_normal + 1]:
        ll = line.split('!')[0].strip()
        ll = [float(f) for f in ll.split(' ') if f]
        chunk.append(ll)
    n, t, d, l = zip(*chunk)

    terms.append(
        dict(type="ResidualHelmholtzPower",
         n=n, d=d, t=t, l=l)
        )

    if N_critical > 0:
        chunk = []
        for line in lines[1 + N_normal:N_normal + 1 + N_critical]:
            ll = line.split('!')[0].strip()
            ll = [float(f) for f in ll.split(' ') if f]
            chunk.append(ll)
        n, t, d, l, dummy1, neg_eta, neg_beta, gamma, epsilon, dummy2, dummy3, dummy4 = zip(*chunk)
        eta, beta = [[-el for el in vec] for vec in [neg_eta, neg_beta]]

        terms.append(
            dict(type="ResidualHelmholtzGaussian",
             n=n, d=d, t=t, eta=eta, beta=beta, gamma=gamma, epsilon=epsilon)
            )

    return terms


def extract_active_EOS(file_name):
    with open(file_name, 'r') as fp:
        lines = fp.readlines()
    del fp

    # Find the #EOS key
    iEOS = -1
    for i, line in enumerate(lines):
        if line.startswith('#EOS'):
            if iEOS < 0:
                iEOS = i
            else:
                raise KeyError("Found more than one #EOS in file")

    # Find the end of the EOS block
    for i in range(iEOS, len(lines)):
        if len(lines[i].strip()) == 0:
            iEOS_end = i
            break

        if i == len(lines) - 1:
            raise KeyError("Could not find empty line at end of EOS block")

    iEOS_start = iEOS
    for i in range(iEOS, iEOS_end + 1):
        if lines[i].startswith('!'):
            iEOS_start = i
            break

    for i in range(iEOS_start + 1, iEOS_end + 1):
        print(lines[i].strip())

    def unpack(line):
        s = line.split('!')[0].strip()
        try:
            return float(s)
        except ValueError:
            return s

    Tmin, Tmax, pmax, rhomax, CPP, MW, Tt, Pt, rhot, NBP, acentric, crit, reducing, R = [unpack(line) for line in lines[iEOS_start + 1:iEOS_start + 1 + 14]]
    reducing = [float(f) for f in reducing.split(' ') if f]
    crit = [float(f) for f in crit.split(' ') if f]

    alphar = extract_coeffs(lines[iEOS_start + 15:iEOS_end + 1])

    del file_name, iEOS, iEOS_start, iEOS_end, unpack, lines, line, i
    return locals().copy()


def build_states():
    return []


def build_alpha0(ALPHA0):
    return []


def build_alphar(EOS):
    return []


def write(EOS, ALPHA0, json_file):

    base = dict(BibTeX_CP0="",
                BibTeX_EOS="XXXX TO BE FILLED IN XXX",
                STATES=build_states(),
                T_max=EOS['Tmax'],
                T_max_units='K',
                Ttriple=EOS['Tt'],
                Ttriple_units='K',
                acentric=EOS['acentric'],
                acentric_units='-',
                alpha0=build_alpha0(ALPHA0),
                alphar=EOS['alphar'],
                gas_constant=EOS['R'],
                gas_constant_units='J/mol/K',
                molar_mass=EOS['MW'] / 1000.0,
                molar_mass_units='kg/mol',
                p_max=EOS['pmax'] * 1000.0,
                p_max_units='Pa',
                pseudo_pure='XXX TO BE FILLED IN XXX'
                )

    with open(json_file, 'w') as fp:
        json.dump(base, fp, indent=2, sort_keys=True)


if __name__ == '__main__':
    EOS = extract_active_EOS('/Users/ian/Downloads/MDM.FLD')
    ALPHA0 = []
    write(EOS, ALPHA0, 'MDM.json')
