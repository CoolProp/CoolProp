from CoolProp.CoolProp import get_fluid_param_string

lines = open('KunzWagner2012_TableA7.txt', 'r').read()

template = """{{"Name1" : "{Name1:s}",
"Name2" : "{Name2:s}",
"CAS1" : "{CAS1:s}",
"CAS2" : "{CAS2:s}",
"d" : {d:s},
"t" : {t:s},
"n" : {n:s},
"eta" : {eta:s},
"epsilon" : {epsilon:s},
"beta": {beta:s},
"gamma": {gamma:s}
}},"""

chunks = lines.split('\n\n')

for chunk in chunks:
    lines = chunk.split('\n')

    D, T, N, ETA, EPSILON, BETA, GAMMA = [0], [0], [0], [0], [0], [0], [0]
    names = lines.pop(0)

    for line in lines:
        vals = line.strip().split(' ')

        if len(vals) == 4:
            i, d, t, n = vals
            eta = 0
            epsilon = 0
            beta = 0
            gamma = 0
        else:
            i, d, t, n, eta, epsilon, beta, gamma = vals

        D.append(int(d))
        T.append(float(t))
        N.append(float(n))
        ETA.append(float(eta))
        EPSILON.append(float(epsilon))
        BETA.append(float(beta))
        GAMMA.append(float(gamma))

    name1, name2 = names.split('/')

    CAS1 = get_fluid_param_string(name1, 'CAS')
    CAS2 = get_fluid_param_string(name2, 'CAS')

    print(template.format(Name1=name1,
                          Name2=name2,
                          CAS1=CAS1,
                          CAS2=CAS2,
                          d=str(D),
                          t=str(T),
                          n=str(N),
                          eta=str(ETA),
                          epsilon=str(EPSILON),
                          beta=str(BETA),
                          gamma=str(GAMMA)
                          ))
