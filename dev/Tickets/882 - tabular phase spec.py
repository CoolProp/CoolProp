# coding: utf-8

import msgpack, zlib, StringIO
import matplotlib.pyplot as plt, numpy as np, CoolProp

# Modify path as necessary
user_paths = [r'C:\Users\Belli', r'C:\Users\jowr']
with open(user_paths[1] + r'\.CoolProp\Tables\HelmholtzEOSBackend(Water[1.0000000000])/single_phase_logpT.bin.z', 'rb') as fp:
    ph = zlib.decompress(fp.read())
    values = msgpack.load(StringIO.StringIO(ph))
    revision, matrices = values[0:2]
    T, h, p, rho = np.array(matrices['T']), np.array(matrices['hmolar']), np.array(matrices['p']), np.array(matrices['rhomolar'])
    Tt = CoolProp.CoolProp.PropsSI('Ttriple', 'Water')
    Tc = CoolProp.CoolProp.PropsSI('Tcrit', 'Water')
    Ts = np.linspace(Tt, Tc)
    ps = CoolProp.CoolProp.PropsSI('P', 'T', Ts, 'Q', 0, 'Water')

    plt.plot(T, p, '.', color='gray')
    plt.plot(Ts, ps, 'k', lw=2)
    plt.yscale('log')
    plt.xlabel('Temperature / K')
    plt.ylabel('Pressure / Pa')
    plt.show()


# coding: utf-8

# In[ ]:

import CoolProp, numpy as np
AS = CoolProp.AbstractState('BICUBIC&HEOS', 'CO2')
print(AS.rhomolar_critical())
pt = AS.keyed_output(CoolProp.iP_triple)
pc = AS.p_critical()

verbose = False
dTs = np.power(10.0, -np.arange(-1, 4))

for p in np.logspace(np.log10(pt * 1.05), np.log10(pc * 0.95)):
    AS.update(CoolProp.PQ_INPUTS, p, 0)
    rhoL = AS.rhomolar()
    Ts = AS.T()
    # Liquid side
    for specify_phase in [False, True]:
        if specify_phase:
            AS.specify_phase(CoolProp.iphase_liquid)
        else:
            AS.unspecify_phase()
        for dT in dTs:
            if verbose: print("%s %s" % (p, Ts - dT))
            try:
                AS.update(CoolProp.PT_INPUTS, p, Ts - dT)
                if verbose: print("%s %s %s" % (p, Ts - dT, AS.rhomolar()))
            except BaseException as BE:
                if specify_phase: print('Liquid error: %s' % BE)
                else: print('Liquid issue: %s' % BE)
        if verbose: print(p, Ts, rhoL)
    # Gaseous side
    for specify_phase in [False, True]:
        if specify_phase:
            AS.specify_phase(CoolProp.iphase_gas)
        else:
            AS.unspecify_phase()
        for dT in dTs:
            if verbose: print("%s %s" % (p, Ts + dT))
            try:
                AS.update(CoolProp.PT_INPUTS, p, Ts + dT)
                if verbose: print(p, Ts + dT, AS.rhomolar())
            except BaseException as BE:
                if specify_phase: print('   Gas error: %s' % BE)
                else: print('   Gas issue: %s' % BE)
        if verbose: print("%s %s %s" % (p, Ts, rhoL))
