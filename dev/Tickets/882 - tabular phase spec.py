import msgpack, zlib, StringIO
import matplotlib.pyplot as plt, numpy as np, CoolProp

# Modify path as necessary
with open(r'C:\Users\Belli\.CoolProp\Tables\HelmholtzEOSBackend(Water[1.0000000000])/single_phase_logpT.bin.z','rb') as fp:
    ph = zlib.decompress(fp.read())
    values = msgpack.load(StringIO.StringIO(ph))
    revision, matrices = values[0:2]
    T,h,p,rho = np.array(matrices['T']), np.array(matrices['hmolar']), np.array(matrices['p']), np.array(matrices['rhomolar'])
    Tt = CoolProp.CoolProp.PropsSI('Ttriple','Water')
    Tc = CoolProp.CoolProp.PropsSI('Tcrit','Water')
    Ts = np.linspace(Tt, Tc)
    ps = CoolProp.CoolProp.PropsSI('P','T',Ts,'Q',0,'Water')
    
    plt.plot(T, p, '.')
    plt.plot(Ts, ps, 'k', lw = 2)
    plt.yscale('log')
    plt.show()
    
    
# coding: utf-8

# In[ ]:

import CoolProp, numpy as np
AS = CoolProp.AbstractState('BICUBIC&HEOS','CO2')
print AS.rhomolar_critical()
pt = AS.keyed_output(CoolProp.iP_triple)
pc = AS.p_critical()

for p in np.logspace(np.log10(pt*1.05),np.log10(pc*0.95)):
    AS.update(CoolProp.PQ_INPUTS, p, 0)
    rhoL = AS.rhomolar()
    Ts = AS.T()    
    # Liquid side
    for specify_phase in [True]:
        if specify_phase:
            AS.specify_phase(CoolProp.iphase_liquid)
        else:
            AS.unspecify_phase()
        for dT in [10,1,0.01,0.001]:
            print p, Ts-dT
            try:
                AS.update(CoolProp.PT_INPUTS, p, Ts-dT)
                print p, Ts-dT, AS.rhomolar()
            except BaseException as BE:
                print 'ERROR:', BE
        print p, Ts, rhoL
    print ' '
    
                


# In[ ]:



