.. _tabular_interpolation:

**********************
Tabular Interpolation
**********************

Especially when evaluating inputs as a function of pressure and enthalpy (common in many engineering applications), evaluation of the full equation of state is simply too slow, and it is necessary to come up with some means to speed up the calculations.  

As of version 5.1 of CoolProp, the tabular interpolation methods of CoolProp v4 have been brought back from the dead, and significantly improved.  They are approximately 4 times faster than the equivalent methods in v4 of CoolProp due to a more optimized structure.  In order to make the most effective use of the tabular interpolation methods, you must be using the :ref:`low-level interface <low_level_api>`, otherwise significant overhead and slowdown will be experienced.  Thus, this method is best suited to C++, python, and the SWIG wrappers.

There are two backends implemented for tabular interpolation, ``BICUBIC`` and ``TTSE``.  Both consume the same gridded tabular data that is stored to your user home directory in the folder ``HOME/.CoolProp/Tables``.  If you want to find the directory that CoolProp is using as your home directory (``HOME``), you can do something like 

.. ipython::

    In [0]: import CoolProp.CoolProp as CP
    
    In [1]: CP.get_global_param_string("HOME")

This directory is used because the user should under almost all circumstances have read/write access to this folder.  Alternatively, if for some reason you want to use a different directory, you can set the configuration variable ``ALTERNATIVE_TABLES_DIRECTORY`` (see :ref:`configuration`).

.. warning::

    Constructing the tables generates approximately 20 MB of data per fluid. The size of the directory of tabular data can get to be quite large if you use tables with lots of fluids.  By default, CoolProp will warn when this directory is greater than 1 GB in size, and error out at 1.5 x (1.0 GB).  This cap can be lifted by setting a configuration variable (see :ref:`configuration`) as shown here
    
.. ipython::

    In [0]: import CoolProp.CoolProp as CP, json
    
    In [1]: jj = json.loads(CP.get_config_as_json_string())
    
    In [2]: jj['MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB'] = 1.0
    
    In [3]: jj = CP.set_config_as_json_string(json.dumps(jj))

General Information
-------------------

Two sets of tables are generated and written to disk: pressure-enthalpy and pressure-temperature tables.  The tables are generated over the full range of temperatures between the triple point temperature and the maximum temperature, with pressures from the minimum pressure to the maximum pressure.  Pressure-enthalpy and pressure-temperature will be the fastest inputs with these methods because they are the so-called native inputs of the tables.  Other inputs, such as pressure-entropy, require an additional iteration, though the overhead is not as severe as for the full equation of state.

It is critical that you try to only initialize one AbstractState instance and then call its methods. The overhead for generating an AbstractState instance when using TTSE or BICUBIC is not too punitive, but you should try to only do it once.  Each time an instance is generated, all the tabular data is loaded into it.

TTSE Interpolation
------------------

The basic concept behind Tabular Taylor Series Extrapolation (TTSE) extrapolation is that you cache the value and its derivatives with respect to a set of two variables over a regularly (linearly or logarithmically) spaced grid.  It is therefore easy to look up the index of the independent variables ``x`` and ``y`` by interval bisection.  Once the coordinates of the given grid point :math:`i,j` are known, you can used the cached derivatives of the desired variable :math:`z` to obtain it in terms of the given independent variables :math:`x` and :math:`y`:

.. math::

    z = z_{i,j}+\Delta x\left(\frac{\partial z}{\partial x}\right)_{y}+\Delta y\left(\frac{\partial z}{\partial y}\right)_{x}+\frac{1}{2}\Delta x^2\left(\frac{\partial^2 z}{\partial x^2}\right)_{y}+\frac{1}{2}\Delta y^2\left(\frac{\partial^2z}{\partial y^2}\right)_{y}+\Delta x\Delta y\left(\frac{\partial^2z}{\partial y\partial x}\right)
       
.. math::

    \Delta x = x-x_i
    
    \Delta x = y-y_j
    
See the `IAPWS TTSE report <http://www.iapws.org/relguide/TTSE.pdf>`_ for a description of the method.  Analytic derivatives are used to build the tables

Bicubic Interpolation
---------------------

In bicubic interpolation, the values of the output parameter as well as its derivatives are evaluated at the corners of a unit square, and these values are used to fit a bicubic surface over the unit square. `Wikipedia <http://en.wikipedia.org/wiki/Bicubic_interpolation>`_ has excellent coverage of bicubic interpolation, and the method implemented in CoolProp is exactly this method.

Normalized cell values are generated from

.. math::

    \hat x = \frac{x-x_i}{x_{i+1}-x_{i}}
    
    \hat y = \frac{y-y_j}{y_{j+1}-y_{j}}
    
And derivatives must be scaled to be in terms of unit cell values, or 

.. math::

    \frac{\partial z}{\partial \hat x} = \frac{\partial z}{\partial x}\frac{\partial x}{\partial \hat x}
    
    \frac{\partial z}{\partial \hat y} = \frac{\partial z}{\partial y}\frac{\partial y}{\partial \hat y}
    
In CoolProp, after loading the tabular data, the coefficients for all cells are calculated in one shot.

Accuracy comparison
-------------------

Here is a simple comparison of accuracy, the density is obtained for R245fa using the EOS, TTSE extrapolation, and Bicubic interpolation

.. ipython::

    In [0]: import CoolProp
    
    In [1]: HEOS = CoolProp.AbstractState("HEOS", "R245fa")
    
    In [2]: TTSE = CoolProp.AbstractState("TTSE&HEOS", "R245fa")
    
    In [3]: BICU = CoolProp.AbstractState("BICUBIC&HEOS", "R245fa")
    
    In [4]: HEOS.update(CoolProp.PT_INPUTS, 101325, 300); BICU.update(CoolProp.PT_INPUTS, 101325, 300); TTSE.update(CoolProp.PT_INPUTS, 101325, 300)
    
    In [5]: print(HEOS.rhomolar(), TTSE.rhomolar(), BICU.rhomolar())
    
A more complete comparison of the accuracy of these methods can be obtained by studying the following figure for refrigerant R245fa.  You can download the script and change the fluid name to another fluid to investigate the behavior

.. plot::

    import CoolProp
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.ticker
    import numpy as np
    import random

    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_axes((0.08,0.1,0.32,0.83))
    ax2 = fig.add_axes((0.50,0.1,0.32,0.83))

    Ref = 'R245fa'

    BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS',Ref)
    TTSE = CoolProp.AbstractState('TTSE&HEOS',Ref)
    EOS = CoolProp.AbstractState('HEOS',Ref)

    T = np.linspace(CP.PropsSI(Ref,'Tmin')+0.1, CP.PropsSI(Ref,'Tcrit')-0.01, 300)
    pV = CP.PropsSI('P','T',T,'Q',1,Ref)
    hL = CP.PropsSI('Hmass','T',T,'Q',0,Ref)
    hV = CP.PropsSI('Hmass','T',T,'Q',1,Ref)

    HHH1, PPP1, EEE1 = [], [], []
    HHH2, PPP2, EEE2 = [], [], []

    cNorm  = colors.LogNorm(vmin=1e-12, vmax=10)
    scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = plt.get_cmap('jet'))

    for a_useless_counter in range(40000):
            
        h = random.uniform(150000,590000)
        p = 10**random.uniform(np.log10(100000),np.log10(7000000))
        CP.set_debug_level(0)
        try:
            
            EOS.update(CoolProp.HmassP_INPUTS, h, p)
            rhoEOS = EOS.rhomolar(); TEOS = EOS.T()
            
            TTSE.update(CoolProp.HmassP_INPUTS, h, p)
            rhoTTSE = TTSE.rhomolar(); TTTSE = TTSE.T()
            
            BICUBIC.update(CoolProp.HmassP_INPUTS, h, p)
            rhoBICUBIC = BICUBIC.rhomolar(); TBICUBIC = BICUBIC.T()
            
            errorTTSE = abs(rhoTTSE/rhoEOS-1)*100
            errorBICUBIC = abs(rhoBICUBIC/rhoEOS-1)*100
            if errorTTSE > 100 or errorTTSE < 1e-12:
                print(h, p, errorTTSE)

            HHH1.append(h)
            PPP1.append(p)
            EEE1.append(errorTTSE)
            
            HHH2.append(h)
            PPP2.append(p)
            EEE2.append(errorBICUBIC)
            
        except ValueError as VE:
            print('ERROR', VE)
            pass
        
    SC1 = ax1.scatter(HHH1, PPP1, s = 8, c = EEE1, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)
    SC2 = ax2.scatter(HHH2, PPP2, s = 8, c = EEE2, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)

    ax1.set_title('Error in Density from TTSE')
    ax2.set_title('Error in Density from Bicubic')

    for ax in [ax1, ax2]:
        
        ax.set_xlim(250000, 550000)
        ax.set_ylim(100000, 7000000)

        ax.set_yscale('log')
        
        ticks = [100000,200000,400000,600000,800000,1000000,2000000, 4000000, 6000000]
        labels = [str(tick) for tick in ticks]
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        
        ticks = [150000, 250000,350000,450000,550000]
        labels = [str(tick) for tick in ticks]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)

        ax.tick_params(axis='y',which='minor', left='off')

        ax.set_xticklabels(ax.get_xticks()/1e3)
        ax.set_xlabel('Enthalpy [kJ/kg]')
        ax.set_yticklabels(ax.get_yticks()/10**3)
        ax.set_ylabel('Pressure [kPa]')

        ax.plot(hL,pV,'k',lw = 4)
        ax.plot(hV,pV,'k',lw = 4)

    cbar_ax = fig.add_axes([0.85, 0.15, 0.06, 0.7])
    CB = fig.colorbar(SC1, cax=cbar_ax)
    CB.set_label(r'$(\rho/\rho_{EOS}-1)\times 100$ [%]')

Speed comparison
----------------

The primary motivation for the use of tabular interpolation is the improvement in computational speed.  Thus a small summary could be useful.  This tabular data was obtained by this python script : :download:`(link to script) <speed_script.py>`.

.. include :: tabular_data.rst.in

More Information
----------------

The tables are stored in a zipped format using the msgpack package and miniz.  If you want to see what data is serialized in the tabular data, you can unzip and unpack into python (or other high-level languages) using something roughly like::

    import msgpack, zlib, io, numpy as np, matplotlib.pyplot as plt

    root = r'C:\Users\ian\.CoolProp\Tables\REFPROPMixtureBackend(R32[0.8292500000]&R1234yf[0.1707500000])'
    with open(root+'/single_phase_logph.bin.z','rb') as fp:
        values = msgpack.load(io.BytesIO(zlib.decompress(fp.read())))
        revision, matrices = values[0:2]
        T,h,p,rho = np.array(matrices['T']), np.array(matrices['hmolar']), np.array(matrices['p']), np.array(matrices['rhomolar'])
        plt.plot(np.array(matrices['p']),np.array(matrices['hmolar']),'x')
    with open(root+'/phase_envelope.bin.z','rb') as fp:
        values = msgpack.load(io.BytesIO(zlib.decompress(fp.read())))
        revision, matrices = values[0:2]
        plt.plot(np.array(matrices['p']),np.array(matrices['hmolar_vap']),'-')
    plt.show()

    with open(root+'/single_phase_logpT.bin.z','rb') as fp:
        values = msgpack.load(io.BytesIO(zlib.decompress(fp.read())))
        revision, matrices = values[0:2]
        T,h,p,rho = np.array(matrices['T']), np.array(matrices['hmolar']), np.array(matrices['p']), np.array(matrices['rhomolar'])
        plt.plot(np.array(matrices['p']),np.array(matrices['T']),'x')
    with open(root+'/phase_envelope.bin.z','rb') as fp:
        values = msgpack.load(io.BytesIO(zlib.decompress(fp.read())))
        revision, matrices = values[0:2]
        plt.plot(np.array(matrices['p']),np.array(matrices['T']),'-')
    plt.show()

You'll need msgpack wrapper for your target language.        
