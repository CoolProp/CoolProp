.. _python-plotting:

Python Plotting
===============

The following example can be used to create a temperature-entropy (T,s) plot for
propane (R290) with isolines for the vapour quality in steps of 0.2:

.. plot::
    :include-source:

    import CoolProp
    from CoolProp.Plots import PropertyPlot
    ts_plot = PropertyPlot('R290','Ts')
    ts_plot.calc_isolines(CoolProp.iQ, num=6)
    ts_plot.show()

The following example can be used to create a pressure-enthalpy (log p,h) plot for 
R-134a with a couple of isotherms and isentropic lines:

.. plot::
    :include-source:
  
    import CoolProp
    from CoolProp.Plots import PropertyPlot
    plot = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR', tp_limits='ACHP')
    plot.calc_isolines(CoolProp.iQ, num=11)
    plot.calc_isolines(CoolProp.iT, num=25)
    plot.calc_isolines(CoolProp.iSmass, num=15)
    plot.show()

Here is an example for R-245fa using two different ways to generate the isobars:
    
.. plot::
    :include-source:
    import CoolProp
    from CoolProp.Plots import PropertyPlot
    plot = PropertyPlot('HEOS::R245fa', 'TS', unit_system='EUR', tp_limits='ORC')
    plot.calc_isolines(CoolProp.iQ, num=11)
    plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=True)
    plot.draw()
    plot.isolines.clear()
    plot.props[CoolProp.iP]['color'] = 'green'
    plot.props[CoolProp.iP]['lw'] = '0.5'
    plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=False)
    plot.show()
    
The available plots are listed in the :py:class:`CoolProp.Plots.Common.Base2DObject`  
and at this time (May 2016) this list contains :math:`T,s`, :math:`p,h`, 
:math:`h,s`, :math:`p,s`, :math:`p,\rho`, :math:`T,\rho`, :math:`p,T` and 
:math:`p,u` plots. 

Some of the commonly used `Matplotlib <http://www.matplotlib.org>`_ functions,
such as :func:`title`, :func:`grid`, :func:`xlabel` and :func:`ylabel` have been wrapped in
the :py:class:`CoolProp.Plots.Common.BasePlot` class to make the plotting of
graphs a little simpler, for example:

.. plot::
    :include-source:

    from CoolProp.Plots import PropertyPlot
    ts_plot = PropertyPlot('Water', 'Ts')
    ts_plot.calc_isolines(CoolProp.iQ, num=11)
    ts_plot.title(r'$T,s$ Graph for Water')
    ts_plot.xlabel(r'$s$ [{kJ}/{kg K}]')
    ts_plot.ylabel(r'$T$ [K]')
    ts_plot.grid()
    ts_plot.show()


    
Cycle Calculations
==================

It is also possible to carry out simple thermodynamic cycle calculations with the 
CoolProp classes. These calculations are based on the utility classes 
:py:class:`CoolProp.Plots.SimpleCycles.StatePoint` and 
:py:class:`CoolProp.Plots.SimpleCycles.StateContainer`, which can be used on their 
own as demonstrated below. Note that the utility classes support numerous notations
to access their members and you can chose the one you like best or mix them:

.. ipython::

    In [0]: from __future__ import print_function
    In [0]: import CoolProp
    In [0]: from CoolProp.Plots.SimpleCycles import StateContainer
    In [0]: T0 = 300.000; p0 = 200000.000; h0 = 112745.749; s0 = 393.035
    In [0]: cycle_states = StateContainer()
    In [0]: cycle_states[0,'H'] = h0
    In [0]: cycle_states[0]['S'] = s0
    In [0]: cycle_states[0][CoolProp.iP] = p0
    In [0]: cycle_states[0,CoolProp.iT] = T0
    In [0]: cycle_states[1,"T"] = 300.064
    In [0]: print(cycle_states)

    
The utility classes were designed to work well with the plotting objects described above
and this example illustrates how a simple Rankine cycle can be added to to a :math:`T,s` 
graph:

.. plot::
    :include-source:
    import CoolProp
    from CoolProp.Plots.Plots import PropertyPlot
    from CoolProp.Plots.SimpleCycles import SimpleRankineCycle
    pp = PropertyPlot('HEOS::Water', 'TS', unit_system='EUR')
    pp.calc_isolines(CoolProp.iQ, num=11)
    cycle = SimpleRankineCycle('HEOS::Water', 'TS', unit_system='EUR')
    T0 = 300
    pp.state.update(CoolProp.QT_INPUTS,0.0,T0+15)
    p0 = pp.state.keyed_output(CoolProp.iP)
    T2 = 700
    pp.state.update(CoolProp.QT_INPUTS,1.0,T2-150)
    p2 = pp.state.keyed_output(CoolProp.iP)
    cycle.simple_solve(T0, p0, T2, p2, 0.7, 0.8, SI=True)
    cycle.steps = 50
    sc = cycle.get_state_changes()
    pp.draw_process(sc) 

