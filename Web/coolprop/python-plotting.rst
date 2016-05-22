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
