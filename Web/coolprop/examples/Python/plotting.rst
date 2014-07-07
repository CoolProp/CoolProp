.. _python-plotting:

Python Plotting
===============

.. note::
    The python plotting API has been changed as of v4.0. The examples shown
    on this page use the new python plotting API. Examples using the old
    python plotting API can be found here :ref:`python-plotting-old`.

The following example can be used to create a Temperature-Entropy plot for
propane (R290):

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ts_plot = PropsPlot('R290', 'Ts')
    ts_plot.show()


The following example can be used to create a Pressure-Enthalpy plot for R410A:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ph_plot = PropsPlot('R410A', 'Ph')
    ph_plot.show()

The available plots are:

== ====================
PT Pressure-Temperature
PD Pressure-Density
PH Pressure-Enthalpy
PS Pressure-Entropy
TD Temperature-Density
TS Temperatre-Entropy
HS Enthalpy-Entropy
== ====================


The following, more advanced example, can be used to draw lines of constant
properties for n-Pentane. Note the different ways to invoke the
:py:func:`CoolProp.Plots.Plots.PropsPlot.draw_isolines` function draw:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ref_fluid = 'n-Pentane'
    ts_plot = PropsPlot(ref_fluid, 'Ts')
    ts_plot.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
    ts_plot.draw_isolines('P', [100, 2000], num=5)
    ts_plot.draw_isolines('D', [2, 600], num=7)
    ts_plot.set_axis_limits([-2, 1.5, 200, 500])
    ts_plot.show()

Some of the commonly used `Matplotlib <http://www.matplotlib.org>`_ functions,
such as :func:`title`, :func:`xlabel` and :func:`ylabel` have been wrapped in
the :py:class:`CoolProp.Plots.Plots.PropsPlot` class to make the plotting of
graphs a little simpler, for example:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ts_plot = PropsPlot('Water', 'Ts')
    ts_plot.title('Ts Graph for Water')
    ts_plot.xlabel(r's $[{kJ}/{kg K}]$')
    ts_plot.ylabel(r'T $[K]$')
    ts_plot.grid()
    ts_plot.show()

The following two examples show how the :class:`matplotlib.pyplot` functions
and :class:`matplotlib.pyplot.axes` functions can also be used along side
the :py:class:`CoolProp.Plots.Plots.PropsPlot` class

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ph_plot = PropsPlot('Water', 'Ph')
    ax = ph_plot.axis
    ax.set_yscale('log')
    ax.text(400, 5500, 'Saturated Liquid', fontsize=15, rotation=40)
    ax.text(2700, 3500, 'Saturated Vapour', fontsize=15, rotation=-100)
    ph_plot.show()

.. plot::
    :include-source:

    from matplotlib import pyplot
    from CoolProp.Plots import PropsPlot

    ref_fluid = 'R600a'
    fig = pyplot.figure(1, figsize=(10, 10), dpi=100)
    for i, gtype in enumerate(['PT', 'PD', 'PS', 'PH', 'TD', 'TS', 'HS']):
        ax = pyplot.subplot(4, 2, i+1)
        if gtype.startswith('P'):
            ax.set_yscale('log')
        props_plot = PropsPlot(ref_fluid, gtype, axis=ax)
        props_plot.title(gtype)
        props_plot._draw_graph()
    pyplot.tight_layout()
    pyplot.show()

