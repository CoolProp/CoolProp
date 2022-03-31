# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

import numpy as np
import warnings

import CoolProp
from CoolProp.Plots.Common import IsoLine, BasePlot, interpolate_values_1d
from CoolProp.Plots.SimpleCycles import StateContainer


class PropertyPlot(BasePlot):
    def __init__(self, fluid_name, graph_type, **kwargs):
        """
        Create graph for the specified fluid properties

        Parameters
        ----------
        fluid_name : string or AbstractState
            The name of the fluid to be plotted or a state instance
        graph_type : string
            The graph type to be plotted, like \"PH\" or \"TS\"
        axis : :func:`matplotlib.pyplot.gca()`, Optional
            The current axis system to be plotted to.
            Default: create a new axis system
        fig : :func:`matplotlib.pyplot.figure()`, Optional
            The current figure to be plotted to.
            Default: create a new figure
        unit_system : string, ['EUR','KSI','SI']
            Select the units used for the plotting.  'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
        tp_limits : string, ['NONE','DEF','ACHP','ORC']
            Select the limits in T and p.
        reciprocal_density : bool
            NOT IMPLEMENTED: If True, 1/rho will be plotted instead of rho

        Examples
        --------
        >>> from CoolProp.Plots import PropertyPlot
        >>> plot = PropertyPlot('HEOS::Water', 'TS')
        >>> plot.calc_isolines()
        >>> plot.show()

        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> plot = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR', tp_limits='ACHP')
        >>> plot.calc_isolines(CoolProp.iQ, num=11)
        >>> plot.calc_isolines(CoolProp.iT, num=25)
        >>> plot.calc_isolines(CoolProp.iSmass, num=15)
        >>> plot.show()

        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> plot = PropertyPlot('HEOS::R245fa', 'TS', unit_system='EUR', tp_limits='ORC')
        >>> plot.calc_isolines(CoolProp.iQ, num=11)
        >>> plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=True)
        >>> plot.draw()
        >>> plot.isolines.clear()
        >>> plot.props[CoolProp.iP]['color'] = 'green'
        >>> plot.props[CoolProp.iP]['lw'] = '0.5'
        >>> plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=False)
        >>> plot.show()

        .. note::

            See the online documentation for a list of the available fluids and
            graph types
        """
        super(PropertyPlot, self).__init__(fluid_name, graph_type, **kwargs)
        self._isolines = {}
        #self._plines = {}
        #self._ppoints = {}
        self.get_axis_limits()
        self._plot_default_annotations()

    @property
    def isolines(self): return self._isolines
    # @property
    #def plines(self): return self._plines
    # @property
    #def ppoints(self): return self._ppoints

    def show(self):
        self.draw()
        super(PropertyPlot, self).show()

    def savefig(self, *args, **kwargs):
        self.draw()
        super(PropertyPlot, self).savefig(*args, **kwargs)

    def _plotRound(self, values):
        """
        A function round an array-like object while maintaining the
        amount of entries. This is needed for the isolines since we
        want the labels to look pretty (=rounding), but we do not
        know the spacing of the lines. A fixed number of digits after
        rounding might lead to reduced array size.
        """
        inVal = np.unique(np.sort(np.array(values)))
        output = inVal[1:] * 0.0
        digits = -1
        limit = 10
        lim = inVal * 0.0 + 10
        # remove less from the numbers until same length,
        # more than 10 significant digits does not really
        # make sense, does it?
        while len(inVal) > len(output) and digits < limit:
            digits += 1
            val = (np.around(np.log10(np.abs(inVal))) * -1) + digits + 1
            val = np.where(val < lim, val, lim)
            val = np.where(val > -lim, val, -lim)
            output = np.zeros(inVal.shape)
            for i in range(len(inVal)):
                output[i] = np.around(inVal[i], decimals=int(val[i]))
            output = np.unique(output)
        return output

    def calc_isolines(self, iso_type=None, iso_range=None, num=15, rounding=False, points=250):
        """Calculate lines with constant values of type 'iso_type' in terms of x and y as
        defined by the plot object. 'iso_range' either is a collection of values or
        simply the minimum and maximum value between which 'num' lines get calculated.
        The 'rounding' parameter can be used to generate prettier labels if needed.
        """

        if iso_type is None or iso_type == 'all':
            for i_type in IsoLine.XY_SWITCH:
                if IsoLine.XY_SWITCH[i_type].get(self.y_index * 10 + self.x_index, None) is not None:
                    self.calc_isolines(i_type, None, num, rounding, points)
            return

        if iso_range is None:
            if iso_type is CoolProp.iQ:
                iso_range = [0.0, 1.0]
            else:
                limits = self.get_axis_limits(iso_type, CoolProp.iT)
                iso_range = [limits[0], limits[1]]

        if len(iso_range) <= 1 and num != 1:
            raise ValueError('You have to provide two values for the iso_range, {0} is not valid.'.format(iso_range))

        if len(iso_range) == 2 and (num is None or num < 2):
            raise ValueError('Please specify the number of isoline you want e.g. num=10.')

        iso_range = np.sort(np.unique(iso_range))
        # Generate iso ranges
        if len(iso_range) == 2:
            iso_range = self.generate_ranges(iso_type, iso_range[0], iso_range[1], num)
        if rounding:
            iso_range = self._plotRound(iso_range)

        # Limits are already in SI units
        limits = self._get_axis_limits()

        ixrange = self.generate_ranges(self._x_index, limits[0], limits[1], points)
        iyrange = self.generate_ranges(self._y_index, limits[2], limits[3], points)

        dim = self._system[iso_type]

        lines = self.isolines.get(iso_type, [])
        for i in range(num):
            lines.append(IsoLine(iso_type, self._x_index, self._y_index, value=dim.to_SI(iso_range[i]), state=self._state))
            lines[-1].calc_range(ixrange, iyrange)
            lines[-1].sanitize_data()
        self.isolines[iso_type] = lines
        return

    def draw_isolines(self):
        dimx = self._system[self._x_index]
        dimy = self._system[self._y_index]

        sat_props = self.props[CoolProp.iQ].copy()
        if 'lw' in sat_props: sat_props['lw'] *= 2.0
        else: sat_props['lw'] = 1.0
        if 'alpha' in sat_props: min([sat_props['alpha'] * 2.0, 1.0])
        else: sat_props['alpha'] = 1.0

        for i in self.isolines:
            props = self.props[i]
            dew = None; bub = None
            xcrit = None; ycrit = None
            if i == CoolProp.iQ:
                for line in self.isolines[i]:
                    if line.value == 0.0: bub = line
                    elif line.value == 1.0: dew = line
                if dew is not None and bub is not None:
                    xmin, xmax, ymin, ymax = self.get_axis_limits()
                    xmin = dimx.to_SI(xmin)
                    xmax = dimx.to_SI(xmax)
                    ymin = dimy.to_SI(ymin)
                    ymax = dimy.to_SI(ymax)
                    dx = xmax - xmin
                    dy = ymax - ymin
                    dew_filter = np.logical_and(np.isfinite(dew.x), np.isfinite(dew.y))
                    #dew_filter = np.logical_and(dew_filter,dew.x>dew.x[-1])
                    stp = min([dew_filter.size, 10])
                    dew_filter[0:-stp] = False
                    bub_filter = np.logical_and(np.isfinite(bub.x), np.isfinite(bub.y))

                    if self._x_index == CoolProp.iP or self._x_index == CoolProp.iDmass:
                        filter_x = lambda x: np.log10(x)
                    else:
                        filter_x = lambda x: x
                    if self._y_index == CoolProp.iP or self._y_index == CoolProp.iDmass:
                        filter_y = lambda y: np.log10(y)
                    else:
                        filter_y = lambda y: y

                    if (  # (filter_x(dew.x[dew_filter][-1])-filter_x(bub.x[bub_filter][-1])) > 0.010*filter_x(dx) and
                        (filter_x(dew.x[dew_filter][-1]) - filter_x(bub.x[bub_filter][-1])) < 0.050 * filter_x(dx) or
                        (filter_y(dew.y[dew_filter][-1]) - filter_y(bub.y[bub_filter][-1])) < 0.010 * filter_y(dy)):
                        x = np.linspace(bub.x[bub_filter][-1], dew.x[dew_filter][-1], 11)
                        y = interpolate_values_1d(
                          np.append(bub.x[bub_filter], dew.x[dew_filter][::-1]),
                          np.append(bub.y[bub_filter], dew.y[dew_filter][::-1]),
                          x_points=x,
                          kind='cubic')
                        self.axis.plot(dimx.from_SI(x), dimy.from_SI(y), **sat_props)
                        warnings.warn("Detected an incomplete phase envelope, fixing it numerically.")
                        xcrit = x[5]; ycrit = y[5]
                        #Tcrit = self.state.trivial_keyed_output(CoolProp.iT_critical)
                        #Dcrit = self.state.trivial_keyed_output(CoolProp.irhomass_critical)
                        # try:
                        #    self.state.update(CoolProp.DmassT_INPUTS, Dcrit, Tcrit)
                        #    xcrit = self.state.keyed_output(self._x_index)
                        #    ycrit = self.state.keyed_output(self._y_index)
                        # except:
                        #    xcrit = x[5]; ycrit = y[5]
                        #    pass
                        #self.axis.plot(dimx.from_SI(np.array([bub.x[bub_filter][-1], dew.x[dew_filter][-1]])),dimy.from_SI(np.array([bub.y[bub_filter][-1], dew.y[dew_filter][-1]])),'o')
            for line in self.isolines[i]:
                if line.i_index == CoolProp.iQ:
                    if line.value == 0.0 or line.value == 1.0:
                        self.axis.plot(dimx.from_SI(line.x), dimy.from_SI(line.y), **sat_props)
                    else:
                        if xcrit is not None and ycrit is not None:
                            self.axis.plot(dimx.from_SI(np.append(line.x, xcrit)), dimy.from_SI(np.append(line.y, ycrit)), **props)
                            # try:
                            #    x = np.append(line.x,[xcrit])
                            #    y = np.append(line.y,[ycrit])
                            #    fltr = np.logical_and(np.isfinite(x),np.isfinite(y))
                            #    f = interp1d(x[fltr][-3:],y[fltr][-3:],kind='linear') # could also be quadratic
                            #    x = np.linspace(x[fltr][-2], x[fltr][-1], 5)
                            #    y = f(x)
                            #    #f = interp1d(y[fltr][-5:],x[fltr][-5:],kind='cubic')
                            #    #y = np.linspace(y[fltr][-2], y[fltr][-1], 5)
                            #    #x = f(y)
                            #    self.axis.plot(dimx.from_SI(np.append(line.x,x)),dimy.from_SI(np.append(line.y,y)),**props)
                            # except:
                            #    self.axis.plot(dimx.from_SI(np.append(line.x,xcrit)),dimy.from_SI(np.append(line.y,ycrit)),**props)
                            #    pass
                else:
                    self.axis.plot(dimx.from_SI(line.x), dimy.from_SI(line.y), **props)

    def draw(self):
        self.get_axis_limits()
        self.draw_isolines()

    # def label_isolines(self, dx=0.075, dy=0.100):
    #    [xmin, xmax, ymin, ymax] = self.get_axis_limits()
    #    for i in self.isolines:
    #         for line in self.isolines[i]:
    #             if self.get_x_y_dydx(xv, yv, x)

    def draw_process(self, statecontainer, points=None, line_opts=None):
        """ Draw process or cycle from x and y values in axis units

        Parameters
        ----------
        statecontainer : CoolProp.Plots.SimpleCycles.StateContainer()
            A state container object that contains all the information required to draw the process.
            Note that points that appear several times get added to a special of highlighted points.
        line_opts : dict
            Line options (please see :func:`matplotlib.pyplot.plot`), optional
            Use this parameter to pass a label for the legend.

        Examples
        --------
        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> pp = PropertyPlot('HEOS::Water', 'TS', unit_system='EUR')
        >>> pp.calc_isolines(CoolProp.iP        )
        >>> pp.calc_isolines(CoolProp.iHmass    )
        >>> pp.calc_isolines(CoolProp.iQ, num=11)
        >>> cycle = SimpleRankineCycle('HEOS::Water', 'TS', unit_system='EUR')
        >>> T0 = 300
        >>> pp.state.update(CoolProp.QT_INPUTS,0.0,T0+15)
        >>> p0 = pp.state.keyed_output(CoolProp.iP)
        >>> T2 = 700
        >>> pp.state.update(CoolProp.QT_INPUTS,1.0,T2-150)
        >>> p2 = pp.state.keyed_output(CoolProp.iP)
        >>> cycle.simple_solve(T0, p0, T2, p2, 0.7, 0.8, SI=True)
        >>> cycle.steps = 50
        >>> sc = cycle.get_state_changes()
        >>> pp.draw_process(sc)
        >>> # The same calculation can be carried out in another unit system:
        >>> cycle.simple_solve(T0-273.15-10, p0/1e5, T2-273.15+50, p2/1e5-5, 0.7, 0.8, SI=False)
        >>> sc2 = cycle.get_state_changes()
        >>> pp.draw_process(sc2, line_opts={'color':'blue', 'lw':1.5})
        >>> pp.show()

        """
        warnings.warn("You called the function \"draw_process\", which is not tested.", UserWarning)

        # Default values
        line_opts = line_opts or {'color': 'r', 'lw': 1.5}

        dimx = self.system[self.x_index]
        dimy = self.system[self.y_index]

        marker = line_opts.pop('marker', 'o')
        style = line_opts.pop('linestyle', 'solid')
        style = line_opts.pop('ls', style)

        if points is None: points = StateContainer()

        xdata = []
        ydata = []
        old = statecontainer[len(statecontainer) - 1]
        for i in statecontainer:
            point = statecontainer[i]
            if point == old:
                points.append(point)
                old = point
                continue
            xdata.append(point[self.x_index])
            ydata.append(point[self.y_index])
            old = point
        xdata = dimx.from_SI(np.asarray(xdata))
        ydata = dimy.from_SI(np.asarray(ydata))
        self.axis.plot(xdata, ydata, marker='None', linestyle=style, **line_opts)

        xdata = np.empty(len(points))
        ydata = np.empty(len(points))
        for i in points:
            point = points[i]
            xdata[i] = point[self.x_index]
            ydata[i] = point[self.y_index]
        xdata = dimx.from_SI(np.asarray(xdata))
        ydata = dimy.from_SI(np.asarray(ydata))
        line_opts['label'] = ''
        self.axis.plot(xdata, ydata, marker=marker, linestyle='None', **line_opts)


def InlineLabel(xv, yv, x=None, y=None, axis=None, fig=None):
    warnings.warn("You called the deprecated function \"InlineLabel\", use \"BasePlot.inline_label\".", DeprecationWarning)
    plot = PropertyPlot("water", "TS", figure=fig, axis=axis)
    return plot.inline_label(xv, yv, x, y)


class PropsPlot(PropertyPlot):
    def __init__(self, fluid_name, graph_type, units='KSI', reciprocal_density=False, **kwargs):
        super(PropsPlot, self).__init__(fluid_name, graph_type, unit_system=units, reciprocal_density=reciprocal_density, **kwargs)
        warnings.warn("You called the deprecated class \"PropsPlot\", use \"PropertyPlot\".", DeprecationWarning)


if __name__ == "__main__":
    plot = PropertyPlot('HEOS::n-Pentane', 'PD', unit_system='EUR')  # , reciprocal_density=True)
    plot.calc_isolines(CoolProp.iT)
    plot.calc_isolines(CoolProp.iQ, num=11)
    # plot.calc_isolines(CoolProp.iSmass)
    # plot.calc_isolines(CoolProp.iHmass)
    plot.show()
