# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

import numpy, matplotlib, matplotlib.pyplot, math, re
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP

from .Common import BasePlot
from scipy import interpolate
from scipy.spatial.kdtree import KDTree

class IsoLine(object):
    def __init__(self):
        self.DEBUG = False

        # direct geometry
        self.X     = None #
        self.Y     = None #
        self.type  = None #
        self.value = None #
        self.unit  = None #
        self.opts  = None #


def InlineLabel(xv,yv,x = None, y= None, axis = None, fig = None):
    """
    This will give the coordinates and rotation required to align a label with
    a line on a plot
    """

    def ToPixelCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin

        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=width*(Fxmax-Fxmin)
        DELTAY_fig=height*(Fymax-Fymin)

        #Convert coords to pixels
        x=(xv-Axmin)/DELTAX_axis*DELTAX_fig+Fxmin
        y=(yv-Aymin)/DELTAY_axis*DELTAY_fig+Fymin

        return x,y

    def ToDataCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin

        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=(Fxmax-Fxmin)*width
        DELTAY_fig=(Fymax-Fymin)*height

        #Convert back to measurements
        x=(xv-Fxmin)/DELTAX_fig*DELTAX_axis+Axmin
        y=(yv-Fymin)/DELTAY_fig*DELTAY_axis+Aymin

        return x,y

    def get_x_y_dydx(xv,yv,x):
        """Get x and y coordinates and the linear interpolation derivative"""
        # Old implementation:
        ##Get the rotation angle
        #f = interp1d(xv, yv)
        #y = f(x)
        #h = 0.00001*x
        #dy_dx = (f(x+h)-f(x-h))/(2*h)
        #return x,y,dy_dx
        if len(xv)==len(yv)>1: # assure same length
            if len(xv)==len(yv)==2: # only two points
                if numpy.min(xv)<x<numpy.max(xv):
                    dx    = xv[1] - xv[0]
                    dy    = yv[1] - yv[0]
                    dydx  = dy/dx
                    y     = yv[0] + dydx * (x-xv[0])
                    return x,y,dydx
                else:
                    raise ValueError("Your coordinate has to be between the input values.")
            else:
                limit = 1e-10                    # avoid hitting a point directly
                diff  = numpy.array(xv)-x        # get differences
                index = numpy.argmin(diff*diff)  # nearest neighbour
                if (xv[index]<x<xv[index+1]      # nearest below, positive inclination
                  or xv[index]>x>xv[index+1]):   # nearest above, negative inclination
                    if diff[index]<limit:
                        index = [index-1,index+1]
                    else:
                        index = [index,  index+1]
                elif (xv[index-1]<x<xv[index]    # nearest above, positive inclination
                  or xv[index-1]>x>xv[index]):   # nearest below, negative inclination
                    if diff[index]<limit:
                        index = [index-1,index+1]
                    else:
                        index = [index-1,index]
                xvnew = xv[index]
                yvnew = yv[index]
                return get_x_y_dydx(xvnew,yvnew,x) # Allow for a single recursion
        else:
            raise ValueError("You have to provide the same amount of x- and y-pairs with at least two entries each.")

    if axis is None:
        axis=matplotlib.pyplot.gca()

    if fig is None:
        fig=matplotlib.pyplot.gcf()

    if y is None and x is not None:
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        #x is provided but y isn't
        (x,trash)=ToPixelCoords(x,trash,axis,fig)

        #Get the rotation angle and y-value
        x,y,dy_dx = get_x_y_dydx(xv,yv,x)
        rot = numpy.arctan(dy_dx)/numpy.pi*180.

    elif x is None and y is not None:
        #y is provided, but x isn't

        _xv = xv[::-1]
        _yv = yv[::-1]
        #Find x by interpolation
        x = interp1d(yv, xv)(y)
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        (x,trash)=ToPixelCoords(x,trash,axis,fig)

        #Get the rotation angle and y-value
        x,y,dy_dx = get_x_y_dydx(xv,yv,x)
        rot = numpy.arctan(dy_dx)/numpy.pi*180.

    (x,y)=ToDataCoords(x,y,axis,fig)
    return (x,y,rot)


def drawLines(Ref,lines,axis,plt_kwargs=None):
    """
    Just an internal method to systematically plot values from
    the generated 'line' dicts, method is able to cover the whole
    saturation curve. Closes the gap at the critical point and
    adds a marker between the two last points of bubble and
    dew line if they reach up to critical point.
    Returns an array of line objects that can be used to change
    the colour or style afterwards.
    """
    if not plt_kwargs is None:
        for line in lines:
            line['opts'] = plt_kwargs
    plottedLines = []
    if len(lines)==2 and (
      'q' in str(lines[0]['type']).lower() and 'q' in str(lines[1]['type']).lower()
      ) and (
      ( 0 == lines[0]['value'] and 1 == lines[1]['type'] ) or ( 1 == lines[0]['value'] and 0 == lines[1]['type'] ) ):
        # We plot the saturation curve
        bubble = lines[0]
        dew = lines[1]
        line, = axis.plot(bubble['x'],bubble['y'],**bubble['opts'])
        plottedLines.extend([line])
        line, = axis.plot(dew['x'], dew['y'], **dew['opts'])
        plottedLines.extend([line])
        # Do we need to test if this is T or p?
        Tmax = min(bubble['kmax'],dew['kmax'])
        if Tmax>CP.PropsSI(Ref,'Tcrit')-2e-5:
            axis.plot(numpy.r_[bubble['x'][-1],dew['x'][-1]],numpy.r_[bubble['y'][-1],dew['y'][-1]],**bubble['opts'])
            #axis.plot((bubble['x'][-1]+dew['x'][-1])/2.,(bubble['y'][-1]+dew['y'][-1])/2.,'o',color='Tomato')
    else:
        for line in lines:
            line, = axis.plot(line['x'],line['y'],**line['opts'])
            plottedLines.extend([line])

    return plottedLines


class IsoLines(BasePlot):
    def __init__(self, fluid_ref, graph_type, iso_type, unit_system='SI', **kwargs):
        BasePlot.__init__(self, fluid_ref, graph_type, unit_system=unit_system,**kwargs)

        if not isinstance(iso_type, str):
            raise TypeError("Invalid iso_type input, expected a string")

        iso_type = iso_type.upper()
        if iso_type not in self.COLOR_MAP.keys() and iso_type != 'Q':
            raise ValueError('This kind of isoline is not supported for a ' \
                             + str(graph_type) + \
                             ' plot. Please choose from '\
                             + str(self.COLOR_MAP.keys()) + ' or Q.')

        self.iso_type = iso_type

    def __set_axis_limits(self, swap_xy):
        """
        Generates limits for the axes in terms of x,y defined by 'plot'
        based on temperature and pressure.

        Returns a tuple containing ((xmin, xmax), (ymin, ymax))
        """
        # Get current axis limits, be sure to set those before drawing isolines
        # if no limits are set, use triple point and critical conditions
        X = [CP.PropsSI(self.graph_type[1],
                      'T', 1.5*CP.PropsSI(self.fluid_ref, 'Tcrit'),
                      'P', CP.PropsSI(self.fluid_ref, 'ptriple'),
                      self.fluid_ref),
             CP.PropsSI(self.graph_type[1],
                      'T', 1.1*CP.PropsSI(self.fluid_ref, 'Tmin'),
                      'P', 1.5*CP.PropsSI(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.PropsSI(self.graph_type[1],
                      'T', 1.5*CP.PropsSI(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.PropsSI(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.PropsSI(self.graph_type[1],
                      'T', 1.1*CP.PropsSI(self.fluid_ref, 'Tmin'),
                      'P', CP.PropsSI(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        Y = [CP.PropsSI(self.graph_type[0],
                      'T', 1.5*CP.PropsSI(self.fluid_ref, 'Tcrit'),
                      'P', CP.PropsSI(self.fluid_ref, 'ptriple'),
                       self.fluid_ref),
             CP.PropsSI(self.graph_type[0],
                      'T', 1.1*CP.PropsSI(self.fluid_ref, 'Tmin') ,
                      'P', 1.5*CP.PropsSI(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.PropsSI(self.graph_type[0],
                      'T', 1.1*CP.PropsSI(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.PropsSI(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.PropsSI(self.graph_type[0],
                      'T', 1.5*CP.PropsSI(self.fluid_ref, 'Tmin') ,
                      'P', CP.PropsSI(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        limits = [[min(X), max(X)], [min(Y), max(Y)]]
        if not self.axis.get_autoscalex_on():
            limits[0][0] = max([limits[0][0], min(self.axis.get_xlim())])
            limits[0][1] = min([limits[0][1], max(self.axis.get_xlim())])
            limits[1][0] = max([limits[1][0], min(self.axis.get_ylim())])
            limits[1][1] = min([limits[1][1], max(self.axis.get_ylim())])
            
        # Limits correction in case of KSI unit_system
        if self.unit_system == 'KSI':
            limits[0] = [l*self.KSI_SCALE_FACTOR[self.graph_type[1]] for l in limits[0]]
            limits[1] = [l*self.KSI_SCALE_FACTOR[self.graph_type[0]] for l in limits[1]]

        self.axis.set_xlim(limits[0])
        self.axis.set_ylim(limits[1])
        return limits

    def __plotRound(self, values):
        """
        A function round an array-like object while maintaining the
        amount of entries. This is needed for the isolines since we
        want the labels to look pretty (=rounding), but we do not
        know the spacing of the lines. A fixed number of digits after
        rounding might lead to reduced array size.
        """
        inVal   = numpy.unique(numpy.sort(numpy.array(values)))
        output  = inVal[1:] * 0.0
        digits  = -1
        limit   = 10
        lim     = inVal * 0.0 + 10
        # remove less from the numbers until same length,
        # more than 10 significant digits does not really
        # make sense, does it?
        while len(inVal) > len(output) and digits < limit:
            digits += 1
            val     = ( numpy.around(numpy.log10(numpy.abs(inVal))) * -1) + digits + 1
            val     = numpy.where(val < lim, val,  lim)
            val     = numpy.where(val >-lim, val, -lim)
            output  = numpy.zeros(inVal.shape)
            for i in range(len(inVal)):
                output[i] = numpy.around(inVal[i],decimals=int(val[i]))
            output = numpy.unique(output)
        return output

    def get_isolines(self, iso_range=[], num=None, rounding=False):
        """
        This is the core method to obtain lines in the dimensions defined
        by 'plot' that describe the behaviour of fluid 'Ref'. The constant
        value is determined by 'iName' and has the values of 'iValues'.

        'iValues' is an array-like object holding at least one element. Lines
        are calculated for every entry in 'iValues'. If the input 'num' is
        larger than the amount of entries in 'iValues', an internally defined
        pattern is used to calculate an appropriate line spacing between the maximum
        and minimum values provided in 'iValues'.

        Returns lines[num] - an array of dicts containing 'x' and 'y'
        coordinates for bubble and dew line. Additionally, the dict holds
        the keys 'label' and 'opts', those can be used for plotting as well.
        """
        if iso_range is None or (len(iso_range) == 1 and num != 1):
            raise ValueError('Automatic interval detection for isoline \
                              boundaries is not supported yet, use the \
                              iso_range=[min, max] parameter.')

        if len(iso_range) == 2 and num is None:
            raise ValueError('Please specify the number of isoline you want \
                              e.g. num=10')

        iso_range = numpy.sort(numpy.unique(iso_range))

        def generate_ranges(xmin, xmax, num):
            if self.iso_type in ['P', 'D']:
                return numpy.logspace(math.log(xmin, 2.),
                                      math.log(xmax, 2.),
                                      num=num,
                                      base=2.)
            return numpy.linspace(xmin, xmax, num=num)

        # Generate iso ranges
        if len(iso_range) == 2:
            iso_range = generate_ranges(iso_range[0], iso_range[1], num)
            #iso_range = plotRound(iso_range)
        #else:
        #    TODO: Automatic interval detection
        #    iVal = [CP.PropsSI(iName,'T',T_c[i],'D',rho_c[i],Ref) for i in range(len(T_c))]
        #    iVal = patterns[iName]([numpy.min(iVal),numpy.max(iVal),num])

        if rounding:
            iso_range = self.__plotRound(iso_range)

        switch_xy_map = {'D': ['TS', 'PH', 'PS'],
                         'S': ['PH', 'PD', 'PT'],
                         'T': ['PH', 'PS'],
                         'H': ['PD']}
        #TS: TD is defined, SD is not
        #PH: PD is defined, HD is not
        #PS: PD is defined, SD is not
        #PH: PS is more stable than HS
        #PD: PS is defined, DS is not
        #PT: PS is defined, TS is not
        #PH: PT is defined, HT is not
        #PS: PT is defined, ST is not
        #PD: PH is defined, DH is not

        iso_error_map = {'TD': ['S', 'H'],
                         'HS': ['T', 'D'],}

        switch_xy = False
        if self.iso_type in ['D', 'S', 'T', 'H']:
            if self.graph_type in switch_xy_map[self.iso_type]:
                switch_xy = True

        if self.graph_type in ['TD', 'HS']:
            if self.iso_type in iso_error_map[self.graph_type]:
                raise ValueError('You should not reach this point!')

        axis_limits = self.__set_axis_limits(switch_xy)
        req_prop = self.graph_type[0]
        prop2_name = self.graph_type[1]
        if switch_xy:
            axis_limits.reverse()
            req_prop = self.graph_type[1]
            prop2_name = self.graph_type[0]

        # Calculate the points
        if self.iso_type == 'Q':
            lines = self._get_sat_lines(x=iso_range)
            return lines

        # TODO: Determine saturation state if two phase region present
        x_range = numpy.linspace(axis_limits[0][0], axis_limits[0][1], 1000.)
        x_mesh = [x_range for i in iso_range]

        plot_data = self._get_fluid_data(req_prop,
                                         self.iso_type, iso_range,
                                         prop2_name, x_mesh)

        if switch_xy:
            plot_data = plot_data[::-1]

        lines = []
        for j in range(len(plot_data[0])):
            line = {
              'x': plot_data[0][j],
              'y': plot_data[1][j],
              # TODO
              'label': "", #_getIsoLineLabel(self.iso_type, iso_range[j]),
              'type': self.iso_type,
              'opts': {'color': self.COLOR_MAP[self.iso_type], 'lw':0.75, 'alpha':0.5 }
              }
            lines.append(line)

        return lines

    def draw_isolines(self, iso_range, num=None, rounding=False):
        """
        Draw lines with constant values of type 'which' in terms of x and y as
        defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value between
        which 'num' get drawn.

        There should also be helpful error messages...
        """
        if iso_range is None or (len(iso_range) == 1 and num != 1):
            raise ValueError('Automatic interval detection for isoline \
                              boundaries is not supported yet, use the \
                              iso_range=[min, max] parameter.')

        if len(iso_range) == 2 and num is None:
            raise ValueError('Please specify the number of isoline you want \
                              e.g. num=10')

        if self.iso_type == 'all':
            raise ValueError('Plotting all lines automatically is not \
                              supported, yet..')

        if self.iso_type != 'all':
            lines = self.get_isolines(iso_range, num, rounding)
            drawn_lines = drawLines(self.fluid_ref, lines, self.axis)
            self._plot_default_annotations()
            return drawn_lines
        #else:
        #    # TODO: assign limits to values automatically
        #    ll = _getIsoLineIds(plot)
        #    if not len(ll)==len(iValues):
        #        raise ValueError('Please provide a properly sized array of bounds.')
        #    for c,l in enumerate(ll):
        #        drawIsoLines(Ref, plot, l, iValues=iValues[c], num=num, axis=axis, fig=fig)


class PropsPlot(BasePlot):
    def __init__(self, fluid_name, graph_type, units = 'KSI', reciprocal_density = False, **kwargs):
        """
        Create graph for the specified fluid properties

        Parameters
        ----------
        fluid_ref : string
            The name of the fluid to be plotted
        graph_type : string
            The graph type to be plotted
        axis : :func:`matplotlib.pyplot.gca()`, Optional
            The current axis system to be plotted to.
            Default: create a new axis system
        fig : :func:`matplotlib.pyplot.figure()`, Optional
            The current figure to be plotted to.
            Default: create a new figure
        units : string, ['KSI','SI']
            Select the units used for the plotting.  'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
        reciprocal_density : bool
            If True, 1/rho will be plotted instead of rho

        Examples
        --------
        >>> from CoolProp.Plots import PropsPlot
        >>> plt = PropsPlot('Water', 'Ph')
        >>> plt.show()

        >>> plt = PropsPlot('n-Pentane', 'Ts')
        >>> plt.set_axis_limits([-0.5, 1.5, 300, 530])
        >>> plt.draw_isolines('Q', [0.1, 0.9])
        >>> plt.draw_isolines('P', [100, 2000])
        >>> plt.draw_isolines('D', [2, 600])
        >>> plt.show()

        .. note::

            See the online documentation for a list of the available fluids and
            graph types
        """
        BasePlot.__init__(self, fluid_name, graph_type, unit_system=units, **kwargs)

        self.smin = kwargs.get('smin', None)
        self.smax = kwargs.get('smax', None)
        
        self._draw_graph()

    def __draw_region_lines(self):
        lines = self._get_sat_lines(kind='T',
                                    smin=self.smin,
                                    smax=self.smax)
        drawLines(self.fluid_ref, lines, self.axis)

    def _draw_graph(self):
        self.__draw_region_lines()
        self._plot_default_annotations()

    def draw_isolines(self, iso_type, iso_range, num=10, rounding=False):
        iso_lines = IsoLines(self.fluid_ref,
                             self.graph_type,
                             iso_type, unit_system = self.unit_system,
                             axis=self.axis)
        iso_lines.draw_isolines(iso_range, num, rounding)


def Ts(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'Ts', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    plt._draw_graph()
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis


def Ph(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'Ph', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis


def Ps(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'Ps', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def PT(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'PT', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def Prho(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    """
    plt = PropsPlot(Ref, 'PD', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def Trho(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'TD', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def hs(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Deprecated.  Use :py:func:`CoolProps.Plots.PropsPlot`
    """
    plt = PropsPlot(Ref, 'hs', smin=Tmin, smax=Tmax, axis=axis, *args, **kwargs)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def drawIsoLines(Ref, plot, which, iValues=[], num=0, show=False, axis=None):
    """
    Draw lines with constant values of type 'which' in terms of x and y as
    defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value
    between which 'num' get drawn.

    :Note:
        :func:`CoolProps.Plots.drawIsoLines` will be depreciated in future
        releases and replaced with :func:`CoolProps.Plots.IsoLines`

    Parameters
    ----------
    Ref : str
        The given reference fluid
    plot : str
        The plot type used
    which : str
        The iso line type
    iValues : list
        The list of constant iso line values
    num : int, Optional
        The number of iso lines
        (Default: 0 - Use iValues list only)
    show : bool, Optional
        Show the current plot
        (Default: False)
    axis : :func:`matplotlib.pyplot.gca()`, Optional
        The current axis system to be plotted to.
        (Default: create a new axis system)

    Examples
    --------
    >>> from matplotlib import pyplot
    >>> from CoolProp.Plots import Ts, drawIsoLines
    >>>
    >>> Ref = 'n-Pentane'
    >>> ax = Ts(Ref)
    >>> ax.set_xlim([-0.5, 1.5])
    >>> ax.set_ylim([300, 530])
    >>> quality = drawIsoLines(Ref, 'Ts', 'Q', [0.3, 0.5, 0.7, 0.8], axis=ax)
    >>> isobars = drawIsoLines(Ref, 'Ts', 'P', [100, 2000], num=5, axis=ax)
    >>> isochores = drawIsoLines(Ref, 'Ts', 'D', [2, 600], num=7, axis=ax)
    >>> pyplot.show()
    """
    isolines = IsoLines(Ref, plot, which, axis=axis)
    lines = isolines.draw_isolines(iValues, num)
    if show:
        isolines.show()
    return lines
