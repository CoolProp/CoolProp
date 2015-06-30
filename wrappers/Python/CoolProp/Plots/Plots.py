# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

import numpy, matplotlib, matplotlib.pyplot, math, re
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP

from scipy import interpolate
from scipy.spatial.kdtree import KDTree
import warnings
from CoolProp.Plots.Common import IsoLine,BasePlot
import CoolProp
import sys



class PropertyPlot(BasePlot):
    def __init__(self, fluid_name, graph_type, units = 'KSI', **kwargs):
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
        units : string, ['EUR','KSI','SI']
            Select the units used for the plotting.  'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
        reciprocal_density : bool
            NOT IMPLEMENTED: If True, 1/rho will be plotted instead of rho 

        Examples
        --------
        >>> from CoolProp.Plots import PropertyPlot
        >>> plot = PropertyPlot('Water', 'Ts')
        >>> plot.show()

        >>> plot = PropertyPlot('HEOS::n-Pentane', 'ph')
        >>> plot.calc_isolines(CoolProp.iQ,[0.0,1.0],num=11)
        >>> Ts_lim = plot.get_axis_limits(CoolProp.iT, CoolProp.iSmass)
        >>> plot.calc_isolines(CoolProp.iT,Ts_lim[0:2])
        >>> plot.calc_isolines(CoolProp.iSmass,Ts_lim[2:4])
        >>> plot.savefig('pentane_ph.pdf')

        .. note::

            See the online documentation for a list of the available fluids and
            graph types
        """
        super(PropertyPlot, self).__init__(fluid_name, graph_type, unit_system=units, **kwargs)
        self._isolines = {} 
        #self._plines = {}
        #self._ppoints = {}
        self.get_axis_limits()
        self._plot_default_annotations()
        
    @property
    def isolines(self): return self._isolines
    #@property
    #def plines(self): return self._plines
    #@property
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
        
    def calc_isolines(self, iso_type, iso_range, num=15, rounding=False, points=200):
        """Calculate lines with constant values of type 'iso_type' in terms of x and y as
        defined by the plot object. 'iso_range' either is a collection of values or 
        simply the minimum and maximum value between which 'num' lines get calculated.
        The 'rounding' parameter can be used to generate prettier labels if needed.
        """
        
        if iso_range is None or (len(iso_range) == 1 and num != 1):
            raise ValueError('Automatic interval detection for isoline \
                              boundaries is not supported yet, use the \
                              iso_range=[min, max] parameter.')
 
        if len(iso_range) == 2 and num is None:
            raise ValueError('Please specify the number of isoline you want \
                              e.g. num=10')

        if iso_type == 'all':
            for i_type in IsoLine.XY_SWITCH:
                if IsoLine.XY_SWITCH[i_type].get(self.y_index*10+self.x_index,None) is not None:
                    # TODO implement the automatic interval detection.
                    limits = self._get_axis_limits(i_type, CoolProp.iT)
                    self.calc_isolines(i_type, [limits[0],limits[1]], num, rounding, points)
            return 
                    
        iso_range = numpy.sort(numpy.unique(iso_range))
        # Generate iso ranges
        if len(iso_range) == 2:
            iso_range = self.generate_ranges(iso_type, iso_range[0], iso_range[1], num)
        if rounding:
            iso_range = self._plotRound(iso_range)
        
        # Limits are alreadyin SI units
        limits = self._get_axis_limits()
        
        ixrange = self.generate_ranges(self._x_index,limits[0],limits[1],points)
        iyrange = self.generate_ranges(self._y_index,limits[2],limits[3],points)
        
        dim = self._system.dimensions[iso_type]
        
        lines  = self.isolines.get(iso_type, [])
        for i in range(num):
            lines.append(IsoLine(iso_type,self._x_index,self._y_index, value=dim.to_SI(iso_range[i]), state=self._state))
            lines[-1].calc_range(ixrange,iyrange)
            lines[-1].sanitize_data()
        self.isolines[iso_type] = lines 
        return 
    
    
    def draw_isolines(self):
        for i in self.isolines:
            props = self.props[i]
            dimx = self._system.dimensions[self._x_index]
            dimy = self._system.dimensions[self._y_index]
            for line in self.isolines[i]:
                if line.i_index == CoolProp.iQ and \
                  (line.value == 0.0 or line.value == 1.0):
                    plot_props = props.copy()
                    if 'lw' in plot_props: plot_props['lw'] *= 2.0
                    else: plot_props['lw'] = 1.0
                    if 'alpha' in plot_props: plot_props['alpha'] *= 2.0
                    else: plot_props['alpha'] = 1.0
                else:
                    plot_props = props
                self.axis.plot(dimx.from_SI(line.x),dimy.from_SI(line.y),**plot_props)
    
    def draw(self):
        self.draw_isolines()
        
    #def label_isolines(self, dx=0.075, dy=0.100):
    #    [xmin, xmax, ymin, ymax] = self.get_axis_limits()
    #    for i in self.isolines:
    #         for line in self.isolines[i]:
    #             if self.get_x_y_dydx(xv, yv, x)
             
         
                
                
    def draw_process(self, states, iso_types=None, line_opts={'color' : 'r', 'lw' : 1.5}):
        """ Draw process or cycle from x and y values in axis units

        Parameters
        ----------
        states : list of (x,y) tuples, required
        iso_types : list 
            isobars that should be used to illustrate the processes, one element less than states, optional
        line_opts : dict
            Line options (please see :func:`matplotlib.pyplot.plot`), optional
        """
        warnings.warn("You called the function \"draw_process\", which is not tested.",UserWarning)

        # plot above other lines
        line_opts['zorder'] = 10
        
        if iso_types is not None and len(states)!=len(iso_types)+1:
            raise ValueError("If you specifiy the isotypes, they have to have the length of the state list - 1.")
        
        X = []
        Y = []

        for i in range(len(states)):
            if i == 0: continue
            (x2, y2) = states[i]
            (x1, y1) = states[i-1]
            
            iso_type = None
            if iso_types is not None and iso_types[i-1] is not None:
                iso_type = self._get_index(iso_types[i-1])
            else: # TODO: detect it!
                iso_type = None
                
            iso_line = None 
            if iso_type is not None:
                switch = IsoLine.XY_SWITCH[iso_type].get(self.y_index*10+self.x_index,None)
                if switch is not None:
                    try: 
                        dimx = self.system.dimensions[self.x_index]
                        dimy = self.system.dimensions[self.y_index]
                        dimi = self.system.dimensions[iso_type]
                        pair, out1, out2 = CP.generate_update_pair(self.x_index,dimx.to_SI(x1),self.y_index,dimy.to_SI(y1))
                        self.state.update(pair, out1, out2)
                        i_val1 = self.state.keyed_output(iso_type)
                        pair, out1, out2 = CP.generate_update_pair(self.x_index,dimx.to_SI(x2),self.y_index,dimy.to_SI(y2))
                        self.state.update(pair, out1, out2)
                        i_val2 = self.state.keyed_output(iso_type)
                        i_val = dimi.from_SI((i_val1 + i_val2)/2.0) 
                        self.calc_isolines(iso_type, [i_val], num=1)
                        iso_line = self.isolines[iso_type].pop()
                        idx1 = numpy.argmin(numpy.abs(iso_line.x - x1))
                        idx2 = numpy.argmin(numpy.abs(iso_line.x - x2))
                        if idx1>idx2:
                            iso_line.x = iso_line.x[idx2+1:idx1]
                            iso_line.y = iso_line.y[idx2+1:idx1]
                        else:
                            iso_line.x = iso_line.x[idx1+1:idx2]
                            iso_line.y = iso_line.y[idx1+1:idx2]
                    except Exception as e:
                        warnings.warn(
                          "There was a problem with the isolines: {0:s}".format(str(e)),
                          UserWarning)
            
            if iso_line is None:
                iso_line = IsoLine(CoolProp.iT, self.x_index, self.y_index) # Just a dummy
                iso_line.x = [x1,x2]
                iso_line.y = [y1,y2]
            
            self.axis.plot(iso_line.x,iso_line.y,**line_opts)

def InlineLabel(xv,yv,x=None,y=None,axis=None,fig=None):
    warnings.warn("You called the deprecated function \"InlineLabel\", use \"BasePlot.inline_label\".",DeprecationWarning)
    plot = PropertyPlot("water","TS",figure=fig,axis=axis)
    return plot.inline_label(xv,yv,x,y)

class PropsPlot(PropertyPlot):
    def __init__(self, fluid_name, graph_type, units = 'KSI', reciprocal_density = False, **kwargs):
        super(PropsPlot, self).__init__(fluid_name, graph_type, units=units, reciprocal_density=reciprocal_density, **kwargs)
        warnings.warn("You called the deprecated class \"PropsPlot\", use \"PropertyPlot\".",DeprecationWarning)


if __name__ == "__main__":
    plot = PropertyPlot('HEOS::n-Pentane', 'PH', units='EUR')
    Ts = plot.get_axis_limits(CoolProp.iT, CoolProp.iSmass)
    TD = plot.get_axis_limits(CoolProp.iT, CoolProp.iDmass)
    plot.calc_isolines(CoolProp.iT,     Ts[0:2])
    plot.calc_isolines(CoolProp.iQ,     [0.0,1.0], num=11)
    plot.calc_isolines(CoolProp.iSmass, Ts[2:4])
    plot.calc_isolines(CoolProp.iDmass, TD[2:4])
#     plot.calc_isolines('all', None)
    plot.draw_isolines()
    #
    Tcrit = plot.state.trivial_keyed_output(CoolProp.iT_critical)
    Dcrit = plot.state.trivial_keyed_output(CoolProp.irhomass_critical)
    plot.state.update(CoolProp.DmassT_INPUTS, Dcrit, Tcrit)
    p1 = plot.state.keyed_output(CoolProp.iP)/1e5 / 2.00
    h1 = plot.state.keyed_output(CoolProp.iHmass)/1e3 * 1.25
    p2 = plot.state.keyed_output(CoolProp.iP)/1e5 / 2.25 
    h2 = plot.state.keyed_output(CoolProp.iHmass)/1e3 * 1.50
    plot.draw_process(zip([h1,h2],[p1,p2]))
    #
    
    #
    plot.savefig("Plots.pdf")
    #for i in plot.isolines:
    #    print(plot.isolines[i][0].x,plot.isolines[i][0].y)
