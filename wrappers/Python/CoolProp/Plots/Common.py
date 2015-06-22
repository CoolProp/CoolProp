# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals

import matplotlib
import numpy

import CoolProp.CoolProp as CP
from abc import ABCMeta


class BaseQuantity(object):
    """A very basic property that can convert an input to and from a 
    given unit system, note that the conversion from SI units starts 
    with a multiplication. If you need to remove an offset, use the 
    off_SI property. 
    Examples with temperature:
    celsius = BaseQuantity(add_SI=-273.15)
    fahrenheit = BaseQuantity(add_SI=32.0, mul_SI=1.8, off_SI=-273.15)
    Examples with pressure:
    bar = BaseQuantity(mul_SI=1e-5)
    psi = BaseQuantity(mul_SI=0.000145037738)    
    """
    def __init__(self, add_SI=0.0, mul_SI=1.0, off_SI=0.0):
        self._add_SI = add_SI
        self._mul_SI = mul_SI
        self._off_SI = off_SI

    @property
    def add_SI(self): return self._add_SI
    @add_SI.setter
    def add_SI(self, value): self._add_SI = value
    @property
    def mul_SI(self): return self._mul_SI
    @mul_SI.setter
    def mul_SI(self, value): self._mul_SI = value
    @property
    def off_SI(self): return self._off_SI
    @off_SI.setter
    def off_SI(self, value): self._off_SI = value
    
    def from_SI(self, value): return ((value+self.off_SI)*self.mul_SI)+self.add_SI
    def to_SI(self, value): return (value-self.add_SI)/self.mul_SI-self.off_SI
    

class BaseDimension(BaseQuantity):
    """A dimension is a class that extends the BaseQuantity and adds a label, a symbol and a unit label"""
    def __init__(self, add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='', symbol='', unit=''):
        self._label = label
        self._symbol = symbol
        self._unit = unit
        super(BaseDimension, self).__init__(add_SI=add_SI, mul_SI=mul_SI, off_SI=off_SI)

    @property
    def label(self): return self._label
    @label.setter
    def label(self, value): self._label = value
    @property
    def symbol(self): return self._symbol
    @symbol.setter
    def symbol(self, value): self._symbol = value
    @property
    def unit(self): return self._unit
    @unit.setter
    def unit(self, value): self._unit = value


class UnitSystem(object):
    """A collection of dimensions for all the required quantities"""
    __metaclass__ = ABCMeta
    @property
    def D(self): return self._D
    @D.setter
    def D(self, value): self._D = value
    @property
    def H(self): return self._H
    @H.setter
    def H(self, value): self._H = value
    @property
    def P(self): return self._P
    @P.setter
    def P(self, value): self._P = value
    @property
    def S(self): return self._S
    @S.setter
    def S(self, value): self._S = value
    @property
    def T(self): return self._T
    @T.setter
    def T(self, value): self._T = value
    @property
    def U(self): return self._U
    @U.setter
    def U(self, value): self._U = value
    
    def dimensions(self): 
        return {
      'D' : self._D,
      'H' : self._H,
      'P' : self._P,
      'S' : self._S,
      'T' : self._T,
      'U' : self._U
    }


class SIunits(UnitSystem):
    def __init__(self):
        self._D = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Density',                  symbol=r'$\rho$', unit=r'kg/m$^3$')
        self._H = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Enthalpy',        symbol=r'$h$',    unit=r'J/kg')
        self._P = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Pressure',                 symbol=r'$p$',    unit=r'Pa')
        self._S = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Entropy',         symbol=r'$s$',    unit=r'J/kg/K')
        self._T = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Temperature',              symbol=r'$T$',    unit=r'K')
        self._U = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Internal Energy', symbol=r'$u$',    unit=r'J/kg')

class KSIunits(SIunits):
    def __init__(self):
        super(KSIunits, self).__init__()
        self.H.mul_SI=1e-3
        self.H.unit=r'kJ/kg'
        self.P.mul_SI=1e-3
        self.P.unit=r'kPa'
        self.S.mul_SI=1e-3
        self.S.unit=r'kJ/kg/K'
        self.U.mul_SI=1e-3
        self.U.unit=r'kJ/kg'

class EURunits(KSIunits):
    def __init__(self):
        super(EURunits, self).__init__()
        self.P.mul_SI=1e-5
        self.P.unit=r'bar'
        self.T.add_SI=-273.15
        self.T.unit=r'C'



if __name__ == "__main__":
    for syst in [SIunits(), KSIunits(), EURunits()]:
        print(syst.P.label)
        print(syst.P.to_SI(20))
        print(syst.T.label)
        print(syst.T.to_SI(20))
    
    



SMALL = 1E-5

class BasePlot(object):
    
    
    #TODO: Simplify / Consolidate dictionary maps
    AXIS_LABELS = {'KSI': {'T': ["Temperature", r"[K]"],
                           'P': ["Pressure", r"[kPa]"],
                           'S': ["Entropy", r"[kJ/kg/K]"],
                           'H': ["Enthalpy", r"[kJ/kg]"],
                           'U': ["Internal Energy", r"[kJ/kg]"],
                           'D': ["Density", r"[kg/m$^3$]"]
                          },
                    'SI': {'T': ["Temperature", r"[K]"],
                           'P': ["Pressure", r"[Pa]"],
                           'S': ["Entropy", r"[J/kg/K]"],
                           'H': ["Enthalpy", r"[J/kg]"],
                           'U': ["Internal Energy", r"[J/kg]"],
                           'D': ["Density", r"[kg/m$^3$]"]
                          }
                   }

    COLOR_MAP = {'T': 'Darkred',
                 'P': 'DarkCyan',
                 'H': 'DarkGreen',
                 'D': 'DarkBlue',
                 'S': 'DarkOrange',
                 'Q': 'black'}

    #: Scale factors to multiply SI units by in order to obtain kSI units
    KSI_SCALE_FACTOR = {'T' : 1.0,
                        'P' : 0.001,
                        'H' : 0.001,
                        'U' : 0.001,
                        'D' : 1,
                        'S' : 0.001,
                        'Q' : 1.0}
                      
    SYMBOL_MAP_KSI = {'T' : [r'$T = ', r'$ K'],
                      'P' : [r'$p = ', r'$ kPa'],
                      'H' : [r'$h = ', r'$ kJ/kg'],
                      'U' : [r'$h = ', r'$ kJ/kg'],
                      'D' : [r'$\rho = ', r'$ kg/m$^3$'],
                      'S' : [r'$s = ', r'$ kJ/kg-K'],
                      'Q' : [r'$x = ', r'$']}
                  
    SYMBOL_MAP_SI = {'T' : [r'$T = ', r'$ K'],
                     'P' : [r'$p = ', r'$ Pa'],
                     'H' : [r'$h = ', r'$ J/kg'],
                     'U' : [r'$h = ', r'$ J/kg'],
                     'D' : [r'$\rho = ', r'$ kg/m$^3$'],
                     'S' : [r'$s = ', r'$ J/kg-K'],
                     'Q' : [r'$x = ', r'$']}

    LINE_IDS = {'TS': ['P', 'D'], #'H'],
                'PH': ['S', 'T', 'D'],
                'HS': ['P'], #'T', 'D'],
                'PS': ['H', 'T', 'D'],
                'PD': ['T', 'S', 'H'],
                'TD': ['P'], #'S', 'H'],
                'PT': ['D', 'P', 'S'],
                'PU': []}

    def __init__(self, fluid_ref, graph_type, unit_system = 'KSI', **kwargs):
        if not isinstance(graph_type, str):
            raise TypeError("Invalid graph_type input, expected a string")

        graph_type = graph_type.upper()
        if len(graph_type) >= 2 and graph_type[1:len(graph_type)] == 'RHO':
            graph_type = graph_type[0] + graph_type[1:len(graph_type)]

        if graph_type.upper() not in self.LINE_IDS.keys():
            raise ValueError(''.join(["You have to specify the kind of ",
                                      "plot, use one of",
                                      str(self.LINE_IDS.keys())]))

        self.graph_drawn = False
        self.fluid_ref = fluid_ref
        self.graph_type = graph_type.upper()
        
        self.unit_system = unit_system
#         if unit_system == 'KSI':
#             self.unit_system = KSIunits()
#         elif unit_system == 'EUR':
#             self.unit_system = EURunits()
#         else:
#             self.unit_system = SIunits()
            

        self.axis = kwargs.get('axis', None)
        if self.axis is None:
            self.axis = matplotlib.pyplot.gca()
            

    def __sat_bounds(self, kind, smin=None, smax=None):
        """
        Generates limits for the saturation line in either T or p determined
        by 'kind'. If xmin or xmax are provided, values will be checked
        against the allowable range for the EOS and an error might be
        generated.

        Returns a tuple containing (xmin, xmax)
        """
        if kind == 'P':
            name = 'pressure'
            min_key = 'ptriple'
        elif kind == 'T':
            name = 'temperature'
            min_key = 'Tmin'

        fluid_min = CP.PropsSI(self.fluid_ref, min_key)
        fluid_crit = CP.PropsSI(self.fluid_ref, ''.join([kind, 'crit']))

        if smin is None:
            smin = fluid_min + SMALL
        elif smin > fluid_crit:
            raise ValueError(''.join(['Minimum ', name,
                             ' cannot be greater than fluid critical ',
                             name, '.']))

        if smax is None:
            smax = fluid_crit - SMALL
        elif smax > fluid_crit:
            raise ValueError(''.join(['Maximum ', name,
                             ' cannot be greater than fluid critical ',
                             name, '.']))

        smin = max(smin, fluid_min + SMALL)
        smax = min(smax, fluid_crit - SMALL)

        return (smin, smax)

    def _get_fluid_data(self, req_prop,
                        prop1_name, prop1_vals,
                        prop2_name, prop2_vals):
        """
        Calculates lines for constant iName (iVal) over an interval of xName
        (xVal). Returns (x[],y[]) - a tuple of arrays containing the values
        in x and y dimensions.
        """
        if len(prop1_vals) != len(prop2_vals):
            raise ValueError(''.join(['We need the same number of x value ',
                                      'arrays as iso quantities.']))

        y_vals = []
        x_vals = []

        # Calculate the values in SI units
        for i, p1_val in enumerate(prop1_vals):
            x_vals.append(prop2_vals[i])
            y_vals.append(CP.PropsSI(req_prop,
                                     prop1_name, [p1_val]*len(prop2_vals[i]), # Convert to an iterable the same size as second input array
                                     prop2_name, prop2_vals[i],
                                     self.fluid_ref))
        return numpy.array([x_vals, y_vals])

    def _get_sat_lines(self, kind='T', smin=None,
                       smax=None, num=500, x=[0., 1.]):
        """
        Calculates bubble and dew line in the quantities for your plot.
        You can specify if you need evenly spaced entries in either
        pressure or temperature by supplying kind='p' and kind='T'
        (default), respectively.
        Limits can be set with kmin (default: minimum from EOS) and
        kmax (default: critical value).
        Returns lines[] - a 2D array of dicts containing 'x' and 'y'
        coordinates for bubble and dew line. Additionally, the dict holds
        the keys 'kmax', 'label' and 'opts', those can be used for plotting
        as well.
        """
        if not kind.upper() in ['T', 'P']:
            raise ValueError(''.join(["Invalid input for determining the ",
                                      "saturation lines... Expected either ",
                                      "'T' or 'P'"]))

        smin, smax = self.__sat_bounds(kind, smin=smin, smax=smax)
        sat_range = numpy.linspace(smin, smax, num)
        sat_mesh = numpy.array([sat_range for i in x])

        x_vals = sat_mesh
        y_vals = sat_mesh
        if self.graph_type[1] != kind:
            _, x_vals = self._get_fluid_data(self.graph_type[1],
                                             'Q', x,
                                             kind, sat_mesh)

        if self.graph_type[0] != kind:
            _, y_vals = self._get_fluid_data(self.graph_type[0],
                                             'Q', x,
                                             kind, sat_mesh)
                                             
        if self.unit_system == 'KSI':
            x_vals *= self.KSI_SCALE_FACTOR[self.graph_type[1]]
            y_vals *= self.KSI_SCALE_FACTOR[self.graph_type[0]]

        # Merge the two lines, capital Y holds important information.
        # We merge on X values
        # Every entry, eg. Xy, contains two arrays of values.
        sat_lines = []
        for i in range(len(x_vals)): # two dimensions: i = {0,1}
            line = {'x': x_vals[i],
                    'y': y_vals[i],
                    'smax': smax}

            line['label'] = self.SYMBOL_MAP_KSI['Q'][0] + str(x[i])
            line['type'] = 'Q'
            line['value'] = x[i]
            line['unit'] = self.SYMBOL_MAP_KSI['Q'][1]
            line['opts'] = {'color': self.COLOR_MAP['Q'],
                            'lw': 1.0}

            if x[i] == 0.:
                line['label'] = 'bubble line'
            elif x[i] == 1.:
                line['label'] = 'dew line'
            else:
                line['opts']['lw'] = 0.75
                line['opts']['alpha'] = 0.5

            sat_lines.append(line)

        return sat_lines

    def _plot_default_annotations(self):
        def filter_fluid_ref(fluid_ref):
            fluid_ref_string = fluid_ref
            if fluid_ref.startswith('REFPROP-MIX'):
                end = 0
                fluid_ref_string = ''
                while fluid_ref.find('[', end + 1) != -1:
                    start = fluid_ref.find('&', end + 1)
                    if end == 0:
                        start = fluid_ref.find(':', end + 1)
                    end = fluid_ref.find('[', end + 1)
                    fluid_ref_string = ' '.join([fluid_ref_string,
                                                fluid_ref[start+1:end], '+'])
                fluid_ref_string = fluid_ref_string[0:len(fluid_ref_string)-2]
            return fluid_ref_string

        if len(self.graph_type) == 2:
            y_axis_id = self.graph_type[0]
            x_axis_id = self.graph_type[1]
        else:
            y_axis_id = self.graph_type[0]
            x_axis_id = self.graph_type[1:len(self.graph_type)]

        tl_str = "%s - %s Graph for %s"
        if not self.axis.get_title():
            self.axis.set_title(tl_str % (self.AXIS_LABELS[self.unit_system][y_axis_id][0],
                                          self.AXIS_LABELS[self.unit_system][x_axis_id][0],
                                          filter_fluid_ref(self.fluid_ref)))
        if not self.axis.get_xlabel():
            self.axis.set_xlabel(' '.join(self.AXIS_LABELS[self.unit_system][x_axis_id]))
        if not self.axis.get_ylabel():
            self.axis.set_ylabel(' '.join(self.AXIS_LABELS[self.unit_system][y_axis_id]))

    def _draw_graph(self):
        return

    def title(self, title):
        self.axis.set_title(title)

    def xlabel(self, xlabel):
        self.axis.set_xlabel(xlabel)

    def ylabel(self, ylabel):
        self.axis.set_ylabel(ylabel)

    def grid(self, b=None, **kwargs):
        g_map = {'on': True, 'off': False}
        if b is not None:
            b = g_map[b.lower()]
        if len(kwargs) == 0:
            self.axis.grid(b)
        else:
            self.axis.grid(kwargs)

    def set_axis_limits(self, limits):
        self.axis.set_xlim([limits[0], limits[1]])
        self.axis.set_ylim([limits[2], limits[3]])

    def show(self):
        self._draw_graph()
        matplotlib.pyplot.show()
        
    def savefig(self, *args, **kwargs):
        self._draw_graph()
        matplotlib.pyplot.savefig(*args, **kwargs)
