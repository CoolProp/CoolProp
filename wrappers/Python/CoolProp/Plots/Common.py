# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals

import matplotlib
import numpy as np

import CoolProp.CoolProp as CP
from abc import ABCMeta
from CoolProp import AbstractState
import CoolProp
import warnings


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
        self.T.unit=ur'\u00B0 C'


class IsoLine(object):
    """An object that holds the functions to calculate a line of 
    a constant property in the dimensions of a property plot. This 
    class only uses SI units."""
    
    # A list of supported plot
    PLOTS = ['TS','PH','HS','PS','PD','TD','PT','PU']
    
    # Normally we calculate a sweep in x-dimensions, but 
    # sometimes a sweep in y-dimensions is better.
    XY_SWITCH = {
      'D': { 'TS':True , 'PH':True , 'HS':False, 'PS':True , 'PD':False, 'TD':False, 'PT':False, 'PU':False},
      'H': { 'TS':False, 'PH':False, 'HS':False, 'PS':False, 'PD':False, 'TD':False, 'PT':False, 'PU':False},
      'P': { 'TS':False, 'PH':False, 'HS':False, 'PS':False, 'PD':False, 'TD':False, 'PT':False, 'PU':False},
      'S': { 'TS':False, 'PH':False, 'HS':False, 'PS':False, 'PD':False, 'TD':False, 'PT':False, 'PU':False},
      'T': { 'TS':False, 'PH':False, 'HS':False, 'PS':False, 'PD':False, 'TD':False, 'PT':False, 'PU':False},
      'U': { 'TS':False, 'PH':False, 'HS':False, 'PS':False, 'PD':False, 'TD':False, 'PT':False, 'PU':False},  
    }
    
    @classmethod
    def get_update_pair(cls,i,x,y):
        """Processes the values for the isoproperty and the graph dimensions
        to figure which should be used as inputs to the state update. Returns
        a tuple with the indices for the update call and the property constant.
        For an isobar in a Ts-diagram it returns the default order and the 
        correct constant for the update pair:
        get_update_pair('P','S','T') -> (0,1,2,CoolProp.PSmass_INPUTS)
        other values require switching and swapping
        get_update_pair('S','P','H') -> (1,0,2,CoolProp.PSmass_INPUTS)
        """
        
        #ii = CP.get_parameter_index(i)
        
        # Figure out if x or y-dimension should be used
        switch = cls.XY_SWITCH[i]
        if switch is None:
            raise ValueError("This isoline cannot be calculated!")
        elif switch is False:
            oo = CP.get_parameter_index(i)
            tt = CP.get_parameter_index(x)
            second = None
            third  = 2
        elif switch is True:
            oo = CP.get_parameter_index(i)
            tt = CP.get_parameter_index(y)
            second = 2
            third  = None
        else:
            raise ValueError("Unknown error!")
        
        pair, out1, _ = CP.generate_update_pair(oo,0.0,tt,1.0)
        if out1==0.0: # Correct order
            first = 0
            if second is None:
                second = 1
            else:
                third  = 1 
        else: # Wrong order
            first = 1
            if second is None:
                second = 0
            else:
                third  = 0
        return first,second,third,pair 

    @property
    def x(self): return self._x
    @x.setter
    def x(self, value): self._x = np.array(value)
    @property
    def y(self): return self._y
    @y.setter
    def y(self, value): self._y = np.array(value)
    @property
    def line_type(self): return self._line_type
    @line_type.setter
    def line_type(self, value): self._line_type = str(value).upper()
    @property
    def graph_type(self): return self._graph_type
    @graph_type.setter
    def graph_type(self, value): self._graph_type = str(value).upper()
    @property
    def value(self): return self._value
    @value.setter
    def value(self, value): self._value = float(value)


class IsoLineCalculator(object):
    
    
    def __init__(self, state):
        self.DEBUG = False

        # direct geometry
        self.X     = None #
        self.Y     = None #
        self.type  = None #
        self.value = None #
        self.unit  = None #
        self.opts  = None #


class BasePlot(object):
    """The base class for all plots. It can be instantiated itself, but provides many 
    general facilities to be used in the different plots. """
    __metaclass__ = ABCMeta
    
    # Define the iteration keys
    PROPERTIES = {
      'D': 'density',
      'H': 'specific enthalpy',
      'P': 'pressure',
      'S': 'specific entropy',
      'T': 'temperature',
      'U': 'specific internal energy'
    }
    
    # A list of supported plot
    PLOTS = ['TS','PH','HS','PS','PD','TD','PT','PU']
    
    # Define the unit systems
    UNIT_SYSTEMS = {
      'SI' : SIunits(),
      'KSI': KSIunits(),
      'EUR': EURunits()
    }
    
    LINE_COLORS = {
      'T': 'Darkred',
      'P': 'DarkCyan',
      'H': 'DarkGreen',
      'D': 'DarkBlue',
      'S': 'DarkOrange',
      'Q': 'black'
    }

    def __init__(self, fluid_ref, graph_type, unit_system = 'KSI', **kwargs):
        
        # Process the fluid and set self._state
        if isinstance(fluid_ref, str):
            # TODO: Fix the backend extraction etc
            fluid_def = fluid_ref.split('::')
            if len(fluid_def)==2:
                backend = fluid_def[0]
                fluid = fluid_def[1]
            elif len(fluid_def)==1:
                backend = "HEOS"
                fluid = fluid_def[0]
            else: 
                raise ValueError("This is not a valid fluid_ref string: {0:s}".format(str(fluid_ref)))
            self._state = AbstractState(backend, fluid)
        elif isinstance(fluid_ref, AbstractState):
            self._state = fluid_ref
        else:
            raise TypeError("Invalid fluid_ref input, expected a string or an abstract state instance")
        
        # Process the graph_type and set self._x_type and self._y_type
        graph_type = graph_type.upper()
        graph_type = graph_type.replace(r'RHO',r'D')
        if graph_type in self.PLOTS:
            self._y_type = graph_type[0]
            self._x_type = graph_type[1]
        else:
            raise ValueError("Invalid graph_type input, expected a string from {0:s}".format(str(self.PLOTS)))
        
        # Process the unit_system and set self._system
        unit_system = unit_system.upper()
        if unit_system in self.UNIT_SYSTEMS:
            self._system = self.UNIT_SYSTEMS[unit_system]
        else:
            raise ValueError("Invalid unit_system input, expected a string from {0:s}".format(str(self.UNIT_SYSTEMS.keys())))
        
        self._axis = kwargs.get('axis', matplotlib.pyplot.gca())        
        self.small = kwargs.get('small', 1e-5)

        self._colors = self.LINE_COLORS.copy()
        colors = kwargs.get('colors', None)
        if colors is not None:
            self._colors.update(colors)
        
        self._graph_drawn = False

    @property
    def small(self): return self._small
    @small.setter
    def small(self, value):
        self._T_small = self._state.trivial_keyed_output(CoolProp.iT_critical)*value
        self._P_small = self._state.trivial_keyed_output(CoolProp.iP_critical)*value
        self._small   = value

    def __sat_bounds(self, kind, smin=None, smax=None):
        """Generates limits for the saturation line in either T or p determined
        by 'kind'. If xmin or xmax are provided, values will be checked
        against the allowable range for the EOS and a warning might be
        generated. Returns a tuple containing (xmin, xmax)"""

        # TODO: REFPROP backend does not have ptriple.
        T_triple = self._state.trivial_keyed_output(CoolProp.iT_triple)
        T_min    = self._state.trivial_keyed_output(CoolProp.iT_min)        
        self._state.update(CoolProp.QT_INPUTS, 0, max([T_triple,T_min])+self._T_small)
        kind = kind.upper()
        if kind == 'P':
            fluid_min = self._state.keyed_output(CoolProp.iP)
            fluid_max = self._state.trivial_keyed_output(CoolProp.iP_critical)-self._P_small
        elif kind == 'T':
            fluid_min = self._state.keyed_output(CoolProp.iT)
            fluid_max = self._state.trivial_keyed_output(CoolProp.iT_critical)-self._T_small
        else:
            raise ValueError("Saturation boundaries have to be defined in T or P, but not in {0:s}".format(str(kind)))
                
        if smin is not None: 
            if fluid_min < smin < fluid_max:
                sat_min = smin
            else:
                warnings.warn(
                  "Your minimum {0:s} has been ignored, {1:f} is not between {2:f} and {3:f}".format(self.PROPERTIES[kind],smin,fluid_min,fluid_max),
                  UserWarning)
                sat_min = fluid_min
        else:
            sat_min = fluid_min
            
        if smax is not None: 
            if fluid_min < smax < fluid_max:
                sat_max = smax
            else:
                warnings.warn(
                  "Your maximum {0:s} has been ignored, {1:f} is not between {2:f} and {3:f}".format(self.PROPERTIES[kind],smax,fluid_min,fluid_max),
                  UserWarning)
                sat_max = fluid_max
        else:
            sat_max = fluid_max

        return (sat_min, sat_max)
    
    
    def _get_sat_lines(self, kind='T', smin=None,
                       smax=None, num=500, x=[0., 1.]):
        """
        Calculates bubble and dew line in the quantities for your plot.
        You can specify if you need evenly spaced entries in either
        pressure or temperature by supplying kind='p' and kind='T'
        (default), respectively.
        Limits can be set with kmin (default: triple point or EOS minimum) and
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
