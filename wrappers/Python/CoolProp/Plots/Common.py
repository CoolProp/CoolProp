# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np

import CoolProp.CoolProp as CP
from abc import ABCMeta
from CoolProp import AbstractState
from CoolProp.CoolProp import PropsSI
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
    @property
    def Q(self): return self._Q
    @Q.setter
    def Q(self, value): self._Q = value
    
    @property
    def dimensions(self): 
        return {
      CoolProp.iDmass : self._D,
      CoolProp.iHmass : self._H,
      CoolProp.iP     : self._P,
      CoolProp.iSmass : self._S,
      CoolProp.iT     : self._T,
      CoolProp.iUmass : self._U,
      CoolProp.iQ     : self._Q
    }


class SIunits(UnitSystem):
    def __init__(self):
        self._D = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Density',                  symbol=r'\rho', unit=r'kg/m$^3$')
        self._H = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Enthalpy',        symbol=r'h',    unit=r'J/kg')
        self._P = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Pressure',                 symbol=r'p',    unit=r'Pa')
        self._S = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Entropy',         symbol=r's',    unit=r'J/kg/K')
        self._T = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Temperature',              symbol=r'T',    unit=r'K')
        self._U = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Internal Energy', symbol=r'u',    unit=r'J/kg')
        self._Q = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Vapour Quality',           symbol=r'x',    unit=r'')

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


class Base2DObject(object):
    """A container for shared settings and constants for the 
    isolines and the property plots."""
    
    __metaclass__ = ABCMeta
    
    # A list of supported plot
    TS = CoolProp.iT*10     + CoolProp.iSmass
    PH = CoolProp.iP*10     + CoolProp.iHmass
    HS = CoolProp.iHmass*10 + CoolProp.iSmass
    PS = CoolProp.iP*10     + CoolProp.iSmass
    PD = CoolProp.iP*10     + CoolProp.iDmass
    TD = CoolProp.iT*10     + CoolProp.iDmass
    PT = CoolProp.iP*10     + CoolProp.iT
    PU = CoolProp.iP*10     + CoolProp.iUmass
    
    PLOTS = {
      'TS': TS,
      'PH': PH,
      'HS': HS,
      'PS': PS,
      'PD': PD,
      'TD': TD,
      'PT': PT,
    }
    
    PLOTS_INV = {v: k for k, v in PLOTS.items()}
    
#     # A list of supported plot
#     @property
#     def TS(self): return type(self).TS
#     @property
#     def PH(self): return CoolProp.iP*10     + CoolProp.iHmass
#     @property
#     def HS(self): return CoolProp.iHmass*10 + CoolProp.iSmass
#     @property
#     def PS(self): return CoolProp.iP*10     + CoolProp.iSmass
#     @property
#     def PD(self): return CoolProp.iP*10     + CoolProp.iDmass
#     @property
#     def TD(self): return CoolProp.iT*10     + CoolProp.iDmass
#     @property
#     def PT(self): return CoolProp.iP*10     + CoolProp.iT
#     @property
#     def PU(self): return CoolProp.iP*10     + CoolProp.iUmass

    def __init__(self, x_type, y_type, state=None):
        self._x_index = self._get_index(x_type)
        self._y_index = self._get_index(y_type)        
        if state is not None: self.state = state
        else: self._state = None 

    # A list of supported plot
    @property
    def x_index(self): return self._x_index
    @property
    def y_index(self): return self._y_index
    @property
    def state(self): return self._state
    @state.setter
    def state(self, value):
        if isinstance(value, AbstractState): self._state = value
        else: raise TypeError("Invalid state input, expected an AbstractState instance.")
    
    def _get_index(self,prop):
        if isinstance(prop, basestring):
            return CP.get_parameter_index(prop)
        elif isinstance(prop, int):
            return prop
        else:
            raise ValueError("Invalid input, expected a string or an int, not {0:s}.".format(str(prop)))
    

class IsoLine(Base2DObject):
    """An object that holds the functions to calculate a line of 
    a constant property in the dimensions of a property plot. This 
    class only uses SI units."""
    
    # Normally we calculate a sweep in x-dimensions, but 
    # sometimes a sweep in y-dimensions is better.
    XY_SWITCH = {
      CoolProp.iDmass: { Base2DObject.TS:True , Base2DObject.PH:True , Base2DObject.HS:False, Base2DObject.PS:True , Base2DObject.PD:None , Base2DObject.TD:None , Base2DObject.PT:False},
      CoolProp.iHmass: { Base2DObject.TS:False, Base2DObject.PH:None , Base2DObject.HS:None , Base2DObject.PS:True , Base2DObject.PD:True , Base2DObject.TD:False, Base2DObject.PT:False},
      CoolProp.iP    : { Base2DObject.TS:False, Base2DObject.PH:None , Base2DObject.HS:False, Base2DObject.PS:None , Base2DObject.PD:None , Base2DObject.TD:False, Base2DObject.PT:None },
      CoolProp.iSmass: { Base2DObject.TS:None , Base2DObject.PH:True , Base2DObject.HS:None , Base2DObject.PS:None , Base2DObject.PD:True , Base2DObject.TD:False, Base2DObject.PT:True },
      CoolProp.iT    : { Base2DObject.TS:None , Base2DObject.PH:True , Base2DObject.HS:False, Base2DObject.PS:False, Base2DObject.PD:False, Base2DObject.TD:None , Base2DObject.PT:None }
    }
    
    def __init__(self, i_index, x_index, y_index, value=0.0, state=None):
        super(IsoLine, self).__init__(x_index, y_index, state)
        self._i_index = self._get_index(i_index)
        if value is not None: self.value = value
        else: self._value = None 
        self._x       = None
        self._y       = None
        
    @property
    def i_index(self): return self._i_index
    @property
    def value(self): return self._value
    @value.setter
    def value(self, value): self._value = float(value)
    @property
    def x(self): return self._x
    @x.setter
    def x(self, value): self._x = np.array(value)
    @property
    def y(self): return self._y
    @y.setter
    def y(self, value): self._y = np.array(value)
    
    def _get_update_pair(self):
        """Processes the values for the isoproperty and the graph dimensions
        to figure which should be used as inputs to the state update. Returns
        a tuple with the indices for the update call and the property constant.
        For an isobar in a Ts-diagram it returns the default order and the 
        correct constant for the update pair:
        get_update_pair(CoolProp.iP,CoolProp.iSmass,CoolProp.iT) -> (0,1,2,CoolProp.PSmass_INPUTS)
        other values require switching and swapping.
        """
        # Figure out if x or y-dimension should be used
        switch = self.XY_SWITCH[self.i_index][self.y_index*10+self.x_index]
        if switch is None:
            raise ValueError("This isoline cannot be calculated!")
        elif switch is False:
            pair, out1, _ = CP.generate_update_pair(self.i_index,0.0,self.x_index,1.0)
        elif switch is True:
            pair, out1, _ = CP.generate_update_pair(self.i_index,0.0,self.y_index,1.0)
        else:
            raise ValueError("Unknown error!")

        if out1==0.0: # Correct order
            swap = False
        else: # Wrong order
            swap = True
        
        if not switch and not swap:
            return 0,1,2,pair
        elif switch and not swap:
            return 0,2,1,pair
        elif not switch and swap:
            return 1,0,2,pair
        elif switch and swap:
            return 1,2,0,pair
        else:
            raise ValueError("Check the code, this should not happen!")
    
    def calc_range(self,xvals=None,yvals=None):
        ipos,xpos,ypos,pair = self._get_update_pair()
        
        order = [ipos,xpos,ypos]
        idxs  = [v for (_,v) in sorted(zip(order,[self.i_index        , self.x_index, self.y_index]))]
        vals  = [v for (_,v) in sorted(zip(order,[np.array(self.value), xvals       , yvals       ]))]
        
        if vals[0] is None or vals[1] is None:
            raise ValueError("One required input is missing, make sure to supply the correct xvals or yvals: {0:s} - {1:s}".format(str(xvals),str(yvals)))
         
        if vals[0].size > vals[1].size:
            vals[1] = np.resize(vals[1],vals[0].shape)
        elif vals[0].size < vals[1].size:
            vals[0] = np.resize(vals[0],vals[1].shape)
            
        it = np.nditer([vals[0], vals[1], None])
        for x, y, z in it:
            self.state.update(pair, x, y)
            z[...] = self.state.keyed_output(idxs[2])
            
        for i,v in enumerate(idxs):
            if v == self.x_index: self.x = it.operands[i]
            if v == self.y_index: self.y = it.operands[i]
        
            

        
        
        
        
                        
        
        
        



class BasePlot(Base2DObject):
    """The base class for all plots. It can be instantiated itself, but provides many 
    general facilities to be used in the different plots. """
    #__metaclass__ = ABCMeta
    
    # Define the iteration keys
    PROPERTIES = {
      'D': 'density',
      'H': 'specific enthalpy',
      'P': 'pressure',
      'S': 'specific entropy',
      'T': 'temperature',
      'U': 'specific internal energy'
    }
    
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
        if isinstance(fluid_ref, basestring):
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
            state = AbstractState(backend, fluid)
        elif isinstance(fluid_ref, AbstractState):
            state = fluid_ref
        else:
            raise TypeError("Invalid fluid_ref input, expected a string or an abstract state instance")
        
        # Process the graph_type and set self._x_type and self._y_type
        graph_type = graph_type.upper()
        graph_type = graph_type.replace(r'RHO',r'D')
        if graph_type not in self.PLOTS:
            raise ValueError("Invalid graph_type input, expected a string from {0:s}".format(str(self.PLOTS)))
        
        # call the base class
        super(BasePlot, self).__init__(graph_type[1], graph_type[0], state)
        
        # Process the unit_system and set self._system
        unit_system = unit_system.upper()
        if unit_system in self.UNIT_SYSTEMS:
            self._system = self.UNIT_SYSTEMS[unit_system]
        else:
            raise ValueError("Invalid unit_system input, expected a string from {0:s}".format(str(self.UNIT_SYSTEMS.keys())))
        
        self._axis = kwargs.get('axis', plt.gca())        
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

    def _get_sat_bounds(self, kind, smin=None, smax=None):
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
    
    
    def _get_iso_label(self, isoline, unit=True):
        if self._system is not None:
            dim = self._system.dimensions[isoline.i_index]
            return str(r"$"+dim.symbol+"="+str(dim.from_SI(isoline.value))+ "$ "+dim.unit if unit else "$").strip()
        return str(isoline.value).strip()
        
        
    
    
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
        sat_range = np.linspace(smin, smax, num)
        sat_mesh = np.array([sat_range for i in x])

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
        plt.show()
        
    def savefig(self, *args, **kwargs):
        self._draw_graph()
        plt.savefig(*args, **kwargs)



if __name__ == "__main__":
    for sys in [SIunits(), KSIunits(), EURunits()]:
        print(sys.H.label)
        print(sys.H.to_SI(20))
        print(sys.P.label)
        print(sys.P.to_SI(20))
        
        #i_index, x_index, y_index, value=None, state=None)
        iso = IsoLine('T','H','P')
        print(iso._get_update_pair())
        
        state = AbstractState("HEOS","water")
        iso = IsoLine('T','H','P', 300.0, state)
        hr = PropsSI("H","T",[290,310],"P",[1e5,1e5],"water")
        iso.calc_range(hr,np.array([0.9e5,1.1e5]))
        print(iso.x,iso.y)        
        
        
        #bp = BasePlot(fluid_ref, graph_type, unit_system = 'KSI', **kwargs):
        bp = BasePlot('n-Pentane', 'PH')
        print(bp._get_sat_bounds('P'))
        print(bp._get_iso_label(iso))
        
        
        # get_update_pair(CoolProp.iP,CoolProp.iSmass,CoolProp.iT) -> (0,1,2,CoolProp.PSmass_INPUTS)
        #other values require switching and swapping
        #get_update_pair(CoolProp.iSmass,CoolProp.iP,CoolProp.iHmass) -> (1,0,2,CoolProp.PSmass_INPUTS)