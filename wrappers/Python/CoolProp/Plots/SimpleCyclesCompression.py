# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

import numpy as np

import CoolProp
from .Common import process_fluid_state
from .SimpleCycles import BaseCycle, StateContainer


class BaseCompressionCycle(BaseCycle):
    """A thermodynamic cycle for vapour compression processes.

    Defines the basic properties and methods to unify access to
    compression cycle-related quantities.
    """

    def __init__(self, fluid_ref='HEOS::Water', graph_type='PH', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BaseCycle` for details."""
        BaseCycle.__init__(self, fluid_ref, graph_type, **kwargs)

    def eta_carnot_heating(self):
        """Carnot efficiency

        Calculates the Carnot efficiency for a heating process, :math:`\eta_c = \frac{T_h}{T_h-T_c}`.

        Returns
        -------
        float
        """
        Tvector = self._cycle_states.T
        return np.max(Tvector) / (np.max(Tvector) - np.min(Tvector))

    def eta_carnot_cooling(self):
        """Carnot efficiency

        Calculates the Carnot efficiency for a cooling process, :math:`\eta_c = \frac{T_c}{T_h-T_c}`.

        Returns
        -------
        float
        """
        Tvector = self._cycle_states.T
        return np.min(Tvector) / (np.max(Tvector) - np.min(Tvector))


class SimpleCompressionCycle(BaseCompressionCycle):
    """A simple vapour compression cycle"""
    STATECOUNT = 4
    STATECHANGE = [
      lambda inp: BaseCycle.state_change(inp, 'S', 'P', 0, ty1='log', ty2='log'),  # Compression process
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 1, ty1='lin', ty2='lin'),  # Heat removal
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 2, ty1='log', ty2='log'),  # Expansion
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 3, ty1='lin', ty2='lin')  # Heat addition
      ]

    def __init__(self, fluid_ref='HEOS::Water', graph_type='PH', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCyclesCompression.BaseCompressionCycle` for details."""
        BaseCompressionCycle.__init__(self, fluid_ref, graph_type, **kwargs)

    def simple_solve(self, T0, p0, T2, p2, eta_com, fluid=None, SI=True):
        """"
        A simple vapour compression cycle calculation

        Parameters
        ----------
        T0 : float
            The evaporated fluid, before the compressor
        p0 : float
            The evaporated fluid, before the compressor
        T2 : float
            The condensed fluid, before the expansion valve
        p2 : float
            The condensed fluid, before the expansion valve
        eta_com : float
            Isentropic compressor efficiency

        Examples
        --------
        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> from CoolProp.Plots import SimpleCompressionCycle
        >>> pp = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR')
        >>> pp.calc_isolines(CoolProp.iQ, num=11)
        >>> cycle = SimpleCompressionCycle('HEOS::R134a', 'PH', unit_system='EUR')
        >>> T0 = 280
        >>> pp.state.update(CoolProp.QT_INPUTS,0.0,T0-15)
        >>> p0 = pp.state.keyed_output(CoolProp.iP)
        >>> T2 = 310
        >>> pp.state.update(CoolProp.QT_INPUTS,1.0,T2+10)
        >>> p2 = pp.state.keyed_output(CoolProp.iP)
        >>> cycle.simple_solve(T0, p0, T2, p2, 0.7, SI=True)
        >>> cycle.steps = 50
        >>> sc = cycle.get_state_changes()
        >>> import matplotlib.pyplot as plt
        >>> plt.close(cycle.figure)
        >>> pp.draw_process(sc)

        """
        if fluid is not None:
            self.state = process_fluid_state(fluid)
        if self._state is None:
            raise ValueError("You have to specify a fluid before you can calculate.")

        cycle_states = StateContainer(unit_system=self._system)

        if not SI:
            conv_t = self._system[CoolProp.iT].to_SI
            conv_p = self._system[CoolProp.iP].to_SI
            T0 = conv_t(T0)
            p0 = conv_p(p0)
            T2 = conv_t(T2)
            p2 = conv_p(p2)

        # Gas from evaporator
        self.state.update(CoolProp.PT_INPUTS, p0, T0)
        h0 = self.state.hmass()
        s0 = self.state.smass()
        # Just a showcase for the different accessor methods
        cycle_states[0, 'H'] = h0
        cycle_states[0]['S'] = s0
        cycle_states[0][CoolProp.iP] = p0
        cycle_states[0, CoolProp.iT] = T0

        # Pressurised vapour
        p1 = p2
        self.state.update(CoolProp.PSmass_INPUTS, p1, s0)
        h1 = h0 + (self.state.hmass() - h0) / eta_com
        self.state.update(CoolProp.HmassP_INPUTS, h1, p1)
        s1 = self.state.smass()
        T1 = self.state.T()
        cycle_states[1, 'H'] = h1
        cycle_states[1, 'S'] = s1
        cycle_states[1, 'P'] = p1
        cycle_states[1, 'T'] = T1

        # Condensed vapour
        self.state.update(CoolProp.PT_INPUTS, p2, T2)
        h2 = self.state.hmass()
        s2 = self.state.smass()
        cycle_states[2, 'H'] = h2
        cycle_states[2, 'S'] = s2
        cycle_states[2, 'P'] = p2
        cycle_states[2, 'T'] = T2

        # Expanded fluid, 2-phase
        p3 = p0
        h3 = h2
        self.state.update(CoolProp.HmassP_INPUTS, h3, p3)
        s3 = self.state.smass()
        T3 = self.state.T()
        cycle_states[3, 'H'] = h3
        cycle_states[3, 'S'] = s3
        cycle_states[3, 'P'] = p3
        cycle_states[3, 'T'] = T3

        w_net = h0 - h1
        q_evap = h0 - h3

        self.cycle_states = cycle_states
        self.fill_states()

    def simple_solve_dt(self, Te, Tc, dT_sh, dT_sc, eta_com, fluid=None, SI=True):
        """"
        A simple vapour compression cycle calculation based on
        superheat, subcooling and temperatures.

        Parameters
        ----------
        Te : float
            The evaporation temperature
        Tc : float
            The condensation temperature
        dT_sh : float
            The superheat after the evaporator
        dT_sc : float
            The subcooling after the condenser
        eta_com : float
            Isentropic compressor efficiency

        Examples
        --------
        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> from CoolProp.Plots import SimpleCompressionCycle
        >>> pp = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR')
        >>> pp.calc_isolines(CoolProp.iQ, num=11)
        >>> cycle = SimpleCompressionCycle('HEOS::R134a', 'PH', unit_system='EUR')
        >>> Te = 265
        >>> Tc = 300
        >>> cycle.simple_solve_dt(Te, Tc, 10, 15, 0.7, SI=True)
        >>> cycle.steps = 50
        >>> sc = cycle.get_state_changes()
        >>> import matplotlib.pyplot as plt
        >>> plt.close(cycle.figure)
        >>> pp.draw_process(sc)
        """
        if fluid is not None:
            self.state = process_fluid_state(fluid)
        if self._state is None:
            raise ValueError("You have to specify a fluid before you can calculate.")

        if not SI:
            conv_t = self._system[CoolProp.iT].to_SI
            Te = conv_p(Te)
            Tc = conv_p(Tc)

        # Get the saturation conditions
        self.state.update(CoolProp.QT_INPUTS, 1.0, Te)
        p0 = self.state.p()
        self.state.update(CoolProp.QT_INPUTS, 0.0, Tc)
        p2 = self.state.p()

        T0 = Te + dT_sh
        T2 = Tc - dT_sc

        self.simple_solve(T0, p0, T2, p2, eta_com, fluid=None, SI=True)

    def COP_heating(self):
        """COP for a heating process

        Calculates the coefficient of performance for a heating process, :math:`COP_h = \frac{q_{con}}{w_{comp}}`.

        Returns
        -------
        float
        """
        return (self.cycle_states[1, 'H'] - self.cycle_states[2, 'H']) / (self.cycle_states[1, 'H'] - self.cycle_states[0, 'H'])

    def COP_cooling(self):
        """COP for a cooling process

        Calculates the coefficient of performance for a cooling process, :math:`COP_c = \frac{q_{eva}}{w_{comp}}`.

        Returns
        -------
        float
        """
        return (self.cycle_states[0, 'H'] - self.cycle_states[3, 'H']) / (self.cycle_states[1, 'H'] - self.cycle_states[0, 'H'])
