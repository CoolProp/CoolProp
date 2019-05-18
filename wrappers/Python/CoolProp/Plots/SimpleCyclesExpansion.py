# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

import numpy as np

import CoolProp
from .Common import process_fluid_state
from .SimpleCycles import BaseCycle, StateContainer


class BasePowerCycle(BaseCycle):
    """A thermodynamic cycle for power producing processes.

    Defines the basic properties and methods to unify access to
    power cycle-related quantities.
    """

    def __init__(self, fluid_ref='HEOS::Water', graph_type='TS', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BaseCycle` for details."""
        BaseCycle.__init__(self, fluid_ref, graph_type, **kwargs)

    def eta_carnot(self):
        """Carnot efficiency

        Calculates the Carnot efficiency for the specified process, :math:`\eta_c = 1 - \frac{T_c}{T_h}`.

        Returns
        -------
        float
        """
        Tvector = self._cycle_states.T
        return 1. - np.min(Tvector) / np.max(Tvector)

    def eta_thermal(self):
        """Thermal efficiency

        The thermal efficiency for the specified process(es), :math:`\eta_{th} = \frac{\dot{W}_{exp} - \dot{W}_{pum}}{\dot{Q}_{in}}`.

        Returns
        -------
        float
        """
        raise NotImplementedError("Implement it in the subclass.")


class SimpleRankineCycle(BasePowerCycle):
    """A simple Rankine cycle *without* regeneration"""
    STATECOUNT = 4
    STATECHANGE = [
      lambda inp: BaseCycle.state_change(inp, 'S', 'P', 0, ty1='log', ty2='log'),  # Pumping process
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 1, ty1='lin', ty2='lin'),  # Heat addition
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 2, ty1='log', ty2='log'),  # Expansion
      lambda inp: BaseCycle.state_change(inp, 'H', 'P', 3, ty1='lin', ty2='lin')  # Heat removal
      ]

    def __init__(self, fluid_ref='HEOS::Water', graph_type='TS', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BasePowerCycle` for details."""
        BasePowerCycle.__init__(self, fluid_ref, graph_type, **kwargs)

    def simple_solve(self, T0, p0, T2, p2, eta_exp, eta_pum, fluid=None, SI=True):
        """"
        A simple Rankine cycle calculation

        Parameters
        ----------
        T0 : float
            The coldest point, before the pump
        p0 : float
            The lowest pressure, before the pump
        T2 : float
            The hottest point, before the expander
        p2 : float
            The highest pressure, before the expander
        eta_exp : float
            Isentropic expander efficiency
        eta_pum : float
            Isentropic pump efficiency

        Examples
        --------
        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> from CoolProp.Plots import SimpleRankineCycle
        >>> pp = PropertyPlot('HEOS::Water', 'TS', unit_system='EUR')
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
        >>> import matplotlib.pyplot as plt
        >>> plt.close(cycle.figure)
        >>> pp.draw_process(sc)

        """
        if fluid is not None: self.state = process_fluid_state(fluid)
        if self._state is None:
            raise ValueError("You have to specify a fluid before you can calculate.")

        cycle_states = StateContainer(unit_system=self._system)

        if not SI:
            Tc = self._system[CoolProp.iT].to_SI
            pc = self._system[CoolProp.iP].to_SI
            T0 = Tc(T0)
            p0 = pc(p0)
            T2 = Tc(T2)
            p2 = pc(p2)

        # Subcooled liquid
        self.state.update(CoolProp.PT_INPUTS, p0, T0)
        h0 = self.state.hmass()
        s0 = self.state.smass()
        # Just a showcase for the different accessor methods
        cycle_states[0, 'H'] = h0
        cycle_states[0]['S'] = s0
        cycle_states[0][CoolProp.iP] = p0
        cycle_states[0, CoolProp.iT] = T0

        # Pressurised liquid
        p1 = p2
        self.state.update(CoolProp.PSmass_INPUTS, p1, s0)
        h1 = h0 + (self.state.hmass() - h0) / eta_pum
        self.state.update(CoolProp.HmassP_INPUTS, h1, p1)
        s1 = self.state.smass()
        T1 = self.state.T()
        cycle_states[1, 'H'] = h1
        cycle_states[1, 'S'] = s1
        cycle_states[1, 'P'] = p1
        cycle_states[1, 'T'] = T1

        # Evaporated vapour
        self.state.update(CoolProp.PT_INPUTS, p2, T2)
        h2 = self.state.hmass()
        s2 = self.state.smass()
        cycle_states[2, 'H'] = h2
        cycle_states[2, 'S'] = s2
        cycle_states[2, 'P'] = p2
        cycle_states[2, 'T'] = T2

        # Expanded gas
        p3 = p0
        self.state.update(CoolProp.PSmass_INPUTS, p3, s2)
        h3 = h2 - eta_exp * (h2 - self.state.hmass())
        self.state.update(CoolProp.HmassP_INPUTS, h3, p3)
        s3 = self.state.smass()
        T3 = self.state.T()
        cycle_states[3, 'H'] = h3
        cycle_states[3, 'S'] = s3
        cycle_states[3, 'P'] = p3
        cycle_states[3, 'T'] = T3

        w_net = h2 - h3
        q_boiler = h2 - h1
        eta_th = w_net / q_boiler

        self.cycle_states = cycle_states
        self.fill_states()

    def eta_thermal(self):
        """Thermal efficiency

        The thermal efficiency for the specified process(es), :math:`\eta_{th} = \frac{\dot{W}_{exp} - \dot{W}_{pum}}{\dot{Q}_{in}}`.

        Returns
        -------
        float
        """
        w_net = self.cycle_states[2].H - self.cycle_states[3].H - (self.cycle_states[1].H - self.cycle_states[0].H)
        q_boiler = self.cycle_states[2].H - self.cycle_states[1].H
        return w_net / q_boiler
