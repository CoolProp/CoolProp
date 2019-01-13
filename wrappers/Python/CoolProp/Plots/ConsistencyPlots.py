# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import


import matplotlib.pyplot as plt
import numpy as np
import time, timeit
import six
import pandas
import CoolProp as CP

CP.CoolProp.set_debug_level(00)
from matplotlib.backends.backend_pdf import PdfPages

all_solvers = ['PT', 'DmolarT', 'HmolarP', 'PSmolar', 'SmolarT', 'PUmolar', 'DmolarP', 'DmolarHmolar', 'DmolarSmolar', 'DmolarUmolar', 'HmolarSmolar', 'HmolarT', 'TUmolar', 'SmolarUmolar', 'HmolarUmolar']
not_implemented_solvers = ['SmolarUmolar', 'HmolarUmolar', 'TUmolar', 'HmolarT']
no_two_phase_solvers = ['PT']

implemented_solvers = [pair for pair in all_solvers if pair not in not_implemented_solvers]

param_labels = dict(Hmolar='Enthalpy [J/mol]/1000',
                    Smolar='Entropy [J/mol/K]/1000',
                    Umolar='Int. Ener. [J/mol]/1000',
                    T='Temperature [K]',
                    Dmolar='Density [mol/m3]/1000',
                    P='Pressure [Pa]/1000')


def split_pair(pair):
    for key in ['Dmolar', 'Hmolar', 'Smolar', 'P', 'T', 'Umolar']:
        if pair.startswith(key):
            return key, pair.replace(key, '')


def split_pair_xy(pair):
    if pair == 'HmolarP':
        return 'Hmolar', 'P'
    elif pair == 'PSmolar':
        return 'Smolar', 'P'
    elif pair == 'PUmolar':
        return 'Umolar', 'P'
    elif pair == 'PT':
        return 'T', 'P'
    elif pair == 'DmolarT':
        return 'Dmolar', 'T'
    elif pair == 'SmolarT':
        return 'Smolar', 'T'
    elif pair == 'TUmolar':
        return 'Umolar', 'T'
    elif pair == 'HmolarT':
        return 'Hmolar', 'T'
    elif pair == 'DmolarP':
        return 'Dmolar', 'P'
    elif pair == 'DmolarHmolar':
        return 'Dmolar', 'Hmolar'
    elif pair == 'DmolarSmolar':
        return 'Dmolar', 'Smolar'
    elif pair == 'DmolarUmolar':
        return 'Dmolar', 'Umolar'
    elif pair == 'HmolarSmolar':
        return 'Smolar', 'Hmolar'
    elif pair == 'SmolarUmolar':
        return 'Smolar', 'Umolar'
    elif pair == 'HmolarUmolar':
        return 'Hmolar', 'Umolar'
    else:
        raise ValueError(pair)


DEBUG_LEVEL = 1


def myprint(level, *args, **kwargs):
    if level > DEBUG_LEVEL:
        print(*args, **kwargs)


class ConsistencyFigure(object):
    def __init__(self, fluid, figsize=(15, 23), backend='HEOS', additional_skips=[], mole_fractions=None, p_limits_1phase=None, T_limits_1phase=None, NT_1phase=40, Np_1phase=40, NT_2phase=20, NQ_2phase=20):

        self.fluid = fluid
        self.backend = backend
        print('***********************************************************************************')
        print('*************** ' + backend + '::' + fluid + ' ************************')
        print('***********************************************************************************')
        self.fig, self.axes = plt.subplots(nrows=5, ncols=3, figsize=figsize)
        self.pairs = all_solvers
        pairs_generator = iter(self.pairs)

        states = [CP.AbstractState(backend, fluid) for _ in range(3)]
        if mole_fractions is not None:
            for state in states:
                state.set_mole_fractions(mole_fractions)
        self.axes_list = []
        for row in self.axes:
            for ax in row:
                pair = six.next(pairs_generator)
                kwargs = dict(p_limits_1phase=p_limits_1phase, T_limits_1phase=T_limits_1phase, NT_1phase=NT_1phase, Np_1phase=Np_1phase,
                              NT_2phase=NT_2phase, NQ_2phase=NQ_2phase)
                self.axes_list.append(ConsistencyAxis(ax, self, pair, self.fluid, self.backend, *states, **kwargs))
                ax.set_title(pair)

        self.calc_saturation_curves()
        self.plot_saturation_curves()

        self.calc_Tmax_curve()
        self.plot_Tmax_curve()

        self.calc_melting_curve()
        self.plot_melting_curve()

        self.tight_layout()

        self.fig.subplots_adjust(top=0.95)
        self.fig.suptitle('Consistency plots for ' + self.fluid, size=14)

        errors = []
        for i, (ax, pair) in enumerate(zip(self.axes_list, self.pairs)):
            if pair not in not_implemented_solvers and pair not in additional_skips:
                errors.append(ax.consistency_check_singlephase())
                if pair not in no_two_phase_solvers:
                    ax.consistency_check_twophase()
            else:
                ax.cross_out_axis()

        self.errors = pandas.concat(errors, sort=True)

    def calc_saturation_curves(self):
        """
        Calculate all the saturation curves in one shot using the state class to save computational time
        """
        HEOS = CP.AbstractState(self.backend, self.fluid)
        self.dictL, self.dictV = {}, {}
        for Q, dic in zip([0, 1], [self.dictL, self.dictV]):
            rhomolar, smolar, hmolar, T, p, umolar = [], [], [], [], [], []
            for _T in np.logspace(np.log10(HEOS.keyed_output(CP.iT_triple)), np.log10(HEOS.keyed_output(CP.iT_critical)), 500):
                try:
                    HEOS.update(CP.QT_INPUTS, Q, _T)
                    if (HEOS.p() < 0): raise ValueError('P is negative:' + str(HEOS.p()))
                    HEOS.T(), HEOS.p(), HEOS.rhomolar(), HEOS.hmolar(), HEOS.smolar()
                    HEOS.umolar()

                    T.append(HEOS.T())
                    p.append(HEOS.p())
                    rhomolar.append(HEOS.rhomolar())
                    hmolar.append(HEOS.hmolar())
                    smolar.append(HEOS.smolar())
                    umolar.append(HEOS.umolar())
                except ValueError as VE:
                    myprint(1, 'satT error:', VE, '; T:', '{T:0.16g}'.format(T=_T), 'T/Tc:', _T / HEOS.keyed_output(CP.iT_critical))

            dic.update(dict(T=np.array(T),
                            P=np.array(p),
                            Dmolar=np.array(rhomolar),
                            Hmolar=np.array(hmolar),
                            Smolar=np.array(smolar),
                            Umolar=np.array(umolar)))

    def plot_saturation_curves(self):
        for ax in self.axes_list:
            ax.label_axes()
            ax.plot_saturation_curves()

    def calc_Tmax_curve(self):
        HEOS = CP.AbstractState(self.backend, self.fluid)
        rhomolar, smolar, hmolar, T, p, umolar = [], [], [], [], [], []

        for _p in np.logspace(np.log10(HEOS.keyed_output(CP.iP_min) * 1.01), np.log10(HEOS.keyed_output(CP.iP_max)), 300):
            try:
                HEOS.update(CP.PT_INPUTS, _p, HEOS.keyed_output(CP.iT_max))
            except ValueError as VE:
                myprint(1, 'Tmax', _p, VE)
                continue

            try:
                T.append(HEOS.T())
                p.append(HEOS.p())
                rhomolar.append(HEOS.rhomolar())
                hmolar.append(HEOS.hmolar())
                smolar.append(HEOS.smolar())
                umolar.append(HEOS.umolar())
            except ValueError as VE:
                myprint(1, 'Tmax access', VE)

        self.Tmax = dict(T=np.array(T),
                         P=np.array(p),
                         Dmolar=np.array(rhomolar),
                         Hmolar=np.array(hmolar),
                         Smolar=np.array(smolar),
                         Umolar=np.array(umolar))

    def plot_Tmax_curve(self):
        for ax in self.axes_list:
            ax.plot_Tmax_curve()

    def calc_melting_curve(self):
        state = CP.AbstractState('HEOS', self.fluid)
        rhomolar, smolar, hmolar, T, p, umolar = [], [], [], [], [], []

        # Melting line if it has it
        if state.has_melting_line():
            pmelt_min = max(state.melting_line(CP.iP_min, -1, -1), state.keyed_output(CP.iP_triple)) * 1.01
            pmelt_max = min(state.melting_line(CP.iP_max, -1, -1), state.keyed_output(CP.iP_max)) * 0.99

            for _p in np.logspace(np.log10(pmelt_min), np.log10(pmelt_max), 100):
                try:
                    Tm = state.melting_line(CP.iT, CP.iP, _p)
                    state.update(CP.PT_INPUTS, _p, Tm)
                    T.append(state.T())
                    p.append(state.p())
                    rhomolar.append(state.rhomolar())
                    hmolar.append(state.hmolar())
                    smolar.append(state.smolar())
                    umolar.append(state.umolar())
                except ValueError as VE:
                    myprint(1, 'melting', VE)

        self.melt = dict(T=np.array(T),
                         P=np.array(p),
                         Dmolar=np.array(rhomolar),
                         Hmolar=np.array(hmolar),
                         Smolar=np.array(smolar),
                         Umolar=np.array(umolar))

    def plot_melting_curve(self):
        for ax in self.axes_list:
            ax.plot_melting_curve()

    def tight_layout(self):
        self.fig.tight_layout()

    def add_to_pdf(self, pdf):
        """ Add this figure to the pdf instance """
        pdf.savefig(self.fig)

    def savefig(self, fname, **kwargs):
        self.fig.savefig(fname, **kwargs)


class ConsistencyAxis(object):
    def __init__(self, axis, fig, pair, fluid, backend, state1, state2, state3,
                 p_limits_1phase=None, T_limits_1phase=None, NT_1phase=40, Np_1phase=40,
                 NT_2phase=20, NQ_2phase=20
                 ):
        self.ax = axis
        self.fig = fig
        self.pair = pair
        self.fluid = fluid
        self.backend = backend
        self.state = state1
        self.state_PT = state2
        self.state_QT = state3
        self.p_limits_1phase = p_limits_1phase
        self.T_limits_1phase = T_limits_1phase
        self.NT_1phase = NT_1phase
        self.Np_1phase = Np_1phase
        self.NQ_2phase = NQ_2phase
        self.NT_2phase = NT_2phase
        # self.saturation_curves()

    def label_axes(self):
        """ Label the axes for the given pair """
        xparam, yparam = split_pair_xy(self.pair)
        self.ax.set_xlabel(param_labels[xparam])
        self.ax.set_ylabel(param_labels[yparam])

        if xparam in ['P', 'Dmolar']:
            self.ax.set_xscale('log')
        if yparam in ['P', 'Dmolar']:
            self.ax.set_yscale('log')

    def plot_saturation_curves(self):
        xparam, yparam = split_pair_xy(self.pair)
        xL = self.to_axis_units(xparam, self.fig.dictL[xparam])
        yL = self.to_axis_units(yparam, self.fig.dictL[yparam])
        xV = self.to_axis_units(xparam, self.fig.dictV[xparam])
        yV = self.to_axis_units(yparam, self.fig.dictV[yparam])
        self.ax.plot(xL, yL, 'k', lw=1)
        self.ax.plot(xV, yV, 'k', lw=1)

    def plot_Tmax_curve(self):
        xparam, yparam = split_pair_xy(self.pair)
        x = self.to_axis_units(xparam, self.fig.Tmax[xparam])
        y = self.to_axis_units(yparam, self.fig.Tmax[yparam])
        self.ax.plot(x, y, 'r', lw=1)

    def plot_melting_curve(self):
        xparam, yparam = split_pair_xy(self.pair)
        x = self.to_axis_units(xparam, self.fig.melt[xparam])
        y = self.to_axis_units(yparam, self.fig.melt[yparam])
        self.ax.plot(x, y, 'b', lw=1)

    def to_axis_units(self, label, vals):
        """ Convert to the units used in the plot """
        if label in ['Hmolar', 'Smolar', 'Umolar', 'Dmolar', 'P']:
            return vals / 1000
        elif label in ['T']:
            return vals
        else:
            raise ValueError(label)

    def consistency_check_singlephase(self):

        tic = time.time()

        # Update the state given the desired set of inputs
        param1, param2 = split_pair(self.pair)
        key1 = getattr(CP, 'i' + param1)
        key2 = getattr(CP, 'i' + param2)
        pairkey = getattr(CP, self.pair + '_INPUTS')

        # Get the keys and indices and values for the inputs needed
        xparam, yparam = split_pair_xy(self.pair)
        xkey = getattr(CP, 'i' + xparam)
        ykey = getattr(CP, 'i' + yparam)

        data = []

        if self.p_limits_1phase is not None:
            # User-specified limits were provided, use them
            p_min, p_max = self.p_limits_1phase
        else:
            # No user-specified limits were provided, use the defaults
            p_min = self.state.keyed_output(CP.iP_min) * 1.01
            p_max = self.state.keyed_output(CP.iP_max)

        for p in np.logspace(np.log10(p_min), np.log10(p_max), self.Np_1phase):

            if self.T_limits_1phase is None:
                # No user-specified limits were provided, using the defaults
                Tmin = self.state.keyed_output(CP.iT_triple)
                if self.state.has_melting_line():
                    try:
                        pmelt_min = self.state.melting_line(CP.iP_min, -1, -1)
                        if p < pmelt_min:
                            T0 = Tmin
                        else:
                            T0 = self.state.melting_line(CP.iT, CP.iP, p)
                    except Exception as E:
                        T0 = Tmin + 1.1
                        data.append(dict(err=str(E), type="melting", input=p))
                        myprint(1, 'MeltingLine:', E)
                else:
                    T0 = Tmin + 1.1
                Tvec = np.linspace(T0, self.state.keyed_output(CP.iT_max), self.NT_1phase)
            else:
                # Use the provided limits for T
                Tvec = np.linspace(self.T_limits_1phase[0], self.T_limits_1phase[1], self.NT_1phase)

            for T in Tvec:

                try:
                    # Update the state using PT inputs in order to calculate all the remaining inputs
                    self.state_PT.update(CP.PT_INPUTS, p, T)
                except ValueError as VE:
                    data.append(dict(err=str(VE), cls="EXCEPTION", type="update", in1="P", val1=p, in2="T", val2=T))
                    myprint(1, 'consistency', VE)
                    continue

                _exception = False
                tic2 = timeit.default_timer()
                try:
                    val1, val2 = self.state_PT.keyed_output(key1), self.state_PT.keyed_output(key2)
                    self.state.update(pairkey, val1, val2)
                    toc2 = timeit.default_timer()
                except ValueError as VE:
                    data.append(dict(err=str(VE), cls="EXCEPTION", type="update", in1=param1, val1=val1, in2=param2, val2=val2))
                    myprint(1, 'update(1p)', self.pair, 'P', p, 'T', T, 'D', self.state_PT.keyed_output(CP.iDmolar), '{0:18.16g}, {1:18.16g}'.format(self.state_PT.keyed_output(key1), self.state_PT.keyed_output(key2)), VE)
                    _exception = True

                x = self.to_axis_units(xparam, self.state_PT.keyed_output(xkey))
                y = self.to_axis_units(yparam, self.state_PT.keyed_output(ykey))

                if not _exception:
                    # Check the error on the density
                    if abs(self.state_PT.rhomolar() / self.state.rhomolar() - 1) < 1e-3 and abs(self.state_PT.p() / self.state.p() - 1) < 1e-3 and abs(self.state_PT.T() - self.state.T()) < 1e-3:
                        data.append(dict(cls="GOOD", x=x, y=y, elapsed=toc2 - tic2))
                        if 'REFPROP' not in self.backend:
                            if self.state_PT.phase() != self.state.phase():
                                myprint(1, 'bad phase', self.pair, '{0:18.16g}, {1:18.16g}'.format(self.state_PT.keyed_output(key1), self.state_PT.keyed_output(key2)), self.state.phase(), 'instead of', self.state_PT.phase())
                    else:
                        data.append(dict(cls="INCONSISTENT", type="update", in1=param1, val1=val1, in2=param2, val2=val2, x=x, y=y))
                        myprint(1, 'bad', self.pair, '{0:18.16g}, {1:18.16g}'.format(self.state_PT.keyed_output(key1), self.state_PT.keyed_output(key2)), 'T:', self.state_PT.T(), 'Drho:', abs(self.state_PT.rhomolar() / self.state.rhomolar() - 1), abs(self.state_PT.p() / self.state.p() - 1), 'DT:', abs(self.state_PT.T() - self.state.T()))

        toc = time.time()
        df = pandas.DataFrame(data)
        bad = df[df.cls == 'INCONSISTENT']
        good = df[df.cls == 'GOOD']
        slowgood = good[good.elapsed > 0.01]
        excep = df[df.cls == 'EXCEPTION']
        badphase = df[df.cls == 'BAD_PHASE']

        self.ax.plot(bad.x, bad.y, 'r+', ms=3)
        self.ax.plot(good.x, good.y, 'k.', ms=1)
        self.ax.plot(excep.x, excep.y, 'rx', ms=3)
        self.ax.plot(slowgood.x, slowgood.y, 'b*', ms=6)
        self.ax.plot(badphase.x, badphase.y, 'o', ms=3, mfc='none')

        print('1-phase took ' + str(toc - tic) + ' s for ' + self.pair)

        if self.pair == 'HmolarSmolar':
            # plt.plot(good.elapsed)
            # plt.title(self.pair)
            # plt.show()

            good.to_excel('times_water.xlsx')
        return df[df.cls != 'GOOD']

    def consistency_check_twophase(self):

        tic = time.time()
        state = self.state

        try:
            if state.fluid_param_string('pure') == 'false':
                print("Not a pure-fluid, skipping two-phase evaluation")
                return
        except:
            pass

        # Update the state given the desired set of inputs
        param1, param2 = split_pair(self.pair)
        key1 = getattr(CP, 'i' + param1)
        key2 = getattr(CP, 'i' + param2)
        pairkey = getattr(CP, self.pair + '_INPUTS')

        # Get the keys and indices and values for the inputs needed
        xparam, yparam = split_pair_xy(self.pair)
        xkey = getattr(CP, 'i' + xparam)
        ykey = getattr(CP, 'i' + yparam)

        data = []
        for q in np.linspace(0, 1, self.NQ_2phase):

            Tmin = state.keyed_output(CP.iT_triple) + 1

            for T in np.linspace(Tmin, state.keyed_output(CP.iT_critical) - 1, self.NT_2phase):

                try:
                    # Update the state using QT inputs in order to calculate all the remaining inputs
                    self.state_QT.update(CP.QT_INPUTS, q, T)
                except ValueError as VE:
                    data.append(dict(err=str(VE), cls="EXCEPTION", type="update", in1="Q", val1=q, in2="T", val2=T))
                    myprint(1, 'consistency', VE)
                    continue

                _exception = False
                try:
                    val1, val2 = self.state_QT.keyed_output(key1), self.state_QT.keyed_output(key2)
                    state.update(pairkey, val1, val2)
                except ValueError as VE:
                    data.append(dict(err=str(VE), cls="EXCEPTION", type="update", in1=param1, val1=val1, in2=param2, val2=val2))
                    myprint(1, 'update_QT', T, q)
                    myprint(1, 'update', param1, self.state_QT.keyed_output(key1), param2, self.state_QT.keyed_output(key2), VE)
                    _exception = True

                x = self.to_axis_units(xparam, self.state_QT.keyed_output(xkey))
                y = self.to_axis_units(yparam, self.state_QT.keyed_output(ykey))

                if not _exception:
                    # Check the error on the density
                    if abs(self.state_QT.rhomolar() / self.state.rhomolar() - 1) < 1e-3 and abs(self.state_QT.p() / self.state.p() - 1) < 1e-3 and abs(self.state_QT.T() - self.state.T()) < 1e-3:
                        data.append(dict(cls="GOOD", x=x, y=y))
                        if 'REFPROP' not in self.backend:
                            if self.state_QT.phase() != self.state.phase():
                                myprint(1, 'bad phase (2phase)', self.pair, '{0:18.16g}, {1:18.16g}'.format(self.state_QT.keyed_output(key1), self.state_QT.keyed_output(key2)), self.state.phase(), 'instead of', self.state_QT.phase())

                    else:
                        myprint(1, 'Q', q)
                        myprint(1, 'bad(2phase)', self.pair, '{0:18.16g}, {1:18.16g}'.format(self.state_QT.keyed_output(key1), self.state_QT.keyed_output(key2)), 'pnew:', self.state.p(), 'pold:', self.state_QT.p(), 'Tnew:', self.state.T(), 'T:', self.state_QT.T(), 'Drho:', abs(self.state_QT.rhomolar() / self.state.rhomolar() - 1), 'DP', abs(self.state_QT.p() / self.state.p() - 1), 'DT:', abs(self.state_QT.T() - self.state.T()))
                        data.append(dict(cls="INCONSISTENT", type="update", in1=param1, val1=val1, in2=param2, val2=val2, x=x, y=y))

        toc = time.time()
        df = pandas.DataFrame(data)
        bad = df[df.cls == 'INCONSISTENT']
        good = df[df.cls == 'GOOD']
        excep = df[df.cls == 'EXCEPTION']
        badphase = df[df.cls == 'BAD_PHASE']

        self.ax.plot(bad.x, bad.y, 'r+', ms=3)
        self.ax.plot(good.x, good.y, 'k.', ms=1)
        self.ax.plot(excep.x, excep.y, 'rx', ms=3)
        self.ax.plot(badphase.x, badphase.y, 'o', ms=3, mfc='none')
        print('2-phase took ' + str(toc - tic) + ' s for ' + self.pair)

    def cross_out_axis(self):
        xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()
        self.ax.plot([xlims[0], xlims[1]], [ylims[0], ylims[1]], lw=3, c='r')
        self.ax.plot([xlims[0], xlims[1]], [ylims[1], ylims[0]], lw=3, c='r')

        xparam, yparam = split_pair_xy(self.pair)
        x = 0.5 * xlims[0] + 0.5 * xlims[1]
        y = 0.5 * ylims[0] + 0.5 * ylims[1]
        if xparam in ['P', 'Dmolar']:
            x = (xlims[0] * xlims[1])**0.5
        if yparam in ['P', 'Dmolar']:
            y = (ylims[0] * ylims[1])**0.5

        self.ax.text(x, y, 'Not\nImplemented', ha='center', va='center', bbox=dict(fc='white'))


if __name__ == '__main__':
    PVT = PdfPages('Consistency.pdf')
    CP.CoolProp.set_debug_level(0)
    open('timelog.txt', 'w')
    with open('timelog.txt', 'a+', buffering=1) as fp:
        for fluid in ['METHANOL']:  # CP.__fluids__:
            tic = timeit.default_timer()
            skips = ['DmolarHmolar', 'DmolarSmolar', 'DmolarUmolar', 'HmolarSmolar']
            skips = []
            ff = ConsistencyFigure(fluid, backend='HEOS', additional_skips=skips)  # , NT_1phase = 10, Np_1phase = 10, NT_2phase = 100, NQ_2phase = 0)
            ff.errors.to_excel('Errors' + fluid + '.xlsx')
            toc = timeit.default_timer()
            print('Time to build:', toc - tic, 'seconds')
            ff.add_to_pdf(PVT)
            ff.savefig(fluid + '.png')
            ff.savefig(fluid + '.pdf')
            plt.close()
            fp.write('Time to build: {0} seconds for {1}\n'.format(toc - tic, fluid))
            del ff
    PVT.close()
