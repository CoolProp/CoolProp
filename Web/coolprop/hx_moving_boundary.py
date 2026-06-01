"""Provider-based moving-boundary counterflow heat-exchanger solver.

Faithful Python port of the original paper supplemental script
(``dev/reference/HX.py``) for

    I.H. Bell et al., "A generalized moving-boundary algorithm to predict the
    heat transfer rate of counterflow heat exchangers for any phase
    configuration," Applied Thermal Engineering 79 (2015) 192-201.

Every thermophysical property call routes through a :class:`PropertyProvider`
so the identical solver runs against either the full Helmholtz EOS (``HEOS``)
or the SVD-compressed tabular backend (``SVDSBTL&HEOS``).  Used by
``SVDSBTLHeatExchangerDemo.ipynb``.
"""
from __future__ import division
from math import log

import numpy as np
import scipy.optimize
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PT_INPUTS, HmassP_INPUTS, PQ_INPUTS, QT_INPUTS

# Headline Q [W] for the Table-3 evaporator at A = 4 m^2, from the recovered
# oracle dev/reference/HX.py on CoolProp 7.x.  Regression anchor for the port.
ORACLE_Q = 4577.242219


class PropertyProvider(object):
    """Property access for one fluid through a chosen backend.

    The ``p,h <-> T,rho`` transforms (:meth:`h_pT`, :meth:`Ts_ph`) go through
    the chosen backend; saturation (:meth:`Tsat`, :meth:`hsat_TQ`) always goes
    through a HEOS ancillary state, so saturation sourcing is identical for both
    backends and the only timed difference is the table lookup.
    """

    backend = None  # set by subclasses

    def __init__(self, fluid):
        self.AS = CP.AbstractState(self.backend, fluid)
        self.sat = CP.AbstractState('HEOS', fluid)

    def h_pT(self, p, T):
        self.AS.update(PT_INPUTS, p, T)
        return self.AS.hmass()

    def Ts_ph(self, p, h):
        self.AS.update(HmassP_INPUTS, h, p)
        return self.AS.T(), self.AS.smass()

    def Tsat(self, p, Q):
        self.sat.update(PQ_INPUTS, p, Q)
        return self.sat.T()

    def hsat_TQ(self, T, Q):
        self.sat.update(QT_INPUTS, Q, T)
        return self.sat.hmass()


class HEOSProvider(PropertyProvider):
    backend = 'HEOS'


class SVDSBTLProvider(PropertyProvider):
    backend = 'SVDSBTL&HEOS'


class HeatExchanger(object):
    """Moving-boundary counterflow HX, faithful to dev/reference/HX.py."""

    def __init__(self, prov_h, prov_c, mdot_h, p_hi, h_hi, mdot_c, p_ci, h_ci):
        self.ph, self.pc = prov_h, prov_c
        self.mdot_h, self.p_hi, self.h_hi = mdot_h, p_hi, h_hi
        self.mdot_c, self.p_ci, self.h_ci = mdot_c, p_ci, h_ci
        self.T_ci, _ = self.pc.Ts_ph(p_ci, h_ci)
        self.T_hi, _ = self.ph.Ts_ph(p_hi, h_hi)
        self.T_cbubble = self.pc.Tsat(p_ci, 0)
        self.T_cdew = self.pc.Tsat(p_ci, 1)
        self.T_hbubble = self.ph.Tsat(p_hi, 0)
        self.T_hdew = self.ph.Tsat(p_hi, 1)
        self.h_cbubble = self.pc.hsat_TQ(self.T_cbubble, 0)
        self.h_cdew = self.pc.hsat_TQ(self.T_cdew, 1)
        self.h_hbubble = self.ph.hsat_TQ(self.T_hbubble, 0)
        self.h_hdew = self.ph.hsat_TQ(self.T_hdew, 1)

    def external_pinching(self):
        self.h_ho = self.ph.h_pT(self.p_hi, self.T_ci)        # Eq 5
        Qmaxh = self.mdot_h * (self.h_hi - self.h_ho)         # Eq 4
        self.h_co = self.pc.h_pT(self.p_ci, self.T_hi)        # Eq 7
        Qmaxc = self.mdot_c * (self.h_co - self.h_ci)         # Eq 6
        Qmax = min(Qmaxh, Qmaxc)
        self.calculate_cell_boundaries(Qmax)
        return Qmax

    def calculate_cell_boundaries(self, Q):
        self.h_co = self.h_ci + Q / self.mdot_c
        self.h_ho = self.h_hi - Q / self.mdot_h
        self.hvec_c = [self.h_ci, self.h_co]
        self.hvec_h = [self.h_ho, self.h_hi]
        if self.h_hi > self.h_hdew > self.h_ho:
            self.hvec_h.insert(-1, self.h_hdew)
        if self.h_hi > self.h_hbubble > self.h_ho:
            self.hvec_h.insert(1, self.h_hbubble)
        if self.h_ci < self.h_cdew < self.h_co:
            self.hvec_c.insert(-1, self.h_cdew)
        if self.h_ci < self.h_cbubble < self.h_co:
            self.hvec_c.insert(1, self.h_cbubble)
        k = 0
        while k < len(self.hvec_c) - 1 or k < len(self.hvec_h) - 1:
            if len(self.hvec_c) == 2 and len(self.hvec_h) == 2:
                break
            Qcell_hk = self.mdot_h * (self.hvec_h[k + 1] - self.hvec_h[k])
            Qcell_ck = self.mdot_c * (self.hvec_c[k + 1] - self.hvec_c[k])
            if abs(Qcell_hk / Qcell_ck - 1) < 1e-6:
                k += 1
                break
            elif Qcell_hk > Qcell_ck:
                self.hvec_h.insert(k + 1, self.hvec_h[k] + Qcell_ck / self.mdot_h)
            else:
                self.hvec_c.insert(k + 1, self.hvec_c[k] + Qcell_hk / self.mdot_c)
            k += 1
        assert len(self.hvec_h) == len(self.hvec_c)
        self.Tvec_c = np.array([self.pc.Ts_ph(self.p_ci, h)[0] for h in self.hvec_c])
        self.Tvec_h = np.array([self.ph.Ts_ph(self.p_hi, h)[0] for h in self.hvec_h])
        self.phases_h = self._phases(self.hvec_h, self.h_hbubble, self.h_hdew)
        self.phases_c = self._phases(self.hvec_c, self.h_cbubble, self.h_cdew)

    @staticmethod
    def _phases(hvec, hbubble, hdew):
        out = []
        for i in range(len(hvec) - 1):
            havg = (hvec[i] + hvec[i + 1]) / 2.0
            if havg < hbubble:
                out.append('liquid')
            elif havg > hdew:
                out.append('vapor')
            else:
                out.append('two-phase')
        return out

    def internal_pinching(self, stream):
        if stream == 'hot':
            for i in range(1, len(self.hvec_h) - 1):
                if abs(self.hvec_h[i] - self.h_hdew) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]:
                    h_c_pinch = self.pc.h_pT(self.p_ci, self.T_hdew)       # Eq 10
                    Qright = self.mdot_h * (self.h_hi - self.h_hdew)       # Eq 9
                    Qmax = self.mdot_c * (h_c_pinch - self.h_ci) + Qright  # Eq 12
                    self.calculate_cell_boundaries(Qmax)
                    return Qmax
        elif stream == 'cold':
            for i in range(1, len(self.hvec_c) - 1):
                if abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]:
                    h_h_pinch = self.ph.h_pT(self.p_hi, self.T_cbubble)    # Eq 14
                    Qleft = self.mdot_c * (self.h_cbubble - self.h_ci)     # Eq 13
                    Qmax = Qleft + self.mdot_h * (self.h_hi - h_h_pinch)   # Eq 16
                    self.calculate_cell_boundaries(Qmax)
                    return Qmax
        else:
            raise ValueError('stream must be "hot" or "cold"')
        return None

    def objective_function(self, Q):
        self.calculate_cell_boundaries(Q)
        w = []
        for k in range(len(self.hvec_c) - 1):
            DTA = self.Tvec_h[k + 1] - self.Tvec_c[k + 1]
            DTB = self.Tvec_h[k] - self.Tvec_c[k]
            LMTD = DTA if DTA == DTB else (DTA - DTB) / log(abs(DTA / DTB))
            UA_req = self.mdot_h * (self.hvec_h[k + 1] - self.hvec_h[k]) / LMTD
            alpha_c = 100 if self.phases_c[k] in ('liquid', 'vapor') else 1000
            alpha_h = 100 if self.phases_h[k] in ('liquid', 'vapor') else 1000
            UA_avail = 1 / (1 / (alpha_h * self.A_h) + 1 / (alpha_c * self.A_c))
            w.append(UA_req / UA_avail)
        return 1 - sum(w)

    def solve(self):
        self.Q = scipy.optimize.brentq(
            self.objective_function, 1e-5, self.Qmax - 1e-10, rtol=1e-14, xtol=1e-10)
        return self.Q

    def run(self, and_solve=True):
        self.Qmax = self.external_pinching()
        for stream in ('hot', 'cold'):
            qi = self.internal_pinching(stream)
            if qi is not None:
                self.Qmax = qi
        return self.solve() if and_solve else None


def make_evaporator(backend, A=4.0):
    """Construct the water-heated n-Propane evaporator of the paper's Table 3.

    Matches the recovered oracle: the HOT stream (Water) carries mdot = 0.1 kg/s
    and the COLD stream (n-Propane) mdot = 0.01 kg/s -- transposed relative to
    the paper's printed Table 3 (see the design spec and dev/reference/HX.py).
    ``backend`` is 'HEOS' or 'SVDSBTL&HEOS' ('SVDSBTL' is accepted as an alias).
    """
    if backend == 'HEOS':
        ph, pc = HEOSProvider('Water'), HEOSProvider('n-Propane')
    elif backend in ('SVDSBTL', 'SVDSBTL&HEOS'):
        ph, pc = SVDSBTLProvider('Water'), SVDSBTLProvider('n-Propane')
    else:
        raise ValueError("backend must be 'HEOS' or 'SVDSBTL&HEOS', got %r" % (backend,))
    p_w = 101325.0
    h_w = ph.h_pT(p_w, 330.0)
    pc.sat.update(QT_INPUTS, 1, 300.0)
    p_r = pc.sat.p()
    h_r = pc.h_pT(p_r, 275.0)
    hx = HeatExchanger(ph, pc, mdot_h=0.1, p_hi=p_w, h_hi=h_w,
                       mdot_c=0.01, p_ci=p_r, h_ci=h_r)
    hx.A_h = hx.A_c = A
    return hx
