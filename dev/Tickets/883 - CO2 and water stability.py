#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI, AbstractState
import sys

back = "HEOS"
fluids = ["water", "CO2"]

num_p = 10
num_h = 500

msgs = ""

for fluid in fluids:
    state = AbstractState(back, fluid)
    fluid = back + "::" + fluid
    p_lo = 1.05 * PropsSI("pmin", "", 0, "", 0, fluid)
    p_hi = 0.95 * PropsSI("pcrit", "", 0, "", 0, fluid)
    p_r = np.logspace(np.log10(p_lo), np.log10(p_hi), num_p)
    h_r = []; d_r = []; d_o = []
    for p in p_r:
        if state.has_melting_line():
            T_m = state.melting_line(CoolProp.constants.iT, CoolProp.constants.iP, p)
            h_m = -np.Inf
            dT = 0.0
            while dT < 2.0 and not np.isfinite(h_m):
                dT += 0.1
                T_t = T_m + dT
                try:
                    h_m = PropsSI("H", "P", p, "T", T_t, fluid)
                except:
                    pass
        else:
            h_m = -np.Inf
        h1 = PropsSI("H", "P", p, "Q", 0, fluid)
        dh = 0.15 * (PropsSI("H", "P", p, "Q", 1, fluid) - h1)
        h0 = np.maximum(h1 - dh, h_m)
        h = np.linspace(h0, h1, num=num_h)
        err = [""] * h.size
        d = np.empty_like(h)
        for i, h_i in enumerate(h):
            try:
                d[i] = PropsSI("D", "P", p, "H", h_i, fluid)
                err[i] = ""
            except Exception as e:
                d[i] = np.NaN
                err[i] = str(e)

        err = np.array(err)
        valid = np.isfinite(d)
        invalid = np.logical_not(valid)
        for i in range(np.sum(invalid)):
            msgs += "[Error] {0} at h={1} J/kg, p={2} Pa: {3}\n".format(fluid, h[invalid][i], p, err[invalid][i])

        print("[Error] {0} at p={1:8.4f} bar: Failed for {2:4.1f}% of the calls.".format(fluid, p / 1e5, np.sum(invalid) / h.size * 1e2))

# print(msgs)

sys.exit()
