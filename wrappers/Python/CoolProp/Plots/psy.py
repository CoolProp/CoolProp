#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import cPickle
from ConfigParser import ConfigParser

import numpy as np

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from pylab import Figure

from CoolProp.HumidAirProp import HAProps, HAProps_Aux

from PyQt4.QtGui import (QDialog, QGridLayout, QProgressBar, QLabel,
                         QDialogButtonBox, QPushButton, QFileDialog, QApplication)

Preferences = ConfigParser()
config_path = os.path.join(os.path.dirname(__file__), "psyrc")
Preferences.read(config_path)
P = Preferences.getfloat("General", "P")


def _Pbar(Z):
    """
    ASHRAE Fundamentals Handbook pag 1.1 eq. 3
    input:
        Z: altitude, m
    return
        standard atmosphere barometric pressure, Pa
    """
    return 101325. * (1 - 2.25577e-5 * Z)**5.256


class PsychroPlot(FigureCanvasQTAgg):
    """
    Plot widget for psychrometric chart
        Add custom margins
        Define a pointer to text state properties, to remove and redraw
    """

    def __init__(self, parent=None, width=15, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        self.axes2D = self.fig.add_subplot(111)
        FigureCanvasQTAgg.updateGeometry(self)
        self.axes2D.figure.subplots_adjust(left=0.01, right=0.92,
                                           bottom=0.05, top=0.98)
        self.notes = []

    def plot(self, *args, **kwargs):
        self.axes2D.plot(*args, **kwargs)

    def config(self):
        self.axes2D.set_autoscale_on(False)
        self.axes2D.set_xlabel(u"Tdb, ºC")
        self.axes2D.set_ylabel("Absolute humidity, kg/kg")
        self.axes2D.yaxis.set_ticks_position("right")
        self.axes2D.yaxis.set_label_position("right")

        tmin = Preferences.getfloat("Psychr", "isotdbStart") - 273.15
        tmax = Preferences.getfloat("Psychr", "isotdbEnd") - 273.15

        self.axes2D.set_xlim(tmin, tmax)
        self.axes2D.set_ylim(0, 0.04)

    def showPointData(self, state):
        self.clearPointData()

        yi = 0.99
        keys = "tdb", "tdp", "twb", "HR", "w", "h", "v", "rho"
        units = u"ºC", u"ºC", u"ºC", "%", "kgw/kgda", "kJ/kg", u"m³/kg", u"kg/m³"
        for key, txt in zip(keys, units):
            self.notes.append(self.axes2D.annotate(
                "%s: %0.4f %s" % (key, state.__getattribute__(key), txt), (0.01, yi),
                xycoords='axes fraction', size="small", va="center"))
            yi -= 0.025
        self.draw()

    def clearPointData(self):
        while self.notes:
            anotation = self.notes.pop()
            anotation.remove()
        self.draw()


class PsyCoolprop(object):
    """
    Psychrometric state using coolprop library

    kwargs definition parameters:
        P: Pressure, Pa
        z: altitude, m

        tdp: dew-point temperature
        tdb: dry-bulb temperature
        twb: web-bulb temperature
        w: Humidity Ratio [kg water/kg dry air]
        HR: Relative humidity
        h: Mixture enthalpy
        v: Mixture specified volume

    P: mandatory input for barometric pressure, z is an alternate pressure input
    it needs other two input parameters:
        0 - tdb, w
        1 - tdb, HR
        2 - tdb, twb
        3 - tdb, tdp
        4 - tdp, HR
        5 - tdp, twb
        6 - twb, w
    """
    kwargs = {"z": 0.0,
              "P": 0.0,

              "tdb": 0.0,
              "tdb": 0.0,
              "twb": 0.0,
              "w": None,
              "HR": None,
              "h": None,
              "v": 0.0}
    status = 0
    msg = "Unknown variables"

    def __init__(self, **kwargs):
        self.kwargs = self.__class__.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()
            self.msg = "Solved"

    @property
    def calculable(self):
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

        self._mode = 0
        if tdb and w is not None:
            self._mode = ("Tdb", "W")
        elif tdb and HR is not None:
            self._mode = ("Tdb", "RH")
        elif tdb and twb:
            self._mode = ("Tdb", "Twb")
        elif tdb and tdp:
            self._mode = ("Tdb", "Tdp")
        elif tdp and HR is not None:
            self._mode = ("Tdp", "RH")

        return bool(self._mode)

    def calculo(self):
        tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h = self._lib()
        self.tdp = tdp - 273.15
        self.tdb = tdb - 273.15
        self.twb = twb - 273.15
        self.P = P
        self.Pvs = Pvs
        self.Pv = Pv
        self.ws = ws
        self.w = w
        self.HR = HR
        self.mu = w / ws * 100
        self.v = v
        self.rho = 1 / v
        self.h = h
        self.Xa = 1 / (1 + self.w / 0.62198)
        self.Xw = 1 - self.Xa

    def args(self):
        # Correct coolprop custom namespace versus pychemqt namespace
        if "Tdb" in self._mode:
            self.kwargs["Tdb"] = self.kwargs["tdb"]
        if "Twb" in self._mode:
            self.kwargs["Twb"] = self.kwargs["twb"]
        if "Tdp" in self._mode:
            self.kwargs["Tdp"] = self.kwargs["tdp"]
        if "RH" in self._mode:
            self.kwargs["RH"] = self.kwargs["HR"]
        if "W" in self._mode:
            self.kwargs["W"] = self.kwargs["w"]

        var1 = self.kwargs[self._mode[0]]
        var2 = self.kwargs[self._mode[1]]

        # units conversion to coolprop expected unit:
        # HR in 0-1, H in kJ/kg, S in kJ/kgK
        if "RH" in self._mode[0]:
            var1 /= 100.
        if "RH" in self._mode[1]:
            var2 /= 100.

        args = ("P", self._P_kPa, self._mode[0], var1, self._mode[1], var2)
        return args

    def _P(self):
        """Barometric pressure calculation, Pa"""
        if self.kwargs["P"]:
            P = self.kwargs["P"]
        elif self.kwargs["z"]:
            P = _Pbar(self.kwargs["z"])
        else:
            P = 101325.
        return P

    @property
    def _P_kPa(self):
        """Property for ease access to pressure in kPa"""
        P = self._P()
        return P / 1000.

    def _lib(self):
        args = self.args()
        P = self._P()

        if "Tdb" in self._mode:
            tdb = self.kwargs["Tdb"]
        else:
            tdb = HAProps("Tdb", *args)
        tdp = HAProps("Tdp", *args)
        twb = HAProps("Twb", *args)
        w = HAProps("W", *args)
        HR = HAProps("RH", *args) * 100
        Pvs = HAProps_Aux("p_ws", tdb, self._P_kPa, w)[0] * 1000
        Pv = Pvs * HR / 100
        ws = HAProps("W", "P", self._P_kPa, "Tdb", tdb, "RH", 1)
        v = HAProps("V", *args)
        h = HAProps("H", *args)

        return tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h

    @classmethod
    def calculatePlot(cls, parent):
        """Function to calculate points in chart"""

        data = {}
        P = Preferences.getfloat("General", "P")
        P_kPa = P / 1000
        t = cls.LineList("isotdb", Preferences)

        # Saturation line
        Hs = []
        for tdb in t:
            Hs.append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", 1))
            parent.progressBar.setValue(5 * len(Hs) / len(t))
        data["t"] = t
        data["Hs"] = Hs

        # left limit of isow lines
        H = cls.LineList("isow", Preferences)
        th = []
        for w in H:
            if w:
                tdp = HAProps("Tdp", "P", 101.325, "W", w, "RH", 1)
                th.append(tdp - 273.15)
            else:
                tmin = Preferences.getfloat("Psychr", "isotdbStart")
                th.append(tmin - 273.15)
        data["H"] = H
        data["th"] = th

        # Humidity ratio lines
        hr = cls.LineList("isohr", Preferences)
        Hr = {}
        cont = 0
        for i in hr:
            Hr[i] = []
            for tdb in t:
                Hr[i].append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", i / 100.))
                cont += 1
                parent.progressBar.setValue(5 + 10 * cont / len(hr) / len(Hs))
        data["Hr"] = Hr

        # Twb
        lines = cls.LineList("isotwb", Preferences)
        Twb = {}
        cont = 0
        for T in lines:
            ws = HAProps("W", "P", P_kPa, "RH", 1, "Tdb", T)
            H = [ws, 0]
            Tw = [T - 273.15, HAProps("Tdb", "P", P_kPa, "Twb", T, "RH", 0) - 273.15]
            cont += 1
            parent.progressBar.setValue(15 + 75 * cont / len(lines))
            Twb[T] = (H, Tw)
        data["Twb"] = Twb

        # v
        lines = cls.LineList("isochor", Preferences)
        V = {}
        rh = np.arange(1, -0.05, -0.05)
        for cont, v in enumerate(lines):
            w = []
            Td = []
            for r in rh:
                w.append(HAProps("W", "P", P_kPa, "RH", r, "V", v))
                Td.append(HAProps("Tdb", "P", P_kPa, "RH", r, "V", v) - 273.15)
            parent.progressBar.setValue(90 + 10 * cont / len(lines))
            V[v] = (Td, w)
        data["v"] = V

        return data

    @staticmethod
    def LineList(name, Preferences):
        """Return a list with the values of isoline name to plot"""
        if Preferences.getboolean("Psychr", name + "Custom"):
            t = []
            for i in Preferences.get("Psychr", name + 'List').split(','):
                if i:
                    t.append(float(i))
        else:
            start = Preferences.getfloat("Psychr", name + "Start")
            end = Preferences.getfloat("Psychr", name + "End")
            step = Preferences.getfloat("Psychr", name + "Step")
            t = np.arange(start, end, step)
        return t


class UI_Psychrometry(QDialog):
    """Psychrometric charts tool"""

    def __init__(self, parent=None):
        super(UI_Psychrometry, self).__init__(parent)
        self.showMaximized()
        self.setWindowTitle("Psychrometric chart")

        layout = QGridLayout(self)
        self.diagrama2D = PsychroPlot(self, dpi=90)
        self.diagrama2D.fig.canvas.mpl_connect('motion_notify_event', self.scroll)
        layout.addWidget(self.diagrama2D, 1, 1, 1, 2)
        self.progressBar = QProgressBar()
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar, 2, 1)
        self.status = QLabel()
        layout.addWidget(self.status, 2, 1)

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Close)
        butonPNG = QPushButton("Save as PNG")
        self.buttonBox.addButton(butonPNG, QDialogButtonBox.AcceptRole)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.accepted.connect(self.savePNG)
        layout.addWidget(self.buttonBox, 2, 2)

        self.plot()

    def savePNG(self):
        """Save chart image to png file"""
        fname = unicode(QFileDialog.getSaveFileName(
            self, "Save chart to file",
            "./", "Portable Network Graphics (*.png)"))
        self.diagrama2D.fig.savefig(fname, facecolor='#eeeeee')

    def drawlabel(self, name, Preferences, t, W, label, unit):
        """
        Draw annotation for isolines
            name: name of isoline
            Preferences: Configparse instance of pychemqt preferences
            t: x array of line
            W: y array of line
            label: text value to draw
            unit: text units to draw
        """
        if Preferences.getboolean("Psychr", name + "label"):
            tmin = Preferences.getfloat("Psychr", "isotdbStart") - 273.15
            tmax = Preferences.getfloat("Psychr", "isotdbEnd") - 273.15
            x = tmax - tmin
            wmin = Preferences.getfloat("Psychr", "isowStart")
            wmax = Preferences.getfloat("Psychr", "isowEnd")
            y = wmax - wmin

            i = 0
            for ti, wi in zip(t, W):
                if tmin <= ti <= tmax and wmin <= wi <= wmax:
                    i += 1
            label = str(label)
            if Preferences.getboolean("Psychr", name + "units"):
                label += unit
            pos = Preferences.getfloat("Psychr", name + "position")
            p = int(i * pos / 100 - 1)
            rot = np.arctan((W[p] - W[p - 1]) / y / (t[p] - t[p - 1]) * x) * 360.0 / 2.0 / np.pi
            self.diagrama2D.axes2D.annotate(label, (t[p], W[p]),
                rotation=rot, size="small", ha="center", va="center")

    def plot(self):
        """Plot chart"""
        Preferences = ConfigParser()
        Preferences.read("psyrc")

        self.diagrama2D.axes2D.clear()
        self.diagrama2D.config()
        filename = "%i.pkl" % P
        if os.path.isfile(filename):
            with open(filename, "r") as archivo:
                data = cPickle.load(archivo)
                self.status.setText("Loading cached data...")
                QApplication.processEvents()
        else:
            self.progressBar.setVisible(True)
            self.status.setText("Calculating data, be patient...")
            QApplication.processEvents()
            data = PsyCoolprop.calculatePlot(self)
            cPickle.dump(data, open(filename, "w"))
            self.progressBar.setVisible(False)
        self.status.setText("Plotting...")
        QApplication.processEvents()

        tmax = Preferences.getfloat("Psychr", "isotdbEnd") - 273.15

        t = [ti - 273.15 for ti in data["t"]]
        Hs = data["Hs"]
        format = {}
        format["ls"] = Preferences.get("Psychr", "saturationlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "saturationlineWidth")
        format["color"] = Preferences.get("Psychr", "saturationColor")
        format["marker"] = Preferences.get("Psychr", "saturationmarker")
        format["markersize"] = 3
        self.diagrama2D.plot(t, Hs, **format)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isotdblineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isotdblineWidth")
        format["color"] = Preferences.get("Psychr", "isotdbColor")
        format["marker"] = Preferences.get("Psychr", "isotdbmarker")
        format["markersize"] = 3
        for i, T in enumerate(t):
            self.diagrama2D.plot([T, T], [0, Hs[i]], **format)

        H = data["H"]
        th = data["th"]
        format = {}
        format["ls"] = Preferences.get("Psychr", "isowlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isowlineWidth")
        format["color"] = Preferences.get("Psychr", "isowColor")
        format["marker"] = Preferences.get("Psychr", "isowmarker")
        format["markersize"] = 3
        for i, H in enumerate(H):
            self.diagrama2D.plot([th[i], tmax], [H, H], **format)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isohrlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isohrlineWidth")
        format["color"] = Preferences.get("Psychr", "isohrColor")
        format["marker"] = Preferences.get("Psychr", "isohrmarker")
        format["markersize"] = 3
        for Hr, H0 in data["Hr"].iteritems():
            self.diagrama2D.plot(t, H0, **format)
            self.drawlabel("isohr", Preferences, t, H0, Hr, "%")

        format = {}
        format["ls"] = Preferences.get("Psychr", "isotwblineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isotwblineWidth")
        format["color"] = Preferences.get("Psychr", "isotwbColor")
        format["marker"] = Preferences.get("Psychr", "isotwbmarker")
        format["markersize"] = 3
        for T, (H, Tw) in data["Twb"].iteritems():
            self.diagrama2D.plot(Tw, H, **format)
            value = T - 273.15
            txt = u"ºC"
            self.drawlabel("isotwb", Preferences, Tw, H, value, txt)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isochorlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isochorlineWidth")
        format["color"] = Preferences.get("Psychr", "isochorColor")
        format["marker"] = Preferences.get("Psychr", "isochormarker")
        format["markersize"] = 3
        for v, (Td, H) in data["v"].iteritems():
            self.diagrama2D.plot(Td, H, **format)
            value = v
            txt = u"m³/kg"
            self.drawlabel("isochor", Preferences, Td, H, value, txt)

        self.diagrama2D.draw()
        self.status.setText("P = %i Pa" % P)

    def scroll(self, event):
        """Update graph annotate when mouse scroll over chart"""
        if event.xdata and event.ydata:
            punto = self.createState(event.xdata, event.ydata)
            if event.ydata < punto.ws:
                self.diagrama2D.showPointData(punto)
            else:
                self.diagrama2D.clearPointData()

    def createState(self, x, y):
        """Create psychrometric state from click or mouse position"""
        tdb = x + 273.15
        punto = PsyCoolprop(P=P, tdb=tdb, w=y)
        return punto


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    aireHumedo = UI_Psychrometry()
    aireHumedo.show()
    sys.exit(app.exec_())
