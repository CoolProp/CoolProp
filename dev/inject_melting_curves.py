import os, json


Simon_curves = {
    "n-Propane": {
        "BibTeX": "Reeves-JCP-1964", "T_m": -187.75 + 273.15, "parts": [{"T_0": 85.3, "a": 7.180e8, "c": 1.283, "p_0": 0.0, "T_max": 168.63}]
    },
    "n-Pentane": {
        "BibTeX": "Reeves-JCP-1964", "T_m": -129.89 + 273.15, "parts": [{"T_0": 143.5, "a": 6.600e8, "c": 1.649, "p_0": 0.0, "T_max": 156.2}]
    },
    "Isopentane": {
        "BibTeX": "Reeves-JCP-1964", "T_m": -159.92 + 273.15, "parts": [{"T_0": 112.5, "a": 5.916e8, "c": 1.563, "p_0": 0, "T_max": 212.16}]
    },
    "Propylene": {
        "BibTeX": "Reeves-JCP-1964", "T_m": -185.09 + 273.15, "parts": [{"T_0": 86.0, "a": 3.196e8, "c": 2.821, "p_0": 0, "T_min": 86.0, "T_max": 129},
                                                                         {"T_0": 109.6, "a": 3.064e8, "c": 3.871, "p_0": 4.450e8, "T_min": 129, "T_max": 145.3}]
    },
    "Cyclohexane": {
        "BibTeX": "Penoncello-IJT-1995", "T_m": 6.81 + 273.15, "parts": [{"T_0": 279.7, "a": 383.4e6, "c": 1.41, "p_0": 0, "T_max": 401.7}]
    },
    "Krypton": {
        "BibTeX": "Michels-PHYSICA-1962", "T_m": 115.95, "parts": [{"T_0": 1, "a": 109479.2307, "c": 1.6169841, "p_0": -237497645.7, "T_max": 168.7}]
    },
    "Xenon": {
        "BibTeX": "Michels-PHYSICA-1962", "T_m": 165.02, "parts": [{"T_0": 1, "a": 80890.5544859, "c": 1.5891650, "p_0": -260932309.446, "T_max": 366.4}]
    },
    "CarbonMonoxide": {
        "BibTeX": "Barreiros-JCT-1982", "T_m": 68.3, "parts": [{"T_0": 1, "a": 19560.8, "c": 2.10747, "p_0": -142921439.2, "T_max": 87.5}]
    },
    "Oxygen": {
        "BibTeX": "Younglove-NIST-1982", "T_m": 54.75, "parts": [{"T_0": 1, "a": 227606.348, "c": 1.769, "p_0": -266999247.652, "T_max": 63.1}]
    },
    "ParaHydrogen": {
        "BibTeX": "Younglove-NIST-1982", "T_m": 18.9, "parts": [{"T_0": 1, "a": 125746.643, "c": 1.955, "p_0": -21155737.752, "T_min": 13.8033, "T_max": 22},
                                                    {"T_0": 1, "a": 248578.596, "c": 1.764739, "p_0": -26280332.904, "T_min": 22, "T_max": 164.5}]
    },
    "Methane": {
        "BibTeX": "Abramson-HPR-2011", "T_m": 90.7, "parts": [{"T_0": 90.6941, "a": 0.208e9, "c": 1.698, "p_0": 1.17e4, "T_max": 600}]
    },
    "Helium": {
        "BibTeX": "Datchi-PRB-2000", "T_m": 1.15, "parts": [{"T_0": 1, "a": 1.6067e6, "c": 1.565, "p_0": -1.6067e6, "T_max": 700}]
    },
    "Neon": {
        "BibTeX": "SantamariaPerez-PRB-2010", "T_m": -1, "parts": [{"T_0": 24.4, "a": 1.7e9, "c": 1 / 0.77, "p_0": 101325, "T_max": 700}]
    },
    "Hydrogen": {
        "BibTeX": "Datchi-PRB-2000", "T_m": 14.009985, "parts": [{"T_0": 1, "a": 2.31e5, "c": 1.7627, "p_0": -0.0052e6 - 2.31e5, "T_max": 700}]
    }
}

polynomial_in_Tr = {
    "Argon": {
        "BibTeX": "Tegeler-JPCRD-1999", "T_m": 87.28, "parts": [{"T_0": 83.8058, "a": [-7476.2665, 9959.0613], "t": [1.05, 1.275], "p_0": 68891, "T_max": 254.0}]
    },
    "Fluorine": {
        "BibTeX": "deReuck-BOOK-1990", "T_m": 53.15, "parts": [{"T_0": 53.4811, "a": [988043.478261], "t": [2.1845], "p_0": 252, "T_max": 55.4}]
    },
    "Nitrogen": {
        "BibTeX": "Span-JPCRD-2000", "T_m": 77.34, "parts": [{"T_0": 63.151, "a": [12798.61], "t": [1.78963], "p_0": 12523, "T_max": 283.8}]
    },
    "Ethane": {
        "BibTeX": "Buecker-JCRD-2006", "T_m": 90.4, "parts": [{"T_0": 90.368, "a": [2.23626315e8, 1.05262374e8], "t": [1.0, 2.55], "p_0": 1.14, "T_max": 110.2}]
    },
    "Isobutane": {
        "BibTeX": "Buecker-JPCRD-2006B", "T_m": 113.55, "parts": [{"T_0": 113.73, "a": [1.9536371309e9], "t": [6.12], "p_0": 0.0219, "T_max": 124.9}]
    },
    "Ethylene": {
        "BibTeX": "Smukala-JPCRD-2000", "T_m": 169, "parts": [{"T_0": 103.989, "a": [2947001.84], "t": [2.045], "p_0": 122.65, "T_min": 103.989, "T_max": 110.369},
                                                               {"T_0": 110.369, "a": [6.82693421], "t": [1.089], "p_0": 46.8e6, "T_min": 110.369, "T_max": 188}]
    },
    "n-Butane": {
        "BibTeX": "Buecker-JPCRD-2006B", "T_m": -137.92 + 273.15, "parts": [{"T_0": 134.895, "a": [5.585582364e8], "t": [2.206], "p_0": 0.653, "T_max": 163.9}]
    },
    "Water": {
        "BibTeX": "IAPWS", "T_m": -1, "parts": [{"T_0": 273.16, "a": [-0.119539337e7, -0.808183159e5, -0.333826860e4], "t": [0.3000000e1, 0.257500e2, 0.103750e3], "p_0": 611.657, "T_min": 273.16, "T_max": 251.165},
                                                 {"T_0": 251.165, "a": [0.299948], "t": [60], "p_0": 208.566e6, "T_min": 251.165, "T_max": 256.164},
                                                 {"T_0": 256.164, "a": [1.18721], "t": [8], "p_0": 350.1e6, "T_min": 256.164, "T_max": 273.31},
                                                 {"T_0": 273.31, "a": [1.07476], "t": [4.6], "p_0": 623.4e6, "T_min": 273.31, "T_max": 355}
        ]
    }
}

polynomial_in_theta = {
    "Methanol": {
        "BibTeX": "deReuck-BOOK-1993", "T_m": 337.8, "parts": [{"T_0": 175.61, "a": [5.330770e9, 4.524780e9, 3.888861e10], "t": [1, 1.5, 4], "p_0": 0.187, "T_max": 245.9}]
    },
    "CarbonDioxide": {
        "BibTeX": "Span-JPCRD-1996", "T_m": 216.58, "parts": [{"T_0": 216.592, "a": [1955.5390, 2055.4593], "t": [1, 2], "p_0": 517950, "T_max": 327.6}]
    }
}

import CoolProp
__ = 0
for fluid in CoolProp.__fluids__:
    if fluid not in Simon_curves and fluid not in polynomial_in_Tr and fluid not in polynomial_in_theta:
        print(fluid)
        __ += 1
    else:
        print(' ' * 30, fluid)
print(__)

import CoolProp.CoolProp as CP
import json, numpy as np, matplotlib.pyplot as plt, pandas

ip = 1
irho = 1
Nrow, Ncol = 5, 5

figp = plt.figure(figsize=(20, 20))
figrho = plt.figure(figsize=(20, 20))


def plot_rho(T, rho, fit=False):
    x, y = (T - T[0]) / (T[len(T) - 1] - T[0]), (rho - rho[0]) / (rho[len(rho) - 1] - rho[0])

    c = np.polyfit(x, y, 3)
    yfit = np.polyval(c, x)
    err = yfit - y
    rms = np.sqrt(np.mean(np.power(err, 2)))

    rhofit = yfit * (rho[len(rho) - 1] - rho[0]) + rho[0]
    if fit:
        return T, (rhofit / rho - 1) * 100
    else:
        return x, y


def simon():
    global ip, irho

    for fluid, values in Simon_curves.iteritems():
        axp = figp.add_subplot(Nrow, Ncol, ip); ip += 1
        axrho = figrho.add_subplot(Nrow, Ncol, irho); irho += 1
        axp.set_xlabel('T [K]')
        axp.set_ylabel('p [Pa]')
        axrho.set_xlabel('T [K]')
        axrho.set_ylabel('rho [mol/m$^3$]')
        axp.set_title(fluid + ' - ' + str(round(CP.Props(fluid, "molemass"), 2)))
        axrho.set_title(fluid)

        fname = os.path.join('fluids', fluid + '.json')
        j = json.load(open(fname, 'r'))
        for part in values['parts']:
            if 'T_min' not in part:
                part['T_min'] = round(CP.Props(fluid, "Tmin"), 4)
        values['type'] = 'Simon'

        j['ANCILLARIES']['melting_line'] = values

        fp = open(fname, 'w')
        from package_json import json_options
        fp.write(json.dumps(j, **json_options))
        fp.close()

#         if not isinstance(values, list):
#             values = [values]
#             df = pandas.read_csv('melting_curves/'+fluid+'.mlt',names=['T','p','rho'])
#             axp.plot(df['T'], df['p'], 'o', mfc='none')
#             x,y = plot_rho(df['T'],df['rho'],fit = True)
#             axrho.plot(x,y, 'o', mfc='none')
#         else:
#             for i in ['I','II']:
#                 df = pandas.read_csv('melting_curves/'+fluid+'-'+i+'.mlt',names=['T','p','rho'])
#                 axp.plot(df['T'], df['p'], 'o', mfc='none')
#                 x,y = plot_rho(df['T'],df['rho'],fit = True)
#                 axrho.plot(x,y, 'o', mfc='none')

        T_m = values['T_m']
        for i, value in enumerate(values['parts']):

            Tmin = value.get('T_min', CP.Props(fluid, "Tmin"))
            Tmax = value['T_max']

            T = np.linspace(Tmin, Tmax, 200)
            T_0 = value['T_0']
            p_0 = value['p_0']
            a = value['a']
            c = value['c']

            p = p_0 + a * ((T / T_0)**c - 1)

            axp.plot(T, p)

            cc = 1.75
            aa = 3e8  # (101325-p_0)/((T_m/T_0)**cc-1)
            pt = CP.Props(fluid, 'ptriple')
            pp = pt + aa * ((T / Tmin)**cc - 1)
            axp.plot(T_m, 101325, '*')
            axp.plot(T, pp, '--')

            print("%s %s %s %s" % (fluid, CP.Props(fluid, "molemass"), CP.Props(fluid, 'accentric'), pp[-1] / p[-1] - 1))

#             if fluid == 'Helium':
#                 T = np.array([326.2,345.1,362.8,385.1,419.4,459,499,535.7,570,608])
#                 p = p_0 + a*((T/T_0)**c - 1)
#                 print p


def Tr():
    global ip, irho
    for fluid, values in polynomial_in_Tr.iteritems():

        axp = figp.add_subplot(Nrow, Ncol, ip); ip += 1
        axrho = figrho.add_subplot(Nrow, Ncol, irho); irho += 1
        axp.set_xlabel('T [K]')
        axp.set_ylabel('p [Pa]')
        axrho.set_xlabel('T [K]')
        axrho.set_ylabel('rho [mol/m$^3$]')
        axp.set_title(fluid + ' - ' + str(round(CP.Props(fluid, "molemass"), 2)))
        axrho.set_title(fluid)

        fname = os.path.join('fluids', fluid + '.json')
        j = json.load(open(fname, 'r'))
        for part in values['parts']:
            if 'T_min' not in part:
                part['T_min'] = round(CP.Props(fluid, "Tmin"), 4)
        values['type'] = 'polynomial_in_Tr'

        j['ANCILLARIES']['melting_line'] = values

        fp = open(fname, 'w')
        from package_json import json_options
        fp.write(json.dumps(j, **json_options))
        fp.close()

        if fluid == 'Ethylene':
            T = [104.003, 104.059, 104.13, 104.2, 104.27, 104.41, 104.55, 104.69, 104.83, 104.969, 105.108, 105.386, 106.077, 106.764, 107.446, 111.384, 119.283, 127.136, 158.146, 188.621]
            p = np.array([0.1, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 25, 50, 75, 100, 200, 300]) * 1e6

            axp.plot(T, p, '*')

#         if not isinstance(values, list):
#             values = [values]
#             df = pandas.read_csv('melting_curves/'+fluid+'.mlt',names=['T','p','rho'])
#             axp.plot(df['T'], df['p'], 'o', mfc='none')
#             x,y = plot_rho(df['T'],df['rho'],fit = True)
#             axrho.plot(x,y, 'o', mfc='none')
#
#         else:
#             for i in ['I','II']:
#                 df = pandas.read_csv('melting_curves/'+fluid+'-'+i+'.mlt',names=['T','p','rho'])
#                 axp.plot(df['T'], df['p'], 'o', mfc='none')
#                 x,y = plot_rho(df['T'],df['rho'],fit = True)
#                 axrho.plot(x,y, 'o', mfc='none')

        T_m = values['T_m']
        for i, value in enumerate(values['parts']):

            Tmin = value.get('T_min', CP.Props(fluid, "Tmin"))
            Tmax = value['T_max']
            T = np.linspace(Tmin, Tmax, 200)

            a = value['a']
            t = value['t']
            T_t = value['T_0']
            p_t = value['p_0']

            RHS = 0
            for i in range(len(a)):
                RHS += a[i] * ((T / T_t)**t[i] - 1)

            p = p_t * (RHS + 1)

            axp.plot(T, p)

            cc = 1.75
            aa = 3e8  # (101325-p_0)/((T_m/T_0)**cc-1)
            pt = CP.Props(fluid, 'ptriple')
            pp = pt + aa * ((T / Tmin)**cc - 1)
            axp.plot(T_m, 101325, '*')
            axp.plot(T, pp, '--')

            print("%s %s %s %s" % (fluid, CP.Props(fluid, "molemass"), CP.Props(fluid, 'accentric'), pp[-1] / p[-1] - 1))


def theta():
    global ip, irho
    for fluid, values in polynomial_in_theta.iteritems():

        axp = figp.add_subplot(Nrow, Ncol, ip); ip += 1
        axrho = figrho.add_subplot(Nrow, Ncol, irho); irho += 1
        axp.set_xlabel('T [K]')
        axp.set_ylabel('p [Pa]')
        axrho.set_xlabel('T [K]')
        axrho.set_ylabel('rho [mol/m$^3$]')
        axp.set_title(fluid + ' - ' + str(round(CP.Props(fluid, "molemass"), 2)))
        axrho.set_title(fluid)

        fname = os.path.join('fluids', fluid + '.json')
        j = json.load(open(fname, 'r'))
        for part in values['parts']:
            if 'T_min' not in part:
                part['T_min'] = round(CP.Props(fluid, "Tmin"), 4)
        values['type'] = 'polynomial_in_Theta'

        j['ANCILLARIES']['melting_line'] = values

        fp = open(fname, 'w')
        from package_json import json_options
        fp.write(json.dumps(j, **json_options))
        fp.close()

        T_m = values['T_m']
        for value in values['parts']:

            a = value['a']
            t = value['t']
            T_t = value['T_0']
            p_t = value['p_0']

            Tmin = T_t
            Tmax = value['T_max']
            T = np.linspace(Tmin, Tmax, 200)

            RHS = 0
            for i in range(len(a)):
                RHS += a[i] * (T / T_t - 1)**t[i]

            p = p_t * (RHS + 1)

            #df = pandas.read_csv('melting_curves/' + fluid + '.mlt', names=['T','p','rho'])

            #axp.plot(df['T'], df['p'], 'o', mfc='none')

            axp.plot(T, p)

            #x,y = plot_rho(df['T'],df['rho'],fit = True)
            #axrho.plot(x,y, 'o', mfc='none')

            cc = 1.75
            aa = 3e8  # (101325-p_0)/((T_m/T_0)**cc-1)
            pt = CP.Props(fluid, 'ptriple')
            pp = pt + aa * ((T / Tmin)**cc - 1)
            axp.plot(T_m, 101325, '*')
            axp.plot(T, pp, '--')

            print("%s %s %s %s" % (fluid, CP.Props(fluid, "molemass"), CP.Props(fluid, 'accentric'), pp[-1] / p[-1] - 1))


if __name__ == '__main__':
    simon()
    Tr()
    theta()

    figp.tight_layout()
    figrho.tight_layout()

    figp.savefig('p.pdf')
    figrho.savefig('rho.pdf')

    plt.close()
