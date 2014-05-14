Simon_curves = {
    "n-Propane" : {
        "T_0" : 85.3, "a" : 7.180e8, "c" : 1.283, "p_0" : 0.0, "T_max" : 168.63, "BibTeX" : "Reeves-JCP-1964"
    },
    "n-Pentane" : {
        "T_0" : 143.5, "a" : 6.600e8, "c" : 1.649, "p_0" : 0.0, "T_max" : 156.2, "BibTeX" : "Reeves-JCP-1964"
    },
    "Isopentane" : {
        "T_0" : 112.5, "a" : 5.916e8, "c" : 1.563, "p_0" : 0, "T_max" : 212.16, "BibTeX" : "Reeves-JCP-1964"
    },
    "Propylene" : [
        {
            "T_0" : 86.0, "a" : 3.196e8, "c" : 2.821, "p_0" : 0, "T_min": 86.0, "T_max" : 129, "BibTeX" : "Reeves-JCP-1964"
        },
        {
            "T_0" : 109.6, "a" : 3.064e8, "c" : 3.871, "p_0" : 4.450e8, "T_min": 129, "T_max" : 145.3, "BibTeX" : "Reeves-JCP-1964"
        }
    ],
    "Cyclohexane" : {
        "T_0" : 279.7, "a" : 383.4e6, "c" : 1.41, "p_0" : 0, "T_max" : 401.7, "BibTeX" : "Penoncello-IJT-1995"
    },
    "Krypton" : {
        "T_0" : 1, "a" : 109479.2307, "c" : 1.6169841, "p_0" : -237497645.7, "T_max" : 168.7, "BibTeX" : "Michels-PHYSICA-1962"
    },
    "Xenon" : {
        "T_0" : 1, "a" : 80890.5544859, "c" : 1.5891650, "p_0" : -260932309.446, "T_max" : 366.4, "BibTeX" : "Michels-PHYSICA-1962"
    },
    "CarbonMonoxide" : {
        "T_0" : 1, "a" : 19560.8, "c" : 2.10747, "p_0" : -142921439.2, "T_max" : 87.5, "BibTeX" : "Barreiros-JCT-1982"
    },
    "Oxygen": {
        "T_0" : 1, "a" : 227606.348, "c" : 1.769, "p_0" : -266999247.652, "T_max" : 63.1, "BibTeX" : "Younglove-NIST-1982"
    },
    "ParaHydrogen": [
    {
        "T_0" : 1, "a" : 125746.643, "c" : 1.955, "p_0" : -21155737.752, "T_min" : 13.8033, "T_max" : 22, "BibTeX" : "Younglove-NIST-1982"
    },
    {
        "T_0" : 1, "a" : 248578.596, "c" : 1.764739, "p_0" : -26280332.904, "T_min" : 22, "T_max" : 164.5, "BibTeX" : "Younglove-NIST-1982"
    }
    ]
}

polynomial_in_Tr = {
    "Argon" : {
        "T_t" : 83.8058, "a" : [-7476.2665, 9959.0613], "t" : [1.05,1.275], "p_t" : 68891, "T_max" : 254.0, "BibTeX" : "Tegeler-JPCRD-1999"
    },
    "Fluorine" : {
        "T_t" : 53.4811, "a" : [988043.478261], "t" : [2.1845], "p_t" : 252, "T_max" : 55.4, "BibTeX" : "deReuck-BOOK-1990"
    },
    "Nitrogen" : {
        "T_t" : 63.151, "a" : [12798.61], "t" : [1.78963], "p_t" : 12523, "T_max" : 283.8, "BibTeX" : "Span-JPCRD-2000"
    },
    "Ethane" : {
        "T_t" : 90.368, "a" : [2.23626315e8, 1.05262374e8], "t" : [1.0, 2.55], "p_t" : 1.14, "T_max" : 110.2, "BibTeX" : "Buecker-JCRD-2006"
    },
    "Isobutane" : {
        "T_t" : 113.73, "a" : [1.9536371309e9], "t" : [6.12], "p_t" : 0.0219, "T_max" : 124.9, "BibTeX" : "Buecker-JPCRD-2006B"
    },
    "Ethylene" : [
    {
        "T_t" : 103.989, "a" : [2947001.84], "t" : [2.045], "p_t" : 122.65, "T_min" : 103.989, "T_max" : 110.369, "BibTeX" : "Smukala-JPCRD-2000"
    },
    {
        "T_t" : 110.369, "a" : [6.82693421], "t" : [1.089], "p_t" : 46.8e6, "T_min" : 110.369, "T_max" : 188, "BibTeX" : "Smukala-JPCRD-2000"
    }
    ],
    "n-Butane" : {
        "T_t" : 134.895, "a" : [5.585582364e8], "t" : [2.206], "p_t" : 0.653, "T_max" : 163.9, "BibTeX" : "Buecker-JPCRD-2006B"
    }
}

polynomial_in_theta = {
    "Methanol" : {
        "T_t" : 175.61, "a" : [5.330770e9, 4.524780e9, 3.888861e10], "t" : [1, 1.5, 4], "p_t" : 0.187, "T_max" : 245.9, "BibTeX" : "deReuck-BOOK-1993"
    },
    "CarbonDioxide" : {
        "T_t" : 216.592, "a" : [1955.5390, 2055.4593], "t" : [1, 2], "p_t" : 517950, "T_max" : 327.6, "BibTeX" : "Span-JPCRD-1996"
    }
}

import json, numpy as np, matplotlib.pyplot as plt, pandas

ip = 1
irho = 1
Nrow,Ncol = 5,5

figp = plt.figure(figsize = (20,20))
figrho = plt.figure(figsize = (20,20))

def plot_rho(T, rho, fit = False):
    x, y = (T-T[0])/(T[len(T)-1]-T[0]), (rho-rho[0])/(rho[len(rho)-1]-rho[0])
    
    c = np.polyfit(x, y, 3)
    yfit = np.polyval(c, x)
    err = yfit - y
    rms = np.sqrt(np.mean(np.power(err,2)))

    rhofit = yfit*(rho[len(rho)-1]-rho[0])+rho[0]
    if fit:
        return T, (rhofit/rho-1)*100
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
        axp.set_title(fluid)
        axrho.set_title(fluid)
        
        if not isinstance(values, list):
            values = [values]
            df = pandas.read_csv('melting_curves/'+fluid+'.mlt',names=['T','p','rho'])
            axp.plot(df['T'], df['p'], 'o', mfc='none')
            x,y = plot_rho(df['T'],df['rho'],fit = True)
            axrho.plot(x,y, 'o', mfc='none')
        else:
            for i in ['I','II']:
                df = pandas.read_csv('melting_curves/'+fluid+'-'+i+'.mlt',names=['T','p','rho'])
                axp.plot(df['T'], df['p'], 'o', mfc='none')
                x,y = plot_rho(df['T'],df['rho'],fit = True)
                axrho.plot(x,y, 'o', mfc='none')
        
        for i,value in enumerate(values):
            Tmin = value.get('T_min',min(df['T']))
            Tmax = value['T_max']
            
            T = np.linspace(Tmin, Tmax, 200)
            T_0 = value['T_0']
            p_0 = value['p_0']
            a = value['a']
            c = value['c']
            
            p = p_0 + a*((T/T_0)**c - 1)
            
            axp.plot(T, p)

def Tr():
    global ip, irho
    for fluid, values in polynomial_in_Tr.iteritems():
        
        axp = figp.add_subplot(Nrow, Ncol, ip); ip += 1
        axrho = figrho.add_subplot(Nrow, Ncol, irho); irho += 1
        axp.set_xlabel('T [K]')
        axp.set_ylabel('p [Pa]')
        axrho.set_xlabel('T [K]')
        axrho.set_ylabel('rho [mol/m$^3$]')
        axp.set_title(fluid)
        axrho.set_title(fluid)
        
        if fluid == 'Ethylene':
            T = [104.003, 104.059, 104.13, 104.2, 104.27, 104.41, 104.55, 104.69, 104.83, 104.969, 105.108, 105.386, 106.077, 106.764, 107.446, 111.384, 119.283, 127.136, 158.146, 188.621]
            p = np.array([0.1, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 25, 50, 75, 100, 200, 300])*1e6
            
            axp.plot(T,p,'*')
            
        if not isinstance(values, list):
            values = [values]
            df = pandas.read_csv('melting_curves/'+fluid+'.mlt',names=['T','p','rho'])
            axp.plot(df['T'], df['p'], 'o', mfc='none')
            x,y = plot_rho(df['T'],df['rho'],fit = True)
            axrho.plot(x,y, 'o', mfc='none')
        
        else:
            for i in ['I','II']:
                df = pandas.read_csv('melting_curves/'+fluid+'-'+i+'.mlt',names=['T','p','rho'])
                axp.plot(df['T'], df['p'], 'o', mfc='none')
                x,y = plot_rho(df['T'],df['rho'],fit = True)
                axrho.plot(x,y, 'o', mfc='none')
                
        for i,value in enumerate(values):
            Tmin = value.get('T_min',min(df['T']))
            Tmax = value['T_max']
            T = np.linspace(Tmin, Tmax, 200)
                
            a = value['a']
            t = value['t']
            T_t = value['T_t']
            p_t = value['p_t']
            
            RHS = 0
            for i in range(len(a)):
                RHS += a[i]*((T/T_t)**t[i] - 1)
        
            p = p_t*(RHS + 1)
            
            axp.plot(T, p)
        
    
def theta():
    global ip, irho
    for fluid, value in polynomial_in_theta.iteritems():
        
        axp = figp.add_subplot(Nrow, Ncol, ip); ip += 1
        axrho = figrho.add_subplot(Nrow, Ncol, irho); irho += 1
        axp.set_xlabel('T [K]')
        axp.set_ylabel('p [Pa]')
        axrho.set_xlabel('T [K]')
        axrho.set_ylabel('rho [mol/m$^3$]')
        axp.set_title(fluid)
        axrho.set_title(fluid)
        
        a = value['a']
        t = value['t']
        T_t = value['T_t']
        p_t = value['p_t']
        
        Tmin = T_t
        Tmax = value['T_max']
        T = np.linspace(Tmin, Tmax, 200)
        
        RHS = 0
        for i in range(len(a)):
            RHS += a[i]*(T/T_t - 1)**t[i]
        
        p = p_t*(RHS + 1)
        
        df = pandas.read_csv('melting_curves/' + fluid + '.mlt', names=['T','p','rho'])
        
        axp.plot(df['T'], df['p'], 'o', mfc='none')
        
        axp.plot(T, p)
        
        x,y = plot_rho(df['T'],df['rho'],fit = True)
        axrho.plot(x,y, 'o', mfc='none')
        
if __name__=='__main__':
    simon()
    Tr()
    theta()

    figp.tight_layout()
    figrho.tight_layout()
    
    figp.savefig('p.pdf')
    figrho.savefig('rho.pdf')
    
    plt.close()