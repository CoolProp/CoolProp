from __future__ import division
import numpy as np
from scipy.odr import *
import scipy.optimize, random
import matplotlib.pyplot as plt
import textwrap
import random
from CoolProp.CoolProp import Props, get_REFPROPname


def rsquared(x, y):
    """
    Return R^2 where x and y are array-like.

    from http://stackoverflow.com/questions/893657/how-do-i-calculate-r-squared-using-python-and-numpy
    """
    import scipy.stats
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2


def saturation_density(Ref, ClassName, form='A', LV='L', perc_error_allowed=0.3, fName=None, add_critical=True):
    """

    Parameters
    ----------
    Ref : string
        The fluid name for the fluid that will be used to generate the saturation data
    ClassName : The name of the class that will be used in the C++ code
    form : string
        If ``'A'``, use a term of the form
    """

    if fName is None:
        Tc = Props(Ref, 'Tcrit')
        pc = Props(Ref, 'pcrit')
        rhoc = Props(Ref, 'rhocrit')
        Tmin = Props(Ref, 'Tmin')
        print("%s %s %s" % (Ref, Tmin, Props(Ref, 'Ttriple')))

        TT = np.linspace(Tmin, Tc - 1, 1000)
        TT = list(TT)
        # Add temperatures around the critical temperature
        if add_critical:
            for dT in [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
                TT.append(Tc - dT)
        TT = np.array(sorted(TT))

        p = Props('P', 'T', TT, 'Q', 0, Ref)
        rhoL = Props('D', 'T', TT, 'Q', 0, Ref)
        rhoV = Props('D', 'T', TT, 'Q', 1, Ref)
    else:
        Tc = 423.27
        pc = 3533
        rhoc = 470
        Tmin = 273
        lines = open(fName, 'r').readlines()
        TT, p, rhoL, rhoV = [], [], [], []
        for line in lines:
            _T, _p, _rhoL, _rhoV = line.split(' ')
            TT.append(float(_T))
            p.append(float(_p))
            rhoL.append(float(_rhoL))
            rhoV.append(float(_rhoV))

        TT = np.array(TT)

    # Start with a large library of potential powers
    n = [i / 6.0 for i in range(1, 200)]  # +[0.35+i/200.0 for i in range(1,70)]+[0.05+0.01*i for i in range(1,70)]

    x = 1.0 - TT / Tc

    if LV == 'L':
        rho_EOS = rhoL
    elif LV == 'V':
        rho_EOS = rhoV
    else:
        raise ValueError

    if form == 'A':
        y = np.array(rho_EOS) / rhoc - 1
    elif form == 'B':
        y = (np.log(rho_EOS) - np.log(rhoc)) * TT / Tc
    else:
        raise ValueError

    max_abserror = 0
    while len(n) > 3:
        print("%s %s" % (max_abserror, len(n)))

        def f_p(B, x):
            # B is a vector of the parameters.
            # x is an array of the current x values.
            return sum([_B * x**(_n) for _B, _n in zip(B, n)])

        linear = Model(f_p)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[0] * len(n))
        myoutput = myodr.run()

        beta = myoutput.beta
        sd = myoutput.sd_beta

        if form == 'A':
            rho_fit = (f_p(myoutput.beta, x) + 1) * rhoc
        elif form == 'B':
            rho_fit = np.exp(f_p(myoutput.beta, x) * Tc / TT) * rhoc
        else:
            raise ValueError

        print('first,last %s %s %s %s %s %s' % (TT[0], TT[-1], rho_fit[0], rho_fit[-1], rho_EOS[0], rho_EOS[-1]))

        max_abserror = np.max(np.abs(rho_fit / rho_EOS - 1)) * 100

        dropped_indices = [i for i in range(len(n)) if abs(sd[i]) < 1e-15]
        if dropped_indices:
            for i in reversed(sorted(dropped_indices)):
                n.pop(i)
            print('popping... %s terms remaining' % len(n))
            continue

        if max_abserror > perc_error_allowed:
            break  # The last good run will be used
        else:
            print(max_abserror)
            Ncoeffs = str(list(myoutput.beta)).lstrip('[').rstrip(']')
            tcoeffs = str(n).lstrip('[').rstrip(']')
            maxerror = max_abserror
            if form == 'A':
                code_template = textwrap.dedent(
                """
                for (int i=1; i<={count:d}; i++)
                {{
                    summer += N[i]*pow(theta,t[i]);
                }}
                return reduce.rho*(summer+1);
                """.format(count=len(n))
                )
            elif form == 'B':
                code_template = textwrap.dedent(
                """
                for (int i=1; i<={count:d}; i++)
                {{
                    summer += N[i]*pow(theta,t[i]);
                }}
                return reduce.rho*exp(reduce.T/T*summer);
                """.format(count=len(n))
                )
            else:
                raise ValueError

        # Find the least significant entry (the one with the largest relative standard error)
        # and remove it
        n.pop(np.argmax(np.abs(sd / beta)))

        # Remove elements that are not
    template = textwrap.dedent(
    """
    double {name:s}Class::rhosat{LV:s}(double T)
    {{
        // Maximum absolute error is {error:f} % between {Tmin:f} K and {Tmax:f} K
        const double t[] = {{0, {tcoeffs:s}}};
        const double N[] = {{0, {Ncoeffs:s}}};
        double summer=0,theta;
        theta=1-T/reduce.T;
        \t{code:s}
    }}
    """)
    the_string = template.format(tcoeffs=tcoeffs,
                            Ncoeffs=Ncoeffs,
                            name=ClassName,
                            Tmin=Tmin,
                            Tmax=TT[-1],
                            error=maxerror,
                            code=code_template,
                            LV=LV
                            )
    f = open('anc.txt', 'a')
    f.write(the_string)
    f.close()
    return the_string


def saturation_pressure_brute(Ref, ClassName):

    Tc = Props(Ref, 'Tcrit')
    pc = Props(Ref, 'pcrit')
    rhoc = Props(Ref, 'rhocrit')
    Tmin = Props(Ref, 'Tmin')

    TT = np.linspace(Tmin + 1e-6, Tc - 0.00001, 300)
    p = np.array([Props('P', 'T', T, 'Q', 0, Ref) for T in TT])
    rhoL = np.array([Props('D', 'T', T, 'Q', 0, Ref) for T in TT])
    rhoV = np.array([Props('D', 'T', T, 'Q', 1, Ref) for T in TT])

    Np = 10

    max_abserror = 99999
    bbest = []

    x = 1.0 - TT / Tc
    y = (np.log(rhoL) - np.log(rhoc)) * TT / Tc

    def f_p(B, x):
        # B is a vector of the parameters.
        # x is an array of the current x values.
        return sum([_B * x**(_n) for _B, _n in zip(B, b)])

    linear = Model(f_p)
    mydata = Data(x, y)

    for attempt in range(300):

        n = [i / 6.0 for i in range(1, 100)] + [0.35 + i / 200.0 for i in range(1, 70)] + [0.05 + 0.01 * i for i in range(1, 70)]
        b = []
        for _ in range(6):
            i = random.randint(0, len(n) - 1)
            b.append(n.pop(i))

        myodr = ODR(mydata, linear, beta0=[1] * len(b))
        myoutput = myodr.run()

        b = np.array(b)

        keepers = np.abs(myoutput.sd_beta / myoutput.beta) < 0.1
        if any(keepers):
            b = b[keepers]

        myodr = ODR(mydata, linear, beta0=[1] * len(b))
        myoutput = myodr.run()

        rho_fit = np.exp(f_p(myoutput.beta, x) * Tc / TT) * rhoc
        abserror = np.max(np.abs(rho_fit / rhoL - 1)) * 100
        print('.')
        if abserror < max_abserror:
            max_abserror = abserror
            bbest = b
            betabest = myoutput.beta
            print("%s %s %s" % (abserror, myoutput.sum_square, myoutput.sd_beta / myoutput.beta))


def saturation_pressure(Ref, ClassName, fName=None, LV=None):

    if fName is None:
        Tc = Props(Ref, 'Tcrit')
        pc = Props(Ref, 'pcrit')
        rhoc = Props(Ref, 'rhocrit')
        Tmin = Props(Ref, 'Tmin')

        TT = np.linspace(Tmin + 1e-6, Tc - 0.00001, 300)
        pL = Props('P', 'T', TT, 'Q', 0, Ref)
        pV = Props('P', 'T', TT, 'Q', 1, Ref)
        rhoL = Props('D', 'T', TT, 'Q', 0, Ref)
        rhoV = Props('D', 'T', TT, 'Q', 1, Ref)
    else:
        Tc = 423.27
        pc = 3533
        rhoc = 470
        Tmin = 273
        lines = open(fName, 'r').readlines()
        TT, p, rhoL, rhoV = [], [], [], []
        for line in lines:
            _T, _p, _rhoL, _rhoV = line.split(' ')
            TT.append(float(_T))
            p.append(float(_p))
            rhoL.append(float(_rhoL))
            rhoV.append(float(_rhoV))

        TT = np.array(TT)

    Np = 60
    n = range(1, Np)
    max_abserror = 0
    while len(n) > 3:

        def f_p(B, x):
            # B is a vector of the parameters.
            # x is an array of the current x values.
            return sum([_B * x**(_n / 2.0) for _B, _n in zip(B, n)])

        x = 1.0 - TT / Tc
        if LV == 'L':
            y = (np.log(pL) - np.log(pc)) * TT / Tc
        elif LV == 'V' or LV is None:
            y = (np.log(pV) - np.log(pc)) * TT / Tc

        linear = Model(f_p)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[0] * len(n))
        myoutput = myodr.run()

        beta = myoutput.beta
        sd = myoutput.sd_beta

        p_fit = np.exp(f_p(myoutput.beta, x) * Tc / TT) * pc
        if LV == 'L':
            max_abserror = np.max(np.abs((p_fit / pL) - 1) * 100)
        elif LV == 'V' or LV is None:
            max_abserror = np.max(np.abs((p_fit / pV) - 1) * 100)

        print(max_abserror)
        psat_error = max_abserror

        dropped_indices = [i for i in range(len(n)) if abs(sd[i]) < 1e-15]
        if dropped_indices:
            # for i in reversed(dropped_indices):
            # randomly drop one of them
            n.pop(random.choice(dropped_indices))
            continue

        if max_abserror < 0.5:  # Max error is 0.5%
            Ncoeffs = str(list(myoutput.beta)).lstrip('[').rstrip(']')
            tcoeffs = str(n).lstrip('[').rstrip(']')
            maxerror = max_abserror
            N = len(n)
        else:
            break

        # Find the least significant entry (the one with the largest standard error)
        # and remove it
        n.pop(np.argmax(sd))

        # Remove elements that are not
    import textwrap
    template = textwrap.dedent(
    """
    double {name:s}Class::psat{LV:s}(double T)
    {{
        // Maximum absolute error is {psat_error:f} % between {Tmin:f} K and {Tmax:f} K
        const double t[]={{0, {tcoeffs:s}}};
        const double N[]={{0, {Ncoeffs:s}}};
        double summer=0,theta;
        theta=1-T/reduce.T;
        for (int i=1;i<={N:d};i++)
        {{
            summer += N[i]*pow(theta,t[i]/2);
        }}
        return reduce.p.Pa*exp(reduce.T/T*summer);
    }}
    """)
    the_string = template.format(N=len(n) + 1,
                            tcoeffs=tcoeffs,
                            Ncoeffs=Ncoeffs,
                            name=ClassName,
                            Tmin=Tmin,
                            Tmax=TT[-1],
                            psat_error=maxerror,
                            LV=LV if LV in ['L', 'V'] else ''
                            )

    f = open('anc.txt', 'a')
    f.write(the_string)
    f.close()
    return the_string


if __name__ == '__main__':
    for RPFluid, Fluid in [('REFPROP-MIX:R32[0.47319469]&R125[0.2051091]&R134a[0.32169621]', 'R407F'),
    # for RPFluid,Fluid in [('R11','R11'),
                        ]:
        # saturation_pressure_brute(RPFluid, Fluid
        # saturation_pressure(RPFluid, Fluid, LV = 'L')
        # saturation_pressure(RPFluid, Fluid, LV = 'V')
        saturation_density(RPFluid, Fluid, form='A', LV='L')
        saturation_density(RPFluid, Fluid, form='B', LV='V', perc_error_allowed=0.4)
