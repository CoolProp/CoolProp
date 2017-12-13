# A template the for the .h file
PPF_h_template = """
#ifndef {RefUpper:s}_H
#define {RefUpper:s}_H

    class {Ref:s}Class : public Fluid{{

    public:
        {Ref:s}Class();
        ~{Ref:s}Class(){{}};
        double psatL(double);
        double psatV(double);
        double rhosatL(double);
        double rhosatV(double);
    }};
#endif
"""

PPF_cpp_template = """
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "{Ref:s}.h"

{Ref:s}Class::{Ref:s}Class()
{{
    // Constants for the ideal-gas contribution
    static double a[]={{{acoeffs:s}}};
    static double b[]={{{bcoeffs:s}}};

    // Constants for the residual contribution
    static double N[]={{{Ncoeffs:s}}};
    static double t[]={{{tcoeffs:s}}};
    static double d[]={{{dcoeffs:s}}};
    static double l[]={{{Lcoeffs:s}}};

    // Other fluid parameters
    params.molemass = {molemass:g}; //[kg/kmol]
    params.Ttriple = {Ttriple:g}; //[K]
    params.accentricfactor = {accentric:g}; //[-]
    params.R_u = 8.314472;
    isPure = false;

    // Critical parameters
    crit.rho = {rhocrit:g};
    crit.p = PressureUnit({pcrit:g},UNIT_KPA);
    crit.T = {Tcrit:g};
    crit.v = 1.0/crit.rho;

    phirlist.push_back(new phir_power(N,d,t,l,1,{N_phir:d}-1,{N_phir:d}));

    phi0list.push_back(new phi0_lead(0, 0));
    phi0list.push_back(new phi0_logtau(-1.0));
    phi0list.push_back(new phi0_cp0_poly(a[1],b[1],crit.T,298.15));
    phi0list.push_back(new phi0_Planck_Einstein(a,b,2,{N_cp0:d},{N_cp0:d}+1));

    // Adjust to the IIR reference state (h=200 kJ/kg, s = 1 kJ/kg for sat. liq at 0C)
    params.HSReferenceState = "IIR";

    // Limits of EOS
    limits.Tmin = params.Ttriple;

    name.assign("{Ref:s}");
}}
{pL:s}
{pV:s}
{rhoL:s}
{rhoV:s}
"""
