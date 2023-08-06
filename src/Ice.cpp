
#ifndef __powerpc__
#    include <math.h>
#    include <complex>
#    include <iostream>
static std::complex<double> t1(0.368017112855051e-1, 0.510878114959572e-1);
static std::complex<double> r1(0.447050716285388e2, 0.656876847463481e2);
static std::complex<double> t2(0.337315741065416, 0.335449415919309);
static std::complex<double> r20(-0.725974574329220e2, -0.781008427112870e2);
static std::complex<double> r21(-0.557107698030123e-4, 0.464578634580806e-4);
static std::complex<double> r22(0.234801409215913e-10, -0.285651142904972e-10);
#endif

#include "Ice.h"

static double T_t = 273.16,  ///< Triple point temperature in K
  p_t = 611.657,             ///< Triple point pressure in Pa
  p_0 = 101325;              ///< Ambient pressure in Pa

// Complex Constants for EOS
static double g00 = -0.632020233449497e6;
static double g01 = 0.655022213658955;
static double g02 = -0.189369929326131e-7;
static double g03 = 0.339746123271053e-14;
static double g04 = -0.556464869058991e-21;
static double s0 = -0.332733756492168e4;

double IsothermCompress_Ice(double T, double p) {
#ifndef __powerpc__
    // Inputs in K, Pa, Output in 1/Pa
    return -dg2_dp2_Ice(T, p) / dg_dp_Ice(T, p);
#else
    return 1e99;
#endif
}
double psub_Ice(double T) {
#ifndef __powerpc__
    double a[] = {0, -0.212144006e2, 0.273203819e2, -0.610598130e1};
    double b[] = {0, 0.333333333e-2, 0.120666667e1, 0.170333333e1};
    double summer = 0, theta;
    theta = T / T_t;
    for (int i = 1; i <= 3; i++) {
        summer += a[i] * pow(theta, b[i]);
    }
    return p_t * exp(1 / theta * summer);
#else
    return 1e99;
#endif
}

double g_Ice(double T, double p) {
#ifndef __powerpc__
    std::complex<double> r2, term1, term2;
    double g0, theta, pi, pi_0;
    theta = T / T_t;
    pi = p / p_t;
    pi_0 = p_0 / p_t;
    g0 = g00 * pow(pi - pi_0, 0.0) + g01 * pow(pi - pi_0, 1.0) + g02 * pow(pi - pi_0, 2.0) + g03 * pow(pi - pi_0, 3.0) + g04 * pow(pi - pi_0, 4.0);
    r2 = r20 * pow(pi - pi_0, 0.0) + r21 * pow(pi - pi_0, 1.0) + r22 * pow(pi - pi_0, 2.0);
    // The two terms of the summation
    term1 = r1 * ((t1 - theta) * log(t1 - theta) + (t1 + theta) * log(t1 + theta) - 2.0 * t1 * log(t1) - theta * theta / t1);
    term2 = r2 * ((t2 - theta) * log(t2 - theta) + (t2 + theta) * log(t2 + theta) - 2.0 * t2 * log(t2) - theta * theta / t2);
    return g0 - s0 * T_t * theta + T_t * real(term1 + term2);
#else
    return 1e99;
#endif
}

double dg_dp_Ice(double T, double p) {
#ifndef __powerpc__
    std::complex<double> r2_p;
    double g0_p, theta, pi, pi_0;
    theta = T / T_t;
    pi = p / p_t;
    pi_0 = p_0 / p_t;
    g0_p = g01 * 1.0 / p_t * pow(pi - pi_0, 1 - 1.0) + g02 * 2.0 / p_t * pow(pi - pi_0, 2 - 1.0) + g03 * 3.0 / p_t * pow(pi - pi_0, 3 - 1.0)
           + g04 * 4.0 / p_t * pow(pi - pi_0, 4 - 1.0);
    r2_p = r21 * 1.0 / p_t * pow(pi - pi_0, 1 - 1.0) + r22 * 2.0 / p_t * pow(pi - pi_0, 2 - 1.0);
    return g0_p + T_t * real(r2_p * ((t2 - theta) * log(t2 - theta) + (t2 + theta) * log(t2 + theta) - 2.0 * t2 * log(t2) - theta * theta / t2));
#else
    return 1e99;
#endif
}

double dg2_dp2_Ice(double T, double p) {
#ifndef __powerpc__
    std::complex<double> r2_pp;
    double g0_pp, theta, pi, pi_0;
    theta = T / T_t;
    pi = p / p_t;
    pi_0 = p_0 / p_t;
    g0_pp = g02 * 2.0 * (2.0 - 1.0) / p_t / p_t * pow(pi - pi_0, 2.0 - 2.0) + g03 * 3.0 * (3.0 - 1.0) / p_t / p_t * pow(pi - pi_0, 3.0 - 2.0)
            + g04 * 4.0 * (4.0 - 1.0) / p_t / p_t * pow(pi - pi_0, 4 - 2.0);
    r2_pp = r22 * 2.0 / p_t / p_t;
    return g0_pp + T_t * real(r2_pp * ((t2 - theta) * log(t2 - theta) + (t2 + theta) * log(t2 + theta) - 2.0 * t2 * log(t2) - theta * theta / t2));
#else
    return 1e99;
#endif
}

double dg_dT_Ice(double T, double p) {
#ifndef __powerpc__
    std::complex<double> r2, term1, term2;
    double theta, pi, pi_0;
    theta = T / T_t;
    pi = p / p_t;
    pi_0 = p_0 / p_t;
    r2 = r20 * pow(pi - pi_0, 0.0) + r21 * pow(pi - pi_0, 1.0) + r22 * pow(pi - pi_0, 2.0);
    // The two terms of the summation
    term1 = r1 * (-log(t1 - theta) + log(t1 + theta) - 2.0 * theta / t1);
    term2 = r2 * (-log(t2 - theta) + log(t2 + theta) - 2.0 * theta / t2);
    return -s0 + real(term1 + term2);
#else
    return 1e99;
#endif
}

double h_Ice(double T, double p) {
#ifndef __powerpc__
    // Returned value is in units of J/kg
    return g_Ice(T, p) - T * dg_dT_Ice(T, p);
#else
    return 1e99;
#endif
}

double rho_Ice(double T, double p) {
#ifndef __powerpc__
    // Returned value is in units of kg/m3
    return 1 / dg_dp_Ice(T, p);
#else
    return 1e99;
#endif
}

double s_Ice(double T, double p) {
#ifndef __powerpc__
    // Returned value is in units of J/kg/K
    return -dg_dT_Ice(T, p);
#else
    return 1e99;
#endif
}
