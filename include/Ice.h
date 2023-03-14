#ifndef ICE_H
#define ICE_H

double psub_Ice(double T);
double g_Ice(double T, double p);
double dg_dp_Ice(double T, double p);
double dg2_dp2_Ice(double T, double p);
double IsothermCompress_Ice(double T, double p);
double dg_dT_Ice(double T, double p);
double h_Ice(double T, double p);
double s_Ice(double T, double p);
double rho_Ice(double T, double p);

#endif