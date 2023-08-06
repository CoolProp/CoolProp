#include <vector>
#include <iostream>
class EOSFitter;
#include "DataTypes.h"
#include "Fitter.h"

double NonlinearExperimentalDataPoint::residual(const std::vector<double>& n) {
    double summer = a_0(n);
    for (unsigned int i = 1; i < n.size(); i++) {
        summer -= n[i] * a_i(i);
    }
    return summer;
}

double LinearExperimentalDataPoint::residual(const std::vector<double>& n) {
    double summer = a_0();
    for (unsigned int i = 1; i < n.size(); i++) {
        summer -= n[i] * a_i(i);
    }
    return summer;
}

PressureDataPoint::PressureDataPoint(EOSFitter* EOS, double T, double rho, double p, double variance) {
    this->EOS = EOS;
    this->T = T;
    this->rho = rho;
    this->p = p;
    this->tau = EOS->Tr / this->T;
    this->delta = this->rho / EOS->rhor;
    this->log_tau = log(tau);
    this->log_delta = log(delta);
    this->variance = variance;
}
/// The part that does not depend on the coefficients
double PressureDataPoint::a_0() {
    double rhoRT = this->rho * EOS->R * this->T;
    return this->p / rhoRT - 1;
}
/// The part that multiplies the coefficients
double PressureDataPoint::a_i(int i) {
    return delta * EOS->dA_dDelta(log_tau, log_delta, delta, i);
}

SpecificHeatCPDataPoint::SpecificHeatCPDataPoint(EOSFitter* EOS, double T, double rho, double cp, double variance) {
    this->EOS = EOS;
    this->T = T;
    this->rho = rho;
    this->cp = cp;
    this->cp_over_R = this->cp / EOS->R;
    this->tau = EOS->Tr / this->T;
    this->delta = this->rho / EOS->rhor;
    this->log_tau = log(tau);
    this->log_delta = log(delta);
    this->variance = variance;
}
/// The part that does not depend on the coefficients
/// Here it requires that the coefficients be passed in to calculate the precorrelation factor
double SpecificHeatCPDataPoint::a_0(const std::vector<double>& n) {
    // Only calculate this function once to save on calls
    double _dalpha_ddelta = EOS->dalphar_dDelta(log_tau, log_delta, delta);
    // The precorrelation factor
    double e_cp = (pow(1 + delta * _dalpha_ddelta - delta * tau * EOS->d2alphar_dDelta_dTau(log_tau, log_delta, delta), (int)2)
                   / (1 + 2 * delta * _dalpha_ddelta + delta * delta * EOS->d2alphar_dDelta2(log_tau, log_delta, delta)));
    // The a_0 term
    return cp_over_R + tau * tau * EOS->d2alpha0_dTau2(tau, delta) - e_cp;
};
/// The part that multiplies the coefficients
double SpecificHeatCPDataPoint::a_i(int i) {
    return -pow(tau, (int)2) * EOS->d2A_dTau2(log_tau, log_delta, delta, i);
}