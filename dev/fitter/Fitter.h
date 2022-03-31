#ifndef FITTER_FITTER_H
#define FITTER_FITTER_H

#include <vector>
#include "DataTypes.h"
#include "Helmholtz.h"

class EOSFitter
{
   public:
    double Tr,  /// The reducing temperature for tau [K]
      rhor,     /// The reducing density for delta [mol/m^3]
      R;        /// The universal gas constant [J/mol/K]
    std::vector<LinearExperimentalDataPoint*> linear_data_points;
    std::vector<NonlinearExperimentalDataPoint*> nonlinear_data_points;
    phir_power alphar;            // Temporary for now
    std::vector<phi_BC*> alpha0;  /// A vector of instances of the phi_BC classes for the ideal-gas Helmholtz energy contribution
    EOSFitter();

    double dA_dDelta(double log_tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double log_delta, double delta, int i);
    double dalphar_dDelta(double log_tau, double log_delta, double delta);
    double d2alphar_dDelta2(double log_tau, double log_delta, double delta);
    double d2alphar_dTau2(double log_tau, double log_delta, double delta);
    double d2alpha0_dTau2(double tau, double delta);
    double d2alphar_dDelta_dTau(double log_tau, double log_delta, double delta);
    /// Set the coefficients in the EOS
    void set_n(const std::vector<double>& n);

    void solve_for_n(std::vector<double>& n, bool non_linear_terms_enabled);
    double sum_squares(std::vector<double>& n, bool non_linear_terms_enabled);
};

class EOSFitterFixedForm : public EOSFitter
{
   public:
    EOSFitterFixedForm(double Tr, double rhor, double R) {
        this->Tr = Tr;
        this->rhor = rhor;
        this->R = R;
    };
};

#endif