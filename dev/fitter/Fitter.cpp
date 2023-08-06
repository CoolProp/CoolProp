#include <iostream>
#include "Eigen/Dense"
#include "Helmholtz.h"

class EOSFitter;
#include "DataTypes.h"
#include "Fitter.h"

EOSFitter::EOSFitter() {
    this->Tr = Tr;
    this->rhor = rhor;
    this->R = R;
};
double EOSFitter::dA_dDelta(double log_tau, double log_delta, double delta, int i) {
    return alphar.dA_dDelta(log_tau, log_delta, delta, i);
};
double EOSFitter::d2A_dTau2(double log_tau, double log_delta, double delta, int i) {
    return alphar.d2A_dTau2(log_tau, log_delta, delta, i);
};
double EOSFitter::dalphar_dDelta(double log_tau, double log_delta, double delta) {
    double summer = 0;
    for (unsigned int i = 0; i < alphar.n.size(); i++) {
        summer += alphar.n[i] * alphar.dA_dDelta(log_tau, log_delta, delta, i);
    }
    return summer;
};
double EOSFitter::d2alphar_dDelta2(double log_tau, double log_delta, double delta) {
    double summer = 0;
    for (unsigned int i = 0; i < alphar.n.size(); i++) {
        summer += alphar.n[i] * alphar.d2A_dDelta2(log_tau, log_delta, delta, i);
    }
    return summer;
};
double EOSFitter::d2alphar_dTau2(double log_tau, double log_delta, double delta) {
    double summer = 0;
    for (unsigned int i = 0; i < alphar.n.size(); i++) {
        summer += alphar.n[i] * alphar.d2A_dTau2(log_tau, log_delta, delta, i);
    }
    return summer;
};
double EOSFitter::d2alpha0_dTau2(double tau, double delta) {
    double summer = 0;
    for (std::vector<phi_BC*>::iterator it = alpha0.begin(); it != alpha0.end(); it++) {
        summer += (*it)->dTau2(tau, delta);
    }
    return summer;
};
double EOSFitter::d2alphar_dDelta_dTau(double log_tau, double log_delta, double delta) {
    double summer = 0;
    for (unsigned int i = 0; i < alphar.n.size(); i++) {
        summer += alphar.n[i] * alphar.d2A_dDelta_dTau(log_tau, log_delta, delta, i);
    }
    return summer;
};
/// Set the coefficients in the EOS
void EOSFitter::set_n(const std::vector<double>& n) {
    alphar.n = n;
};

void EOSFitter::solve_for_n(std::vector<double>& n, bool non_linear_terms_enabled) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(21, 21);
    Eigen::VectorXd Q = Eigen::VectorXd::Random(21);

    // Build the A matrix and the Q vector
    for (int i = 1; i <= A.rows(); i++) {
        // The i-th row of the A matrix (Span 2000 Eq. 4.9)
        for (int j = 1; j <= A.cols(); j++) {
            // The entry for the j-th column and i-th row
            double summer = 0;
            for (unsigned int m = 0; m < linear_data_points.size(); m++) {
                LinearExperimentalDataPoint& pt = *linear_data_points[m];
                summer += (pt.a_i(i) * pt.a_i(j)) / pow(pt.variance, (int)2);
            }
            if (non_linear_terms_enabled) {
                for (unsigned int m = 0; m < nonlinear_data_points.size(); m++) {
                    NonlinearExperimentalDataPoint& pt = *nonlinear_data_points[m];
                    summer += (pt.a_i(i) * pt.a_i(j)) / pow(pt.variance, (int)2);
                }
            }
            A(i - 1, j - 1) = summer;
        }
        // The i-th entry in the Q column vector
        double summer = 0;
        for (unsigned int m = 0; m < linear_data_points.size(); m++) {
            LinearExperimentalDataPoint& pt = *linear_data_points[m];
            summer += (pt.a_i(i) * pt.a_0()) / pow(pt.variance, (int)2);
        }
        if (non_linear_terms_enabled) {
            for (unsigned int m = 0; m < nonlinear_data_points.size(); m++) {
                NonlinearExperimentalDataPoint& pt = *nonlinear_data_points[m];
                summer += (pt.a_i(i) * pt.a_0(n)) / pow(pt.variance, (int)2);
            }
        }
        Q(i - 1) = summer;
    }

    Eigen::VectorXd N = A.colPivHouseholderQr().solve(Q);
    for (unsigned int i = 0; i < n.size() - 1; i++) {
        n[i + 1] = N(i);
    }

    double relative_error = (A * N - Q).norm() / Q.norm();
};

double EOSFitter::sum_squares(std::vector<double>& n, bool non_linear_terms_enabled) {
    double summer = 0;
    for (unsigned int m = 0; m < linear_data_points.size(); m++) {
        LinearExperimentalDataPoint& pt = *linear_data_points[m];
        summer += pt.sum_squares(n);
    }
    if (non_linear_terms_enabled) {
        for (unsigned int m = 0; m < nonlinear_data_points.size(); m++) {
            NonlinearExperimentalDataPoint& pt = *nonlinear_data_points[m];
            summer += pt.sum_squares(n);
        }
    }
    return summer;
}