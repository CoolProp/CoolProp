//
//  VTPRCubic.h
//  CoolProp
//
//  Created by Ian on 7/17/16.
//
//

#include "GeneralizedCubic.h"
#include "Exceptions.h"

#ifndef VTPRCubic_h
#    define VTPRCubic_h

class VTPRCubic : public PengRobinson
{
   private:
    UNIFAC::UNIFACMixture unifaq;

   public:
    VTPRCubic(std::vector<double> Tc, std::vector<double> pc, std::vector<double> acentric, double R_u,
              const UNIFACLibrary::UNIFACParameterLibrary& lib)
      : PengRobinson(Tc, pc, acentric, R_u), unifaq(lib, T_r){};

    VTPRCubic(double Tc, double pc, double acentric, double R_u, const UNIFACLibrary::UNIFACParameterLibrary& lib)
      : PengRobinson(std::vector<double>(1, Tc), std::vector<double>(1, pc), std::vector<double>(1, acentric), R_u), unifaq(lib, T_r){};

    /// Get a reference to the managed UNIFAC instance
    UNIFAC::UNIFACMixture& get_unifaq() {
        return unifaq;
    }

    /// Calculate the non-dimensionalized gE/RT term
    double gE_R_RT(double tau, const std::vector<double>& x, std::size_t itau) {
        double summer = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            summer += x[i] * unifaq.ln_gamma_R(tau, i, itau);
        }
        return summer;
    }
    double d_gE_R_RT_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
        if (xN_independent) {
            return unifaq.ln_gamma_R(tau, i, itau);
        } else {
            return unifaq.ln_gamma_R(tau, i, itau) - unifaq.ln_gamma_R(tau, N - 1, itau);
        }
    }
    double gE_R(double tau, const std::vector<double>& x, std::size_t itau) {
        if (x.size() == 1) {
            return 0.;
        } else {
            switch (itau) {
                case 0: {
                    return R_u * T_r / tau * gE_R_RT(tau, x, 0);
                }
                case 1: {
                    return R_u * T_r / tau * (-gE_R_RT(tau, x, 0) / tau + gE_R_RT(tau, x, 1));
                }
                case 2: {
                    return R_u * T_r / tau * (2 * gE_R_RT(tau, x, 0) / powInt(tau, 2) - 2 * gE_R_RT(tau, x, 1) / tau + gE_R_RT(tau, x, 2));
                }
                case 3: {
                    return R_u * T_r / tau
                           * (-6 * gE_R_RT(tau, x, 0) / powInt(tau, 3) + 6 * gE_R_RT(tau, x, 1) / powInt(tau, 2) - 3 * gE_R_RT(tau, x, 2) / tau
                              + gE_R_RT(tau, x, 3));
                }
                case 4: {
                    return R_u * T_r / tau
                           * (24 * gE_R_RT(tau, x, 0) / powInt(tau, 4) - 24 * gE_R_RT(tau, x, 1) / powInt(tau, 3)
                              + 12 * gE_R_RT(tau, x, 2) / powInt(tau, 2) - 4 * gE_R_RT(tau, x, 3) / tau + gE_R_RT(tau, x, 4));
                }
                default:
                    throw CoolProp::ValueError(format("itau (%d) is invalid", itau));
            }
        }
    }
    double d_gE_R_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
        if (x.size() == 1) {
            return 0.;
        } else {
            switch (itau) {
                case 0: {
                    return R_u * T_r / tau * d_gE_R_RT_dxi(tau, x, 0, i, xN_independent);
                }
                case 1: {
                    return R_u * T_r / tau * (-d_gE_R_RT_dxi(tau, x, 0, i, xN_independent) / tau + d_gE_R_RT_dxi(tau, x, 1, i, xN_independent));
                }
                case 2: {
                    return R_u * T_r / tau
                           * (2 * d_gE_R_RT_dxi(tau, x, 0, i, xN_independent) / powInt(tau, 2) - 2 * d_gE_R_RT_dxi(tau, x, 1, i, xN_independent) / tau
                              + d_gE_R_RT_dxi(tau, x, 2, i, xN_independent));
                }
                case 3: {
                    return R_u * T_r / tau
                           * (-6 * d_gE_R_RT_dxi(tau, x, 0, i, xN_independent) / powInt(tau, 3)
                              + 6 * d_gE_R_RT_dxi(tau, x, 1, i, xN_independent) / powInt(tau, 2)
                              - 3 * d_gE_R_RT_dxi(tau, x, 2, i, xN_independent) / tau + d_gE_R_RT_dxi(tau, x, 3, i, xN_independent));
                }
                case 4: {
                    return R_u * T_r / tau
                           * (24 * d_gE_R_RT_dxi(tau, x, 0, i, xN_independent) / powInt(tau, 4)
                              - 24 * d_gE_R_RT_dxi(tau, x, 1, i, xN_independent) / powInt(tau, 3)
                              + 12 * d_gE_R_RT_dxi(tau, x, 2, i, xN_independent) / powInt(tau, 2)
                              - 4 * d_gE_R_RT_dxi(tau, x, 3, i, xN_independent) / tau + d_gE_R_RT_dxi(tau, x, 4, i, xN_independent));
                }
                default:
                    throw CoolProp::ValueError(format("itau (%d) is invalid", itau));
            }
        }
    }
    double am_term(double tau, const std::vector<double>& x, std::size_t itau) {
        return bm_term(x) * (sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087));
    }
    double sum_xi_aii_bii(double tau, const std::vector<double>& x, std::size_t itau) {
        double summeram = 0;
        for (int i = 0; i < N; ++i) {
            summeram += x[i] * aii_term(tau, i, itau) / b0_ii(i);
        }
        return summeram;
    }
    double d_sum_xi_aii_bii_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
        if (xN_independent) {
            return aii_term(tau, i, itau) / b0_ii(i);
        } else {
            return aii_term(tau, i, itau) / b0_ii(i) - aii_term(tau, N - 1, itau) / b0_ii(N - 1);
        }
    }
    double d_am_term_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
        return d_bm_term_dxi(x, i, xN_independent) * (sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087))
               + bm_term(x) * (d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087));
    }
    double d2_am_term_dxidxj(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent) {
        return d2_bm_term_dxidxj(x, i, j, xN_independent) * (sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087))
               + d_bm_term_dxi(x, i, xN_independent)
                   * (d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087))
               + d_bm_term_dxi(x, j, xN_independent)
                   * (d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087));
    }
    double d3_am_term_dxidxjdxk(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k,
                                bool xN_independent) {
        return d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent) * (sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087))
               + d2_bm_term_dxidxj(x, i, k, xN_independent)
                   * (d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087))
               + d2_bm_term_dxidxj(x, j, k, xN_independent)
                   * (d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087));
    }

    double bm_term(const std::vector<double>& x) {
        double summerbm = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                summerbm += x[i] * x[j] * bij_term(i, j);
            }
        }
        return summerbm;
    }
    double bij_term(std::size_t i, std::size_t j) {
        return pow((pow(b0_ii(i), 0.75) + pow(b0_ii(j), 0.75)) / 2.0, 4.0 / 3.0);
    }
    double d_bm_term_dxi(const std::vector<double>& x, std::size_t i, bool xN_independent) {
        double summer = 0;
        if (xN_independent) {
            for (int j = N - 1; j >= 0; --j) {
                summer += x[j] * bij_term(i, j);
            }
            return 2 * summer;
        } else {
            for (int k = N - 2; k >= 0; --k) {
                summer += x[k] * (bij_term(i, k) - bij_term(k, N - 1));
            }
            return 2 * (summer + x[N - 1] * (bij_term(N - 1, i) - bij_term(N - 1, N - 1)));
        }
    }
    double d2_bm_term_dxidxj(const std::vector<double>& x, std::size_t i, std::size_t j, bool xN_independent) {
        if (xN_independent) {
            return 2 * bij_term(i, j);
        } else {
            return 2 * (bij_term(i, j) - bij_term(j, N - 1) - bij_term(N - 1, i) + bij_term(N - 1, N - 1));
        }
    }
    double d3_bm_term_dxidxjdxk(const std::vector<double>& x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent) {
        return 0;
    }
    // Allows to modify the unifac interaction parameters aij, bij and cij
    void set_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter, const double value) {
        unifaq.set_interaction_parameter(mgi1, mgi2, parameter, value);
    }

    // Allows to modify the surface parameter Q_k of the sub group sgi
    void set_Q_k(const size_t sgi, const double value) {
        unifaq.set_Q_k(sgi, value);
    }

    // Get the surface parameter Q_k of the sub group sgi
    double get_Q_k(const size_t sgi) const {
        return unifaq.get_Q_k(sgi);
    }

    // Allows to modify the unifac interaction parameters aij, bij and cij
    double get_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter) {
        return unifaq.get_interaction_parameter(mgi1, mgi2, parameter);
    }
    std::vector<double> ln_fugacity_coefficient(const std::vector<double>& z, double rhomolar, double p, double T) {
        double v = 1 / rhomolar;
        // Common terms for all components
        double tau = get_Tr() / T;
        double b = bm_term(z);
        double c = cm_term();
        double R = get_R_u();
        std::vector<double> ln_phi;
        double bracket = log((v + c + (1 + sqrt(2.0)) * b) / (v + c + (1 - sqrt(2.0)) * b));
        for (std::size_t i = 0; i < z.size(); ++i) {
            double summer1 = 0;
            for (std::size_t j = 0; j < z.size(); ++j) {
                summer1 += z[j] * bij_term(i, j);
            }
            double a_i_over_b_i = aii_term(tau, i, 0) / b0_ii(i);
            double c_i = 0;  // TODO: fix this, allow for volume translation
            double _ln_phi = (2 / b * summer1 - 1) * (p * (v + c) / (R * T) - 1) - p * c_i / (R * T) - log(p * (v + c - b) / (R * T))
                             - 1.0 / (2.0 * sqrt(2.0) * R * T) * (a_i_over_b_i + R * T * unifaq.ln_gamma_R(tau, i, 0) / -0.53087) * bracket;
            ln_phi.push_back(_ln_phi);
        }
        return ln_phi;
    }
};

#endif /* VTPRCubic_h */
