//
//  VTPRCubic.h
//  CoolProp
//
//  Created by Ian on 7/17/16.
//
//

#include "GeneralizedCubic.h"

#ifndef VTPRCubic_h
#define VTPRCubic_h

class VTPRCubic : public PengRobinson
{
private:
    UNIFAQ::UNIFAQMixture unifaq;
public:
    VTPRCubic(std::vector<double> Tc,
        std::vector<double> pc,
        std::vector<double> acentric,
        double R_u,
        const UNIFAQLibrary::UNIFAQParameterLibrary & lib
    )
        : PengRobinson(Tc, pc, acentric, R_u), unifaq(lib) {};

    VTPRCubic(double Tc,
        double pc,
        double acentric,
        double R_u,
        const UNIFAQLibrary::UNIFAQParameterLibrary & lib)
        : PengRobinson(std::vector<double>(1, Tc), std::vector<double>(1, pc), std::vector<double>(1, acentric), R_u), unifaq(lib) {};

    /// Get a reference to the managed UNIFAQ instance
    UNIFAQ::UNIFAQMixture &get_unifaq() { return unifaq; }

    /// Calculate the non-dimensionalized gE/RT term
    double gE_R_RT(const std::vector<double> &x) {
        double summer = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            summer += x[i] * unifaq.ln_gamma_R(i);
        }
        return summer;
    }
    double d_gE_R_RT_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent) {
        if (xN_independent)
        {
            return unifaq.ln_gamma_R(i);
        }
        else {
            return unifaq.ln_gamma_R(i) - unifaq.ln_gamma_R(N - 1);
        }
    }
    double gE_R(double tau, const std::vector<double> &x, std::size_t itau) {
        if (x.size() == 1) {
            return 0.;
        }
        else {
            if (itau == 0) {
                set_temperature(T_r / tau, x);
                return R_u*T_r / tau*gE_R_RT(x);
            }
            else {
                double dtau = 0.01*tau;
                return (gE_R(tau + dtau, x, itau - 1) - gE_R(tau - dtau, x, itau - 1)) / (2 * dtau);
            }
        }
    }
    double d_gE_R_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent) {
        if (x.size() == 1) {
            return 0.;
        }
        else {
            if (itau == 0) {
                set_temperature(T_r / tau, x);
                return R_u*T_r / tau*d_gE_R_RT_dxi(x, i, xN_independent);
            }
            else {
                double dtau = 0.01*tau;
                return (d_gE_R_dxi(tau + dtau, x, itau - 1, i, xN_independent) - d_gE_R_dxi(tau - dtau, x, itau - 1, i, xN_independent)) / (2 * dtau);
            }
        }
    }
    double am_term(double tau, const std::vector<double> &x, std::size_t itau) {
        return bm_term(x)*(sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087));
    }
    double sum_xi_aii_bii(double tau, const std::vector<double> &x, std::size_t itau) {
        double summeram = 0;
        for (std::size_t i = 0; i < N; ++i) {
            summeram += x[i] * aii_term(tau, i, itau) / b0_ii(i);
        }
        return summeram;
    }
    double d_sum_xi_aii_bii_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent) {
        if (xN_independent)
        {
            return aii_term(tau, i, itau) / b0_ii(i);
        }
        else {
            return aii_term(tau, i, itau) / b0_ii(i) - aii_term(tau, N - 1, itau) / b0_ii(N - 1);
        }
    }
    double d_am_term_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent)
    {
        return d_bm_term_dxi(x, i, xN_independent)*(sum_xi_aii_bii(tau,x,itau) + gE_R(tau, x, itau) / (-0.53087))
            + bm_term(x)*(d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau , x, itau, i, xN_independent) / (-0.53087));
    }
    double d2_am_term_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent)
    {
        return d2_bm_term_dxidxj(x, i, j, xN_independent)*(sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087))
            + d_bm_term_dxi(x, i, xN_independent)*(d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087))
            + d_bm_term_dxi(x, j, xN_independent)*(d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087));
    }
    double d3_am_term_dxidxjdxk(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
    {
        return d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent)*(sum_xi_aii_bii(tau, x, itau) + gE_R(tau, x, itau) / (-0.53087))
            + d2_bm_term_dxidxj(x, i, k, xN_independent)*(d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087))
            + d2_bm_term_dxidxj(x, j, k, xN_independent)*(d_sum_xi_aii_bii_dxi(tau, x, itau, i, xN_independent) + d_gE_R_dxi(tau, x, itau, i, xN_independent) / (-0.53087));
    }

    double bm_term(const std::vector<double> &x) {
        double summerbm = 0;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                summerbm += x[i] * x[j] * bij_term(i, j);
            }
        }
        return summerbm;
    }
    double bij_term(std::size_t i, std::size_t j)
    {
        return pow((pow(b0_ii(i), 0.75) + pow(b0_ii(j), 0.75)) / 2.0, 4.0 / 3.0);
    }
    double d_bm_term_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent)
    {
        double summer = 0;
        if (xN_independent)
        {
            for (int j = N - 1; j >= 0; --j)
            {
                summer += x[j] * bij_term(i, j);
            }
            return 2 * summer;
        }
        else {
            for (int k = N - 2; k >= 0; --k)
            {
                summer += x[k] * (bij_term(i, k) - bij_term(k, N - 1));
            }
            return 2 * (summer + x[N - 1] * (bij_term(N - 1, i) - bij_term(N - 1, N - 1)));
        }
    }
    double d2_bm_term_dxidxj(const std::vector<double> &x, std::size_t i, std::size_t j, bool xN_independent)
    {
        if (xN_independent)
        {
            return 2 * bij_term(i, j);
        }
        else {
            return 2 * (bij_term(i, j) - bij_term(j, N - 1) - bij_term(N - 1, i) + bij_term(N - 1, N - 1));
        }
    }
    double d3_bm_term_dxidxjdxk(const std::vector<double> &x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
    {
        return 0;
    }

    void set_temperature(const double T, const std::vector<double> &z) { unifaq.set_temperature(T, z); }
};

#endif /* VTPRCubic_h */
