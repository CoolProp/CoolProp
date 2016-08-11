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
    : PengRobinson(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u), unifaq(lib) {};
    
    /// Get a reference to the managed UNIFAQ instance
    UNIFAQ::UNIFAQMixture &get_unifaq(){ return unifaq; }
    
    /// The attractive part in cubic EOS
    double a_alpha(double T, std::size_t i){
        return pow(1 + m_ii(i)*(1-sqrt(T/Tc[i])), 2);
    }
    /// Calculate the non-dimensionalized gE/RT term
    double gE_R_RT(){
        const std::vector<double> &z = unifaq.get_mole_fractions();
        double summer = 0;
        for (std::size_t i = 0; i < z.size(); ++i) {
            summer += z[i]*unifaq.ln_gamma_R(i);
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
    /// The co-volume for the i-th pure component
    double b_ii(std::size_t i){
        return 0.0778*R_u*Tc[i]/pc[i];
    }
    /// The attractive parameter for the i-th pure component
    double a_ii(std::size_t i){
        double a0 = 0.45724*pow(R_u*Tc[i], 2)/pc[i];
        double alpha = a_alpha(unifaq.get_temperature(), i);
        return a0*alpha;
    }
    double am_term(double tau, const std::vector<double> &x, std::size_t itau){
        if (itau == 0){
            set_temperature(T_r/tau);
			return bm_term(x)*(sum_xi_aii_bii(x) + R_u*unifaq.get_temperature()*gE_R_RT() / (-0.53087));
        }
        else{
			double dtau = 0.01*tau;
			return (am_term(tau + dtau, x, itau - 1) - am_term(tau - dtau, x, itau - 1)) / (2 * dtau);
        }
    }
	double sum_xi_aii_bii(const std::vector<double> &x) {
		double summeram = 0;
		for (std::size_t i = 0; i < N; ++i) {
			summeram += x[i] * a_ii(i) / b_ii(i);
		}
		return summeram;
	}
	double d_sum_xi_aii_bii_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent) {
		if (xN_independent)
		{
			return x[i] * a_ii(i) / b_ii(i);
		}
		else {
			return x[i] * a_ii(i) / b_ii(i) - x[N - 1] * a_ii(N - 1) / b_ii(N - 1);
		}
	}
	double d_am_term_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent)
	{
		set_temperature(T_r / tau);
		return d_bm_term_dxi(x, i, xN_independent)*(sum_xi_aii_bii(x) + R_u*unifaq.get_temperature()*gE_R_RT() / (-0.53087))
			+ bm_term(x)*(d_sum_xi_aii_bii_dxi(x, i, xN_independent) + R_u*unifaq.get_temperature()*d_gE_R_RT_dxi(x, i, xN_independent) / (-0.53087));
	}
	double d2_am_term_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent)
	{
		set_temperature(T_r / tau);
		return d2_bm_term_dxidxj(x, i, j, xN_independent)*(sum_xi_aii_bii(x) + R_u*unifaq.get_temperature()*gE_R_RT() / (-0.53087))
			+ d_bm_term_dxi(x, i, xN_independent)*(d_sum_xi_aii_bii_dxi(x, j, xN_independent) + R_u*unifaq.get_temperature()*d_gE_R_RT_dxi(x, j, xN_independent) / (-0.53087))
			+ d_bm_term_dxi(x, j, xN_independent)*(d_sum_xi_aii_bii_dxi(x, i, xN_independent) + R_u*unifaq.get_temperature()*d_gE_R_RT_dxi(x, i, xN_independent) / (-0.53087));
	}
	double d3_am_term_dxidxjdxk(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
	{
		set_temperature(T_r / tau);
		return d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent)*(sum_xi_aii_bii(x) + R_u*unifaq.get_temperature()*gE_R_RT() / (-0.53087))
			+ d2_bm_term_dxidxj(x, i, k, xN_independent)*(d_sum_xi_aii_bii_dxi(x, j, xN_independent) + R_u*unifaq.get_temperature()*d_gE_R_RT_dxi(x, j, xN_independent) / (-0.53087))
			+ d2_bm_term_dxidxj(x, j, k, xN_independent)*(d_sum_xi_aii_bii_dxi(x, i, xN_independent) + R_u*unifaq.get_temperature()*d_gE_R_RT_dxi(x, i, xN_independent) / (-0.53087));
	}

    double bm_term(const std::vector<double> &x){
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
		return pow((pow(b_ii(i), 0.75) + pow(b_ii(j), 0.75)) / 2.0, 4.0 / 3.0);
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

    void set_temperature(const double T){ unifaq.set_temperature(T); }
};

#endif /* VTPRCubic_h */
