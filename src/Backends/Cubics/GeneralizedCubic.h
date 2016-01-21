#ifndef CUBIC_H
#define CUBIC_H

#include <vector>
#include <cmath>

class AbstractCubic
{
protected:
	std::vector<double> Tc, ///< Vector of critical temperatures (in K)
		                pc, ///< Vector of critical pressures (in Pa) 
						acentric; ///< Vector of acentric factors (unitless)
	double R_u; ///< The universal gas constant  in J/(mol*K)
	double Delta_1, ///< The first cubic constant
		   Delta_2; ///< The second cubic constant
	int N; ///< Number of components in the mixture
	std::vector< std::vector<double> > k;///< The interaction parameters (k_ii = 0)
public:
	static const double rho_r, T_r;
	/**
	 \brief 
	 */
	AbstractCubic(std::vector<double> Tc, 
		          std::vector<double> pc, 
				  std::vector<double> acentric,
				  double R_u,
				  double Delta_1,
				  double Delta_2) 
		: Tc(Tc), pc(pc), acentric(acentric), R_u(R_u), Delta_1(Delta_1), Delta_2(Delta_2) 
		{
			N = static_cast<int>(Tc.size());
			k.resize(N, std::vector<double>(N, 0));
		};
	
	/// Get the leading constant in the expression for the pure fluid attractive energy term 
	/// (must be implemented by derived classes)
	virtual double a0_ii(std::size_t i) = 0;
	/// Get the leading constant in the expression for the pure fluid covolume term 
	/// (must be implemented by derived classes)
	virtual double b0_ii(std::size_t i) = 0;
	
	virtual double m_ii(std::size_t i) = 0;
	
	///
	double alphar(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta);
	double d_alphar_dxi(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
	double d2_alphar_dxidxj(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
	
	/** Take the n-th derivative of \f$a_m\f$
	 * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
	 * \param x The vector of mole fractions
	 * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
	 */
	double am_term(double tau, const std::vector<double> &x, std::size_t itau);
	double d_am_term_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent);
	double d2_am_term_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);


	/** The term \f$b_m\f$
	 * \param x The vector of mole fractions
	 */
	double bm_term(const std::vector<double> &x);
	double d_bm_term_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent);
	double d2_bm_term_dxidxj(const std::vector<double> &x, std::size_t i, std::size_t j, bool xN_independent);

	/** Take the n-th tau derivative of \f$a_{ij}(\tau)\f$
	 * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
	 * \param i The first index
	 * \param j The second index
	 * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
	 */
	double aij_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
	/** Take the n-th tau derivative of \f$u(\tau)\f$, the argument of sqrt in the cross aij term
	 * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
	 * \param i The first index
	 * \param j The first index
	 * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
	 */
	double u_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
	/** Take the n-th tau derivative of the \f$a_{ii}(\tau)\f$ pure fluid contribution
	 * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
	 * \param i The index of the component
	 * \param itau The number of derivatives of \f$u\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
	 */
	double aii_term(double tau, std::size_t i, std::size_t itau);

	double psi_minus(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta);
	double d_psi_minus_dxi(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
	double d2_psi_minus_dxidxj(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);

	double PI_12(double delta, const std::vector<double> &x, std::size_t idelta);
	double d_PI_12_dxi(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, bool xN_independent);
	double d2_PI_12_dxidxj(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);

	double tau_times_a(double tau, const std::vector<double> &x, std::size_t itau);
	double d_tau_times_a_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent);
	double d2_tau_times_a_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);

	double psi_plus(double delta, const std::vector<double> &x, std::size_t idelta);
	double d_psi_plus_dxi(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, bool xN_independent);
	double d2_psi_plus_dxidxj(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
};

class PengRobinson : public AbstractCubic
{
public:
	PengRobinson(std::vector<double> Tc, 
		         std::vector<double> pc, 
		         std::vector<double> acentric,
		         double R_u) 
		: AbstractCubic(Tc, pc, acentric, R_u, 1+sqrt(2), 1-sqrt(2)) {};
	PengRobinson(double Tc, 
		double pc, 
		double acentric,
		double R_u) 
		: AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1+sqrt(2), 1-sqrt(2)) {};

	double a0_ii(std::size_t i);
	double b0_ii(std::size_t i);
	double m_ii(std::size_t i);
};

class SRK : public AbstractCubic
{
public:
	SRK(std::vector<double> Tc, 
		std::vector<double> pc, 
		std::vector<double> acentric,
		double R_u) 
		: AbstractCubic(Tc, pc, acentric, R_u, 1, 0) {};
	SRK(double Tc, 
		double pc, 
		double acentric,
		double R_u) 
		: AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1, 0) {};

	double a0_ii(std::size_t i);
	double b0_ii(std::size_t i);
	double m_ii(std::size_t i);
};


#endif