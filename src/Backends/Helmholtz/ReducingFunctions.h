/** \brief Code for all the binary pairs in the mixture
 * 
 * This includes both binary pair information for the reducing functions as well as the departure 
 * functions for the given binary pair.
 */
 
#ifndef MIXTURE_BINARY_PAIRS_H
#define MIXTURE_BINARY_PAIRS_H

#include <vector>
#include "CoolPropFluid.h"
#include "crossplatform_shared_ptr.h"

namespace CoolProp{

typedef std::vector<std::vector<long double> > STLMatrix;

enum x_N_dependency_flag{XN_INDEPENDENT, ///< x_N is an independent variable, and not calculated by \f$ x_N = 1-\sum_i x_i\f$
                         XN_DEPENDENT ///< x_N is an dependent variable, calculated by \f$ x_N = 1-\sum_i x_i\f$
                         };
                 
std::string get_reducing_function_name(std::string CAS1, std::string CAS2);

/** \brief Abstract base class for reducing function
 * An abstract base class for the reducing function to allow for
 * Lemmon-Jacobsen, GERG, or other reducing function to yield the 
 * reducing parameters \f$\rho_r\f$ and \f$T_r\f$
*/
class ReducingFunction
{
protected:
	std::size_t N;
public:
	ReducingFunction(){};
	virtual ~ReducingFunction(){};

    /// A factory function to generate the requiredreducing function
    static shared_ptr<ReducingFunction> factory(const std::vector<CoolPropFluid*> &components, std::vector< std::vector< long double> > &F);

	/// The reduced temperature
	virtual long double Tr(const std::vector<long double> &x) = 0;
	/// The derivative of reduced temperature with respect to component i mole fraction
	virtual long double dTrdxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag) = 0;
	/// The molar reducing density
	virtual long double rhormolar(const std::vector<long double> &x) = 0;
	///Derivative of the molar reducing density with respect to component i mole fraction
	virtual long double drhormolardxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag) = 0;

	virtual long double d2rhormolardxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag) = 0;
	virtual long double d2rhormolardxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) = 0;
	virtual long double d2Trdxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag) = 0;
	virtual long double d2Trdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) = 0;

	/** \brief GERG 2004 Monograph equation 7.56:
     * 
	 * \f[
	 * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2T_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2T_r}{\partial x_j \partial x_k}\right)
	 * \f]
	 */
	long double d_ndTrdni_dxj__constxi(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);
	/** \brief GERG 2004 Monograph equation 7.55:
     * 
	 * \f[
	 * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2\rho_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2\rho_r}{\partial x_j \partial x_k}\right)
	 * \f]
	 */
	long double d_ndrhorbardni_dxj__constxi(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

	long double ndrhorbardni__constnj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	long double ndTrdni__constnj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
};

/** \brief The reducing function model of GERG-2008
 * 
 * Used by the GERG-2008 formulation to yield the
 * reducing parameters \f$ \rho_r \f$ and \f$ T_r \f$ and derivatives thereof
 */
class GERG2008ReducingFunction : public ReducingFunction
{
private:
    GERG2008ReducingFunction(const GERG2008ReducingFunction& that); // No copying
protected:
	STLMatrix v_c; ///< \f$ v_{c,ij} = \frac{1}{8}\left(v_{c,i}^{1/3}+v_{c,j}^{1/3}\right)^{3}\f$ from GERG-2008
	STLMatrix T_c; ///< \f$ T_{c,ij} = \sqrt{T_{c,i}T_{c,j}} \f$ from GERG=2008
	STLMatrix beta_v; ///< \f$ \beta_{v,ij} \f$ from GERG-2008
	STLMatrix gamma_v; ///< \f$ \gamma_{v,ij} \f$ from GERG-2008
	STLMatrix beta_T; ///< \f$ \beta_{T,ij} \f$ from GERG-2008
	STLMatrix gamma_T; ///< \f$ \gamma_{T,ij} \f$ from GERG-2008
    std::vector<long double> Yc_T; ///< Vector of critical temperatures for all components
    std::vector<long double> Yc_v; ///< Vector of critical molar volumes for all components
	std::vector<CoolPropFluid *> pFluids; ///< List of pointer to fluids

public:
	GERG2008ReducingFunction(std::vector<CoolPropFluid *> pFluids, STLMatrix beta_v, STLMatrix gamma_v, STLMatrix beta_T, STLMatrix gamma_T)
	{
		this->pFluids = pFluids;
		this->beta_v = beta_v;
		this->gamma_v = gamma_v;
		this->beta_T = beta_T;
		this->gamma_T = gamma_T;
		this->N = pFluids.size();
		T_c.resize(N,std::vector<long double>(N,0));
		v_c.resize(N,std::vector<long double>(N,0));
        Yc_T.resize(N);
        Yc_v.resize(N);
		for (std::size_t i = 0; i < N; ++i)
		{
			for (std::size_t j = 0; j < N; j++)
			{
				T_c[i][j] = sqrt(pFluids[i]->pEOS->reduce.T*pFluids[j]->pEOS->reduce.T);
				v_c[i][j] = 1.0/8.0*pow(pow(pFluids[i]->pEOS->reduce.rhomolar, -1.0/3.0)+pow(pFluids[j]->pEOS->reduce.rhomolar, -1.0/3.0),3);
			}
            Yc_T[i] = pFluids[i]->pEOS->reduce.T;
            Yc_v[i] = 1/pFluids[i]->pEOS->reduce.rhomolar;
		}
	};

	/// Default destructor
	~GERG2008ReducingFunction(){};
	/// The reduced temperature
	long double Tr(const std::vector<long double> &x);
	/// The derivative of reduced temperature with respect to component i mole fraction
	long double dTrdxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	/// The molar reducing density
	long double rhormolar(const std::vector<long double> &x);
	///Derivative of the molar reducing density with respect to component i mole fraction
	long double drhormolardxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	long double dvrmolardxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag  xN_fla);

	long double d2vrmolardxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	long double d2rhormolardxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	long double d2vrmolardxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);
	long double d2rhormolardxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);
	long double d2Trdxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag);
	long double d2Trdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    long double Yr(const std::vector<long double> &x, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc);
    long double dYrdxi__constxj(const std::vector<long double> &x, std::size_t i, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag);
    long double d2Yrdxi2__constxj(const std::vector<long double> &x, std::size_t i, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag);
    long double d2Yrdxidxj__constxj(const std::vector<long double> &x, std::size_t i, std::size_t j, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag);
	long double c_Y_ij(std::size_t i, std::size_t j, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c);
	long double f_Y_ij(const std::vector<long double> &x, std::size_t i, std::size_t j, const STLMatrix &beta);

	long double dfYkidxi__constxk(const std::vector<long double> &x, std::size_t k, std::size_t i, const STLMatrix &beta);
	long double dfYikdxi__constxk(const std::vector<long double> &x, std::size_t i, std::size_t k, const STLMatrix &beta);
	long double d2fYkidxi2__constxk(const std::vector<long double> &x, std::size_t k, std::size_t i, const STLMatrix &beta);
	long double d2fYikdxi2__constxk(const std::vector<long double> &x, std::size_t i, std::size_t k, const STLMatrix &beta);
	long double d2fYijdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t k, const STLMatrix &beta);
};

/** \brief Reducing function converter for dry air and HFC blends
 * 
 * From Lemmon, JPCRD, 2000 for the properties of Dry Air, and also from Lemmon, JPCRD, 2004 for the properties of R404A, R410A, etc.
 * \f[
 * \rho_r(\bar x) = \left[ \sum_{i=1}^m\frac{x_i}{\rho_{c_i}}+\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}x_ix_j\zeta_{ij}\right]^{-1}
 * \f]
 * \f[
 * T_r(\bar x) = \sum_{i=1}^mx_iT_{c_i}+\sum_{i=1}^{m-1}\sum_{j=i+1}^mx_ix_j\xi_{ij}
 * \f]
 * 
 * These can be converted to the form of GERG by the following equations:
 * \f[
 * \beta_T = 1\ \ \ \ \beta_v = 1 
 * \f] 
 * and
 * \f[
 * \boxed{\gamma_T = \dfrac{T_{c0}+T_{c1}+\xi_{01}}{2\sqrt{T_{c0}T_{c1}}}}
 * \f]
 * and
 * \f[
 * \boxed{\gamma_v = \dfrac{v_{c0}+v_{c1}+\zeta_{01}}{\frac{1}{4}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}}}
 * \f]
 */
class LemmonAirHFCReducingFunction
{
protected:
    LemmonAirHFCReducingFunction(const LemmonAirHFCReducingFunction &);
public:
	/// Set the coefficients based on reducing parameters loaded from JSON
	static void convert_to_GERG(const std::vector<CoolPropFluid*> &pFluids, std::size_t i, std::size_t j, Dictionary d, long double &beta_T, long double &beta_v, long double &gamma_T, long double &gamma_v)
    {
        long double xi_ij = d.get_number("xi");
        long double zeta_ij = d.get_number("zeta");
	    beta_T = 1;
	    beta_v = 1;
	    gamma_T = (pFluids[i]->pEOS->reduce.T + pFluids[j]->pEOS->reduce.T + xi_ij)/(2*sqrt(pFluids[i]->pEOS->reduce.T*pFluids[j]->pEOS->reduce.T));
	    long double v_i = 1/pFluids[i]->pEOS->reduce.rhomolar;
	    long double v_j = 1/pFluids[j]->pEOS->reduce.rhomolar;
	    long double one_third = 1.0/3.0;
	    gamma_v = (v_i + v_j + zeta_ij)/(0.25*pow(pow(v_i, one_third)+pow(v_j, one_third),3));
    };
};

} /* namespace CoolProp */
#endif