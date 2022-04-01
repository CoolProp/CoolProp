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

namespace CoolProp {

typedef std::vector<std::vector<CoolPropDbl>> STLMatrix;

enum x_N_dependency_flag
{
    XN_INDEPENDENT,  ///< x_N is an independent variable, and not calculated by \f$ x_N = 1-\sum_i x_i\f$
    XN_DEPENDENT     ///< x_N is an dependent variable, calculated by \f$ x_N = 1-\sum_i x_i\f$
};

std::string get_reducing_function_name(const std::string& CAS1, const std::string& CAS2);

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
    ReducingFunction() : N(0){};
    virtual ~ReducingFunction(){};

    virtual ReducingFunction* copy() = 0;

    virtual void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, double value) = 0;

    virtual double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) const = 0;

    /// A factory function to generate the required reducing function
    static shared_ptr<ReducingFunction> factory(const std::vector<CoolPropFluid*>& components, STLMatrix& F);

    /// The reduced temperature
    virtual CoolPropDbl Tr(const std::vector<CoolPropDbl>& x) const = 0;
    /// The derivative of reduced temperature with respect to component i mole fraction
    virtual CoolPropDbl dTrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const = 0;
    /// The molar reducing density
    virtual CoolPropDbl rhormolar(const std::vector<CoolPropDbl>& x) const = 0;
    ///Derivative of the molar reducing density with respect to component i mole fraction
    virtual CoolPropDbl drhormolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl dTr_dgammaT(const std::vector<CoolPropDbl>& x) const {
        throw CoolProp::NotImplementedError("dTr_dgammaT is not implemented for this backend");
    }
    virtual CoolPropDbl dTr_dbetaT(const std::vector<CoolPropDbl>& x) const {
        throw CoolProp::NotImplementedError("dTr_dbetaT is not implemented for this backend");
    }
    virtual CoolPropDbl drhormolar_dgammaV(const std::vector<CoolPropDbl>& x) const {
        throw CoolProp::NotImplementedError("drhormolar_dgammaV is not implemented for this backend");
    }
    virtual CoolPropDbl drhormolar_dbetaV(const std::vector<CoolPropDbl>& x) const {
        throw CoolProp::NotImplementedError("drhormolar_dbetaV is not implemented for this backend");
    }
    virtual CoolPropDbl d2Tr_dxidgammaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        throw CoolProp::NotImplementedError("d2Tr_dxidgammaT is not implemented for this backend");
    }
    virtual CoolPropDbl d2Tr_dxidbetaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        throw CoolProp::NotImplementedError("d2Tr_dxidbetaT is not implemented for this backend");
    }
    virtual CoolPropDbl d2rhormolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        throw CoolProp::NotImplementedError("d2rhormolar_dxidgammaV is not implemented for this backend");
    }
    virtual CoolPropDbl d2rhormolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        throw CoolProp::NotImplementedError("d2rhormolar_dxidbetaV is not implemented for this backend");
    }

    virtual CoolPropDbl d2rhormolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl d2rhormolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl d3rhormolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                             x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl d2Trdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl d2Trdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const = 0;
    virtual CoolPropDbl d3Trdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                      x_N_dependency_flag xN_flag) const = 0;

    /** \brief GERG 2004 Monograph equation 7.56:
     *
     * If the \f$x_i\f$ are all independent
     * \f[
     * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2T_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}-\sum_{k=0}^{N-1}x_k\left(\frac{\partial^2T_r}{\partial x_j \partial x_k}\right)
     * \f]
     * If \f$x_N = 1-\sum x_i\f$:
     * \f[
     * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2T_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}-\sum_{k=0}^{N-1}x_k\left(\frac{\partial^2T_r}{\partial x_j \partial x_k}\right)
     * \f]
     */
    virtual CoolPropDbl d_ndTrdni_dxj__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;

    virtual CoolPropDbl d2_ndTrdni_dxj_dxk__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                    x_N_dependency_flag xN_flag) const;

    /** \brief
     *
     * GERG 2004 Monograph equation 7.55:
     * If the \f$x_i\f$ are all independent
     * \f[
     * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2\rho_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}-\sum_{k=0}^{N-1}x_k\left(\frac{\partial^2\rho_r}{\partial x_j \partial x_k}\right)
     * \f]
     * Gernert, JPCRD, 2014, A28
     * If \f$x_N = 1-\sum x_i\f$:
     * \f[
     * \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2\rho_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}-\sum_{k=0}^{N-2}x_k\left(\frac{\partial^2\rho_r}{\partial x_j \partial x_k}\right)
     * \f]
     */
    virtual CoolPropDbl d_ndrhorbardni_dxj__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                                    x_N_dependency_flag xN_flag) const;

    virtual CoolPropDbl d2_ndrhorbardni_dxj_dxk__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                         x_N_dependency_flag xN_flag) const;

    virtual CoolPropDbl ndrhorbardni__constnj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl ndTrdni__constnj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl PSI_rho(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl d_PSI_rho_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl d2_PSI_rho_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                           x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl PSI_T(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl d_PSI_T_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;
    virtual CoolPropDbl d2_PSI_T_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                         x_N_dependency_flag xN_flag) const;
};

/** \brief The reducing function model of GERG-2008
 *
 * Used by the GERG-2008 formulation to yield the
 * reducing parameters \f$ \rho_r \f$ and \f$ T_r \f$ and derivatives thereof
 */
class GERG2008ReducingFunction : public ReducingFunction
{
   private:
    GERG2008ReducingFunction(const GERG2008ReducingFunction& that);  // No copying
   protected:
    STLMatrix v_c;                       ///< \f$ v_{c,ij} = \frac{1}{8}\left(v_{c,i}^{1/3}+v_{c,j}^{1/3}\right)^{3}\f$ from GERG-2008
    STLMatrix T_c;                       ///< \f$ T_{c,ij} = \sqrt{T_{c,i}T_{c,j}} \f$ from GERG=2008
    STLMatrix beta_v;                    ///< \f$ \beta_{v,ij} \f$ from GERG-2008
    STLMatrix gamma_v;                   ///< \f$ \gamma_{v,ij} \f$ from GERG-2008
    STLMatrix beta_T;                    ///< \f$ \beta_{T,ij} \f$ from GERG-2008
    STLMatrix gamma_T;                   ///< \f$ \gamma_{T,ij} \f$ from GERG-2008
    std::vector<CoolPropDbl> Yc_T;       ///< Vector of critical temperatures for all components
    std::vector<CoolPropDbl> Yc_v;       ///< Vector of critical molar volumes for all components
    std::vector<CoolPropFluid> pFluids;  ///< List of fluids

   public:
    GERG2008ReducingFunction(const std::vector<CoolPropFluid>& pFluids, const STLMatrix& beta_v, const STLMatrix& gamma_v, STLMatrix beta_T,
                             const STLMatrix& gamma_T) {
        this->pFluids = pFluids;
        this->beta_v = beta_v;
        this->gamma_v = gamma_v;
        this->beta_T = beta_T;
        this->gamma_T = gamma_T;
        this->N = pFluids.size();
        T_c.resize(N, std::vector<CoolPropDbl>(N, 0));
        v_c.resize(N, std::vector<CoolPropDbl>(N, 0));
        Yc_T.resize(N);
        Yc_v.resize(N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; j++) {
                T_c[i][j] = sqrt(pFluids[i].EOS().reduce.T * pFluids[j].EOS().reduce.T);
                v_c[i][j] = 1.0 / 8.0 * pow(pow(pFluids[i].EOS().reduce.rhomolar, -1.0 / 3.0) + pow(pFluids[j].EOS().reduce.rhomolar, -1.0 / 3.0), 3);
            }
            Yc_T[i] = pFluids[i].EOS().reduce.T;
            Yc_v[i] = 1 / pFluids[i].EOS().reduce.rhomolar;
        }
    };

    ReducingFunction* copy() {
        return new GERG2008ReducingFunction(pFluids, beta_v, gamma_v, beta_T, gamma_T);
    };

    /// Default destructor
    ~GERG2008ReducingFunction(){};

    /// Set all beta and gamma values in one shot
    void set_binary_interaction_double(const std::size_t i, const std::size_t j, double betaT, double gammaT, double betaV, double gammaV) {
        // bound-check indices
        if (i < 0 || i >= N) {
            if (j < 0 || j >= N) {
                throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
            } else {
                throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
            }
        } else if (j < 0 || j >= N) {
            throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
        }
        beta_T[i][j] = betaT;
        beta_T[j][i] = 1 / betaT;
        gamma_T[i][j] = gammaT;
        gamma_T[j][i] = gammaT;
        beta_v[i][j] = betaV;
        beta_v[j][i] = 1 / betaV;
        gamma_v[i][j] = gammaV;
        gamma_v[j][i] = gammaV;
    }

    /// Set a parameter
    virtual void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, double value) {
        // bound-check indices
        if (i < 0 || i >= N) {
            if (j < 0 || j >= N) {
                throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
            } else {
                throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
            }
        } else if (j < 0 || j >= N) {
            throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
        }
        if (parameter == "betaT") {
            beta_T[i][j] = value;
            beta_T[j][i] = 1 / value;
        } else if (parameter == "gammaT") {
            gamma_T[i][j] = value;
            gamma_T[j][i] = value;
        } else if (parameter == "betaV") {
            beta_v[i][j] = value;
            beta_v[j][i] = 1 / value;
        } else if (parameter == "gammaV") {
            gamma_v[i][j] = value;
            gamma_v[j][i] = value;
        } else {
            throw KeyError(format("This key [%s] is invalid to set_binary_interaction_double", parameter.c_str()));
        }
    }
    /// Get a parameter
    virtual double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) const {
        if (parameter == "betaT") {
            return beta_T[i][j];
        } else if (parameter == "gammaT") {
            return gamma_T[i][j];
        } else if (parameter == "betaV") {
            return beta_v[i][j];
        } else if (parameter == "gammaV") {
            return gamma_v[i][j];
        } else {
            throw KeyError(format("This key [%s] is invalid to get_binary_interaction_double", parameter.c_str()));
        }
    }

    /** \brief The reducing temperature
     * Calculated from \ref Yr with \f$T = Y\f$
     */
    CoolPropDbl Tr(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing temperature with respect to gammaT
     * Calculated from \ref dYr_gamma with \f$T = Y\f$
     */
    CoolPropDbl dTr_dgammaT(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing temperature with respect to betaT
     * Calculated from \ref dYr_beta with \f$T = Y\f$
     */
    CoolPropDbl dTr_dbetaT(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing temperature with respect to gammaT and composition
    */
    CoolPropDbl d2Tr_dxidgammaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing temperature with respect to betaT and composition
    */
    CoolPropDbl d2Tr_dxidbetaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing temperature with respect to component i mole fraction
     *
     * Calculated from \ref dYrdxi__constxj with \f$T = Y\f$
     */
    CoolPropDbl dTrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief The second derivative of reducing temperature with respect to component i mole fraction
     *
     * Calculated from \ref d2Yrdxi2__constxj with \f$T = Y\f$
     */
    CoolPropDbl d2Trdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief The second derivative of reducing temperature with respect to component i and j mole fractions
     *
     * Calculated from \ref d2Yrdxidxj with \f$T = Y\f$
     */
    CoolPropDbl d2Trdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;
    /** \brief The third derivative of reducing temperature with respect to component i, j and k mole fractions
	*
	* Calculated from \ref d3Yrdxidxjdxk with \f$T = Y\f$
	*/
    CoolPropDbl d3Trdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing molar volume with respect to component i mole fraction
     *
     * Calculated from \ref dYrdxi__constxj with \f$v = Y\f$
     */
    CoolPropDbl dvrmolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief The second derivative of reducing molar volume with respect to component i mole fraction
     *
     * Calculated from \ref d2Yrdxi2__constxj with \f$v = Y\f$
     */
    CoolPropDbl d2vrmolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief The second derivative of reducing molar volume with respect to component i and j mole fractions
     *
     * Calculated from \ref d2Yrdxidxj with \f$v = Y\f$
     */
    CoolPropDbl d2vrmolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;
    /** \brief The third derivative of reducing molar volume with respect to component i, j and k mole fractions
	*
	* Calculated from \ref d3Yrdxidxjdxk with \f$v = Y\f$
	*/
    CoolPropDbl d3vrmolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) const;
    /** \brief The molar reducing density
     *
     * Given by \f$ \rho_r = 1/v_r \f$
     */
    CoolPropDbl rhormolar(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing density with respect to gammaV
     * Calculated from \ref dYr_gamma with \f$v = Y\f$
     */
    CoolPropDbl drhormolar_dgammaV(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing density with respect to betaV
     * Calculated from \ref dYr_beta with \f$v = Y\f$
     */
    CoolPropDbl drhormolar_dbetaV(const std::vector<CoolPropDbl>& x) const;

    /** \brief The derivative of reducing volume with respect to gammaV and composition
    */
    CoolPropDbl d2vrmolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing volume with respect to betaV and composition
    */
    CoolPropDbl d2vrmolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing density with respect to betaV and composition
    */
    CoolPropDbl d2rhormolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief The derivative of reducing density with respect to gammaV and composition
    */
    CoolPropDbl d2rhormolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;

    /** \brief Derivative of the molar reducing density with respect to component i mole fraction
     *
     * See also GERG 2004, Eqn. 7.57
     * \f[
     * \left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_{i\neq j}} = -\rho_r^2\left(\frac{\partial v_r}{\partial x_i}\right)_{x_{i\neq j}}
     * \f]
     */
    CoolPropDbl drhormolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief Derivative of the molar reducing density with respect to component i mole fraction
     *
     * See also GERG 2004, Eqn. 7.58
     * \f[
     * \left(\frac{\partial^2 \rho_r}{\partial x_i^2}\right)_{x_{i\neq j}} = 2\rho_r^3\left(\left(\frac{\partial v_r}{\partial x_i}\right)_{x_{i\neq j}}\right)^2-\rho_r\left(\left(\frac{\partial^2 v_r}{\partial x_i^2}\right)_{x_{i\neq j}}\right)
     * \f]
     */
    CoolPropDbl d2rhormolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const;
    /** \brief Derivative of the molar reducing density with respect to component i  and j mole fractions
     *
     * See also GERG 2004, Eqn. 7.59
     * \f[
     * \left(\frac{\partial^2 \rho_r}{\partial x_i\partial x_j}\right) = 2\rho_r^3\left(\left(\frac{\partial v_r}{\partial x_i}\right)_{x_{i\neq j}}\right)\left(\left(\frac{\partial v_r}{\partial x_j}\right)_{x_{i\neq j}}\right)-\rho_r^2\left(\left(\frac{\partial v_r}{\partial x_i\partial x_j}\right)\right)
     * \f]
     */
    CoolPropDbl d2rhormolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const;
    /** \brief Derivative of the molar reducing density with respect to component i, j, and k mole fractions
	*
	*/
    CoolPropDbl d3rhormolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                     x_N_dependency_flag xN_flag) const;

    /** \brief Generalized reducing term \f$Y_r\f$
     *
     * \f[
     * Y_r = \sum_{i=1}^{N}x_iY_{c,i}^2+\sum_{i=1}^{N-1}\sum_{j=i+1}^{N} c_{Y,ij}f_{Y,ij}(x_i,x_j)
     * \f]
     */
    CoolPropDbl Yr(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                   const std::vector<CoolPropDbl>& Yc) const;

    /** \brief Derivative of reducing term \f$Y_r\f$ with respect to \f$\gamma\f$
     *
     * \f[
     * \frac{\partial Y_r}{\partial \gamma} = \sum_{i=1}^{N-1}\sum_{j=i+1}^{N} 2\beta_{ij}Y_{c,ij}f_{Y,ij}(x_i,x_j)
     * \f]
     */
    CoolPropDbl dYr_dgamma(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                           const std::vector<CoolPropDbl>& Yc) const;

    /** \brief Derivative of reducing term \f$Y_r\f$ with respect to \f$\beta\f$
     */
    CoolPropDbl dYr_dbeta(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                          const std::vector<CoolPropDbl>& Yc) const;

    /** \brief First composition derivative of \f$Y_r\f$ with \f$x_i\f$
     *
     * If \f$x_N\f$ is given by \f$ x_N = 1-\sum_{i=1}^{N-1}x_i\f$ (Gernert, FPE, 2014, Table S1):
     * \f{eqnarray*}{
     * \left(\frac{\partial Y_r}{\partial x_i}\right)_{\substack{x_{j\neq i} \\ i<N}} &=& 2x_iY_{c,j}-2x_NY_{c,N} + \sum_{k=1}^{i-1}c_{Y,ki}\frac{\partial f_{Y,ki}(x_k, x_i)}{\partial x_i}+\sum_{k=i+1}^{N-1}c_{Y,ik}\frac{\partial f_{Y,ik}(x_i, x_k)}{\partial x_i} \\
     * &&+c_{Y,iN}\left(\frac{x_{N}(x_i+x_{N})}{\beta_{Y,iN}^2x_i+x_{N}}+(1-\beta_{Y,iN}^2)\frac{x_ix_{N}^2}{(\beta_{Y,i(N)}^2x_i+x_{N})^2}\right) \\
     * &&+\sum_{k=0}^{N-2}c_{Y,kN}\left(-\frac{x_k(x_k+x_N)}{\beta_{Y,kN}^2 x_k+x_N}+(1-\beta_{Y,kN}^2)\frac{x_Nx_k^2}{(\beta_{Y,kN}^2x_k+x_N)^2}\right)
     * \f}
     *
     * Otherwise, if \f$x_i\f$ are all independent:
     * \f[
     * \left(\frac{\partial Y_r}{\partial x_i}\right)_{\substack{x_{j\neq i}}} = 2x_iY_{c,i} + \sum_{k=1}^{i-1}c_{Y,ki}\frac{\partial f_{Y,ki}(x_k,x_i)}{\partial x_i} + \sum_{k=i+1}^{N}c_{Y,ik}\frac{\partial f_{Y,ik}(x_i,x_k)}{\partial x_i}
     * \f]
     *
     */
    CoolPropDbl dYrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                                const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;

    /** \brief Derivative of derivative of reducing term \f$\frac{\partial Y_r}{\partial x_i}\f$ with respect to \f$\beta\f$
    */
    CoolPropDbl d2Yrdxidgamma(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                              const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;
    /** \brief Derivative of derivative of reducing term \f$\frac{\partial Y_r}{\partial x_i}\f$ with respect to \f$\gamma\f$
    */
    CoolPropDbl d2Yrdxidbeta(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                             const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;

    /** \brief Second composition derivative of \f$Y_r\f$ with \f$x_i\f$
     *
     * If \f$x_N\f$ is given by \f$ x_N = 1-\sum_{i=1}^{N-1}x_i\f$ (Gernert, FPE, 2014, Table S1):
     * \f{eqnarray*}{
     * \left(\frac{\partial^2 Y_r}{\partial x_i^2}\right)_{\substack{x_{j\neq i} \\ i<N}} &=& 2(Y_{c,i}+Y_{c,N}) + \sum_{k=1}^{i-1}c_{Y,ki}\frac{\partial^2 f_{Y,ki}(x_k, x_i)}{\partial x_i^2}+\sum_{k=i+1}^{N-1}c_{Y,ik}\frac{\partial^2 f_{Y,ik}(x_i, x_k)}{\partial x_i^2} \\
     * &&+2c_{Y,iN}\left(-\frac{x_i+x_N}{\beta_{Y,iN}^2x_i+x_N}+(1-\beta_{Y,iN}^2)\left( \frac{x_N^2}{(\beta_{Y,iN}^2x_i+x_N)^2} + \frac{(1-\beta_{Y,iN}^2)x_ix_{N}^2-\beta_{Y,iN}^2x_i^2x_N}{(\beta_{Y,iN}^2x_i+x_N)^3}\right)\right) \\
     * &&+\sum_{k=1}^{N-1}2c_{Y,kN}x_k^2\frac{1-\beta_{Y,kN}^2}{(\beta_{Y,kN}^2x_k+x_N)^2} \left(\frac{x_N}{\beta_{Y,kN}^2 x_k+x_N}-1\right)
     * \f}
     *
     * Otherwise, if \f$x_i\f$ are all independent:
     * \f[
     * \left(\frac{\partial^2 Y_r}{\partial x_i^2}\right)_{x_{j\neq i}} = 2Y_{c,i} + \sum_{k=1}^{i-1}c_{Y,ki}\frac{\partial^2 f_{Y,ki}(x_k,x_i)}{\partial x_i^2} + \sum_{k=i+1}^{N}c_{Y,ik}\frac{\partial^2 f_{Y,ik}(x_i,x_k)}{\partial x_i^2}
     * \f]
     */
    CoolPropDbl d2Yrdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                                  const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;
    /** \brief Second mixed composition derivative of \f$Y_r\f$ with \f$x_i\f$ and  \f$x_j\f$
     *
     * If \f$x_N\f$ is given by \f$ x_N = 1-\sum_{i=1}^{N-1}x_i\f$ (Gernert, FPE, 2014, Table S1):
     * \f{eqnarray*}{
     * \left(\frac{\partial^2 Y_r}{\partial x_i\partial x_j}\right)_{\substack{x_{k\neq j\neq i} \\ i<N \\ j < N}} &=& 2Y_{c,N} + c_{Y,ij}\frac{\partial^2 f_{Y,ij}(x_i, x_j)}{\partial x_i\partial x_j}+\sum_{k=1}^{N-1}2c_{Y,kN}x_k^2 \frac{1-\beta_{Y,kN}^2}{(\beta_{Y,kN}^2x_k+x_N)^2} \left(\frac{x_N}{\beta_{Y,kN}^2 x_k+x_N}-1\right) \\
     * &&+c_{Y,iN}\left((1-\beta_{Y,iN}^2)\left(\frac{2x_ix_N^2}{(\beta_{Y,iN}^2x_i+x_N)^3}-\frac{x_ix_N}{(\beta_{Y,iN}^2x_i+x_N)^2}\right) - \frac{x_i+x_N}{\beta_{Y,iN}^2x_i+x_N} \right) \\
     * &&-c_{Y,jN}\left((1-\beta_{Y,jN}^2)\left(\frac{2x_j^2x_N\beta_{Y,jN}^2}{(\beta_{Y,jN}^2x_j+x_N)^3}-\frac{x_jx_N}{(\beta_{Y,jN}^2x_j+x_N)^2}\right) + \frac{x_j+x_N}{\beta_{Y,jN}^2x_j+x_N} \right)
     * \f}
     *
     * Otherwise, if \f$x_i\f$ are all independent:
     * \f[
     * \left(\frac{\partial^2 Y_r}{\partial x_i\partial x_j}\right)_{\substack{x_{k\neq j\neq i}}} = c_{Y,ij}\frac{\partial^2f_{Y,ij}(x_i,x_j)}{\partial x_i\partial x_j}
     * \f]
     */
    CoolPropDbl d2Yrdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta, const STLMatrix& gamma,
                           const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;

    /** \brief Third mixed composition derivative of \f$Y_r\f$ with \f$x_i\f$ and \f$x_j\f$ and \f$x_k\f$
	*
	*/
    CoolPropDbl d3Yrdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k, const STLMatrix& beta,
                              const STLMatrix& gamma, const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const;

    /** \brief The coefficient \f$ c_{Y,ij} \f$
     *
     * \f[
     * c_{Y,ij} = 2\beta_{Y,ij}\gamma_{Y,ij}Y_{c,ij}
     * \f]
     */
    const CoolPropDbl c_Y_ij(const std::size_t i, const std::size_t j, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c) const;

    /** \brief The function \f$ f_{Y,ij}(x_i,x_j) \f$
     *
     * \f[ f_{Y,ij}(x_i,x_j) = x_ix_j\frac{x_i+x_j}{\beta_{Y,ij}^2x_i+x_j} \f]
     */
    CoolPropDbl f_Y_ij(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const;
    /**
     *
     * \f[
     * \left(\frac{\partial f_{Y,ki}(x_k, x_i)}{\partial x_i}\right)_{x_{k\neq i}} = x_k\frac{x_k+x_i}{\beta_{Y,ki}^2x_k+x_i} + \frac{x_kx_i}{\beta_{Y,ki}^2x_k+x_i}\left(1-\frac{x_k+x_i}{\beta_{Y,ki}^2x_k+x_i}\right)
     * \f]
     */
    CoolPropDbl dfYkidxi__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i, const STLMatrix& beta) const;
    /**
     *
     * \f[
     * \left(\frac{\partial f_{Y,ik}(x_i, x_k)}{\partial x_i}\right)_{x_k} = x_k\frac{x_i+x_k}{\beta_{Y,ik}^2x_i+x_k} + \frac{x_ix_k}{\beta_{Y,ik}^2x_i+x_k}\left(1-\beta_{Y,ik}^2\frac{x_i+x_k}{\beta_{Y,ik}^2x_i+x_k}\right)
     * \f]
     */
    CoolPropDbl dfYikdxi__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k, const STLMatrix& beta) const;
    /**
     * \f[
     * \left(\frac{\partial^2 f_{Y,ki}(x_k, x_i)}{\partial x_i^2}\right)_{x_{k\neq i}} = \frac{1}{\beta_{Y,ki}^2x_k+x_i}\left(1-\frac{x_k+x_i}{\beta_{Y,ki}^2x_k+x_i}\right)\left(2x_k-\frac{2x_kx_i}{\beta_{Y,ki}^2x_k+x_i}\right)
     * \f]
     */
    CoolPropDbl d2fYkidxi2__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i, const STLMatrix& beta) const;
    /**
     * \f[
     * \left(\frac{\partial^2 f_{Y,ik}(x_i, x_k)}{\partial x_i^2}\right)_{x_{k}} = \frac{1}{\beta_{Y,ik}^2x_i+x_k}\left(1-\beta_{Y,ik}^2\frac{x_i+x_k}{\beta_{Y,ik}^2x_i+x_k}\right)\left(2x_k-\frac{2x_ix_k\beta_{Y,ik}^2}{\beta_{Y,ik}^2x_i+x_k}\right)
     * \f]
     */
    CoolPropDbl d2fYikdxi2__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k, const STLMatrix& beta) const;
    /**
     * \f{eqnarray*}{
     * \left(\frac{\partial^2 f_{Y,ki}(x_k, x_i)}{\partial x_i\partial x_j}\right)_{x_{k\neq j\neq i}} &=& \frac{x_i+x_j}{\beta_{Y,ij}^2x_i+x_j} + \frac{x_j}{\beta_{Y,ij}^2x_i+x_j}\left(1-\frac{x_i+x_j}{\beta_{Y,ij}^2x_i+x_j}\right) \\
     * &+& \frac{x_i}{\beta_{Y,ij}^2x_i+x_j}\left(1-\beta_{Y,ij}^2\frac{x_i+x_j}{\beta_{Y,ij}^2x_i+x_j}\right) - \frac{x_ix_j}{(\beta_{Y,ij}^2x_i+x_j)^2}\left(1+\beta_{Y,ij}^2-2\beta_{Y,ij}^2\frac{x_i+x_j}{\beta_{Y,ij}^2x_i+x_j}\right)
     * \f}
     */
    CoolPropDbl d2fYijdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k, const STLMatrix& beta) const;

    /* Third order mixed partial derivative of \f$f_{Y,ij}\f$
	  * \f[
	  * \left(\dfrac{\partial^3 f_{Y,ij}(x_i, x_j)}{\partial x_i^2\partial x_j}\right)_{j\neq i} = \dfrac{-6 \beta^{2} x_{i} x_{j}^{2} \left(\beta^{2} - 1\right)}{\beta^{8} x_{i}^{4} + 4 \beta^{6} x_{i}^{3} x_{j} + 6 \beta^{4} x_{i}^{2} x_{j}^{2} + 4 \beta^{2} x_{i} x_{j}^{3} + x_{j}^{4}}
	  * \f]
	  */
    CoolPropDbl d3fYijdxi2dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const;
    /* Third order mixed partial derivative of \f$f_{Y,ij}\f$
	 * \f[
	 * \left(\dfrac{\partial^3 f_{Y,ij}(x_i, x_j)}{\partial x_i\partial x_j^2}\right)_{j\neq i} = \dfrac{6 \beta^{2} x_{i}^{2} x_{j} \left(\beta^{2} - 1\right)}{\beta^{8} x_{i}^{4} + 4 \beta^{6} x_{i}^{3} x_{j} + 6 \beta^{4} x_{i}^{2} x_{j}^{2} + 4 \beta^{2} x_{i} x_{j}^{3} + x_{j}^{4}}
	 * \f]
	 */
    CoolPropDbl d3fYijdxidxj2(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const;
    /* Third order mixed partial derivative of \f$f_{Y,ij}\f$
	 * \f[
	 * \left(\dfrac{ \partial ^ 3 f_{ Y, ki }(x_k, x_i) }{\partial x_i ^ 3}\right)_{ k\neq i } = \dfrac{ \beta_{ Y ki }^{2} x_{ k }^{3} \left(-6 \beta_{ Y ki }^{2} +6\right) }{\left(\beta_{ Y ki }^{2} x_{ k } +x_{ i }\right) ^ { 4 }}
	 * \f]
	 */
    CoolPropDbl d3fYkidxi3__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i, const STLMatrix& beta) const;
    /* Third order mixed partial derivative of \f$f_{Y,ij}\f$
	 * \f[
	 * \left(\dfrac{\partial^3 f_{Y,ik}(x_i, x_k)}{\partial x_i^3}\right)_{k\neq i} = \dfrac{6 \beta_{Y ik}^{2} x_{k}^{3} \left(\beta_{Y ik}^{2} - 1\right)}{\beta_{Y ik}^{8} x_{i}^{4} + 4 \beta_{Y ik}^{6} x_{i}^{3} x_{k} + 6 \beta_{Y ik}^{4} x_{i}^{2} x_{k}^{2} + 4 \beta_{Y ik}^{2} x_{i} x_{k}^{3} + x_{k}^{4}}
	 * \f]
	 */
    CoolPropDbl d3fYikdxi3__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k, const STLMatrix& beta) const;
};

/** \brief A constant reducing function that does not vary with composition.  Think for instance the
 *  reducing function for the cubic EOS
 *
 * Yields the reducing parameters \f$ \rho_r \f$ and \f$ T_r \f$
 */
class ConstantReducingFunction : public ReducingFunction
{
   private:
    ConstantReducingFunction(const ConstantReducingFunction& that);  // No copying
    double T_c, rhomolar_c;

   public:
    ConstantReducingFunction(const double T_c, const double rhomolar_c) : T_c(T_c), rhomolar_c(rhomolar_c){};

    ReducingFunction* copy() {
        return new ConstantReducingFunction(T_c, rhomolar_c);
    };

    void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, double value) {
        return;
    }
    double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) const {
        return _HUGE;
    }

    /// \brief The reducing temperature
    CoolPropDbl Tr(const std::vector<CoolPropDbl>& x) const {
        return T_c;
    };
    /// \brief The derivative of reducing temperature with respect to component i mole fraction
    CoolPropDbl dTrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief The second derivative of reducing temperature with respect to component i mole fraction
    CoolPropDbl d2Trdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief The second derivative of reducing temperature with respect to component i and j mole fractions
    CoolPropDbl d2Trdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief The third derivative of reducing temperature with respect to component i, j and k mole fractions
    CoolPropDbl d3Trdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) const {
        return 0;
    };

    /// \brief The molar reducing density
    CoolPropDbl rhormolar(const std::vector<CoolPropDbl>& x) const {
        return rhomolar_c;
    };
    /// \brief Derivative of the molar reducing density with respect to component i mole fraction
    CoolPropDbl drhormolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief Derivative of the molar reducing density with respect to component i mole fraction
    CoolPropDbl d2rhormolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief Derivative of the molar reducing density with respect to component i  and j mole fractions
    CoolPropDbl d2rhormolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    /// \brief Derivative of the molar reducing density with respect to component i, j, and k mole fractions
    CoolPropDbl d3rhormolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                     x_N_dependency_flag xN_flag) const {
        return 0;
    };

    CoolPropDbl dTr_dgammaT(const std::vector<CoolPropDbl>& x) const {
        return 0;
    }
    CoolPropDbl dTr_dbetaT(const std::vector<CoolPropDbl>& x) const {
        return 0;
    }
    CoolPropDbl drhormolar_dgammaV(const std::vector<CoolPropDbl>& x) const {
        return 0;
    }
    CoolPropDbl drhormolar_dbetaV(const std::vector<CoolPropDbl>& x) const {
        return 0;
    }

    CoolPropDbl d2Tr_dxidgammaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    }
    CoolPropDbl d2Tr_dxidbetaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    }
    CoolPropDbl d2rhormolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    }
    CoolPropDbl d2rhormolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    }

    //virtual CoolPropDbl d_ndTrdni_dxj__constxi(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){ return 0; };
    //virtual CoolPropDbl d2_ndTrdni_dxj_dxk__constxi(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){ return 0; };
    //virtual CoolPropDbl d_ndrhorbardni_dxj__constxi(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){ return 0; };
    //virtual CoolPropDbl d2_ndrhorbardni_dxj_dxk__constxi(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){ return 0; };
    //virtual CoolPropDbl ndrhorbardni__constnj(const std::vector<CoolPropDbl> &x, std::size_t i, x_N_dependency_flag xN_flag){ return 0; };
    //virtual CoolPropDbl ndTrdni__constnj(const std::vector<CoolPropDbl> &x, std::size_t i, x_N_dependency_flag xN_flag){ return 0; };

    /// Note: this one is one, not zero
    virtual CoolPropDbl PSI_rho(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 1;
    };
    virtual CoolPropDbl d_PSI_rho_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    virtual CoolPropDbl d2_PSI_rho_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                           x_N_dependency_flag xN_flag) const {
        return 0;
    };
    virtual CoolPropDbl PSI_T(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    virtual CoolPropDbl d_PSI_T_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
        return 0;
    };
    virtual CoolPropDbl d2_PSI_T_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                         x_N_dependency_flag xN_flag) const {
        return 0;
    };
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
    LemmonAirHFCReducingFunction(const LemmonAirHFCReducingFunction&);

   public:
    /// Set the coefficients based on reducing parameters loaded from JSON
    static void convert_to_GERG(std::vector<CoolPropFluid>& pFluids, std::size_t i, std::size_t j, const Dictionary& d, CoolPropDbl& beta_T,
                                CoolPropDbl& beta_v, CoolPropDbl& gamma_T, CoolPropDbl& gamma_v) {
        CoolPropDbl xi_ij = d.get_number("xi");
        CoolPropDbl zeta_ij = d.get_number("zeta");
        beta_T = 1;
        beta_v = 1;
        gamma_T = (pFluids[i].EOS().reduce.T + pFluids[j].EOS().reduce.T + xi_ij) / (2 * sqrt(pFluids[i].EOS().reduce.T * pFluids[j].EOS().reduce.T));
        CoolPropDbl v_i = 1 / pFluids[i].EOS().reduce.rhomolar;
        CoolPropDbl v_j = 1 / pFluids[j].EOS().reduce.rhomolar;
        CoolPropDbl one_third = 1.0 / 3.0;
        gamma_v = (v_i + v_j + zeta_ij) / (0.25 * pow(pow(v_i, one_third) + pow(v_j, one_third), 3));
    };
};

} /* namespace CoolProp */
#endif
