/**
 * This file contains derivatives needed in the mixture model.  The derivatives are quite nasty, and there
 * are a lot of them, so they are put in this file for cleanness.  The MixtureDerivatives class is a friend class
 * of the HelmholtzEOSMixtureBackend, so it can access all the private members of the HelmholtzEOSMixtureBackend
 * class
*/

// ***************************************************************
// ***************************************************************
// *****************  MIXTURE DERIVATIVES  ***********************
// ***************************************************************
// ***************************************************************

#ifndef MIXTURE_DERIVATIVES_H
#define MIXTURE_DERIVATIVES_H

#include <Eigen/Core>
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "ReducingFunctions.h"

namespace CoolProp {

class HelmholtzEOSMixtureBackend;

/**
This class is a friend class of HelmholtzEOSMixtureBackend, therefore the
static methods contained in it have access to the private and
protected variables in the HelmholtzEOSMixtureBackend instance.

In this way the derivative routines can be kept in their own separate file
and not pollute the HelmholtzEOSMixtureBackend namespace
*/
class MixtureDerivatives
{
   public:
    /** \brief GERG 2004 Monograph equation 7.62
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial p}{\partial V} \right)_{T,\bar n} = -\rho^2 RT(1+2\delta \alpha_{\delta}^r+\delta^2\alpha^r_{\delta\delta})
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     */
    static CoolPropDbl ndpdV__constT_n(HelmholtzEOSMixtureBackend& HEOS);

    /** \brief GERG 2004 Monograph equation 7.61
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial p}{\partial T} \right)_{V,\bar n} = \rho R(1+\delta \alpha_{\delta}^r-\delta \tau \alpha^r_{\delta\tau})
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     */
    static CoolPropDbl dpdT__constV_n(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl dpdrho__constT_n(HelmholtzEOSMixtureBackend& HEOS);

    /** \brief GERG 2004 Monograph equation 7.63
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial p}{\partial n_i} \right)_{T,V,n_j} = \rho RT\left[1+\delta\alpha_{\delta}^r\left[2- \frac{1}{\rho_r}\cdot n\left( \frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] +\delta\cdot n\left(\frac{\partial\alpha_{\delta}^r}{\partial n_i}\right)_{T,V,n_j}\right]
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl ndpdni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 monograph Eqn. 7.32
     *
     * The partial molar volume
     * \f[
     * \hat v_i = \left( \frac{\partial V}{\partial n_i}\right)_{T,p,n_j} = \frac{-\left(\dfrac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\dfrac{\partial p}{\partial V}\right)_{T,\bar n}}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl partial_molar_volume(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief Fugacity of the i-th component
     *
     * Given by the equation
     * \f[
     * f_i(\delta, \tau, \bar x) = x_i\rho R T \exp\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_{j \neq i}}
     * \f]
     */
    static CoolPropDbl fugacity_i(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief Natural logarithm of the fugacity coefficient
     *
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl ln_fugacity_coefficient(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief Derivative of the natural logarithm of the fugacity with respect to T
     *
     * From Witzke, Eqn. 3.14
     * \f[
     *  \left(\frac{\partial \ln(f_i)}{\partial T} \right)_{\rho,x} = -\frac{1}{T}\left(1-\tau\alpha^r_{\tau}-\tau n\left(\frac{\partial\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}}{\partial \tau}\right)_{\delta,\bar x}    \right)
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl dln_fugacity_i_dT__constrho_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /**    \brief Derivative of the natural logarithm of the fugacity with respect to T
     *
     * From Witzke, Eqn. 3.15
     * \f[
     *  \left(\frac{\partial \ln(f_i)}{\partial \rho} \right)_{T, x} = \frac{1}{\rho}\left(1+\delta\alpha^r_{\delta}+\delta n\left(\frac{\partial\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}}{\partial \delta}\right)_{\tau,\bar x}    \right)
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl dln_fugacity_i_drho__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph Eqn. 7.29
     *
     * The derivative term
     *
     * \f[
     * \left(\frac{\partial \ln \phi_i}{\partial T} \right)_{p,\bar n} = \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} + \frac{1}{T}-\frac{\hat v}{RT}\left(\frac{\partial p}{\partial T}\right)_{V,\bar n}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl dln_fugacity_coefficient_dT__constp_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl dln_fugacity_i_dT__constp_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl dln_fugacity_i_dp__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl dln_fugacity_dxj__constT_p_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);
    static CoolPropDbl dln_fugacity_i_dtau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl dln_fugacity_i_ddelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl dln_fugacity_dxj__constT_rho_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    static CoolPropDbl ndln_fugacity_i_dnj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);
    static CoolPropDbl d_ndln_fugacity_i_dnj_dtau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                  x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                  x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_ndln_fugacity_i_dnj_ddelta2__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                  x_N_dependency_flag xN_flag);
    static CoolPropDbl d_ndln_fugacity_i_dnj_ddelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                x_N_dependency_flag xN_flag);
    static CoolPropDbl d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                  x_N_dependency_flag xN_flag);

    static CoolPropDbl d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                   x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                   x_N_dependency_flag xN_flag);

    static CoolPropDbl nd_ndln_fugacity_i_dnj_dnk__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                               x_N_dependency_flag xN_flag);

    static CoolPropDbl nAij(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ndln_fugacity_i_dnj__constT_V_xi(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl n2Aijk(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) {
        return nd_ndln_fugacity_i_dnj_dnk__constT_V_xi(HEOS, i, j, k, xN_flag) - nAij(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d_nAij_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return d_nAij_dX(HEOS, i, j, xN_flag, iTau);
    }
    static CoolPropDbl d_nAij_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return d_nAij_dX(HEOS, i, j, xN_flag, iDelta);
    }
    static CoolPropDbl d_nAij_dX(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag, parameters WRT) {
        if (WRT == iTau) {
            return MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(HEOS, i, j, xN_flag);
        } else if (WRT == iDelta) {
            return MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(HEOS, i, j, xN_flag);
        } else {
            throw ValueError(format("wrong WRT"));
        }
    }
    static CoolPropDbl d_n2Aijk_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) {
        return d_n2Aijk_dX(HEOS, i, j, k, xN_flag, iTau);
    }
    static CoolPropDbl d_n2Aijk_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag) {
        return d_n2Aijk_dX(HEOS, i, j, k, xN_flag, iDelta);
    }
    static CoolPropDbl d_n2Aijk_dX(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag,
                                   parameters WRT) {
        double summer = 0;
        if (WRT == iTau) {
            summer += d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(HEOS, i, j, xN_flag) * ndtaudni__constT_V_nj(HEOS, k, xN_flag);
            summer += d_ndln_fugacity_i_dnj_dtau__constdelta_x(HEOS, i, j, xN_flag) * d_ndtaudni_dTau(HEOS, k, xN_flag);
            summer += d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(HEOS, i, j, xN_flag) * nddeltadni__constT_V_nj(HEOS, k, xN_flag);
            summer += d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(HEOS, i, j, k, xN_flag);
            std::size_t mmax = HEOS.mole_fractions.size();
            if (xN_flag == XN_DEPENDENT) {
                mmax--;
            }
            for (std::size_t m = 0; m < mmax; ++m) {
                summer -= HEOS.mole_fractions[m] * d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(HEOS, i, j, m, xN_flag);
            }
        } else if (WRT == iDelta) {
            summer += d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(HEOS, i, j, xN_flag) * ndtaudni__constT_V_nj(HEOS, k, xN_flag);
            summer += d2_ndln_fugacity_i_dnj_ddelta2__consttau_x(HEOS, i, j, xN_flag) * nddeltadni__constT_V_nj(HEOS, k, xN_flag);
            summer += d_ndln_fugacity_i_dnj_ddelta__consttau_x(HEOS, i, j, xN_flag) * d_nddeltadni_dDelta(HEOS, k, xN_flag);
            summer += d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(HEOS, i, j, k, xN_flag);
            std::size_t mmax = HEOS.mole_fractions.size();
            if (xN_flag == XN_DEPENDENT) {
                mmax--;
            }
            for (std::size_t m = 0; m < mmax; ++m) {
                summer -= HEOS.mole_fractions[m] * d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(HEOS, i, j, m, xN_flag);
            }
        } else {
            return _HUGE;
        }
        return summer - d_nAij_dX(HEOS, i, j, xN_flag, WRT);
    }
    static Eigen::MatrixXd Lstar(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
        std::size_t N = HEOS.mole_fractions.size();
        Eigen::MatrixXd L;
        L.resize(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i; j < N; ++j) {
                L(i, j) = nAij(HEOS, i, j, xN_flag);
            }
        }
        // Fill in the symmetric elements
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                L(i, j) = L(j, i);
            }
        }
        return L;
    }
    static Eigen::MatrixXd dLstar_dX(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag, parameters WRT) {

        std::size_t N = HEOS.mole_fractions.size();
        Eigen::MatrixXd dLstar_dX(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i; j < N; ++j) {
                dLstar_dX(i, j) = d_nAij_dX(HEOS, i, j, xN_flag, WRT);
            }
        }
        // Fill in the symmetric elements
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                dLstar_dX(i, j) = dLstar_dX(j, i);
            }
        }
        return dLstar_dX;
    }
    static Eigen::MatrixXd d2Lstar_dX2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag, parameters WRT1, parameters WRT2) {

        std::size_t N = HEOS.mole_fractions.size();
        Eigen::MatrixXd d2Lstar_dX2(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i; j < N; ++j) {
                if (WRT1 == iTau && WRT2 == iTau) {
                    d2Lstar_dX2(i, j) = MixtureDerivatives::d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(HEOS, i, j, xN_flag);
                } else {
                    throw ValueError(format("d2Lstar_dX2 invalid WRT"));
                }
            }
        }
        // Fill in the symmetric elements
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                d2Lstar_dX2(i, j) = d2Lstar_dX2(j, i);
            }
        }
        return d2Lstar_dX2;
    }
    static Eigen::MatrixXd Mstar(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag, Eigen::MatrixXd& L) {

        std::size_t N = HEOS.mole_fractions.size();
        Eigen::MatrixXd M = L, adjL = adjugate(L);

        // Last row
        for (std::size_t i = 0; i < N; ++i) {
            Eigen::MatrixXd n2dLdni(N, N);
            for (std::size_t j = 0; j < N; ++j) {
                for (std::size_t k = j; k < N; ++k) {
                    n2dLdni(j, k) = n2Aijk(HEOS, j, k, i, xN_flag);
                    // Fill in the symmetric elements
                    n2dLdni(k, j) = n2dLdni(j, k);
                }
            }
            M(N - 1, i) = (adjL * n2dLdni).trace();
        }
        return M;
    }
    static Eigen::MatrixXd dMstar_dX(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag, parameters WRT, Eigen::MatrixXd& L,
                                     Eigen::MatrixXd& dL_dX) {

        std::size_t N = HEOS.mole_fractions.size();
        Eigen::MatrixXd dMstar = dL_dX, adjL = adjugate(L), d_adjL_dX = adjugate_derivative(L, dL_dX);

        // Last row in the d(Mstar)/d(X) requires derivatives of
        for (std::size_t i = 0; i < N; ++i) {
            Eigen::MatrixXd n2dLdni(N, N), d_n2dLdni_dX(N, N);
            for (std::size_t j = 0; j < N; ++j) {
                for (std::size_t k = j; k < N; ++k) {
                    n2dLdni(j, k) = n2Aijk(HEOS, j, k, i, xN_flag);
                    d_n2dLdni_dX(j, k) = d_n2Aijk_dX(HEOS, j, k, i, xN_flag, WRT);
                    // Fill in the symmetric elements
                    n2dLdni(k, j) = n2dLdni(j, k);
                    d_n2dLdni_dX(k, j) = d_n2dLdni_dX(j, k);
                }
            }
            dMstar(N - 1, i) = (n2dLdni * d_adjL_dX + adjL * d_n2dLdni_dX).trace();
        }
        return dMstar;
    }

    /** \brief Table B4, Kunz, JCED, 2012 for the original term and the subsequent substitutions
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j}
     * \f]
     * which is equal to
     * \f{eqnarray*}{
     * n\left(\frac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} &=& \delta \alpha^r_{\delta}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
     * && +\tau \alpha^r_{\tau}\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
     * && +\alpha^r_{x_i}-\sum_{k=1}^{N}x_k\alpha^r_{x_k}
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl ndalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /// GERG Equation 7.42
    static CoolPropDbl dnalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph Eqn. 7.30
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial \ln \phi_i}{\partial p} \right)_{T,\bar n} = \frac{\hat v_i}{RT}-\frac{1}{p}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl dln_fugacity_coefficient_dp__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph Equation 7.31
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1+\frac{n}{RT}\frac{\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\left(\frac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\frac{\partial p}{\partial V}\right)_{T,\bar n}}
     * \f]
     * which is also equal to
     * \f[
     * n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1-\frac{\hat v_i}{RT}\left[n\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\right]
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param j The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl ndln_fugacity_coefficient_dnj__constT_p(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);

    /** \brief Gernert Equation 3.115
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial \ln \phi_i}{\partial x_j}\right)_{T,p,x_{k\neq j}} = \left(\frac{\partial^2n\alpha^r}{\partial x_j \partial n_i} \right)_{T,V}+\frac{1}{RT}\frac{\left(\frac{\partial p}{\partial n_i}\right)_{T,V,n_{k\neq i}}\left(\frac{\partial p}{\partial x_j}\right)_{T,V,x_{k\neq j}}}{\left(\frac{\partial p}{\partial V}\right)_{T,\bar n}}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The first index of interest
     * @param j The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     *
    */
    static CoolPropDbl dln_fugacity_coefficient_dxj__constT_p_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag);

    /** \brief Gernert Equation 3.130
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial p}{\partial x_j} \right)_{T,V,x_{k\neq j}} = \rho RT\left(-\frac{1}{\rho_r}\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_{k\neq j}} \delta\alpha_{\delta}^r + \delta\left(\frac{\partial}{\partial x_j}\left(\left( \frac{\partial \alpha^r}{\partial \delta}\right)_{\tau,\bar x}\right)\right)_{T,V,x_{k\neq j}}\right)
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param j The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl dpdxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag);

    /** \brief Gernert Equation 3.117
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial^2n\alpha^r}{\partial x_j\partial n_i} \right)_{T,V} = \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_{j\neq i}}\right)\right)_{T,V,x_{k\neq j}} +\left(\frac{\partial \alpha^r}{\partial x_j}\right)_{T,V,x_{k\neq j}}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param j The first index of interest
     * @param i The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl d2nalphar_dxj_dni__constT_V(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, std::size_t i, x_N_dependency_flag xN_flag) {
        return MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(HEOS, i, j, xN_flag)
               + MixtureDerivatives::dalphar_dxj__constT_V_xi(HEOS, j, xN_flag);
    };

    /// Gernert Equation 3.119
    /// Catch test provided
    static CoolPropDbl dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag);

    /// Gernert Equation 3.118
    /// Catch test provided
    static CoolPropDbl d_ndalphardni_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /// Gernert Equation 3.134
    /// Catch test provided
    static CoolPropDbl d_dalpharddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag);

    /// Gernert Equation 3.121
    /// Catch test provided
    static CoolPropDbl ddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag);

    /// Gernert Equation 3.122
    /// Catch test provided
    static CoolPropDbl dtau_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph, equations 7.44 and 7.51
     *
     * The derivative term
     * \f[
     * \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = \left( \frac{\partial}{\partial T}\left(\frac{\partial n \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right)_{V,\bar n}
     * \f]
     * \f[
     * \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = -\frac{\tau}{T}\left[\alpha_{\tau}^r +\left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\bar x}\right]
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl d2nalphar_dni_dT(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph Equation 7.51 and Table B4, Kunz, JCED, 2012
     *
     * The derivative term
     * \f{eqnarray*}{
     * \frac{\partial }{\partial \tau} \left( n\left(\frac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right) &=& \delta \alpha^r_{\delta\tau}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
     * && +(\tau \alpha^r_{\tau\tau}+\alpha^r_{\tau})\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
     * && +\alpha^r_{x_i\tau}-\sum_{k=1}^{N}x_k\alpha^r_{x_k\tau}
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl d_ndalphardni_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief \f$\tau\f$ derivartive of GERG 2004 Monograph Equation 7.51
     *
     * The derivative term
     * \f[
     * \dfrac{\partial^2 }{\partial \tau^2} \left( n\left(\dfrac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right) = \delta \alpha^r_{\delta\tau\tau}\Psi_{\rho}+ (2\alpha^r_{\tau\tau}+\tau\alpha^r_{\tau\tau\tau})\Psi_T+\alpha^r_{x_i\tau\tau}-\sum_{k=1}^{N}x_k\alpha^r_{x_k\tau\tau}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl d2_ndalphardni_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d3_ndalphardni_dTau3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph Equation 7.50 and Table B4, Kunz, JCED, 2012
     *
     * The derivative term
     * \f{eqnarray*}{
     * \left(\frac{\partial }{\partial \delta} \left( n\left(\frac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\bar x} &=& (\alpha_{\delta}^r+\delta\alpha_{\delta\delta}^r)\left[1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j} \right] \\
     * &+&\tau\alpha^r_{\delta\tau}\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
     * &+&\alpha^r_{\delta x_i}-\sum_{k=1}^{N}x_k\alpha^r_{\delta x_k}
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief \f$\delta\f$ derivative of GERG 2004 Monograph Equation 7.50
    *
    * The derivative term
    * \f[
    * \begin{array}{ccl}
    * \left(\dfrac{\partial^2 }{\partial \delta^2} \left( n\left(\dfrac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\bar x} &=& (2\alpha_{\delta\delta}^r+\delta\alpha_{\delta\delta\delta}^r)\Psi_{\rho} +\tau\alpha^r_{\delta\delta\tau}\Psi_T+\alpha^r_{\delta\delta x_i}-\sum_{k=1}^{N}x_k\alpha^r_{\delta\delta x_k}
    * \end{array}
    * \f]
    * @param HEOS The HelmholtzEOSMixtureBackend to be used
    * @param i The index of interest
    * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl d2_ndalphardni_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d3_ndalphardni_dDelta3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief \f$\tau\f$ derivative of GERG 2004 Monograph Equation 7.50
    *
    * The derivative term
    * \f[
    * \begin{array}{ccl}
    * \left(\dfrac{\partial^2 }{\partial \delta\partial \tau} \left( n\left(\dfrac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\bar x} &=& (\alpha_{\delta\tau}^r+\delta\alpha_{\delta\delta\tau}^r)\Psi_{\rho} +(\tau\alpha^r_{\delta\tau\tau} + \alpha^r_{\delta\tau})\Psi_T
+\alpha^r_{\delta\tau x_i}-\sum_{k=1}^{N}x_k\alpha^r_{\delta\tau x_k}
    * \end{array}
    * \f]
    * @param HEOS The HelmholtzEOSMixtureBackend to be used
    * @param i The index of interest
    * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl d2_ndalphardni_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d3_ndalphardni_dDelta2_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl d3_ndalphardni_dDelta_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /** \brief \f$x_j\f$ derivative of GERG 2004 Monograph Equation 7.50
    *
    * The derivative term
    * \f[
    * \begin{array}{ccl}
    * \dfrac{\partial}{\partial x_j}\left(\dfrac{\partial }{\partial \delta} \left( n\left(\dfrac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right)_{\tau,\bar x}\right)_{\tau, \delta, x_i} &=& (\alpha_{\delta}^r+\delta\alpha_{\delta\delta}^r)\dfrac{\partial\Psi_{\rho}}{\partial x_j} + (\alpha_{\delta x_j}^r+\delta\alpha_{\delta\delta x_j}^r)\Psi_{\rho}\\
    * &+&\tau\alpha^r_{\delta\tau}\dfrac{\partial\Psi_{T}}{\partial x_j} + \tau\alpha^r_{\delta\tau x_j}\Psi_T\\
    * &+&\alpha^r_{\delta x_i x_j}-\sum_{k=1}^{N}\left[x_k\alpha^r_{\delta x_k x_j} + \dfrac{d x_k}{d x_j}\alpha^r_{\delta x_k}\right]
    * \end{array}
    * \f]
    * @param HEOS The HelmholtzEOSMixtureBackend to be used
    * @param i The index of interest
    * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl d2_ndalphardni_dxj_dDelta__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                              x_N_dependency_flag xN_flag);

    /** \brief \f$x_j\f$ derivative of GERG 2004 Monograph Equation 7.51
    *
    * The derivative term
    * \f[
    * \begin{array}{ccl}
    * \dfrac{\partial}{\partial x_j}\left(\dfrac{\partial }{\partial \tau} \left( n\left(\dfrac{\partial \alpha^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\delta,x_i} &=& \delta \alpha^r_{\delta\tau}\dfrac{\partial \Psi_{\rho}}{\partial x_j} + \delta \alpha^r_{\delta\tau x_j}\Psi_{\rho}\\
    * &+& (\tau \alpha^r_{\tau\tau}+\alpha^r_{\tau})\dfrac{\partial \Psi_T}{\partial x_j}+ \left(\tau \alpha^r_{\tau\tau x_j} +\alpha^r_{\tau x_j}\right)\Psi_T\\
    * &+&\alpha^r_{\tau x_i x_j}-\sum_{k=1}^{N}\left[x_k\alpha^r_{\tau x_k x_j} + \dfrac{d x_k}{d x_j}\alpha^r_{\tau x_k}\right]
    * \end{array}
    * \f]
    * @param HEOS The HelmholtzEOSMixtureBackend to be used
    * @param i The index of interest
    * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
    */
    static CoolPropDbl d2_ndalphardni_dxj_dTau__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                              x_N_dependency_flag xN_flag);
    static CoolPropDbl d3_ndalphardni_dxj_dTau2__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);
    static CoolPropDbl d3_ndalphardni_dxj_dDelta2__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);
    static CoolPropDbl d3_ndalphardni_dxj_dDelta_dTau__constxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph equation 7.41
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial^2n\alpha^r}{\partial n_i \partial n_j} \right)_{T,V} = n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i}
     * \f]
     * and
     * GERG 2004 Monograph equation 7.46:
     * \f[
     * n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i} = n\left( \frac{\partial \alpha^r}{\partial n_j}\right)_{T,V,n_i} + n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i}
     * \f]
     * GERG 2004 Monograph equation 7.47:
     * \f{eqnarray*}{
     * n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i} &=& \left( \frac{\partial}{\partial \delta}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\delta}{\partial n_j}\right)_{T,V,n_i}\\
     * &+& \left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\tau}{\partial n_j}\right)_{T,V,n_i}\\
     * &+& \left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}-\sum_{k=1}^{N}x_k \left( \frac{\partial}{\partial x_k}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}\\
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param j The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /* \brief GERG 2004 Eqn. 7.47
     * \f{eqnarray*}{
     * n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i} &=& \left( \frac{\partial}{\partial \delta}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\delta}{\partial n_j}\right)_{T,V,n_i}\\
     * &+& \left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\tau}{\partial n_j}\right)_{T,V,n_i}\\
     * &+& \left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}-\sum_{k=1}^{N}x_k \left( \frac{\partial}{\partial x_k}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}\\
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param j The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl nd_ndalphardni_dnj__constT_V(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /* \brief \f$\tau\f$ derivative of GERG 2004 7.47
     *
     */
    static CoolPropDbl d_nd_ndalphardni_dnj_dTau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_nd_ndalphardni_dnj_dTau2__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_nd_ndalphardni_dnj_dDelta2__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_nd_ndalphardni_dnj_dDelta_dTau__constx(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag);

    /* \brief \f$\delta\f$ derivative of GERG 2004 7.47
     *
     */
    static CoolPropDbl d_nd_ndalphardni_dnj_dDelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag);

    /* \brief \f$x_k\f$ derivative of GERG 2004 7.47
    *
    */
    static CoolPropDbl d_nd_ndalphardni_dnj_dxk__consttau_delta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                x_N_dependency_flag xN_flag);

    static CoolPropDbl d2_nd_ndalphardni_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                  x_N_dependency_flag xN_flag);
    static CoolPropDbl d2_nd_ndalphardni_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                  x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph equation 7.48
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial \delta}{\partial n_i} \right)_{T,V,n_j} = \delta - \frac{\delta}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl nddeltadni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d_nddeltadni_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d_nddeltadni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    static CoolPropDbl d2_nddeltadni_dxj_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                          x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph equation 7.49
     *
     * The derivative term
     * \f[
     * n\left(\frac{\partial \tau}{\partial n_i} \right)_{T,V,n_j} = \frac{\tau}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}
     * \f]
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d_ndtaudni_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    static CoolPropDbl d_ndtaudni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    static CoolPropDbl d2_ndtaudni_dxj_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /** \brief GERG 2004 Monograph equation 7.52
     *
     * The derivative term
     * \f{eqnarray*}{
     * \left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i} &=& \delta\alpha_{\delta x_j}^{r}\left[ 1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
     * &-& \delta\alpha_{\delta}^{r}\frac{1}{\rho_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{\rho_r}\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
     * &+& \tau\alpha_{\tau x_j}^r\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
     * &+& \tau\alpha_{\tau}^{r}\frac{1}{T_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{T_r}\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right] \\
     * &+& \alpha_{x_ix_j}^r-\alpha_{x_j}^r-\sum_{m=1}^Nx_m\alpha_{x_jx_m}^r
     * \f}
     * @param HEOS The HelmholtzEOSMixtureBackend to be used
     * @param i The index of interest
     * @param j The second index of interest
     * @param xN_flag A flag specifying whether the all mole fractions are independent or only the first N-1
     */
    static CoolPropDbl d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                            x_N_dependency_flag xN_flag);

    /* \brief \f$x_k\f$ derivative of GERG 2004 7.52
    *
    */
    static CoolPropDbl d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                 x_N_dependency_flag xN_flag);

    static CoolPropDbl d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                  x_N_dependency_flag xN_flag);
    static CoolPropDbl d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                  x_N_dependency_flag xN_flag);

    static CoolPropDbl dalpha0_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl d2alpha0_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl d2alpha0_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);
    static CoolPropDbl d2alpha0dxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /// Return the Helmholtz energy density \f$\psi = \rho a\f$
    static CoolPropDbl psi(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return (HEOS.alphar() + HEOS.alpha0()) * HEOS.rhomolar() * HEOS.gas_constant() * HEOS.T();
    }
    /// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\delta\f$
    static CoolPropDbl dpsi_dDelta(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\tau\f$
    static CoolPropDbl dpsi_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the first partial of residual Helmholtz energy density \f$\psi^{\rm r} = \rho a^{\rm r}\f$ with respect to \f$\tau\f$
    static CoolPropDbl dpsir_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$x_i\f$
    static CoolPropDbl dpsi_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /// Return the first partial of Helmholtz energy density \f$\psi^{\rm r} = \rho a^{\rm r}\f$ with respect to \f$x_i\f$
    static CoolPropDbl dpsir_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /// Return the first partial of the product \f$\rho_{\rm r}T_{\rm r}\f$ with respect to \f$x_i\f$
    static CoolPropDbl d_rhorTr_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag);

    /// Return the second partial of the product \f$\rho_{\rm r}T_{\rm r}\f$ with respect to \f$x_i\f$ and \f$x_j\f$
    static CoolPropDbl d2_rhorTr_dxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag);

    /// Return the second partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\delta\f$
    static CoolPropDbl d2psi_dDelta2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second cross partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\delta\f$ and \f$\tau\f$
    static CoolPropDbl d2psi_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second cross partial of Helmholtz energy density \f$\psi^{\rm r}  = \rho a^{\rm r} \f$ with respect to \f$\delta\f$ and \f$\tau\f$
    static CoolPropDbl d2psir_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\tau\f$
    static CoolPropDbl d2psi_dTau2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of residual Helmholtz energy density \f$\psi^{\rm r} = \rho a^{\rm r}\f$ with respect to \f$\tau\f$
    static CoolPropDbl d2psir_dTau2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\tau\f$ and \f$x_i\f$
    static CoolPropDbl d2psi_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of residual Helmholtz energy density \f$\psi^{\rm r} = \rho a^{\rm r}\f$ with respect to \f$\tau\f$ and \f$x_i\f$
    static CoolPropDbl d2psir_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\delta\f$ and \f$x_i\f$
    static CoolPropDbl d2psi_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$x_i\f$ and \f$x_j\f$
    static CoolPropDbl d2psi_dxi_dxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// Return the second partial of residual Helmholtz energy density \f$\psi^{\rm r} = \rho a^{\rm r}\f$ with respect to \f$x_i\f$ and \f$x_j\f$
    static CoolPropDbl d2psir_dxi_dxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag = XN_INDEPENDENT);

    /// ****************************************************************************
    /// ****************************************************************************
    /// ****************************************************************************
    ///                 Shim functions for testing of derivatives
    /// ****************************************************************************
    /// ****************************************************************************
    /// ****************************************************************************
    /// (these are needed because this class is a friend of HelmholtzEOSMixtureBackend and therefore has access to class private variables)

    static CoolPropDbl PSI_rho(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    }
    static CoolPropDbl d_PSI_rho_dxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                     CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    }
    static CoolPropDbl d2_PSI_rho_dxj_dxk(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                          CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    }
    static CoolPropDbl PSI_T(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    }
    static CoolPropDbl d_PSI_T_dxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    }
    static CoolPropDbl d2_PSI_T_dxj_dxk(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                        CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    }

    static CoolPropDbl alpha0(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.calc_alpha0_deriv_nocache(0, 0, HEOS.mole_fractions, HEOS.tau(), HEOS.delta(), HEOS.T_reducing(), HEOS.rhomolar_reducing());
    }

    static CoolPropDbl alphar(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        bool cache_values = false;
        HelmholtzDerivatives HD = HEOS.residual_helmholtz->all(HEOS, HEOS.mole_fractions, HEOS.tau(), HEOS.delta(), cache_values);
        return HD.alphar;
    }
    static CoolPropDbl dalphar_dxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->dalphar_dxi(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alphar_dxi_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, i, xN_flag);
    }
    static CoolPropDbl d3alphar_dxi_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, i, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dTau3(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alphar_dxi_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, i, xN_flag);
    }
    static CoolPropDbl d3alphar_dxi_dDelta2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, i, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dDelta3(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, i, xN_flag);
    }
    static CoolPropDbl d3alphar_dxi_dDelta_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, i, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dDelta2_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, i, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dDelta_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alphardxidxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                      CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d2alphardxidxj(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dxj_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                              CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dxj_dDelta2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d4alphar_dxi_dxj_dDelta_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                    CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d3alphar_dxi_dxj_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                               CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d3alphar_dxi_dxj_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                             CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl d3alphardxidxjdxk(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                         CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.residual_helmholtz->d3alphardxidxjdxk(HEOS, i, j, k, xN_flag);
    }
    static CoolPropDbl tau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.tau();
    }
    static CoolPropDbl Tr(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->Tr(HEOS.get_mole_fractions());
    }
    static CoolPropDbl dTrdxi__constxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->dTrdxi__constxj(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d2Trdxi2__constxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2Trdxi2__constxj(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d2Trdxidxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2Trdxidxj(HEOS.get_mole_fractions(), i, j, xN_flag);
    }
    static CoolPropDbl d3Trdxidxjdxk(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                     CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d3Trdxidxjdxk(HEOS.get_mole_fractions(), i, j, k, xN_flag);
    }
    static CoolPropDbl ndTrdni__constnj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->ndTrdni__constnj(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d_ndTrdni_dxj__constxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                              CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.get_mole_fractions(), i, j, xN_flag);
    }
    static CoolPropDbl d2_ndTrdni_dxj_dxk__constxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                   CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2_ndTrdni_dxj_dxk__constxi(HEOS.get_mole_fractions(), i, j, k, xN_flag);
    }
    static CoolPropDbl dTr_dgammaT(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->dTr_dgammaT(HEOS.get_mole_fractions());
    }
    static CoolPropDbl dTr_dbetaT(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->dTr_dbetaT(HEOS.get_mole_fractions());
    }
    static CoolPropDbl d2Tr_dxidgammaT(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2Tr_dxidgammaT(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d2Tr_dxidbetaT(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2Tr_dxidbetaT(HEOS.get_mole_fractions(), i, xN_flag);
    }

    static CoolPropDbl delta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.delta();
    }
    static CoolPropDbl rhormolar(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->rhormolar(HEOS.get_mole_fractions());
    }
    static CoolPropDbl drhormolardxi__constxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->drhormolardxi__constxj(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d2rhormolardxidxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                         CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2rhormolardxidxj(HEOS.get_mole_fractions(), i, j, xN_flag);
    }
    static CoolPropDbl d3rhormolardxidxjdxk(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                            CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d3rhormolardxidxjdxk(HEOS.get_mole_fractions(), i, j, k, xN_flag);
    }
    static CoolPropDbl drhormolar_dgammaV(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->drhormolar_dgammaV(HEOS.get_mole_fractions());
    }
    static CoolPropDbl drhormolar_dbetaV(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->drhormolar_dbetaV(HEOS.get_mole_fractions());
    }
    static CoolPropDbl d2rhormolar_dxidgammaV(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2rhormolar_dxidgammaV(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d2rhormolar_dxidbetaV(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2rhormolar_dxidbetaV(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl ndrhorbardni__constnj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->ndrhorbardni__constnj(HEOS.get_mole_fractions(), i, xN_flag);
    }
    static CoolPropDbl d_ndrhorbardni_dxj__constxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                   CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.get_mole_fractions(), i, j, xN_flag);
    }
    static CoolPropDbl d2_ndrhorbardni_dxj_dxk__constxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                        CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.Reducing->d2_ndrhorbardni_dxj_dxk__constxi(HEOS.get_mole_fractions(), i, j, k, xN_flag);
    }
    static CoolPropDbl p(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag) {
        return HEOS.p();
    }

    static CoolPropDbl dalpha0_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalpha0_dDelta();
    }
    static CoolPropDbl d2alpha0_dDelta2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alpha0_dDelta2();
    }
    static CoolPropDbl dalpha0_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalpha0_dTau();
    }
    static CoolPropDbl d2alpha0_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alpha0_dTau2();
    }

    static CoolPropDbl dalphar_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalphar_dDelta();
    }
    static CoolPropDbl d2alphar_dDelta2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dDelta2();
    }
    static CoolPropDbl d2alphar_dDelta_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dDelta_dTau();
    }
    static CoolPropDbl dalphar_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalphar_dTau();
    }
    static CoolPropDbl d2alphar_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dTau2();
    }

    static CoolPropDbl alpha(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.alphar() + alpha0(HEOS, xN_flag);
    }
    static CoolPropDbl dalpha_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalphar_dDelta() + HEOS.dalpha0_dDelta();
    }
    static CoolPropDbl dalpha_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.dalphar_dTau() + HEOS.dalpha0_dTau();
    }
    static CoolPropDbl d2alpha_dDelta2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dDelta2() + HEOS.d2alpha0_dDelta2();
    }
    static CoolPropDbl d2alpha_dDelta_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dDelta_dTau() + HEOS.d2alpha0_dDelta_dTau();
    }
    static CoolPropDbl d2alpha_dTau2(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag = XN_INDEPENDENT) {
        return HEOS.d2alphar_dTau2() + HEOS.d2alpha0_dTau2();
    }
    static CoolPropDbl dalpha_dxi(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return dalphar_dxi(HEOS, i, xN_flag) + dalpha0_dxi(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alpha_dxi_dDelta(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return d2alphar_dxi_dDelta(HEOS, i, xN_flag) + d2alpha0_dxi_dDelta(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alpha_dxi_dTau(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        return d2alphar_dxi_dTau(HEOS, i, xN_flag) + d2alpha0_dxi_dTau(HEOS, i, xN_flag);
    }
    static CoolPropDbl d2alphadxidxj(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                     CoolProp::x_N_dependency_flag xN_flag) {
        return d2alphardxidxj(HEOS, i, j, xN_flag) + d2alpha0dxidxj(HEOS, i, j, xN_flag);
    }
    static CoolPropDbl ln_fugacity(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag) {
        double f_i = MixtureDerivatives::fugacity_i(HEOS, i, xN_flag);
        return log(f_i);
    }

}; /* class MixtureDerivatives */

} /* namespace CoolProp*/
#endif