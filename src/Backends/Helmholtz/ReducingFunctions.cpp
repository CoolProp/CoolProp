#include "ReducingFunctions.h"

namespace CoolProp {

CoolPropDbl ReducingFunction::d_ndTrdni_dxj__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                                     x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        CoolPropDbl s = 0;
        for (std::size_t k = 0; k < N; k++) {
            s += x[k] * d2Trdxidxj(x, j, k, xN_flag);
        }
        return d2Trdxidxj(x, i, j, xN_flag) - dTrdxi__constxj(x, j, xN_flag) - s;
    } else if (xN_flag == XN_DEPENDENT) {
        if (j == N - 1) {
            return 0;
        }
        if (N == 0) {
            return 0;
        }
        CoolPropDbl s = 0;
        for (std::size_t k = 0; k < N - 1; k++) {
            s += x[k] * d2Trdxidxj(x, k, j, xN_flag);
        }
        CoolPropDbl val = d2Trdxidxj(x, j, i, xN_flag) - dTrdxi__constxj(x, j, xN_flag) - s;
        return val;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl ReducingFunction::d2_ndTrdni_dxj_dxk__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                          x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        CoolPropDbl s = 0;
        for (std::size_t m = 0; m < N; m++) {
            s += x[m] * d3Trdxidxjdxk(x, j, k, m, xN_flag);
        }
        return d3Trdxidxjdxk(x, i, j, k, xN_flag) - 2 * d2Trdxidxj(x, j, k, xN_flag) - s;
    } else if (xN_flag == XN_DEPENDENT) {
        if (N == 0) {
            return 0;
        }
        if (j == N - 1) {
            return 0;
        }
        CoolPropDbl s = 0;
        for (std::size_t m = 0; m < N - 1; m++) {
            s += x[m] * d3Trdxidxjdxk(x, k, j, m, xN_flag);
        }
        CoolPropDbl val = d3Trdxidxjdxk(x, i, j, k, xN_flag) - d2Trdxidxj(x, j, k, xN_flag) - s;
        return val;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl ReducingFunction::d_ndrhorbardni_dxj__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                                          x_N_dependency_flag xN_flag) const {
    CoolPropDbl s = 0;
    for (std::size_t k = 0; k < N; k++) {
        s += x[k] * d2rhormolardxidxj(x, j, k, xN_flag);
    }
    return d2rhormolardxidxj(x, j, i, xN_flag) - drhormolardxi__constxj(x, j, xN_flag) - s;
}
CoolPropDbl ReducingFunction::d2_ndrhorbardni_dxj_dxk__constxi(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                               x_N_dependency_flag xN_flag) const {
    CoolPropDbl s = 0;
    for (std::size_t m = 0; m < N; m++) {
        s += x[m] * d3rhormolardxidxjdxk(x, j, k, m, xN_flag);
    }
    return d3rhormolardxidxjdxk(x, i, j, k, xN_flag) - 2 * d2rhormolardxidxj(x, j, k, xN_flag) - s;
}
CoolPropDbl ReducingFunction::ndrhorbardni__constnj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        CoolPropDbl summer_term1 = 0;
        for (std::size_t j = 0; j < N; j++) {
            summer_term1 += x[j] * drhormolardxi__constxj(x, j, xN_flag);
        }
        return drhormolardxi__constxj(x, i, xN_flag) - summer_term1;
    } else if (xN_flag == XN_DEPENDENT) {
        CoolPropDbl summer_term1 = 0;
        if (N == 0) {
            return 0;
        }
        for (std::size_t k = 0; k < N - 1; ++k) {
            summer_term1 += x[k] * drhormolardxi__constxj(x, k, xN_flag);
        }
        return drhormolardxi__constxj(x, i, xN_flag) - summer_term1;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl ReducingFunction::ndTrdni__constnj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        // GERG Equation 7.54
        CoolPropDbl summer_term1 = 0;
        for (std::size_t j = 0; j < N; j++) {
            summer_term1 += x[j] * dTrdxi__constxj(x, j, xN_flag);
        }
        return dTrdxi__constxj(x, i, xN_flag) - summer_term1;
    } else if (xN_flag == XN_DEPENDENT) {
        CoolPropDbl summer_term1 = 0;
        if (N == 0) {
            return 0;
        }
        for (std::size_t k = 0; k < N - 1; ++k) {
            summer_term1 += x[k] * dTrdxi__constxj(x, k, xN_flag);
        }
        return dTrdxi__constxj(x, i, xN_flag) - summer_term1;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl ReducingFunction::PSI_rho(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return 1 - 1 / rhormolar(x) * ndrhorbardni__constnj(x, i, xN_flag);
}
CoolPropDbl ReducingFunction::d_PSI_rho_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
    return -1 / rhormolar(x) * (d_ndrhorbardni_dxj__constxi(x, i, j, xN_flag) - drhormolardxi__constxj(x, j, xN_flag) * (1 - PSI_rho(x, i, xN_flag)));
}
CoolPropDbl ReducingFunction::d2_PSI_rho_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                 x_N_dependency_flag xN_flag) const {
    double line1 = d2_ndrhorbardni_dxj_dxk__constxi(x, i, j, k, xN_flag);
    double line2 = -1 / rhormolar(x) * drhormolardxi__constxj(x, k, xN_flag) * d_ndrhorbardni_dxj__constxi(x, i, j, xN_flag);
    double line3 = drhormolardxi__constxj(x, j, xN_flag) * d_PSI_rho_dxj(x, i, k, xN_flag);
    double line4 =
      -(d2rhormolardxidxj(x, j, k, xN_flag) - 1 / rhormolar(x) * drhormolardxi__constxj(x, k, xN_flag) * drhormolardxi__constxj(x, j, xN_flag))
      * (1 - PSI_rho(x, i, xN_flag));
    return -1 / rhormolar(x) * (line1 + line2 + line3 + line4);
}
CoolPropDbl ReducingFunction::PSI_T(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return 1 / Tr(x) * ndTrdni__constnj(x, i, xN_flag);
}
CoolPropDbl ReducingFunction::d_PSI_T_dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
    return 1 / Tr(x) * (d_ndTrdni_dxj__constxi(x, i, j, xN_flag) - dTrdxi__constxj(x, j, xN_flag) * PSI_T(x, i, xN_flag));
}
CoolPropDbl ReducingFunction::d2_PSI_T_dxj_dxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                               x_N_dependency_flag xN_flag) const {
    double line1 = d2_ndTrdni_dxj_dxk__constxi(x, i, j, k, xN_flag);
    double line2 = -1 / Tr(x) * dTrdxi__constxj(x, k, xN_flag) * d_ndTrdni_dxj__constxi(x, i, j, xN_flag);
    double line3 = -dTrdxi__constxj(x, j, xN_flag) * d_PSI_T_dxj(x, i, k, xN_flag);
    double line4 =
      -(d2Trdxidxj(x, j, k, xN_flag) - 1 / Tr(x) * dTrdxi__constxj(x, k, xN_flag) * dTrdxi__constxj(x, j, xN_flag)) * PSI_T(x, i, xN_flag);
    return 1 / Tr(x) * (line1 + line2 + line3 + line4);
}

CoolPropDbl GERG2008ReducingFunction::Tr(const std::vector<CoolPropDbl>& x) const {
    return Yr(x, beta_T, gamma_T, T_c, Yc_T);
}
CoolPropDbl GERG2008ReducingFunction::dTr_dbetaT(const std::vector<CoolPropDbl>& x) const {
    return dYr_dbeta(x, beta_T, gamma_T, T_c, Yc_T);
}
CoolPropDbl GERG2008ReducingFunction::dTr_dgammaT(const std::vector<CoolPropDbl>& x) const {
    return dYr_dgamma(x, beta_T, gamma_T, T_c, Yc_T);
}
CoolPropDbl GERG2008ReducingFunction::d2Tr_dxidgammaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxidgamma(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
};
CoolPropDbl GERG2008ReducingFunction::d2Tr_dxidbetaT(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxidbeta(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
};

CoolPropDbl GERG2008ReducingFunction::dTrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return dYrdxi__constxj(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2Trdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxi2__constxj(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2Trdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) const {
    return d2Yrdxidxj(x, i, j, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d3Trdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                    x_N_dependency_flag xN_flag) const {
    return d3Yrdxidxjdxk(x, i, j, k, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::rhormolar(const std::vector<CoolPropDbl>& x) const {
    return 1 / Yr(x, beta_v, gamma_v, v_c, Yc_v);
}

CoolPropDbl GERG2008ReducingFunction::d2rhormolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    CoolPropDbl rhor = rhormolar(x);
    return -rhor * rhor * d2vrmolar_dxidgammaV(x, i, xN_flag)
           + 2 * POW3(rhor) * dvrmolardxi__constxj(x, i, xN_flag) * dYr_dgamma(x, beta_v, gamma_v, v_c, Yc_v);
}
CoolPropDbl GERG2008ReducingFunction::d2rhormolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    CoolPropDbl rhor = rhormolar(x);
    return -rhor * rhor * d2vrmolar_dxidbetaV(x, i, xN_flag)
           + 2 * POW3(rhor) * dvrmolardxi__constxj(x, i, xN_flag) * dYr_dbeta(x, beta_v, gamma_v, v_c, Yc_v);
}
CoolPropDbl GERG2008ReducingFunction::d2vrmolar_dxidgammaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxidgamma(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
};
CoolPropDbl GERG2008ReducingFunction::d2vrmolar_dxidbetaV(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxidbeta(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
};
CoolPropDbl GERG2008ReducingFunction::drhormolar_dgammaV(const std::vector<CoolPropDbl>& x) const {
    CoolPropDbl rhor = rhormolar(x);
    return -rhor * rhor * dYr_dgamma(x, beta_v, gamma_v, v_c, Yc_v);
}
CoolPropDbl GERG2008ReducingFunction::drhormolar_dbetaV(const std::vector<CoolPropDbl>& x) const {
    CoolPropDbl rhor = rhormolar(x);
    return -rhor * rhor * dYr_dbeta(x, beta_v, gamma_v, v_c, Yc_v);
}
CoolPropDbl GERG2008ReducingFunction::drhormolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    CoolPropDbl rhor = rhormolar(x);
    return -rhor * rhor * dvrmolardxi__constxj(x, i, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::dvrmolardxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return dYrdxi__constxj(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2vrmolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    return d2Yrdxi2__constxj(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2vrmolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                                      x_N_dependency_flag xN_flag) const {
    return d2Yrdxidxj(x, i, j, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d3vrmolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                         x_N_dependency_flag xN_flag) const {
    return d3Yrdxidxjdxk(x, i, j, k, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2rhormolardxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) const {
    CoolPropDbl rhor = this->rhormolar(x);
    CoolPropDbl dvrbardxi = this->dvrmolardxi__constxj(x, i, xN_flag);
    return 2 * pow(rhor, 3) * pow(dvrbardxi, 2) - pow(rhor, 2) * this->d2vrmolardxi2__constxj(x, i, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d2rhormolardxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                                        x_N_dependency_flag xN_flag) const {
    double rhor = this->rhormolar(x);
    double dvrbardxi = this->dvrmolardxi__constxj(x, i, xN_flag);
    double dvrbardxj = this->dvrmolardxi__constxj(x, j, xN_flag);
    return 2 * pow(rhor, 3) * dvrbardxi * dvrbardxj - pow(rhor, 2) * this->d2vrmolardxidxj(x, i, j, xN_flag);
}
CoolPropDbl GERG2008ReducingFunction::d3rhormolardxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                           x_N_dependency_flag xN_flag) const {
    double rhor = this->rhormolar(x);
    double line1 = -pow(rhor, 2) * d3vrmolardxidxjdxk(x, i, j, k, xN_flag);
    double line2 = 2 * pow(rhor, 3)
                   * (dvrmolardxi__constxj(x, k, xN_flag) * d2vrmolardxidxj(x, i, j, xN_flag)
                      + dvrmolardxi__constxj(x, j, xN_flag) * d2vrmolardxidxj(x, i, k, xN_flag)
                      + dvrmolardxi__constxj(x, i, xN_flag) * d2vrmolardxidxj(x, j, k, xN_flag));
    double line3 =
      -6 * pow(rhor, 4) * dvrmolardxi__constxj(x, i, xN_flag) * dvrmolardxi__constxj(x, j, xN_flag) * dvrmolardxi__constxj(x, k, xN_flag);
    return line1 + line2 + line3;
}

CoolPropDbl GERG2008ReducingFunction::Yr(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                                         const std::vector<CoolPropDbl>& Yc) const {
    CoolPropDbl Yr = 0;
    for (std::size_t i = 0; i < N; i++) {
        double xi = x[i];
        Yr += xi * xi * Yc[i];

        // The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
        if (i == N - 1) {
            break;
        }

        for (std::size_t j = i + 1; j < N; j++) {
            Yr += c_Y_ij(i, j, beta, gamma, Y_c_ij) * f_Y_ij(x, i, j, beta);
        }
    }
    return Yr;
}

CoolPropDbl GERG2008ReducingFunction::dYr_dgamma(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma,
                                                 const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc) const {
    CoolPropDbl dYr_dgamma = 0;
    for (std::size_t i = 0; i < N; i++) {
        // The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
        if (i == N - 1) {
            break;
        }
        for (std::size_t j = i + 1; j < N; j++) {
            dYr_dgamma += 2 * beta[i][j] * Y_c_ij[i][j] * f_Y_ij(x, i, j, beta);
        }
    }
    return dYr_dgamma;
}
CoolPropDbl GERG2008ReducingFunction::dYr_dbeta(const std::vector<CoolPropDbl>& x, const STLMatrix& beta, const STLMatrix& gamma,
                                                const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc) const {
    CoolPropDbl dYr_dbeta = 0;
    for (std::size_t i = 0; i < N; i++) {
        // The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
        if (i == N - 1) {
            break;
        }

        for (std::size_t j = i + 1; j < N; j++) {
            double xj = x[j], xi = x[i], beta_Y = beta[i][j], beta_Y_squared = beta_Y * beta_Y;
            if (std::abs(xi) < 10 * DBL_EPSILON && std::abs(xj) < 10 * DBL_EPSILON) {
                return 0;
            }
            double dfYij_dbeta = xi * xj * (-(xi + xj) * (2 * beta_Y * xi)) / pow(beta_Y_squared * xi + xj, 2);
            dYr_dbeta += c_Y_ij(i, j, beta, gamma, Y_c_ij) * dfYij_dbeta + f_Y_ij(x, i, j, beta) * 2 * gamma[i][j] * Y_c_ij[i][j];
        }
    }
    return dYr_dbeta;
}

CoolPropDbl GERG2008ReducingFunction::dYrdxi__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                                                      const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc,
                                                      x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        CoolPropDbl xi = x[i];
        CoolPropDbl dYr_dxi = 2 * xi * Yc[i];
        for (std::size_t k = 0; k < i; k++) {
            dYr_dxi += c_Y_ij(k, i, beta, gamma, Y_c_ij) * dfYkidxi__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N; k++) {
            dYr_dxi += c_Y_ij(i, k, beta, gamma, Y_c_ij) * dfYikdxi__constxk(x, i, k, beta);
        }
        return dYr_dxi;
    } else if (xN_flag == XN_DEPENDENT) {
        // Table S1 from Gernert, 2014, supplemental information
        if (i == N - 1) {
            return 0.0;
        }
        CoolPropDbl dYr_dxi = 2 * x[i] * Yc[i] - 2 * x[N - 1] * Yc[N - 1];
        for (std::size_t k = 0; k < i; k++) {
            dYr_dxi += c_Y_ij(k, i, beta, gamma, Y_c_ij) * dfYkidxi__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N - 1; k++) {
            dYr_dxi += c_Y_ij(i, k, beta, gamma, Y_c_ij) * dfYikdxi__constxk(x, i, k, beta);
        }
        double beta_Y_iN = beta[i][N - 1], xN = x[N - 1];
        dYr_dxi += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij)
                   * (xN * (x[i] + xN) / (pow(beta_Y_iN, 2) * x[i] + xN)
                      + (1 - beta_Y_iN * beta_Y_iN) * x[i] * xN * xN / POW2(beta_Y_iN * beta_Y_iN * x[i] + xN));
        for (std::size_t k = 0; k < N - 1; ++k) {
            double beta_Y_kN = beta[k][N - 1], xk = x[k], beta_Y_kN_squared = beta_Y_kN * beta_Y_kN;
            dYr_dxi +=
              c_Y_ij(k, N - 1, beta, gamma, Y_c_ij)
              * (-xk * (xk + xN) / (beta_Y_kN_squared * xk + xN) + (1 - beta_Y_kN_squared) * xN * xk * xk / POW2(beta_Y_kN_squared * xk + xN));
        }
        return dYr_dxi;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl GERG2008ReducingFunction::d2Yrdxidbeta(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                                                   const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        CoolPropDbl xi = x[i];
        CoolPropDbl deriv = 0;
        for (std::size_t k = 0; k < i; k++) {
            /*
            Sympy code:
            x_k,x_i,beta_Y = symbols('x_k,x_i,beta_Y')
            dfkidxi = x_k*(x_k+x_i)/(beta_Y**2*x_k+x_i) + x_k*x_i/(beta_Y**2*x_k+x_i)*(1-(x_k+x_i)/(beta_Y**2*x_k+x_i))
            simplify(diff(dfkidxi, beta_Y))
            */
            double xk = x[k], beta_Y = beta[k][i], beta_Y_squared = beta_Y * beta_Y;
            double d2fYkidxidbeta = 2 * beta_Y * pow(xk, 2) * (xi * (xi + xk * (-beta_Y_squared + 1) + xk) - (xi + xk) * (beta_Y_squared * xk + xi))
                                    / pow(beta_Y_squared * xk + xi, 3);
            deriv += c_Y_ij(k, i, beta, gamma, Y_c_ij) * d2fYkidxidbeta + dfYkidxi__constxk(x, k, i, beta) * 2 * gamma[k][i] * Y_c_ij[k][i];
        }
        for (std::size_t k = i + 1; k < N; k++) {
            /*
            x_k,x_i,beta_Y = symbols('x_k,x_i,beta_Y')
            dfikdxi = x_k*(x_i+x_k)/(beta_Y**2*x_i+x_k) + x_i*x_k/(beta_Y**2*x_i+x_k)*(1-beta_Y**2*(x_i+x_k)/(beta_Y**2*x_i+x_k))
            print(ccode(simplify(diff(dfikdxi, beta_Y))))
            */
            double xk = x[k], beta_Y = beta[i][k], beta_Y_squared = beta_Y * beta_Y;
            double d2fYikdxidbeta =
              2 * beta_Y * xi * xk
              * (xi * (-beta_Y_squared * xi + beta_Y_squared * (xi + xk) - xk) - xk * (xi + xk) - (xi + xk) * (beta_Y_squared * xi + xk))
              / pow(beta_Y_squared * xi + xk, 3);
            deriv += c_Y_ij(i, k, beta, gamma, Y_c_ij) * d2fYikdxidbeta + dfYikdxi__constxk(x, i, k, beta) * 2 * gamma[i][k] * Y_c_ij[i][k];
        }
        return deriv;
    } else if (xN_flag == XN_DEPENDENT) {
        throw NotImplementedError("Not yet implemented for xN_dependent");
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl GERG2008ReducingFunction::d2Yrdxidgamma(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta, const STLMatrix& gamma,
                                                    const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        CoolPropDbl deriv = 0;
        for (std::size_t k = 0; k < i; k++) {
            deriv += 2 * beta[k][i] * Y_c_ij[k][i] * dfYkidxi__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N; k++) {
            deriv += 2 * beta[i][k] * Y_c_ij[i][k] * dfYikdxi__constxk(x, i, k, beta);
        }
        return deriv;
    } else if (xN_flag == XN_DEPENDENT) {
        // Table S1 from Gernert, 2014, supplemental information
        if (i == N - 1) {
            return 0.0;
        }
        CoolPropDbl deriv = 0;
        for (std::size_t k = 0; k < i; k++) {
            deriv += 2 * beta[k][i] * Y_c_ij[k][i] * dfYkidxi__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N - 1; k++) {
            deriv += 2 * beta[i][k] * Y_c_ij[i][k] * dfYikdxi__constxk(x, i, k, beta);
        }
        double beta_Y_iN = beta[i][N - 1], xN = x[N - 1];
        deriv += 2 * beta[i][N - 1] * Y_c_ij[i][N - 1]
                 * (xN * (x[i] + xN) / (pow(beta_Y_iN, 2) * x[i] + xN)
                    + (1 - beta_Y_iN * beta_Y_iN) * x[i] * xN * xN / POW2(beta_Y_iN * beta_Y_iN * x[i] + xN));
        for (std::size_t k = 0; k < N - 1; ++k) {
            double beta_Y_kN = beta[k][N - 1], xk = x[k], beta_Y_kN_squared = beta_Y_kN * beta_Y_kN;
            deriv += 2 * beta[k][N - 1] * Y_c_ij[k][N - 1]
                     * (-xk * (xk + xN) / (beta_Y_kN_squared * xk + xN) + (1 - beta_Y_kN_squared) * xN * xk * xk / POW2(beta_Y_kN_squared * xk + xN));
        }
        return deriv;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}
CoolPropDbl GERG2008ReducingFunction::d2Yrdxi2__constxj(const std::vector<CoolPropDbl>& x, std::size_t i, const STLMatrix& beta,
                                                        const STLMatrix& gamma, const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc,
                                                        x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        CoolPropDbl d2Yr_dxi2 = 2 * Yc[i];
        for (std::size_t k = 0; k < i; k++) {
            d2Yr_dxi2 += c_Y_ij(k, i, beta, gamma, Y_c_ij) * d2fYkidxi2__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N; k++) {
            d2Yr_dxi2 += c_Y_ij(i, k, beta, gamma, Y_c_ij) * d2fYikdxi2__constxk(x, i, k, beta);
        }
        return d2Yr_dxi2;
    } else if (xN_flag == XN_DEPENDENT) {
        // Table S1 from Gernert, 2014, supplemental information
        if (i == N - 1) {
            return 0.0;
        }
        CoolPropDbl d2Yr_dxi2 = 2 * Yc[i] + 2 * Yc[N - 1];
        for (std::size_t k = 0; k < i; k++) {
            d2Yr_dxi2 += c_Y_ij(k, i, beta, gamma, Y_c_ij) * d2fYkidxi2__constxk(x, k, i, beta);
        }
        for (std::size_t k = i + 1; k < N - 1; k++) {
            d2Yr_dxi2 += c_Y_ij(i, k, beta, gamma, Y_c_ij) * d2fYikdxi2__constxk(x, i, k, beta);
        }
        double beta_Y_iN = beta[i][N - 1], xN = x[N - 1];
        d2Yr_dxi2 += 2 * c_Y_ij(i, N - 1, beta, gamma, Y_c_ij)
                     * (-(x[i] + xN) / (pow(beta_Y_iN, 2) * x[i] + xN)
                        + (1 - beta_Y_iN * beta_Y_iN)
                            * (xN * xN / pow(beta_Y_iN * beta_Y_iN * x[i] + xN, 2)
                               + ((1 - beta_Y_iN * beta_Y_iN) * x[i] * xN * xN - beta_Y_iN * beta_Y_iN * x[i] * x[i] * xN)
                                   / pow(beta_Y_iN * beta_Y_iN * x[i] + xN, 3)));
        for (std::size_t k = 0; k < N - 1; ++k) {
            double beta_Y_kN = beta[k][N - 1], xk = x[k], beta_Y_kN_squared = beta_Y_kN * beta_Y_kN;
            d2Yr_dxi2 += 2 * c_Y_ij(k, N - 1, beta, gamma, Y_c_ij) * xk * xk * (1 - beta_Y_kN_squared) / pow(beta_Y_kN_squared * xk + xN, 2)
                         * (xN / (beta_Y_kN_squared * xk + xN) - 1);
        }
        return d2Yr_dxi2;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}

CoolPropDbl GERG2008ReducingFunction::d2Yrdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta,
                                                 const STLMatrix& gamma, const STLMatrix& Y_c_ij, const std::vector<CoolPropDbl>& Yc,
                                                 x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        if (i == j) {
            return d2Yrdxi2__constxj(x, i, beta, gamma, Y_c_ij, Yc, xN_flag);
        } else {
            // See Table B9 from Kunz Wagner 2012 (GERG 2008)
            return c_Y_ij(i, j, beta, gamma, Y_c_ij) * d2fYijdxidxj(x, i, j, beta);
        }
    } else if (xN_flag == XN_DEPENDENT) {
        // Table S1 from Gernert, 2014, supplemental information
        if (j == N - 1 || i == N - 1) {
            return 0.0;
        }
        if (i == j) {
            return d2Yrdxi2__constxj(x, i, beta, gamma, Y_c_ij, Yc, xN_flag);
        }
        CoolPropDbl d2Yr_dxidxj = 2 * Yc[N - 1];
        d2Yr_dxidxj += c_Y_ij(i, j, beta, gamma, Y_c_ij) * d2fYijdxidxj(x, i, j, beta);

        for (std::size_t k = 0; k < N - 1; k++) {
            d2Yr_dxidxj += c_Y_ij(k, N - 1, beta, gamma, Y_c_ij) * d2fYkidxi2__constxk(x, k, N - 1, beta);
        }
        d2Yr_dxidxj -= c_Y_ij(i, N - 1, beta, gamma, Y_c_ij) * d2fYijdxidxj(x, i, N - 1, beta);
        d2Yr_dxidxj -= c_Y_ij(j, N - 1, beta, gamma, Y_c_ij) * d2fYijdxidxj(x, j, N - 1, beta);
        return d2Yr_dxidxj;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}

CoolPropDbl GERG2008ReducingFunction::d3Yrdxidxjdxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                                    const STLMatrix& beta, const STLMatrix& gamma, const STLMatrix& Y_c_ij,
                                                    const std::vector<CoolPropDbl>& Yc, x_N_dependency_flag xN_flag) const {
    if (xN_flag == XN_INDEPENDENT) {
        if (i != j && j != k && k != i) {
            return 0;
        } else if (k == i && i != j) {
            return c_Y_ij(i, j, beta, gamma, Y_c_ij) * d3fYijdxi2dxj(x, i, j, beta);
        } else if (k == j && i != j) {
            return c_Y_ij(i, j, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, i, j, beta);
        } else if (i == j && i != k) {
            return c_Y_ij(i, k, beta, gamma, Y_c_ij) * d3fYijdxi2dxj(x, i, k, beta);
        } else {
            CoolPropDbl d3Yr_dxi3 = 0;
            for (std::size_t m = 0; m < i; m++) {
                d3Yr_dxi3 += c_Y_ij(m, i, beta, gamma, Y_c_ij) * d3fYkidxi3__constxk(x, m, i, beta);
            }
            for (std::size_t m = i + 1; m < N; m++) {
                d3Yr_dxi3 += c_Y_ij(i, m, beta, gamma, Y_c_ij) * d3fYikdxi3__constxk(x, i, m, beta);
            }
            return d3Yr_dxi3;
        }
    } else if (xN_flag == XN_DEPENDENT) {
        CoolPropDbl summer = 0;
        // Needed for all third partials
        for (std::size_t m = 0; m < N - 1; m++) {
            summer -= c_Y_ij(m, N - 1, beta, gamma, Y_c_ij) * d3fYkidxi3__constxk(x, m, N - 1, beta);
        }
        if (i != j && j != k && k != i) {
            summer += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, i, N - 1, beta);
            summer += c_Y_ij(j, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, j, N - 1, beta);
            summer += c_Y_ij(k, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, k, N - 1, beta);
        } else if (k == i && i != j) {  // two i, one j
            summer += c_Y_ij(i, j, beta, gamma, Y_c_ij) * d3fYijdxi2dxj(x, i, j, beta);
            summer += c_Y_ij(j, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, j, N - 1, beta);
            summer += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij) * (2 * d3fYijdxidxj2(x, i, N - 1, beta) - d3fYijdxi2dxj(x, i, N - 1, beta));
        } else if (k == j && i != j) {  // two j, one i
            summer += c_Y_ij(i, j, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, i, j, beta);
            summer += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, i, N - 1, beta);
            summer += c_Y_ij(j, N - 1, beta, gamma, Y_c_ij) * (2 * d3fYijdxidxj2(x, j, N - 1, beta) - d3fYijdxi2dxj(x, j, N - 1, beta));
        } else if (i == j && i != k) {  // two i, one k
            summer += c_Y_ij(i, k, beta, gamma, Y_c_ij) * d3fYijdxi2dxj(x, i, k, beta);
            summer += c_Y_ij(k, N - 1, beta, gamma, Y_c_ij) * d3fYijdxidxj2(x, k, N - 1, beta);
            summer += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij) * (2 * d3fYijdxidxj2(x, i, N - 1, beta) - d3fYijdxi2dxj(x, i, N - 1, beta));
        } else {
            for (std::size_t m = 0; m < i; m++) {
                summer += c_Y_ij(m, i, beta, gamma, Y_c_ij) * d3fYkidxi3__constxk(x, m, i, beta);
            }
            for (std::size_t m = i + 1; m < N - 1; m++) {
                summer += c_Y_ij(i, m, beta, gamma, Y_c_ij) * d3fYikdxi3__constxk(x, i, m, beta);
            }
            summer += c_Y_ij(i, N - 1, beta, gamma, Y_c_ij)
                      * (3 * d3fYijdxidxj2(x, i, N - 1, beta) - 3 * d3fYijdxi2dxj(x, i, N - 1, beta) + d3fYikdxi3__constxk(x, i, N - 1, beta));
        }
        return summer;
    } else {
        throw ValueError(format("xN dependency flag invalid"));
    }
}

CoolPropDbl GERG2008ReducingFunction::dfYkidxi__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i,
                                                        const STLMatrix& beta) const {
    double xk = x[k], xi = x[i], beta_Y = beta[k][i], beta_Y_squared = beta_Y * beta_Y;
    return xk * (xk + xi) / (beta_Y_squared * xk + xi) + xk * xi / (beta_Y_squared * xk + xi) * (1 - (xk + xi) / (beta_Y_squared * xk + xi));
}
CoolPropDbl GERG2008ReducingFunction::dfYikdxi__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k,
                                                        const STLMatrix& beta) const {
    double xk = x[k], xi = x[i], beta_Y = beta[i][k], beta_Y_squared = beta_Y * beta_Y;
    return xk * (xi + xk) / (beta_Y_squared * xi + xk)
           + xi * xk / (beta_Y_squared * xi + xk) * (1 - beta_Y_squared * (xi + xk) / (beta_Y_squared * xi + xk));
}
const CoolPropDbl GERG2008ReducingFunction::c_Y_ij(const std::size_t i, const std::size_t j, const STLMatrix& beta, const STLMatrix& gamma,
                                                   const STLMatrix& Y_c) const {
    return 2 * beta[i][j] * gamma[i][j] * Y_c[i][j];
}
CoolPropDbl GERG2008ReducingFunction::f_Y_ij(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const {
    double xi = x[i], xj = x[j], beta_Y = beta[i][j];
    return xi * xj * (xi + xj) / (beta_Y * beta_Y * xi + xj);
}
CoolPropDbl GERG2008ReducingFunction::d2fYikdxi2__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k,
                                                          const STLMatrix& beta) const {
    double xi = x[i], xk = x[k], beta_Y = beta[i][k];
    return 1 / (beta_Y * beta_Y * xi + xk) * (1 - beta_Y * beta_Y * (xi + xk) / (beta_Y * beta_Y * xi + xk))
           * (2 * xk - xi * xk * 2 * beta_Y * beta_Y / (beta_Y * beta_Y * xi + xk));
}
CoolPropDbl GERG2008ReducingFunction::d2fYkidxi2__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i,
                                                          const STLMatrix& beta) const {
    double xi = x[i], xk = x[k], beta_Y = beta[k][i];
    return 1 / (beta_Y * beta_Y * xk + xi) * (1 - (xk + xi) / (beta_Y * beta_Y * xk + xi)) * (2 * xk - xk * xi * 2 / (beta_Y * beta_Y * xk + xi));
}
CoolPropDbl GERG2008ReducingFunction::d2fYijdxidxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const {
    double xi = x[i], xj = x[j], beta_Y = beta[i][j], beta_Y2 = beta_Y * beta_Y;
    return (xi + xj) / (beta_Y2 * xi + xj) + xj / (beta_Y2 * xi + xj) * (1 - (xi + xj) / (beta_Y2 * xi + xj))
           + xi / (beta_Y2 * xi + xj) * (1 - beta_Y2 * (xi + xj) / (beta_Y2 * xi + xj))
           - xi * xj / pow(beta_Y2 * xi + xj, 2) * (1 + beta_Y2 - 2 * beta_Y2 * (xi + xj) / (beta_Y2 * xi + xj));
}
CoolPropDbl GERG2008ReducingFunction::d3fYijdxi2dxj(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const {
    double x_i = x[i], x_j = x[j], beta_Y = beta[i][j], beta_Y2 = beta_Y * beta_Y;
    double den = pow(beta_Y, 8) * pow(x_i, 4) + 4 * pow(beta_Y, 6) * pow(x_i, 3) * x_j + 6 * pow(beta_Y, 4) * pow(x_i * x_j, 2)
                 + 4 * beta_Y2 * x_i * pow(x_j, 3) + pow(x_j, 4);
    return -6 * beta_Y2 * x_i * x_j * x_j * (beta_Y2 - 1) / den;
}
CoolPropDbl GERG2008ReducingFunction::d3fYijdxidxj2(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, const STLMatrix& beta) const {
    double x_i = x[i], x_j = x[j], beta_Y = beta[i][j], beta_Y2 = beta_Y * beta_Y;
    double den = pow(beta_Y, 8) * pow(x_i, 4) + 4 * pow(beta_Y, 6) * pow(x_i, 3) * x_j + 6 * pow(beta_Y, 4) * pow(x_i * x_j, 2)
                 + 4 * beta_Y2 * x_i * pow(x_j, 3) + pow(x_j, 4);
    return 6 * beta_Y2 * x_i * x_i * x_j * (beta_Y2 - 1) / den;
}
CoolPropDbl GERG2008ReducingFunction::d3fYikdxi3__constxk(const std::vector<CoolPropDbl>& x, std::size_t i, std::size_t k,
                                                          const STLMatrix& beta) const {
    double x_i = x[i], x_k = x[k], beta_Y = beta[i][k], beta_Y2 = beta_Y * beta_Y;
    double den = pow(beta_Y, 8) * pow(x_i, 4) + 4 * pow(beta_Y, 6) * pow(x_i, 3) * x_k + 6 * pow(beta_Y, 4) * pow(x_i * x_k, 2)
                 + 4 * beta_Y2 * x_i * pow(x_k, 3) + pow(x_k, 4);
    return 6 * beta_Y2 * x_k * x_k * x_k * (beta_Y2 - 1) / den;
}
CoolPropDbl GERG2008ReducingFunction::d3fYkidxi3__constxk(const std::vector<CoolPropDbl>& x, std::size_t k, std::size_t i,
                                                          const STLMatrix& beta) const {
    double x_i = x[i], x_k = x[k], beta_Y = beta[k][i], beta_Y2 = beta_Y * beta_Y;
    return 6 * beta_Y2 * x_k * x_k * x_k * (1 - beta_Y2) / pow(beta_Y2 * x_k + x_i, 4);
}

} /* namespace CoolProp */
