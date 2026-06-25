#include "GeneralizedCubic.h"
#include "CoolProp/numerics/numerics.h"
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

double BasicMathiasCopemanAlphaFunction::term(double tau, std::size_t itau) {

    // If we are not using the full Mathias-Copeman formulation for a_ii,
    // we just use the simple results from the supplemental information because
    // they are much more computationally efficient

    // All derivatives have a common bracketed term, so we factor it out
    // and calculate it here
    double B = 1 + m * (1 - sqrt_Tr_Tci * sqrt(1 / tau));

    switch (itau) {
        case 0:
            return a0 * B * B;
        case 1:
            return a0 * m * B / pow(tau, 3.0 / 2.0) * sqrt_Tr_Tci;
        case 2:
            return a0 * m / 2.0 * (m / pow(tau, 3) * Tr_over_Tci - 3 * B / pow(tau, 5.0 / 2.0) * sqrt_Tr_Tci);
        case 3:
            return (3.0 / 4.0) * a0 * m * (-3.0 * m / pow(tau, 4) * Tr_over_Tci + 5 * B / pow(tau, 7.0 / 2.0) * sqrt_Tr_Tci);
        case 4:
            return (3.0 / 8.0) * a0 * m * (29.0 * m / pow(tau, 5) * Tr_over_Tci - 35 * B / pow(tau, 9.0 / 2.0) * sqrt_Tr_Tci);
        default:
            throw -1;
    }
}

double MathiasCopemanAlphaFunction::term(double tau, std::size_t itau) {
    // Here we are using the full Mathias-Copeman formulation, introducing
    // some additional computational effort, so we only evaluate the parameters that
    // we actually need to evaluate, otherwise we just set their value to zero
    // See info on the conditional (ternary) operator : http://www.cplusplus.com/articles/1AUq5Di1/
    // Furthermore, this should help with branch prediction
    double Di = 1 - sqrt_Tr_Tci / sqrt(tau);
    double dDi_dtau = (itau >= 1) ? (1.0 / 2.0) * sqrt_Tr_Tci / (pow(tau, 1.5)) : 0;
    double d2Di_dtau2 = (itau >= 2) ? -(3.0 / 4.0) * sqrt_Tr_Tci / (pow(tau, 2.5)) : 0;
    double d3Di_dtau3 = (itau >= 3) ? (15.0 / 8.0) * sqrt_Tr_Tci / (pow(tau, 3.5)) : 0;
    double d4Di_dtau4 = (itau >= 4) ? -(105.0 / 16.0) * sqrt_Tr_Tci / (pow(tau, 4.5)) : 0;

    double Bi = 1, dBi_dtau = 0, d2Bi_dtau2 = 0, d3Bi_dtau3 = 0, d4Bi_dtau4 = 0;
    for (int n = 1; n <= 3; ++n) {
        Bi += c[n - 1] * pow(Di, n);
        dBi_dtau += (itau < 1) ? 0 : (n * c[n - 1] * pow(Di, n - 1) * dDi_dtau);
        d2Bi_dtau2 += (itau < 2) ? 0 : n * c[n - 1] * ((n - 1) * pow(dDi_dtau, 2) + Di * d2Di_dtau2) * pow(Di, n - 2);
        d3Bi_dtau3 += (itau < 3)
                        ? 0
                        : n * c[n - 1] * (3 * (n - 1) * Di * dDi_dtau * d2Di_dtau2 + (n * n - 3 * n + 2) * pow(dDi_dtau, 3) + pow(Di, 2) * d3Di_dtau3)
                            * pow(Di, n - 3);
        d4Bi_dtau4 +=
          (itau < 4)
            ? 0
            : n * c[n - 1]
                * (6 * (n * n - 3 * n + 2) * Di * pow(dDi_dtau, 2) * d2Di_dtau2 + (n * n * n - 6 * n * n + 11 * n - 6) * pow(dDi_dtau, 4)
                   + (4 * n * dDi_dtau * d3Di_dtau3 + 3 * n * pow(d2Di_dtau2, 2) - 4 * dDi_dtau * d3Di_dtau3 - 3 * pow(d2Di_dtau2, 2)) * pow(Di, 2)
                   + pow(Di, 3) * d4Di_dtau4)
                * pow(Di, n - 4);
    }
    switch (itau) {
        case 0:
            return a0 * Bi * Bi;
        case 1:
            return 2 * a0 * Bi * dBi_dtau;
        case 2:
            return 2 * a0 * (Bi * d2Bi_dtau2 + dBi_dtau * dBi_dtau);
        case 3:
            return 2 * a0 * (Bi * d3Bi_dtau3 + 3 * dBi_dtau * d2Bi_dtau2);
        case 4:
            return 2 * a0 * (Bi * d4Bi_dtau4 + 4 * dBi_dtau * d3Bi_dtau3 + 3 * pow(d2Bi_dtau2, 2));
        default:
            throw -1;
    }
}

double TwuAlphaFunction::term(double tau, std::size_t itau) {
    // Here we are using the Twu formulation, introducing
    // some additional computational effort, so we only evaluate the parameters that
    // we actually need to evaluate, otherwise we just set their value to zero
    // See info on the conditional (ternary) operator : http://www.cplusplus.com/articles/1AUq5Di1/
    // Furthermore, this should help with branch prediction

    const double L = c[0], M = c[1], N = c[2];
    double A = pow(Tr_over_Tci / tau, M * N);
    double B1 = (itau < 1) ? 0 : N / tau * (L * M * A - M + 1);
    double dB1_dtau = (itau < 2) ? 0 : N / powInt(tau, 2) * (-L * M * M * N * A - L * M * A + M - 1);
    double d2B1_dtau2 = (itau < 3) ? 0 : N / powInt(tau, 3) * (L * M * M * M * N * N * A + 3 * L * M * M * N * A + 2 * L * M * A - 2 * M + 2);
    double d3B1_dtau3 =
      (itau < 4) ? 0
                 : -N / powInt(tau, 4)
                     * (L * powInt(M, 4) * powInt(N, 3) * A + 6 * L * M * M * M * N * N * A + 11 * L * M * M * N * A + 6 * L * M * A - 6 * M + 6);

    double dam_dtau = NAN, d2am_dtau2 = NAN, d3am_dtau3 = NAN, d4am_dtau4 = NAN;
    double am = a0 * pow(Tr_over_Tci / tau, N * (M - 1)) * exp(L * (1 - A));

    if (itau == 0) {
        return am;
    } else {
        // Calculate terms as needed
        dam_dtau = am * B1;
        d2am_dtau2 = (itau < 2) ? 0 : B1 * dam_dtau + am * dB1_dtau;
        d3am_dtau3 = (itau < 3) ? 0 : B1 * d2am_dtau2 + am * d2B1_dtau2 + 2 * dB1_dtau * dam_dtau;
        d4am_dtau4 = (itau < 4) ? 0 : B1 * d3am_dtau3 + am * d3B1_dtau3 + 3 * dB1_dtau * d2am_dtau2 + 3 * d2B1_dtau2 * dam_dtau;
    }
    switch (itau) {
        case 1:
            return dam_dtau;
        case 2:
            return d2am_dtau2;
        case 3:
            return d3am_dtau3;
        case 4:
            return d4am_dtau4;
        default:
            throw -1;
    }
}

void BasicMathiasCopemanAlphaFunction::calc_all_terms(double tau, std::array<double, 5>& terms) {
    // Compute B and shared powers of tau once, then fill all 5 derivatives
    const double sq = sqrt_Tr_Tci / sqrt(tau);  // sqrt(Tr/Tci / tau)
    const double B = 1.0 + m * (1.0 - sq);
    const double t15 = sqrt_Tr_Tci / (tau * sqrt(tau));  // sqrt_Tr_Tci / tau^1.5
    const double t25 = t15 / tau;                        // sqrt_Tr_Tci / tau^2.5
    const double t35 = t25 / tau;                        // sqrt_Tr_Tci / tau^3.5
    const double t45 = t35 / tau;                        // sqrt_Tr_Tci / tau^4.5
    const double t3 = Tr_over_Tci / (tau * tau * tau);
    const double t4 = t3 / tau;
    const double t5 = t4 / tau;
    terms[0] = a0 * B * B;
    terms[1] = a0 * m * B * t15;
    terms[2] = a0 * m * 0.5 * (m * t3 - 3.0 * B * t25);
    terms[3] = (3.0 / 4.0) * a0 * m * (-3.0 * m * t4 + 5.0 * B * t35);
    terms[4] = (3.0 / 8.0) * a0 * m * (29.0 * m * t5 - 35.0 * B * t45);
}

void TwuAlphaFunction::calc_all_terms(double tau, std::array<double, 5>& terms) {
    const double L = c[0], M = c[1], N = c[2];
    // Compute the two pow() and one exp() calls exactly once
    const double A = pow(Tr_over_Tci / tau, M * N);
    const double am = a0 * pow(Tr_over_Tci / tau, N * (M - 1)) * exp(L * (1.0 - A));
    terms[0] = am;
    const double B1 = N / tau * (L * M * A - M + 1.0);
    const double dam = am * B1;
    terms[1] = dam;
    const double dB1 = N / powInt(tau, 2) * (-L * M * M * N * A - L * M * A + M - 1.0);
    const double d2am = B1 * dam + am * dB1;
    terms[2] = d2am;
    const double d2B1 = N / powInt(tau, 3) * (L * M * M * M * N * N * A + 3.0 * L * M * M * N * A + 2.0 * L * M * A - 2.0 * M + 2.0);
    const double d3am = B1 * d2am + am * d2B1 + 2.0 * dB1 * dam;
    terms[3] = d3am;
    const double d3B1 =
      -N / powInt(tau, 4)
      * (L * powInt(M, 4) * powInt(N, 3) * A + 6.0 * L * M * M * M * N * N * A + 11.0 * L * M * M * N * A + 6.0 * L * M * A - 6.0 * M + 6.0);
    terms[4] = B1 * d3am + am * d3B1 + 3.0 * dB1 * d2am + 3.0 * d2B1 * dam;
}

void MathiasCopemanAlphaFunction::calc_all_terms(double tau, std::array<double, 5>& terms) {
    // Compute Di = 1 - sqrt(Tr/Tci / tau) and its tau-derivatives once, reusing
    // a common power chain.  All five derivatives share these values.
    const double sq = sqrt_Tr_Tci / sqrt(tau);  // sqrt(Tr_over_Tci / tau)
    const double Di = 1.0 - sq;
    const double p = sq / tau;  // sqrt_Tr_Tci / tau^1.5
    const double d1Di = 0.5 * p;
    const double d2Di = -0.75 * p / tau;
    const double d3Di = 1.875 * p / (tau * tau);
    const double d4Di = -6.5625 * p / (tau * tau * tau);

    const double Di2 = Di * Di;
    const double Di3 = Di2 * Di;
    const double c0 = c[0], c1 = c[1], c2 = c[2];

    // Expand the loop over n=1,2,3 analytically.
    // Many polynomial factors in (n^2-3n+2) and (n^3-6n^2+11n-6) vanish for small n.
    const double Bi = 1.0 + c0 * Di + c1 * Di2 + c2 * Di3;
    const double dBi = d1Di * (c0 + 2.0 * c1 * Di + 3.0 * c2 * Di2);
    const double d2Bi = c0 * d2Di + 2.0 * c1 * (d1Di * d1Di + Di * d2Di) + 3.0 * c2 * Di * (2.0 * d1Di * d1Di + Di * d2Di);
    const double d3Bi =
      c0 * d3Di + 2.0 * c1 * (3.0 * d1Di * d2Di + Di * d3Di) + 3.0 * c2 * (6.0 * Di * d1Di * d2Di + 2.0 * d1Di * d1Di * d1Di + Di2 * d3Di);
    const double d4Bi = c0 * d4Di + 2.0 * c1 * (4.0 * d1Di * d3Di + 3.0 * d2Di * d2Di + Di * d4Di)
                        + 3.0 * c2 * (12.0 * d1Di * d1Di * d2Di + Di * (8.0 * d1Di * d3Di + 6.0 * d2Di * d2Di) + Di2 * d4Di);

    terms[0] = a0 * Bi * Bi;
    terms[1] = 2.0 * a0 * Bi * dBi;
    terms[2] = 2.0 * a0 * (Bi * d2Bi + dBi * dBi);
    terms[3] = 2.0 * a0 * (Bi * d3Bi + 3.0 * dBi * d2Bi);
    terms[4] = 2.0 * a0 * (Bi * d4Bi + 4.0 * dBi * d3Bi + 3.0 * d2Bi * d2Bi);
}

AbstractCubic::AbstractCubic(const std::vector<double>& Tc, std::vector<double> pc, std::vector<double> acentric, double R_u, double Delta_1,
                             double Delta_2, const std::vector<double>& C1, const std::vector<double>& C2, const std::vector<double>& C3)
  : T_r(1.0),
    rho_r(1.0),
    Tc(Tc),
    pc(std::move(pc)),
    acentric(std::move(acentric)),
    R_u(R_u),
    Delta_1(Delta_1),
    Delta_2(Delta_2),
    N(static_cast<int>(Tc.size())),
    cm(0.) {

    k.resize(N, std::vector<double>(N, 0));

    alpha.resize(N);
    m_tau_cache = std::numeric_limits<double>::quiet_NaN();
};

void AbstractCubic::set_alpha(const std::vector<double>& C1, const std::vector<double>& C2, const std::vector<double>& C3) {
    /// Resize the vector of alpha functions
    alpha.resize(Tc.size());
    m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    /// If no Mathias-Copeman coefficients are passed in (all empty vectors), use the predictive scheme for m_ii
    if (C1.empty() && C2.empty() && C3.empty()) {
        for (std::size_t i = 0; i < Tc.size(); ++i) {
            alpha[i] = std::make_shared<BasicMathiasCopemanAlphaFunction>(a0_ii(i), m_ii(i), T_r / Tc[i]);
        }
    } else {
        /// Use the Mathias-Copeman constants passed in to initialize Mathias-Copeman alpha functions
        for (std::size_t i = 0; i < Tc.size(); ++i) {
            alpha[i] = std::make_shared<MathiasCopemanAlphaFunction>(a0_ii(i), C1[i], C2[i], C3[i], T_r / Tc[i]);
        }
    }
}

double AbstractCubic::am_term(double tau, const std::vector<double>& x, std::size_t itau) {
    // am_term and its composition-derivatives are the only callers of aij_term/u_term, so we
    // validate/populate the aii cache for this tau once here.  The inner u_term reads then go
    // straight to m_aii_cache instead of re-validating the cache on every aii lookup, which was
    // the dominant cost (~38% of self-time) in the pure-fluid hot path.
    _ensure_aii_cache(tau);
    double summer = 0;
    for (int i = N - 1; i >= 0; --i) {
        for (int j = N - 1; j >= 0; --j) {
            summer += x[i] * x[j] * aij_term(tau, i, j, itau);
        }
    }
    return summer;
}
double AbstractCubic::d_am_term_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
    _ensure_aii_cache(tau);  // see am_term: guarantees u_term can read m_aii_cache directly
    if (xN_independent) {
        double summer = 0;
        for (int j = N - 1; j >= 0; --j) {
            summer += x[j] * aij_term(tau, i, j, itau);
        }
        return 2 * summer;
    } else {
        double summer = 0;
        for (int k = N - 2; k >= 0; --k) {
            summer += x[k] * (aij_term(tau, i, k, itau) - aij_term(tau, k, N - 1, itau));
        }
        return 2 * (summer + x[N - 1] * (aij_term(tau, N - 1, i, itau) - aij_term(tau, N - 1, N - 1, itau)));
    }
}
double AbstractCubic::d2_am_term_dxidxj(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j,
                                        bool xN_independent) {
    _ensure_aii_cache(tau);  // see am_term: guarantees u_term can read m_aii_cache directly
    if (xN_independent) {
        return 2 * aij_term(tau, i, j, itau);
    } else {
        return 2 * (aij_term(tau, i, j, itau) - aij_term(tau, j, N - 1, itau) - aij_term(tau, N - 1, i, itau) + aij_term(tau, N - 1, N - 1, itau));
    }
}

double AbstractCubic::d3_am_term_dxidxjdxk(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k,
                                           bool xN_independent) {
    return 0;
}

double AbstractCubic::bm_term(const std::vector<double>& x) {
    double summer = 0;
    for (int i = N - 1; i >= 0; --i) {
        summer += x[i] * b0_ii(i);
    }
    return summer;
}
double AbstractCubic::d_bm_term_dxi(const std::vector<double>& x, std::size_t i, bool xN_independent) {
    if (xN_independent) {
        return b0_ii(i);
    } else {
        return b0_ii(i) - b0_ii(N - 1);
    }
}
double AbstractCubic::d2_bm_term_dxidxj(const std::vector<double>& x, std::size_t i, std::size_t j, bool xN_independent) {
    return 0;
}
double AbstractCubic::d3_bm_term_dxidxjdxk(const std::vector<double>& x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent) {
    return 0;
}

double AbstractCubic::cm_term() {
    return cm;
}

double AbstractCubic::aii_term(double tau, std::size_t i, std::size_t itau) {
    if (itau > 4) {
        // m_aii_cache[i] only stores derivatives for itau = 0..4 (see calc_all_terms()); higher
        // orders are unsupported, matching the contract of AbstractCubicAlphaFunction::term().
        throw -1;
    }
    _ensure_aii_cache(tau);
    return m_aii_cache[i][itau];
}
double AbstractCubic::u_term(double tau, std::size_t i, std::size_t j, std::size_t itau) {
    // Read the aii derivatives straight from the cache.  The only callers (the am_term family)
    // populate it for this tau first; the assert pins that invariant in debug builds so a future
    // caller that forgets to do so fails loudly instead of silently returning stale values.
    assert(std::memcmp(&tau, &m_tau_cache, sizeof(double)) == 0 && "u_term: aii cache not populated for this tau");
    (void)tau;
    const std::array<double, 5>& ai = m_aii_cache[i];
    const std::array<double, 5>& aj = m_aii_cache[j];
    const double aii = ai[0], ajj = aj[0];
    switch (itau) {
        case 0:
            return aii * ajj;
        case 1:
            return aii * aj[1] + ajj * ai[1];
        case 2:
            return (aii * aj[2] + 2 * ai[1] * aj[1] + ajj * ai[2]);
        case 3:
            return (aii * aj[3] + 3 * ai[1] * aj[2] + 3 * ai[2] * aj[1] + ajj * ai[3]);
        case 4:
            return (aii * aj[4] + 4 * ai[1] * aj[3] + 6 * ai[2] * aj[2] + 4 * ai[3] * aj[1] + ajj * ai[4]);
        default:
            throw -1;
    }
}
double AbstractCubic::aij_term(double tau, std::size_t i, std::size_t j, std::size_t itau) {
    double u = u_term(tau, i, j, 0);

    switch (itau) {
        case 0:
            return (1 - k[i][j]) * sqrt(u);
        case 1:
            return (1 - k[i][j]) / (2.0 * sqrt(u)) * u_term(tau, i, j, 1);
        case 2:
            return (1 - k[i][j]) / (4.0 * pow(u, 3.0 / 2.0)) * (2 * u * u_term(tau, i, j, 2) - pow(u_term(tau, i, j, 1), 2));
        case 3:
            return (1 - k[i][j]) / (8.0 * pow(u, 5.0 / 2.0))
                   * (4 * pow(u, 2) * u_term(tau, i, j, 3) - 6 * u * u_term(tau, i, j, 1) * u_term(tau, i, j, 2) + 3 * pow(u_term(tau, i, j, 1), 3));
        case 4:
            return (1 - k[i][j]) / (16.0 * pow(u, 7.0 / 2.0))
                   * (-4 * pow(u, 2) * (4 * u_term(tau, i, j, 1) * u_term(tau, i, j, 3) + 3 * pow(u_term(tau, i, j, 2), 2))
                      + 8 * pow(u, 3) * u_term(tau, i, j, 4) + 36 * u * pow(u_term(tau, i, j, 1), 2) * u_term(tau, i, j, 2)
                      - 15 * pow(u_term(tau, i, j, 1), 4));
        default:
            throw -1;
    }
}
double AbstractCubic::psi_minus(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta) {
    if (itau > 0) return 0.0;
    double bmc = bm_term(x) - cm_term();  // appears only in the form (b-c) in the equations
    double bracket = 1 - bmc * delta * rho_r;

    switch (idelta) {
        case 0:
            return -log(bracket);
        case 1:
            return bmc * rho_r / bracket;
        case 2:
            return pow(bmc * rho_r / bracket, 2);
        case 3:
            return 2 * pow(bmc * rho_r / bracket, 3);
        case 4:
            return 6 * pow(bmc * rho_r / bracket, 4);
        default:
            throw -1;
    }
}
double AbstractCubic::d_psi_minus_dxi(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                      bool xN_independent) {
    if (itau > 0) return 0.0;
    double bmc = bm_term(x) - cm_term();  // appears only in the form (b-c) in the equations
    double db_dxi = d_bm_term_dxi(x, i, xN_independent);
    double bracket = 1 - bmc * delta * rho_r;

    switch (idelta) {
        case 0:
            return delta * rho_r * db_dxi / bracket;
        case 1:
            return rho_r * db_dxi / pow(bracket, 2);
        case 2:
            return 2 * pow(rho_r, 2) * bmc * db_dxi / pow(bracket, 3);
        case 3:
            return 6 * pow(rho_r, 3) * pow(bmc, 2) * db_dxi / pow(bracket, 4);
        case 4:
            return 24 * pow(rho_r, 4) * pow(bmc, 3) * db_dxi / pow(bracket, 5);
        default:
            throw -1;
    }
}
double AbstractCubic::d2_psi_minus_dxidxj(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                          std::size_t j, bool xN_independent) {
    if (itau > 0) return 0.0;
    double bmc = bm_term(x) - cm_term();  // appears only in the form (b-c) in the equations
    double db_dxi = d_bm_term_dxi(x, i, xN_independent), db_dxj = d_bm_term_dxi(x, j, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(x, i, j, xN_independent);
    double bracket = 1 - bmc * delta * rho_r;

    switch (idelta) {
        case 0:
            return pow(delta * rho_r, 2) * db_dxi * db_dxj / pow(bracket, 2) + delta * rho_r * d2b_dxidxj / bracket;
        case 1:
            return 2 * delta * pow(rho_r, 2) * db_dxi * db_dxj / pow(bracket, 3) + rho_r * d2b_dxidxj / pow(bracket, 2);
        case 2:
            return 2 * pow(rho_r, 2) * db_dxi * db_dxj / pow(bracket, 4) * (2 * delta * rho_r * bmc + 1)
                   + 2 * pow(rho_r, 2) * bmc * d2b_dxidxj / pow(bracket, 3);
        case 3:
            return 12 * pow(rho_r, 3) * bmc * db_dxi * db_dxj / pow(bracket, 5) * (delta * rho_r * bmc + 1)
                   + 6 * pow(rho_r, 3) * pow(bmc, 2) * d2b_dxidxj / pow(bracket, 4);
        case 4:
            return 24 * pow(rho_r, 4) * pow(bmc, 2) * db_dxi * db_dxj / pow(bracket, 6) * (2 * delta * rho_r * bmc + 3)
                   + 24 * pow(rho_r, 4) * pow(bmc, 3) * d2b_dxidxj / pow(bracket, 5);
        default:
            throw -1;
    }
}
double AbstractCubic::d3_psi_minus_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                             std::size_t j, std::size_t k, bool xN_independent) {
    if (itau > 0) return 0.0;
    double bmc = bm_term(x) - cm_term();  // appears only in the form (b-c) in the equations
    double db_dxi = d_bm_term_dxi(x, i, xN_independent), db_dxj = d_bm_term_dxi(x, j, xN_independent), db_dxk = d_bm_term_dxi(x, k, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(x, i, j, xN_independent), d2b_dxidxk = d2_bm_term_dxidxj(x, i, k, xN_independent),
           d2b_dxjdxk = d2_bm_term_dxidxj(x, j, k, xN_independent), d3b_dxidxjdxk = d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent);
    double bracket = 1 - bmc * delta * rho_r;

    switch (idelta) {
        case 0:
            return delta * rho_r * d3b_dxidxjdxk / bracket + 2 * pow(delta * rho_r, 3) * db_dxi * db_dxj * db_dxk / pow(bracket, 3)
                   + pow(delta * rho_r, 2) / pow(bracket, 2) * (db_dxi * d2b_dxjdxk + db_dxj * d2b_dxidxk + db_dxk * d2b_dxidxj);
        case 1:
            return rho_r * d3b_dxidxjdxk / pow(bracket, 2) + 6 * pow(delta, 2) * pow(rho_r, 3) * db_dxi * db_dxj * db_dxk / pow(bracket, 4)
                   + 2 * delta * pow(rho_r, 2) / pow(bracket, 3) * (db_dxi * d2b_dxjdxk + db_dxj * d2b_dxidxk + db_dxk * d2b_dxidxj);
        default:
            throw -1;
    }
}
double AbstractCubic::PI_12(double delta, const std::vector<double>& x, std::size_t idelta) {
    double bm = bm_term(x);
    double cm = cm_term();
    switch (idelta) {
        case 0:
            return (1 + (Delta_1 * bm + cm) * rho_r * delta) * (1 + (Delta_2 * bm + cm) * rho_r * delta);
        case 1:
            return rho_r * (2 * (bm * Delta_1 + cm) * (bm * Delta_2 + cm) * delta * rho_r + (Delta_1 + Delta_2) * bm + 2 * cm);
        case 2:
            return 2 * (Delta_1 * bm + cm) * (Delta_2 * bm + cm) * pow(rho_r, 2);
        case 3:
            return 0;
        case 4:
            return 0;
        default:
            throw -1;
    }
}
double AbstractCubic::d_PI_12_dxi(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, bool xN_independent) {
    double bm = bm_term(x);
    double cm = cm_term();
    double db_dxi = d_bm_term_dxi(x, i, xN_independent);
    switch (idelta) {
        case 0:
            return delta * rho_r * db_dxi * (2 * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + cm * delta * rho_r));
        case 1:
            return rho_r * db_dxi * (4 * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + 2 * cm * delta * rho_r));
        case 2:
            return 2 * pow(rho_r, 2) * (2 * Delta_1 * Delta_2 * bm + (Delta_1 + Delta_2) * cm) * db_dxi;
        case 3:
            return 0;
        case 4:
            return 0;
        default:
            throw -1;
    }
}
double AbstractCubic::d2_PI_12_dxidxj(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j,
                                      bool xN_independent) {
    double bm = bm_term(x);
    double cm = cm_term();
    double db_dxi = d_bm_term_dxi(x, i, xN_independent), db_dxj = d_bm_term_dxi(x, j, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(x, i, j, xN_independent);
    switch (idelta) {
        case 0:
            return delta * rho_r
                   * (2 * Delta_1 * Delta_2 * delta * rho_r * db_dxi * db_dxj
                      + (2 * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + cm * delta * rho_r)) * d2b_dxidxj);
        case 1:
            return rho_r
                   * (4 * Delta_1 * Delta_2 * delta * rho_r * db_dxi * db_dxj
                      + (4 * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + 2 * cm * delta * rho_r)) * d2b_dxidxj);
        case 2:
            return 2 * pow(rho_r, 2)
                   * (2 * Delta_1 * Delta_2 * db_dxi * db_dxj + (2 * Delta_1 * Delta_2 * bm + (Delta_1 + Delta_2) * cm) * d2b_dxidxj);
        case 3:
            return 0;
        case 4:
            return 0;
        default:
            throw -1;
    }
}
double AbstractCubic::d3_PI_12_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k,
                                         bool xN_independent) {
    double bm = bm_term(x);
    double cm = cm_term();
    double db_dxi = d_bm_term_dxi(x, i, xN_independent), db_dxj = d_bm_term_dxi(x, j, xN_independent), db_dxk = d_bm_term_dxi(x, k, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(x, i, j, xN_independent), d2b_dxidxk = d2_bm_term_dxidxj(x, i, k, xN_independent),
           d2b_dxjdxk = d2_bm_term_dxidxj(x, j, k, xN_independent), d3b_dxidxjdxk = d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent);
    switch (idelta) {
        case 0:
            return delta * rho_r
                   * ((2 * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + cm * delta * rho_r)) * d3b_dxidxjdxk
                      + 2 * Delta_1 * Delta_2 * delta * rho_r * (db_dxi * d2b_dxjdxk + db_dxj * d2b_dxidxk + db_dxk * d2b_dxidxj));
        case 1:
            return rho_r
                   * ((4. * Delta_1 * Delta_2 * bm * delta * rho_r + (Delta_1 + Delta_2) * (1 + 2 * cm * delta * rho_r)) * d3b_dxidxjdxk
                      + 4 * Delta_1 * Delta_2 * delta * rho_r * (db_dxi * d2b_dxjdxk + db_dxj * d2b_dxidxk + db_dxk * d2b_dxidxj));
        default:
            throw -1;
    }
}
double AbstractCubic::psi_plus(double delta, const std::vector<double>& x, std::size_t idelta) {
    switch (idelta) {
        case 0:
            return A_term(delta, x) * c_term(x) / (Delta_1 - Delta_2);
        case 1:
            return rho_r / PI_12(delta, x, 0);
        case 2:
            return -rho_r / pow(PI_12(delta, x, 0), 2) * PI_12(delta, x, 1);
        case 3:
            return rho_r * (-PI_12(delta, x, 0) * PI_12(delta, x, 2) + 2 * pow(PI_12(delta, x, 1), 2)) / pow(PI_12(delta, x, 0), 3);
        case 4:
            // Term -PI_12(delta,x,0)*PI_12(delta,x,3) in the numerator is zero (and removed) since PI_12(delta,x,3) = 0
            return rho_r * (6 * PI_12(delta, x, 0) * PI_12(delta, x, 1) * PI_12(delta, x, 2) - 6 * pow(PI_12(delta, x, 1), 3))
                   / pow(PI_12(delta, x, 0), 4);
        default:
            throw -1;
    }
}
double AbstractCubic::d_psi_plus_dxi(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, bool xN_independent) {
    double bracket = 0;
    if (idelta == 0) {
        return (A_term(delta, x) * d_c_term_dxi(x, i, xN_independent) + c_term(x) * d_A_term_dxi(delta, x, i, xN_independent)) / (Delta_1 - Delta_2);
    }
    // All the terms with at least one delta derivative are multiplied by a common term of -rhor/PI12^2
    // So we just evaluate the bracketed term and then multiply by the common factor in the front
    switch (idelta) {
        case 1:
            bracket = d_PI_12_dxi(delta, x, 0, i, xN_independent);
            break;
        case 2:
            bracket = (d_PI_12_dxi(delta, x, 1, i, xN_independent)
                       + 2 / rho_r * PI_12(delta, x, 0) * PI_12(delta, x, 1) * d_psi_plus_dxi(delta, x, 1, i, xN_independent));
            break;
        case 3: {
            bracket =
              (d_PI_12_dxi(delta, x, 2, i, xN_independent)
               + 2 / rho_r * (pow(PI_12(delta, x, 1), 2) + PI_12(delta, x, 0) * PI_12(delta, x, 2)) * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
               + 4 / rho_r * PI_12(delta, x, 0) * PI_12(delta, x, 1) * d_psi_plus_dxi(delta, x, 2, i, xN_independent));
            break;
        }
        case 4:
            // d_PI_12_dxi(delta, x, 3, i, xN_independent) = 0, and PI_12(delta,x,0)*PI_12(delta,x,3) = 0, so removed from sum
            bracket =
              (6 / rho_r * PI_12(delta, x, 1) * PI_12(delta, x, 2) * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
               + 6 / rho_r * (pow(PI_12(delta, x, 1), 2) + PI_12(delta, x, 0) * PI_12(delta, x, 2)) * d_psi_plus_dxi(delta, x, 2, i, xN_independent)
               + 6 / rho_r * PI_12(delta, x, 0) * PI_12(delta, x, 1) * d_psi_plus_dxi(delta, x, 3, i, xN_independent));
            break;
        default:
            throw -1;
    }
    return -rho_r / pow(PI_12(delta, x, 0), 2) * bracket;
}
double AbstractCubic::d2_psi_plus_dxidxj(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j,
                                         bool xN_independent) {
    double bracket = 0;
    double PI12 = PI_12(delta, x, 0);
    if (idelta == 0) {
        return (A_term(delta, x) * d2_c_term_dxidxj(x, i, j, xN_independent) + c_term(x) * d2_A_term_dxidxj(delta, x, i, j, xN_independent)
                + d_A_term_dxi(delta, x, i, xN_independent) * d_c_term_dxi(x, j, xN_independent)
                + d_A_term_dxi(delta, x, j, xN_independent) * d_c_term_dxi(x, i, xN_independent))
               / (Delta_1 - Delta_2);
    }
    // All the terms with at least one delta derivative have a common factor of -1/PI_12^2 out front
    // so we just calculate the bracketed term and then multiply later on
    switch (idelta) {
        case 1:
            bracket = (rho_r * d2_PI_12_dxidxj(delta, x, 0, i, j, xN_independent)
                       + 2 * PI12 * d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_psi_plus_dxi(delta, x, 1, i, xN_independent));
            break;
        case 2:
            bracket = (rho_r * d2_PI_12_dxidxj(delta, x, 1, i, j, xN_independent)
                       + 2 * (PI12 * d_PI_12_dxi(delta, x, 1, j, xN_independent) + PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 0, j, xN_independent))
                           * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
                       + 2 * PI12 * PI_12(delta, x, 1) * d2_psi_plus_dxidxj(delta, x, 1, i, j, xN_independent)
                       + 2 * PI12 * d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_psi_plus_dxi(delta, x, 2, i, xN_independent));
            break;
        case 3: {
            bracket =
              (rho_r * d2_PI_12_dxidxj(delta, x, 2, i, j, xN_independent)
               + 2 * (PI12 * PI_12(delta, x, 2) + pow(PI_12(delta, x, 1), 2)) * d2_psi_plus_dxidxj(delta, x, 1, i, j, xN_independent)
               + 4 * (PI12 * d_PI_12_dxi(delta, x, 1, j, xN_independent) + PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 0, j, xN_independent))
                   * d_psi_plus_dxi(delta, x, 2, i, xN_independent)
               + 2
                   * (PI12 * d_PI_12_dxi(delta, x, 2, j, xN_independent) + 2 * PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 1, j, xN_independent)
                      + d_PI_12_dxi(delta, x, 0, j, xN_independent) * PI_12(delta, x, 2))
                   * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
               + 4 * PI12 * PI_12(delta, x, 1) * d2_psi_plus_dxidxj(delta, x, 2, i, j, xN_independent)
               + 2 * PI12 * d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_psi_plus_dxi(delta, x, 3, i, xN_independent));
            break;
        }
        case 4:
            // rho_r*d2_PI_12_dxidxj(delta, x, 3, i, j, xN_independent)  = 0
            // PI_12(delta, x, 3) = 0
            // PI12*d_PI_12_dxi(delta, x, 3, j, xN_independent) = 0
            // d_PI_12_dxi(delta, x, 0, j, xN_independent)*PI_12(delta, x, 3) = 0
            bracket =
              (+6 * (PI12 * PI_12(delta, x, 2) + pow(PI_12(delta, x, 1), 2)) * d2_psi_plus_dxidxj(delta, x, 2, i, j, xN_independent)
               + 6 * PI_12(delta, x, 1) * PI_12(delta, x, 2) * d2_psi_plus_dxidxj(delta, x, 1, i, j, xN_independent)
               + 6 * (PI12 * d_PI_12_dxi(delta, x, 1, j, xN_independent) + PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 0, j, xN_independent))
                   * d_psi_plus_dxi(delta, x, 3, i, xN_independent)
               + 6
                   * (PI12 * d_PI_12_dxi(delta, x, 2, j, xN_independent) + 2 * PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 1, j, xN_independent)
                      + d_PI_12_dxi(delta, x, 0, j, xN_independent) * PI_12(delta, x, 2))
                   * d_psi_plus_dxi(delta, x, 2, i, xN_independent)
               + 6
                   * (PI_12(delta, x, 1) * d_PI_12_dxi(delta, x, 2, j, xN_independent)
                      + PI_12(delta, x, 2) * d_PI_12_dxi(delta, x, 1, j, xN_independent))
                   * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
               + 6 * PI12 * PI_12(delta, x, 1) * d2_psi_plus_dxidxj(delta, x, 3, i, j, xN_independent)
               + 2 * PI12 * d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_psi_plus_dxi(delta, x, 4, i, xN_independent));
            break;
        default:
            throw -1;
    }
    return -1 / pow(PI12, 2) * bracket;
}
double AbstractCubic::d3_psi_plus_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j,
                                            std::size_t k, bool xN_independent) {
    double PI12 = PI_12(delta, x, 0);
    switch (idelta) {
        case 0:
            return (A_term(delta, x) * d3_c_term_dxidxjdxk(x, i, j, k, xN_independent)
                    + c_term(x) * d3_A_term_dxidxjdxk(delta, x, i, j, k, xN_independent)
                    + d_A_term_dxi(delta, x, i, xN_independent) * d2_c_term_dxidxj(x, j, k, xN_independent)
                    + d_A_term_dxi(delta, x, j, xN_independent) * d2_c_term_dxidxj(x, i, k, xN_independent)
                    + d_A_term_dxi(delta, x, k, xN_independent) * d2_c_term_dxidxj(x, i, j, xN_independent)
                    + d_c_term_dxi(x, i, xN_independent) * d2_A_term_dxidxj(delta, x, j, k, xN_independent)
                    + d_c_term_dxi(x, j, xN_independent) * d2_A_term_dxidxj(delta, x, i, k, xN_independent)
                    + d_c_term_dxi(x, k, xN_independent) * d2_A_term_dxidxj(delta, x, i, j, xN_independent))
                   / (Delta_1 - Delta_2);
        case 1:
            return -1 / pow(PI12, 2)
                   * (rho_r * d3_PI_12_dxidxjdxk(delta, x, 0, i, j, k, xN_independent)
                      + 2
                          * (PI12 * d2_PI_12_dxidxj(delta, x, 0, j, k, xN_independent)
                             + d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_PI_12_dxi(delta, x, 0, k, xN_independent))
                          * d_psi_plus_dxi(delta, x, 1, i, xN_independent)
                      + 2 * PI12 * d_PI_12_dxi(delta, x, 0, j, xN_independent) * d2_psi_plus_dxidxj(delta, x, 1, i, k, xN_independent)
                      + 2 * PI12 * d_PI_12_dxi(delta, x, 0, k, xN_independent) * d2_psi_plus_dxidxj(delta, x, 1, i, j, xN_independent));
        default:
            throw -1;
    }
}

double AbstractCubic::tau_times_a(double tau, const std::vector<double>& x, std::size_t itau) {
    if (itau == 0) {
        return tau * am_term(tau, x, 0);
    } else {
        return tau * am_term(tau, x, itau) + itau * am_term(tau, x, itau - 1);
    }
}
double AbstractCubic::d_tau_times_a_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent) {
    if (itau == 0) {
        return tau * d_am_term_dxi(tau, x, 0, i, xN_independent);
    } else {
        return tau * d_am_term_dxi(tau, x, itau, i, xN_independent) + itau * d_am_term_dxi(tau, x, itau - 1, i, xN_independent);
    }
}
double AbstractCubic::d2_tau_times_a_dxidxj(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j,
                                            bool xN_independent) {
    if (itau == 0) {
        return tau * d2_am_term_dxidxj(tau, x, 0, i, j, xN_independent);
    } else {
        return tau * d2_am_term_dxidxj(tau, x, itau, i, j, xN_independent) + itau * d2_am_term_dxidxj(tau, x, itau - 1, i, j, xN_independent);
    }
}
double AbstractCubic::d3_tau_times_a_dxidxjdxk(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j,
                                               std::size_t k, bool xN_independent) {
    if (itau == 0) {
        return tau * d3_am_term_dxidxjdxk(tau, x, 0, i, j, k, xN_independent);
    } else {
        return tau * d3_am_term_dxidxjdxk(tau, x, itau, i, j, k, xN_independent)
               + itau * d3_am_term_dxidxjdxk(tau, x, itau - 1, i, j, k, xN_independent);
    }
}
double AbstractCubic::alphar(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta) {
    return psi_minus(delta, x, itau, idelta) - tau_times_a(tau, x, itau) / (R_u * T_r) * psi_plus(delta, x, idelta);
}
double AbstractCubic::d_alphar_dxi(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                   bool xN_independent) {
    return (d_psi_minus_dxi(delta, x, itau, idelta, i, xN_independent)
            - 1 / (R_u * T_r)
                * (d_tau_times_a_dxi(tau, x, itau, i, xN_independent) * psi_plus(delta, x, idelta)
                   + tau_times_a(tau, x, itau) * d_psi_plus_dxi(delta, x, idelta, i, xN_independent)));
}
double AbstractCubic::d2_alphar_dxidxj(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                       std::size_t j, bool xN_independent) {
    return (d2_psi_minus_dxidxj(delta, x, itau, idelta, i, j, xN_independent)
            - 1 / (R_u * T_r)
                * (d2_tau_times_a_dxidxj(tau, x, itau, i, j, xN_independent) * psi_plus(delta, x, idelta)
                   + d_tau_times_a_dxi(tau, x, itau, i, xN_independent) * d_psi_plus_dxi(delta, x, idelta, j, xN_independent)
                   + d_tau_times_a_dxi(tau, x, itau, j, xN_independent) * d_psi_plus_dxi(delta, x, idelta, i, xN_independent)
                   + tau_times_a(tau, x, itau) * d2_psi_plus_dxidxj(delta, x, idelta, i, j, xN_independent)));
}
double AbstractCubic::d3_alphar_dxidxjdxk(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                          std::size_t j, std::size_t k, bool xN_independent) {
    return (d3_psi_minus_dxidxjdxk(delta, x, itau, idelta, i, j, k, xN_independent)
            - 1 / (R_u * T_r)
                * (d2_tau_times_a_dxidxj(tau, x, itau, i, j, xN_independent) * d_psi_plus_dxi(delta, x, idelta, k, xN_independent)
                   + d3_tau_times_a_dxidxjdxk(tau, x, itau, i, j, k, xN_independent) * psi_plus(delta, x, idelta)

                   + d_tau_times_a_dxi(tau, x, itau, i, xN_independent) * d2_psi_plus_dxidxj(delta, x, idelta, j, k, xN_independent)
                   + d2_tau_times_a_dxidxj(tau, x, itau, i, k, xN_independent) * d_psi_plus_dxi(delta, x, idelta, j, xN_independent)

                   + d_tau_times_a_dxi(tau, x, itau, j, xN_independent) * d2_psi_plus_dxidxj(delta, x, idelta, i, k, xN_independent)
                   + d2_tau_times_a_dxidxj(tau, x, itau, j, k, xN_independent) * d_psi_plus_dxi(delta, x, idelta, i, xN_independent)

                   + tau_times_a(tau, x, itau) * d3_psi_plus_dxidxjdxk(delta, x, idelta, i, j, k, xN_independent)
                   + d_tau_times_a_dxi(tau, x, itau, k, xN_independent) * d2_psi_plus_dxidxj(delta, x, idelta, i, j, xN_independent)));
}

double SRK::a0_ii(std::size_t i) {
    // Exact value: 1/(9*(2^(1/3)-1)); see Bell and Deiters, IECR, 2021
    double a = 0.42748023335403414043900347952220 * R_u * R_u * Tc[i] * Tc[i] / pc[i];
    return a;
}
double SRK::b0_ii(std::size_t i) {
    // Exact value: (2^(1/3)-1)/3; see Bell and Deiters, IECR, 2021
    double b = 0.08664034999649577215890158147700 * R_u * Tc[i] / pc[i];
    return b;
}
double SRK::m_ii(std::size_t i) {
    // Values from Soave, 1972 (Equilibrium constants from a ..)
    double omega = acentric[i];
    double m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    return m;
}

double PengRobinson::a0_ii(std::size_t i) {
    // Exact value; see Bell and Deiters, IECR, 2021
    double a = 0.45723552892138218938000849856422 * R_u * R_u * Tc[i] * Tc[i] / pc[i];
    return a;
}
double PengRobinson::b0_ii(std::size_t i) {
    // Exact value; see Bell and Deiters, IECR, 2021
    double b = 0.07779607390388455972148597969400 * R_u * Tc[i] / pc[i];
    return b;
}
double PengRobinson::m_ii(std::size_t i) {
    double omega = acentric[i];
    double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    return m;
}
