// Structure-aware PC-SAFT: the whole A_ij triangle with the expensive arithmetic done in
// univariate Taylor arithmetic.  With eta = zeta_3 = rho q(T) and zeta_n = eta r_n(T),
// at fixed composition
//
//   alphar = mbar [A1(T) eta/(1-eta) + A2(T) eta/(1-eta)^2 + A3(T) ln(1-eta)]
//          - sum_i x_i (m_i - 1) ln[ 1/(1-eta) + a_i(T) eta/(1-eta)^2 + b_i(T) eta^2/(1-eta)^3 ]
//          + B1(T) eta I1(eta) + B2(T) eta C1(eta) I2(eta)
//
//   A1 = 3 r1 r2 / r0,  A2 = r2^3 / r0,  A3 = r2^3 / r0 - 1,  a_i = 3 D_i r2,  b_i = 2 D_i^2 r2^2
//   B1 = -12 E1 (1/T) / s3,  B2 = -6 mbar E2 (1/T)^2 / s3,  D_i = d_i / 2
//
// so every eta-function phi_k is a univariate Taylor series in eta, every coefficient is a
// univariate Taylor series in u (1/T = (1/T0)(1+u)), and the only bivariate work is:
// the powers of H = eta(u,v) - eta0 (N-1 products, shared), one axpy per phi_k, one
// univariate-by-bivariate product per coefficient, and one bivariate log per component.
#include <cmath>
#include "api_all.hpp"
#include "taylor1.hpp"
using namespace spike;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    template <int N>
    static void all(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        using U = nd::Taylor1<N>;  // univariate in u (temperature) or in eta
        using B = nd::Taylor2<N>;  // bivariate in (u, v)
        constexpr double pi = 3.14159265358979323846, NA = 6.02214076e23;
        const std::size_t Nc = p.N;

        // ---- temperature side, univariate in u
        const double ti0 = 1.0 / T;
        const U ti = U::variable(ti0, ti0);  // 1/T = ti0 (1 + u)
        std::array<U, NMAX> Dh{};            // d_i / 2
        U s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        double mbar = 0, E1 = 0, E2 = 0;
        for (std::size_t i = 0; i < Nc; ++i) {
            const U di = p.sigma_A[i] * (1.0 - 0.12 * exp(ti * (-3.0 * p.eps_k[i])));
            Dh[i] = di * 0.5;
            const double xm = x[i] * p.m[i];
            const U d2 = di * di;
            s0 += U(xm);
            s1 += di * xm;
            s2 += d2 * xm;
            s3 += d2 * di * xm;
            mbar += xm;
            for (std::size_t j = 0; j < Nc; ++j) {
                const double sij = 0.5 * (p.sigma_A[i] + p.sigma_A[j]);
                const double pre = x[i] * x[j] * p.m[i] * p.m[j] * sij * sij * sij;
                const double eij = std::sqrt(p.eps_k[i] * p.eps_k[j]);
                E1 += pre * eij;
                E2 += pre * eij * eij;
            }
        }
        const U r0 = s0 / s3, r1 = s1 / s3, r2 = s2 / s3;
        const U r2c = r2 * r2 * r2;
        const U A1 = 3.0 * r1 * r2 / r0, A2 = r2c / r0, A3 = A2 + (-1.0);
        const U B1 = (-12.0 * E1) * ti / s3, B2 = (-6.0 * mbar * E2) * ti * ti / s3;
        const U q = (pi / 6.0) * NA * 1e-30 * s3;  // eta = rho q

        // ---- density side, univariate in eta about eta0
        const double eta0 = rho * q.t[0];
        const U eta = U::variable(eta0);
        const U om = 1.0 - eta, iom = 1.0 / om;
        const U phi1 = eta * iom, phi2 = phi1 * iom, phi3 = log(om), phi4 = iom, psi = phi2 * eta * iom;
        const double c1m = (mbar - 1.0) / mbar, c2m = c1m * (mbar - 2.0) / mbar;
        U I1 = 0.0, I2 = 0.0, ek = 1.0;
        for (int k = 0; k < 7; ++k) {
            I1 += ek * (gs2001::a[0][k] + c1m * gs2001::a[1][k] + c2m * gs2001::a[2][k]);
            I2 += ek * (gs2001::b[0][k] + c1m * gs2001::b[1][k] + c2m * gs2001::b[2][k]);
            ek = ek * eta;
        }
        const U om2 = om * om, tme = 2.0 - eta, e2 = eta * eta;
        const U C1 = 1.0
                     / (1.0 + mbar * (8.0 * eta - 2.0 * e2) / (om2 * om2)
                        + (1.0 - mbar) * (20.0 * eta - 27.0 * e2 + 12.0 * e2 * eta - 2.0 * e2 * e2) / (om * tme * om * tme));
        const U phia = eta * I1, phib = eta * C1 * I2;

        // ---- the only bivariate part: eta(u, v) = rho0 (1 + v) q(u)
        B etaB = nd::embed_u(q) * rho;
        for (int i = 0; i + 1 <= N; ++i)
            etaB.c[B::idx(i, 1)] = rho * q.t[i];  // the v * q(u) part
        B H = etaB;
        H.c[0] = 0.0;
        const nd::PowerBasis<N> P(H);
        const B P1 = P.compose(phi1), P2 = P.compose(phi2), P3 = P.compose(phi3), P4 = P.compose(phi4), Ps = P.compose(psi);
        B a = (nd::mul_u(A1, P1) + nd::mul_u(A2, P2) + nd::mul_u(A3, P3)) * mbar;
        a += nd::mul_u(B1, P.compose(phia)) + nd::mul_u(B2, P.compose(phib));
        for (std::size_t i = 0; i < Nc; ++i) {
            const B g = P4 + nd::mul_u(3.0 * Dh[i] * r2, P2) + nd::mul_u(2.0 * Dh[i] * Dh[i] * r2 * r2, Ps);
            a -= log(g) * (x[i] * (p.m[i] - 1.0));
        }

        fill_nan(o);
        double fi = 1;
        for (int i = 0; i <= N; ++i, fi *= i) {
            double fj = 1;
            for (int j = 0; i + j <= N; ++j, fj *= j)
                o[idx2(i, j)] = fi * fj * a(i, j);
        }
    }
};
}  // namespace
SPIKE_EXPORT_ALL(structured)
