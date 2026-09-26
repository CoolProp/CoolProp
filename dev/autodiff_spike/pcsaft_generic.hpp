// Minimal, dependency-free PC-SAFT (hard chain + dispersion, Gross & Sadowski 2001,
// doi:10.1021/ie0003887), written once and generic over the number type so that every
// differentiation method in this spike runs *identical* model code.
//
// The Tag parameter only exists to mint distinct model types for the "K models"
// compile-time scaling experiment; it does not change the math.
#pragma once
#include <cmath>
#include <type_traits>
#include <array>
#include <cstddef>

namespace spike {

// Pick the "active" type of two operands.  In every call in this spike at most one AD
// type is live at a time, so "whichever is not double" is sufficient and keeps the
// T-only (d_i) and x-only (mbar) intermediates as plain doubles when they can be.
template <class A, class B>
using promote_t = std::conditional_t<std::is_same_v<A, double>, B, A>;

inline constexpr std::size_t NMAX = 3;  // fixed capacity -> no heap traffic in the timed path

struct PCSAFTParams
{
    std::size_t N = 1;
    std::array<double, NMAX> m{}, sigma_A{}, eps_k{};  // segments, sigma [Angstrom], epsilon/k [K]
};

namespace gs2001 {
inline constexpr std::array<std::array<double, 7>, 3> a = {
  {{{0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084}},
   {{-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930}},
   {{-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368}}}};
inline constexpr std::array<std::array<double, 7>, 3> b = {
  {{{0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612}},
   {{-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346}},
   {{0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585}}}};
}  // namespace gs2001

template <int Tag = 0>
struct PCSAFTModel
{
    // alphar(T [K], rho [mol/m^3], x[]) with independent number types for each input.
    template <class TT, class TR, class TX>
    static auto alphar(const PCSAFTParams& p, const TT& T, const TR& rho, const std::array<TX, NMAX>& x) {
        using std::exp;
        using std::log;
        using TTX = promote_t<TT, TX>;
        using R = promote_t<TTX, TR>;
        constexpr double pi = 3.14159265358979323846;
        constexpr double NA = 6.02214076e23;
        const std::size_t N = p.N;

        // Temperature-dependent segment diameter
        std::array<TT, NMAX> d{};
        for (std::size_t i = 0; i < N; ++i) {
            d[i] = p.sigma_A[i] * (1.0 - 0.12 * exp(-3.0 * p.eps_k[i] / T));
        }
        const R rhoN = rho * (NA * 1e-30);  // number density in 1/A^3

        // zeta_n = pi/6 rhoN sum x_i m_i d_i^n
        TTX s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        TX mbar = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            TTX xm = x[i] * p.m[i];
            s0 += xm;
            s1 += xm * d[i];
            s2 += xm * d[i] * d[i];
            s3 += xm * d[i] * d[i] * d[i];
            mbar += x[i] * p.m[i];
        }
        const R z0 = (pi / 6.0) * rhoN * s0, z1 = (pi / 6.0) * rhoN * s1, z2 = (pi / 6.0) * rhoN * s2, z3 = (pi / 6.0) * rhoN * s3;
        const R om = 1.0 - z3;  // 1 - zeta3

        // Hard sphere + hard chain
        const R ahs = (3.0 * z1 * z2 / om + z2 * z2 * z2 / (z3 * om * om) + (z2 * z2 * z2 / (z3 * z3) - z0) * log(om)) / z0;
        R ahc = mbar * ahs;
        for (std::size_t i = 0; i < N; ++i) {
            const TT dij = d[i] / 2.0;  // d_i d_i/(d_i+d_i)
            const R gii = 1.0 / om + dij * 3.0 * z2 / (om * om) + dij * dij * 2.0 * z2 * z2 / (om * om * om);
            ahc -= x[i] * (p.m[i] - 1.0) * log(gii);
        }

        // Dispersion
        const R& eta = z3;
        const TX c1m = (mbar - 1.0) / mbar, c2m = c1m * (mbar - 2.0) / mbar;
        R I1 = 0.0, I2 = 0.0, etai = 1.0;
        for (int i = 0; i < 7; ++i) {
            const TX ai = gs2001::a[0][i] + c1m * gs2001::a[1][i] + c2m * gs2001::a[2][i];
            const TX bi = gs2001::b[0][i] + c1m * gs2001::b[1][i] + c2m * gs2001::b[2][i];
            I1 += ai * etai;
            I2 += bi * etai;
            etai *= eta;
        }
        TTX m2es3 = 0.0, m2e2s3 = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                const double sij = 0.5 * (p.sigma_A[i] + p.sigma_A[j]);
                const TT eij_T = std::sqrt(p.eps_k[i] * p.eps_k[j]) / T;  // k_ij = 0
                const TTX pre = x[i] * x[j] * (p.m[i] * p.m[j] * sij * sij * sij);
                m2es3 += pre * eij_T;
                m2e2s3 += pre * eij_T * eij_T;
            }
        }
        const R om2 = om * om;
        const R C1 = 1.0
                     / (1.0 + mbar * (8.0 * eta - 2.0 * eta * eta) / (om2 * om2)
                        + (1.0 - mbar) * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta * eta * eta - 2.0 * eta * eta * eta * eta)
                            / (om * (2.0 - eta) * om * (2.0 - eta)));
        const R adisp = -2.0 * pi * rhoN * I1 * m2es3 - pi * rhoN * mbar * C1 * I2 * m2e2s3;
        return R(ahc + adisp);
    }
};

}  // namespace spike
