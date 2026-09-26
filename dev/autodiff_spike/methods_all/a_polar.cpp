// Polarization with autodiff's univariate Taylor type: N+1 passes of Real<N> along
// directions (u, v) = (1, s_m) t, then per degree k a Vandermonde solve in s recovers
// c_{k-j, j}.  Needs nothing beyond the Real<N> type teqp already uses.
#include "api_all.hpp"
#include <autodiff/forward/real.hpp>
#include <cmath>
using namespace spike;
using namespace autodiff;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int N>
struct Polar
{
    std::array<double, N + 1> s;
    std::array<std::array<double, N + 1>, N + 1> Vinv;  // inverse of V[m][j] = s_m^j
    Polar() {
        const double pi = 3.14159265358979323846;
        for (int m = 0; m <= N; ++m)
            s[m] = std::cos((2 * m + 1) * pi / (2 * (N + 1)));  // Chebyshev nodes
        // Gauss-Jordan on [V | I]
        std::array<std::array<double, 2 * (N + 1)>, N + 1> A{};
        for (int m = 0; m <= N; ++m) {
            for (int j = 0; j <= N; ++j)
                A[m][j] = std::pow(s[m], j);
            A[m][N + 1 + m] = 1.0;
        }
        for (int c = 0; c <= N; ++c) {
            int piv = c;
            for (int r = c + 1; r <= N; ++r)
                if (std::abs(A[r][c]) > std::abs(A[piv][c])) piv = r;
            std::swap(A[c], A[piv]);
            const double d = A[c][c];
            for (auto& v : A[c])
                v /= d;
            for (int r = 0; r <= N; ++r)
                if (r != c) {
                    const double f = A[r][c];
                    for (int k = 0; k < 2 * (N + 1); ++k)
                        A[r][k] -= f * A[c][k];
                }
        }
        for (int j = 0; j <= N; ++j)
            for (int m = 0; m <= N; ++m)
                Vinv[j][m] = A[j][N + 1 + m];
    }
};
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    template <int N>
    static void all(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        static const Polar<N> P;
        using R = Real<N, double>;
        const double ti0 = 1.0 / T;
        std::array<std::array<double, N + 1>, N + 1> g;  // g[m][k] = k-th normalized Taylor coeff along direction m
        for (int m = 0; m <= N; ++m) {
            const double s = P.s[m];
            R t = 0.0;
            auto f = [&](const R& tt) {
                const R ti = ti0 * (1.0 + tt), r = rho * (1.0 + s * tt);
                return M::alphar(p, R(1.0 / ti), r, x);
            };
            const auto d = derivatives(f, along(1), at(t));
            double kf = 1;
            for (int k = 0; k <= N; ++k, kf *= k)
                g[m][k] = d[k] / kf;
        }
        fill_nan(o);
        for (int k = 0; k <= N; ++k) {
            double fk = 1;  // factorials for i!, j! below
            for (int j = 0; j <= k; ++j) {
                double cj = 0;
                for (int m = 0; m <= N; ++m)
                    cj += P.Vinv[j][m] * g[m][k];
                const int i = k - j;
                double fi = 1, fj = 1;
                for (int q = 2; q <= i; ++q)
                    fi *= q;
                for (int q = 2; q <= j; ++q)
                    fj *= q;
                o[idx2(i, j)] = fi * fj * cj;
            }
            (void)fk;
        }
    }
};
}  // namespace
SPIKE_EXPORT_ALL(polar)
