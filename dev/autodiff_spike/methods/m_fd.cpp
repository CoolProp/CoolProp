// Central finite differences on the double model, textbook step sizes
// h ~ eps^(1/(n+2)) * |x| for an n-th derivative with an O(h^2) stencil.
#include "api.hpp"
#include <cmath>
using namespace spike;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static constexpr double eps = 2.220446049250313e-16;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        auto f = [&](double r) { return M::alphar(p, T, r, x); };
        const double f0 = f(rho);
        const double h1 = std::cbrt(eps) * rho, h2 = std::pow(eps, 0.25) * rho, h3 = std::pow(eps, 0.2) * rho;
        o[0] = f0;
        o[1] = rho * (f(rho + h1) - f(rho - h1)) / (2 * h1);
        o[2] = rho * rho * (f(rho + h2) - 2 * f0 + f(rho - h2)) / (h2 * h2);
        o[3] = rho * rho * rho * (f(rho + 2 * h3) - 2 * f(rho + h3) + 2 * f(rho - h3) - f(rho - 2 * h3)) / (2 * h3 * h3 * h3);
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const double ti = 1.0 / T;
        auto f = [&](double t) { return M::alphar(p, 1.0 / t, rho, x); };
        const double f0 = f(ti), h1 = std::cbrt(eps) * ti, h2 = std::pow(eps, 0.25) * ti;
        o[0] = f0;
        o[1] = ti * (f(ti + h1) - f(ti - h1)) / (2 * h1);
        o[2] = ti * ti * (f(ti + h2) - 2 * f0 + f(ti - h2)) / (h2 * h2);
    }
    static void Ar11(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const double ti = 1.0 / T;
        auto f = [&](double t, double r) { return M::alphar(p, 1.0 / t, r, x); };
        const double ht = std::pow(eps, 0.25) * ti, hr = std::pow(eps, 0.25) * rho;
        const double d = (f(ti + ht, rho + hr) - f(ti + ht, rho - hr) - f(ti - ht, rho + hr) + f(ti - ht, rho - hr)) / (4 * ht * hr);
        o[0] = ti * rho * d;
    }
    static void gradx(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        for (std::size_t i = 0; i < p.N; ++i) {
            const double h = std::cbrt(eps) * std::max(std::abs(x[i]), 1e-3);
            XArr xp = x, xm = x;
            xp[i] += h;
            xm[i] -= h;
            o[i] = (M::alphar(p, T, rho, xp) - M::alphar(p, T, rho, xm)) / (2 * h);
        }
    }
};
}  // namespace
SPIKE_EXPORT(fd)
