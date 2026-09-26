// Multicomplex step (Bell, Deiters & Leal), the teqp alternative backend.
#include "api.hpp"
#include "MultiComplex.hpp"
#include <functional>
using namespace spike;
using MC = mcx::MultiComplex<double>;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        std::function<MC(const MC&)> f = [&](const MC& r) { return M::alphar(p, T, r, x); };
        auto d = mcx::diff_mcx1(f, rho, 3, true);
        for (int n = 0; n <= 3; ++n)
            o[n] = std::pow(rho, n) * d[n];
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        std::function<MC(const MC&)> f = [&](const MC& t) { return M::alphar(p, MC(1.0 / t), rho, x); };
        auto d = mcx::diff_mcx1(f, 1.0 / T, 2, true);
        for (int n = 0; n <= 2; ++n)
            o[n] = std::pow(1.0 / T, n) * d[n];
    }
    static void Ar11(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        using fcn_t = std::function<MC(const std::valarray<MC>&)>;
        const fcn_t f = [&](const std::valarray<MC>& z) { return M::alphar(p, MC(1.0 / z[0]), z[1], x); };
        const std::vector<double> xs = {1.0 / T, rho};
        const std::vector<int> order = {1, 1};
        o[0] = (1.0 / T) * rho * mcx::diff_mcxN(f, xs, order);
    }
    static void gradx(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        for (std::size_t i = 0; i < p.N; ++i) {
            std::function<MC(const MC&)> f = [&](const MC& xi) {
                std::array<MC, NMAX> xm;
                for (std::size_t j = 0; j < NMAX; ++j)
                    xm[j] = (j == i) ? xi : MC(x[j]);
                return M::alphar(p, T, rho, xm);
            };
            o[i] = mcx::diff_mcx1(f, x[i], 1)[0];
        }
    }
};
}  // namespace
SPIKE_EXPORT(mcx)
