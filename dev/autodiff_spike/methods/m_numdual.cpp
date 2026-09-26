// num-dual style: one purpose-built eager type per derivative request.
#include "api.hpp"
#include "numdual.hpp"
using namespace spike;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const auto a = M::alphar(p, T, nd::Dual3<double>(rho, 1.0, 0.0, 0.0), x);
        o[0] = a.re;
        o[1] = rho * a.v1;
        o[2] = rho * rho * a.v2;
        o[3] = rho * rho * rho * a.v3;
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const nd::Dual2<double> ti(1.0 / T, 1.0, 0.0);
        const auto a = M::alphar(p, 1.0 / ti, rho, x);
        o[0] = a.re;
        o[1] = ti.re * a.v1;
        o[2] = ti.re * ti.re * a.v2;
    }
    static void Ar11(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        using H = nd::HyperDual<double>;
        const H ti(1.0 / T, 1.0, 0.0, 0.0), r(rho, 0.0, 1.0, 0.0);
        o[0] = ti.re * rho * M::alphar(p, 1.0 / ti, r, x).e12;
    }
    static void gradx(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        using V = nd::DualVec<double, NMAX>;
        std::array<V, NMAX> xd{};
        for (std::size_t i = 0; i < NMAX; ++i) {
            xd[i].re = x[i];
            xd[i].eps[i] = 1.0;
        }
        const auto a = M::alphar(p, T, rho, xd);
        for (std::size_t i = 0; i < p.N; ++i)
            o[i] = a.eps[i];
    }
};
}  // namespace
SPIKE_EXPORT(numdual)
