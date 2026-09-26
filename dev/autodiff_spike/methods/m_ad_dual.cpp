// autodiff, nested (expression-template) duals only: dual3rd / dual2nd everywhere.
// This is the configuration that stresses template depth the most.
#include "api.hpp"
#include <autodiff/forward/dual.hpp>
using namespace spike;
using namespace autodiff;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        dual3rd r = rho;
        auto f = [&](const dual3rd& rr) { return M::alphar(p, T, rr, x); };
        auto d = derivatives(f, wrt(r, r, r), at(r));
        for (int n = 0; n <= 3; ++n)
            o[n] = std::pow(rho, n) * d[n];
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        dual2nd ti = 1.0 / T;
        auto f = [&](const dual2nd& t) { return M::alphar(p, dual2nd(1.0 / t), rho, x); };
        auto d = derivatives(f, wrt(ti, ti), at(ti));
        for (int n = 0; n <= 2; ++n)
            o[n] = std::pow(1.0 / T, n) * d[n];
    }
    static void Ar11(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        dual2nd ti = 1.0 / T, r = rho;
        auto f = [&](const dual2nd& t, const dual2nd& rr) { return M::alphar(p, dual2nd(1.0 / t), rr, x); };
        auto d = derivatives(f, wrt(ti, r), at(ti, r));
        o[0] = (1.0 / T) * rho * d[2];
    }
    static void gradx(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        std::array<dual, NMAX> xd{};
        for (std::size_t i = 0; i < p.N; ++i)
            xd[i] = x[i];
        for (std::size_t i = 0; i < p.N; ++i) {
            xd[i].grad = 1.0;
            o[i] = derivative<1>(M::alphar(p, T, rho, xd));
            xd[i].grad = 0.0;
        }
    }
};
}  // namespace
SPIKE_EXPORT(ad_dual)
