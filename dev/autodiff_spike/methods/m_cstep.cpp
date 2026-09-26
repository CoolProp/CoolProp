// Complex step: exact-to-roundoff first derivatives only, std::complex, no library.
#include "api.hpp"
#include <complex>
using namespace spike;
using C = std::complex<double>;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static constexpr double h = 1e-100;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const C a = M::alphar(p, T, C(rho, h), x);
        o[0] = a.real();
        o[1] = rho * a.imag() / h;
        o[2] = o[3] = NaN;
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        const C a = M::alphar(p, 1.0 / C(1.0 / T, h), rho, x);
        o[0] = a.real();
        o[1] = (1.0 / T) * a.imag() / h;
        o[2] = NaN;
    }
    static void Ar11(const PCSAFTParams&, double, double, const XArr&, double* o) {
        o[0] = NaN;
    }
    static void gradx(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        for (std::size_t i = 0; i < p.N; ++i) {
            std::array<C, NMAX> xc;
            for (std::size_t j = 0; j < NMAX; ++j)
                xc[j] = x[j];
            xc[i] = C(x[i], h);
            o[i] = M::alphar(p, T, rho, xc).imag() / h;
        }
    }
};
}  // namespace
SPIKE_EXPORT(cstep)
