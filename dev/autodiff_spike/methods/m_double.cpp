// Baseline: plain double, value only.  Sets the floor for compile time and run time.
#include "api.hpp"
using namespace spike;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    static void Ar0n(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        o[0] = M::alphar(p, T, rho, x);
        o[1] = o[2] = o[3] = NaN;
    }
    static void Arn0(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        o[0] = M::alphar(p, T, rho, x);
        o[1] = o[2] = NaN;
    }
    static void Ar11(const PCSAFTParams&, double, double, const XArr&, double* o) {
        o[0] = NaN;
    }
    static void gradx(const PCSAFTParams& p, double, double, const XArr&, double* o) {
        for (std::size_t i = 0; i < p.N; ++i)
            o[i] = NaN;
    }
};
}  // namespace
SPIKE_EXPORT(double)
