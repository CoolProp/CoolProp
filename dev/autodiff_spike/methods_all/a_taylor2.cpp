// One pass of a bivariate truncated Taylor type in scaled variables
// 1/T = (1/T0)(1+u), rho = rho0(1+v)  =>  A_ij = i! j! c_ij directly.
#include "api_all.hpp"
#include "taylor2.hpp"
using namespace spike;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    template <int N>
    static void all(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        using TT = nd::Taylor2<N>;
        TT ti(1.0 / T), r(rho);
        ti.c[TT::idx(1, 0)] = 1.0 / T;
        r.c[TT::idx(0, 1)] = rho;
        const TT a = M::alphar(p, TT(1.0 / ti), r, x);
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
SPIKE_EXPORT_ALL(taylor2)
