// teqp's scheme applied to the whole triangle: Real<N> in rho, Real<N> in 1/T, and one
// HigherOrderDual<i+j> pass per mixed derivative (i, j >= 1), as get_Agenxy does.
#include "api_all.hpp"
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/real.hpp>
#include <cmath>
#include <tuple>
using namespace spike;
using namespace autodiff;
namespace {  // internal linkage: every TU has its own Impl<Tag>
template <typename T, size_t... I>
auto dup_impl(const T& v, std::index_sequence<I...>) {
    return std::make_tuple((static_cast<void>(I), v)...);
}
template <int n, typename T>
auto dup(const T& v) {
    return dup_impl(v, std::make_index_sequence<n>());
}
struct wrt_helper
{
    template <typename... Args>
    auto operator()(Args&&... args) const {
        return Wrt<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
    }
};
template <int Tag>
struct Impl
{
    using M = PCSAFTModel<Tag>;
    template <int iT, int iD>
    static double mixed(const PCSAFTParams& p, double T, double rho, const XArr& x) {
        using ad = HigherOrderDual<iT + iD, double>;
        ad ti = 1.0 / T, r = rho;
        auto f = [&](const ad& t, const ad& rr) { return M::alphar(p, ad(1.0 / t), rr, x); };
        auto w = std::tuple_cat(dup<iT>(std::ref(ti)), dup<iD>(std::ref(r)));
        auto d = derivatives(f, std::apply(wrt_helper(), w), at(ti, r));
        return std::pow(1.0 / T, iT) * std::pow(rho, iD) * d[d.size() - 1];
    }
    template <int N, int i, int j>
    static void mixed_loop(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        if constexpr (i + j <= N) o[idx2(i, j)] = mixed<i, j>(p, T, rho, x);
        if constexpr (j + 1 < N)
            mixed_loop<N, i, j + 1>(p, T, rho, x, o);
        else if constexpr (i + 1 < N)
            mixed_loop<N, i + 1, 1>(p, T, rho, x, o);
    }
    template <int N>
    static void all(const PCSAFTParams& p, double T, double rho, const XArr& x, double* o) {
        fill_nan(o);
        {
            Real<N, double> r = rho;
            auto f = [&](const Real<N, double>& rr) { return M::alphar(p, T, rr, x); };
            const auto d = derivatives(f, along(1), at(r));
            for (int j = 0; j <= N; ++j)
                o[idx2(0, j)] = std::pow(rho, j) * d[j];
        }
        {
            Real<N, double> ti = 1.0 / T;
            auto f = [&](const Real<N, double>& t) { return M::alphar(p, Real<N, double>(1.0 / t), rho, x); };
            const auto d = derivatives(f, along(1), at(ti));
            for (int i = 1; i <= N; ++i)
                o[idx2(i, 0)] = std::pow(1.0 / T, i) * d[i];
        }
        if constexpr (N >= 2) mixed_loop<N, 1, 1>(p, T, rho, x, o);
    }
};
}  // namespace
SPIKE_EXPORT_ALL(teqp)
