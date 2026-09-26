// Interface for the "whole half-matrix in one call" experiment.
//   out[idx(i,j)] = A_ij = (1/T)^i rho^j d^{i+j} alphar / d(1/T)^i d rho^j,  i + j <= N,
// packed by total degree as in nd::Taylor2::idx; 15 slots for N = 4.
#pragma once
#include "api.hpp"

namespace spike {
inline constexpr int NMAXORDER = 4;
inline constexpr int NSLOTS = (NMAXORDER + 1) * (NMAXORDER + 2) / 2;
constexpr int idx2(int i, int j) {
    return (i + j) * (i + j + 1) / 2 + j;
}
struct MethodAll
{
    const char* name;
    std::array<Fn, NMAXORDER + 1> byN;  // byN[N] fills orders <= N (NaN above); nullptr if not built
};
inline void fill_nan(double* o) {
    for (int m = 0; m < NSLOTS; ++m)
        o[m] = NaN;
}
}  // namespace spike

// ONLY_N restricts instantiation to one order (used by the compile-time sweep)
#ifdef ONLY_N
#    define SPIKE_HAS_N(n) ((n) == ONLY_N)
#else
#    define SPIKE_HAS_N(n) true
#endif

#define SPIKE_EXPORT_ALL(NAME)                                                                                            \
    template <std::size_t I, int N>                                                                                       \
    static constexpr spike::Fn pick_() {                                                                                  \
        if constexpr (N >= 1 && SPIKE_HAS_N(N))                                                                           \
            return &Impl<I>::template all<N>;                                                                             \
        else                                                                                                              \
            return nullptr;                                                                                               \
    }                                                                                                                     \
    template <std::size_t... I>                                                                                           \
    static std::array<spike::MethodAll, KMODELS> make_table_(std::index_sequence<I...>) {                                 \
        return {spike::MethodAll{#NAME, {pick_<I, 0>(), pick_<I, 1>(), pick_<I, 2>(), pick_<I, 3>(), pick_<I, 4>()}}...}; \
    }                                                                                                                     \
    spike::MethodAll methodall_##NAME(int k) {                                                                            \
        static const auto t = make_table_(std::make_index_sequence<KMODELS>{});                                           \
        return t[k];                                                                                                      \
    }
