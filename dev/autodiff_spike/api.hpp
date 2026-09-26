// Uniform interface every method TU exports.  Quantities (teqp conventions):
//   Ar0n : out = {A00, A01, A02, A03},  A0n = rho^n  d^n(alphar)/d rho^n
//   Arn0 : out = {A00, A10, A20},       An0 = (1/T)^n d^n(alphar)/d(1/T)^n
//   Ar11 : out = {A11}
//   gradx: out = {d alphar/d x_i}, i < N, compositions treated as independent
// A slot a method cannot produce is written as NaN.
#pragma once
#include <array>
#include <cstddef>
#include <limits>
#include <utility>
#include "pcsaft_generic.hpp"

namespace spike {
using XArr = std::array<double, NMAX>;
using Fn = void (*)(const PCSAFTParams&, double T, double rho, const XArr& x, double* out);
struct MethodTable
{
    const char* name;
    Fn Ar0n, Arn0, Ar11, gradx;
};
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
}  // namespace spike

#ifndef KMODELS
#    define KMODELS 1
#endif

// Instantiate Impl<Tag> for Tag = 0..KMODELS-1 and export the table.  Taking the
// addresses forces every instantiation to be emitted, which is what the compile-time
// scaling experiment measures.
#define SPIKE_EXPORT(NAME)                                                                                      \
    template <std::size_t... I>                                                                                 \
    static std::array<spike::MethodTable, KMODELS> make_table_(std::index_sequence<I...>) {                     \
        return {spike::MethodTable{#NAME, &Impl<I>::Ar0n, &Impl<I>::Arn0, &Impl<I>::Ar11, &Impl<I>::gradx}...}; \
    }                                                                                                           \
    spike::MethodTable method_##NAME(int k) {                                                                   \
        static const auto t = make_table_(std::make_index_sequence<KMODELS>{});                                 \
        return t[k];                                                                                            \
    }
