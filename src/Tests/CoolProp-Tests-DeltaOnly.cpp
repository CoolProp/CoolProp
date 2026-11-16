// Equivalence guard for the delta-only residual-alphar evaluation (all_deltaonly).
//
// The HEOS mixture-flash density-root solvers (spinodal search + pressure/dp-drho residuals) never
// need tau-derivatives, so they call ResidualHelmholtz::all_deltaonly(), which evaluates only the
// delta-derivatives of alpha^r (orders 0..4) and skips the tau/mixed/4th-order-tau work.  The whole
// speedup rests on those delta-derivatives agreeing with the full all().  This test locks that in:
// over a (tau, delta) grid, for a pure fluid (generalized-exponential term only), a GERG mixture
// (adds the departure/excess term) and a CO2/Water pair (adds the non-analytic term).
//
// ALL FIVE ORDERS ARE COMPARED TO 1e-12, SCALED PER FIELD.  No order is asserted bit-for-bit.
// The check is |got - ref| <= 1e-12 * max(|ref|, 1e-3 * M), where M is that field's largest
// magnitude over the grid for that fluid.
//
// Bit-exactness is not an invariant of this code; it was a property of particular builds.  The two
// functions compute every delta-field from textually identical expressions (compare
// src/Helmholtz.cpp ResidualHelmholtzGeneralizedExponential::all() / ::all_deltaonly() and the
// NonAnalytic pair), so there is no algebraic difference to find.  But they are two separately
// compiled function bodies, and neither the language nor the compilers promise that they round
// identically.  History:
//
//   * Orders 3-4 first differed on arm64 / Apple clang -O3: 1-3 ulp at tau = 0.96 for the
//     non-analytic fluids.  Rebuilding Helmholtz.cpp with -ffp-contract=off removed it, so FMA
//     contraction was the mechanism there.  (That build is bit-exact again on all orders today.)
//   * Order 2 first differed when GitHub's ubuntu-latest moved to ubuntu-26.04 (COO-182): 243
//     order-2 and 15 order-3/4 comparisons failed, across all five cases; orders 0 and 1 stayed
//     bit-exact.  The mechanism is FMA contraction again.  The runner image's GCC 15.2 defaults to
//     -march=x86-64-v3, so -mfma is on and Helmholtz.cpp.o carries ~570 FMA instructions under
//     GCC's default -ffp-contract=fast (gnu++17).  The gcc in a stock ubuntu:26.04 container, with
//     the same version string, defaults to plain x86-64 (no FMA) and gives bit-exact agreement,
//     which is why it did not reproduce there; adding -march=x86-64-v3 there reproduces it.  Measured on the runner (AMD EPYC 9V74) with the check below: worst
//     scaled error 4.1e-15 (CarbonDioxide, order 2); worst distance ~700 ulp (orders 3-4).
//
// None of the callers need bit equality.  calc_alphar_deriv_nocache() and the density-root
// residuals use all_deltaonly() for pressure and its delta-derivatives inside iterative solvers.
// Their convergence tolerances are many orders above 1e-12.  So the test asserts the property
// they rely on: the same function, to rounding.
//
// Why the floor 1e-3 * M.  A plain relative test has no fixed bound at a grid point where a
// derivative passes close to zero, because there the value is a small difference of large terms
// and its relative rounding error is unbounded.  The floor keeps the tolerance proportional to the
// size of the terms in that field instead.  1e-12 is ~4500 ulp, about 250x above the worst scaled
// error measured on the ubuntu-26.04 runner.  A real formula change still fails by many orders: perturbing one coefficient in
// all_deltaonly by 1e-7 relative fails 350 of 2865 assertions.  On failure the test captures the values,
// their ulp distance and the scaled error, so rounding noise and a genuine regression are easy to
// tell apart.
//
// Mirrors "Cubic shared-intermediate all() is bit-exact ..." in CoolProp-Tests-CubicU.cpp.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

#    include <algorithm>
#    include <array>
#    include <cmath>
#    include <cstddef>
#    include <cstdint>
#    include <cstring>
#    include <limits>
#    include <memory>
#    include <string>
#    include <utility>
#    include <vector>

using namespace CoolProp;

namespace {

// Agreement to within a few ULP, relative to the larger magnitude.  For two evaluations of the same quantity at the same
// inputs: an inlined and an out-of-line copy of the same code may still contract to FMA differently (COO-182), so the
// check is not bit-exact.  A state whose cached fields belong to another density than the one it reports (e.g. the
// iterate before the last solver step) is off by the relative size of that step, which this bound catches down to 1e-13.
bool same_to_rounding(CoolPropDbl a, CoolPropDbl b) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 1e-13 * std::max(std::abs(a), std::abs(b));
}

// Distance in units in the last place between two doubles, for the failure message only: the
// assertion itself is scaled-relative (see the header).  Maps the IEEE bit pattern onto a monotone
// integer line so that the distance is well-defined across zero and between adjacent binades.
// Non-finite inputs return UINT64_MAX; all_delta_fields_finite() fails those points separately.
std::uint64_t ulp_distance(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    auto ordered = [](double d) {
        std::uint64_t u = 0;
        std::memcpy(&u, &d, sizeof(double));
        // Negative values: flip all bits; non-negative: set the sign bit.  The result orders
        // like the doubles themselves, with -0.0 and +0.0 adjacent.
        return (u & 0x8000000000000000ULL) ? ~u : (u | 0x8000000000000000ULL);
    };
    const std::uint64_t ua = ordered(a), ub = ordered(b);
    return ua > ub ? ua - ub : ub - ua;
}

// The closeness check below fails on NaN, and on inf vs inf (inf - inf is NaN), but it compares
// after narrowing to double.  Confirm the values are real numbers separately, before narrowing.
bool all_delta_fields_finite(const HelmholtzDerivatives& d) {
    // No narrowing to double first: std::isfinite overloads on long double, and narrowing would
    // report a finite long-double value above DBL_MAX as non-finite if CoolPropDbl ever changes.
    return std::isfinite(d.alphar) && std::isfinite(d.dalphar_ddelta) && std::isfinite(d.d2alphar_ddelta2) && std::isfinite(d.d3alphar_ddelta3)
           && std::isfinite(d.d4alphar_ddelta4);
}

struct MixCase
{
    std::string fluids;
    std::vector<double> z;
};

}  // namespace

TEST_CASE("HEOS all_deltaonly agrees with full all() on the delta-derivatives", "[flash][mixture][deltaonly]") {
    const std::vector<MixCase> cases = {
      {"Methane", {1.0}},                                                            // generalized-exponential term only
      {"Methane&Ethane&Propane&n-Butane&Nitrogen", {0.80, 0.10, 0.05, 0.03, 0.02}},  // + GERG departure (excess term)
      {"Water", {1.0}},                                                              // non-analytic term (pure)
      {"CarbonDioxide", {1.0}},                                                      // non-analytic term (pure)
      {"CarbonDioxide&Water", {0.98, 0.02}},                                         // non-analytic term (mixture)
    };
    for (const auto& c : cases) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", c.fluids));
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        heos->set_mole_fractions(c.z);
        const std::vector<CoolPropDbl> z(c.z.begin(), c.z.end());

        // Every comparison below is a closeness test, and 0 is close to 0, so a regression that
        // zeroed a field in BOTH paths -- an
        // accumulator that is never written, a term container that stops dispatching -- would leave
        // every assertion for that field green.  Track the largest magnitude produced PER FIELD; one
        // maximum across fields would be dominated by d4 (by 2 to 5 orders of magnitude here) and
        // would not notice alphar, d1, d2 or d3 collapsing.
        std::array<CoolPropDbl, 5> max_abs{};

        // Pass 1: evaluate both paths over the grid.  The tolerance in pass 2 is scaled per field by
        // that field's largest magnitude on the grid, so it needs the whole grid first.
        struct Point
        {
            double tau, delta;
            std::array<double, 5> got, ref;
        };
        std::vector<Point> points;
        for (int it_tau = 0; it_tau <= 8; ++it_tau) {
            const double tau = 0.6 + 0.18 * it_tau;  // 0.6 .. 2.04
            for (int it_del = 0; it_del <= 8; ++it_del) {
                const double delta = 0.05 + 0.30 * it_del;  // 0.05 .. 2.45 (all > 0)
                CAPTURE(c.fluids, tau, delta);
                const HelmholtzDerivatives full = heos->residual_helmholtz->all(*heos, z, tau, delta, /*cache_values=*/false);
                const HelmholtzDerivatives donly = heos->residual_helmholtz->all_deltaonly(*heos, z, tau, delta);
                CHECK(all_delta_fields_finite(full));
                CHECK(all_delta_fields_finite(donly));
                const Point pt{
                  tau,
                  delta,
                  {static_cast<double>(donly.alphar), static_cast<double>(donly.dalphar_ddelta), static_cast<double>(donly.d2alphar_ddelta2),
                   static_cast<double>(donly.d3alphar_ddelta3), static_cast<double>(donly.d4alphar_ddelta4)},
                  {static_cast<double>(full.alphar), static_cast<double>(full.dalphar_ddelta), static_cast<double>(full.d2alphar_ddelta2),
                   static_cast<double>(full.d3alphar_ddelta3), static_cast<double>(full.d4alphar_ddelta4)}};
                for (std::size_t q = 0; q < 5; ++q) {
                    max_abs[q] = std::max(max_abs[q], static_cast<CoolPropDbl>(std::abs(pt.ref[q])));
                }
                points.push_back(pt);
            }
        }

        // Pass 2: |got - ref| <= 1e-12 * max(|ref|, 1e-3 * max_abs[field]).  See the header comment.
        // The floor stops a point where a derivative crosses zero (and so cancels heavily) from
        // demanding a relative accuracy that rounding cannot deliver.  It also covers ref == 0 exactly.
        const std::array<const char*, 5> names = {"alphar", "dalphar_ddelta", "d2alphar_ddelta2", "d3alphar_ddelta3", "d4alphar_ddelta4"};
        for (const Point& pt : points) {
            for (std::size_t q = 0; q < 5; ++q) {
                const double g = pt.got[q], r = pt.ref[q];
                const double scale = std::max(std::abs(r), 1e-3 * static_cast<double>(max_abs[q]));
                const double err = std::abs(g - r) / scale;
                const std::uint64_t ulps = ulp_distance(g, r);
                CAPTURE(c.fluids, pt.tau, pt.delta, names[q], g, r, ulps, err);
                // Written as !(err > tol) would pass a NaN; this form fails it.
                CHECK(err <= 1e-12);
            }
        }
        // Smallest per-field maximum measured over this grid is ~3e4 (Methane alphar), so 1e-12 is
        // 16 orders of margin: this cannot flake.  It is deliberately coarse -- it fires only when a
        // field is zero at EVERY one of the 81 grid points, not at a subset.
        // Non-finite values need no check here -- all_delta_fields_finite() rejects them at every
        // grid point, and std::max keeps the old maximum when the candidate is NaN.
        for (std::size_t q = 0; q < 5; ++q) {
            CAPTURE(c.fluids, q, max_abs[q]);
            REQUIRE(max_abs[q] > 1e-12);
        }
    }
}

// The warm-started density solve (solver_rho_Tp's Householder path, reached through update_TP_guessrho and
// PT_flash_with_guesses) iterates on the delta-only residual, which writes no state, and then puts the state on the
// density it returns; update_TP_guessrho skips its own update_DmolarT_direct when the state is already there (COO-126).
// Lock in what the callers see: the state is the one update_DmolarT_direct gives at that density -- p and all fifteen cached
// alphar derivatives bit for bit, so nothing is left cached from an iterate -- and that density is the root the cold PT
// flash finds.  (update_DmolarT_direct sets no phase, so phase-checked outputs such as hmolar of a pure fluid or the
// fugacity coefficients are not defined on this state, before or after the change; at fixed composition they are built
// from these alphar fields.)
TEST_CASE("Warm-started density solve leaves the state of update_DmolarT_direct at the root", "[flash][mixture][deltaonly]") {
    struct WarmCase
    {
        std::string backend, fluids;
        std::vector<double> z;
        double T, p;
    };
    const std::vector<WarmCase> cases = {
      {"HEOS", "Methane", {1.0}, 150, 5e6},                                                                  // pure, compressed liquid
      {"HEOS", "Methane", {1.0}, 300, 1e6},                                                                  // pure, gas
      {"GERG2008", "Methane&Ethane&Propane&n-Butane&Nitrogen", {0.80, 0.10, 0.05, 0.03, 0.02}, 300, 10e6},   // dense supercritical
      {"GERG2008", "Methane&Ethane&Propane&n-Butane&Nitrogen", {0.80, 0.10, 0.05, 0.03, 0.02}, 250, 0.5e6},  // gas
      {"HEOS", "CarbonDioxide&Water", {0.99, 0.01}, 450, 20e6},  // non-analytic terms, well clear of water dropout
    };
    using Acc = CoolPropDbl (AbstractState::*)();
    const std::vector<std::pair<const char*, Acc>> fields = {
      {"alphar", &AbstractState::alphar},
      {"dalphar_dDelta", &AbstractState::dalphar_dDelta},
      {"dalphar_dTau", &AbstractState::dalphar_dTau},
      {"d2alphar_dDelta2", &AbstractState::d2alphar_dDelta2},
      {"d2alphar_dDelta_dTau", &AbstractState::d2alphar_dDelta_dTau},
      {"d2alphar_dTau2", &AbstractState::d2alphar_dTau2},
      {"d3alphar_dDelta3", &AbstractState::d3alphar_dDelta3},
      {"d3alphar_dDelta2_dTau", &AbstractState::d3alphar_dDelta2_dTau},
      {"d3alphar_dDelta_dTau2", &AbstractState::d3alphar_dDelta_dTau2},
      {"d3alphar_dTau3", &AbstractState::d3alphar_dTau3},
      {"d4alphar_dDelta4", &AbstractState::d4alphar_dDelta4},
      {"d4alphar_dDelta3_dTau", &AbstractState::d4alphar_dDelta3_dTau},
      {"d4alphar_dDelta2_dTau2", &AbstractState::d4alphar_dDelta2_dTau2},
      {"d4alphar_dDelta_dTau3", &AbstractState::d4alphar_dDelta_dTau3},
      {"d4alphar_dTau4", &AbstractState::d4alphar_dTau4},
    };
    for (const auto& c : cases) {
        for (double guess_factor : {0.97, 1.03}) {
            CAPTURE(c.fluids, c.T, c.p, guess_factor);
            std::shared_ptr<AbstractState> ref(AbstractState::factory(c.backend, c.fluids));
            ref->set_mole_fractions(c.z);
            ref->update(PT_INPUTS, c.p, c.T);
            const double rho_cold = ref->rhomolar();
            REQUIRE(ref->phase() != iphase_twophase);

            auto check_state = [&](AbstractState& AS) {
                const double rho = AS.rhomolar();
                CHECK(std::abs(rho / rho_cold - 1) < 1e-9);  // the same root as the cold flash
                std::shared_ptr<AbstractState> fresh(AbstractState::factory(c.backend, c.fluids));
                fresh->set_mole_fractions(c.z);
                auto* fresh_heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(fresh.get());
                REQUIRE(fresh_heos != nullptr);
                fresh_heos->update_DmolarT_direct(rho, c.T);
                CHECK(std::abs(fresh->p() / c.p - 1) < 1e-8);  // and it is a root
                // Both states are written by update_DmolarT_direct at the same (T, rho, z); see same_to_rounding.
                CHECK(same_to_rounding(AS.p(), fresh->p()));
                for (const auto& f : fields) {
                    CAPTURE(f.first);
                    const CoolPropDbl a = (AS.*f.second)(), b = ((*fresh).*f.second)();
                    CHECK(same_to_rounding(a, b));
                }
            };
            {  // update_TP_guessrho
                std::shared_ptr<AbstractState> AS(AbstractState::factory(c.backend, c.fluids));
                AS->set_mole_fractions(c.z);
                auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
                REQUIRE(heos != nullptr);
                heos->update_TP_guessrho(c.T, c.p, guess_factor * rho_cold);
                check_state(*AS);
            }
            {  // update_with_guesses: PT_flash_with_guesses reads the density off the state
                std::shared_ptr<AbstractState> AS(AbstractState::factory(c.backend, c.fluids));
                AS->set_mole_fractions(c.z);
                GuessesStructure g;
                g.rhomolar = guess_factor * rho_cold;
                AS->update_with_guesses(PT_INPUTS, c.p, c.T, g);
                check_state(*AS);
            }
        }
    }
}

#endif
