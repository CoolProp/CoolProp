// Equivalence guard for the delta-only residual-alphar evaluation (all_deltaonly).
//
// The HEOS mixture-flash density-root solvers (spinodal search + pressure/dp-drho residuals) never
// need tau-derivatives, so they call ResidualHelmholtz::all_deltaonly(), which evaluates only the
// delta-derivatives of alpha^r (orders 0..4) and skips the tau/mixed/4th-order-tau work.  The whole
// speedup rests on those delta-derivatives agreeing with the full all().  This test locks that in:
// over a (tau, delta) grid, for a pure fluid (generalized-exponential term only), a GERG mixture
// (adds the departure/excess term) and a CO2/Water pair (adds the non-analytic term).
//
// ORDERS 0-2 ARE COMPARED BIT-FOR-BIT.  ORDERS 3-4 ARE COMPARED TO 1e-14 RELATIVE, and that is a
// deliberate weakening of an assertion that was previously bit-exact and FAILED ON MASTER for the
// three fluids carrying the non-analytic term (Water, CO2, CO2&Water), at tau = 0.96 only.
//
// Why the weakening is right rather than a papering-over.  The two functions compute these two
// fields from *textually identical* expressions -- compare src/Helmholtz.cpp
// ResidualHelmholtzNonAnalytic::all() and ::all_deltaonly(): every intermediate (theta, PSI, DELTA
// and their delta-derivatives) and both final accumulations are the same source text.  There is no
// algebraic difference to find.  What differs is code generation: all() also computes the tau and
// mixed derivatives, so it carries far more live values, and the compiler makes different
// FMA-contraction and register-allocation choices in the two bodies.  Measured discrepancy:
//
//     Water        tau=0.96 delta=0.65   d3: 1 ulp (2.1e-16 rel)
//     CarbonDioxide tau=0.96 delta=0.05  d3: 3 ulp (3.9e-16 rel)
//     ... 15 failing comparisons in all, every one of them 1-3 ulp
//
// i.e. the last bit of a double, showing through only where the non-analytic term's high-order
// delta-derivatives cancel heavily near the critical point (tau = 0.96 is the only tau on the grid
// where it happens).  Bit-exactness across two separately-compiled function bodies is not something
// the language guarantees, so asserting it was asserting a property of one compiler's output.
//
// 1e-14 is ~45 ulp: far above the observed 1-3, and far below anything a genuine formula change
// could produce -- editing an expression in all_deltaonly moves these values by parts in 1e3 or
// more, not parts in 1e14.  The regression this test exists to catch is still caught.
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
#    include <cstring>
#    include <memory>
#    include <string>
#    include <utility>
#    include <vector>

using namespace CoolProp;

namespace {

// Reinterpret a value's full storage as a byte array so two values compare bit-for-bit.  Keyed on
// CoolPropDbl (not double): when CoolPropDbl is long double, truncating to double would discard the
// extended-precision bits and hide a real difference.  memcpy avoids the suspicious-memory-comparison
// pitfall of memcmp on floating point; std::array gives a well-defined element-wise operator==.
std::array<unsigned char, sizeof(CoolPropDbl)> bits_of(CoolPropDbl d) {
    std::array<unsigned char, sizeof(CoolPropDbl)> u{};
    std::memcpy(u.data(), &d, sizeof(CoolPropDbl));
    return u;
}

// Both comparisons below treat "identical on each side" as agreement, so they can pass on garbage:
// two NaNs with the same payload are bitwise equal (orders 0-2), and WithinAbs(NaN, 1e-300) is
// false but bits_of(NaN) == bits_of(NaN) is true.  Confirm the values are real numbers separately.
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

        // Every comparison below is an equality or a relative-closeness test, and WithinAbs(0,
        // 1e-300) makes 0 == 0 pass, so a regression that zeroed a field in BOTH paths -- an
        // accumulator that is never written, a term container that stops dispatching -- would leave
        // every assertion for that field green.  Track the largest magnitude produced PER FIELD; one
        // maximum across fields would be dominated by d4 (by 2 to 5 orders of magnitude here) and
        // would not notice alphar, d1, d2 or d3 collapsing.
        std::array<CoolPropDbl, 5> max_abs{};

        for (int it_tau = 0; it_tau <= 8; ++it_tau) {
            const double tau = 0.6 + 0.18 * it_tau;  // 0.6 .. 2.04
            for (int it_del = 0; it_del <= 8; ++it_del) {
                const double delta = 0.05 + 0.30 * it_del;  // 0.05 .. 2.45 (all > 0)
                CAPTURE(c.fluids, tau, delta);
                HelmholtzDerivatives full = heos->residual_helmholtz->all(*heos, z, tau, delta, /*cache_values=*/false);
                HelmholtzDerivatives donly = heos->residual_helmholtz->all_deltaonly(*heos, z, tau, delta);
                CHECK(all_delta_fields_finite(full));
                CHECK(all_delta_fields_finite(donly));
                // Orders 0-2: bit-for-bit, which holds for every case on the grid.
                CHECK(bits_of(donly.alphar) == bits_of(full.alphar));
                CHECK(bits_of(donly.dalphar_ddelta) == bits_of(full.dalphar_ddelta));
                CHECK(bits_of(donly.d2alphar_ddelta2) == bits_of(full.d2alphar_ddelta2));
                // Orders 3-4: 1e-14 relative.  See the header comment for the measurement and for
                // why bit-exactness is not assertable here.  WithinAbs is OR-ed in because
                // WithinRel is meaningless when the reference value is exactly zero, which happens
                // on this grid for the terms that vanish at delta = 1.
                CHECK_THAT((double)donly.d3alphar_ddelta3, Catch::Matchers::WithinRel((double)full.d3alphar_ddelta3, 1e-14)
                                                             || Catch::Matchers::WithinAbs((double)full.d3alphar_ddelta3, 1e-300));
                CHECK_THAT((double)donly.d4alphar_ddelta4, Catch::Matchers::WithinRel((double)full.d4alphar_ddelta4, 1e-14)
                                                             || Catch::Matchers::WithinAbs((double)full.d4alphar_ddelta4, 1e-300));

                const std::array<CoolPropDbl, 5> f = {full.alphar, full.dalphar_ddelta, full.d2alphar_ddelta2, full.d3alphar_ddelta3,
                                                      full.d4alphar_ddelta4};
                for (std::size_t q = 0; q < 5; ++q) {
                    max_abs[q] = std::max(max_abs[q], std::abs(f[q]));
                }
            }
        }
        // Smallest per-field maximum measured over this grid is ~3e4 (Methane alphar), so 1e-12 is
        // 16 orders of margin: this cannot flake.  It is deliberately coarse -- it fires only when a
        // field is zero at EVERY one of the 81 grid points, not at a subset.
        // Non-finite values need no check here -- all_delta_fields_finite() rejects them at every
        // grid point, and std::max discards a NaN candidate rather than propagating it.
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
                CHECK(bits_of(AS.p()) == bits_of(fresh->p()));
                for (const auto& f : fields) {
                    CAPTURE(f.first);
                    const CoolPropDbl a = (AS.*f.second)(), b = ((*fresh).*f.second)();
                    CHECK(std::isfinite(a));
                    CHECK(bits_of(a) == bits_of(b));
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
