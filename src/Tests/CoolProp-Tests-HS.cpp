// Binding regression test for the ACTIVE production H,S flash.
//
// This exercises the shipped flash directly via HmolarSmolar_INPUTS (no
// prototype solver) and is the only HS test that runs in the default suite.
// The exploratory prototype solvers and their [.]-hidden A/B characterization
// tests live in CoolProp-Tests-HS-prototypes.cpp; the active production code is
// in src/Backends/Helmholtz/FlashRoutines.cpp.

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <cstdio>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = (n == 1) ? a : a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return v;
}

}  // namespace

// ---------------------------------------------------------------------------
// PRODUCTION two-phase HS round-trip: exercises the new superancillary
// two-phase pre-screen + EOS-exact HS_flash_twophase path in FlashRoutines.
//     CatchTestRunner "[HS_prod2ph]"
// ---------------------------------------------------------------------------
TEST_CASE("HS production two-phase round-trip", "[HS][HS_prod2ph]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "R134a", "MM", "n-Pentane", "Methane");
    CAPTURE(fluid);
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tc = ref->T_critical(), Tmin = ref->Ttriple();
    std::size_t total = 0, ok = 0;
    for (double T : linspace(Tmin + 1.0, Tc - std::max(0.5, 1e-3 * Tc), 25)) {
        for (double Q : {0.05, 0.25, 0.5, 0.75, 0.95}) {
            try {
                ref->update(CoolProp::QT_INPUTS, Q, T);
            } catch (...) {
                continue;
            }
            const double h = ref->hmolar(), s = ref->smolar();
            if (!std::isfinite(h) || !std::isfinite(s)) continue;
            ++total;
            CAPTURE(T, Q, h, s);
            try {
                wrk->update(CoolProp::HmolarSmolar_INPUTS, h, s);
            } catch (const std::exception& e) {
                FAIL_CHECK("two-phase HS threw: " << e.what());
                continue;
            }
            const bool good = std::abs(wrk->T() - T) / T < 1e-5 && std::abs(wrk->Q() - Q) < 1e-4;
            if (good)
                ++ok;
            else {
                CAPTURE(wrk->T(), wrk->Q());
                CHECK(good);
            }
        }
    }
    std::printf("[HS_prod2ph] %s: %zu/%zu two-phase round-trips\n", fluid.c_str(), ok, total);
    REQUIRE(total > 0);  // grid actually exercised the flash
    CHECK(ok == total);  // every two-phase (h,s) round-tripped
}

// ---------------------------------------------------------------------------
// PRODUCTION single-phase HS round-trip across the whole phase diagram (liquid,
// vapor, supercritical).  Each state is established with PT_INPUTS (always
// single-phase), then (h,s) is round-tripped back to (T,rho) through the
// production HmolarSmolar_INPUTS flash; the originating T and rho must return.
// This is the single-phase complement to [HS_prod2ph]; together they cover
// general inputs over the full phase diagram.
//     CatchTestRunner "[HS_prod1ph]"
// ---------------------------------------------------------------------------
TEST_CASE("HS production single-phase round-trip (phase-diagram sweep)", "[HS][HS_prod1ph]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "CarbonDioxide", "R134a", "MM", "n-Pentane", "Methane", "Hydrogen");
    CAPTURE(fluid);
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    double plo;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-6);
    } catch (...) {
        plo = pmax * 1e-6;
    }
    std::size_t total = 0, ok = 0;
    // (log p, linear T) grid: log-p spans the compressed-liquid / dense
    // supercritical region down to dilute vapor; T spans triple to Tmax.
    for (double T : linspace(Tmin + 0.5, Tmax, 20)) {
        for (std::size_t j = 0; j < 20; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / 19.0);
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);  // always single-phase
            } catch (...) {
                continue;  // unreachable (p, T) for this fluid
            }
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar(), Tt = ref->T();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) continue;
            ++total;
            CAPTURE(T, p, rho, h, s);
            try {
                wrk->update(CoolProp::HmolarSmolar_INPUTS, h, s);
            } catch (const std::exception& e) {
                FAIL_CHECK("single-phase HS threw: " << e.what());
                continue;
            }
            const bool good = std::abs(wrk->T() - Tt) / Tt < 1e-5 && std::abs(wrk->rhomolar() - rho) / rho < 1e-5;
            if (good)
                ++ok;
            else {
                CAPTURE(wrk->T(), wrk->rhomolar());
                CHECK(good);
            }
        }
    }
    std::printf("[HS_prod1ph] %s: %zu/%zu single-phase round-trips\n", fluid.c_str(), ok, total);
    REQUIRE(total > 0);  // grid actually exercised the flash
    CHECK(ok == total);  // every single-phase (h,s) round-tripped
}

#endif  // ENABLE_CATCH
