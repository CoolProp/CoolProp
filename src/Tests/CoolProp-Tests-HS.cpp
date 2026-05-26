// Binding regression test for the ACTIVE production H,S flash.
//
// This exercises the shipped flash directly via HmolarSmolar_INPUTS (no
// prototype solver) and is the only HS test that runs in the default suite.
// The exploratory prototype solvers and their [.]-hidden A/B characterization
// tests live in CoolProp-Tests-HS-prototypes.cpp; the active production code is
// in src/Backends/Helmholtz/FlashRoutines.cpp.

#include "AbstractState.h"
#include "CoolProp.h"
#include "DataStructures.h"

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
}

#endif  // ENABLE_CATCH
