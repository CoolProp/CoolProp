// Bit-exactness guard for the delta-only residual-alphar evaluation (all_deltaonly).
//
// The HEOS mixture-flash density-root solvers (spinodal search + pressure/dp-drho residuals) never
// need tau-derivatives, so they call ResidualHelmholtz::all_deltaonly(), which evaluates only the
// delta-derivatives of alpha^r (orders 0..4) and skips the tau/mixed/4th-order-tau work.  The whole
// speedup rests on those delta-derivatives being *identical* to the full all() -- not merely ~equal.
// This test locks that invariant in: over a (tau, delta) grid, for a pure fluid (generalized-
// exponential term only), a GERG mixture (adds the departure/excess term) and a CO2/Water pair
// (adds the non-analytic term), every delta-field must be bit-for-bit equal to the full path.
//
// Mirrors "Cubic shared-intermediate all() is bit-exact ..." in CoolProp-Tests-CubicU.cpp.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

#    include <cstdint>
#    include <cstring>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

// Reinterpret a double's storage as a uint64 so two values compare bit-for-bit (memcpy avoids the
// suspicious-memory-comparison pitfall of memcmp on floating-point).
std::uint64_t bits_of(double d) {
    std::uint64_t u = 0;
    std::memcpy(&u, &d, sizeof(u));
    return u;
}

struct MixCase
{
    std::string fluids;
    std::vector<double> z;
};

}  // namespace

TEST_CASE("HEOS all_deltaonly is bit-exact vs full all() on the delta-derivatives", "[flash][mixture][deltaonly]") {
    const std::vector<MixCase> cases = {
      {"Methane", {1.0}},                                                            // generalized-exponential term only
      {"Methane&Ethane&Propane&n-Butane&Nitrogen", {0.80, 0.10, 0.05, 0.03, 0.02}},  // + GERG departure (excess term)
      {"CarbonDioxide&Water", {0.98, 0.02}},                                         // + non-analytic term
    };
    for (const auto& c : cases) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", c.fluids));
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        heos->set_mole_fractions(c.z);
        const std::vector<CoolPropDbl> z(c.z.begin(), c.z.end());

        for (int it_tau = 0; it_tau <= 8; ++it_tau) {
            const double tau = 0.6 + 0.18 * it_tau;  // 0.6 .. 2.04
            for (int it_del = 0; it_del <= 8; ++it_del) {
                const double delta = 0.05 + 0.30 * it_del;  // 0.05 .. 2.45 (all > 0)
                CAPTURE(c.fluids, tau, delta);
                HelmholtzDerivatives full = heos->residual_helmholtz->all(*heos, z, tau, delta, /*cache_values=*/false);
                HelmholtzDerivatives donly = heos->residual_helmholtz->all_deltaonly(*heos, z, tau, delta);
                CHECK(bits_of(donly.alphar) == bits_of(full.alphar));
                CHECK(bits_of(donly.dalphar_ddelta) == bits_of(full.dalphar_ddelta));
                CHECK(bits_of(donly.d2alphar_ddelta2) == bits_of(full.d2alphar_ddelta2));
                CHECK(bits_of(donly.d3alphar_ddelta3) == bits_of(full.d3alphar_ddelta3));
                CHECK(bits_of(donly.d4alphar_ddelta4) == bits_of(full.d4alphar_ddelta4));
            }
        }
    }
}

#endif
