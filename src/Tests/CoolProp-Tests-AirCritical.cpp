// Dense round-trip consistency coverage of Air's critical region.
//
// Air is modelled as a PSEUDO-PURE fluid, so its property surfaces are not
// perfectly continuous near the critical point (the underlying mixture has a
// temperature glide that the pseudo-pure approximation collapses).  That makes
// the flash routines fragile in this region; in particular the D+{H,S,U}
// back-flashes (HSU_D_flash) used to throw "p is not a valid number" at
// near-critical density because pseudo-pure fluids cannot take the
// is_pure()-gated superancillary happy path and fall onto the legacy
// saturation_D_pure path, which is unstable near rho_c.
//
// This harness sweeps a dense (T, p) grid spanning the critical region
// (T/Tc in [0.96, 2.0], p/pc in [0.5, 5.0]) and, for every single-phase grid
// point, requires that each of the seven non-trivial back-flashes recovers the
// forward (T, rho).  Grid points that land in the two-phase dome are skipped:
// PT inputs there are genuinely unsupported for a pseudo-pure fluid.
//
// Built into CatchTestRunner when COOLPROP_CATCH_MODULE=ON.  Run via
//   ./build_catch/CatchTestRunner "[air_critical]"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <memory>
#    include <cmath>
#    include <cstdio>
#    include <array>

TEST_CASE("Air critical-region round-trip consistency (pseudo-pure)", "[air_critical][HSU_D]") {
    using namespace CoolProp;

    auto fwd = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Air"));
    auto rt = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Air"));

    const double Tc = fwd->T_critical();
    const double pc = fwd->p_critical();

    // Round-trip tolerance.  The critical region is ill-conditioned (dp/drho -> 0
    // near rho_c), so a small relative slack is used; the production solvers stay
    // well inside it everywhere on this grid.
    const double tol = 2e-4;

    const int NT = 60, NP = 40;
    std::size_t single_phase = 0, two_phase_skipped = 0;

    for (int i = 0; i < NT; ++i) {
        const double T = (0.96 + (2.0 - 0.96) * i / (NT - 1)) * Tc;
        for (int j = 0; j < NP; ++j) {
            const double p = (0.5 + (5.0 - 0.5) * j / (NP - 1)) * pc;

            // Forward PT defines the reference single-phase state.  A throw here
            // means the point fell in/at the two-phase dome (PT not supported for
            // a pseudo-pure fluid) -- not a back-flash defect -- so skip it.
            try {
                fwd->update(PT_INPUTS, p, T);
            } catch (...) {
                ++two_phase_skipped;
                continue;
            }
            const double rho = fwd->rhomolar();
            const double h = fwd->hmolar();
            const double s = fwd->smolar();
            const double u = fwd->umolar();
            if (!ValidNumber(rho) || !ValidNumber(h) || !ValidNumber(s) || !ValidNumber(u)) {
                continue;
            }
            ++single_phase;

            // Each back-flash must recover the forward (T, rho).
            struct Case
            {
                const char* name;
                input_pairs pair;
                double v1, v2;
            };
            const std::array<Case, 7> cases = {{
              {"DmolarHmolar", DmolarHmolar_INPUTS, rho, h},
              {"DmolarSmolar", DmolarSmolar_INPUTS, rho, s},
              {"DmolarUmolar", DmolarUmolar_INPUTS, rho, u},
              {"HmolarP", HmolarP_INPUTS, h, p},
              {"PSmolar", PSmolar_INPUTS, p, s},
              {"PUmolar", PUmolar_INPUTS, p, u},
              {"DmolarP", DmolarP_INPUTS, rho, p},
            }};
            for (const auto& c : cases) {
                CAPTURE(c.name, T, p, T / Tc, p / pc, rho);
                bool threw = false;
                try {
                    rt->update(c.pair, c.v1, c.v2);
                } catch (const std::exception& e) {
                    threw = true;
                    FAIL_CHECK("back-flash threw: " << e.what());
                }
                if (threw) {
                    continue;
                }
                CHECK(rt->T() == Catch::Approx(T).epsilon(tol));
                CHECK(rt->rhomolar() == Catch::Approx(rho).epsilon(tol));
            }
        }
    }

    std::printf("[air_critical] Air: %zu single-phase points round-tripped (x7 back-flashes), %zu two-phase PT points skipped\n", single_phase,
                two_phase_skipped);
    // Guard against the grid silently collapsing (e.g. forward PT failing wholesale).
    CHECK(single_phase > 2000);
}

#endif  // ENABLE_CATCH
