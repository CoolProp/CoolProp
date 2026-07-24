// Catch2 benchmark for the cubic-EOS MIXTURE flash forward eval (GitHub #3192 / PR #3224).
//
// The #3192 cubic-mixture perf work targets the residual-Helmholtz evaluations that dominate a
// multicomponent cubic (SRK/PR) flash:
//   * geometric-mean collapse of the attractive mixture term a_m  (O(N^2) -> O(N) per evaluation),
//   * a per-tau a_ij matrix cache read directly by the composition gradient/Hessian
//     (d_am_term_dxi / d2_am_term_dxidxj), removing the per-element cache-validity checks, and
//   * a composition-generation cache for the covolume mixture term bm.
// A PT flash drives the density solve (am_term / bm via all()); a two-phase QT flash additionally
// drives the fugacity-coefficient composition derivatives (d_alphar_dxi -> d_am_term_dxi), which is
// where the a_ij matrix cache pays off most.  Both are timed here.
//
// Opt-in [!benchmark] (run with `CatchTestRunner "[cubic][mixture][!benchmark]"`), mirroring the
// pure-fluid CoolProp-Tests-CubicU.cpp.  Nothing here runs in CI; correctness is covered by the
// [consistency] suite and the kij-invalidation guard in CoolProp-Tests.cpp.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

TEST_CASE("Cubic mixture flash forward-eval timing (#3192)", "[cubic][mixture][!benchmark]") {
    // 5-component natural gas.
    const std::vector<std::string> comps = {"Methane", "Ethane", "Propane", "n-Butane", "Nitrogen"};
    const std::vector<double> z = {0.80, 0.10, 0.05, 0.03, 0.02};

    for (const char* backend : {"SRK", "PR"}) {
        auto make = [&]() {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, comps));
            AS->set_mole_fractions(z);
            return AS;
        };

        // --- Single-phase PT flash: 320 K / 20 bar (blind update(PT_INPUTS)). ---
        {
            const double T = 320.0, p = 20e5;
            auto AS = make();
            AS->update(PT_INPUTS, p, T);  // warm-up: pay one-time lazy-init
            BENCHMARK(std::string("PT flash 5-comp NG 320K/20bar: ") + backend) {
                AS->update(PT_INPUTS, p, T);
                return AS->rhomolar();
            };
        }

        // --- Two-phase flash: bubble/dew at Q = 0.5, T = 220 K (drives the composition derivatives). ---
        {
            const double T = 220.0, Q = 0.5;
            auto AS = make();
            AS->update(QT_INPUTS, Q, T);  // warm-up
            BENCHMARK(std::string("QT flash 5-comp NG Q=0.5 220K: ") + backend) {
                AS->update(QT_INPUTS, Q, T);
                return AS->p();
            };
        }
    }
}

#endif
