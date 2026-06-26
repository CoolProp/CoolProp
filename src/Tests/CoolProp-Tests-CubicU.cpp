// Catch2 benchmarks + bit-exactness guard for the cubic-EOS forward eval u(T, rho).
//
// Motivation: when designing UV tabular inputs for SVDSBTL we want to know whether a
// tabular forward eval is actually faster than calling a cubic EOS (SRK/PR) directly, so
// we time u(T, rho) for the cubics against the full Helmholtz reference (HEOS).  This file
// replaces the old standalone dev/bench_cubic_u.cpp executable: the timing lives in a
// [!benchmark] opt-in test (run with `CatchTestRunner "[cubic][!benchmark]"`) and the
// bit-exactness check is a normal [cubic] test that runs in CI.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"
#    include "../Backends/Cubics/GeneralizedCubic.h"

#    include <cstdint>
#    include <cstring>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

// Reinterpret a double's storage as a uint64 so two values can be compared bit-for-bit.
// The shared-intermediate optimization claims a *term-for-term identical* assembly, not
// merely ~equal results, so the guard below must be exact (and memcpy avoids the
// suspicious-memory-comparison pitfall of memcmp on floating-point).
std::uint64_t bits_of(double d) {
    std::uint64_t u = 0;
    std::memcpy(&u, &d, sizeof(u));
    return u;
}

}  // namespace

// Opt-in via [!benchmark] (mirrors the [SVDComponents][!benchmark] timing tests).  Each cubic
// is timed twice: default (lets update_DmolarT do superancillary phase determination) and with
// the phase imposed (sends update_DmolarT down the direct-pressure branch, skipping that work).
// update(DmolarT) clears the property cache every call, so umolar() genuinely recomputes each
// iteration -- this is the true forward-eval cost, not a cache hit.
TEST_CASE("Cubic u(T,rho) forward-eval timing", "[cubic][!benchmark]") {
    const std::string fluid = "Methane";
    const double T = 300.0;           // K
    const double rho_molar = 5000.0;  // mol/m^3

    struct Case
    {
        const char* backend;
        bool impose_phase;
    };
    for (auto c : {Case{"SRK", false}, Case{"SRK", true}, Case{"PR", false}, Case{"PR", true}, Case{"HEOS", false}}) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(c.backend, fluid));
        if (c.impose_phase) {
            AS->specify_phase(iphase_gas);
        }
        AS->update(DmolarT_INPUTS, rho_molar, T);  // warm-up: pay one-time lazy-init costs
        const std::string label = std::string(c.backend) + (c.impose_phase ? " (phase)" : "");
        BENCHMARK(std::string("u(T,rho): ") + label) {
            AS->update(DmolarT_INPUTS, rho_molar, T);
            return AS->umolar();
        };
    }

    // Floor cost: just the single residual derivative umolar actually needs --
    // alphar(tau, delta, z, itau=1, idelta=0) called directly on the cubic object.
    for (const char* backend : {"SRK", "PR"}) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluid));
        auto* cb = dynamic_cast<AbstractCubicBackend*>(AS.get());
        REQUIRE(cb != nullptr);
        auto& cubic = cb->get_cubic();
        const std::vector<double> z = {1.0};
        const double tau = cubic->get_Tr() / T;
        const double delta = rho_molar / cubic->get_rhor();
        BENCHMARK(std::string("raw alphar itau=1: ") + backend) {
            return cubic->alphar(tau, delta, z, 1, 0);
        };
    }
}

// Correctness guard (runs in CI): the shared-intermediate assembly in the *production*
// CubicResidualHelmholtz::all() must be bit-for-bit identical to 15 independent alphar() calls
// (alphar() itself is unchanged, so it is the reference).  This drives the real all() -- including
// its 15 hand-unrolled derivative assignments -- so a transcription typo in any one line is caught,
// not just the algebraic identity it relies on.  The 15 (itau, idelta) derivatives are exactly the
// pairs with itau + idelta <= 4.
TEST_CASE("Cubic shared-intermediate all() is bit-exact vs 15 alphar() calls", "[cubic]") {
    for (const char* backend : {"SRK", "PR"}) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, "Methane"));
        auto* cb = dynamic_cast<AbstractCubicBackend*>(AS.get());
        REQUIRE(cb != nullptr);
        auto& cubic = cb->get_cubic();
        CubicResidualHelmholtz rh(cb);
        const std::vector<CoolPropDbl> z = {1.0};
        for (int it_tau = 0; it_tau <= 10; ++it_tau) {
            const double tau = 0.7 + 0.13 * it_tau;
            for (int it_del = 0; it_del <= 10; ++it_del) {
                const double delta = 0.2 + 0.27 * it_del;
                // Production assembly under test (the optimized 15-derivative path).
                HelmholtzDerivatives assembled = rh.all(*cb, z, tau, delta, /*cache_values=*/false);
                for (std::size_t it = 0; it <= 4; ++it) {
                    for (std::size_t id = 0; id + it <= 4; ++id) {
                        const double direct = cubic->alphar(tau, delta, z, it, id);
                        CHECK(bits_of(assembled.get(it, id)) == bits_of(direct));
                    }
                }
            }
        }
    }
}

#endif
