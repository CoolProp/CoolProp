// Reconnaissance harness for the resurrected `allEigen` evaluator.
//
// QUESTIONS:
//   1. Does `allEigen` produce numerically equivalent derivatives to scalar `all()`
//      for fluids where it should (no tau_mi_in_u branch, no 4th-order required)?
//   2. By how much, on real per-fluid term tables?
//   3. Is the 4th-order absence + tau_mi_in_u dead-comment fixable in a follow-up,
//      or do they make allEigen a non-starter as a drop-in?
//
// PART 1 — Equivalence: for each test fluid + (tau, delta) point, call both `all()`
//   and `allEigen()` on the same instance, with a zero-initialized HelmholtzDerivatives.
//   Compare every 0th-3rd-order field at tol = 1e-12 relative.  Skip R125/Methanol
//   (known tau_mi_in_u gap — allEigen will miss those terms).
//
// PART 2 — Microbench: for Water, R134a, n-Propane: 1e5 reps, 5 trials each, mean
//   ± stddev for scalar all() vs allEigen(), report ratio.
//
// This file is RECONNAISSANCE — Catch2 hidden ([.]).  Run via:
//   ./build_catch_rel/CatchTestRunner "[alleigen_recon]"

#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Helmholtz.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <chrono>
#    include <cmath>
#    include <cstdio>
#    include <memory>
#    include <numeric>
#    include <string>
#    include <vector>

namespace {

struct EquivResult
{
    bool ok = false;         // false => fluid skipped (factory fail / not HEOS pure)
    double worst_rel = 0.0;  // worst |delta| / max(|scalar|, 1.0) across all 0-3rd-order fields
    std::string worst_field;
    std::size_t N = 0;
};

double rel_diff(double a, double b) {
    const double scale = std::max(std::abs(a), 1.0);
    return std::abs(a - b) / scale;
}

EquivResult compare_at(CoolProp::ResidualHelmholtzGeneralizedExponential& gen, double tau, double delta) {
    EquivResult r;
    CoolProp::HelmholtzDerivatives scalar_d, eigen_d;
    // Zero-init both so the accumulating += inside all/allEigen lands on a known baseline.
    // The default ctor doesn't zero every field per the clang-tidy diagnostic, so do it explicitly.
    scalar_d = CoolProp::HelmholtzDerivatives{};
    eigen_d = CoolProp::HelmholtzDerivatives{};
    gen.all(static_cast<CoolPropDbl>(tau), static_cast<CoolPropDbl>(delta), scalar_d);
    gen.allEigen(static_cast<CoolPropDbl>(tau), static_cast<CoolPropDbl>(delta), eigen_d);

    // Compare 0th-3rd-order fields (allEigen does not write 4th-order, so excluded).
    struct Field
    {
        const char* name;
        double s, e;
    };
    const Field fields[] = {
      {"alphar", scalar_d.alphar, eigen_d.alphar},
      {"dalphar_ddelta", scalar_d.dalphar_ddelta, eigen_d.dalphar_ddelta},
      {"dalphar_dtau", scalar_d.dalphar_dtau, eigen_d.dalphar_dtau},
      {"d2alphar_ddelta2", scalar_d.d2alphar_ddelta2, eigen_d.d2alphar_ddelta2},
      {"d2alphar_dtau2", scalar_d.d2alphar_dtau2, eigen_d.d2alphar_dtau2},
      {"d2alphar_ddelta_dtau", scalar_d.d2alphar_ddelta_dtau, eigen_d.d2alphar_ddelta_dtau},
      {"d3alphar_ddelta3", scalar_d.d3alphar_ddelta3, eigen_d.d3alphar_ddelta3},
      {"d3alphar_dtau3", scalar_d.d3alphar_dtau3, eigen_d.d3alphar_dtau3},
      {"d3alphar_ddelta2_dtau", scalar_d.d3alphar_ddelta2_dtau, eigen_d.d3alphar_ddelta2_dtau},
      {"d3alphar_ddelta_dtau2", scalar_d.d3alphar_ddelta_dtau2, eigen_d.d3alphar_ddelta_dtau2},
    };
    for (const auto& f : fields) {
        const double rd = rel_diff(f.s, f.e);
        if (rd > r.worst_rel) {
            r.worst_rel = rd;
            r.worst_field = f.name;
        }
    }
    r.ok = true;
    return r;
}

CoolProp::ResidualHelmholtzGeneralizedExponential* get_genexp(CoolProp::AbstractState* AS) {
    auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS);
    if (be == nullptr) return nullptr;
    auto& comps = be->get_components();
    if (comps.empty()) return nullptr;
    return &comps[0].EOS().alphar.GenExp;
}

}  // namespace

TEST_CASE("allEigen vs scalar all() numerical equivalence (0th-3rd order, no tau_mi_in_u)", "[alleigen_recon][.]") {
    // Coverage: fluids without tau_mi_in_u branch — should be bit-equivalent within ULP.
    const std::vector<std::string> fluids = {"Water", "CarbonDioxide", "n-Propane", "R134a", "Methane", "Nitrogen", "Hydrogen", "MM"};

    // (tau, delta) test points covering different EOS regimes.
    struct TP
    {
        double tau, delta;
        const char* desc;
    };
    const std::vector<TP> points = {
      {2.0, 2.5, "compressed_liquid"}, {0.8, 0.3, "supercritical"},       {1.5, 1.0, "near_critical"},
      {3.0, 0.01, "dilute_gas"},       {1.05, 1.05, "near_crit_density"},
    };

    std::printf("\n=== allEigen vs all() equivalence (tol=1e-12 relative) ===\n");
    std::printf("%-15s %-20s %15s %s\n", "fluid", "regime", "worst_rel", "worst_field");
    std::printf("%-15s %-20s %15s %s\n", "-----", "------", "---------", "-----------");
    for (const auto& name : fluids) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", name));
        } catch (...) {
            std::printf("%-15s SKIPPED (factory failed)\n", name.c_str());
            continue;
        }
        auto* gen = get_genexp(AS.get());
        if (gen == nullptr) {
            std::printf("%-15s SKIPPED (no GenExp)\n", name.c_str());
            continue;
        }
        for (const auto& p : points) {
            const auto r = compare_at(*gen, p.tau, p.delta);
            if (!r.ok) continue;
            std::printf("%-15s %-20s %15.3e %s\n", name.c_str(), p.desc, r.worst_rel, r.worst_field.c_str());
            // Equivalence requirement: 1e-12 relative.
            CHECK(r.worst_rel < 1e-12);
        }
    }
}

TEST_CASE("allEigen known-gap fluids (R125, Methanol use tau_mi_in_u — divergence expected)", "[alleigen_recon][.]") {
    // R125 and Methanol exercise the tau_mi_in_u branch that allEigen has commented out.
    // We expect noticeable divergence — measure how big so the verdict can quote it.
    const std::vector<std::string> fluids = {"R125", "Methanol"};
    std::printf("\n=== known-gap fluids: tau_mi_in_u dead-coded in allEigen ===\n");
    for (const auto& name : fluids) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", name));
        } catch (...) {
            std::printf("%-15s SKIPPED\n", name.c_str());
            continue;
        }
        auto* gen = get_genexp(AS.get());
        if (gen == nullptr) continue;
        const auto r = compare_at(*gen, 1.5, 1.0);
        std::printf("%-15s near_critical: worst_rel = %.3e on %s\n", name.c_str(), r.worst_rel, r.worst_field.c_str());
    }
}

TEST_CASE("Microbench: scalar all() vs allEigen() per-call wall time", "[alleigen_recon][.]") {
    struct BenchFluid
    {
        std::string name;
        double tau, delta;
    };
    // Representative single-phase points well away from saturation/criticality.
    const std::vector<BenchFluid> benches = {
      {"Water", 1.4, 1.1},
      {"R134a", 1.2, 0.5},
      {"n-Propane", 1.3, 0.6},
      {"MM", 1.5, 0.8},
    };
    constexpr std::size_t reps = 100'000;
    constexpr std::size_t n_trials = 5;

    std::printf("\n=== Microbench: %zu reps x %zu trials ===\n", reps, n_trials);
    std::printf("%-12s %5s   scalar (ns/call)    eigen (ns/call)     ratio\n", "fluid", "N");
    std::printf("%-12s %5s   ---------------     ---------------     -----\n", "-----", "-");

    for (const auto& bf : benches) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bf.name));
        } catch (...) {
            continue;
        }
        auto* gen = get_genexp(AS.get());
        if (gen == nullptr) continue;

        std::vector<double> scalar_ns, eigen_ns;
        for (std::size_t trial = 0; trial < n_trials; ++trial) {
            // scalar
            volatile double scalar_sink = 0;
            const auto t0 = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < reps; ++i) {
                CoolProp::HelmholtzDerivatives d{};
                gen->all(static_cast<CoolPropDbl>(bf.tau), static_cast<CoolPropDbl>(bf.delta), d);
                scalar_sink = scalar_sink + d.alphar;
            }
            const auto t1 = std::chrono::steady_clock::now();
            scalar_ns.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(reps));
            (void)scalar_sink;

            // eigen
            volatile double eigen_sink = 0;
            const auto t2 = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < reps; ++i) {
                CoolProp::HelmholtzDerivatives d{};
                gen->allEigen(static_cast<CoolPropDbl>(bf.tau), static_cast<CoolPropDbl>(bf.delta), d);
                eigen_sink = eigen_sink + d.alphar;
            }
            const auto t3 = std::chrono::steady_clock::now();
            eigen_ns.push_back(std::chrono::duration<double, std::nano>(t3 - t2).count() / static_cast<double>(reps));
            (void)eigen_sink;
        }

        auto stats = [](const std::vector<double>& xs) {
            double m = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
            double v = 0.0;
            for (double x : xs)
                v += (x - m) * (x - m);
            return std::pair<double, double>{m, std::sqrt(v / xs.size())};
        };
        auto [sm, ss] = stats(scalar_ns);
        auto [em, es] = stats(eigen_ns);
        std::printf("%-12s %5zu   %7.1f +- %5.1f    %7.1f +- %5.1f    %.3f\n", bf.name.c_str(), gen->elements.size(), sm, ss, em, es, sm / em);
    }
}

#endif  // ENABLE_CATCH
