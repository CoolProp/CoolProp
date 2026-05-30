// Comprehensive baseline microbench for ResidualHelmholtzGeneralizedExponential::all
// across a representative cross-section of the CoolProp fluid library.
//
// PURPOSE: establish the per-call scalar wall time we're trying to beat.  All
// SIMD/vectorization experiments on this branch compare against the numbers
// produced here.
//
// Run via: ./build_catch_rel/CatchTestRunner "[helmholtz_inner_bench]"
//
// Output CSV at /tmp/helmholtz_inner_bench.csv with one row per fluid:
//   fluid,N,scalar_ns_mean,scalar_ns_std
//
// Hidden ([.]) — does not run in CI default.

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

// One bench point per fluid.  (tau, delta) chosen well away from saturation /
// criticality so we exercise the hot loop without numerical edge cases.
struct BenchPoint
{
    const char* name;
    double tau;
    double delta;
};

const std::vector<BenchPoint> kBenchPoints = {
  // Light + heavy fluids spanning small/medium/large term tables
  {"Water", 1.4, 1.1},    {"CarbonDioxide", 1.2, 0.9}, {"n-Propane", 1.3, 0.6}, {"R134a", 1.2, 0.5},
  {"Methane", 1.5, 0.7},  {"Nitrogen", 1.4, 0.8},      {"Hydrogen", 1.3, 0.6},  {"MM", 1.5, 0.8},
  {"Ammonia", 1.3, 0.7},  {"R32", 1.3, 0.6},           {"R125", 1.3, 0.6},  // R125 uses tau_mi_in_u
  {"Methanol", 1.3, 0.7},                                                   // Methanol uses tau_mi_in_u
  {"R1234yf", 1.3, 0.6},  {"R1234ze(E)", 1.3, 0.6},    {"Ethane", 1.3, 0.7},    {"Propane", 1.3, 0.7},
  {"Argon", 1.4, 0.8},    {"Oxygen", 1.4, 0.8},        {"Toluene", 1.3, 0.6},   {"Benzene", 1.3, 0.7},
  {"R143a", 1.3, 0.6},    {"R152a", 1.3, 0.6},         {"R227ea", 1.3, 0.6},    {"D4", 1.4, 0.8},
};

constexpr std::size_t kReps = 200'000;
constexpr std::size_t kTrials = 7;

struct Stats
{
    double mean = 0.0, std = 0.0;
};

Stats stats_of(const std::vector<double>& xs) {
    Stats s;
    if (xs.empty()) return s;
    s.mean = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    double v = 0.0;
    for (double x : xs)
        v += (x - s.mean) * (x - s.mean);
    s.std = std::sqrt(v / xs.size());
    return s;
}

}  // namespace

// Three callables: scalar all(), vvexp-batched allFastVDSP.
// Add new fast paths here as they're implemented.
enum class WhichVariant
{
    Scalar,
    FastVDSP,
    FastNEON
};

void call_variant(CoolProp::ResidualHelmholtzGeneralizedExponential& gen, WhichVariant which, CoolPropDbl tau, CoolPropDbl delta,
                  CoolProp::HelmholtzDerivatives& d) {
    switch (which) {
        case WhichVariant::Scalar:
            gen.all(tau, delta, d);
            break;
        case WhichVariant::FastVDSP:
            gen.allFastVDSP(tau, delta, d);
            break;
        case WhichVariant::FastNEON:
            gen.allFastNEON(tau, delta, d);
            break;
    }
}

double measure_variant(CoolProp::ResidualHelmholtzGeneralizedExponential& gen, WhichVariant which, CoolPropDbl tau, CoolPropDbl delta) {
    std::vector<double> ns_per_call;
    ns_per_call.reserve(kTrials);
    for (std::size_t trial = 0; trial < kTrials; ++trial) {
        volatile double sink = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < kReps; ++i) {
            CoolProp::HelmholtzDerivatives d{};
            call_variant(gen, which, tau, delta, d);
            sink = sink + d.alphar;
        }
        const auto t1 = std::chrono::steady_clock::now();
        ns_per_call.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(kReps));
        (void)sink;
    }
    return stats_of(ns_per_call).mean;
}

TEST_CASE("Helmholtz inner-loop scalar baseline", "[helmholtz_inner_bench][.]") {
    std::FILE* csv = std::fopen("/tmp/helmholtz_inner_bench.csv", "w");
    if (csv) std::fprintf(csv, "fluid,N,scalar_ns_mean,scalar_ns_std\n");

    std::printf("\n=== Helmholtz inner-loop scalar baseline ===\n");
    std::printf("%-14s %4s   scalar (ns/call)\n", "fluid", "N");
    std::printf("%-14s %4s   ----------------\n", "-----", "-");

    double total_ns_mean = 0.0;
    std::size_t fluid_count = 0;

    for (const auto& bp : kBenchPoints) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bp.name));
        } catch (...) {
            std::printf("%-14s  SKIPPED (factory)\n", bp.name);
            continue;
        }
        auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        if (!be) continue;
        auto& comps = be->get_components();
        if (comps.empty()) continue;
        auto& gen = comps[0].EOS().alphar.GenExp;
        const auto tau = static_cast<CoolPropDbl>(bp.tau);
        const auto delta = static_cast<CoolPropDbl>(bp.delta);

        // Warm up cache + branch predictor
        {
            CoolProp::HelmholtzDerivatives d{};
            for (std::size_t i = 0; i < 1000; ++i)
                gen.all(tau, delta, d);
        }

        std::vector<double> ns_per_call;
        ns_per_call.reserve(kTrials);
        for (std::size_t trial = 0; trial < kTrials; ++trial) {
            volatile double sink = 0;
            const auto t0 = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < kReps; ++i) {
                CoolProp::HelmholtzDerivatives d{};
                gen.all(tau, delta, d);
                sink = sink + d.alphar;
            }
            const auto t1 = std::chrono::steady_clock::now();
            ns_per_call.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(kReps));
            (void)sink;
        }
        const auto s = stats_of(ns_per_call);
        std::printf("%-14s %4zu    %7.1f ± %5.1f\n", bp.name, gen.elements.size(), s.mean, s.std);
        if (csv) std::fprintf(csv, "%s,%zu,%.3f,%.3f\n", bp.name, gen.elements.size(), s.mean, s.std);
        total_ns_mean += s.mean;
        ++fluid_count;
    }
    if (csv) std::fclose(csv);
    if (fluid_count > 0) {
        std::printf("--- avg ns/call across %zu fluids: %.1f ---\n", fluid_count, total_ns_mean / fluid_count);
    }
}

TEST_CASE("allFastVDSP equivalence to scalar all()", "[helmholtz_inner_bench_vdsp][.]") {
    std::printf("\n=== allFastVDSP equivalence (worst |Δ|/max(|s|,1) across all fields) ===\n");
    std::printf("%-14s %4s   %-20s   %12s   %s\n", "fluid", "N", "regime", "worst_rel", "worst_field");
    struct TP
    {
        double tau, delta;
        const char* desc;
    };
    const std::vector<TP> regimes = {
      {1.4, 1.1, "near_critical"},
      {2.0, 2.5, "compressed_liquid"},
      {0.8, 0.3, "supercritical"},
      {3.0, 0.01, "dilute_gas"},
    };
    for (const auto& bp : kBenchPoints) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bp.name));
        } catch (...) {
            continue;
        }
        auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        if (!be) continue;
        auto& comps = be->get_components();
        if (comps.empty()) continue;
        auto& gen = comps[0].EOS().alphar.GenExp;
        for (const auto& r : regimes) {
            CoolProp::HelmholtzDerivatives s{}, v{};
            // all_scalar bypasses the env-toggled dispatcher so this is a
            // genuine scalar-vs-NEON/vvexp comparison, not a self-comparison.
            gen.all_scalar(static_cast<CoolPropDbl>(r.tau), static_cast<CoolPropDbl>(r.delta), s);
            gen.allFastVDSP(static_cast<CoolPropDbl>(r.tau), static_cast<CoolPropDbl>(r.delta), v);
            struct F
            {
                const char* n;
                double s, v;
            };
            const F fs[] = {
              {"alphar", s.alphar, v.alphar},
              {"dar_dd", s.dalphar_ddelta, v.dalphar_ddelta},
              {"dar_dt", s.dalphar_dtau, v.dalphar_dtau},
              {"d2ar_dd2", s.d2alphar_ddelta2, v.d2alphar_ddelta2},
              {"d2ar_dt2", s.d2alphar_dtau2, v.d2alphar_dtau2},
              {"d2ar_ddt", s.d2alphar_ddelta_dtau, v.d2alphar_ddelta_dtau},
              {"d3ar_dd3", s.d3alphar_ddelta3, v.d3alphar_ddelta3},
              {"d3ar_dt3", s.d3alphar_dtau3, v.d3alphar_dtau3},
              {"d4ar_dd4", s.d4alphar_ddelta4, v.d4alphar_ddelta4},
              {"d4ar_dt4", s.d4alphar_dtau4, v.d4alphar_dtau4},
            };
            double worst = 0;
            std::string worst_field;
            for (const auto& f : fs) {
                const double scale = std::max(std::abs(f.s), 1.0);
                const double rd = std::abs(f.s - f.v) / scale;
                if (rd > worst) {
                    worst = rd;
                    worst_field = f.n;
                }
            }
            std::printf("%-14s %4zu   %-20s   %12.3e   %s\n", bp.name, gen.elements.size(), r.desc, worst, worst_field.c_str());
            CHECK(worst < 1e-12);
        }
    }
}

// Equivalence test runs in CI default (no [.] hide).  Code review of the
// SIMD branch flagged that running the bespoke PXcdj_sweep alone (which is
// [.]-hidden) provided no CI coverage of the NEON path's correctness; this
// test now does, against the genuine scalar body (all_scalar) on every
// fluid in kBenchPoints across 4 regimes.
TEST_CASE("allFastNEON equivalence to scalar all()", "[helmholtz_inner_bench_neon][helmholtz_neon_equiv]") {
    std::printf("\n=== allFastNEON equivalence ===\n");
    struct TP
    {
        double tau, delta;
        const char* desc;
    };
    const std::vector<TP> regimes = {
      {1.4, 1.1, "near_critical"},
      {2.0, 2.5, "compressed_liquid"},
      {0.8, 0.3, "supercritical"},
      {3.0, 0.01, "dilute_gas"},
    };
    double worst_seen = 0.0;
    std::string worst_label;
    for (const auto& bp : kBenchPoints) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bp.name));
        } catch (...) {
            continue;
        }
        auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        if (!be) continue;
        auto& comps = be->get_components();
        if (comps.empty()) continue;
        auto& gen = comps[0].EOS().alphar.GenExp;
        for (const auto& r : regimes) {
            CoolProp::HelmholtzDerivatives s{}, v{};
            // all_scalar bypasses the env-toggled dispatcher so this is a
            // genuine scalar-vs-NEON comparison, not a self-comparison.
            gen.all_scalar(static_cast<CoolPropDbl>(r.tau), static_cast<CoolPropDbl>(r.delta), s);
            gen.allFastNEON(static_cast<CoolPropDbl>(r.tau), static_cast<CoolPropDbl>(r.delta), v);
            struct F
            {
                const char* n;
                double s, v;
            };
            const F fs[] = {
              {"alphar", s.alphar, v.alphar},
              {"dar_dd", s.dalphar_ddelta, v.dalphar_ddelta},
              {"dar_dt", s.dalphar_dtau, v.dalphar_dtau},
              {"d2ar_dd2", s.d2alphar_ddelta2, v.d2alphar_ddelta2},
              {"d2ar_dt2", s.d2alphar_dtau2, v.d2alphar_dtau2},
              {"d2ar_ddt", s.d2alphar_ddelta_dtau, v.d2alphar_ddelta_dtau},
              {"d3ar_dd3", s.d3alphar_ddelta3, v.d3alphar_ddelta3},
              {"d3ar_dt3", s.d3alphar_dtau3, v.d3alphar_dtau3},
              {"d4ar_dd4", s.d4alphar_ddelta4, v.d4alphar_ddelta4},
              {"d4ar_dt4", s.d4alphar_dtau4, v.d4alphar_dtau4},
            };
            double worst = 0;
            std::string worst_field;
            for (const auto& f : fs) {
                const double scale = std::max(std::abs(f.s), 1.0);
                const double rd = std::abs(f.s - f.v) / scale;
                if (rd > worst) {
                    worst = rd;
                    worst_field = f.n;
                }
            }
            if (worst > worst_seen) {
                worst_seen = worst;
                worst_label = std::string(bp.name) + " / " + r.desc + " / " + worst_field;
            }
            CHECK(worst < 1e-12);
        }
    }
    std::printf("worst across all fluid×regime×field: %.3e (%s)\n", worst_seen, worst_label.c_str());
}

TEST_CASE("Helmholtz inner-loop: scalar vs allFastNEON", "[helmholtz_inner_bench_neon][.]") {
    std::printf("\n=== scalar all() vs allFastNEON (vvexp + NEON 2-wide SIMD B-chain) ===\n");
    std::printf("%-14s %4s   scalar       neon         scalar/neon\n", "fluid", "N");
    std::printf("%-14s %4s   ------       ----         -----------\n", "-----", "-");
    double total_scalar = 0, total_neon = 0;
    std::size_t fluid_count = 0;
    for (const auto& bp : kBenchPoints) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bp.name));
        } catch (...) {
            continue;
        }
        auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        if (!be) continue;
        auto& comps = be->get_components();
        if (comps.empty()) continue;
        auto& gen = comps[0].EOS().alphar.GenExp;
        const auto tau = static_cast<CoolPropDbl>(bp.tau);
        const auto delta = static_cast<CoolPropDbl>(bp.delta);
        CoolProp::HelmholtzDerivatives d{};
        for (std::size_t i = 0; i < 1000; ++i)
            gen.all(tau, delta, d);
        for (std::size_t i = 0; i < 1000; ++i)
            gen.allFastNEON(tau, delta, d);
        const double scalar_ns = measure_variant(gen, WhichVariant::Scalar, tau, delta);
        const double neon_ns = measure_variant(gen, WhichVariant::FastNEON, tau, delta);
        const double ratio = scalar_ns / neon_ns;
        std::printf("%-14s %4zu   %7.1f      %7.1f      %5.3f%s\n", bp.name, gen.elements.size(), scalar_ns, neon_ns, ratio,
                    (ratio > 1.0 ? "  ✓" : "  ✗"));
        total_scalar += scalar_ns;
        total_neon += neon_ns;
        ++fluid_count;
    }
    if (fluid_count > 0) {
        std::printf("--- aggregate: scalar %.1f ns, neon %.1f ns, ratio %.3f ---\n", total_scalar, total_neon, total_scalar / total_neon);
    }
}

TEST_CASE("Helmholtz inner-loop: scalar vs allFastVDSP", "[helmholtz_inner_bench_vdsp][.]") {
    std::FILE* csv = std::fopen("/tmp/helmholtz_inner_bench_vdsp.csv", "w");
    if (csv) std::fprintf(csv, "fluid,N,scalar_ns,vdsp_ns,ratio\n");

    std::printf("\n=== scalar all() vs allFastVDSP (Accelerate vvexp batched) ===\n");
    std::printf("%-14s %4s   scalar       vdsp         scalar/vdsp\n", "fluid", "N");
    std::printf("%-14s %4s   ------       ----         -----------\n", "-----", "-");

    double total_scalar = 0, total_vdsp = 0;
    std::size_t fluid_count = 0;

    for (const auto& bp : kBenchPoints) {
        std::shared_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory("HEOS", bp.name));
        } catch (...) {
            continue;
        }
        auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        if (!be) continue;
        auto& comps = be->get_components();
        if (comps.empty()) continue;
        auto& gen = comps[0].EOS().alphar.GenExp;
        const auto tau = static_cast<CoolPropDbl>(bp.tau);
        const auto delta = static_cast<CoolPropDbl>(bp.delta);

        // Warm caches
        CoolProp::HelmholtzDerivatives d{};
        for (std::size_t i = 0; i < 1000; ++i)
            gen.all(tau, delta, d);
        for (std::size_t i = 0; i < 1000; ++i)
            gen.allFastVDSP(tau, delta, d);

        const double scalar_ns = measure_variant(gen, WhichVariant::Scalar, tau, delta);
        const double vdsp_ns = measure_variant(gen, WhichVariant::FastVDSP, tau, delta);
        const double ratio = scalar_ns / vdsp_ns;
        std::printf("%-14s %4zu   %7.1f      %7.1f      %5.3f%s\n", bp.name, gen.elements.size(), scalar_ns, vdsp_ns, ratio,
                    (ratio > 1.0 ? "  ✓" : "  ✗"));
        if (csv) std::fprintf(csv, "%s,%zu,%.3f,%.3f,%.4f\n", bp.name, gen.elements.size(), scalar_ns, vdsp_ns, ratio);
        total_scalar += scalar_ns;
        total_vdsp += vdsp_ns;
        ++fluid_count;
    }
    if (csv) std::fclose(csv);
    if (fluid_count > 0) {
        std::printf("--- aggregate (sum of means): scalar %.1f ns, vdsp %.1f ns, ratio %.3f ---\n", total_scalar, total_vdsp,
                    total_scalar / total_vdsp);
    }
}

#endif  // ENABLE_CATCH
