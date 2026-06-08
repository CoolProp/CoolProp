// Emit per-region CSVs of conformance residuals against IAPWS IF97
// truth, for both SVDSBTL&HEOS and SVDSBTL&IF97 source backends.
// Tagged [!benchmark] so it doesn't run on every CI invocation;
// invoke explicitly with one of:
//
//   ./CatchTestRunner '[SVDSBTL][fail_map][if97_src]'         > /tmp/failmap_if97.csv
//   ./CatchTestRunner '[SVDSBTL][fail_map][heos_src_vs_if97]' > /tmp/failmap_heos_vs_if97.csv
//   ./CatchTestRunner '[SVDSBTL][fail_map][heos_src_vs_heos]' > /tmp/failmap_heos_vs_heos.csv

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <cstdio>
#    include <memory>
#    include <random>
#    include <string>

#    include "CoolProp/AbstractState.h"

namespace {

struct Box
{
    double T_lo, T_hi, p_lo_MPa, p_hi_MPa;
};
const std::vector<std::pair<int, Box>> kRegions = {
  {1, {273.15, 623.15, 1.0e-3, 100.0}},
  {2, {273.15, 1073.15, 1.0e-3, 100.0}},
  {3, {623.15, 863.15, 16.5292, 100.0}},
  {5, {1073.15, 2273.15, 1.0e-3, 50.0}},
};

double p_B23_MPa(double T) {
    return 0.10192970039326e-2 * (T - 0.57254459862746e3) * (T - 0.57254459862746e3) + 0.1391883776670e2;
}

void refine_to_forward_h(::CoolProp::AbstractState& s, double p, double T_forward, double h_target) {
    try {
        s.update(::CoolProp::HmassP_INPUTS, h_target, p);
    } catch (...) {  // NOLINT(bugprone-empty-catch)
        s.update(::CoolProp::PT_INPUTS, p, T_forward);
        return;
    }
    const ::CoolProp::phases phase0 = s.phase();
    const bool single_phase =
      (phase0 == ::CoolProp::iphase_liquid || phase0 == ::CoolProp::iphase_gas || phase0 == ::CoolProp::iphase_supercritical_liquid
       || phase0 == ::CoolProp::iphase_supercritical_gas || phase0 == ::CoolProp::iphase_supercritical);
    bool pinned = false;
    if (single_phase) {
        try {
            s.specify_phase(phase0);
            pinned = true;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // HEOS doesn't always honor specify_phase mid-iteration; skip.
        }
    }
    for (int k = 0; k < 5; ++k) {
        const double h_now = s.hmass();
        const double dh = h_now - h_target;
        if (std::abs(dh) < 1e-10 * std::abs(h_target)) break;
        const double cp = s.cpmass();
        if (!std::isfinite(cp) || cp <= 0.0) break;
        const double T_new = s.T() - dh / cp;
        try {
            s.update(::CoolProp::PT_INPUTS, p, T_new);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            break;
        }
    }
    if (pinned) {
        try {
            s.unspecify_phase();
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    }
}

int classify(double T, double p_MPa, CoolProp::AbstractState& heos) {
    if (T < 273.15 || T > 2273.15) return -1;
    if (p_MPa <= 0.0 || p_MPa > 100.0) return -1;
    if (T > 1073.15) return (p_MPa <= 50.0) ? 5 : -1;
    if (T <= 623.15) {
        try {
            heos.update(::CoolProp::QT_INPUTS, 0.0, T);
            const double psat_MPa = heos.p() / 1e6;
            return (p_MPa > psat_MPa) ? 1 : 2;
        } catch (...) {
            return -1;
        }
    }
    if (T <= 863.15 && p_MPa > p_B23_MPa(T)) return 3;
    return 2;
}

void run_failmap(const std::string& source_backend, const std::string& truth_backend, int rank = 0, int NT = 0, int NR = 0) {
    // Optional grid overrides via factory-string options blob.  rank<=0 / NT<=0 / NR<=0
    // keeps the corresponding default (200 / 800 / 20 per SVDSBTLBackend's kDefaultGrid).
    std::string fluid_with_opts = "Water";
    std::vector<std::string> kv;
    if (NT > 0) kv.push_back("\"NT\":" + std::to_string(NT));
    if (NR > 0) kv.push_back("\"NR\":" + std::to_string(NR));
    if (rank > 0) kv.push_back("\"rank\":" + std::to_string(rank));
    if (!kv.empty()) {
        std::string opts = "\"grid\":{";
        for (size_t i = 0; i < kv.size(); ++i) {
            if (i) opts += ",";
            opts += kv[i];
        }
        opts += "}";
        fluid_with_opts = "Water?{" + opts + "}";
    }
    auto svd = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("SVDSBTL&" + source_backend, fluid_with_opts));
    auto truth = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory(truth_backend, "Water"));
    auto classifier = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));

    constexpr int N_per_region = 2000;
    std::mt19937_64 rng(0xC0DECAFE);

    std::printf("region,h_J_kg,p_Pa,T_truth_K,dT_mK,dv_pct,ds_J_kgK,dw_pct\n");
    for (const auto& [region, box] : kRegions) {
        std::uniform_real_distribution<double> uT(box.T_lo, box.T_hi);
        std::uniform_real_distribution<double> uLp(std::log(box.p_lo_MPa), std::log(box.p_hi_MPa));
        for (int i = 0; i < N_per_region; ++i) {
            const double T = uT(rng);
            const double p_MPa = std::exp(uLp(rng));
            if (classify(T, p_MPa, *classifier) != region) continue;
            const double p_Pa = p_MPa * 1e6;
            try {
                truth->update(::CoolProp::PT_INPUTS, p_Pa, T);
                if (truth->Q() > 0.0 && truth->Q() < 1.0) continue;
                const double h_truth = truth->hmass();
                refine_to_forward_h(*truth, p_Pa, T, h_truth);
                const double T_truth_refined = truth->T();
                const double v_truth = 1.0 / truth->rhomass();
                const double s_truth = truth->smass();
                const double w_truth = truth->speed_sound();
                svd->update(::CoolProp::HmassP_INPUTS, h_truth, p_Pa);
                const double T_svd = svd->T();
                const double v_svd = 1.0 / svd->rhomass();
                const double s_svd = svd->smass();
                const double w_svd = svd->speed_sound();
                const double dT_mK = (T_svd - T_truth_refined) * 1000.0;
                const double dv_pct = (v_svd - v_truth) / v_truth * 100.0;
                const double ds = s_svd - s_truth;
                const double dw_pct = (w_svd - w_truth) / w_truth * 100.0;
                std::printf("%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", region, h_truth, p_Pa, T_truth_refined, dT_mK, dv_pct, ds, dw_pct);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                continue;
            }
        }
    }
    std::fflush(stdout);
}

}  // namespace

TEST_CASE("SVDSBTL&IF97 vs IF97 conformance fail-map", "[SVDSBTL][fail_map][if97_src][!benchmark]") {
    run_failmap("IF97", "IF97");
}

TEST_CASE("SVDSBTL&HEOS vs IF97 conformance fail-map", "[SVDSBTL][fail_map][heos_src_vs_if97][!benchmark]") {
    run_failmap("HEOS", "IF97");
}

TEST_CASE("SVDSBTL&HEOS vs HEOS conformance fail-map", "[SVDSBTL][fail_map][heos_src_vs_heos][!benchmark]") {
    run_failmap("HEOS", "HEOS");
}

#endif  // ENABLE_CATCH
