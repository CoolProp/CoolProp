// State-point benchmark for the SVDSBTL backend in (h, p) coordinates.
//
// Walks a (T, p) grid via HEOS to get reference values, then queries
// SVDSBTL at (h_truth, p) and records: rel-err of rho / T / s / u,
// and the average per-call wall time for SVDSBTL update+rhomass and
// HEOS update+rhomass at that point.
//
// Output: a CSV at /tmp/bench_svdsbtl_ph_<fluid>.csv with one row per
// single-phase grid cell.
//
// Build:
//   cmake -B build -DCOOLPROP_BUILD_SVDSBTL_BENCH=ON
//   cmake --build build --target bench_svdsbtl_ph -j
//
// Run:
//   ./build/bench_svdsbtl_ph                  # defaults to Water
//   ./build/bench_svdsbtl_ph Water Propane    # subset
//
// Env:
//   BENCH_NT (default 80) — temperature grid points
//   BENCH_NP (default 60) — pressure grid points
//   BENCH_REPEATS (default 200) — calls per point for timing
//
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractState.h"
#include "Configuration.h"
#include "DataStructures.h"

namespace {

std::size_t env_size_t(const char* name, std::size_t fallback) {
    const char* v = std::getenv(name);
    if (v == nullptr || v[0] == 0) {
        return fallback;
    }
    try {
        return static_cast<std::size_t>(std::stoul(v));
    } catch (...) {
        return fallback;
    }
}

struct Row
{
    double h;
    double p;
    double T;
    double rho_truth;
    double rho_pred;
    double T_pred;
    double s_truth;
    double s_pred;
    double u_truth;
    double u_pred;
    double rel_err_rho;
    double rel_err_T;
    double rel_err_s;
    double rel_err_u;
    double ns_per_call_svd;
    double ns_per_call_heos;
    double ns_per_call_refprop;  // NaN when REFPROP unavailable
    double rho_refprop;          // NaN when REFPROP unavailable
    double rel_err_rho_refprop;  // SVDSBTL vs REFPROP (NaN if no REFPROP)
};

bool has_refprop_at_(const std::string& path) {
    return !path.empty();
}

// Run `repeats` (update + rhomass) cycles and return the average
// wall-clock per call in nanoseconds.  Uses volatile aggregation so
// the compiler can't optimise the work away.
double time_update_rhomass(::CoolProp::AbstractState& AS, ::CoolProp::input_pairs pair, double v1, double v2, std::size_t repeats) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        AS.update(pair, v1, v2);
        sink += AS.rhomass();
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

void bench_one(const std::string& fluid, std::size_t NT, std::size_t NP, std::size_t repeats, const std::string& refprop_path) {
    auto svd = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("SVDSBTL", fluid));
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", fluid));

    std::shared_ptr<::CoolProp::AbstractState> refprop;
    if (has_refprop_at_(refprop_path)) {
        ::CoolProp::set_config_string(ALTERNATIVE_REFPROP_PATH, refprop_path);
        try {
            refprop.reset(::CoolProp::AbstractState::factory("REFPROP", fluid));
        } catch (const std::exception& e) {
            std::printf("  REFPROP unavailable: %s\n", e.what());
            refprop.reset();
        }
    }

    const double T_min = std::max(heos->Ttriple(), heos->Tmin()) + 1.0;
    const double T_max = heos->Tmax() - 1.0;
    const double p_min = heos->p_triple() * 1.5;
    const double p_max = heos->p_critical() * 0.98;

    std::vector<Row> rows;
    rows.reserve(NT * NP);

    for (std::size_t i = 0; i < NT; ++i) {
        const double T = T_min + (T_max - T_min) * static_cast<double>(i) / static_cast<double>(NT - 1);
        for (std::size_t j = 0; j < NP; ++j) {
            const double log_p = std::log(p_min) + (std::log(p_max) - std::log(p_min)) * static_cast<double>(j) / static_cast<double>(NP - 1);
            const double p = std::exp(log_p);
            try {
                heos->update(::CoolProp::PT_INPUTS, p, T);
                if (heos->Q() > 0.0 && heos->Q() < 1.0) {
                    continue;  // two-phase, not covered by SVDSBTL subcritical regions
                }
                Row r{};
                r.h = heos->hmass();
                r.p = p;
                r.T = T;
                r.rho_truth = heos->rhomass();
                r.s_truth = heos->smass();
                r.u_truth = heos->umass();

                svd->update(::CoolProp::HmassP_INPUTS, r.h, r.p);
                r.rho_pred = svd->rhomass();
                r.T_pred = svd->T();
                r.s_pred = svd->smass();
                r.u_pred = svd->umass();

                if (std::isnan(r.rho_pred)) {
                    continue;  // outside any region (e.g. dome edge slip)
                }

                r.rel_err_rho = std::abs(r.rho_pred - r.rho_truth) / r.rho_truth;
                r.rel_err_T = std::abs(r.T_pred - r.T) / r.T;
                r.rel_err_s = (r.s_truth != 0.0) ? std::abs(r.s_pred - r.s_truth) / std::abs(r.s_truth) : 0.0;
                r.rel_err_u = (r.u_truth != 0.0) ? std::abs(r.u_pred - r.u_truth) / std::abs(r.u_truth) : 0.0;

                r.ns_per_call_svd = time_update_rhomass(*svd, ::CoolProp::HmassP_INPUTS, r.h, r.p, repeats);
                r.ns_per_call_heos = time_update_rhomass(*heos, ::CoolProp::PT_INPUTS, p, T, repeats);

                if (refprop) {
                    try {
                        // Use HmassP for a fair shape comparison with SVDSBTL.
                        refprop->update(::CoolProp::HmassP_INPUTS, r.h, p);
                        r.rho_refprop = refprop->rhomass();
                        r.rel_err_rho_refprop = std::abs(r.rho_pred - r.rho_refprop) / r.rho_refprop;
                        r.ns_per_call_refprop = time_update_rhomass(*refprop, ::CoolProp::HmassP_INPUTS, r.h, r.p, repeats);
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        r.rho_refprop = std::nan("");
                        r.rel_err_rho_refprop = std::nan("");
                        r.ns_per_call_refprop = std::nan("");
                    }
                } else {
                    r.rho_refprop = std::nan("");
                    r.rel_err_rho_refprop = std::nan("");
                    r.ns_per_call_refprop = std::nan("");
                }

                rows.push_back(r);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
            }
        }
    }

    const std::string out_path = std::string("/tmp/bench_svdsbtl_ph_") + fluid + ".csv";
    std::ofstream ofs(out_path);
    if (!ofs) {
        throw std::runtime_error("bench_one: failed to open " + out_path);
    }
    ofs << "h,p,T,rho_truth,rho_pred,T_pred,s_truth,s_pred,u_truth,u_pred,"
        << "rel_err_rho,rel_err_T,rel_err_s,rel_err_u,"
        << "ns_per_call_svd,ns_per_call_heos,ns_per_call_refprop,"
        << "rho_refprop,rel_err_rho_refprop\n";
    ofs.precision(17);
    for (const auto& r : rows) {
        ofs << r.h << "," << r.p << "," << r.T << "," << r.rho_truth << "," << r.rho_pred << "," << r.T_pred << "," << r.s_truth << "," << r.s_pred
            << "," << r.u_truth << "," << r.u_pred << "," << r.rel_err_rho << "," << r.rel_err_T << "," << r.rel_err_s << "," << r.rel_err_u << ","
            << r.ns_per_call_svd << "," << r.ns_per_call_heos << "," << r.ns_per_call_refprop << "," << r.rho_refprop << "," << r.rel_err_rho_refprop
            << "\n";
    }

    // Console summary.
    double max_rho = 0;
    double max_T = 0;
    double max_s = 0;
    double max_u = 0;
    double sum_ns_svd = 0;
    double sum_ns_heos = 0;
    double sum_ns_refprop = 0;
    std::size_t n_refprop = 0;
    for (const auto& r : rows) {
        max_rho = std::max(max_rho, r.rel_err_rho);
        max_T = std::max(max_T, r.rel_err_T);
        max_s = std::max(max_s, r.rel_err_s);
        max_u = std::max(max_u, r.rel_err_u);
        sum_ns_svd += r.ns_per_call_svd;
        sum_ns_heos += r.ns_per_call_heos;
        if (!std::isnan(r.ns_per_call_refprop)) {
            sum_ns_refprop += r.ns_per_call_refprop;
            ++n_refprop;
        }
    }
    const auto n = static_cast<double>(rows.size());
    if (n_refprop > 0) {
        const double mean_ns_refprop = sum_ns_refprop / static_cast<double>(n_refprop);
        std::printf("%-12s  N=%5zu  max_rel rho=%.2e T=%.2e s=%.2e u=%.2e  mean ns/call: svd=%.0f heos=%.0f refprop=%.0f"
                    "  speedup: vs heos=%.1fx vs refprop=%.1fx  -> %s\n",
                    fluid.c_str(), rows.size(), max_rho, max_T, max_s, max_u, sum_ns_svd / n, sum_ns_heos / n, mean_ns_refprop,
                    (sum_ns_heos / n) / (sum_ns_svd / n), mean_ns_refprop / (sum_ns_svd / n), out_path.c_str());
    } else {
        std::printf("%-12s  N=%5zu  max_rel: rho=%.2e T=%.2e s=%.2e u=%.2e  mean ns/call: svd=%.0f heos=%.0f  speedup=%.1fx  -> %s\n", fluid.c_str(),
                    rows.size(), max_rho, max_T, max_s, max_u, sum_ns_svd / n, sum_ns_heos / n, (sum_ns_heos / n) / (sum_ns_svd / n),
                    out_path.c_str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t NT = env_size_t("BENCH_NT", 80);
    const std::size_t NP = env_size_t("BENCH_NP", 60);
    const std::size_t repeats = env_size_t("BENCH_REPEATS", 200);
    const char* refprop_env = std::getenv("BENCH_REFPROP_PATH");
    const std::string refprop_path = (refprop_env != nullptr) ? std::string(refprop_env) : std::string();

    std::vector<std::string> fluids;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            fluids.emplace_back(argv[i]);
        }
    } else {
        fluids = {"Water"};
    }

    std::printf("SVDSBTL PH benchmark (NT=%zu  NP=%zu  repeats=%zu  refprop=%s)\n\n", NT, NP, repeats,
                refprop_path.empty() ? "off" : refprop_path.c_str());
    int failed = 0;
    for (const auto& fluid : fluids) {
        try {
            bench_one(fluid, NT, NP, repeats, refprop_path);
        } catch (const std::exception& e) {
            std::printf("%-12s  ERROR: %s\n", fluid.c_str(), e.what());
            ++failed;
        }
    }
    return failed == 0 ? 0 : 1;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)
