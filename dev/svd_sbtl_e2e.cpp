// Phase 2a end-to-end validation: multi-fluid SVD ρ(h, p) accuracy
// using the new CoolProp::svd + CoolProp::region components, ported
// from dev/svd_sbtl_error_plot.py.
//
// Per fluid, this executable:
//   1. Constructs a HEOS AbstractState.
//   2. Determines the subcritical pressure range
//      [max(1.5 * p_triple, 1e-3 * p_crit),  0.999 * p_crit].
//   3. Builds four PiecewiseChebyshevCurve boundaries over log p:
//        h(T_min,p)       ← LIQUID lower bound
//        h_sat,L(p)       ← LIQUID upper bound (and VAPOR-side dome)
//        h_sat,V(p)       ← VAPOR lower bound
//        h(T_max - 0.5,p) ← VAPOR upper bound
//   4. Builds a 2-region atlas (LIQUID, VAPOR) keyed on (p, h).
//   5. For each region: samples log(ρ) on an (xnorm, xi_log_p) grid
//      of size NT × NR by calling HEOS::HmassP_INPUTS at each cell;
//      builds a rank-RANK SVD with OutputTransform::EXP.
//   6. Probes N_PROBES random single-phase (T, p) states, computes
//      h via HEOS::PT_INPUTS, dispatches via the atlas, evaluates
//      the appropriate region's SVD, and compares to HEOS rhomass().
//   7. Writes a CSV of (T, p, h, rho_truth, rho_pred, rel_err) for
//      the companion Python plotter.
//   8. Prints a one-line summary per fluid plus an aggregate.
//
// Build: enable -DCOOLPROP_SVD_E2E=ON.  Run: ./build/SVDSBTL_E2E.
//
// This is a developer tool, not production code; the file is exempt
// from a handful of clang-tidy style checks (printf and empty-catch
// idioms, deterministic-seed RNG for reproducible probes).
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cppcoreguidelines-pro-bounds-pointer-arithmetic,cert-err33-c,cert-msc32-c,cert-msc51-cpp,bugprone-empty-catch,hicpp-vararg)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractState.h"

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/PiecewiseChebyshevCurve.h"
#include "CoolProp/region/Region.h"
#include "CoolProp/region/RegionAtlas.h"
#include "CoolProp/svd/SVDBuilder.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "CoolProp/svd/SVDEvaluator.h"

namespace cp_svd = CoolProp::svd;
namespace cp_region = CoolProp::region;

namespace {

// Grid sizes start small so the smoke-iteration loop is fast.  Once the
// pipeline produces sane numbers we'll ramp NT/NR back to the PoC's
// (200, 800).  Setting via environment variables (SVD_NT, SVD_NR,
// SVD_RANK, SVD_PROBES) overrides the defaults so the production
// validation run doesn't require a rebuild.
constexpr std::size_t NT_DEFAULT = 100;
constexpr std::size_t NR_DEFAULT = 200;
constexpr std::int32_t RANK_DEFAULT = 20;
constexpr std::size_t N_PROBES_DEFAULT = 5000;
constexpr std::size_t SAT_PIECES = 6;   // sat-curve Chebyshev pieces
constexpr std::size_t SAT_DEGREE = 10;  // sat-curve Chebyshev degree

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

std::int32_t env_int(const char* name, std::int32_t fallback) {
    const char* v = std::getenv(name);
    if (v == nullptr || v[0] == 0) {
        return fallback;
    }
    try {
        return static_cast<std::int32_t>(std::stoi(v));
    } catch (...) {
        return fallback;
    }
}

std::size_t NT = NT_DEFAULT;
std::size_t NR = NR_DEFAULT;
std::int32_t RANK = RANK_DEFAULT;
std::size_t N_PROBES = N_PROBES_DEFAULT;

// Materialised per-region build artefacts.  The SVDDecomposition is
// heap-allocated so its address is stable even if RegionData is moved
// (the SVDEvaluator stores a borrowed pointer that would dangle if the
// decomposition lived inline in a vector that gets push_back'd).
struct RegionData
{
    std::unique_ptr<cp_svd::SVDDecomposition> decomp;
    std::unique_ptr<cp_svd::SVDEvaluator> eval;
};

struct FluidStats
{
    std::string name;
    std::size_t kept = 0;
    double max_rel = 0.0;
    double p99_rel = 0.0;
    double p50_rel = 0.0;
    double build_seconds = 0.0;
    double eval_seconds = 0.0;
};

// Compute the percentile of a sorted vector.  pct in [0, 1].
double percentile_sorted(const std::vector<double>& v, double pct) {
    if (v.empty()) {
        return std::nan("");
    }
    const auto n = static_cast<double>(v.size());
    const auto idx = static_cast<std::size_t>(std::min(n - 1.0, std::max(0.0, pct * (n - 1.0))));
    return v[idx];
}

// Build the four boundary curves for one fluid and return the assembled
// 2-region atlas alongside p_min / p_max for the SVD grid generation.
struct AtlasBundle
{
    cp_region::RegionAtlas atlas;
    double p_min;
    double p_max;
};

AtlasBundle build_atlas(CoolProp::AbstractState& heos) {
    const double p_crit = heos.p_critical();
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();
    // The PoC's triple-point pressure: a hair above the actual triple
    // point to avoid PQ flash failures at the boundary.
    heos.update(CoolProp::QT_INPUTS, 0.0, heos.Ttriple() * 1.001);
    const double p_triple = heos.p();
    const double p_min = std::max(p_triple * 1.5, p_crit * 1e-3);
    const double p_max = 0.999 * p_crit;

    using PCC = cp_region::PiecewiseChebyshevCurve;

    // LIQUID lower bound: h on the cold-isotherm floor.  For fluids
    // with steep melting curves (Methane, Propane), T_min_eos can lie
    // below T_melt(p) at high pressure; in that case walk T up in 0.5 K
    // increments until HEOS accepts the (T, p) state.  The resulting
    // floor is non-isothermal but still defines a valid region with
    // monotone `h_lo(p)`.
    auto h_lo_liq = PCC::build(p_min, p_max, SAT_PIECES, SAT_DEGREE, PCC::ParamScale::LOG, [&](double p) {
        for (int k = 0; k < 40; ++k) {
            const double T_try = T_min_eos + 0.5 * k;
            try {
                heos.update(CoolProp::PT_INPUTS, p, T_try);
                return heos.hmass();
            } catch (...) {
                // Below melting line at this p; bump T up and retry.
            }
        }
        throw std::runtime_error("LIQUID floor unreachable within 20 K of T_min_eos");
    });
    // LIQUID upper bound / VAPOR side of dome: h_sat,L(p).
    auto h_sat_L = PCC::build(p_min, p_max, SAT_PIECES, SAT_DEGREE, PCC::ParamScale::LOG, [&](double p) {
        heos.update(CoolProp::PQ_INPUTS, p, 0.0);
        return heos.hmass();
    });
    // VAPOR lower bound: h_sat,V(p).
    auto h_sat_V = PCC::build(p_min, p_max, SAT_PIECES, SAT_DEGREE, PCC::ParamScale::LOG, [&](double p) {
        heos.update(CoolProp::PQ_INPUTS, p, 1.0);
        return heos.hmass();
    });
    // VAPOR upper bound: h on the hot-isotherm ceiling.  Stay safely
    // inside the HEOS validity envelope (Tmax - 0.5 K).
    auto h_hi_vap = PCC::build(p_min, p_max, SAT_PIECES, SAT_DEGREE, PCC::ParamScale::LOG, [&](double p) {
        heos.update(CoolProp::PT_INPUTS, p, T_max_eos - 0.5);
        return heos.hmass();
    });

    // Atlas: region 0 = LIQUID, region 1 = VAPOR.  Both share the same
    // primary AxisTransform (LOG over [p_min, p_max]).
    cp_region::RegionAtlas atlas;
    atlas.add(cp_region::Region(cp_region::AxisTransform::make(cp_region::AxisScale::LOG, p_min, p_max), std::move(h_lo_liq), std::move(h_sat_L)));
    auto axis_v = cp_region::AxisTransform::make(cp_region::AxisScale::LOG, p_min, p_max);
    atlas.add(cp_region::Region(axis_v, std::move(h_sat_V), std::move(h_hi_vap)));
    return AtlasBundle{std::move(atlas), p_min, p_max};
}

// Sample log(ρ) on the (xnorm, xi_log_p) ∈ [0,1]² grid for one region
// and build the rank-truncated SVD with EXP output transform.
cp_svd::SVDDecomposition build_region_svd(CoolProp::AbstractState& heos, const cp_region::Region& region) {
    std::vector<double> xn_grid(NT);
    std::vector<double> xi_grid(NR);
    // Avoid the exact endpoints — Chebyshev curves extrapolate slightly
    // and HmassP flash failures cluster at the corners.
    for (std::size_t i = 0; i < NT; ++i) {
        xn_grid[i] = 0.001 + (0.999 - 0.001) * static_cast<double>(i) / static_cast<double>(NT - 1);
    }
    for (std::size_t j = 0; j < NR; ++j) {
        xi_grid[j] = static_cast<double>(j) / static_cast<double>(NR - 1);
    }

    std::vector<double> M(NT * NR, std::nan(""));
    for (std::size_t j = 0; j < NR; ++j) {
        // Map xi to p, then evaluate both bounds at this p.
        const double p = region.primary().inverse(xi_grid[j]);
        for (std::size_t i = 0; i < NT; ++i) {
            // Map (xnorm, p) to (xnorm, h) via the region's boundaries.
            const auto [_, h] = region.from_normalized(xi_grid[j], xn_grid[i]);
            try {
                heos.update(CoolProp::HmassP_INPUTS, h, p);
                if (heos.Q() > 0.0 && heos.Q() < 1.0) {
                    continue;
                }
                const double rho = heos.rhomass();
                if (rho > 0.0 && std::isfinite(rho)) {
                    M[i * NR + j] = std::log(rho);
                }
            } catch (...) {
                // HEOS occasionally fails very close to the dome / corners;
                // leave the cell as NaN and let the row-fill below patch it.
            }
        }
    }
    // Row-wise linear fill along log p for any holes (typically near the
    // dome corner).  Same trick as the PoC.
    for (std::size_t i = 0; i < NT; ++i) {
        std::size_t first_good = NR, last_good = NR;
        for (std::size_t j = 0; j < NR; ++j) {
            if (std::isfinite(M[i * NR + j])) {
                if (first_good == NR) {
                    first_good = j;
                }
                last_good = j;
            }
        }
        if (first_good == NR) {
            // No data in this row at all — fall back to the median of
            // the whole matrix.  Extremely rare for the PoC's range.
            continue;
        }
        // Fill leading and trailing NaNs with the nearest-good value.
        for (std::size_t j = 0; j < first_good; ++j) {
            M[i * NR + j] = M[i * NR + first_good];
        }
        for (std::size_t j = last_good + 1; j < NR; ++j) {
            M[i * NR + j] = M[i * NR + last_good];
        }
        // Linear-interp interior NaN runs.
        for (std::size_t j = first_good + 1; j < last_good; ++j) {
            if (std::isfinite(M[i * NR + j])) {
                continue;
            }
            std::size_t k_lo = j;
            while (k_lo > 0 && !std::isfinite(M[i * NR + k_lo - 1])) {
                --k_lo;
            }
            --k_lo;  // k_lo now points at the last finite cell before the run
            std::size_t k_hi = j;
            while (k_hi < NR && !std::isfinite(M[i * NR + k_hi])) {
                ++k_hi;
            }
            const double t = (xi_grid[j] - xi_grid[k_lo]) / (xi_grid[k_hi] - xi_grid[k_lo]);
            M[i * NR + j] = (1.0 - t) * M[i * NR + k_lo] + t * M[i * NR + k_hi];
        }
    }
    // Final guard: if any NaNs remain (shouldn't), replace with the
    // matrix median to keep the SVD numerics sane.
    std::vector<double> finite_vals;
    finite_vals.reserve(M.size());
    for (double v : M) {
        if (std::isfinite(v)) {
            finite_vals.push_back(v);
        }
    }
    if (finite_vals.size() < M.size()) {
        std::nth_element(finite_vals.begin(), finite_vals.begin() + finite_vals.size() / 2, finite_vals.end());
        const double median = finite_vals[finite_vals.size() / 2];
        for (double& v : M) {
            if (!std::isfinite(v)) {
                v = median;
            }
        }
    }

    cp_svd::SVDBuildOptions opts;
    opts.rank = std::min<std::int32_t>(RANK, static_cast<std::int32_t>(std::min(NT, NR)));
    opts.out_transform = cp_svd::OutputTransform::EXP;
    opts.slope_source = cp_svd::SlopeSource::NATURAL_CUBIC_SPLINE;
    return cp_svd::build_svd(xn_grid, xi_grid, M, opts);
}

FluidStats run_fluid(const std::string& fluid, const std::string& csv_dir) {
    FluidStats stats;
    stats.name = fluid;

    auto heos_ptr = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    CoolProp::AbstractState& heos = *heos_ptr;

    std::printf("== %s ==  ", fluid.c_str());
    std::fflush(stdout);

    auto t_build_0 = std::chrono::steady_clock::now();
    auto bundle = build_atlas(heos);
    std::vector<RegionData> regions;
    regions.reserve(bundle.atlas.size());
    for (std::size_t r = 0; r < bundle.atlas.size(); ++r) {
        regions.emplace_back();
        regions.back().decomp = std::make_unique<cp_svd::SVDDecomposition>(build_region_svd(heos, bundle.atlas.region(r)));
        regions.back().eval = std::make_unique<cp_svd::SVDEvaluator>(*regions.back().decomp);
    }
    auto t_build_1 = std::chrono::steady_clock::now();
    stats.build_seconds = std::chrono::duration<double>(t_build_1 - t_build_0).count();

    // Single-phase (T, p) probe set.
    std::mt19937 rng(42);
    const double T_min_q = std::max(heos.Ttriple(), heos.Tmin()) + 0.5;
    const double T_max_q = heos.Tmax() - 0.5;
    std::uniform_real_distribution<double> uT(T_min_q, T_max_q);
    std::uniform_real_distribution<double> u_log_p(std::log(bundle.p_min * 1.05), std::log(bundle.p_max * 0.99));

    std::ofstream csv(csv_dir + "/svd_sbtl_e2e_" + fluid + ".csv");
    csv << "T,p,h,rho_truth,rho_pred,rel_err\n";
    csv << std::scientific << std::setprecision(12);

    std::vector<double> errs;
    errs.reserve(N_PROBES);

    auto t_eval_0 = std::chrono::steady_clock::now();
    std::size_t attempt = 0;
    while (errs.size() < N_PROBES && attempt < N_PROBES * 5) {
        ++attempt;
        const double T = uT(rng);
        const double p = std::exp(u_log_p(rng));
        try {
            heos.update(CoolProp::PT_INPUTS, p, T);
            if (heos.Q() > 0.0 && heos.Q() < 1.0) {
                continue;
            }
            const double h = heos.hmass();
            const double rho_truth = heos.rhomass();
            const int region_idx = bundle.atlas.find_region(p, h);
            if (region_idx < 0) {
                continue;
            }
            // Region's to_normalized returns (xi_primary, eta_secondary)
            // = (xi_of_log_p, xnorm).  SVD's x-axis is xnorm and y-axis
            // is xi_of_log_p — so pass (eta, xi).
            const auto [xi, eta] = bundle.atlas.region(static_cast<std::size_t>(region_idx)).to_normalized(p, h);
            const double rho_pred = regions[static_cast<std::size_t>(region_idx)].eval->eval(eta, xi);
            const double rel = std::abs(rho_pred - rho_truth) / rho_truth;
            errs.push_back(rel);
            csv << T << ',' << p << ',' << h << ',' << rho_truth << ',' << rho_pred << ',' << rel << '\n';
        } catch (...) {
            // Treat any HEOS exception as a skipped probe.
        }
    }
    auto t_eval_1 = std::chrono::steady_clock::now();
    stats.eval_seconds = std::chrono::duration<double>(t_eval_1 - t_eval_0).count();

    stats.kept = errs.size();
    std::sort(errs.begin(), errs.end());
    stats.max_rel = errs.empty() ? std::nan("") : errs.back();
    stats.p99_rel = percentile_sorted(errs, 0.99);
    stats.p50_rel = percentile_sorted(errs, 0.50);

    std::printf("kept=%5zu  max=%.3e  p99=%.3e  p50=%.3e  build=%.1fs  eval=%.2fs\n", stats.kept, stats.max_rel, stats.p99_rel, stats.p50_rel,
                stats.build_seconds, stats.eval_seconds);
    return stats;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string csv_dir = (argc > 1) ? argv[1] : "/tmp";
    NT = env_size_t("SVD_NT", NT_DEFAULT);
    NR = env_size_t("SVD_NR", NR_DEFAULT);
    RANK = env_int("SVD_RANK", RANK_DEFAULT);
    N_PROBES = env_size_t("SVD_PROBES", N_PROBES_DEFAULT);

    const std::vector<std::string> fluids = {"Water", "R134a", "Ammonia", "Methane", "Propane"};

    std::printf("Phase 2a — SVD ρ(h,p) end-to-end validation\n");
    std::printf("  NT=%zu  NR=%zu  RANK=%d  N_PROBES=%zu\n", NT, NR, RANK, N_PROBES);
    std::printf("  Writing CSVs to %s/\n\n", csv_dir.c_str());

    std::vector<FluidStats> results;
    results.reserve(fluids.size());
    for (const auto& f : fluids) {
        try {
            results.push_back(run_fluid(f, csv_dir));
        } catch (const std::exception& ex) {
            std::printf("  ERROR (%s): %s\n", f.c_str(), ex.what());
        }
    }

    std::printf("\n%-12s %8s %12s %12s %12s\n", "fluid", "kept", "max", "p99", "p50");
    std::printf("%-12s %8s %12s %12s %12s\n", "-----", "----", "---", "---", "---");
    for (const auto& r : results) {
        std::printf("%-12s %8zu %12.3e %12.3e %12.3e\n", r.name.c_str(), r.kept, r.max_rel, r.p99_rel, r.p50_rel);
    }
    return 0;
}
// NOLINTEND(cppcoreguidelines-pro-type-vararg,cppcoreguidelines-pro-bounds-pointer-arithmetic,cert-err33-c,cert-msc32-c,cert-msc51-cpp,bugprone-empty-catch,hicpp-vararg)
