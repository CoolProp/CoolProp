#include "CoolProp/sbtl/SVDSurfaceFactory.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

#include "AbstractState.h"
#include "Configuration.h"
#include "CoolProp/region/Region.h"
#include "CoolProp/sbtl/SatBoundaryFactory.h"
#include "CoolProp/svd/SVDBuilder.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "DataStructures.h"

namespace CoolProp {
namespace sbtl {

namespace {

// Build the per-region (xn, xi) grids returned to both the sampler and
// the SVD builder.  Single source of truth for the secondary-axis
// spacing — the SVD's stored y_grid MUST match the η positions where
// the source was actually sampled, so sample_grid() and the
// build_surface() loop below both go through here.
//
// Secondary axis (xn) uses Chebyshev-cosine spacing: cells crowd toward
// both η=0 and η=1 boundaries.  In LIQUID/VAPOR regions the saturation
// curve is one of the η boundaries, and properties have sharp curvature
// in a thin η-strip right against it.  Uniform η spacing puts only ~1
// grid cell of resolution into that strip — measured worst-case error
// in SVDSBTL&HEOS R2 lives there (h_truth − h_sat,V ≈ 0.33%, η ≈ 0.005,
// only ~1 grid cell of support at NT=200).  Chebyshev squashes bulk
// cells toward the middle and crowds them at both ends, dropping
// per-region failure counts by 6-9× on the SVDSBTL&HEOS conformance
// fail-map (2026-05-18).
//
// Tiny pad off the literal 0 / 1 corners keeps boundary-curve evaluation
// safe — the cubic-spline / Chebyshev sat boundaries can extrapolate
// epsilon-outside their fit interval at the dome corners, and HmassP
// flash failures cluster there.
//
// Primary axis (xi) stays on uniform spacing — its scale is handled by
// the region's AxisTransform (linear or log) at from_normalized() time.
std::pair<std::vector<double>, std::vector<double>> make_grid_axes(std::size_t NT, std::size_t NR) {
    constexpr double kPadEnd = 0.001;
    std::vector<double> xn_grid(NT);
    std::vector<double> xi_grid(NR);
    for (std::size_t i = 0; i < NT; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(NT - 1);
        const double eta_cheb = 0.5 * (1.0 - std::cos(M_PI * t));
        xn_grid[i] = kPadEnd + (1.0 - 2.0 * kPadEnd) * eta_cheb;
    }
    for (std::size_t j = 0; j < NR; ++j) {
        xi_grid[j] = static_cast<double>(j) / static_cast<double>(NR - 1);
    }
    return {std::move(xn_grid), std::move(xi_grid)};
}

// Walk the (xnorm, log_p_norm) ∈ [0,1] × [0,1] grid for one Region,
// sample HEOS per cell via the spec's update_state / read_property
// callbacks, and fill the per-property matrix dictionary.  Cells
// where HEOS throws (e.g. two-phase) are left NaN and the row-fill
// patches them.  Returns NT × NR row-major matrices indexed by
// property.
//
// Matches dev/svd_sbtl_e2e.cpp's build_region_svd loop layout so
// numerics agree with Phase 2a's validated output.
std::vector<std::vector<double>> sample_grid(::CoolProp::AbstractState& heos, const SurfaceSpec& spec, const region::Region& region) {
    const std::size_t NT = spec.NT;
    const std::size_t NR = spec.NR;
    const std::size_t n_props = spec.properties.size();

    const auto [xn_grid, xi_grid] = make_grid_axes(NT, NR);

    // M[prop][i*NR + j] for each output property.
    std::vector<std::vector<double>> M(n_props);
    for (auto& m : M) {
        m.assign(NT * NR, std::nan(""));
    }

    // Per-cell sampling lambda — captures the grid axes + region by
    // reference (read-only) and the output matrix M by reference
    // (per-cell writes are to disjoint indices so concurrent
    // invocations are safe).  `src` is whichever AbstractState the
    // caller hands in — the original `heos` for the serial path, a
    // per-thread factory-built clone for the parallel path.
    auto sample_cell = [&](std::size_t i, std::size_t j, ::CoolProp::AbstractState& src) {
        // Avoid structured-binding capture (a C++20 feature; the
        // project compiles at C++17).  Use the pair fields directly.
        const auto ab = region.from_normalized(xi_grid[j], xn_grid[i]);
        const double a = ab.first;
        const double b = ab.second;
        try {
            spec.update_state(src, a, b);
            if (src.Q() > 0.0 && src.Q() < 1.0) {
                return;  // two-phase — leave NaN
            }
            for (std::size_t p = 0; p < n_props; ++p) {
                const double v = spec.read_property(src, spec.properties[p].key);
                if (!std::isfinite(v)) {
                    continue;
                }
                // Apply the property's transform (typically LOG
                // for density via OutputTransform::EXP) before
                // storing into M.
                const double stored = (spec.properties[p].transform == svd::OutputTransform::EXP) ? (v > 0.0 ? std::log(v) : std::nan("")) : v;
                M[p][i * NR + j] = stored;
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // HEOS failures (above melting line, OOB, ...) → NaN.
            // Row-fill below patches them.
        }
    };

    // Decide serial vs parallel.  Parallel needs all of:
    //   * config opt-in (SVDSBTL_SAMPLING_THREADS != 1)
    //   * non-REFPROP source (REFPROP's FORTRAN runtime is process-
    //     global; concurrent SETUPdll'd states share state and corrupt
    //     each other under parallel updates)
    //   * a non-empty spec.source_backend + spec.fluid_name (so worker
    //     threads can factory-build their own AbstractState)
    //   * at least 2 effective worker threads after resolving "auto"
    auto resolve_thread_count = []() -> std::size_t {
        const int cfg = ::CoolProp::get_config_int(SVDSBTL_SAMPLING_THREADS);
        if (cfg == 1) return 1;  // serial fast-path
        if (cfg <= 0) {
            // 0 (or negative) → auto = hardware_concurrency
            const unsigned hw = std::thread::hardware_concurrency();
            return hw <= 1 ? 1 : static_cast<std::size_t>(hw);
        }
        return static_cast<std::size_t>(cfg);
    };
    const std::size_t n_threads_cfg = resolve_thread_count();
    const bool can_parallelise = n_threads_cfg > 1 && spec.source_backend != "REFPROP" && !spec.source_backend.empty() && !spec.fluid_name.empty();

    if (!can_parallelise) {
        for (std::size_t j = 0; j < NR; ++j) {
            for (std::size_t i = 0; i < NT; ++i) {
                sample_cell(i, j, heos);
            }
        }
    } else {
        // Cap workers at NT (no point launching more threads than rows
        // — each thread gets one full row at a time).  Each worker
        // factory-builds its own AbstractState once and reuses across
        // its row slice; per-thread construction cost (~ms) amortises
        // across hundreds of cells per row.
        const std::size_t n_threads = std::min(n_threads_cfg, NT);
        auto worker = [&](std::size_t i_lo, std::size_t i_hi) {
            std::unique_ptr<::CoolProp::AbstractState> local_src;
            try {
                local_src.reset(::CoolProp::AbstractState::factory(spec.source_backend, spec.fluid_name));
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // Factory failure leaves local_src null → bail; the
                // serial fallback below would have hit the same error
                // on the shared heos handle.
                return;
            }
            for (std::size_t i = i_lo; i < i_hi; ++i) {
                for (std::size_t j = 0; j < NR; ++j) {
                    sample_cell(i, j, *local_src);
                }
            }
        };
        std::vector<std::thread> workers;
        workers.reserve(n_threads);
        const std::size_t chunk = (NT + n_threads - 1) / n_threads;
        for (std::size_t t = 0; t < n_threads; ++t) {
            const std::size_t i_lo = t * chunk;
            const std::size_t i_hi = std::min(i_lo + chunk, NT);
            if (i_lo >= i_hi) break;
            workers.emplace_back(worker, i_lo, i_hi);
        }
        for (auto& w : workers) {
            w.join();
        }
    }

    // Row-wise linear fill in log_p for each property (mirrors the
    // PoC's NaN-patch logic).  Used to keep the SVD numerics sane on
    // cells where HEOS gave up.
    for (auto& m : M) {
        for (std::size_t i = 0; i < NT; ++i) {
            std::size_t first_good = NR, last_good = NR;
            for (std::size_t j = 0; j < NR; ++j) {
                if (std::isfinite(m[i * NR + j])) {
                    if (first_good == NR) {
                        first_good = j;
                    }
                    last_good = j;
                }
            }
            if (first_good == NR) {
                continue;
            }
            for (std::size_t j = 0; j < first_good; ++j) {
                m[i * NR + j] = m[i * NR + first_good];
            }
            for (std::size_t j = last_good + 1; j < NR; ++j) {
                m[i * NR + j] = m[i * NR + last_good];
            }
            for (std::size_t j = first_good + 1; j < last_good; ++j) {
                if (std::isfinite(m[i * NR + j])) {
                    continue;
                }
                std::size_t k_lo = j;
                while (k_lo > 0 && !std::isfinite(m[i * NR + k_lo - 1])) {
                    --k_lo;
                }
                --k_lo;
                std::size_t k_hi = j;
                while (k_hi < NR && !std::isfinite(m[i * NR + k_hi])) {
                    ++k_hi;
                }
                const double t = (xi_grid[j] - xi_grid[k_lo]) / (xi_grid[k_hi] - xi_grid[k_lo]);
                m[i * NR + j] = (1.0 - t) * m[i * NR + k_lo] + t * m[i * NR + k_hi];
            }
        }
        // Final guard: any remaining NaNs (fully-NaN rows in pathological cases)
        // get filled with the median.  Falling through with NaN would crash the
        // SVD builder via the finite-value precondition.
        std::vector<double> finite_vals;
        finite_vals.reserve(m.size());
        for (double v : m) {
            if (std::isfinite(v)) {
                finite_vals.push_back(v);
            }
        }
        if (finite_vals.size() < m.size()) {
            if (finite_vals.empty()) {
                // Every sample for this property in this region was
                // non-finite — typically means the (T, p) grid walked
                // outside the HEOS validity envelope for this region's
                // axis transform.  Failing fast beats producing a NaN-
                // filled "table" and silently corrupting the SVD math.
                throw std::runtime_error("build_surface: all-NaN property matrix in region (HEOS sampling produced no finite values)");
            }
            std::nth_element(finite_vals.begin(), finite_vals.begin() + static_cast<std::ptrdiff_t>(finite_vals.size() / 2), finite_vals.end());
            const double median = finite_vals[finite_vals.size() / 2];
            for (double& v : m) {
                if (!std::isfinite(v)) {
                    v = median;
                }
            }
        }
    }

    return M;
}

}  // namespace

SVDSurface build_surface(::CoolProp::AbstractState& heos, SurfaceSpec spec, const BuildOptions& opts) {
    if (spec.regions.empty()) {
        throw std::invalid_argument("build_surface: SurfaceSpec.regions is empty");
    }
    if (spec.properties.empty()) {
        throw std::invalid_argument("build_surface: SurfaceSpec.properties is empty");
    }
    if (!spec.update_state || !spec.read_property) {
        throw std::invalid_argument("build_surface: SurfaceSpec is missing update_state / read_property callbacks");
    }
    // Grid axes use (NT - 1) / (NR - 1) as denominators; values < 2 would
    // emit NaN-shaped grids and crash downstream SVD math.  Rank < 1 has
    // no usable decomposition.
    if (spec.NT < 2 || spec.NR < 2) {
        throw std::invalid_argument("build_surface: SurfaceSpec.NT and .NR must each be >= 2");
    }
    if (spec.rank < 1) {
        throw std::invalid_argument("build_surface: SurfaceSpec.rank must be >= 1");
    }

    // Build the SVDSurface, transferring regions in order, then SVDs.
    std::vector<::CoolProp::parameters> prop_keys;
    prop_keys.reserve(spec.properties.size());
    for (const auto& p : spec.properties) {
        prop_keys.push_back(p.key);
    }
    SVDSurface surface(spec.fluid_name, spec.input_pair, prop_keys);

    const auto [xn_grid, xi_grid] = make_grid_axes(spec.NT, spec.NR);

    for (std::size_t r = 0; r < spec.regions.size(); ++r) {
        // Move the RegionSpec's boundary curves into a real Region,
        // then immediately add it to the surface.
        region::Region region(spec.regions[r].primary, std::move(spec.regions[r].b_lo), std::move(spec.regions[r].b_hi));
        // Sample HEOS for this region BEFORE moving it into the
        // surface — we need the (a, b) → from_normalized while we
        // still hold an immediate reference.
        const auto M_per_prop = sample_grid(heos, spec, region);
        const std::size_t region_idx = surface.add_region(std::move(region));

        // Build one SVD per property.
        for (std::size_t p = 0; p < spec.properties.size(); ++p) {
            svd::SVDBuildOptions svdopts;
            svdopts.rank = spec.rank;
            svdopts.out_transform = spec.properties[p].transform;
            svdopts.slope_source = opts.slope_source;
            svd::SVDDecomposition decomp = svd::build_svd(xn_grid, xi_grid, M_per_prop[p], svdopts);
            surface.add_region_property_svd(region_idx, spec.properties[p].key, std::move(decomp));
            if (opts.verbose) {
                // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
                std::printf("  region %zu property %d: SVD built (rank=%d)\n", region_idx, static_cast<int>(spec.properties[p].key), spec.rank);
            }
        }
    }

    surface.seal();
    return surface;
}

}  // namespace sbtl
}  // namespace CoolProp
