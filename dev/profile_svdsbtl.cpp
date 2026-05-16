// Microbenchmark that decomposes the SVDSBTL per-call cost into
// stages, so we can see where the 178 ns/call from
// dev/bench_svdsbtl_ph.cpp is actually going.
//
// Picks one representative warm state for Water (1 MPa, 350 K), runs
// each stage in a tight loop with volatile aggregation, and prints
// the breakdown.
//
// Build:
//   cmake -B build -DCOOLPROP_BUILD_SVDSBTL_PROFILE=ON
//   cmake --build build --target profile_svdsbtl -j
//
// Run:
//   ./build/profile_svdsbtl                    # Water at (1 MPa, 350 K)
//   ./build/profile_svdsbtl 5e5 400            # custom p (Pa), T (K)
//
// Env:
//   PROFILE_REPEATS (default 5,000,000)
//
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

#include "AbstractState.h"
#include "CoolProp/region/RegionAtlas.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SVDSurfaceSerializer.h"
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

// Time a callable that returns a double; returns mean ns/call.
template <typename F>
double timeit_double(F&& f, std::size_t repeats) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        sink += f();
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

// Time a callable that returns an int; returns mean ns/call.
template <typename F>
double timeit_int(F&& f, std::size_t repeats) {
    volatile int sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        sink += f();
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

void report(const char* label, double ns_per_call, double baseline_ns) {
    if (baseline_ns > 0.0) {
        std::printf("  %-58s  %7.1f ns   (%5.1f%% of baseline)\n", label, ns_per_call, 100.0 * ns_per_call / baseline_ns);
    } else {
        std::printf("  %-58s  %7.1f ns\n", label, ns_per_call);
    }
}

}  // namespace

int main(int argc, char** argv) {
    // Clamp to >=1 so the timing math doesn't divide by zero on
    // PROFILE_REPEATS=0.
    const std::size_t repeats = std::max<std::size_t>(1, env_size_t("PROFILE_REPEATS", 5'000'000));

    double p = 1.0e6;
    double T = 350.0;
    if (argc >= 3) {
        // Catch std::stod throws (invalid_argument / out_of_range) on
        // bad numeric input rather than letting the process abort with
        // a confusing terminate().  Same posture as the env-parsing
        // helpers above.
        try {
            p = std::stod(argv[1]);
            T = std::stod(argv[2]);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "profile_svdsbtl: bad numeric argument: %s (expected: ./profile_svdsbtl <p_Pa> <T_K>)\n", e.what());
            return 2;
        }
    }

    // Setup: HEOS for h_truth, SVDSBTL backend for end-to-end timing,
    // raw SVDSurface for per-stage timing.
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    heos->update(::CoolProp::PT_INPUTS, p, T);
    const double h = heos->hmass();
    const double rho_truth = heos->rhomass();

    auto svd_backend = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));

    // Direct surface (skip backend wrapper) -- load the same cache
    // file SVDSBTLBackend just hit.
    const std::string cache_path = ::CoolProp::sbtl::SVDSurfaceSerializer::default_cache_path("Water", ::CoolProp::HmassP_INPUTS);
    auto surface = ::CoolProp::sbtl::SVDSurfaceSerializer::load_from_file(cache_path);

    // Sanity check that the warm state lands in-domain.
    auto resolved = surface.resolve(p, h);
    if (resolved.region_idx < 0) {
        std::printf("ERROR: state (p=%.3g Pa, T=%.3g K, h=%.3g J/kg) is out of SVDSBTL domain\n", p, T, h);
        return 1;
    }

    std::printf("SVDSBTL timing breakdown -- Water at (p=%.3g Pa, T=%.3g K, h=%.3g J/kg)\n", p, T, h);
    std::printf("  rho_truth=%.4f kg/m^3  region_idx=%d  repeats=%zu\n\n", rho_truth, resolved.region_idx, repeats);

    // ===== Stage 1: full backend path =====
    const double ns_backend_full = timeit_double(
      [&] {
          svd_backend->update(::CoolProp::HmassP_INPUTS, h, p);
          return svd_backend->rhomass();
      },
      repeats);

    // Stage 2: surface.eval — skips the backend wrapper but does
    // atlas dispatch + transform + SVD eval.
    const double ns_surface_eval = timeit_double([&] { return surface.eval(::CoolProp::iDmass, p, h); }, repeats);

    // Stage 3: surface.eval_with_region — same as eval but caller
    // supplied region + transformed coords.  Pure SVD eval cost.
    const double ns_eval_with_region =
      timeit_double([&] { return surface.eval_with_region(::CoolProp::iDmass, resolved.region_idx, resolved.svd_x, resolved.svd_y); }, repeats);

    // Stage 4: atlas dispatch alone (no SVD eval, no transform).
    const double ns_atlas_find = timeit_int([&] { return surface.atlas().find_region(p, h); }, repeats);

    // Stage 5: just resolve (atlas + to_normalized) -- no SVD eval.
    const double ns_resolve = timeit_int(
      [&] {
          auto r = surface.resolve(p, h);
          return r.region_idx;
      },
      repeats);

    // ===== Stage 6: backend wrapper-only =====
    // update() only -- no rhomass() call afterwards.  This isolates
    // the "store inputs + resolve" work that SVDSBTLBackend::update()
    // does on top of surface.resolve().
    const double ns_backend_update = timeit_double(
      [&] {
          svd_backend->update(::CoolProp::HmassP_INPUTS, h, p);
          return 0.0;
      },
      repeats);

    // Stage 7: just calc_rhomass() after a single update() (cached path).
    svd_backend->update(::CoolProp::HmassP_INPUTS, h, p);  // prime once
    const double ns_backend_rhomass_cached = timeit_double([&] { return svd_backend->rhomass(); }, repeats);

    // Header
    std::printf("Stage breakdowns (mean ns / call, %zu iters each)\n", repeats);
    std::printf("  %-58s  %7s\n", "                              ", "ns/call");
    std::printf("  %-58s  %7s\n", "------------------------------", "-------");

    report("[A] Backend update + rhomass (the headline number)", ns_backend_full, 0.0);
    report("[B] surface.eval(iDmass, p, h)  -- no backend wrapper", ns_surface_eval, ns_backend_full);
    report("[C] surface.eval_with_region(...)  -- pure SVD math", ns_eval_with_region, ns_backend_full);
    report("[D] surface.atlas().find_region(p, h)  -- AABB + curve", ns_atlas_find, ns_backend_full);
    report("[E] surface.resolve(p, h)  -- atlas + to_normalized", ns_resolve, ns_backend_full);
    report("[F] backend.update() only (no rhomass)", ns_backend_update, ns_backend_full);
    report("[G] backend.rhomass() only (already updated)", ns_backend_rhomass_cached, ns_backend_full);

    std::printf("\nDerived breakdown of the [A] baseline:\n");
    const double backend_overhead = ns_backend_full - ns_surface_eval;
    const double atlas_dispatch = ns_atlas_find;
    const double transform_only = ns_resolve - ns_atlas_find;
    const double svd_eval = ns_eval_with_region;
    const double leftover = ns_surface_eval - ns_resolve - ns_eval_with_region;

    auto show = [&](const char* name, double ns) { std::printf("  %-58s  %7.1f ns   (%5.1f%%)\n", name, ns, 100.0 * ns / ns_backend_full); };
    show("backend wrapper (update + rhomass dispatch + cache)", backend_overhead);
    show("atlas.find_region (AABB scan + curve_contains)", atlas_dispatch);
    show("region.to_normalized (axis + boundary curve eval)", transform_only);
    show("SVDEvaluator inner product (rank-20 x 2D hermite)", svd_eval);
    show("(unaccounted -- inlining / branch / measurement noise)", leftover);

    return 0;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)
