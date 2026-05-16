// Bulk-build SVDSBTL tables for a fluid set and write them to
// ~/.CoolProp/SVDTables/<fluid>.<input_pair>.svd.bin.z.
//
// Build:
//   cmake -B build -DCOOLPROP_BUILD_SVD_TABLES=ON
//   cmake --build build --target build_svd_tables -j
//
// Run:
//   ./build/build_svd_tables                 # builds the 7 Phase 2a fluids
//   ./build/build_svd_tables Water Methane   # build a subset
//
// Per fluid, builds PH and PT surfaces using the presets at production
// resolution (NT=200, NR=800, rank=20).  Reports per-file time, on-disk
// size, and a quick read-back sanity check.
//
// This is a developer tool — exempt from a few clang-tidy style nags
// (printf, deterministic-seed RNG, empty catch on optional HEOS
// failures).
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AbstractState.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SVDSurfaceSerializer.h"

namespace cp_sbtl = CoolProp::sbtl;

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

double file_size_kb(const std::string& path) {
    std::error_code ec;
    const auto sz = std::filesystem::file_size(path, ec);
    if (ec) {
        return 0.0;
    }
    return static_cast<double>(sz) / 1024.0;
}

void build_one(const std::string& fluid, std::size_t NT, std::size_t NR, std::int32_t rank) {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", fluid));

    for (const auto pair : {::CoolProp::HmassP_INPUTS, ::CoolProp::PT_INPUTS}) {
        const std::string label = (pair == ::CoolProp::HmassP_INPUTS) ? "PH" : "PT";
        const std::string out_path = cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid, pair);

        std::printf("  %s %s  building... ", fluid.c_str(), label.c_str());
        std::fflush(stdout);

        const auto t0 = std::chrono::steady_clock::now();
        cp_sbtl::SurfaceSpec spec = (pair == ::CoolProp::HmassP_INPUTS) ? cp_sbtl::presets::ph_subcritical(*heos, NT, NR, rank)
                                                                        : cp_sbtl::presets::pt_subcritical(*heos, NT, NR, rank);
        cp_sbtl::SVDSurface surface = cp_sbtl::build_surface(*heos, std::move(spec));
        const auto t1 = std::chrono::steady_clock::now();

        cp_sbtl::SVDSurfaceSerializer::save_to_file(surface, out_path);
        const auto t2 = std::chrono::steady_clock::now();

        // Quick load-back sanity check.
        const auto reloaded = cp_sbtl::SVDSurfaceSerializer::load_from_file(out_path);
        if (reloaded.region_count() != surface.region_count()) {
            throw std::runtime_error("build_one: load-back region_count mismatch for " + fluid + " " + label);
        }

        const double build_s = std::chrono::duration<double>(t1 - t0).count();
        const double save_s = std::chrono::duration<double>(t2 - t1).count();
        std::printf("build=%.1fs  save=%.2fs  size=%.1f kB  → %s\n", build_s, save_s, file_size_kb(out_path), out_path.c_str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    // Clamp env overrides to valid ranges before kicking off the build
    // loop -- (NT - 1) and (NR - 1) appear in grid denominators downstream,
    // and rank < 1 produces no usable decomposition.
    const std::size_t NT = std::max<std::size_t>(2, env_size_t("SVD_NT", 200));
    const std::size_t NR = std::max<std::size_t>(2, env_size_t("SVD_NR", 800));
    const std::int32_t rank = std::max<std::int32_t>(1, env_int("SVD_RANK", 20));

    std::vector<std::string> fluids;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            fluids.emplace_back(argv[i]);
        }
    } else {
        // Default: the 7 Phase 2a fluids.
        fluids = {"Water", "CarbonDioxide", "R134a", "Ammonia", "Methane", "Propane", "D6"};
    }

    std::printf("SVDSBTL bulk-build  (NT=%zu  NR=%zu  rank=%d)\n", NT, NR, rank);
    std::printf("Cache dir: %s\n\n", cp_sbtl::SVDSurfaceSerializer::default_cache_dir().c_str());

    int ok = 0;
    int failed = 0;
    const auto t_global = std::chrono::steady_clock::now();
    for (const auto& fluid : fluids) {
        try {
            build_one(fluid, NT, NR, rank);
            ++ok;
        } catch (const std::exception& e) {
            std::printf("  %s  ERROR: %s\n", fluid.c_str(), e.what());
            ++failed;
        }
    }
    const double total_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_global).count();
    std::printf("\nDone in %.1fs: %d ok, %d failed.\n", total_s, ok, failed);
    return failed == 0 ? 0 : 1;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)
