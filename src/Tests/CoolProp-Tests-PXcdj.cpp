// Characterization tests for P + {H, S, U} flash inputs.
//
// The legacy P+X solver is HSU_P_flash_singlephase_Brent in
// src/Backends/Helmholtz/FlashRoutines.cpp; PH/PS/PU all dispatch through
// FlashRoutines::HSU_P_flash.  This file is the acceptance harness for the
// rewrite tracked in bd CoolProp-cdj — it characterises both the named
// regressions we expect to keep passing (CI default) and a broad single-phase
// (T, p) sweep that maps the full surface for the production code path
// (opt-in heavy test, writes a CSV fixture).
//
// Two test cases:
//   [PXcdj]              named regressions, must pass on current master
//   [PXcdj_sweep][.]     opt-in dense grid across CoolProp's pure-fluid list
//                        — writes a CSV with timing + round-trip results, used
//                        to produce dev/fixtures/pxcdj_master_baseline.csv.gz
//
// Input pair argument orders (verified against include/DataStructures.h):
//   HmolarP_INPUTS : (value1=h, value2=p)
//   PSmolar_INPUTS : (value1=p, value2=s)
//   PUmolar_INPUTS : (value1=p, value2=u)
//
// Built into CatchTestRunner when COOLPROP_CATCH_MODULE=ON.  Run via
//   ./build_catch_rel/CatchTestRunner "[PXcdj]"        # CI default
//   ./build_catch_rel/CatchTestRunner "[PXcdj_sweep]"  # opt-in heavy

#include "AbstractState.h"
#include "CoolProp.h"
#include "CPstrings.h"
#include "DataStructures.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <chrono>
#    include <cmath>
#    include <cstdio>
#    include <cstdlib>
#    include <fstream>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

// Read an env var as a size_t, returning def if absent or malformed.
std::size_t env_or(const char* name, std::size_t def) {
    const char* v = std::getenv(name);
    if (v == nullptr) return def;
    try {
        const long n = std::stol(v);
        if (n > 1) return static_cast<std::size_t>(n);
    } catch (...) {
        return def;
    }
    return def;
}

std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = (n == 1) ? a : a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return v;
}

std::vector<double> logspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = (n == 1) ? a : a * std::pow(b / a, static_cast<double>(i) / static_cast<double>(n - 1));
    }
    return v;
}

const char* pair_label(CoolProp::input_pairs pr) {
    switch (pr) {
        case CoolProp::HmolarP_INPUTS:
            return "PH";
        case CoolProp::PSmolar_INPUTS:
            return "PS";
        case CoolProp::PUmolar_INPUTS:
            return "PU";
        default:
            return "??";
    }
}

double caloric_molar(CoolProp::AbstractState& AS, CoolProp::input_pairs pr) {
    switch (pr) {
        case CoolProp::HmolarP_INPUTS:
            return AS.hmolar();
        case CoolProp::PSmolar_INPUTS:
            return AS.smolar();
        case CoolProp::PUmolar_INPUTS:
            return AS.umolar();
        default:
            return _HUGE;
    }
}

// Re-flash via the chosen P+X pair and check (T, rho) recover within 1e-5
// relative.  Uses CHECK (not FAIL_CHECK on the recover assertions) so callers
// see WHICH coordinate missed.  Wraps the flash in try/catch and FAIL_CHECKs
// on throw — every named-regression point must complete on current master.
void check_PX_roundtrip(CoolProp::AbstractState& wrk, CoolProp::input_pairs pr, double p, double value, double T_expect, double rho_expect) {
    INFO("pair=" << pair_label(pr) << " p=" << p << " value=" << value << " T_expect=" << T_expect << " rho_expect=" << rho_expect);
    try {
        if (pr == CoolProp::HmolarP_INPUTS) {
            wrk.update(pr, value, p);  // (h, p)
        } else {
            wrk.update(pr, p, value);  // (p, s) or (p, u)
        }
    } catch (const std::exception& e) {
        FAIL_CHECK("update threw: " << e.what());
        return;
    } catch (...) {
        FAIL_CHECK("update threw non-std exception");
        return;
    }
    const double T_out = wrk.T(), rho_out = wrk.rhomolar();
    CHECK(std::isfinite(T_out));
    CHECK(std::isfinite(rho_out));
    CHECK(T_out == Catch::Approx(T_expect).epsilon(1e-5));
    CHECK(rho_out == Catch::Approx(rho_expect).epsilon(1e-5));
}

// One named-regression record.
struct NamedCase
{
    const char* fluid;
    double T_K;
    double p_Pa;
    const char* note;
};

// Round-trip a single (fluid, T, p) point across all three P+X pairs.  The
// reference state is set via PT_INPUTS; rho, h, s, u are read off; then
// PH/PS/PU each get a fresh state and must recover (T, rho).
void run_named_case(const NamedCase& c) {
    CAPTURE(c.fluid, c.note, c.T_K, c.p_Pa);
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
    try {
        ref->update(CoolProp::PT_INPUTS, c.p_Pa, c.T_K);
    } catch (const std::exception& e) {
        FAIL_CHECK("PT reference update threw for " << c.fluid << ": " << e.what());
        return;
    }
    const double T_ref = ref->T(), rho_ref = ref->rhomolar();
    const double h = ref->hmolar(), s = ref->smolar(), u = ref->umolar();
    if (!std::isfinite(rho_ref) || !std::isfinite(h) || !std::isfinite(s) || !std::isfinite(u)) {
        FAIL_CHECK("non-finite reference state for " << c.fluid);
        return;
    }
    {
        INFO("PH");
        check_PX_roundtrip(*wrk, CoolProp::HmolarP_INPUTS, c.p_Pa, h, T_ref, rho_ref);
    }
    {
        INFO("PS");
        check_PX_roundtrip(*wrk, CoolProp::PSmolar_INPUTS, c.p_Pa, s, T_ref, rho_ref);
    }
    {
        INFO("PU");
        check_PX_roundtrip(*wrk, CoolProp::PUmolar_INPUTS, c.p_Pa, u, T_ref, rho_ref);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// CI default: focused PH/PS/PU round-trips that should PASS on current master.
// Every point is established via PT_INPUTS first, then re-flashed through each
// of HmolarP/PSmolar/PUmolar — the recovered (T, rho) must agree to 1e-5
// relative.  Covers:
//   * gas / liquid / supercritical points across a representative fluid set
//     (Water, CO2, R134a, Propane, Nitrogen, MM)
//   * three prototype-era known-hard points that current master DOES handle:
//       - Nitrogen PS at the supercritical-cold wrong-basin trap
//       - R134a PH at the spinodal step
//       - Water near the density anomaly
//
// NOTE: The two Nitrogen PS supercritical-cold failures on master at
// (T=86.35K, p=4.754MPa) and (T=90.20K, p=7.768MPa) are deliberately omitted
// here — they are fixed in PR #3020.  They show up in the [PXcdj_sweep] CSV.
// ---------------------------------------------------------------------------
TEST_CASE("PXcdj named regressions", "[PXcdj]") {
    static const NamedCase kCases[] = {
      // Water: gas, liquid, supercritical
      {"Water", 600.0, 1.0e5, "gas"},
      {"Water", 320.0, 5.0e6, "liquid"},
      {"Water", 700.0, 30.0e6, "supercrit"},
      // CO2: gas, supercritical, liquid
      {"CarbonDioxide", 350.0, 1.0e6, "gas"},
      {"CarbonDioxide", 320.0, 10.0e6, "supercrit"},
      {"CarbonDioxide", 250.0, 5.0e6, "liquid"},
      // R134a: liquid, gas
      {"R134a", 250.0, 3.0e5, "liquid"},
      {"R134a", 350.0, 2.0e5, "gas"},
      // Propane: gas, liquid
      {"n-Propane", 350.0, 1.0e6, "gas"},
      {"n-Propane", 250.0, 5.0e6, "liquid"},
      // Nitrogen: gas, supercritical (NOT the 86.35K / 90.20K failing points)
      {"Nitrogen", 300.0, 1.0e5, "gas"},
      {"Nitrogen", 200.0, 10.0e6, "supercrit"},
      // MM (hexamethyldisiloxane): gas, liquid
      {"MM", 450.0, 1.0e5, "gas"},
      {"MM", 350.0, 5.0e6, "liquid"},
      // Prototype-era known-hard points that master handles correctly.
      {"Nitrogen", 121.0, 3.7e6, "supercritical-cold wrong-basin trap"},
      {"R134a", 250.0, 3.0e5, "spinodal step"},
      {"Water", 277.0, 5.0e6, "density anomaly"},
    };
    for (const auto& c : kCases) {
        DYNAMIC_SECTION(c.fluid << " " << c.note << " T=" << c.T_K << "K p=" << c.p_Pa << "Pa") {
            run_named_case(c);
        }
    }
}

// ---------------------------------------------------------------------------
// Opt-in heavy sweep: iterate CoolProp's pure-fluid list, build a dense
// (log-p, linear-T) grid, establish a reference via PT_INPUTS, then re-flash
// each P+X pair through the production code and check (T, rho) recovery.
// Writes a single CSV with timing + per-point success — the moderate (30x30)
// run is committed as dev/fixtures/pxcdj_master_baseline.csv.gz to serve as
// a regression fingerprint.
//
// Env knobs:
//   PXCDJ_NT      grid size in T (default 40)
//   PXCDJ_NP      grid size in p (default 40)
//   PXCDJ_CSV     output CSV path (default "pxcdj_sweep.csv" in cwd)
//   PXCDJ_REPS    timing reps per point (default 3); best of reps is kept
// ---------------------------------------------------------------------------
TEST_CASE("PXcdj broad sweep", "[PXcdj_sweep][.]") {
    const std::size_t NT = env_or("PXCDJ_NT", 40);
    const std::size_t NP = env_or("PXCDJ_NP", 40);
    const std::size_t reps = env_or("PXCDJ_REPS", 3);
    const char* csv_env = std::getenv("PXCDJ_CSV");
    const std::string csv_path = (csv_env != nullptr) ? csv_env : "pxcdj_sweep.csv";

    std::ofstream out(csv_path);
    if (!out.is_open()) {
        FAIL("Failed to open PXcdj sweep CSV: " << csv_path);
    }
    out << "fluid,pair,T_K,p_Pa,rho_molm3,value,T_solved,rho_solved,success,microseconds\n";

    const std::string fluids_csv = CoolProp::get_global_param_string("FluidsList");
    const std::vector<std::string> fluids = strsplit(fluids_csv, ',');

    static const CoolProp::input_pairs kPairs[] = {CoolProp::HmolarP_INPUTS, CoolProp::PSmolar_INPUTS, CoolProp::PUmolar_INPUTS};

    std::size_t fluids_attempted = 0, fluids_completed = 0;
    std::size_t grand_total = 0, grand_ok = 0;

    for (const auto& fluid : fluids) {
        if (fluid.empty()) continue;
        ++fluids_attempted;

        std::shared_ptr<CoolProp::AbstractState> ref, wrk;
        try {
            ref.reset(CoolProp::AbstractState::factory("HEOS", fluid));
            wrk.reset(CoolProp::AbstractState::factory("HEOS", fluid));
        } catch (const std::exception& e) {
            std::printf("[PXcdj_sweep] %s: factory threw — skipped (%s)\n", fluid.c_str(), e.what());
            continue;
        } catch (...) {
            std::printf("[PXcdj_sweep] %s: factory threw non-std — skipped\n", fluid.c_str());
            continue;
        }

        double Tmin = 0, Tmax = 0, pmax = 0, ptriple = 0;
        try {
            Tmin = ref->Tmin();
            Tmax = ref->Tmax();
            pmax = ref->pmax();
        } catch (const std::exception& e) {
            std::printf("[PXcdj_sweep] %s: limits query threw — skipped (%s)\n", fluid.c_str(), e.what());
            continue;
        }
        try {
            ptriple = ref->p_triple();
        } catch (...) {
            ptriple = 0;
        }
        const double plo = std::max(std::max(ptriple, pmax * 1e-8), 1.0);
        if (!(Tmax > Tmin) || !(pmax > plo)) {
            std::printf("[PXcdj_sweep] %s: degenerate (T,p) bounds — skipped\n", fluid.c_str());
            continue;
        }

        const auto T_vec = linspace(Tmin + 0.1, Tmax, NT);
        const auto p_vec = logspace(plo, pmax, NP);

        // Per-pair counters: total points, successful round-trips.
        std::size_t pair_total[3] = {0, 0, 0};
        std::size_t pair_ok[3] = {0, 0, 0};
        std::size_t fluid_total = 0, fluid_ok = 0;

        for (double T : T_vec) {
            for (double p : p_vec) {
                // Establish reference state; skip out-of-domain or two-phase.
                try {
                    ref->update(CoolProp::PT_INPUTS, p, T);
                } catch (...) {
                    continue;
                }
                if (ref->phase() == CoolProp::iphase_twophase) continue;
                const double rho_ref = ref->rhomolar();
                if (!std::isfinite(rho_ref) || rho_ref <= 0) continue;

                for (std::size_t pi = 0; pi < 3; ++pi) {
                    const CoolProp::input_pairs pr = kPairs[pi];
                    double value = _HUGE;
                    try {
                        value = caloric_molar(*ref, pr);
                    } catch (...) {
                        continue;
                    }
                    if (!std::isfinite(value)) continue;

                    ++pair_total[pi];
                    ++fluid_total;
                    ++grand_total;

                    double best_us = 1e300;
                    bool threw = false;
                    double T_out = std::nan(""), rho_out = std::nan("");
                    for (std::size_t r = 0; r < reps; ++r) {
                        const auto t0 = std::chrono::steady_clock::now();
                        try {
                            if (pr == CoolProp::HmolarP_INPUTS) {
                                wrk->update(pr, value, p);  // (h, p)
                            } else {
                                wrk->update(pr, p, value);  // (p, s/u)
                            }
                        } catch (...) {
                            threw = true;
                        }
                        const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
                        if (us < best_us) best_us = us;
                        if (threw) break;  // no point timing repeats once it threw
                    }

                    bool point_ok = false;
                    if (!threw) {
                        try {
                            T_out = wrk->T();
                            rho_out = wrk->rhomolar();
                        } catch (...) {
                            T_out = std::nan("");
                            rho_out = std::nan("");
                        }
                        if (std::isfinite(T_out) && std::isfinite(rho_out) && std::abs(T_out - T) / T < 1e-5
                            && std::abs(rho_out - rho_ref) / rho_ref < 1e-5) {
                            point_ok = true;
                        }
                    }
                    if (point_ok) {
                        ++pair_ok[pi];
                        ++fluid_ok;
                        ++grand_ok;
                    }

                    // Always write a row, even on failure / throw.  T_solved /
                    // rho_solved are nan on throw; success is 0/1.
                    out << fluid << ',' << pair_label(pr) << ',' << T << ',' << p << ',' << rho_ref << ',' << value << ',' << T_out << ',' << rho_out
                        << ',' << (point_ok ? 1 : 0) << ',' << best_us << '\n';
                }
            }
        }

        ++fluids_completed;
        std::printf("[PXcdj_sweep] %s: %zu/%zu round-trips (PH:%zu/%zu PS:%zu/%zu PU:%zu/%zu)\n", fluid.c_str(), fluid_ok, fluid_total, pair_ok[0],
                    pair_total[0], pair_ok[1], pair_total[1], pair_ok[2], pair_total[2]);
    }

    out.close();
    std::printf("[PXcdj_sweep] wrote %zu rows to %s\n", grand_total, csv_path.c_str());
    std::printf("[PXcdj_sweep] fluids: %zu attempted, %zu completed; round-trips: %zu/%zu (%.2f%%)\n", fluids_attempted, fluids_completed, grand_ok,
                grand_total, (grand_total > 0) ? 100.0 * grand_ok / grand_total : 0.0);
    // The sweep is a characterisation map, not a pass/fail: don't fail the
    // test when some points fail (e.g. the Nitrogen/PS supercritical-cold
    // points fixed by PR #3020).  The CSV is the artifact of interest.
}

#endif  // ENABLE_CATCH
