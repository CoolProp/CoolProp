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
#include "CoolProp/px_preconditioner.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <atomic>
#    include <chrono>
#    include <cmath>
#    include <cstdio>
#    include <cstdlib>
#    include <fstream>
#    include <memory>
#    include <string>
#    include <thread>
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

// ---------------------------------------------------------------------------
// Phase 3a — preconditioner unit tests.
//
// The preconditioner is a standalone unit (include/CoolProp/px_preconditioner.h,
// src/Backends/Helmholtz/PXPreconditioner.cpp).  Phase 3b will wire it into
// HSU_P_flash_singlephase_Brent; this phase just unit-tests the standalone
// pieces:
//   1. ensure_melt_curve_table / MeltCurveTable::rho_at
//   2. regime_classified_rho_seed (supercritical / liquid / vapor / edges)
//   3. lazy-init idempotence + (light) thread safety
// ---------------------------------------------------------------------------

namespace {

// Cast helper: get the HelmholtzEOSMixtureBackend from an AbstractState built
// via the "HEOS" factory.  Returns nullptr if the cast fails (shouldn't on the
// pure fluids we test below, but be defensive).
CoolProp::HelmholtzEOSMixtureBackend* as_heos(CoolProp::AbstractState& AS) {
    return dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(&AS);
}

}  // namespace

TEST_CASE("PXcdj preconditioner: melt-curve table", "[PXcdj][PXcdj_precond]") {
    // Fluids known to ship a melting line in CoolProp (dev/fluids/*.json with a
    // "melting_line" block).  R134a / R32 / R1234yf etc. do NOT have one; the
    // no-melt branch is exercised in the dedicated section below.
    static const char* kMeltFluids[] = {"Water", "CarbonDioxide", "n-Propane", "Methane"};
    for (const char* fluid : kMeltFluids) {
        DYNAMIC_SECTION("melt table: " << fluid) {
            std::shared_ptr<CoolProp::AbstractState> AS;
            try {
                AS.reset(CoolProp::AbstractState::factory("HEOS", fluid));
            } catch (const std::exception& e) {
                FAIL("factory threw for " << fluid << ": " << e.what());
            }
            auto* H = as_heos(*AS);
            REQUIRE(H != nullptr);
            REQUIRE(H->has_melting_line());

            CoolProp::pxprecond::ensure_melt_curve_table(*H);
            const auto& tbl = CoolProp::pxprecond::get_melt_curve_table(*H);
            CAPTURE(tbl.p.size());
            REQUIRE(tbl.has_data);
            REQUIRE(tbl.p.size() >= 5);
            REQUIRE(tbl.p.size() == tbl.T_melt.size());
            REQUIRE(tbl.p.size() == tbl.rho_melt.size());

            // Monotone-increasing in p.
            for (std::size_t i = 1; i < tbl.p.size(); ++i) {
                REQUIRE(tbl.p[i] > tbl.p[i - 1]);
            }

            // Spot-check Water at p=1MPa: melt-side liquid density ~ 55000 mol/m^3
            // (standard water is ~ 55556 mol/m^3 at NTP; near the melting point it
            // is similar to within ~50%).
            if (std::string(fluid) == "Water") {
                const double rhoM = tbl.rho_at(1.0e6);
                CAPTURE(rhoM);
                REQUIRE(std::isfinite(rhoM));
                REQUIRE(rhoM > 0.5 * 55000.0);
                REQUIRE(rhoM < 1.5 * 55000.0);
                // Diagnostic printout — useful for the wake-up report.
                std::printf("[PXcdj_precond] Water melt table (%zu pts):\n", tbl.p.size());
                std::printf("  p_lo=%.4e  T_melt=%.4f K  rho_melt=%.4f mol/m^3\n", tbl.p.front(), tbl.T_melt.front(), tbl.rho_melt.front());
                std::printf("  p_hi=%.4e  T_melt=%.4f K  rho_melt=%.4f mol/m^3\n", tbl.p.back(), tbl.T_melt.back(), tbl.rho_melt.back());
                for (double pq : {1.0e6, 1.0e7, 1.0e8, 5.0e8}) {
                    std::printf("  rho_at(%.2e Pa) = %.4f mol/m^3\n", pq, tbl.rho_at(pq));
                }
            }

            // rho_at across the table range returns finite positive values.
            for (std::size_t i = 0; i < tbl.p.size(); ++i) {
                const double rho = tbl.rho_at(tbl.p[i]);
                REQUIRE(std::isfinite(rho));
                REQUIRE(rho > 0.0);
            }
            // Mid-points too.
            for (std::size_t i = 1; i < tbl.p.size(); ++i) {
                const double pm = 0.5 * (tbl.p[i - 1] + tbl.p[i]);
                const double rho = tbl.rho_at(pm);
                REQUIRE(std::isfinite(rho));
                REQUIRE(rho > 0.0);
            }

            // Clamp on either side — must not throw, must return finite positive.
            const double rho_lo = tbl.rho_at(tbl.p.front() * 0.001);
            const double rho_hi = tbl.rho_at(tbl.p.back() * 1000.0);
            CHECK(std::isfinite(rho_lo));
            CHECK(rho_lo > 0.0);
            CHECK(std::isfinite(rho_hi));
            CHECK(rho_hi > 0.0);
            // Clamped value at p_lo*1e-3 == rho_melt.front(); at p_hi*1e3 == rho_melt.back().
            CHECK(rho_lo == tbl.rho_melt.front());
            CHECK(rho_hi == tbl.rho_melt.back());
        }
    }

    SECTION("fluid without a melting line") {
        // Find the first fluid (from a small candidate list) that has no melting
        // line.  Document which we picked.  R134a / R1234yf / R32 / R125 do not
        // ship a melting line.  Helium, IsoButane, n-Butane DO have one (per
        // dev/fluids/*.json) — keep them out of this candidate list.
        static const char* kNoMeltCandidates[] = {"R134a", "R1234yf", "R32", "R125", "MM", "MDM"};
        const char* picked = nullptr;
        std::shared_ptr<CoolProp::AbstractState> AS;
        for (const char* f : kNoMeltCandidates) {
            try {
                AS.reset(CoolProp::AbstractState::factory("HEOS", f));
            } catch (...) {
                continue;
            }
            if (!AS->has_melting_line()) {
                picked = f;
                break;
            }
        }
        CAPTURE(picked);
        REQUIRE(picked != nullptr);
        auto* H = as_heos(*AS);
        REQUIRE(H != nullptr);
        CoolProp::pxprecond::ensure_melt_curve_table(*H);
        const auto& tbl = CoolProp::pxprecond::get_melt_curve_table(*H);
        REQUIRE_FALSE(tbl.has_data);
        REQUIRE(tbl.p.empty());
        // rho_at on a no-data table returns NaN.
        const double r = tbl.rho_at(1.0e6);
        CHECK(std::isnan(r));
    }
}

TEST_CASE("PXcdj preconditioner: regime-classified seed", "[PXcdj][PXcdj_precond]") {
    using CoolProp::pxprecond::regime_classified_rho_seed;

    struct Point
    {
        const char* fluid;
        double T;
        double p;
        const char* regime;  // "super", "liquid", "vapor", "edge"
        const char* note;
    };
    // Representative points.  Pick points clearly inside each regime
    // (don't sit near T_sat or near Tc — there the seed can legitimately be
    // far from rho_sat,L or rho_sat,V).  Tc's used to pick: Water 647.1K,
    // CO2 304.1K, n-Propane 369.8K, R134a 374.2K.
    //
    // Regime labels:
    //   "super-ideal"  T > 0.8*Tmax or p < 0.01*pc -> seed == p/(RT) exactly
    //   "super-blend"  supercritical blend; seed lies between rho_c and p/(RT)
    //   "liquid"       T < T_sat(p), p > p_sat(T); seed within factor 2 of rho_sat,L
    //   "vapor"        T < Tc, p < p_sat(T); seed within factor 2 of rho_sat,V or p/(RT)
    static const Point kPts[] = {
      // Supercritical-ideal: hit the cutoffs.
      // Water: Tmax=2000K so 0.8*Tmax=1600K; pc=22.064 MPa so 0.01*pc=220.64 kPa.
      // 200 kPa is below the low-p cutoff.
      {"Water", 700.0, 2.0e5, "super-ideal", "low-p super (< 0.01*pc) -> ideal-gas"},
      // CO2: pc=7.37 MPa, 0.01*pc=73.7 kPa; 50 kPa is below.
      {"CarbonDioxide", 400.0, 5.0e4, "super-ideal", "low-p super -> ideal-gas"},
      // Supercritical-blend: moderate p, moderate T above Tc.
      {"Water", 700.0, 30.0e6, "super-blend", "blended"},
      {"CarbonDioxide", 400.0, 1.0e6, "super-blend", "blended"},
      {"n-Propane", 450.0, 5.0e5, "super-blend", "blended"},
      // Subcritical liquid: T well below T_sat(p).  All Ts < Tc, p above p_sat(T).
      {"Water", 320.0, 5.0e6, "liquid", "moderate p"},
      {"CarbonDioxide", 250.0, 5.0e6, "liquid", "moderate p"},
      {"n-Propane", 280.0, 2.0e6, "liquid", "moderate p"},
      {"R134a", 280.0, 5.0e5, "liquid", "moderate p"},
      // Subcritical vapor (T < Tc and p < p_sat(T)).  Choose T well below Tc.
      {"Water", 500.0, 1.0e5, "vapor", "T<Tc(647), p<p_sat(500)"},
      {"CarbonDioxide", 280.0, 1.0e5, "vapor", "T<Tc(304), p<p_sat(280)"},
      {"n-Propane", 320.0, 1.0e5, "vapor", "T<Tc(370), p<p_sat(320)"},
      {"R134a", 350.0, 1.0e5, "vapor", "T<Tc(374), p<p_sat(350)"},
    };

    for (const auto& pt : kPts) {
        DYNAMIC_SECTION(pt.fluid << " " << pt.regime << " T=" << pt.T << "K p=" << pt.p << "Pa (" << pt.note << ")") {
            std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", pt.fluid));
            auto* H = as_heos(*AS);
            REQUIRE(H != nullptr);

            const double seed = regime_classified_rho_seed(*H, pt.T, pt.p);
            CAPTURE(seed);
            REQUIRE(std::isfinite(seed));
            REQUIRE(seed > 0.0);

            const double R = H->gas_constant();
            const double ig = pt.p / (R * pt.T);
            REQUIRE(std::isfinite(ig));
            REQUIRE(ig > 0.0);

            if (std::string(pt.regime) == "super-ideal") {
                // The ideal-gas branches (high T or low p) should give exactly p/(RT).
                CHECK(seed == Catch::Approx(ig).epsilon(1e-12));
            } else if (std::string(pt.regime) == "super-blend") {
                // Blended (1-alpha)*rho_c + alpha*p/(RT).  Seed lies between
                // ideal-gas and rho_c (whichever is larger/smaller); just
                // sanity-check it is bounded by min/max of those two.
                const double rho_c = H->rhomolar_critical();
                REQUIRE(std::isfinite(rho_c));
                REQUIRE(rho_c > 0.0);
                const double lo = std::min(ig, rho_c);
                const double hi = std::max(ig, rho_c);
                CAPTURE(rho_c, lo, hi);
                CHECK(seed >= lo - 1e-9 * hi);
                CHECK(seed <= hi + 1e-9 * hi);
            } else if (std::string(pt.regime) == "liquid") {
                // Compare to rho_sat,L(T) via the same superanc the seed used.
                auto sa = H->get_superanc();
                REQUIRE(sa);
                double rhoL = 0.0;
                REQUIRE_NOTHROW(rhoL = sa->eval_sat(pt.T, 'D', 0));
                REQUIRE(rhoL > 0.0);
                CAPTURE(rhoL);
                // Seed must be within a factor of 2 of rho_sat,L for moderate p.
                CHECK(seed > 0.5 * rhoL);
                CHECK(seed < 2.0 * rhoL);
            } else if (std::string(pt.regime) == "vapor") {
                auto sa = H->get_superanc();
                REQUIRE(sa);
                double rhoV = 0.0;
                REQUIRE_NOTHROW(rhoV = sa->eval_sat(pt.T, 'D', 1));
                REQUIRE(rhoV > 0.0);
                CAPTURE(rhoV);
                // Seed must be within a factor of 2 of rho_sat,V for moderate p.
                // The seed may also legitimately equal p/(RT) for low p — accept that too.
                const bool near_rhoV = (seed > 0.5 * rhoV && seed < 2.0 * rhoV);
                const bool near_ig = (seed > 0.5 * ig && seed < 2.0 * ig);
                CAPTURE(near_rhoV, near_ig, rhoV, ig);
                CHECK((near_rhoV || near_ig));
            }
        }
    }

    SECTION("edge cases") {
        // For each of a few fluids, hit edges and check seed is finite positive.
        for (const char* f : {"Water", "CarbonDioxide", "n-Propane", "R134a"}) {
            DYNAMIC_SECTION(f) {
                std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
                auto* H = as_heos(*AS);
                REQUIRE(H != nullptr);
                const double Tc = H->T_critical();
                const double pc = H->p_critical();
                const double Tmax = H->Tmax();
                double p_triple = 1.0;
                try {
                    p_triple = H->p_triple();
                } catch (...) {
                    p_triple = 1.0;
                }
                struct EdgePt
                {
                    double T;
                    double p;
                    const char* note;
                };
                const EdgePt edges[] = {
                  {Tc, pc, "T=Tc, p=pc"},
                  {Tc, std::max(p_triple, pc * 0.5), "T=Tc"},
                  {std::min(Tmax, 1.2 * Tc), 0.5 * pc, "near Tmax-ish"},
                  {0.5 * Tc, std::max(p_triple, 100.0), "low T low p (~p_triple)"},
                  {0.7 * Tc, pc, "subcrit at pc"},
                };
                for (const auto& e : edges) {
                    CAPTURE(e.note, e.T, e.p);
                    const double seed = CoolProp::pxprecond::regime_classified_rho_seed(*H, e.T, e.p);
                    CAPTURE(seed);
                    CHECK(std::isfinite(seed));
                    CHECK(seed > 0.0);
                }
            }
        }
    }
}

TEST_CASE("PXcdj preconditioner: idempotent lazy-init", "[PXcdj][PXcdj_precond]") {
    using clock_t = std::chrono::steady_clock;
    SECTION("Argon: second ensure_melt_curve_table call is essentially free") {
        // Pick a melt-line fluid that is unlikely to have been queried in earlier
        // test_cases (the table_map is process-wide; Water/CO2/Propane/Methane
        // are all queried by the melt-table test_case above, so reusing them
        // here would measure "cached" both times).  Argon has a melting line.
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Argon"));
        auto* H = as_heos(*AS);
        REQUIRE(H != nullptr);
        REQUIRE(H->has_melting_line());

        const auto t0 = clock_t::now();
        CoolProp::pxprecond::ensure_melt_curve_table(*H);
        const auto t1 = clock_t::now();
        CoolProp::pxprecond::ensure_melt_curve_table(*H);
        const auto t2 = clock_t::now();

        const double first_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        const double second_us = std::chrono::duration<double, std::micro>(t2 - t1).count();
        std::printf("[PXcdj_precond] Argon lazy-init: first=%.1f us, second=%.3f us\n", first_us, second_us);

        // First call should genuinely build (PT_INPUTS updates etc.); on this
        // machine that's O(ms) IF the fluid has not already been seeded
        // elsewhere in the test process.  Second call must be a small
        // constant — the map lookup + initialized check.  Documenting actual
        // numbers via the printf above; the assertion is a loose ceiling.
        CHECK(first_us >= 0.0);     // sanity (steady_clock guarantees non-decreasing)
        CHECK(second_us < 1000.0);  // < 1 ms — loose, would pass even under heavy load

        // Table contents unchanged across the two calls.
        const auto& tbl = CoolProp::pxprecond::get_melt_curve_table(*H);
        REQUIRE(tbl.has_data);
        const auto p_size = tbl.p.size();
        CoolProp::pxprecond::ensure_melt_curve_table(*H);
        const auto& tbl2 = CoolProp::pxprecond::get_melt_curve_table(*H);
        REQUIRE(tbl2.p.size() == p_size);
        for (std::size_t i = 0; i < p_size; ++i) {
            CHECK(tbl.p[i] == tbl2.p[i]);
            CHECK(tbl.rho_melt[i] == tbl2.rho_melt[i]);
        }
    }

    SECTION("concurrent ensure + seed (4 threads, Water)") {
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        auto* H = as_heos(*AS);
        REQUIRE(H != nullptr);

        std::atomic<bool> ok{true};
        std::atomic<int> nonpos{0};
        auto worker = [&]() {
            for (int i = 0; i < 200; ++i) {
                try {
                    CoolProp::pxprecond::ensure_melt_curve_table(*H);
                    const double T = 320.0 + 0.1 * i;
                    const double p = 1.0e5 + 1.0e4 * i;
                    const double seed = CoolProp::pxprecond::regime_classified_rho_seed(*H, T, p);
                    if (!(std::isfinite(seed) && seed > 0.0)) {
                        ++nonpos;
                        ok = false;
                    }
                } catch (...) {
                    ok = false;
                }
            }
        };
        std::vector<std::thread> threads;
        for (int t = 0; t < 4; ++t)
            threads.emplace_back(worker);
        for (auto& t : threads)
            t.join();
        CAPTURE(nonpos.load());
        CHECK(ok.load());
    }
}

#endif  // ENABLE_CATCH
