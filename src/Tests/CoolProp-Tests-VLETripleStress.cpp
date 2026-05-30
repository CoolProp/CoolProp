// VLE precision stress test — measures the achievable noise floor of CoolProp's
// saturation flash + EOS evaluation along the saturation curve, with emphasis
// on the triple point where iterative convergence breaks down first.
//
// METHODOLOGY
//   At each (fluid, T) on the saturation curve:
//     1. Run VLE flash twice via QT_INPUTS: once with Q=0 (liquid), once Q=1
//        (vapor).  This produces (rho_L, p_sat_L) and (rho_V, p_sat_V).
//     2. Independently re-evaluate pressure via DmolarT_INPUTS at the liquid
//        and vapor densities.  This computes p_L = p(T, rho_L) and p_V =
//        p(T, rho_V) directly from the EOS, bypassing the solver-internal
//        pressure.
//     3. The metric of interest is |p_L - p_V| / p_avg — the **noise floor of
//        the (solver + EOS evaluation) combination**.  At true VLE in exact
//        arithmetic this is 0; in double precision it is bounded by ULP scatter
//        through the EOS terms.
//     4. Also report Gibbs energy residual |g_L - g_V| / max(|g|, 1) — at VLE
//        in exact arithmetic g_L = g_V (Maxwell criterion).
//
// HOT REGIONS
//   - Triple point: extremely low p (e.g., Propane p_t ≈ 1.7e-10 Pa).  Tiny
//     absolute errors become huge relative errors here.
//   - Near critical: ρ_L ≈ ρ_V, Jacobian becomes near-singular.  Different
//     failure mode — about ill-conditioning, not magnitude.
//   - Mid-saturation: the well-conditioned reference for what "good" looks like.
//
// CSV at /tmp/vle_triple_stress.csv with columns:
//   fluid, T, T_reduced, p_sat, rho_L, rho_V, p_L_eval, p_V_eval,
//   abs_dp, rel_dp, g_L, g_V, rel_dg, wall_us, success
//
// Hidden ([.]) — does not run in CI by default.  Run via:
//   ./build_catch_rel/CatchTestRunner "[vle_triple_stress]"

#include "AbstractState.h"
#include "CoolProp.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <chrono>
#    include <cmath>
#    include <cstdio>
#    include <limits>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

struct FluidSpec
{
    const char* name;
    int n_temps;  // number of temperatures to sample from triple → critical
};

// Cross-section spanning differing molecular weights, EOS structure (Span-
// Wagner vs Lemmon vs Setzmann vs Wagner-Pruss), and triple-point pressures
// (Propane at 1.7e-10 Pa down at one extreme; Methanol at the other).
const std::vector<FluidSpec> kFluids = {
  {"Propane", 30},        // triple T=85.5 K, p ≈ 1.7e-10 Pa — extreme low-p
  {"Water", 30},          // triple T=273.16 K, p ≈ 611 Pa — moderate
  {"Methane", 30},        // triple T=90.7 K, p ≈ 11.7 kPa
  {"Nitrogen", 30},       // triple T=63.2 K, p ≈ 12.5 kPa
  {"CarbonDioxide", 30},  // triple T=216.6 K, p ≈ 518 kPa
  {"Ammonia", 30},        // triple T=195.5 K, p ≈ 6.1 kPa
  {"R134a", 30},          // triple T=169.85 K
  {"Methanol", 30},       // pseudo-pure check; triple T=175.6 K
  {"R125", 30},           // tau_mi_in_u user; triple T=172.5 K
  {"Ethane", 30},         // triple T=90.4 K, p ≈ 1.1 Pa — very low
};

struct StressRow
{
    std::string fluid;
    double T = 0, T_reduced = 0, p_sat = 0;
    double rho_L = 0, rho_V = 0;
    double p_L_eval = 0, p_V_eval = 0;
    double abs_dp = 0, rel_dp = 0;  // as-converged by CoolProp (default 1e-10 tol)
    double g_L = std::numeric_limits<double>::quiet_NaN();
    double g_V = std::numeric_limits<double>::quiet_NaN();
    double rel_dg = 0;
    // After test-internal refinement: keep iterating Akasaka beyond CoolProp's
    // hardcoded 1e-10 tolerance until the residual stops shrinking (true noise
    // floor of the (solver + EOS) combination).
    double rho_L_refined = 0, rho_V_refined = 0;
    double abs_dp_refined = 0, rel_dp_refined = 0;
    int refine_iters = 0;
    double wall_us = 0;
    bool success = false;
    const char* fail_stage = "";
};

// Returns a vector of N temperatures sampled from T_triple to T_critical,
// emphasizing the triple-point side via a cubic warp.
std::vector<double> sample_temperatures(double T_t, double T_c, int N) {
    std::vector<double> Ts;
    Ts.reserve(N);
    const double t_margin_lo = 1e-4;  // stay infinitesimally inside the triple
    const double t_margin_hi = 5e-2;  // stay safely below critical
    for (int i = 0; i < N; ++i) {
        const double u = static_cast<double>(i) / static_cast<double>(N - 1);
        // Cubic warp clusters samples near u=0 (triple side)
        const double warped = u * u * u;
        const double frac = t_margin_lo + (1.0 - t_margin_lo - t_margin_hi) * warped;
        Ts.push_back(T_t + frac * (T_c - T_t));
    }
    return Ts;
}

StressRow measure_one(const std::string& fluid, double T) {
    StressRow r;
    r.fluid = fluid;
    r.T = T;
    const auto t0 = std::chrono::steady_clock::now();
    try {
        auto AS_L = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        const double T_t = AS_L->Ttriple();
        const double T_c = AS_L->T_critical();
        r.T_reduced = (T - T_t) / (T_c - T_t);

        // Phase L
        AS_L->update(CoolProp::QT_INPUTS, 0.0, T);
        r.rho_L = AS_L->rhomolar();
        r.p_sat = AS_L->p();
        try {
            r.g_L = AS_L->gibbsmolar();
        } catch (...) {
            r.g_L = std::numeric_limits<double>::quiet_NaN();
        }

        // Phase V
        auto AS_V = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        AS_V->update(CoolProp::QT_INPUTS, 1.0, T);
        r.rho_V = AS_V->rhomolar();
        try {
            r.g_V = AS_V->gibbsmolar();
        } catch (...) {
            r.g_V = std::numeric_limits<double>::quiet_NaN();
        }

        // Independent evaluation: p(T, rho) from EOS at each density
        auto AS_eval = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        AS_eval->update(CoolProp::DmolarT_INPUTS, r.rho_L, T);
        r.p_L_eval = AS_eval->p();
        AS_eval->update(CoolProp::DmolarT_INPUTS, r.rho_V, T);
        r.p_V_eval = AS_eval->p();

        r.abs_dp = std::abs(r.p_L_eval - r.p_V_eval);
        const double p_scale = std::max(std::abs(r.p_L_eval + r.p_V_eval) * 0.5, 1.0);
        r.rel_dp = r.abs_dp / p_scale;

        if (std::isfinite(r.g_L) && std::isfinite(r.g_V)) {
            const double g_scale = std::max(std::abs(r.g_L + r.g_V) * 0.5, 1.0);
            r.rel_dg = std::abs(r.g_L - r.g_V) / g_scale;
        } else {
            r.rel_dg = std::numeric_limits<double>::quiet_NaN();
        }

        // ---- Refinement: keep pushing past CoolProp's 1e-10 cutoff ---------
        // Production saturation_T_pure_Akasaka stops when error > 1e-10 (see
        // src/Backends/Helmholtz/VLERoutines.cpp:910).  Take its converged
        // state as a starting point, then iterate the same Akasaka step in
        // double precision until the residual stops shrinking — that's the
        // true noise floor.  Reuses AS_L / AS_V via DmolarT_INPUTS updates.
        //
        // Akasaka formulation (per Akasaka 2008):
        //   delta_L = rho_L / rho_r,  delta_V = rho_V / rho_r
        //   J(δ) = δ (1 + δ · ∂α_r/∂δ)
        //   K(δ) = δ · ∂α_r/∂δ + α_r + log(δ)
        //   residual = sqrt((K_L - K_V)^2 + (J_L - J_V)^2)
        //   Newton step solves linearized 2x2 system in (δ_L, δ_V)
        const double rho_r = AS_L->rhomolar_reducing();
        double rL = r.rho_L, rV = r.rho_V;
        double prev_res = std::numeric_limits<double>::infinity();
        for (int it = 0; it < 80; ++it) {
            try {
                AS_L->update(CoolProp::DmolarT_INPUTS, rL, T);
                AS_V->update(CoolProp::DmolarT_INPUTS, rV, T);
            } catch (...) {
                break;
            }
            const double dL = rL / rho_r, dV = rV / rho_r;
            const double aL = AS_L->alphar(), aV = AS_V->alphar();
            const double dad_L = AS_L->dalphar_dDelta(), dad_V = AS_V->dalphar_dDelta();
            const double d2ad_L = AS_L->d2alphar_dDelta2(), d2ad_V = AS_V->d2alphar_dDelta2();

            const double JL = dL * (1 + dL * dad_L);
            const double JV = dV * (1 + dV * dad_V);
            const double KL = dL * dad_L + aL + std::log(dL);
            const double KV = dV * dad_V + aV + std::log(dV);

            const double res = std::sqrt((KL - KV) * (KL - KV) + (JL - JV) * (JL - JV));
            // Convergence: residual stopped shrinking (hit noise floor).  Use
            // a tight ratio (0.999) to allow finite-precision tail.
            if (it > 2 && res >= prev_res * 0.999) {
                r.refine_iters = it;
                break;
            }
            prev_res = res;

            const double dJL = 1 + 2 * dL * dad_L + dL * dL * d2ad_L;
            const double dJV = 1 + 2 * dV * dad_V + dV * dV * d2ad_V;
            const double dKL = 2 * dad_L + dL * d2ad_L + 1.0 / dL;
            const double dKV = 2 * dad_V + dV * d2ad_V + 1.0 / dV;
            const double DELTA = dJV * dKL - dJL * dKV;
            if (std::abs(DELTA) < 1e-30) break;

            const double stepL = ((KV - KL) * dJV - (JV - JL) * dKV) / DELTA;
            const double stepV = ((KV - KL) * dJL - (JV - JL) * dKL) / DELTA;

            const double dL_new = dL + stepL;
            const double dV_new = dV + stepV;
            if (!(dL_new > 1.0 && dV_new < 1.0 && dV_new > 0)) {
                // Step would leave physical regime — accept current state as
                // the noise-floor converged value.
                r.refine_iters = it;
                break;
            }
            rL = dL_new * rho_r;
            rV = dV_new * rho_r;
        }
        r.rho_L_refined = rL;
        r.rho_V_refined = rV;

        // Final p evaluation at refined densities
        try {
            AS_eval->update(CoolProp::DmolarT_INPUTS, r.rho_L_refined, T);
            const double pL_ref = AS_eval->p();
            AS_eval->update(CoolProp::DmolarT_INPUTS, r.rho_V_refined, T);
            const double pV_ref = AS_eval->p();
            r.abs_dp_refined = std::abs(pL_ref - pV_ref);
            const double p_scale_ref = std::max(std::abs(pL_ref + pV_ref) * 0.5, 1.0);
            r.rel_dp_refined = r.abs_dp_refined / p_scale_ref;
        } catch (...) {
            r.abs_dp_refined = r.abs_dp;
            r.rel_dp_refined = r.rel_dp;
        }

        r.success = true;
    } catch (CoolProp::CoolPropBaseError& e) {
        r.success = false;
        r.fail_stage = "CoolPropBaseError";
    } catch (std::exception& e) {
        r.success = false;
        r.fail_stage = "std::exception";
    } catch (...) {
        r.success = false;
        r.fail_stage = "unknown";
    }
    const auto t1 = std::chrono::steady_clock::now();
    r.wall_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    return r;
}

}  // namespace

TEST_CASE("VLE saturation-curve precision stress (triple → critical)", "[vle_triple_stress][.]") {
    std::FILE* csv = std::fopen("/tmp/vle_triple_stress.csv", "w");
    if (csv) {
        std::fprintf(csv, "fluid,T,T_reduced,p_sat,rho_L,rho_V,p_L_eval,p_V_eval,"
                          "abs_dp_Pa,rel_dp,g_L,g_V,rel_dg,"
                          "rho_L_refined,rho_V_refined,abs_dp_refined_Pa,rel_dp_refined,refine_iters,"
                          "wall_us,success,fail_stage\n");
    }

    std::printf("\n=== VLE precision stress test (CoolProp default tol vs noise-floor refined) ===\n");
    std::printf("Worst rel_dp per fluid:  as-converged (CoolProp 1e-10 cutoff) vs noise-floor (refined)\n");
    std::printf("%-15s %12s %14s %14s %8s\n", "fluid", "p_sat_min", "rel_dp (default)", "rel_dp (refined)", "max_it");
    std::printf("%-15s %12s %14s %14s %8s\n", "-----", "---------", "----------------", "----------------", "------");

    struct FluidSummary
    {
        std::string name;
        double worst_rel_dp = 0;
        double worst_rel_dp_refined = 0;
        double worst_rel_dg = 0;
        double worst_T_reduced_dp = 0;  // T_reduced at which worst rel_dp occurs
        int max_refine_iters = 0;
        std::size_t n_success = 0, n_fail = 0;
        double min_p_sat = std::numeric_limits<double>::infinity();
    };
    std::vector<FluidSummary> summaries;

    for (const auto& fs : kFluids) {
        FluidSummary fsum;
        fsum.name = fs.name;

        double T_t = 0, T_c = 0;
        try {
            auto probe = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fs.name));
            T_t = probe->Ttriple();
            T_c = probe->T_critical();
        } catch (...) {
            std::printf("%-15s SKIPPED (factory or Ttriple/T_critical failed)\n", fs.name);
            continue;
        }

        const auto Ts = sample_temperatures(T_t, T_c, fs.n_temps);
        for (double T : Ts) {
            const auto r = measure_one(fs.name, T);
            if (csv) {
                std::fprintf(csv,
                             "%s,%.6e,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,"
                             "%.6e,%.6e,%.6e,%.6e,%d,%.2f,%d,%s\n",
                             r.fluid.c_str(), r.T, r.T_reduced, r.p_sat, r.rho_L, r.rho_V, r.p_L_eval, r.p_V_eval, r.abs_dp, r.rel_dp, r.g_L, r.g_V,
                             r.rel_dg, r.rho_L_refined, r.rho_V_refined, r.abs_dp_refined, r.rel_dp_refined, r.refine_iters, r.wall_us,
                             r.success ? 1 : 0, r.fail_stage);
            }
            if (r.success) {
                ++fsum.n_success;
                if (r.rel_dp > fsum.worst_rel_dp) {
                    fsum.worst_rel_dp = r.rel_dp;
                    fsum.worst_T_reduced_dp = r.T_reduced;
                }
                if (r.rel_dp_refined > fsum.worst_rel_dp_refined) fsum.worst_rel_dp_refined = r.rel_dp_refined;
                if (r.refine_iters > fsum.max_refine_iters) fsum.max_refine_iters = r.refine_iters;
                if (std::isfinite(r.rel_dg) && r.rel_dg > fsum.worst_rel_dg) fsum.worst_rel_dg = r.rel_dg;
                if (r.p_sat > 0 && r.p_sat < fsum.min_p_sat) fsum.min_p_sat = r.p_sat;
            } else {
                ++fsum.n_fail;
            }
        }
        std::printf("%-15s %12.3e %14.3e %14.3e %8d%s\n", fsum.name.c_str(), fsum.min_p_sat, fsum.worst_rel_dp, fsum.worst_rel_dp_refined,
                    fsum.max_refine_iters, fsum.n_fail > 0 ? "  FAILS" : "");
        if (fsum.n_fail > 0) {
            std::printf("  ↳ %zu/%zu VLE flashes FAILED for this fluid\n", fsum.n_fail, fsum.n_success + fsum.n_fail);
        }
        summaries.push_back(fsum);
    }
    if (csv) std::fclose(csv);

    // Print aggregate summary
    std::printf("\n=== Aggregate summary (worst rel_dp across all (fluid, T)) ===\n");
    double overall_worst_dp = 0, overall_worst_dg = 0;
    std::string overall_dp_at, overall_dg_at;
    std::size_t total_success = 0, total_fail = 0;
    for (const auto& s : summaries) {
        total_success += s.n_success;
        total_fail += s.n_fail;
        if (s.worst_rel_dp > overall_worst_dp) {
            overall_worst_dp = s.worst_rel_dp;
            overall_dp_at = s.name + std::string(" at T_red=") + std::to_string(s.worst_T_reduced_dp);
        }
        if (s.worst_rel_dg > overall_worst_dg) {
            overall_worst_dg = s.worst_rel_dg;
            overall_dg_at = s.name;
        }
    }
    std::printf("Worst |Δp|/p across %zu (fluid, T) probes: %.3e (%s)\n", total_success + total_fail, overall_worst_dp, overall_dp_at.c_str());
    std::printf("Worst |Δg|/g across same: %.3e (%s)\n", overall_worst_dg, overall_dg_at.c_str());
    std::printf("VLE flash success rate: %zu/%zu = %.2f%%\n", total_success, total_success + total_fail,
                100.0 * static_cast<double>(total_success) / static_cast<double>(total_success + total_fail));
}

#endif  // ENABLE_CATCH
