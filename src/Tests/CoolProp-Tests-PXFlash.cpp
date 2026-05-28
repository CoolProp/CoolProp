// Focused reliability regression test for the P+{H,S,U} single-phase flash
// (HSU_P_flash_singlephase_Brent in src/Backends/Helmholtz/FlashRoutines.cpp).
//
// Runs in the default suite (tag [PXflash], NOT [.]-hidden).  Run explicitly:
//   ./CatchTestRunner "[PXflash]"
//
// Background (bd CoolProp-0nx): the solver left HEOS._p corrupted when the inner
// PT density solve threw mid-Newton; the exception-fallback guard then branched
// on the garbage pressure and re-threw instead of running the 2-D Newton
// fallback.  Two supercritical-cold (p > p_crit) Nitrogen/PS states failed.
// The fix saves p_target at entry and lets the fallback engage for p > p_crit.
//
// Input-pair argument order (verified against include/DataStructures.h):
//   HmolarP_INPUTS : update(HmolarP_INPUTS, h_J_mol, p_Pa)   (H first, P second)
//   PSmolar_INPUTS : update(PSmolar_INPUTS, p_Pa, s_J_mol_K) (P first, S second)
//   PUmolar_INPUTS : update(PUmolar_INPUTS, p_Pa, u_J_mol)   (P first, U second)

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <chrono>
#    include <cmath>
#    include <cstdio>
#    include <memory>
#    include <string>
#    include <vector>

#    include "AbstractState.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "DataStructures.h"

// ---------------------------------------------------------------------------
// The two previously-failing Nitrogen states.  Both are supercritical-cold
// (p > p_crit ~ 3.396 MPa), the compressed-liquid corner.  Before the fix the
// P,S round-trip threw; the fix must make them solve and recover T and rho.
// ---------------------------------------------------------------------------
TEST_CASE("P+X supercritical-cold Nitrogen reliability (was failing)", "[PXflash]") {
    struct State
    {
        double T;  // K
        double p;  // Pa
    };
    const State pt = GENERATE(State{86.35, 4.754e6}, State{90.20, 7.768e6});
    CAPTURE(pt.T, pt.p);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));

    // Establish the reference state at (T, p) and read its entropy.
    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, pt.p, pt.T));
    const double s = ref->smolar();
    const double rho = ref->rhomolar();
    REQUIRE(std::isfinite(s));
    REQUIRE(std::isfinite(rho));

    // The P,S flash that used to throw must now solve and round-trip.
    REQUIRE_NOTHROW(wrk->update(CoolProp::PSmolar_INPUTS, pt.p, s));
    CHECK(wrk->T() == Catch::Approx(pt.T).epsilon(1e-5));
    CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
}

// ---------------------------------------------------------------------------
// Small, fast P+X round-trip guard across a few fluids/phases.  Sets a known
// (T, p) state, reads h/s/u, then runs each of the PH / PS / PU flashes and
// confirms T and rho recover.  Keeps CI cost low (a handful of points).
// ---------------------------------------------------------------------------
TEST_CASE("P+X single-phase round-trips (PH/PS/PU)", "[PXflash]") {
    struct Case
    {
        const char* fluid;
        double T;  // K
        double p;  // Pa
        const char* label;
    };
    const Case c = GENERATE(Case{"Water", 600.0, 1.0e5, "Water gas"},              // superheated vapor
                            Case{"Water", 320.0, 5.0e6, "Water liquid"},           // compressed liquid
                            Case{"CarbonDioxide", 320.0, 1.0e7, "CO2 supercrit"},  // p,T > critical
                            Case{"R134a", 250.0, 3.0e5, "R134a liquid"},           // subcooled liquid
                            Case{"R134a", 350.0, 2.0e5, "R134a gas"});             // superheated vapor
    CAPTURE(c.label, c.fluid, c.T, c.p);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));

    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, c.p, c.T));
    const double h = ref->hmolar();
    const double s = ref->smolar();
    const double u = ref->umolar();
    const double rho = ref->rhomolar();
    REQUIRE(std::isfinite(h));
    REQUIRE(std::isfinite(s));
    REQUIRE(std::isfinite(u));

    SECTION("PH") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::HmolarP_INPUTS, h, c.p));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
    SECTION("PS") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PSmolar_INPUTS, c.p, s));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
    SECTION("PU") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PUmolar_INPUTS, c.p, u));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
}

// ---------------------------------------------------------------------------
// Non-EOS overhead profile of a single Propane P+H flash (bd CoolProp-0nx).
// Tagged [PXFLASH_PROFILE][.] so it is hidden from the default suite; invoke
// explicitly:   ./CatchTestRunner "[PXFLASH_PROFILE]"
//
// Times each component of one outer step of HSU_P_flash_singlephase_Brent in
// isolation, plus the full P+H solve, at the representative Propane gas state
// p = 1 MPa, T = 350 K (round-tripped via HmolarP_INPUTS).  Goal: localize the
// ~15 us per-solve overhead beyond raw EOS to specific levers.  See
// CoolProp-0nx for the analysis output.
// ---------------------------------------------------------------------------
TEST_CASE("PXflash profile - Propane gas (1 MPa, 350 K)", "[PXFLASH_PROFILE][.]") {
    using clk = std::chrono::steady_clock;
    using ns = std::chrono::duration<double, std::nano>;

    // Establish the (T, p) reference state and read h.
    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Propane"));
    auto* ref = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    REQUIRE(ref);
    const double p = 1.0e6, T = 350.0;
    ref->update(CoolProp::PT_INPUTS, p, T);
    const double h = ref->hmolar();
    const double rho = ref->rhomolar();
    const auto mole_fractions = ref->get_mole_fractions();
    const double Tr = ref->get_reducing_state().T;
    const double rhor = ref->get_reducing_state().rhomolar;
    const double tau = Tr / T;
    const double delta = rho / rhor;

    INFO("p = " << p << " Pa, T = " << T << " K, rho = " << rho << " mol/m^3");
    INFO("h = " << h << " J/mol, tau = " << tau << ", delta = " << delta);

    // Working state for the timed kernels (separate from `ref`).
    auto wrksp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Propane"));
    auto* wrk = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(wrksp.get());
    REQUIRE(wrk);

    // Number of inner reps per timed block.  Outer trials average over noise.
    constexpr int REPS = 100000;
    constexpr int TRIALS = 5;

    // ---- (1) Raw Helmholtz residual kernel at fixed (tau, delta) -----------
    // calc_alphar_deriv_nocache calls residual_helmholtz->all internally.
    auto time_block = [&](const char* name, auto&& thunk) {
        double best = 1e300;
        double sum = 0.0;
        for (int t = 0; t < TRIALS; ++t) {
            const auto t0 = clk::now();
            for (int i = 0; i < REPS; ++i) thunk(i);
            const auto t1 = clk::now();
            const double per = ns(t1 - t0).count() / REPS;
            sum += per;
            if (per < best) best = per;
        }
        const double mean = sum / TRIALS;
        std::printf("  %-50s  best %8.1f ns   mean %8.1f ns\n", name, best, mean);
        return mean;
    };

    std::printf("\n[PXflash profile] Propane gas, p=1 MPa, T=350 K   (%d reps x %d trials, Release -O3)\n", REPS, TRIALS);

    // Perturb tau slightly to defeat any per-call cache (calc_alphar_deriv_nocache
    // is documented as bypassing the cache, but we still avoid bit-identical input).
    volatile double sink = 0.0;
    const double raw = time_block("raw alphar all (residual_helmholtz->all)", [&](int i) {
        const double t_eps = 1.0 + 1e-12 * (i & 0xff);
        sink += wrk->calc_alphar_deriv_nocache(0, 0, mole_fractions, tau * t_eps, delta);
    });

    // ---- (2) update_DmolarT_direct(rho, T) ---------------------------------
    const double upd_DT = time_block("update_DmolarT_direct(rho, T)", [&](int i) {
        const double rho_eps = rho * (1.0 + 1e-12 * (i & 0xff));
        wrk->update_DmolarT_direct(rho_eps, T);
        sink += wrk->p();
    });

    // ---- (3) update_TP_guessrho (warm) -------------------------------------
    // Seed the inner Newton with a rho already very close to the root.
    const double upd_TPwarm = time_block("update_TP_guessrho(T, p, rho) [WARM]", [&](int i) {
        const double rho_guess = rho * (1.0 + 1e-10 * (i & 0xff));
        wrk->update_TP_guessrho(T, p, rho_guess);
        sink += wrk->p();
    });

    // ---- (4) update(PT_INPUTS, p, T)  [COLD: full inner Newton] ------------
    const double upd_PT = time_block("update(PT_INPUTS, p, T)        [COLD]", [&](int i) {
        const double T_eps = T * (1.0 + 1e-12 * (i & 0xff));
        wrk->update(CoolProp::PT_INPUTS, p, T_eps);
        sink += wrk->rhomolar();
    });

    // ---- (5) keyed_output(iHmolar) alone after a state-set -----------------
    wrk->update_DmolarT_direct(rho, T);
    const double key_h = time_block("keyed_output(iHmolar) after state set", [&](int /*i*/) {
        // The first call hits an uncached cache slot; subsequent reps in the
        // same wrk should hit it.  We re-set the state every K reps to mimic
        // the per-outer-step cache miss the solver actually pays.
        // Cache misses are what the solver experiences:
        wrk->update_DmolarT_direct(rho, T);
        sink += wrk->keyed_output(CoolProp::iHmolar);
    });

    // ---- (6) Solver inner step: update_TP_guessrho + keyed_output ---------
    const double inner_step = time_block("resid.call(T) approx: TP_warm + keyed_h", [&](int i) {
        const double T_eps = T * (1.0 + 1e-12 * (i & 0xff));
        wrk->update_TP_guessrho(T_eps, p, rho);
        sink += wrk->keyed_output(CoolProp::iHmolar);
    });

    // Cold variant: full PT cold + keyed_h
    const double inner_step_cold = time_block("resid.call(T) approx: PT_cold + keyed_h", [&](int i) {
        const double T_eps = T * (1.0 + 1e-12 * (i & 0xff));
        wrk->update(CoolProp::PT_INPUTS, p, T_eps);
        sink += wrk->keyed_output(CoolProp::iHmolar);
    });

    // ---- (7) Full P+H solve -----------------------------------------------
    const double full = time_block("FULL update(HmolarP_INPUTS, h, p)", [&](int i) {
        const double h_eps = h * (1.0 + 1e-12 * (i & 0xff));
        wrk->update(CoolProp::HmolarP_INPUTS, h_eps, p);
        sink += wrk->T();
    });

    // ---- (7b) Just update() with no flash work for comparison: PT cold ----
    // (Already measured above, but include here as the outer-frame baseline.)
    // (7c) Full PT solve including phase determination via standard update().
    // (7d) get_T_from_p superancillary inversion alone (called twice per solve:
    //      once in p_phase_determination_pure_or_pseudopure, once in HSU_P_flash
    //      for the gas-phase Tmin.).
    auto superanc_ptr = wrk->get_superanc();
    REQUIRE(superanc_ptr);
    auto& superanc = *superanc_ptr;
    const double sa_TfromP = time_block("superanc.get_T_from_p(p)              ", [&](int i) {
        const double p_eps = p * (1.0 + 1e-12 * (i & 0xff));
        sink += superanc.get_T_from_p(p_eps);
    });
    const double sa_evalsat = time_block("superanc.eval_sat(Tsat, 'D', 1)       ", [&](int i) {
        const double T_eps = T * (1.0 + 1e-12 * (i & 0xff));
        sink += superanc.eval_sat(T_eps, 'D', 1);
    });

    // ---- (7e) Cost of clear() alone ---------------------------------------
    const double t_clear = time_block("wrk->clear() (cache + bulk fields)     ", [&](int /*i*/) {
        wrk->clear();
        sink += 1.0;
    });

    // ---- (7f) keyed_output(iHmolar) WITHOUT re-state-set (cache HIT) ------
    wrk->update_DmolarT_direct(rho, T);
    sink += wrk->keyed_output(CoolProp::iHmolar);  // warm cache
    const double key_h_warm = time_block("keyed_output(iHmolar) WARM cache hit   ", [&](int /*i*/) {
        sink += wrk->keyed_output(CoolProp::iHmolar);
    });

    std::printf("  keyed_output(iHmolar) WARM (no re-set)  = %8.1f\n", key_h_warm);
    std::printf("  -> implied keyed_output dispatch cost   = %8.1f   (key_h_warm - 0)\n", key_h_warm);
    std::printf("  clear() alone                           = %8.1f\n", t_clear);

    // ---- (8) solver_rho_Tp_SRK alone (initial-guess cost) ------------------
    // Need a fresh state (the SRK helper relies on EOS().reduce on components).
    wrk->update(CoolProp::PT_INPUTS, p, T);  // prime reducing-state cache
    const double srk = time_block("solver_rho_Tp_SRK(T, p, iphase_gas)", [&](int i) {
        const double T_eps = T * (1.0 + 1e-12 * (i & 0xff));
        sink += wrk->solver_rho_Tp_SRK(T_eps, p, CoolProp::iphase_gas);
    });

    // ---- Cross-checks and decomposition ------------------------------------
    std::printf("\n[PXflash profile] decomposition (ns, mean):\n");
    std::printf("  raw alphar kernel                       = %8.1f\n", raw);
    std::printf("  update_DmolarT_direct - raw             = %8.1f   (pre/post_update bookkeeping)\n", upd_DT - raw);
    std::printf("  update_TP_guessrho [WARM] - upd_DT      = %8.1f   (inner Newton: ~K*upd_DT + Householder setup)\n", upd_TPwarm - upd_DT);
    std::printf("  update(PT) [COLD] - update_TP_guessrho  = %8.1f   (SRK guess + ancillary + extra Newton iters)\n", upd_PT - upd_TPwarm);
    std::printf("  solver_rho_Tp_SRK                       = %8.1f   (initial-guess cubic root)\n", srk);
    std::printf("  keyed_output(iHmolar) (incl. state set) = %8.1f\n", key_h);
    std::printf("  resid.call WARM = TP_warm + keyed_h     = %8.1f   (target: ~%0.1f)\n", inner_step, upd_TPwarm + (key_h - upd_DT));
    std::printf("  resid.call COLD = PT_cold + keyed_h     = %8.1f   (target: ~%0.1f)\n", inner_step_cold, upd_PT + (key_h - upd_DT));
    std::printf("  FULL solve                              = %8.1f\n", full);

    // Approximate iteration-count back-calculation
    // Prior data: ~9 outer iters; first 2 cold, rest warm (per HSU_P_flash_singlephase_Brent).
    const double approx_9outer = 2 * inner_step_cold + 7 * inner_step;
    const double toms_overhead = full - approx_9outer;
    std::printf("  ~2*cold + 7*warm                        = %8.1f\n", approx_9outer);
    std::printf("  TOMS748 outer bookkeeping (full - calc) = %8.1f\n", toms_overhead);

    std::printf("\n[PXflash profile] full-solve stacked decomposition (estimated):\n");
    const double full_raw_eos = full > 0 ? (21.0 * raw) : 0.0;  // ~21 EOS evals per solve
    std::printf("  ~21 raw EOS evals                       = %8.1f  (%.0f%%)\n", full_raw_eos, 100.0 * full_raw_eos / full);
    const double full_pre_post = full > 0 ? (21.0 * (upd_DT - raw)) : 0.0;
    std::printf("  ~21 pre/post_update bookkeeping         = %8.1f  (%.0f%%)\n", full_pre_post, 100.0 * full_pre_post / full);
    const double full_srk = 2.0 * srk;  // SRK guess per cold call (first 2 iters cold)
    std::printf("  2 SRK initial guesses (2 cold iters)    = %8.1f  (%.0f%%)\n", full_srk, 100.0 * full_srk / full);
    const double full_keyed = 9.0 * (key_h - upd_DT);
    std::printf("  9 keyed_output dispatches               = %8.1f  (%.0f%%)\n", full_keyed, 100.0 * full_keyed / full);
    const double accounted = full_raw_eos + full_pre_post + full_srk + full_keyed;
    std::printf("  unaccounted (TOMS748 outer + misc)      = %8.1f  (%.0f%%)\n", full - accounted, 100.0 * (full - accounted) / full);

    // Defeat sink elision.
    REQUIRE(std::isfinite(static_cast<double>(sink)) == std::isfinite(static_cast<double>(sink)));

    // Sanity: full solve should land in the ballpark of the prior baseline (~22 us
    // per outer P+H solve at this state).  Loose bounds; this is a profile, not a
    // pass/fail.
    CHECK(full > 5e3);   // > 5 us
    CHECK(full < 60e3);  // < 60 us
}

#endif  // ENABLE_CATCH
