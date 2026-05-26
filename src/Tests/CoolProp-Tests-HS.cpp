// Characterization grid for SINGLE-PHASE H,S (HmolarSmolar_INPUTS) flashes.
//
// HS is the last input family not superancillary-accelerated.  Two-phase HS is
// a 1D-in-T Qh-Qs rootfind (mirror of the D+X happy path); single-phase HS is a
// genuine 2D (T, rho) inversion whose current solver seeds from a RANDOM guess
// (HS_flash_generate_TP_singlephase_guess), the root cause of #2022.
//
// This builds a dense (log p, linear T) grid over hard fluids -- every PT point
// is single-phase, so (h, s) <-> (T, rho) is a bijection and T/rho recovery is a
// clean, unambiguous check (unlike D+X, no non-uniqueness).  It maps the current
// HS_flash failure surface as the baseline to beat with the continuation method
// (bd CoolProp-j3n HS child).
//
// Opt-in (the current solver fails on much of this):
//     CatchTestRunner "[HS_grid]"
//     CatchTestRunner "[HS_grid]" -c Water
// Grid is env-tunable: HS_GRID_NT / HS_GRID_NP (default 80 x 80).

#include "AbstractState.h"
#include "CoolProp.h"
#include "DataStructures.h"
#include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "../Backends/Helmholtz/FlashRoutines.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include "boost/math/tools/toms748_solve.hpp"
#    include "CPstrings.h"
#    include <array>
#    include <chrono>
#    include <cmath>
#    include <cstdlib>
#    include <fstream>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = (n == 1) ? a : a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return v;
}

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

// Outcome of the prototype continuation solver.
struct HSResult
{
    bool ok = false;
    double T = 0, rho = 0;
    int nsteps = 0;  // continuation steps used
    int nevals = 0;  // EOS state evaluations (update_DmolarT_direct calls)
    int method = 0;  // cascade leg that produced the accepted result: 1=sat 2=isentrope 3=departure
    double anchorT = 0, anchorRho = 0;
};

// Impose single-phase so update_DmolarT_direct + h/s/derivative queries take the
// single-phase EOS path everywhere (including metastable (T,rho) on a continuation
// path); RAII-restore on every exit.
struct PhaseGuard
{
    CoolProp::HelmholtzEOSMixtureBackend& AS;
    explicit PhaseGuard(CoolProp::HelmholtzEOSMixtureBackend& AS) : AS(AS) {
        AS.specify_phase(CoolProp::iphase_gas);
    }
    ~PhaseGuard() {
        AS.unspecify_phase();
    }
    PhaseGuard(const PhaseGuard&) = delete;
    PhaseGuard& operator=(const PhaseGuard&) = delete;
    PhaseGuard(PhaseGuard&&) = delete;
    PhaseGuard& operator=(PhaseGuard&&) = delete;
};

// Core (T, w=ln rho) homotopy corrector shared by the saturation-anchor and
// supercritical-isentrope solvers.  Homotopes (h,s) linearly from the anchor's
// values to the target with adaptive subdivision; the (dp/drho)_T>0 guard keeps
// the corrector on the mechanically stable branch (det J single-signed there, so
// the stable (h,s)->(T,rho) root is unique).  Caller must have imposed single
// phase (PhaseGuard).  Returns true and fills R.{T,rho,nsteps} on success.
bool continue_HS(CoolProp::HelmholtzEOSMixtureBackend& AS, double T0, double rho0, double h_t, double s_t, HSResult& R) {
    using namespace CoolProp;
    const double Rgas = AS.gas_constant();
    const double Tsc = AS.T_critical();
    auto eval = [&](double T, double rho) {
        AS.update_DmolarT_direct(rho, T);
        ++R.nevals;
    };
    eval(T0, rho0);
    const double h0 = AS.hmolar(), s0 = AS.smolar();
    // Dimensional scales (NOT max(|h|,1): water's h passes through 0 at its
    // reference state, which would make a relative tolerance unreachable).
    const double hscale = Rgas * Tsc;  // molar energy scale
    const double sscale = Rgas;        // molar entropy scale

    const char* nmax_env = std::getenv("HS_CONT_NMAX");
    const int Nmax = (nmax_env != nullptr) ? std::atoi(nmax_env) : 128;
    for (int N = 1; N <= Nmax; N *= 2) {
        double T = T0, w = std::log(rho0);
        bool failed = false;
        for (int k = 1; k <= N && !failed; ++k) {
            const double lam = static_cast<double>(k) / N;
            const double ht = h0 + lam * (h_t - h0), st = s0 + lam * (s_t - s0);
            bool conv = false;
            double best_norm = 1e300, best_T = T, best_w = w;
            for (int iter = 0; iter < 40; ++iter) {
                const double rho = std::exp(w);
                eval(T, rho);
                const double rh = AS.hmolar() - ht, rs = AS.smolar() - st;
                const double norm = std::abs(rh) / hscale + std::abs(rs) / sscale;
                if (norm < best_norm) {
                    best_norm = norm;
                    best_T = T;
                    best_w = w;
                }
                if (norm < 1e-11) {
                    conv = true;
                    break;
                }
                // Jacobian in (T, w): columns d/dT, d/dw = rho*d/drho.
                const double hT = AS.first_partial_deriv(iHmolar, iT, iDmolar);
                const double hr = AS.first_partial_deriv(iHmolar, iDmolar, iT);
                const double sT = AS.first_partial_deriv(iSmolar, iT, iDmolar);
                const double sr = AS.first_partial_deriv(iSmolar, iDmolar, iT);
                const double a11 = hT, a12 = rho * hr, a21 = sT, a22 = rho * sr;
                const double det = a11 * a22 - a12 * a21;
                if (!std::isfinite(det) || std::abs(det) < 1e-300) break;
                double dT = -(a22 * rh - a12 * rs) / det;
                double dw = -(-a21 * rh + a11 * rs) / det;
                // Damp so each Newton iterate stays near the EOS T-range (and |dw|
                // not explosive); rho<=0 is impossible because w = ln(rho) is the
                // variable.  The lower bound is Tmin with a small relative slack: a
                // hard clamp at Tmin blocks two opposite needs.  Too loose, and a
                // cold-cryogen corrector slides far below Tmin into a SPURIOUS
                // (h,s)-equal root (dense ortho-hydrogen, ~30% below Tmin).  Too
                // tight, and a compressed-liquid homotopy from the triple point
                // cannot follow its (dome-free) path, which folds a few hundredths of
                // a percent below Tmin near the density anomaly (water at ~273 K) and
                // then climbs back to an in-range target.  The 2% slack admits the
                // shallow fold while still blocking the deep spurious slide; accept()
                // guarantees the FINAL root is in [Tmin, Tmax] regardless.
                const double Tlo = AS.Tmin() * (1.0 - 2e-2);
                double f = 1.0;
                while ((T + f * dT < Tlo || T + f * dT > 1.5 * AS.Tmax() || std::abs(f * dw) > 2.0) && f > 1e-6)
                    f *= 0.5;
                // Mechanical-stability guard: keep the iterate where (dp/drho)_T > 0.
                for (int g = 0; g < 8 && f > 1e-3; ++g) {
                    eval(T + f * dT, std::exp(w + f * dw));
                    if (AS.first_partial_deriv(iP, iDmolar, iT) > 0) break;  // mechanically stable
                    f *= 0.5;
                }
                T += f * dT;
                w += f * dw;
            }
            // Accept the best iterate if it reached EOS-level precision even without
            // hitting the tight flag (e.g. the cold saturated-liquid boundary).
            if (!conv && best_norm < 1e-8) {
                T = best_T;
                w = best_w;
                conv = true;
            }
            if (!conv) failed = true;
        }
        if (!failed) {
            R.ok = true;
            R.T = T;
            R.rho = std::exp(w);
            R.nsteps = N;
            return true;
        }
    }
    return false;
}

// Prototype: single-phase HS solve via homotopy continuation in (T, ln rho)
// from a superancillary saturation anchor.  Working in w = ln(rho) makes rho<0
// impossible and conditions the entropy direction (s = s_ig(T) - R*ln rho +
// dep).  The anchor is the saturated state whose entropy equals the target's
// (liquid branch if the target sits at/below critical entropy, else vapor),
// which keeps the straight (h,s) homotopy on one side of the dome.  Adaptive
// subdivision (double the step count on a failed Newton) provides the guarantee.
HSResult solve_HS_continuation(CoolProp::HelmholtzEOSMixtureBackend& AS, CoolProp::superancillary::SuperAncillary<std::vector<double>>& sa,
                               double h_t, double s_t) {
    using namespace CoolProp;
    HSResult R;
    const double Tc = sa.get_Tcrit_num();
    const double Tmin = sa.get_Tmin();
    const double Tcap = Tc - std::max(1e-4, 1e-6 * Tc);  // stay off the exact critical point

    PhaseGuard phase_guard(AS);

    auto eval = [&](double T, double rho) {
        AS.update_DmolarT_direct(rho, T);
        ++R.nevals;
    };
    // Saturated-branch entropy as a function of T (caller frame): density from
    // the (reference-free) D superancillary, entropy from the EOS.
    auto sat_s = [&](double T, short Q) {
        eval(T, sa.eval_sat(T, 'D', Q));
        return AS.smolar();
    };

    // --- 1. Anchor selection ---
    // Enumerate every saturation temperature whose saturated entropy equals the
    // target's via the superancillary's native S-inversion (get_all_intersections
    // on the caloric 'S' superancillary) -- Chebyshev rootfinding, no EOS calls,
    // returns ALL roots (handles the retrograde / non-monotonic saturated entropy
    // of siloxanes).  Each root may be on either branch, so disambiguate with the
    // SECOND input: keep the candidate whose saturated enthalpy is closest to h_t
    // -- a same-phase, near-target anchor, so the straight (h,s) homotopy stays on
    // one side of the dome.
    const double Rgas = AS.gas_constant();
    AS.ensure_caloric_superancillaries();  // build H/S/U so get_all_intersections('S') works
    const double a = Tmin + 1e-3, b = Tcap;
    const double s_crit = sat_s(Tcap, 0);  // ~critical entropy (both branches meet)
    // Shift s_t into the superancillary's stamped reference frame (see
    // resolve_T_via_superancillary / #2773): s_cache = s_t - R*(a1_cache - a1_caller).
    double s_cache = s_t;
    {
        const auto stamp = sa.get_caloric_alpha0_stamp();
        if (stamp.has_value()) {
            const auto [caller_a1, caller_a2] = CoolProp::FlashRoutines::alpha0_offset_total(AS);
            s_cache = s_t - Rgas * (stamp->first - caller_a1);
        }
    }
    double T0 = 0, rho0 = 0, best_hgap = 1e300;
    bool found = false, near_critical = false;
    auto consider = [&](double Tcand, short Q) {
        const double rho = sa.eval_sat(Tcand, 'D', Q);
        eval(Tcand, rho);
        const double hgap = std::abs(AS.hmolar() - h_t);
        if (std::isfinite(hgap) && hgap < best_hgap) {
            best_hgap = hgap;
            T0 = Tcand;
            rho0 = rho;
            found = true;
        }
    };
    for (const auto& pr : sa.get_all_intersections('S', s_cache, 48, 100, 1e-12)) {
        const double Tcand = pr.first;
        // A matching saturation root that sits in the [Tcap, Tc] critical-exclusion
        // sliver means s_t is ~ s_crit: the target is effectively a critical-isentrope
        // state.  Flag it so the fallback below hands it to the supercritical leg
        // rather than the (wrong) dilute-gas anchor.
        if (Tcand >= b) near_critical = true;
        if (Tcand <= a || Tcand >= b) continue;
        // get_all_intersections concatenates both branches without labelling, so
        // recover the branch this root belongs to (the one whose saturated entropy
        // equals s_cache there) and anchor on THAT branch only.  Trying both
        // branches blindly lets the h-disambiguation pick a wrong-branch anchor.
        const short Qroot = (std::abs(sa.eval_sat(Tcand, 'S', 0) - s_cache) <= std::abs(sa.eval_sat(Tcand, 'S', 1) - s_cache)) ? 0 : 1;
        consider(Tcand, Qroot);
    }
    if (!found) {
        if (near_critical) {
            // The only matching saturation root sits against the critical point
            // (get_all_intersections DOES return it -- e.g. water at s ~ s_crit
            // anchors at T = Tc - 6e-4 K -- but it lies in the [Tcap, Tc] exclusion
            // sliver).  The target is a critical-isentrope state; the dilute-gas
            // anchor below would be wrong (the target is dense, not dilute), making
            // the homotopy thrash adaptive subdivision to its ceiling (~4000 evals)
            // and fail anyway.  Hand it straight to the supercritical isentrope leg,
            // which anchors at T = Tmax (far from the singular critical region) and
            // solves these in ~30 evals.
            return R;  // R.ok == false -> cascade falls through to solve_HS_isentrope
        }
        // s_t genuinely outside both saturation branches (far from critical).
        if (s_t > s_crit) {
            // High entropy -> dilute gas: tiny density (nearly ideal), T from the
            // 1D (monotone) enthalpy inversion h(rho0, T) = h_t.
            rho0 = sa.get_rhocrit_num() * 1e-4;
            auto hres = [&](double T) {
                eval(T, rho0);
                return AS.hmolar() - h_t;
            };
            double Ta = Tmin + 1e-3, Tb = 1.5 * AS.Tmax();
            double fa = hres(Ta), fb = hres(Tb);
            if (fa * fb <= 0) {
                boost::math::uintmax_t it = 60;
                auto [l, r] = boost::math::tools::toms748_solve(hres, Ta, Tb, fa, fb, boost::math::tools::eps_tolerance<double>(30), it);
                T0 = 0.5 * (l + r);
            } else {
                T0 = (std::abs(fa) <= std::abs(fb)) ? Ta : Tb;
            }
        } else {
            // Low entropy -> cold compressed liquid: triple-point saturated liquid.
            T0 = a;
            rho0 = sa.eval_sat(a, 'D', 0);
        }
    }
    R.anchorT = T0;
    R.anchorRho = rho0;

    // --- 2. Continuation in (T, w=ln rho), adaptive subdivision ---
    continue_HS(AS, T0, rho0, h_t, s_t, R);
    return R;
}

// Supercritical-target solver: anchor on the dome-free isotherm at T = Tmax with
// a density chosen so the anchor entropy already equals s_t, then run the shared
// (T, ln rho) corrector.  Because s_anchor == s_t, the (h,s) homotopy holds
// entropy fixed -- the path is the isentrope s = s_t, which for a supercritical
// target stays at T > Tc the whole way down and so never touches the dome.  This
// is the third leg of the cascade: it owns the supercritical points where the
// saturation-anchor line and the ideal-gas departure path both clip the dome.
HSResult solve_HS_isentrope(CoolProp::HelmholtzEOSMixtureBackend& AS, double h_t, double s_t) {
    using namespace CoolProp;
    HSResult R;
    PhaseGuard phase_guard(AS);
    const double Ta = AS.Tmax();
    const double rhoc = AS.rhomolar_critical();
    auto sres = [&](double rho) {
        AS.update_DmolarT_direct(rho, Ta);
        ++R.nevals;
        return AS.smolar() - s_t;
    };
    // Entropy decreases monotonically with density along an isotherm; bracket the
    // anchor density between a near-ideal gas and a dense (liquid-like) state.  Low
    // target entropy (cold, very compressed states) needs a denser anchor than a
    // fixed multiple of rho_crit, so grow the upper bound until it brackets s_t (or
    // the EOS refuses the density) rather than clamping -- a clamped anchor leaves
    // the corrector seeded too far off for the extreme high-pressure corner.
    const double ra = 1e-6 * rhoc;
    const double ga = sres(ra);
    double rb = 6.0 * rhoc, gb = sres(rb);
    for (int e = 0; e < 12 && ga * gb > 0; ++e) {
        const double rb_next = rb * 2.0;
        double gb_next;
        try {
            gb_next = sres(rb_next);
        } catch (...) {
            break;  // past the EOS density ceiling; keep the last valid bound
        }
        if (!std::isfinite(gb_next)) break;
        rb = rb_next;
        gb = gb_next;
    }
    double rho0;
    if (ga * gb > 0) {
        rho0 = (std::abs(ga) <= std::abs(gb)) ? ra : rb;  // clamp to the closer bound
    } else {
        boost::math::uintmax_t it = 80;
        auto [l, r] = boost::math::tools::toms748_solve(sres, ra, rb, ga, gb, boost::math::tools::eps_tolerance<double>(40), it);
        rho0 = 0.5 * (l + r);
    }
    R.anchorT = Ta;
    R.anchorRho = rho0;
    continue_HS(AS, Ta, rho0, h_t, s_t, R);
    return R;
}

// ===========================================================================
// DEPARTURE (ideal-gas -> real) homotopy.  Scale the residual Helmholtz energy
// by lambda: a(tau,delta;lambda) = a0(tau,delta) + lambda*ar(tau,delta).  At
// lambda=0 the fluid is an IDEAL GAS -- it has NO two-phase dome at all, and the
// (h,s)->(T,rho) inversion is two trivial monotone 1D rootfinds.  As lambda->1
// the dome grows from nothing while we track the stable branch (p_rho>0 guard, so
// det J != 0).  Because the anchor problem vanishes (ideal gas) and the early path
// is dome-free, this structurally avoids the saturation-anchor dome-crossing that
// the (h,s)-line homotopy suffers for supercritical targets.  Self-contained: no
// superancillary needed.
// ===========================================================================
HSResult solve_HS_departure(CoolProp::HelmholtzEOSMixtureBackend& AS, double h_t, double s_t) {
    using namespace CoolProp;
    HSResult R;
    const double Rg = AS.gas_constant();
    const auto& red = AS.get_reducing_state();
    const double Tr = red.T, rhor = red.rhomolar;
    const auto& x = AS.get_mole_fractions_ref();
    const double Tmin = AS.Tmin(), Tmax = AS.Tmax();

    // h, s, p_rho and the (T,rho) Jacobian of the lambda-scaled model.
    struct LP
    {
        double h, s, prho, hT, hr, sT, sr;
    };
    auto lprops = [&](double T, double rho, double lam) -> LP {
        const double tau = Tr / T, delta = rho / rhor;
        ++R.nevals;
        const double a0 = AS.calc_alpha0_deriv_nocache(0, 0, x, tau, delta, Tr, rhor);
        const double a0t = AS.calc_alpha0_deriv_nocache(1, 0, x, tau, delta, Tr, rhor);
        const double a0tt = AS.calc_alpha0_deriv_nocache(2, 0, x, tau, delta, Tr, rhor);
        const double ar = AS.calc_alphar_deriv_nocache(0, 0, x, tau, delta);
        const double art = AS.calc_alphar_deriv_nocache(1, 0, x, tau, delta);
        const double ard = AS.calc_alphar_deriv_nocache(0, 1, x, tau, delta);
        const double artt = AS.calc_alphar_deriv_nocache(2, 0, x, tau, delta);
        const double ardd = AS.calc_alphar_deriv_nocache(0, 2, x, tau, delta);
        const double artd = AS.calc_alphar_deriv_nocache(1, 1, x, tau, delta);
        const double at = a0t + lam * art, att = a0tt + lam * artt;
        const double H = 1.0 + tau * at + delta * lam * ard;
        LP L;
        L.h = Rg * T * H;
        L.s = Rg * (tau * at - (a0 + lam * ar));
        L.prho = Rg * T * (1.0 + 2.0 * delta * lam * ard + delta * delta * lam * ardd);
        const double dHdtau = at + tau * att + delta * lam * artd;
        const double dHddelta = lam * ard + tau * lam * artd + delta * lam * ardd;
        L.hT = Rg * (H - tau * dHdtau);
        L.hr = Rg * T * dHddelta / rhor;
        L.sT = -Rg * tau * tau * att / T;                                 // = cv/T
        L.sr = Rg * (tau * lam * artd - 1.0 / delta - lam * ard) / rhor;  // ideal -1/delta term is lambda-free
        return L;
    };

    // --- Anchor at lambda = 0 (ideal gas): two monotone 1D inversions ---
    // h_ig(T) increases with T (cp_ig > 0); s_ig decreases with rho (-R ln rho).
    double T0 = 0, rho0 = 0;
    {
        auto hig = [&](double T) { return lprops(T, rhor, 0.0).h - h_t; };  // rho cancels in h_ig
        double Ta = 0.3 * Tmin, Tb = 3.0 * Tmax, fa = hig(Ta), fb = hig(Tb);
        if (fa * fb > 0) return R;  // h_t unreachable as ideal gas (should not happen)
        boost::math::uintmax_t it = 80;
        auto [l, r] = boost::math::tools::toms748_solve(hig, Ta, Tb, fa, fb, boost::math::tools::eps_tolerance<double>(40), it);
        T0 = 0.5 * (l + r);
        auto sig = [&](double rho) { return lprops(T0, rho, 0.0).s - s_t; };
        double ra = 1e-8 * rhor, rb = 50.0 * rhor, ga = sig(ra), gb = sig(rb);
        if (ga * gb > 0) {
            rho0 = (std::abs(ga) <= std::abs(gb)) ? ra : rb;
        } else {
            boost::math::uintmax_t it2 = 80;
            auto [l2, r2] = boost::math::tools::toms748_solve(sig, ra, rb, ga, gb, boost::math::tools::eps_tolerance<double>(40), it2);
            rho0 = 0.5 * (l2 + r2);
        }
    }
    R.anchorT = T0;
    R.anchorRho = rho0;

    const double hscale = Rg * Tr, sscale = Rg;
    for (int N = 1; N <= 256; N *= 2) {
        double T = T0, w = std::log(rho0);
        bool failed = false;
        for (int k = 1; k <= N && !failed; ++k) {
            const double lam = static_cast<double>(k) / N;
            bool conv = false;
            double best_norm = 1e300, best_T = T, best_w = w;
            for (int iter = 0; iter < 40; ++iter) {
                const double rho = std::exp(w);
                LP L = lprops(T, rho, lam);
                const double rh = L.h - h_t, rs = L.s - s_t;
                const double norm = std::abs(rh) / hscale + std::abs(rs) / sscale;
                if (norm < best_norm) {
                    best_norm = norm;
                    best_T = T;
                    best_w = w;
                }
                if (norm < 1e-11) {
                    conv = true;
                    break;
                }
                const double a11 = L.hT, a12 = rho * L.hr, a21 = L.sT, a22 = rho * L.sr;
                const double det = a11 * a22 - a12 * a21;
                if (!std::isfinite(det) || std::abs(det) < 1e-300) break;
                double dT = -(a22 * rh - a12 * rs) / det;
                double dw = -(-a21 * rh + a11 * rs) / det;
                double f = 1.0;
                while ((T + f * dT <= 0 || T + f * dT > 3.0 * Tmax || std::abs(f * dw) > 2.0) && f > 1e-6)
                    f *= 0.5;
                for (int g = 0; g < 8 && f > 1e-3; ++g) {  // stay on the mechanically stable branch
                    if (lprops(T + f * dT, std::exp(w + f * dw), lam).prho > 0) break;
                    f *= 0.5;
                }
                T += f * dT;
                w += f * dw;
            }
            if (!conv && best_norm < 1e-8) {
                T = best_T;
                w = best_w;
                conv = true;
            }
            if (!conv) failed = true;
        }
        if (!failed) {
            R.ok = true;
            R.T = T;
            R.rho = std::exp(w);
            R.nsteps = N;
            return R;
        }
    }
    return R;
}

// PRODUCTION cascade: try each method in turn and ACCEPT the first result that is
// both (a) a faithful (h,s) reproduction and (b) mechanically stable (dp/drho>0).
// The det J = -(1/rho)[(dp/dT)_rho^2/rho^2 + cv*(dp/drho)_T/T] argument makes the
// STABLE (h,s)->(T,rho) root unique, so a stable+matching result is THE physical
// root -- no knowledge of the true (T,rho) is needed (unlike the A/B scorer).  A
// method that "converges" by crossing the dome lands on a metastable/unstable
// branch and is rejected, falling through to the next leg.  Order is fastest-
// first: saturation-anchor line, then supercritical isentrope, then ideal-gas
// departure as the universal-but-slow backstop.
HSResult solve_HS_cascade(CoolProp::HelmholtzEOSMixtureBackend& AS, CoolProp::superancillary::SuperAncillary<std::vector<double>>* sa, double h_t,
                          double s_t) {
    using namespace CoolProp;
    const double Rg = AS.gas_constant(), Tc = AS.T_critical();
    const double htol = 1e-6 * Rg * Tc, stol = 1e-6 * Rg;
    const double Tlo = AS.Tmin() * (1 - 1e-6), Thi = AS.Tmax() * (1 + 1e-6);
    auto accept = [&](const HSResult& r) -> bool {
        if (!r.ok || !std::isfinite(r.T) || !std::isfinite(r.rho) || r.rho <= 0) return false;
        // Reject solutions outside the EOS validity range (the small relative slack
        // admits roots that land exactly on Tmin/Tmax): below Tmin the single-phase
        // extrapolation grows a SPURIOUS (h,s)-matching root (e.g. dense cryogenic
        // hydrogen pairs solve to a fake T ~ 7 K state).
        if (r.T < Tlo || r.T > Thi) return false;
        PhaseGuard g(AS);
        AS.update_DmolarT_direct(r.rho, r.T);
        if (std::abs(AS.hmolar() - h_t) > htol || std::abs(AS.smolar() - s_t) > stol) return false;
        // Full intrinsic stability: mechanical (dp/drho|_T > 0) AND thermal (cv > 0).
        // The (h,s) -> (T,rho) map is injective only over the FULLY stable region; a
        // mechanically-stable-but-thermally-unstable (cv < 0) extrapolated state can
        // share the target (h,s) (confirmed for dense supercritical hydrogen: a
        // T~16 K cv=-38 root duplicates a T~101 K cv=+19 state).  Requiring cv > 0
        // rejects that artifact so the cascade falls through to the unique physical
        // root.  cv = -R*tau^2*(a0_tautau + ar_tautau); query via cvmolar().
        double cv = -1;
        try {
            cv = AS.cvmolar();
        } catch (...) {
            return false;
        }
        if (!std::isfinite(cv) || cv <= 0) return false;
        return AS.first_partial_deriv(iP, iDmolar, iT) > 0;  // mechanically stable
    };
    int evals = 0;  // accumulate work across rejected legs for honest cost accounting
    // Each leg is guarded: a leg whose rootfind/Newton wanders out of the EOS domain
    // can throw (e.g. T driven just below Tmin); on a throw we simply fall through to
    // the next leg rather than failing the whole solve.  accept() is evaluated inside
    // the same guard because its (h,s) re-check also touches the EOS.
    auto try_leg = [&](HSResult (*fn)(CoolProp::HelmholtzEOSMixtureBackend&, double, double), int methodnum, HSResult& out) -> bool {
        try {
            out = fn(AS, h_t, s_t);
            evals += out.nevals;
            if (accept(out)) {
                out.method = methodnum;
                out.nevals = evals;
                return true;
            }
        } catch (...) {
        }
        return false;
    };
    HSResult r;
    if (sa != nullptr) {
        try {
            r = solve_HS_continuation(AS, *sa, h_t, s_t);
            evals += r.nevals;
            if (accept(r)) {
                r.method = 1;
                r.nevals = evals;
                return r;
            }
        } catch (...) {
        }
    }
    if (try_leg(&solve_HS_isentrope, 2, r)) return r;
    if (try_leg(&solve_HS_departure, 3, r)) return r;
    r.ok = false;
    r.nevals = evals;
    return r;
}

}  // namespace

// ---------------------------------------------------------------------------
// Baseline: round-trip every single-phase (p, T) grid point through
// HmolarSmolar_INPUTS and record where the current HS_flash fails to recover
// the originating state.
// ---------------------------------------------------------------------------
TEST_CASE("HS single-phase grid round-trip (baseline)", "[HS][HS_grid][.]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "Hydrogen", "Methane", "Oxygen", "R134a", "R1234yf", "R32", "R245fa",
                                       "MM", "MDM", "D4", "D6", "n-Pentane", "n-Decane");
    CAPTURE(fluid);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit"), pc = CoolProp::Props1SI(fluid, "pcrit");
    double plo = pmax * 1e-7;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-7);
    } catch (...) {
        plo = pmax * 1e-7;
    }

    const std::size_t NT = env_or("HS_GRID_NT", 80), NP = env_or("HS_GRID_NP", 80);

    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);  // always single-phase
            } catch (...) {
                continue;  // outside the EOS domain (e.g. solid) -> skip
            }
            const double h = ref->hmolar(), s = ref->smolar();
            const double rho = ref->rhomolar();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) continue;
            // Region label for the failure map.
            const char* region = (T > Tc) ? (p > pc ? "supercrit" : "supercrit_gas") : (rho > ref->rhomolar_critical() ? "liquid" : "vapor");
            CAPTURE(T, p, rho, region);
            try {
                rt->update(CoolProp::HmolarSmolar_INPUTS, h, s);
            } catch (const std::exception& e) {
                FAIL_CHECK("HS update threw: " << e.what());
                continue;
            }
            // Single-phase (h,s) is a bijection to (T,rho): both must recover.
            CHECK(rt->T() == Catch::Approx(T).epsilon(1e-5));
            CHECK(rt->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
        }
    }
}

// ---------------------------------------------------------------------------
// Prototype A/B: run the (T, ln rho) continuation solver over the same
// single-phase grid and report success vs the baseline, plus step/eval counts.
//     CatchTestRunner "[HS_cont]"
// ---------------------------------------------------------------------------
TEST_CASE("HS single-phase continuation prototype A/B", "[HS][HS_cont][.]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "Hydrogen", "Methane", "Oxygen", "R134a", "R1234yf", "R32", "R245fa",
                                       "MM", "MDM", "D4", "D6", "n-Pentane", "n-Decane");
    CAPTURE(fluid);

    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* ref = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* wrk = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    REQUIRE(sa);

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    double plo = pmax * 1e-7;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-7);
    } catch (...) {
        plo = pmax * 1e-7;
    }
    const std::size_t NT = env_or("HS_GRID_NT", 80), NP = env_or("HS_GRID_NP", 80);

    std::size_t total = 0, cont_ok = 0, cont_fail = 0;
    std::size_t step_hist[8] = {0};  // N = 1,2,4,...,128
    std::size_t eval_sum = 0;

    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);
            } catch (...) {
                continue;
            }
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) continue;
            ++total;
            CAPTURE(T, p, rho);
            HSResult r = solve_HS_continuation(*wrk, *sa, h, s);
            if (r.ok && std::abs(r.T - T) / T < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5) {
                ++cont_ok;
                eval_sum += r.nevals;
                int b = 0, n = r.nsteps;
                while (n > 1) {
                    n /= 2;
                    ++b;
                }
                if (b < 8) ++step_hist[b];
            } else {
                ++cont_fail;
                CAPTURE(r.ok, r.T, r.rho, r.nsteps, r.anchorT, r.anchorRho);
                CHECK(r.ok);  // record the point that the prototype missed
                if (r.ok) CHECK(std::abs(r.T - T) / T < 1e-5);
            }
        }
    }
    std::printf("\n[HS_cont] %s: %zu single-phase points; continuation OK=%zu FAIL=%zu; mean EOS evals=%.1f\n", fluid.c_str(), total, cont_ok,
                cont_fail, cont_ok ? static_cast<double>(eval_sum) / cont_ok : 0.0);
    std::printf("[HS_cont] %s: continuation-step histogram N=1:%zu 2:%zu 4:%zu 8:%zu 16:%zu 32:%zu 64:%zu 128:%zu\n", fluid.c_str(), step_hist[0],
                step_hist[1], step_hist[2], step_hist[3], step_hist[4], step_hist[5], step_hist[6], step_hist[7]);
}

// ---------------------------------------------------------------------------
// Speed: continuation solver vs the legacy HmolarSmolar_INPUTS flash.
//     CatchTestRunner "[HS_bench]"
// ---------------------------------------------------------------------------
TEST_CASE("HS single-phase: continuation vs legacy speed", "[HS][HS_bench][.]") {
    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    auto leg = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    struct Case
    {
        const char* name;
        double p, T;
    };
    for (auto c : {Case{"compressed liquid (50 MPa, 320 K)", 50e6, 320.0}, Case{"superheated vapor (0.1 MPa, 600 K)", 1e5, 600.0},
                   Case{"supercritical (30 MPa, 700 K)", 30e6, 700.0}}) {
        ref->update(CoolProp::PT_INPUTS, c.p, c.T);
        const double h = ref->hmolar(), s = ref->smolar();
        BENCHMARK(std::string("continuation: ") + c.name) {
            return solve_HS_continuation(*wrk, *sa, h, s).T;
        };
        BENCHMARK(std::string("legacy HS_flash: ") + c.name) {
            leg->update(CoolProp::HmolarSmolar_INPUTS, h, s);
            return leg->T();
        };
    }
}

// ---------------------------------------------------------------------------
// Diagnostic: why did get_all_intersections('S') miss anchor roots that the EOS
// saturated-entropy scan finds (retrograde MDM)?  Compare the two root sets and
// dump the reference-frame shift.   CatchTestRunner "[HS_dbg]"
// ---------------------------------------------------------------------------
TEST_CASE("HS anchor diagnostic: S-superancillary vs EOS scan (MDM)", "[HS_dbg]") {
    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "MDM"));
    auto* AS = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    auto sa = AS->get_superanc();
    REQUIRE(sa);
    AS->specify_phase(CoolProp::iphase_gas);
    AS->ensure_caloric_superancillaries();
    const double Tmin = sa->get_Tmin(), Tcap = sa->get_Tcrit_num() - 1e-3, Rgas = AS->gas_constant();

    // A representative single-phase MDM state (compressed liquid).
    {
        auto pt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "MDM"));
        pt->update(CoolProp::PT_INPUTS, 1e6, 400.0);
        const double s_t = pt->smolar(), h_t = pt->hmolar();
        std::printf("\nMDM target: T=400 p=1e6  h=%.6g  s=%.6g  has_S=%d\n", h_t, s_t, (int)sa->has_variable('S'));

        // frame shift
        double s_cache = s_t;
        const auto stamp = sa->get_caloric_alpha0_stamp();
        if (stamp.has_value()) {
            const auto [a1, a2] = CoolProp::FlashRoutines::alpha0_offset_total(*AS);
            s_cache = s_t - Rgas * (stamp->first - a1);
            std::printf("  stamp.a1=%.6g caller.a1=%.6g  shift=%.6g  s_cache=%.6g\n", stamp->first, a1, -Rgas * (stamp->first - a1), s_cache);
        }

        std::printf("  get_all_intersections('S', s_cache) roots:\n");
        for (auto& pr : sa->get_all_intersections('S', s_cache, 48, 100, 1e-12))
            std::printf("    T=%.4f\n", pr.first);

        std::printf("  EOS-scan sign changes of sat_s(T)-s_t (caller frame), both branches:\n");
        for (short Q = 0; Q <= 1; ++Q) {
            double prevT = Tmin + 1e-3;
            AS->update_DmolarT_direct(sa->eval_sat(prevT, 'D', Q), prevT);
            double prevg = AS->smolar() - s_t;
            for (int i = 1; i <= 60; ++i) {
                double Ti = (Tmin + 1e-3) + (Tcap - Tmin - 1e-3) * i / 60.0;
                AS->update_DmolarT_direct(sa->eval_sat(Ti, 'D', Q), Ti);
                double gi = AS->smolar() - s_t;
                if (std::isfinite(prevg) && std::isfinite(gi) && prevg * gi <= 0)
                    std::printf("    Q=%d bracket [%.2f, %.2f]  (EOS sat_s crosses s_t)\n", Q, prevT, Ti);
                prevT = Ti;
                prevg = gi;
            }
        }
        // Also: at one mid-T, compare EOS sat_s vs SA eval_sat('S') to check the frame.
        double Tm = 0.6 * sa->get_Tcrit_num();
        AS->update_DmolarT_direct(sa->eval_sat(Tm, 'D', 0), Tm);
        std::printf("  frame check @T=%.2f: EOS sL=%.6g  SA eval_sat('S',0)=%.6g  (SA should ~= EOS sL + (s_cache-s_t))\n", Tm, AS->smolar(),
                    sa->eval_sat(Tm, 'S', 0));
    }

    // Low-pressure superheated vapor (s_t likely above the saturated-vapor range).
    {
        auto pt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "MDM"));
        pt->update(CoolProp::PT_INPUTS, 200.0, 450.0);
        const double s_t = pt->smolar(), h_t = pt->hmolar();
        const double sV_lo = [&] {
            AS->update_DmolarT_direct(sa->eval_sat(Tmin + 1e-3, 'D', 1), Tmin + 1e-3);
            return AS->smolar();
        }();
        const double sV_hi = [&] {
            AS->update_DmolarT_direct(sa->eval_sat(Tcap, 'D', 1), Tcap);
            return AS->smolar();
        }();
        std::printf("\nMDM vapor: T=450 p=200  h=%.6g  s=%.6g | sat-vapor s range [%.6g, %.6g]\n", h_t, s_t, sV_hi, sV_lo);
        auto roots = sa->get_all_intersections('S', s_t, 48, 100, 1e-12);
        std::printf("  get_all_intersections('S') -> %zu roots%s\n", roots.size(), roots.empty() ? "  (none -> needs ideal-gas anchor)" : "");
        for (auto& pr : roots)
            std::printf("    T=%.4f\n", pr.first);
    }
    AS->unspecify_phase();
}

// ---------------------------------------------------------------------------
// Departure (ideal-gas -> real) homotopy over the single-phase grid.
//     CatchTestRunner "[HS_dep]"
// ---------------------------------------------------------------------------
TEST_CASE("HS single-phase departure-homotopy A/B", "[HS][HS_dep][.]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "Hydrogen", "Methane", "Oxygen", "R134a", "R1234yf", "R32", "R245fa",
                                       "MM", "MDM", "D4", "D6", "n-Pentane", "n-Decane");
    CAPTURE(fluid);
    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* ref = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* wrk = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(wrksp.get());

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    double plo = pmax * 1e-7;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-7);
    } catch (...) {
        plo = pmax * 1e-7;
    }
    auto sa = wrk->get_superanc();
    auto w2sp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* w2 = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(w2sp.get());
    auto w3sp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* w3 = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(w3sp.get());
    const double Tc = ref->T_critical();
    const std::size_t NT = env_or("HS_GRID_NT", 80), NP = env_or("HS_GRID_NP", 80);
    std::size_t total = 0, sat_ok = 0, dep_ok = 0, isen_ok = 0, casc_ok = 0, casc_fail = 0;

    auto good = [](const HSResult& r, double T, double rho) { return r.ok && std::abs(r.T - T) / T < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5; };
    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);
            } catch (...) {
                continue;
            }
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) continue;
            ++total;
            const bool sgood = sa && good(solve_HS_continuation(*wrk, *sa, h, s), T, rho);
            const bool dgood = good(solve_HS_departure(*w2, h, s), T, rho);
            const bool igood = good(solve_HS_isentrope(*w3, h, s), T, rho);
            sat_ok += sgood;
            dep_ok += dgood;
            isen_ok += igood;
            const bool cascade = sgood || igood || dgood;
            casc_ok += cascade;
            if (!cascade) {
                ++casc_fail;
                const bool supercrit = T > Tc;
                CAPTURE(T, p, rho, supercrit);
                CHECK(cascade);  // record points NO method solved
            }
        }
    }
    std::printf("[HS_dep] %s: %zu pts; sat=%zu dep=%zu isen=%zu CASCADE_OK=%zu CASCADE_FAIL=%zu\n", fluid.c_str(), total, sat_ok, dep_ok, isen_ok,
                casc_ok, casc_fail);
}

// ---------------------------------------------------------------------------
// Production cascade: selection uses ONLY (h,s) + stability (no truth), scoring
// uses the true (T,rho).  This is the real end-to-end test of the HS solver.
//     CatchTestRunner "[HS_casc]"
// ---------------------------------------------------------------------------
TEST_CASE("HS single-phase production cascade", "[HS][HS_casc][.]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "Hydrogen", "Methane", "Oxygen", "R134a", "R1234yf", "R32", "R245fa",
                                       "MM", "MDM", "D4", "D6", "n-Pentane", "n-Decane");
    CAPTURE(fluid);
    auto refsp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* ref = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto* wrk = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    double plo = pmax * 1e-7;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-7);
    } catch (...) {
        plo = pmax * 1e-7;
    }
    const std::size_t NT = env_or("HS_GRID_NT", 80), NP = env_or("HS_GRID_NP", 80);
    std::size_t total = 0, ok = 0, fail = 0, eval_sum = 0;

    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);
            } catch (...) {
                continue;
            }
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) continue;
            ++total;
            HSResult r = solve_HS_cascade(*wrk, sa.get(), h, s);
            const bool correct = r.ok && std::abs(r.T - T) / T < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5;
            if (correct) {
                ++ok;
                eval_sum += r.nevals;
            } else {
                ++fail;
                CAPTURE(T, p, rho, r.ok, r.T, r.rho);
                CHECK(correct);
            }
        }
    }
    std::printf("[HS_casc] %s: %zu pts; OK=%zu FAIL=%zu; mean evals=%.1f\n", fluid.c_str(), total, ok, fail,
                ok ? static_cast<double>(eval_sum) / ok : 0.0);
}

// ---------------------------------------------------------------------------
// AGGRESSIVE stress + no-regression sweep over EVERY pure fluid in the library.
// For each single-phase (log p, lin T) grid point: (a) the cascade must recover
// (T,rho) AND (b) it must never FAIL where the CURRENT production HS_flash
// succeeds (the no-regression invariant).  Also probes (c) two-phase (h,s) inputs
// -- the cascade must NOT report a single-phase solution for them -- and (d) a
// near-critical refinement band.  Opt-in / slow; build Release for sane runtime:
//     ./build_rel/CatchTestRunner "[HS_stress]"
// Env: HS_STRESS_NT/NP (grid, default 35), HS_STRESS_NFLUIDS (cap fluid count).
// ---------------------------------------------------------------------------
TEST_CASE("HS cascade stress + no-regression, all pure fluids", "[HS][HS_stress][.]") {
    using namespace CoolProp;
    std::vector<std::string> fluids = strsplit(get_global_param_string("fluids_list"), ',');
    const std::size_t cap = env_or("HS_STRESS_NFLUIDS", fluids.size());
    if (cap < fluids.size()) fluids.resize(cap);
    const std::size_t NT = env_or("HS_STRESS_NT", 35), NP = env_or("HS_STRESS_NP", 35);

    std::size_t g_total = 0, g_ok = 0, g_fail = 0, g_regress = 0, g_improve = 0, g_tp = 0, g_tp_bad = 0, g_tp_stable = 0;
    std::array<std::size_t, 4> g_method = {0, 0, 0, 0};
    std::size_t fluids_done = 0, fluids_skipped = 0;

    auto good = [](const HSResult& r, double T, double rho) { return r.ok && std::abs(r.T - T) / T < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5; };

    for (const std::string& fluid : fluids) {
        std::shared_ptr<AbstractState> refsp, wrksp, prodsp;
        HelmholtzEOSMixtureBackend *ref = nullptr, *wrk = nullptr;
        try {
            refsp.reset(AbstractState::factory("HEOS", fluid));
            wrksp.reset(AbstractState::factory("HEOS", fluid));
            prodsp.reset(AbstractState::factory("HEOS", fluid));
            ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
            wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
        } catch (...) {
            ++fluids_skipped;
            continue;
        }
        if (ref == nullptr || wrk == nullptr || !wrk->is_pure()) {
            ++fluids_skipped;  // the happy-path cascade targets pure fluids (mirror of HSU_D is_pure gate)
            continue;
        }
        double Tmin = 0, Tmax = 0, pmax = 0, Tc = 0, pc = 0;
        try {
            Tmin = ref->Tmin();
            Tmax = ref->Tmax();
            pmax = ref->pmax();
            Tc = ref->T_critical();  // throws for pseudo-pure -> skip (out of scope)
            pc = ref->p_critical();
        } catch (...) {
            ++fluids_skipped;
            continue;
        }
        ++fluids_done;
        std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> sa;
        try {
            sa = wrk->get_superanc();  // pure-only; pseudo-pure -> null, cascade uses legs 2/3
        } catch (...) {
        }
        double plo = pmax * 1e-7;
        try {
            plo = std::max(ref->p_triple(), pmax * 1e-7);
        } catch (...) {
        }

        std::size_t f_total = 0, f_fail = 0;

        // Evaluate whatever single-phase state `ref` currently holds: cascade-solve
        // from its (h,s), score against the true (T,rho), and compare no-regression
        // against the current production HS_flash.
        auto process = [&]() {
            if (ref->phase() == iphase_twophase) return;
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar(), Ttrue = ref->T();
            if (!std::isfinite(h) || !std::isfinite(s) || !std::isfinite(rho)) return;
            // Only count states that are a VALID single-phase flash target: fully
            // intrinsically stable (cv>0 AND dp/drho|_T>0).  The near-critical
            // density band can place a sample inside the spinodal (cv<0) -- an
            // unstable extrapolation that ref->phase() does not flag as two-phase;
            // such an (h,s) has no physical single-phase solution, so feeding it
            // would be a bogus test of the solver (and of "no-regression").
            try {
                const double cvref = ref->cvmolar();
                if (!std::isfinite(cvref) || cvref <= 0) return;
                if (ref->first_partial_deriv(iP, iDmolar, iT) <= 0) return;
            } catch (...) {
                return;
            }
            ++f_total;
            ++g_total;

            HSResult r = solve_HS_cascade(*wrk, sa.get(), h, s);
            const bool cgood = good(r, Ttrue, rho);
            if (cgood) {
                ++g_ok;
                if (r.method >= 1 && r.method <= 3) ++g_method[r.method];
            } else {
                ++g_fail;
                ++f_fail;
                CAPTURE(fluid, Ttrue, rho, r.ok, r.T, r.rho, r.method);
                CHECK(cgood);
            }

            // No-regression: does the CURRENT production HS_flash get this point?
            bool prod_ok = false;
            try {
                prodsp->update(HmolarSmolar_INPUTS, h, s);
                prod_ok = std::abs(prodsp->T() - Ttrue) / Ttrue < 1e-4 && std::abs(prodsp->rhomolar() - rho) / rho < 1e-4;
            } catch (...) {
                prod_ok = false;
            }
            const bool regressed = prod_ok && !cgood;  // cascade must not lose what production had
            if (regressed) {
                ++g_regress;
                CAPTURE(fluid, Ttrue, rho);
                CHECK_FALSE(regressed);
            }
            if (cgood && !prod_ok) ++g_improve;
        };

        // (1) uniform (log p, lin T) grid.
        for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
            for (std::size_t j = 0; j < NP; ++j) {
                try {
                    ref->update(PT_INPUTS, plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1)), T);
                } catch (...) {
                    continue;
                }
                process();
            }
        }
        // (2) near-critical refinement, also in (p,T): a denser (log p, lin T) band
        // bracketing the critical point.  Sampling in (p,T) -- as the ConsistencyPlot
        // suite does -- guarantees every point is a single-phase, intrinsically
        // stable state, so the truth is always a valid HS flash target (density-based
        // sampling can land inside the spinodal, an unstable non-target).
        if (std::isfinite(Tc) && Tc > Tmin && std::isfinite(pc) && pc > 0) {
            const double Tlo = std::max(0.85 * Tc, Tmin + 0.1), Thi = std::min(1.25 * Tc, Tmax);
            const double pclo = std::max(0.3 * pc, plo), pchi = std::min(4.0 * pc, pmax);
            for (double T : linspace(Tlo, Thi, 16)) {
                for (std::size_t j = 0; j < 16; ++j) {
                    try {
                        ref->update(PT_INPUTS, pclo * std::pow(pchi / pclo, static_cast<double>(j) / 15.0), T);
                    } catch (...) {
                        continue;
                    }
                    process();
                }
            }
        }

        // (c) Two-phase inputs: a two-phase (h,s) lies inside the (h,s) lens, where
        // the only single-phase states sharing that (h,s) are METASTABLE extensions
        // (inside the dome, between binodal and spinodal -- still cv>0, dp/drho|_T>0).
        // The single-phase cascade is allowed to return such a root; phase
        // determination is the DISPATCHER's job (must run before this solver).  What
        // would be a genuine BUG is returning a state OUTSIDE the dome (a truly
        // stable single-phase state with a two-phase state's (h,s) -- impossible).
        // So classify any returned root by dome membership and only fail on the
        // outside-dome (stable) case.
        if (std::isfinite(Tc) && Tc > Tmin && sa) {
            for (double T : linspace(Tmin + 1.0, Tc - std::max(1e-3, 1e-4 * Tc), 8)) {
                for (double Q : {0.2, 0.5, 0.8}) {
                    try {
                        ref->update(QT_INPUTS, Q, T);
                    } catch (...) {
                        continue;
                    }
                    const double h = ref->hmolar(), s = ref->smolar();
                    ++g_tp;
                    HSResult r = solve_HS_cascade(*wrk, sa.get(), h, s);
                    if (!r.ok) continue;
                    bool inside_dome = false;
                    if (r.T < Tc) {
                        const double rhoV = sa->eval_sat(r.T, 'D', 1), rhoL = sa->eval_sat(r.T, 'D', 0);
                        inside_dome = (r.rho > rhoV && r.rho < rhoL);
                    }
                    if (inside_dome) {
                        ++g_tp_bad;  // metastable extension: expected, not a defect
                    } else {
                        ++g_tp_stable;  // INVARIANT VIOLATION: a stable single-phase root for a two-phase (h,s)
                        CAPTURE(fluid, T, Q, r.T, r.rho, r.method);
                        CHECK(inside_dome);
                    }
                }
            }
        }

        std::printf("[HS_stress] %-28s pts=%zu fail=%zu\n", fluid.c_str(), f_total, f_fail);
    }

    std::printf("\n[HS_stress] ============ SUMMARY ============\n");
    std::printf("[HS_stress] fluids: %zu tested, %zu skipped\n", fluids_done, fluids_skipped);
    std::printf("[HS_stress] single-phase points: %zu  OK=%zu  FAIL=%zu  (%.5f%%)\n", g_total, g_ok, g_fail, g_total ? 100.0 * g_ok / g_total : 0.0);
    std::printf("[HS_stress] no-regression vs production HS_flash: REGRESSIONS=%zu  improvements=%zu\n", g_regress, g_improve);
    std::printf("[HS_stress] method wins: sat=%zu isentrope=%zu departure=%zu\n", g_method[1], g_method[2], g_method[3]);
    std::printf("[HS_stress] two-phase inputs: %zu probed -> %zu metastable-extension roots (expected), %zu STABLE roots (BUG)\n", g_tp, g_tp_bad,
                g_tp_stable);

    CHECK(g_fail == 0);       // cascade solves every single-phase point
    CHECK(g_regress == 0);    // and never regresses against production
    CHECK(g_tp_stable == 0);  // and never returns a stable single-phase root for a two-phase (h,s)
}

// ---------------------------------------------------------------------------
// Alternate-reference-state robustness: the saturation-anchor leg shifts the
// target entropy into the caloric superancillary's stamped frame (#2773).  Re-run
// the production cascade for water under several reference states to prove the
// frame handling is correct (a wrong shift would desync method 1 silently).
//     CatchTestRunner "[HS_refstate]"
// ---------------------------------------------------------------------------
TEST_CASE("HS cascade under alternate reference states (water)", "[HS][HS_refstate][.]") {
    using namespace CoolProp;
    for (const std::string& refstate : {std::string("DEF"), std::string("NBP"), std::string("IIR"), std::string("ASHRAE")}) {
        try {
            set_reference_stateS("Water", refstate);
        } catch (...) {
            continue;
        }
        auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
        auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
        auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
        auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
        auto sa = wrk->get_superanc();
        const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
        const double plo = std::max(ref->p_triple(), pmax * 1e-7);
        std::size_t total = 0, ok = 0;
        for (double T : linspace(Tmin + 0.1, Tmax, 30)) {
            for (std::size_t j = 0; j < 30; ++j) {
                const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / 29.0);
                try {
                    ref->update(PT_INPUTS, p, T);
                } catch (...) {
                    continue;
                }
                if (ref->phase() == iphase_twophase) continue;
                const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar(), Tt = ref->T();
                if (!std::isfinite(h) || !std::isfinite(s)) continue;
                ++total;
                HSResult r = solve_HS_cascade(*wrk, sa.get(), h, s);
                if (r.ok && std::abs(r.T - Tt) / Tt < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5)
                    ++ok;
                else {
                    CAPTURE(refstate, Tt, rho, r.ok, r.T, r.rho);
                    CHECK(false);
                }
            }
        }
        std::printf("[HS_refstate] Water/%-7s: %zu/%zu\n", refstate.c_str(), ok, total);
    }
    set_reference_stateS("Water", "DEF");  // restore default for other tests
}

// ---------------------------------------------------------------------------
// Water domain map: dump per-point cascade cost (EOS evals + wall time + winning
// leg) over the full single-phase (p,T) domain to CSV for the interactive Python
// plot (wrappers/Python/CoolProp/Plots/hs_water_timing.py).  Build Release for
// meaningful timings:  ./build_rel/CatchTestRunner "[HS_watermap]"
// Env: HS_MAP_NT/NP (default 220), HS_MAP_CSV (output path).
// ---------------------------------------------------------------------------
TEST_CASE("HS water domain timing map", "[HS][HS_watermap][.]") {
    using namespace CoolProp;
    auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    // Legacy production HS_flash, timed on the SAME points for a head-to-head.
    auto legsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    const double plo = std::max(ref->p_triple(), pmax * 1e-7);
    const std::size_t NT = env_or("HS_MAP_NT", 220), NP = env_or("HS_MAP_NP", 220);
    const char* path_env = std::getenv("HS_MAP_CSV");
    const std::string path = (path_env != nullptr) ? path_env : "water_hs_map.csv";
    std::ofstream out(path);
    out << "T_K,p_Pa,rho_moldm3,h_Jmol,s_JmolK,ok,T_solved,rho_solved,evals,method,microseconds,legacy_ok,legacy_us\n";

    // Warm up: the FIRST solve lazily builds the (caloric) superancillaries, which
    // costs milliseconds and would otherwise pollute the first timed cell.
    try {
        ref->update(PT_INPUTS, std::sqrt(plo * pmax), 0.5 * (Tmin + Tmax));
        solve_HS_cascade(*wrk, sa.get(), ref->hmolar(), ref->smolar());
        legsp->update(HmolarSmolar_INPUTS, ref->hmolar(), ref->smolar());
    } catch (...) {
    }

    // The legacy production HS_flash is ~1000x slower (heavy tail; one near-critical
    // point can take ~0.5 s), so the full-grid head-to-head is opt-in: set
    // HS_MAP_LEGACY=1.  Default off keeps the map fast (cascade-only timing).
    const bool do_legacy = std::getenv("HS_MAP_LEGACY") != nullptr;
    std::size_t total = 0, ok = 0, leg_ok = 0;
    double casc_us_sum = 0, leg_us_sum = 0;
    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(PT_INPUTS, p, T);
            } catch (...) {
                continue;
            }
            if (ref->phase() == iphase_twophase) continue;
            const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar(), Tt = ref->T();
            if (!std::isfinite(h) || !std::isfinite(s)) continue;
            ++total;
            // New cascade.
            const auto t0 = std::chrono::steady_clock::now();
            HSResult r = solve_HS_cascade(*wrk, sa.get(), h, s);
            const auto t1 = std::chrono::steady_clock::now();
            const double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
            const bool correct = r.ok && std::abs(r.T - Tt) / Tt < 1e-5 && std::abs(r.rho - rho) / rho < 1e-5;
            ok += correct;
            casc_us_sum += us;
            // Legacy production HS_flash on the same (h,s) (timed incl. any throw).
            bool leg_correct = false;
            double leg_us = 0;
            if (do_legacy) {
                const auto l0 = std::chrono::steady_clock::now();
                try {
                    legsp->update(HmolarSmolar_INPUTS, h, s);
                    leg_correct = std::abs(legsp->T() - Tt) / Tt < 1e-5 && std::abs(legsp->rhomolar() - rho) / rho < 1e-5;
                } catch (...) {
                    leg_correct = false;
                }
                leg_us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - l0).count();
                leg_ok += leg_correct;
                leg_us_sum += leg_us;
            }
            out << Tt << ',' << p << ',' << rho / 1000.0 << ',' << h << ',' << s << ',' << (correct ? 1 : 0) << ',' << r.T << ',' << r.rho << ','
                << r.nevals << ',' << r.method << ',' << us << ',' << (leg_correct ? 1 : 0) << ',' << leg_us << '\n';
        }
    }
    out.close();
    std::printf("[HS_watermap] wrote %zu rows to %s\n", total, path.c_str());
    std::printf("[HS_watermap] cascade: %zu/%zu correct (%zu fail), total %.1f ms (mean %.2f us)\n", ok, total, total - ok, casc_us_sum / 1e3,
                casc_us_sum / total);
    if (do_legacy) {
        std::printf("[HS_watermap] legacy : %zu/%zu correct (%zu fail), total %.1f ms (mean %.2f us)  -> cascade is %.1fx faster\n", leg_ok, total,
                    total - leg_ok, leg_us_sum / 1e3, leg_us_sum / total, leg_us_sum / casc_us_sum);
    }
    CHECK(ok == total);
}

// Diagnostic probe: walk the (h,s) chord from the triple-point saturated-liquid
// anchor to a cold compressed-liquid target and print the (T,rho) trajectory --
// is the obstruction a fold (det J sign flip), the Tmin clamp, or a dome exit?
//     CatchTestRunner "[HS_probe6]"
TEST_CASE("HS probe cold-liquid homotopy trajectory", "[HS][HS_probe6][.]") {
    using namespace CoolProp;
    auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    const double Tmin = ref->Tmin();
    ref->update(PT_INPUTS, 6.943e7, 273.26);
    const double h_t = ref->hmolar(), s_t = ref->smolar();
    const double T0 = Tmin + 1e-3, rho0 = sa->eval_sat(T0, 'D', 0);
    wrk->specify_phase(iphase_gas);
    wrk->update_DmolarT_direct(rho0, T0);
    const double h0 = wrk->hmolar(), s0 = wrk->smolar();
    std::printf("[HS_probe6] anchor T=%.4f rho=%.1f h0=%.2f s0=%.4f | target h=%.2f s=%.4f (Tmin=%.4f)\n", T0, rho0, h0, s0, h_t, s_t, Tmin);
    double T = T0, w = std::log(rho0);
    const int N = 2000;
    double prev_det = 0;
    for (int k = 1; k <= N; ++k) {
        const double lam = static_cast<double>(k) / N;
        const double ht = h0 + lam * (h_t - h0), st = s0 + lam * (s_t - s0);
        for (int iter = 0; iter < 60; ++iter) {
            const double rho = std::exp(w);
            wrk->update_DmolarT_direct(rho, T);
            const double rh = wrk->hmolar() - ht, rs = wrk->smolar() - st;
            const double hT = wrk->first_partial_deriv(iHmolar, iT, iDmolar), hr = wrk->first_partial_deriv(iHmolar, iDmolar, iT);
            const double sT = wrk->first_partial_deriv(iSmolar, iT, iDmolar), sr = wrk->first_partial_deriv(iSmolar, iDmolar, iT);
            const double a11 = hT, a12 = rho * hr, a21 = sT, a22 = rho * sr;
            const double det = a11 * a22 - a12 * a21;
            if (k == 1 && iter == 0) prev_det = det;
            if (det * prev_det < 0) {
                std::printf("[HS_probe6] *** det J sign flip (FOLD) at lam=%.4f T=%.4f rho=%.1f ***\n", lam, T, std::exp(w));
                prev_det = det;
            }
            prev_det = det;
            if (std::abs(rh) / (8.314 * 647) + std::abs(rs) / 8.314 < 1e-11) break;
            double dT = -(a22 * rh - a12 * rs) / det, dw = -(-a21 * rh + a11 * rs) / det;
            double f = 1.0;
            const double Tlo = 0.90 * Tmin;  // relaxed: allow a small sub-Tmin excursion
            while ((T + f * dT < Tlo || std::abs(f * dw) > 2.0) && f > 1e-6)
                f *= 0.5;
            T += f * dT;
            w += f * dw;
        }
        if (k % 200 == 0 || k == 1) std::printf("[HS_probe6] lam=%.3f T=%.4f rho=%.1f p=%.3g\n", lam, T, std::exp(w), wrk->p());
    }
    wrk->unspecify_phase();
    std::printf("[HS_probe6] final T=%.4f rho=%.1f (target 273.26 / %.1f)\n", T, std::exp(w), ref->rhomolar());
}

// Diagnostic probe: why do cold s<0 compressed-liquid points fall to leg 2?
//     CatchTestRunner "[HS_probe5]"
TEST_CASE("HS probe cold compressed-liquid cost", "[HS][HS_probe5][.]") {
    using namespace CoolProp;
    auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    std::printf("[HS_probe5] Tmin=%.4f  sat-liq @Tmin: rho=%.1f\n", ref->Tmin(), sa->eval_sat(ref->Tmin() + 1e-3, 'D', 0));
    for (auto tp : {std::make_pair(273.26, 6.943e7), std::make_pair(273.26, 1.437e8)}) {
        ref->update(PT_INPUTS, tp.second, tp.first);
        const double h = ref->hmolar(), s = ref->smolar(), rho = ref->rhomolar();
        HSResult c = solve_HS_continuation(*wrk, *sa, h, s);
        HSResult i = solve_HS_isentrope(*wrk, h, s);
        std::printf("[HS_probe5] T=%.2f p=%.3g s=%+.4f rho_true=%.1f | leg1 ok=%d nsteps=%d nevals=%d anchorT=%.3f anchorRho=%.1f T=%.3f rho=%.1f"
                    " | leg2 ok=%d nevals=%d\n",
                    tp.first, tp.second, s, rho, (int)c.ok, c.nsteps, c.nevals, c.anchorT, c.anchorRho, c.T, c.rho, (int)i.ok, i.nevals);
    }
}

// Diagnostic probe: why do critical-isentrope supercritical points cost ~4000
// evals?  Walk the s = s_crit isentrope and report each leg's nsteps/nevals.
//     CatchTestRunner "[HS_probe3]"
TEST_CASE("HS probe critical-isentrope cost", "[HS][HS_probe3][.]") {
    using namespace CoolProp;
    auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    const double pc = ref->p_critical();
    // The exact worst-cost points from the water map (s ~ s_crit, supercritical).
    for (auto tp : {std::make_pair(833.07, 10.087), std::make_pair(699.03, 1.970), std::make_pair(722.68, 2.731)}) {
        ref->update(PT_INPUTS, tp.second * pc, tp.first);
        const double h = ref->hmolar(), s = ref->smolar();
        HSResult c = solve_HS_continuation(*wrk, *sa, h, s);  // leg 1 (saturation anchor)
        HSResult i = solve_HS_isentrope(*wrk, h, s);          // leg 2 (isentrope from Tmax)
        std::printf(
          "[HS_probe3] T=%.1f p/pc=%5.2f: leg1 ok=%d nsteps=%d nevals=%d T=%.2f rho=%.1f | leg2 ok=%d nsteps=%d nevals=%d T=%.2f rho=%.1f\n",
          tp.first, tp.second, (int)c.ok, c.nsteps, c.nevals, c.T, c.rho, (int)i.ok, i.nsteps, i.nevals, i.T, i.rho);
    }
}

// Diagnostic probe: is there a coverage GAP in the entropy superancillary near
// s_crit?  For the worst-cost supercritical points, does a saturated state at the
// target entropy exist (EOS + SA scan), and does get_all_intersections find it?
//     CatchTestRunner "[HS_probe4]"
TEST_CASE("HS probe entropy-superancillary gap near s_crit", "[HS][HS_probe4][.]") {
    using namespace CoolProp;
    auto refsp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* ref = dynamic_cast<HelmholtzEOSMixtureBackend*>(refsp.get());
    auto wrksp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
    auto* wrk = dynamic_cast<HelmholtzEOSMixtureBackend*>(wrksp.get());
    auto sa = wrk->get_superanc();
    wrk->ensure_caloric_superancillaries();
    const double Rgas = wrk->gas_constant();
    const double Tc_sa = sa->get_Tcrit_num(), Tmin_sa = sa->get_Tmin();
    std::printf("[HS_probe4] SA: Tmin=%.4f Tcrit=%.6f\n", Tmin_sa, Tc_sa);
    // SA saturated entropy on both branches as T -> Tcrit (look for a gap/truncation).
    std::printf("[HS_probe4] SA sat entropy near critical (Tc-dT):\n");
    for (double dT : {5.0, 1.0, 0.1, 0.01, 0.001}) {
        const double T = Tc_sa - dT;
        std::printf("           dT=%-7.3g sL=%.5f sV=%.5f\n", dT, sa->eval_sat(T, 'S', 0), sa->eval_sat(T, 'S', 1));
    }
    for (auto tp : {std::make_pair(833.07, 10.087), std::make_pair(699.03, 1.970)}) {
        ref->update(PT_INPUTS, tp.second * ref->p_critical(), tp.first);
        const double s_t = ref->smolar();
        // Shift into the SA stamped frame (same as solve_HS_continuation).
        double s_cache = s_t;
        const auto stamp = sa->get_caloric_alpha0_stamp();
        if (stamp.has_value()) {
            const auto [a1, a2] = CoolProp::FlashRoutines::alpha0_offset_total(*wrk);
            s_cache = s_t - Rgas * (stamp->first - a1);
        }
        // Does a saturated VAPOR state at s_cache exist?  Scan the SA vapor branch.
        bool vapor_crosses = false;
        double prev = sa->eval_sat(Tmin_sa + 1e-3, 'S', 1);
        for (int i = 1; i <= 400; ++i) {
            const double T = (Tmin_sa + 1e-3) + (Tc_sa - 1e-4 - (Tmin_sa + 1e-3)) * i / 400.0;
            const double sv = sa->eval_sat(T, 'S', 1);
            if ((prev - s_cache) * (sv - s_cache) <= 0) {
                vapor_crosses = true;
            }
            prev = sv;
        }
        auto roots = sa->get_all_intersections('S', s_cache, 48, 100, 1e-12);
        std::printf("[HS_probe4] T=%.1f: s_t=%.5f s_cache=%.5f | vapor-branch crosses s_cache? %d | get_all_intersections -> %zu roots\n", tp.first,
                    s_t, s_cache, (int)vapor_crosses, roots.size());
        for (auto& r : roots)
            std::printf("              root T=%.4f\n", r.first);
    }
}

// Diagnostic probe: why does accept() reject the (correct) dense near-critical
// Fluorine state?   CatchTestRunner "[HS_probe2]"
TEST_CASE("HS probe Fluorine near-critical reject", "[HS][HS_probe2][.]") {
    using namespace CoolProp;
    auto sp = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Fluorine"));
    auto* AS = dynamic_cast<HelmholtzEOSMixtureBackend*>(sp.get());
    std::printf("[HS_probe2] Tc=%.4f rhoc=%.2f Tmax=%.2f\n", AS->T_critical(), AS->rhomolar_critical(), AS->Tmax());
    for (auto pr : {std::make_pair(137.19370617114313404, 43568.35336238020681776), std::make_pair(137.19370617114313404, 46808.97468685475905659)}) {
        const double T = pr.first, rho = pr.second;
        AS->specify_phase(iphase_gas);
        AS->update_DmolarT_direct(rho, T);
        const double pr_T = AS->first_partial_deriv(iP, iDmolar, iT);
        double cv = std::nan("");
        try {
            cv = AS->cvmolar();
        } catch (const std::exception& e) {
            std::printf("[HS_probe2] cvmolar threw: %s\n", e.what());
        }
        std::printf("[HS_probe2] T=%.6g rho=%.6g (%.3f*rhoc): dp/drho|T=%.5g cv=%.5g p=%.6g\n", T, rho, rho / AS->rhomolar_critical(), pr_T, cv,
                    AS->p());
        AS->unspecify_phase();
    }
}

// Diagnostic probe: are the dense-supercritical-hydrogen "wrong" roots genuine
// co-(h,s) stable states or pathologies?   CatchTestRunner "[HS_probe]"
TEST_CASE("HS probe hydrogen dual root", "[HS][HS_probe][.]") {
    auto sp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Hydrogen"));
    auto* AS = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(sp.get());
    AS->update(CoolProp::PT_INPUTS, 659950683.68263840675354004, 101.41903797468354753);
    const double h_t = AS->hmolar(), s_t = AS->smolar();
    std::printf("[HS_probe] target h=%.10g s=%.10g (T=101.419 rho=69734.1)\n", h_t, s_t);
    auto report = [&](double T, double rho, const char* tag) -> double {
        AS->specify_phase(CoolProp::iphase_gas);
        AS->update_DmolarT_direct(rho, T);
        const double h = AS->hmolar(), s = AS->smolar();
        const double pr = AS->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT);
        const double prs = AS->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iSmolar);
        double cv = std::nan("");
        try {
            cv = AS->cvmolar();
        } catch (...) {
        }
        AS->unspecify_phase();
        std::printf("[HS_probe] %-6s T=%.6g rho=%.6g: dh=%.3e ds=%.3e  dp/drho|T=%.4g  dp/drho|s=%.4g  cv=%.4g\n", tag, T, rho, h - h_t, s - s_t, pr,
                    prs, cv);
        return cv;
    };
    const double cv_true = report(101.41903797468354753, 69734.11745331903512124, "true");
    const double cv_alt = report(16.00826076445090251, 72872.58412607655918691, "alt");
    // The crux of HS non-uniqueness: both states reproduce the target (h,s) and are
    // mechanically stable, but the spurious one is THERMALLY unstable (cv < 0).  The
    // cv > 0 gate in solve_HS_cascade::accept relies on exactly this separation.
    CHECK(cv_true > 0);
    CHECK(cv_alt < 0);
}

// ---------------------------------------------------------------------------
// PRODUCTION two-phase HS round-trip: exercises the new superancillary
// two-phase pre-screen + EOS-exact HS_flash_twophase path in FlashRoutines.
//     CatchTestRunner "[HS_prod2ph]"
// ---------------------------------------------------------------------------
TEST_CASE("HS production two-phase round-trip", "[HS][HS_prod2ph]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "Nitrogen", "R134a", "MM", "n-Pentane", "Methane");
    CAPTURE(fluid);
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tc = ref->T_critical(), Tmin = ref->Ttriple();
    std::size_t total = 0, ok = 0;
    for (double T : linspace(Tmin + 1.0, Tc - std::max(0.5, 1e-3 * Tc), 25)) {
        for (double Q : {0.05, 0.25, 0.5, 0.75, 0.95}) {
            try {
                ref->update(CoolProp::QT_INPUTS, Q, T);
            } catch (...) {
                continue;
            }
            const double h = ref->hmolar(), s = ref->smolar();
            if (!std::isfinite(h) || !std::isfinite(s)) continue;
            ++total;
            CAPTURE(T, Q, h, s);
            try {
                wrk->update(CoolProp::HmolarSmolar_INPUTS, h, s);
            } catch (const std::exception& e) {
                FAIL_CHECK("two-phase HS threw: " << e.what());
                continue;
            }
            const bool good = std::abs(wrk->T() - T) / T < 1e-5 && std::abs(wrk->Q() - Q) < 1e-4;
            if (good)
                ++ok;
            else {
                CAPTURE(wrk->T(), wrk->Q());
                CHECK(good);
            }
        }
    }
    std::printf("[HS_prod2ph] %s: %zu/%zu two-phase round-trips\n", fluid.c_str(), ok, total);
}

#endif  // ENABLE_CATCH
