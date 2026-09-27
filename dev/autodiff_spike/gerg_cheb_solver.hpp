// Chebyshev all-roots density solver for the GERG-2008 multi-fluid model (spike).
// See README.md, Experiments 6-7.  Header-only so the benchmark (gerg_cheb.cpp) and the
// validation harness (gerg_validate.cpp) share one implementation.
//
//   G(delta) = delta Z(delta) - t,   t = p / (rho_r(x) R T),   Z = 1 + sum_g W_g(T, x) chi_g(delta)
//
// Tiers: component set -> build(); (T, x) -> assemble(); p -> roots().
// Thread safety: a built Solver is read-only; State is per thread; counters are thread_local.
#pragma once
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Configuration.h"
#include <stdexcept>
#include "cheb2bern.hpp"

#ifndef NQ_DEG
#    define NQ_DEG 16
#endif
#if NQ_DEG == 16
#    define GERGCHEB_C2B CHEB2BERN_17
#elif NQ_DEG == 20
#    define GERGCHEB_C2B CHEB2BERN_21
#elif NQ_DEG == 24
#    define GERGCHEB_C2B CHEB2BERN_25
#else
#    error "NQ_DEG must be 16, 20 or 24 (Chebyshev-to-Bernstein matrices are generated for those)"
#endif

namespace gergcheb {
constexpr double PI = 3.14159265358979323846;
constexpr int NQ = NQ_DEG, NG = NQ + 1;
using VecQ = std::array<double, NQ + 1>;
using VecG = std::array<double, NG + 1>;

// instrumentation (per thread)
inline thread_local long g_na_active = 0, g_dropped = 0, g_na_deg8 = 0;
inline thread_local long g_nodes = 0, g_calls = 0, g_pieces_open = 0, g_ref_its = 0, g_ref_calls = 0, g_uncertain = 0, g_pol_its = 0, g_pol_calls = 0,
                         g_pol_uncert = 0;

// ------------------------------------------------------------------ terms
struct Term
{
    int i, j;  // component (j = -1 for pure-fluid terms)
    double n, d, t;
    bool has_cl, has_om, has_e1, has_e2, has_b1, has_b2;
    double c, l, om, m, e1, eps1, e2, eps2, b1, g1, b2, g2;
    int di = -1, li = -1;  // d and l as small non-negative integers (powers from a table), else -1
    // chi = delta d/d(delta) [delta^d e^{u(delta)}] = delta^d e^u (d + delta u')
    double chi(double D) const {
        double f, df;
        chi_d(D, f, df);
        return f;
    }
    // chi and chi' = delta^(d-1) e^u [a^2 + delta u' + delta^2 u''],  a = d + delta u'
    void chi_d(double D, double& f, double& df) const {
        double u = 0, du = 0, d2u = 0;
        if (has_cl) {
            const double dl = std::pow(D, l);
            u -= c * dl;
            du -= c * l * dl / D;
            d2u -= c * l * (l - 1) * dl / (D * D);
        }
        if (has_e1) {
            u -= e1 * (D - eps1);
            du -= e1;
        }
        if (has_e2) {
            u -= e2 * (D - eps2) * (D - eps2);
            du -= 2 * e2 * (D - eps2);
            d2u -= 2 * e2;
        }
        const double base = std::pow(D, d) * std::exp(u), a = d + D * du;
        f = base * a;
        df = base / D * (a * a + D * du + D * D * d2u);
    }
    // chi and chi' with delta^k taken from pw[k] = delta^k (k <= 24) when the exponents are integers
    void chi_d_pw(double D, const double* pw, double& f, double& df) const {
        if (di < 0 || (has_cl && li < 0)) {
            chi_d(D, f, df);
            return;
        }
        double u = 0, du = 0, d2u = 0;
        const double iD = 1 / D;
        if (has_cl) {
            const double dl = pw[li];
            u -= c * dl;
            du -= c * l * dl * iD;
            d2u -= c * l * (l - 1) * dl * iD * iD;
        }
        if (has_e1) {
            u -= e1 * (D - eps1);
            du -= e1;
        }
        if (has_e2) {
            u -= e2 * (D - eps2) * (D - eps2);
            du -= 2 * e2 * (D - eps2);
            d2u -= 2 * e2;
        }
        const double base = pw[di] * std::exp(u), a = d + D * du;
        f = base * a;
        df = base * iD * (a * a + D * du + D * D * d2u);
    }
    // n tau^t exp(u_tau)
    double kappa(double tau, double ltau) const {
        double u = t * ltau;
        if (has_om) u -= om * std::exp(m * ltau);
        if (has_b1) u -= b1 * (tau - g1);
        if (has_b2) u -= b2 * (tau - g2) * (tau - g2);
        return n * std::exp(u);
    }
    bool same_delta(const Term& a) const {
        return a.d == d && a.has_cl == has_cl && (!has_cl || (a.c == c && a.l == l)) && a.has_e1 == has_e1
               && (!has_e1 || (a.e1 == e1 && a.eps1 == eps1)) && a.has_e2 == has_e2 && (!has_e2 || (a.e2 == e2 && a.eps2 == eps2));
    }
};
inline bool finite_(double v) {
    return std::isfinite(v);
}
inline void collect(const CoolProp::ResidualHelmholtzGeneralizedExponential& g, int i, int j, std::vector<Term>& out) {
    for (const auto& el : g.elements) {
        Term T{};
        T.i = i;
        T.j = j;
        T.n = el.n;
        T.d = el.d;
        T.t = el.t;
        T.has_cl = g.delta_li_in_u && finite_(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON;
        T.c = el.c;
        T.l = el.l_double;
        T.has_om = g.tau_mi_in_u && std::abs(el.m_double) > 0;
        T.om = el.omega;
        T.m = el.m_double;
        T.has_e1 = g.eta1_in_u && finite_(el.eta1);
        T.e1 = el.eta1;
        T.eps1 = el.epsilon1;
        T.has_e2 = g.eta2_in_u && finite_(el.eta2);
        T.e2 = el.eta2;
        T.eps2 = el.epsilon2;
        T.has_b1 = g.beta1_in_u && finite_(el.beta1);
        T.b1 = el.beta1;
        T.g1 = el.gamma1;
        T.has_b2 = g.beta2_in_u && finite_(el.beta2);
        T.b2 = el.beta2;
        T.g2 = el.gamma2;
        if (T.d >= 0 && T.d <= 24 && T.d == std::floor(T.d)) T.di = static_cast<int>(T.d);
        if (T.has_cl && T.l >= 0 && T.l <= 24 && T.l == std::floor(T.l)) T.li = static_cast<int>(T.l);
        out.push_back(T);
    }
}

// ------------------------------------------------------------------ Chebyshev utilities
template <int M, class F>
std::array<double, M + 1> chebfit_n(F f, double lo, double hi) {
    std::array<double, M + 1> fv{}, c{};
    for (int j = 0; j <= M; ++j)
        fv[j] = f(lo + (hi - lo) * (std::cos(PI * (j + 0.5) / (M + 1)) + 1) / 2);
    for (int k = 0; k <= M; ++k) {
        double s = 0;
        for (int j = 0; j <= M; ++j)
            s += fv[j] * std::cos(PI * k * (j + 0.5) / (M + 1));
        c[k] = s * (k == 0 ? 1.0 : 2.0) / (M + 1);
    }
    return c;
}
template <class F>
VecQ chebfit(F f, double lo, double hi) {
    VecQ fv{}, c{};
    for (int j = 0; j <= NQ; ++j)
        fv[j] = f(lo + (hi - lo) * (std::cos(PI * (j + 0.5) / (NQ + 1)) + 1) / 2);
    for (int k = 0; k <= NQ; ++k) {
        double s = 0;
        for (int j = 0; j <= NQ; ++j)
            s += fv[j] * std::cos(PI * k * (j + 0.5) / (NQ + 1));
        c[k] = s * (k == 0 ? 1.0 : 2.0) / (NQ + 1);
    }
    return c;
}
inline VecG times_x(const VecQ& c, double lo, double hi) {  // (m + h u) * Q exactly
    const double m = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
    VecG r{};
    for (int k = 0; k <= NQ; ++k) {
        r[k] += m * c[k];
        if (k == 0)
            r[1] += h * c[0];
        else {
            r[k + 1] += 0.5 * h * c[k];
            r[k - 1] += 0.5 * h * c[k];
        }
    }
    return r;
}
inline double clenshaw_n(const VecG& c, double u, int n = NG) {
    double b1 = 0, b2 = 0;
    for (int k = n; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + u * b1 - b2;
}
inline void clenshaw_fd(const VecG& c, double u, int n, double& f, double& df) {
    double b1 = 0, b2 = 0, d1 = 0, d2 = 0;
    for (int k = n; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2, d0 = 2 * b1 + 2 * u * d1 - d2;
        b2 = b1;
        b1 = b0;
        d2 = d1;
        d1 = d0;
    }
    f = c[0] + u * b1 - b2;
    df = b1 + u * d1 - d2;
}
inline double l1tail(const VecG& c, int n) {
    double s = 0;
    for (int k = 1; k <= n; ++k)
        s += std::abs(c[k]);
    return s;
}
inline VecG chebder(const VecG& c, int n) {
    VecG d{};
    if (n == 0) return d;
    d[n - 1] = 2 * n * c[n];
    if (n >= 2) d[n - 2] = 2 * (n - 1) * c[n - 1];
    for (int k = n - 3; k >= 0; --k)
        d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
    d[0] *= 0.5;
    return d;
}
// Safeguarded Newton on a Chebyshev series inside a sign-change bracket [a, b] (u coordinates).
inline double refine_cheb(const VecG& c, double a, double b, double fa) {
    ++g_ref_calls;
    double x = 0.5 * (a + b);
    for (int it = 0; it < 60; ++it) {
        ++g_ref_its;
        double f, df;
        clenshaw_fd(c, x, NG, f, df);
        if (f == 0) return x;
        if ((f < 0) == (fa < 0))
            a = x;
        else
            b = x;
        double xn = x - f / df;
        if (!(xn > a && xn < b)) xn = 0.5 * (a + b);
        const double tolx = 4 * DBL_EPSILON * std::max(1.0, std::abs(x));
        if (std::abs(xn - x) <= tolx || b - a <= tolx) return xn;
        x = xn;
    }
    return x;
}

// ------------------------------------------------------------------ Bernstein / Descartes
// Descartes' rule in the Bernstein basis on de Casteljau subintervals.  tol = margin (fit
// error bound of the tables vs the true equation + roundoff) + a rigorous bound on the
// Chebyshev-to-Bernstein conversion roundoff; coefficients within tol have unknown sign.
struct RootOut
{
    double u, ua, ub;  // root and bracket in the piece variable
    bool certified;    // bracket end signs certified for the true equation, and exactly one root
};
constexpr int MAXDEPTH_AMB = 16;  // ambiguous signs (|G| within tolerance): subdividing cannot resolve it; bound the work
constexpr int MAXDEPTH = 48;      // unambiguous but >= 2 sign changes: a close pair, separate it down to ~ulp
constexpr int MAXROOTS = 64;
inline void bern_rec(const VecG& b, const VecG& ctop, double ua, double ub, double tol, int depth, RootOut* out, int& nr) {
    ++g_nodes;
    bool amb = false, anypos = false, anyneg = false;
    int V = 0, last = 0;
    for (int j = 0; j <= NG; ++j) {
        const int sgn = b[j] > tol ? 1 : b[j] < -tol ? -1 : 0;
        if (sgn == 0) {
            amb = true;
            continue;
        }
        (sgn > 0 ? anypos : anyneg) = true;
        if (last != 0 && sgn != last) ++V;
        last = sgn;
    }
    if (!amb && !(anypos && anyneg)) return;  // convex hull excludes zero
    if (!amb && V == 0) return;
    const bool ends_opposite = (b[0] < -tol && b[NG] > tol) || (b[0] > tol && b[NG] < -tol);
    if (!amb && V == 1 && ends_opposite) {  // exactly one root, certified
        if (nr >= MAXROOTS) return;
        out[nr++] = {refine_cheb(ctop, ua, ub, clenshaw_n(ctop, ua)), ua, ub, true};
        return;
    }
    if ((amb && depth >= MAXDEPTH_AMB) || depth >= MAXDEPTH) {  // unresolved: ambiguous signs, or a pair not yet separated
        ++g_uncertain;
        const double fa = clenshaw_n(ctop, ua), fb = clenshaw_n(ctop, ub);
        if ((fa < 0) != (fb < 0) && nr < MAXROOTS) out[nr++] = {refine_cheb(ctop, ua, ub, fa), ua, ub, false};
        return;
    }
    VecG L{}, R{}, t = b;
    L[0] = b[0];
    R[NG] = b[NG];
    for (int r = 1; r <= NG; ++r) {
        for (int i = 0; i <= NG - r; ++i)
            t[i] = 0.5 * (t[i] + t[i + 1]);
        L[r] = t[0];
        R[NG - r] = t[NG - r];
    }
    const double um = 0.5 * (ua + ub);
    bern_rec(L, ctop, ua, um, tol, depth + 1, out, nr);
    bern_rec(R, ctop, um, ub, tol, depth + 1, out, nr);
}
inline int bern_roots(const VecG& c, double margin, RootOut* out) {
    VecG b{};
    double err = 0;
    for (int j = 0; j <= NG; ++j) {
        double s = 0, sa = 0;
        for (int k = 0; k <= NG; ++k) {
            s += GERGCHEB_C2B[j][k] * c[k];
            sa += std::abs(GERGCHEB_C2B[j][k] * c[k]);
        }
        b[j] = s;
        err = std::max(err, 2.3e-16 * (NG + 2) * sa);
    }
    int nr = 0;
    bern_rec(b, c, -1, 1, margin + err, 0, out, nr);
    return nr;
}

// ------------------------------------------------------------------ the solver
struct Solver
{
    std::vector<Term> terms;
    std::vector<int> group;  // term -> distinct delta-function
    std::vector<Term> reps;  // one representative per group
    std::vector<double> edges;
    std::vector<VecQ> C;     // [p * NGRP + g], piece-major
    std::vector<double> Cn;  // l1 norm of each fit (roundoff scale)
    std::vector<double> Ct;  // tail estimate of each fit (fit-error scale)
    std::vector<std::vector<double>> F;
    shared_ptr<CoolProp::ReducingFunction> Red;
    double R = 0, delta_max = 0, tau_max = 0;
    int N = 0, P = 0;
    // Non-analytic (critical-region) terms, e.g. IAPWS-95 water and Span-Wagner CO2.  Not separable
    // in (tau, delta): added per (T, x) by fitting at the Chebyshev nodes when active.
    struct NAComp
    {
        int i;
        const CoolProp::ResidualHelmholtzNonAnalytic* na;
        double Dmin;  // smallest D: the terms carry exp(-D (tau-1)^2)
        double Cmin;  // smallest C: ... and exp(-C (delta-1)^2)
        double nsum;  // sum |n|
    };
    double table_tol = 1e-6;
    std::vector<NAComp> nacomps;
    std::vector<double> Ri;  // component gas constants (mole-fraction average unless normalized)
    bool R_normalized = true;
    // delta d(alpha_NA)/d(delta) and its delta-derivative for one component.  Closed-form first and
    // second delta-derivatives of alpha = n Delta^b delta psi (IAPWS-95 Table 5 form), instead of
    // CoolProp's all(), which computes everything to 4th order.  Checked against all() in Exp. 7.
    static void na_z(const NAComp& c, double tau, double D, double& f, double& df) {
        double ad = 0, add = 0;
        const double dm = (std::abs(D - 1) < 10 * DBL_EPSILON) ? 10 * DBL_EPSILON : D - 1, d2 = dm * dm;  // as CoolProp
        const double tm = (std::abs(tau - 1) < 10 * DBL_EPSILON) ? 1 + 10 * DBL_EPSILON : tau;
        const double L = std::log(d2);  // all ((d-1)^2)^s as exp(s L)
        // IAPWS-95 and Span-Wagner share a = 3.5, beta = 0.3 (and often A, B) across their terms: reuse
        double la = NAN, lbeta = NAN, lA = NAN, lB = NAN, pw = 0, pa = 0, theta = 0, Delta = 0, dDelta = 0, d2Delta = 0, lnDelta = 0;
        for (const auto& el : c.na->elements) {
            const double ib = 1 / (2 * el.beta);
            const bool same_ab = el.a == la && el.beta == lbeta;
            if (!same_ab) {
                pw = std::exp((ib - 1) * L);    // ((d-1)^2)^(1/(2 beta) - 1)
                pa = std::exp((el.a - 1) * L);  // ((d-1)^2)^(a-1)
            }
            if (!same_ab || el.A != lA || el.B != lB) {
                theta = (1 - tm) + el.A * pw * d2;
                Delta = theta * theta + el.B * pa * d2;
                dDelta = dm * (el.A * theta * (2 / el.beta) * pw + 2 * el.B * el.a * pa);
                d2Delta = dDelta / dm
                          + d2
                              * (4 * el.B * el.a * (el.a - 1) * pa / d2 + 2 * el.A * el.A / (el.beta * el.beta) * pw * pw
                                 + el.A * theta * 4 / el.beta * (ib - 1) * pw / d2);
                lnDelta = std::log(Delta);
            }
            la = el.a;
            lbeta = el.beta;
            lA = el.A;
            lB = el.B;
            const double Db = std::exp(el.b * lnDelta), dDb = el.b * Db / Delta * dDelta,
                         d2Db = el.b * (Db / Delta * d2Delta + (el.b - 1) * Db / (Delta * Delta) * dDelta * dDelta);
            const double psi = std::exp(-el.C * d2 - el.D * (tm - 1) * (tm - 1)), dpsi = -2 * el.C * dm * psi,
                         d2psi = (2 * el.C * d2 - 1) * 2 * el.C * psi;
            ad += el.n * (Db * (psi + D * dpsi) + dDb * D * psi);
            add += el.n * (Db * (2 * dpsi + D * d2psi) + 2 * dDb * (psi + D * dpsi) + d2Db * D * psi);
        }
        f = D * ad;
        df = ad + D * add;
    }
    static void na_z_coolprop(const NAComp& c, double tau, double D, double& f, double& df) {  // reference
        CoolProp::HelmholtzDerivatives d;
        const_cast<CoolProp::ResidualHelmholtzNonAnalytic*>(c.na)->all(tau, D, d);  // reads only its own members
        f = D * d.dalphar_ddelta;
        df = d.dalphar_ddelta + D * d.d2alphar_ddelta2;
    }

    // Tier 1: component set.  tau_max bounds Tr(x)/T over the intended use.
    void build(CoolProp::HelmholtzEOSMixtureBackend* HEOS, double dmax, double taumax, double tol) {
        delta_max = dmax;
        tau_max = taumax;
        table_tol = tol;
        Red = HEOS->Reducing;
        R = HEOS->gas_constant();
        auto& comps = HEOS->get_components();
        N = static_cast<int>(comps.size());
        auto& Ex = HEOS->residual_helmholtz->Excess;
        F = Ex.F;
        for (int i = 0; i < N; ++i)
            collect(comps[i].EOS().alphar.GenExp, i, -1, terms);
        R_normalized = CoolProp::get_config_bool(NORMALIZE_GAS_CONSTANTS);
        for (int i = 0; i < N; ++i) {
            auto& ar = comps[i].EOS().alphar;
            Ri.push_back(comps[i].gas_constant());
            if (ar.NonAnalytic.N > 0) {
                double Dmin = 1e300, Cmin = 1e300, nsum = 0;
                for (const auto& el : ar.NonAnalytic.elements) {
                    Dmin = std::min(Dmin, el.D);
                    Cmin = std::min(Cmin, el.C);
                    nsum += std::abs(el.n);
                }
                nacomps.push_back({i, &ar.NonAnalytic, Dmin, Cmin, nsum});
            }
            // every residual contribution must be GenExp or NonAnalytic, else refuse
            for (double tau : {0.7, 1.3})
                for (double D : {0.3, 1.7}) {
                    CoolProp::HelmholtzDerivatives full = ar.all(tau, D, false), ge, na;
                    ar.GenExp.all(tau, D, ge);
                    const_cast<CoolProp::ResidualHelmholtzNonAnalytic&>(ar.NonAnalytic).all(tau, D, na);
                    if (std::abs(full.alphar - ge.alphar - na.alphar) > 1e-12 * (1 + std::abs(full.alphar)))
                        throw std::runtime_error("component " + std::to_string(i) + " has residual terms other than GenExp/NonAnalytic");
                }
        }
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                if (Ex.F[i][j] != 0 && Ex.DepartureFunctionMatrix[i][j]) collect(Ex.DepartureFunctionMatrix[i][j]->phi, i, j, terms);
        for (const auto& T : terms) {
            const auto it = std::find_if(reps.begin(), reps.end(), [&](const Term& a) { return T.same_delta(a); });
            if (it == reps.end()) {
                group.push_back(static_cast<int>(reps.size()));
                reps.push_back(T);
            } else
                group.push_back(static_cast<int>(it - reps.begin()));
        }
        std::vector<double> gmax(reps.size(), 0.0);
        for (std::size_t k = 0; k < terms.size(); ++k) {
            double km = 0;
            for (int s = 1; s <= 400; ++s) {
                const double tau = tau_max * s / 400.0;
                km = std::max(km, std::abs(terms[k].kappa(tau, std::log(tau))));
            }
            gmax[group[k]] += km;  // |x factor| <= 1 and |F_ij| <= 1 in GERG
        }
        // Shared pieces.  Allowed fit error of each weighted group: tol in absolute units of Z,
        // relaxed only to the roundoff floor of the group's own largest value on the piece (no floor
        // on the piece touching delta = 0).  At low T the GERG terms cancel by up to ~1e9 relative to
        // Z, so a budget relative to the term size lets the tables err by more than Z (Exp. 7).
        std::vector<std::pair<double, double>> todo = {{0.0, dmax}}, done;
        while (!todo.empty()) {
            auto [lo, hi] = todo.back();
            todo.pop_back();
            bool ok = true;
            for (std::size_t g = 0; g < reps.size() && ok; ++g) {
                const auto& T = reps[g];
                const VecQ c = chebfit([&](double D) { return T.chi(D); }, lo, hi);
                double big = 0;
                for (double v : c)
                    big = std::max(big, std::abs(v));
#ifndef GERGCHEB_ABSOLUTE_BUDGET  // default: tol relative to the local term magnitude: tol relative to the local term magnitude (fails when terms cancel, see README Exp. 7)
                const double loc = std::min(std::abs(T.chi(lo > 0 ? lo : 1e-300)), std::abs(T.chi(hi)));
                const double allowed = std::max(tol * std::max(1.0, gmax[g] * loc), lo > 0 ? 4 * DBL_EPSILON * gmax[g] * big : 0.0);
#else  // absolute budget in units of Z, relaxed only to the roundoff floor (274 MB of tables in Exp. 7: not viable)
                const double allowed = std::max(tol, lo > 0 ? 4 * DBL_EPSILON * gmax[g] * big : 0.0);
#endif
                ok = (std::abs(c[NQ]) + std::abs(c[NQ - 1])) * gmax[g] <= allowed;
            }
            if (ok || hi - lo < 1e-3)
                done.push_back({lo, hi});
            else {
                todo.push_back({0.5 * (lo + hi), hi});
                todo.push_back({lo, 0.5 * (lo + hi)});
            }
        }
        std::sort(done.begin(), done.end());
        if (!nacomps.empty())  // non-analytic at delta = 1: make it a piece edge
            for (std::size_t k = 0; k < done.size(); ++k)
                if (done[k].first < 1.0 && done[k].second > 1.0) {
                    const auto pr = done[k];
                    done[k] = {pr.first, 1.0};
                    done.insert(done.begin() + k + 1, {1.0, pr.second});
                    break;
                }
        edges = {0.0};
        for (auto& pr : done)
            edges.push_back(pr.second);
        P = static_cast<int>(done.size());
        const std::size_t NGRP = reps.size();
        C.resize(NGRP * P);
        Cn.resize(NGRP * P);
        Ct.resize(NGRP * P);
        for (std::size_t g = 0; g < NGRP; ++g)
            for (int p = 0; p < P; ++p) {
                const VecQ c = chebfit([&](double D) { return reps[g].chi(D); }, edges[p], edges[p + 1]);
                double s = 0;
                for (double v : c)
                    s += std::abs(v);
                C[p * NGRP + g] = c;
                Cn[p * NGRP + g] = s;
                // fit-error bound: measured max |fit - chi| on 64 interior points, never below the
                // tail estimate (the tail alone underestimated it: unflagged spinodal misses, Exp. 7)
                double meas = 0;
                VecG cg{};
                for (int m = 0; m <= NQ; ++m)
                    cg[m] = c[m];
                for (int j = 0; j < 64; ++j) {
                    const double u = -1 + 2 * (j + 0.5) / 64, D = edges[p] + (edges[p + 1] - edges[p]) * (u + 1) / 2;
                    meas = std::max(meas, std::abs(clenshaw_n(cg, u, NQ) - reps[g].chi(D)));
                }
                Ct[p * NGRP + g] = std::max(meas, std::abs(c[NQ]) + std::abs(c[NQ - 1]));
            }
    }

    // Tier 2: per (T, x)
    struct State
    {
        std::vector<VecG> G;
        std::vector<double> margin;
        std::vector<double> W;
        std::vector<double> x;
        std::vector<int> na_on;  // indices into nacomps active at this (T, x)
        double T = 0, tau = 0, rhor = 0, t_scale = 0;
    };
    void assemble(double T, const std::vector<double>& x, State& S) const {
        const double Tr = Red->Tr(x), rhor = Red->rhormolar(x), tau = Tr / T, lt = std::log(tau);
        const std::size_t NGRP = reps.size();
        S.W.assign(NGRP, 0.0);
        for (std::size_t k = 0; k < terms.size(); ++k) {
            const Term& tm = terms[k];
            const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * F[tm.i][tm.j];
            S.W[group[k]] += X * tm.kappa(tau, lt);
        }
        S.x = x;
        S.tau = tau;
        S.na_on.clear();
        for (std::size_t ci = 0; ci < nacomps.size(); ++ci)
            if (x[nacomps[ci].i] > 0 && nacomps[ci].Dmin * (tau - 1) * (tau - 1) < 50) S.na_on.push_back(static_cast<int>(ci));
        if (!S.na_on.empty()) ++g_na_active;
        S.G.resize(P);
        S.margin.resize(P);
        for (int p = 0; p < P; ++p) {
            VecQ q{};
            q[0] = 1.0;
            double scale = 1.0, fiterr = 0.0;
            const VecQ* cp = &C[p * NGRP];
            const double* cn = &Cn[p * NGRP];
            const double* ct = &Ct[p * NGRP];
            for (std::size_t g = 0; g < NGRP; ++g) {
                const double w = S.W[g];
                for (int m = 0; m <= NQ; ++m)
                    q[m] += w * cp[g][m];
                scale += std::abs(w) * cn[g];
                fiterr += std::abs(w) * ct[g];
            }
            for (int ci : S.na_on) {  // non-analytic add-in: fit at the nodes of this piece
                const NAComp& c = nacomps[ci];
                const double xi = x[c.i];
                {  // skip where negligible: psi_max times a generous polynomial factor, into the margin instead
                    const double lo = edges[p], hi = edges[p + 1];
                    const double dmin = (lo <= 1 && hi >= 1) ? 0.0 : std::min(std::abs(lo - 1), std::abs(hi - 1));
                    const double psimax = std::exp(-c.Cmin * dmin * dmin - c.Dmin * (tau - 1) * (tau - 1));
                    const double bound = 1e3 * xi * c.nsum * psimax * hi * (1 + hi) * (1 + 2 * c.Cmin * (1 + hi));
                    if (bound < 1e-3 * table_tol) {
                        fiterr += bound;
                        continue;
                    }
                }
                auto fna = [&](double D) {
                    double f, df;
                    na_z(c, tau, D, f, df);
                    return xi * f;
                };
                VecQ cf{};
                {  // degree 8 first; degree NQ only if its tail is not negligible against the table tolerance
                    const auto c8 = chebfit_n<8>(fna, edges[p], edges[p + 1]);
                    const double tail8 = std::abs(c8[8]) + std::abs(c8[7]);
                    if (tail8 < 1e-3 * table_tol) {
                        for (int m = 0; m <= 8; ++m)
                            cf[m] = c8[m];
                        ++g_na_deg8;
                    } else
                        cf = chebfit(fna, edges[p], edges[p + 1]);
                }
                double s1 = 0;
                for (int m = 0; m <= NQ; ++m) {
                    q[m] += cf[m];
                    s1 += std::abs(cf[m]);
                }
                scale += s1;
                // Fit-error bound: MEASURED at points between the fitting nodes, not a tail estimate alone -- the
                // slowly decaying CO2 terms (C = 10) on wide pieces beat the tail estimate by ~100x (Exp. 10).
                double meas = 0;
                {
                    VecG cg{};
                    for (int m = 0; m <= NQ; ++m)
                        cg[m] = cf[m];
                    for (double u : {-0.97, -0.62, -0.21, 0.21, 0.62, 0.97}) {
                        const double D = edges[p] + (edges[p + 1] - edges[p]) * (u + 1) / 2;
                        meas = std::max(meas, std::abs(clenshaw_n(cg, u, NQ) - fna(D)));
                    }
                }
                fiterr += std::max(10 * meas, 5 * (std::abs(cf[NQ]) + std::abs(cf[NQ - 1])));
            }
            S.G[p] = times_x(q, edges[p], edges[p + 1]);
            // bound on |G_tables - G_true| on the piece: 2 x fit error bound + roundoff, times delta <= hi
            S.margin[p] = edges[p + 1] * (2 * fiterr + 1e-14 * scale);
        }
        double Rmix = R;
        if (!R_normalized) {
            Rmix = 0;
            for (int i = 0; i < N; ++i)
                Rmix += x[i] * Ri[i];
        }
        S.T = T;
        S.rhor = rhor;
        S.t_scale = 1.0 / (rhor * Rmix * T);
    }

    // true equation from the grouped terms (exact regrouping of the model)
    void true_G(const State& S, double D, double t, double& G, double& dG, double& scale) const {
        if (D <= 0) {  // limit: G(0) = -t exactly (the chi_g vanish like delta^d; direct evaluation is 0/0)
            G = -t;
            dG = 1;
            scale = t;
            return;
        }
        double Z = 1, dZ = 0, a = 0;
        double pw[25];
        pw[0] = 1;
        for (int k = 1; k < 25; ++k)
            pw[k] = pw[k - 1] * D;
        for (std::size_t g = 0; g < reps.size(); ++g) {
            double f, df;
            reps[g].chi_d_pw(D, pw, f, df);
            Z += S.W[g] * f;
            dZ += S.W[g] * df;
            a += std::abs(S.W[g] * f);
        }
        for (int ci : S.na_on) {
            double f, df;
            na_z(nacomps[ci], S.tau, D, f, df);
            const double xi = S.x[nacomps[ci].i];
            Z += xi * f;
            dZ += xi * df;
            a += std::abs(xi * f);
        }
        G = D * Z - t;
        dG = Z + D * dZ;
        scale = D * (1 + a) + t;
    }

    struct Root
    {
        double rho;       // polished if polished == true, else the table root
        double rho_cheb;  // from the tables
        double resid;     // backward error in delta after polish (0 if not polished)
        bool certified;   // bracket certified by the subdivision
        bool bracket_ok;  // polish stayed inside a sign-checked bracket
        bool polished;
        double Da, Db, sa;  // bracket in delta and the table sign at Da, kept for a later polish()
    };
    enum class Polish
    {
        None,
        All
    };
    // Tier 3: per p.  Returns the number of roots in (0, delta_max].  With Polish::None the roots are
    // the table roots (accuracy ~ table tolerance); polish(S, p, root) refines a chosen one later.
    int roots(const State& S, double p, Root* out, Polish mode = Polish::All) const {
        const double t = p * S.t_scale;
        ++g_calls;
        int n = 0;
        for (int pc = 0; pc < P; ++pc) {
            VecG g = S.G[pc];
            g[0] -= t;
            ++g_nodes;
            if (std::abs(g[0]) > l1tail(g, NG) + S.margin[pc]) continue;
            ++g_pieces_open;
            RootOut ro[MAXROOTS];
            const int nr = bern_roots(g, S.margin[pc], ro);
            const double lo = edges[pc], w = edges[pc + 1] - edges[pc];
            for (int k = 0; k < nr; ++k) {
                const double D = lo + w * (ro[k].u + 1) / 2;
                if (n > 0 && std::abs(D * S.rhor - out[n - 1].rho_cheb) < 1e-10 * out[n - 1].rho_cheb) continue;
                Root r{D * S.rhor,
                       D * S.rhor,
                       0,
                       ro[k].certified,
                       false,
                       false,
                       lo + w * (ro[k].ua + 1) / 2,
                       lo + w * (ro[k].ub + 1) / 2,
                       clenshaw_n(g, ro[k].ua)};
                if (!r.certified) {  // not certified: keep it only if the TRUE G changes sign across its bracket
                    double Ga, Gb, dG, sc;
                    true_G(S, r.Da, t, Ga, dG, sc);
                    true_G(S, r.Db, t, Gb, dG, sc);
                    if ((Ga < 0) == (Gb < 0)) {
                        ++g_dropped;
                        continue;
                    }
                }
                if (mode == Polish::All) polish(S, p, r);
                if (n < MAXROOTS) out[n++] = r;
            }
        }
        return n;
    }
    void polish(const State& S, double p, Root& r) const {
        if (r.polished) return;
        r.rho = S.rhor * polish_root(S, p * S.t_scale, r.rho_cheb / S.rhor, r.Da, r.Db, r.sa, r.resid, r.bracket_ok, r.certified);
        r.polished = true;
    }
    // Bracketed Newton on the true equation: the bracket shrinks with every true-G sign, a step
    // leaving it becomes a bisection.  bracket_ok = the true G has the expected sign at Da.
    double polish_root(const State& S, double t, double D, double Da, double Db, double sa, double& resid, bool& bracket_ok,
                       bool certified = false) const {
        ++g_pol_calls;
        const double D0 = D;
        double G, dG, sc;
        if (certified)  // the subdivision's margin bounds |G_table - G_true|, so the sign at Da is already known
            bracket_ok = true;
        else {
            true_G(S, Da, t, G, dG, sc);
            bracket_ok = (G < 0) == (sa < 0);
        }
        if (!bracket_ok) ++g_pol_uncert;
        double a = Da, b = Db;
        for (int it = 0; it < 20; ++it) {
            ++g_pol_its;
            true_G(S, D, t, G, dG, sc);
            resid = std::abs(G) / (std::abs(dG) * D + 1e-300);  // backward error: relative change in delta that zeroes G
            if (G == 0) break;
            if ((G < 0) == (sa < 0))
                a = D;
            else
                b = D;
            const double step = G / dG;
            if (std::abs(step) <= 2 * DBL_EPSILON * D) {  // converged: test before the bracket check, since
                D -= step;                                // D itself was just made a bracket end
                break;
            }
            double Dn = D - step;
            if (bracket_ok && !(Dn > a && Dn < b)) Dn = 0.5 * (a + b);
            D = Dn;
        }
        true_G(S, D, t, G, dG, sc);
        resid = std::abs(G) / (std::abs(dG) * D + 1e-300);      // backward error: relative change in delta that zeroes G
        if (!bracket_ok && !(std::abs(D - D0) <= 1e-3 * D0)) {  // unbracketed Newton wandered: keep the table root
            D = D0;
            true_G(S, D, t, G, dG, sc);
            resid = std::abs(G) / (std::abs(dG) * D + 1e-300);  // backward error: relative change in delta that zeroes G
        }
        return D;
    }

    // term-by-term evaluation (independent of the grouping), for spot checks
    double G_terms(double T, const std::vector<double>& x, double D, double p) const {
        const double Tr = Red->Tr(x), rhor = Red->rhormolar(x), tau = Tr / T, lt = std::log(tau);
        double Z = 1;
        for (const auto& tm : terms) {
            const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * F[tm.i][tm.j];
            Z += X * tm.kappa(tau, lt) * tm.chi(D);
        }
        return D * Z - p / (rhor * R * T);
    }
};

}  // namespace gergcheb
