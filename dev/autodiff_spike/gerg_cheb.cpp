// Experiment 6: Chebyshev all-roots density solver for the GERG-2008 multi-fluid model.
//
// Every GERG-2008 residual term (pure-fluid and departure) is a generalized-exponential
// element that factors exactly as  n tau^t e^{u_tau(tau)} * delta^d e^{u_delta(delta)}.
// At fixed (T, x) the pressure equation p = rho R T Z with rho = delta rho_r(x) becomes
//     G(delta) = delta Z(delta) - t = 0,   t = p / (rho_r(x) R T),
//     Z = 1 + sum_k X_k(x) kappa_k(tau) chi_k(delta),
//     chi_k = delta d/d(delta) [delta^d e^{u_delta}] = delta^d e^{u_delta} (d + delta u_delta')
// with X_k = x_i (pure-fluid term) or x_i x_j F_ij (departure term).  Tiers:
//   component set : collect the terms; choose shared delta-pieces; Chebyshev-fit chi_k on each
//   per (T, x)    : tau, rho_r; weights W_k = X_k kappa_k(tau); Q = 1 + sum W_k chi_k per piece
//   per p         : G = delta Q - t (delta applied exactly); certified subdivision for roots
//
// Build against a Release static CoolProp (see README, Experiment 3).
#include <algorithm>
#include <array>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include "AbstractState.h"
#include "cheb2bern.hpp"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

namespace {
constexpr double PI = 3.14159265358979323846;
#ifndef NQ_DEG
#    define NQ_DEG 12
#endif
constexpr int NQ = NQ_DEG, NG = NQ + 1;
using VecQ = std::array<double, NQ + 1>;
using VecG = std::array<double, NG + 1>;

// ------------------------------------------------------------------ terms
struct Term
{
    int i, j;  // component (j = -1 for pure-fluid terms)
    double n, d, t;
    bool has_cl, has_om, has_e1, has_e2, has_b1, has_b2;
    double c, l, om, m, e1, eps1, e2, eps2, b1, g1, b2, g2;
    // delta part
    double chi(double D) const {
        double u = 0, du = 0;
        if (has_cl) {
            const double dl = std::pow(D, l);
            u -= c * dl;
            du -= c * l * dl / D;
        }
        if (has_e1) {
            u -= e1 * (D - eps1);
            du -= e1;
        }
        if (has_e2) {
            u -= e2 * (D - eps2) * (D - eps2);
            du -= 2 * e2 * (D - eps2);
        }
        return std::pow(D, d) * std::exp(u) * (d + D * du);
    }
    // chi and d(chi)/d(delta): chi = delta^d e^u (d + delta u'),  chi' = delta^(d-1) e^u [a^2 + delta u' + delta^2 u''], a = d + delta u'
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
    // tau part
    double kappa(double tau, double ltau) const {
        double u = t * ltau;
        if (has_om) u -= om * std::exp(m * ltau);
        if (has_b1) u -= b1 * (tau - g1);
        if (has_b2) u -= b2 * (tau - g2) * (tau - g2);
        return n * std::exp(u);
    }
};
static bool finite(double v) {
    return std::isfinite(v);
}
void collect(const CoolProp::ResidualHelmholtzGeneralizedExponential& g, int i, int j, std::vector<Term>& out) {
    for (const auto& el : g.elements) {
        Term T{};
        T.i = i;
        T.j = j;
        T.n = el.n;
        T.d = el.d;
        T.t = el.t;
        T.has_cl = g.delta_li_in_u && finite(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON;
        T.c = el.c;
        T.l = el.l_double;
        T.has_om = g.tau_mi_in_u && std::abs(el.m_double) > 0;
        T.om = el.omega;
        T.m = el.m_double;
        T.has_e1 = g.eta1_in_u && finite(el.eta1);
        T.e1 = el.eta1;
        T.eps1 = el.epsilon1;
        T.has_e2 = g.eta2_in_u && finite(el.eta2);
        T.e2 = el.eta2;
        T.eps2 = el.epsilon2;
        T.has_b1 = g.beta1_in_u && finite(el.beta1);
        T.b1 = el.beta1;
        T.g1 = el.gamma1;
        T.has_b2 = g.beta2_in_u && finite(el.beta2);
        T.b2 = el.beta2;
        T.g2 = el.gamma2;
        out.push_back(T);
    }
}

// ------------------------------------------------------------------ Chebyshev utilities
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
VecG times_x(const VecQ& c, double lo, double hi) {  // (m + h u) * Q exactly
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
double clenshaw(const VecG& c, double u) {
    double b1 = 0, b2 = 0;
    for (int k = NG; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + u * b1 - b2;
}
double l1tail(const VecG& c, int n) {
    double s = 0;
    for (int k = 1; k <= n; ++k)
        s += std::abs(c[k]);
    return s;
}

// ------------------------------------------------------------------ certified subdivision
struct HalfMaps
{
    double L[NG + 1][NG + 1], Rm[NG + 1][NG + 1];
    HalfMaps() {
        for (int side = 0; side < 2; ++side)
            for (int col = 0; col <= NG; ++col) {
                double f[NG + 1];
                for (int j = 0; j <= NG; ++j) {
                    const double x = std::cos(PI * (j + 0.5) / (NG + 1)), u = side == 0 ? (x - 1) / 2 : (x + 1) / 2;
                    f[j] = std::cos(col * std::acos(u));
                }
                for (int k = 0; k <= NG; ++k) {
                    double s = 0;
                    for (int j = 0; j <= NG; ++j)
                        s += f[j] * std::cos(PI * k * (j + 0.5) / (NG + 1));
                    (side == 0 ? L : Rm)[k][col] = s * (k == 0 ? 1.0 : 2.0) / (NG + 1);
                }
            }
    }
};
const HalfMaps HM;
// re-expansion of a degree-n series onto a half; the image of T_j only involves T_0..T_j
VecG half(const VecG& c, bool right, int n) {
    const auto& A = right ? HM.Rm : HM.L;
    VecG r{};
    for (int k = 0; k <= n; ++k) {
        double s = 0;
        for (int j = k; j <= n; ++j)
            s += A[k][j] * c[j];
        r[k] = s;
    }
    return r;
}
VecG chebder(const VecG& c, int n) {
    VecG d{};
    if (n == 0) return d;
    d[n - 1] = 2 * n * c[n];
    if (n >= 2) d[n - 2] = 2 * (n - 1) * c[n - 1];
    for (int k = n - 3; k >= 0; --k)
        d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
    d[0] *= 0.5;
    return d;
}
double clenshaw_n(const VecG& c, double u, int n) {
    double b1 = 0, b2 = 0;
    for (int k = n; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + u * b1 - b2;
}
long g_ill_its = 0, g_ill_calls = 0;
// value and derivative of a Chebyshev series (Clenshaw for f and f')
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
// Safeguarded Newton inside a sign-change bracket [a, b]: take the Newton step when it stays
// inside the bracket, else bisect; shrink the bracket every iteration.  Stops at ~4 ulp.
double illinois(const VecG& c, int n, double a, double b, double fa, double fb) {
    ++g_ill_calls;
    (void)fb;
    double x = 0.5 * (a + b);
    for (int it = 0; it < 60; ++it) {
        ++g_ill_its;
        double f, df;
        clenshaw_fd(c, x, n, f, df);
        if (f == 0) return x;
        if ((f < 0) == (fa < 0)) {
            a = x;
            fa = f;
        } else
            b = x;
        double xn = x - f / df;
        if (!(xn > a && xn < b)) xn = 0.5 * (a + b);
        if (std::abs(xn - x) <= 4 * DBL_EPSILON * std::max(1.0, std::abs(x)) || b - a <= 4 * DBL_EPSILON * std::max(1.0, std::abs(x))) return xn;
        x = xn;
    }
    return x;
}
constexpr int MAXDEPTH = 16;
long g_nodes = 0, g_calls = 0, g_pieces_open = 0;
// c has degree n on this subinterval; margin bounds (fit error + roundoff + trimmed mass).
void cert_rec(const VecG& c, int n, double ua, double ub, double margin, int depth, double* out, int& nr) {
    ++g_nodes;
    if (std::abs(c[0]) > l1tail(c, n) + margin) return;  // no root, with margin
    if (depth == 0) ++g_pieces_open;
    const VecG d = chebder(c, n);
    if (n <= 1 || std::abs(d[0]) > l1tail(d, n - 1) || depth >= MAXDEPTH) {  // monotone (or give up refining)
        double fa = 0, fb = 0;
        for (int k = 0; k <= n; ++k) {
            fb += c[k];
            fa += (k % 2 ? -c[k] : c[k]);
        }
        if ((fa < 0) != (fb < 0) || fa == 0) out[nr++] = ua + (ub - ua) * ((fa == 0 ? -1.0 : illinois(c, n, -1, 1, fa, fb)) + 1) / 2;
        return;
    }
    const double um = 0.5 * (ua + ub);
    for (int side = 0; side < 2; ++side) {
        VecG h = half(c, side == 1, n);
        int m = n;
        double trimmed = 0;  // drop trailing coefficients worth < margin/10 in total; carry them in the margin
        while (m > 1 && trimmed + std::abs(h[m]) < 0.1 * margin)
            trimmed += std::abs(h[m--]);
        for (int k = m + 1; k <= n; ++k)
            h[k] = 0;
        cert_rec(h, m, side ? um : ua, side ? ub : um, margin + trimmed, depth + 1, out, nr);
    }
}

// ------------------------------------------------------------------ Bernstein / Descartes rootfinder
// Convert a piece to the Bernstein basis once (exact matrix, rounded), then subdivide by
// de Casteljau.  Descartes' rule in the Bernstein basis: #roots <= #sign changes of the
// coefficients, with equal parity, so 0 changes -> no root and 1 change -> exactly one.
// Coefficients within tol (fit margin + a rigorous bound on the conversion roundoff) have
// unknown sign; a subinterval containing one is subdivided.  Roots are refined by Illinois
// on the original Chebyshev series, which is well conditioned.
#if NQ_DEG == 16
#    define C2B CHEB2BERN_17
#elif NQ_DEG == 20
#    define C2B CHEB2BERN_21
#elif NQ_DEG == 24
#    define C2B CHEB2BERN_25
#endif
#ifdef C2B
void bern_rec(const VecG& b, const VecG& ctop, double ua, double ub, double tol, int depth, double* out, int& nr) {
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
    if ((!amb && V == 1 && ends_opposite) || depth >= MAXDEPTH) {
        if ((b[0] < 0) != (b[NG] < 0) || b[0] == 0) {
            const double fa = clenshaw_n(ctop, ua, NG), fb = clenshaw_n(ctop, ub, NG);
            if ((fa < 0) != (fb < 0) || fa == 0) out[nr++] = fa == 0 ? ua : illinois(ctop, NG, ua, ub, fa, fb);
        }
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
int bern_roots(const VecG& c, double margin, double* out) {
    VecG b{};
    double err = 0;
    for (int j = 0; j <= NG; ++j) {
        double s = 0, sa = 0;
        for (int k = 0; k <= NG; ++k) {
            s += C2B[j][k] * c[k];
            sa += std::abs(C2B[j][k] * c[k]);
        }
        b[j] = s;
        err = std::max(err, 2.3e-16 * (NG + 2) * sa);
    }
    int nr = 0;
    bern_rec(b, c, -1, 1, margin + err, 0, out, nr);
    return nr;
}
#endif

// ------------------------------------------------------------------ the solver
struct Solver
{
    std::vector<Term> terms;
    std::vector<double> edges;
    std::vector<int> group;  // term -> distinct delta-function
    std::vector<Term> reps;  // one representative term per group
    std::vector<VecQ> C;     // [p * NGRP + g]  (piece-major: assembly streams contiguously)
    std::vector<double> Cn;  // l1 norms of C, same layout (for the margin)
    int N = 0, P = 0;
    double delta_max = 0;
    CoolProp::HelmholtzEOSMixtureBackend* H = nullptr;

    // Tier 1: component set
    // tau_max: largest reduced inverse temperature the tables will be used at (lowest T).
    // Acceptance per term and piece: tail(chi_k) * kmax_k <= tol * max(1, kmax_k * max|chi_k|),
    // i.e. the fit error of the weighted term is below tol on the O(1) scale of Z, or below tol
    // relative to the term itself where it is large.  A criterion relative to the term alone
    // never converges for exp(-delta^6) terms, which are ~1e-24 (irrelevant) by delta ~ 2.
    void build(CoolProp::HelmholtzEOSMixtureBackend* HEOS, double dmax, double tau_max, double tol = 1e-14) {
        H = HEOS;
        delta_max = dmax;
        auto& comps = HEOS->get_components();
        N = static_cast<int>(comps.size());
        for (int i = 0; i < N; ++i)
            collect(comps[i].EOS().alphar.GenExp, i, -1, terms);
        auto& Ex = HEOS->residual_helmholtz->Excess;
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                if (Ex.F[i][j] != 0 && Ex.DepartureFunctionMatrix[i][j]) collect(Ex.DepartureFunctionMatrix[i][j]->phi, i, j, terms);
        // Terms whose delta-parameters coincide share chi(delta); tabulate each distinct one once.
        for (const auto& T : terms) {
            auto same = [&](const Term& a) {
                return a.d == T.d && a.has_cl == T.has_cl && (!T.has_cl || (a.c == T.c && a.l == T.l)) && a.has_e1 == T.has_e1
                       && (!T.has_e1 || (a.e1 == T.e1 && a.eps1 == T.eps1)) && a.has_e2 == T.has_e2
                       && (!T.has_e2 || (a.e2 == T.e2 && a.eps2 == T.eps2));
            };
            const auto it = std::find_if(reps.begin(), reps.end(), same);
            if (it == reps.end()) {
                group.push_back(static_cast<int>(reps.size()));
                reps.push_back(T);
            } else
                group.push_back(static_cast<int>(it - reps.begin()));
        }
        std::vector<double> kmax(terms.size(), 0.0);
        for (std::size_t k = 0; k < terms.size(); ++k)
            for (int s = 1; s <= 400; ++s) {  // x-factor <= 1 (|F_ij| <= 1 in GERG), so bound |kappa| over tau
                const double tau = tau_max * s / 400.0;
                kmax[k] = std::max(kmax[k], std::abs(terms[k].kappa(tau, std::log(tau))));
            }
        std::vector<double> gmax(reps.size(), 0.0);
        for (std::size_t k = 0; k < terms.size(); ++k)
            gmax[group[k]] += kmax[k];
        // shared pieces: split until every group meets the acceptance test above
        std::vector<std::pair<double, double>> todo = {{0.0, dmax}}, done;
        while (!todo.empty()) {
            auto [lo, hi] = todo.back();
            todo.pop_back();
            bool ok = true;
            for (std::size_t g = 0; g < reps.size(); ++g) {
                const auto& T = reps[g];
                const VecQ c = chebfit([&](double D) { return T.chi(D); }, lo, hi);
                double big = 0;
                for (double v : c)
                    big = std::max(big, std::abs(v));
                // allowed fit error, in units of Z: tol relative to the *smallest* local magnitude of the
                // weighted group on the piece (ends), never below tol absolute (the O(1) scale of Z near
                // delta = 0, where the gas root lives), with a floor at the unavoidable roundoff of the
                // largest value on the piece.
                const double loc = std::min(std::abs(T.chi(lo > 0 ? lo : 1e-300)), std::abs(T.chi(hi)));
                // (no roundoff floor on the piece touching delta = 0: the terms shrink like delta^d there,
                // so splitting terminates, and the gas root needs the absolute accuracy)
                const double allowed = std::max(tol * std::max(1.0, gmax[g] * loc), lo > 0 ? 4 * DBL_EPSILON * gmax[g] * big : 0.0);
                if ((std::abs(c[NQ]) + std::abs(c[NQ - 1])) * gmax[g] > allowed) {
                    ok = false;
                    break;
                }
            }
            if (ok || hi - lo < 1e-3)
                done.push_back({lo, hi});
            else {
                const double mid = 0.5 * (lo + hi);
                todo.push_back({mid, hi});
                todo.push_back({lo, mid});
            }
        }
        std::sort(done.begin(), done.end());
        edges = {0.0};
        for (auto& pr : done)
            edges.push_back(pr.second);
        P = static_cast<int>(done.size());
        C.resize(reps.size() * P);
        Cn.resize(reps.size() * P);
        const std::size_t NGRP = reps.size();
        for (std::size_t g = 0; g < NGRP; ++g)
            for (int p = 0; p < P; ++p) {
                C[p * NGRP + g] = chebfit([&](double D) { return reps[g].chi(D); }, edges[p], edges[p + 1]);
                double s = 0;
                for (double v : C[p * NGRP + g])
                    s += std::abs(v);
                Cn[p * NGRP + g] = s;
            }
    }

    // Tier 2: per (T, x)
    struct State
    {
        std::vector<VecG> G;
        std::vector<double> margin;
        std::vector<double> W;  // group weights, reused by the polish
        double rhor, t_scale;   // t = p * t_scale
    };
    void assemble(double T, const std::vector<double>& x, State& S, std::vector<double>& W) const {
        const double Tr = H->Reducing->Tr(x), rhor = H->Reducing->rhormolar(x), tau = Tr / T, lt = std::log(tau);
        W.assign(reps.size(), 0.0);
        for (std::size_t k = 0; k < terms.size(); ++k) {
            const Term& tm = terms[k];
            const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * H->residual_helmholtz->Excess.F[tm.i][tm.j];
            W[group[k]] += X * tm.kappa(tau, lt);
        }
        S.G.resize(P);
        S.margin.resize(P);
        for (int p = 0; p < P; ++p) {
            VecQ q{};
            q[0] = 1.0;
            double scale = 1.0;
            const std::size_t NGRP = reps.size();
            const VecQ* cp = &C[p * NGRP];
            const double* cn = &Cn[p * NGRP];
            for (std::size_t g = 0; g < NGRP; ++g) {
                const double w = W[g];
                for (int m = 0; m <= NQ; ++m)
                    q[m] += w * cp[g][m];
                scale += std::abs(w) * cn[g];
            }
            S.G[p] = times_x(q, edges[p], edges[p + 1]);
            S.margin[p] = 1e-13 * scale * edges[p + 1];  // fit tolerance + roundoff, times |delta| <= hi
        }
        S.W = W;
        S.rhor = rhor;
        S.t_scale = 1.0 / (rhor * H->gas_constant() * T);
    }
    // Tier 3: per p; returns molar densities
    int roots(const State& S, double p, double* rho) const {
        const double t = p * S.t_scale;
        ++g_calls;
        int n = 0;
        for (int pc = 0; pc < P; ++pc) {
            VecG g = S.G[pc];
            g[0] -= t;
            double u[64];
            int nr = 0;
#if defined(C2B) && !defined(USE_CHEB_TESTS)
            ++g_nodes;
            if (std::abs(g[0]) > l1tail(g, NG) + S.margin[pc]) continue;  // cheap exclusion before converting
            ++g_pieces_open;
            nr = bern_roots(g, S.margin[pc], u);
#else
            cert_rec(g, NG, -1, 1, S.margin[pc], 0, u, nr);
#endif
            for (int k = 0; k < nr; ++k) {
                const double D = edges[pc] + (edges[pc + 1] - edges[pc]) * (u[k] + 1) / 2;
                if (n > 0 && std::abs(D * S.rhor - rho[n - 1]) < 1e-10 * rho[n - 1]) continue;
                rho[n++] = D * S.rhor;
            }
        }
        if (polish)
            for (int k = 0; k < n; ++k)
                rho[k] = S.rhor * polish_root(S, t, rho[k] / S.rhor);
        return n;
    }
    // Newton on the true equation (grouped term sum with the assembly's weights)
    bool polish = false;
    mutable long pol_its = 0, pol_calls = 0, pol_fail = 0;
    double polish_root(const State& S, double t, double D) const {
        ++pol_calls;
        const double D0 = D;
        for (int it = 0; it < 8; ++it) {
            ++pol_its;
            double Z = 1, dZ = 0;
            for (std::size_t g = 0; g < reps.size(); ++g) {
                double f, df;
                reps[g].chi_d(D, f, df);
                Z += S.W[g] * f;
                dZ += S.W[g] * df;
            }
            const double G = D * Z - t, dG = Z + D * dZ, step = G / dG;
            D -= step;
            if (std::abs(step) <= 2 * DBL_EPSILON * D) break;
        }
        if (!(std::abs(D - D0) <= 1e-3 * D0)) {  // wandered off: keep the Chebyshev root
            ++pol_fail;
            return D0;
        }
        return D;
    }
    // direct evaluation of G(delta) from the terms (reference)
    double G_direct(double T, const std::vector<double>& x, double D, double p) const {
        const double Tr = H->Reducing->Tr(x), rhor = H->Reducing->rhormolar(x), tau = Tr / T, lt = std::log(tau);
        double Z = 1;
        for (const auto& tm : terms) {
            const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * H->residual_helmholtz->Excess.F[tm.i][tm.j];
            Z += X * tm.kappa(tau, lt) * tm.chi(D);
        }
        return D * Z - p / (rhor * H->gas_constant() * T);
    }
};

int ref_roots(const Solver& Sv, double T, const std::vector<double>& x, double p, double dscan, double* rho, int M = 100000) {
    const double rhor = Sv.H->Reducing->rhormolar(x);
    auto g = [&](double D) { return Sv.G_direct(T, x, D, p); };
    int n = 0;
    double d0 = 1e-12, g0 = g(d0);
    for (int j = 1; j <= M; ++j) {
        const double d1 = dscan * j / M, g1 = g(d1);
        if ((g0 < 0) != (g1 < 0)) {
            double a = d0, b = d1, fa = g0;
            for (int it = 0; it < 70; ++it) {
                const double m = 0.5 * (a + b), fm = g(m);
                if ((fm < 0) == (fa < 0)) {
                    a = m;
                    fa = fm;
                } else
                    b = m;
            }
            rho[n++] = 0.5 * (a + b) * rhor;
        }
        d0 = d1;
        g0 = g1;
    }
    return n;
}
}  // namespace

int main(int argc, char** argv) {
    struct Mix
    {
        std::string name, fluids;
        std::vector<double> x;
        std::vector<double> Ts;
    };
    const std::vector<double> Tstd = {100, 120, 150, 180, 200, 250, 300, 400, 500};
    const std::vector<Mix> mixes = {
      {"C1/C2 50/50", "Methane&Ethane", {0.5, 0.5}, Tstd},
      {"C1/C2/C3 50/30/20", "Methane&Ethane&Propane", {0.5, 0.3, 0.2}, Tstd},
      {"natural gas (5)", "Methane&Nitrogen&CarbonDioxide&Ethane&Propane", {0.85, 0.03, 0.02, 0.07, 0.03}, Tstd},
      {"C1/H2S 50/50 (type III)", "Methane&HydrogenSulfide", {0.5, 0.5}, {150, 180, 200, 220, 250, 300, 350, 400, 500}},
      {"humid air",
       "Nitrogen&Oxygen&Argon&CarbonDioxide&Water",
       {0.7710, 0.2069, 0.0092, 0.0004, 0.0125},
       {200, 250, 273.15, 300, 330, 373.15, 450, 550, 700}},
    };
    const std::vector<double> ps = {1e3, 1e4, 1e5, 1e6, 3e6, 1e7, 3e7, 1e8};
    const double DMAX = argc > 1 ? std::atof(argv[1]) : 4.0, DSCAN = 6.0, TOL = argc > 2 ? std::atof(argv[2]) : 1e-14;
    const bool FAST = argc > 3;  // skip the brute-force validation

    std::printf("degree %d, tol %.0e, delta_max = %.1f; reference scan to delta = %.1f\n\n", NQ, TOL, DMAX, DSCAN);
    std::printf("%-26s %5s %6s %6s %6s %9s %9s %8s %9s %9s %11s\n", "mixture", "grps", "pieces", "kB", "roots", "max_rel", "count_ok", ">dmax",
                "asm_us", "root_us", "cp_rhoTp_us");
    for (const auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("GERG2008", mx.fluids));
        auto* HEOS = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        AS->set_mole_fractions(mx.x);
        Solver Sv;
        Sv.polish = std::getenv("POLISH") != nullptr;
        const auto tb = std::chrono::steady_clock::now();
        double Tcmax = 0;
        for (const auto& c : HEOS->get_components())
            Tcmax = std::max(Tcmax, c.EOS().reduce.T);
        const double tau_max = Tcmax / *std::min_element(mx.Ts.begin(), mx.Ts.end());  // generous bound on Tr(x)/T
        Sv.build(HEOS, DMAX, tau_max, TOL);
        const double t_build = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tb).count();

        // 0. the term sum must reproduce CoolProp's pressure
        double worst_p = 0;
        AS->specify_phase(CoolProp::iphase_gas);  // skip mixture phase determination; just evaluate the EOS
        for (double T : {150.0, 300.0, 500.0})
            for (double D : {0.01, 0.5, 1.5, 2.5}) {
                const double rho = D * HEOS->Reducing->rhormolar(mx.x);
                AS->update(CoolProp::DmolarT_INPUTS, rho, T);
                const double p_cp = AS->p(), p_terms = (Sv.G_direct(T, mx.x, D, 0.0)) * HEOS->Reducing->rhormolar(mx.x) * HEOS->gas_constant() * T;
                worst_p = std::max(worst_p, std::abs(p_terms - p_cp) / std::max(std::abs(p_cp), 1.0));
            }
        AS->unspecify_phase();

        // 1. accuracy and root counts vs brute-force scan of the direct term sum
        double worst = 0;
        char wdesc[400] = "";
        int nroots = 0, bad = 0, beyond = 0;
        Solver::State S;
        std::vector<double> W;
        for (double T : mx.Ts) {
            if (FAST) break;
            Sv.assemble(T, mx.x, S, W);
            for (double p : ps) {
                double r1[64], r2[64];
                const int n1 = Sv.roots(S, p, r1), n2all = ref_roots(Sv, T, mx.x, p, DSCAN, r2);
                int n2 = 0;
                for (int k = 0; k < n2all; ++k) {
                    if (r2[k] / S.rhor <= DMAX)
                        r2[n2++] = r2[k];
                    else
                        ++beyond;
                }
                nroots += n2;
                if (n1 != n2) {
                    ++bad;
                    std::printf("   count mismatch T=%g p=%g: cheb %d ref %d  [", T, p, n1, n2);
                    for (int k = 0; k < n1; ++k)
                        std::printf(" %.6g", r1[k] / S.rhor);
                    std::printf(" ] vs [");
                    for (int k = 0; k < n2; ++k)
                        std::printf(" %.6g", r2[k] / S.rhor);
                    std::printf(" ] (delta)\n");
                    continue;
                }
                for (int k = 0; k < n1; ++k) {
                    const double e = std::abs(r1[k] - r2[k]) / r2[k];
                    if (e > worst) {
                        worst = e;
                        std::snprintf(wdesc, sizeof wdesc,
                                      "T=%g p=%g root %d/%d delta_ref=%.10g delta_cheb=%.10g Gdirect(cheb)=%.2e Gdirect(ref)=%.2e", T, p, k + 1, n1,
                                      r2[k] / S.rhor, r1[k] / S.rhor, Sv.G_direct(T, mx.x, r1[k] / S.rhor, p),
                                      Sv.G_direct(T, mx.x, r2[k] / S.rhor, p));
                    }
                }
            }
        }

        if (!FAST) std::printf("   worst: %s\n", wdesc);
        if (Sv.polish && Sv.pol_calls)
            std::printf("   polish: %.2f Newton its per root, %ld of %ld rejected\n", double(Sv.pol_its) / Sv.pol_calls, Sv.pol_fail, Sv.pol_calls);
        // 2. timing
        volatile double sink = 0;
        auto time = [&](auto&& fn, int REP) {
            double best = 1e300;
            for (int rep = 0; rep < 5; ++rep) {
                const auto a = std::chrono::steady_clock::now();
                for (int it = 0; it < REP; ++it)
                    fn(it);
                best = std::min(best, std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - a).count() / REP);
            }
            return best;
        };
        const int nT = static_cast<int>(mx.Ts.size()), nP = static_cast<int>(ps.size());
        const double t_asm = time(
          [&](int it) {
              Sv.assemble(mx.Ts[it % nT] * (1 + 1e-12 * it), mx.x, S, W);
              sink = sink + S.G[0][0];
          },
          20000);
        std::vector<Solver::State> pre(nT);
        for (int i = 0; i < nT; ++i)
            Sv.assemble(mx.Ts[i], mx.x, pre[i], W);
        const double t_root = time(
          [&](int it) {
              double r[64];
              sink = sink + Sv.roots(pre[it % nT], ps[(it / nT) % nP], r);
          },
          20000);
        // CoolProp's own single-root density solver, for context (default guess)
        const double t_cp = time(
          [&](int it) {
              try {
                  sink = sink + HEOS->solver_rho_Tp(mx.Ts[it % nT], ps[(it / nT) % nP]);
              } catch (...) {
              }
          },
          2000);
        std::printf(
          "   per call: %.1f recursion nodes, %.2f pieces not excluded at top level, %d pieces; %.2f roots refined, %.1f Illinois its each\n",
          double(g_nodes) / g_calls, double(g_pieces_open) / g_calls, Sv.P, double(g_ill_calls) / g_calls,
          double(g_ill_its) / std::max(1L, g_ill_calls));
        g_ill_its = g_ill_calls = 0;
        g_nodes = g_calls = g_pieces_open = 0;
        const double kB = Sv.C.size() * sizeof(VecQ) / 1024.0;
        std::printf("%-26s %5zu %6d %6.0f %6d %9.1e %9s %8d %9.2f %9.2f %11.1f   (build %.0f ms, p-check %.0e)\n", mx.name.c_str(), Sv.reps.size(),
                    Sv.P, kB, nroots, worst, bad ? "NO" : "yes", beyond, t_asm, t_root, t_cp, t_build, worst_p);
    }
}
