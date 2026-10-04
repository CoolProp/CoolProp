#include "ChebDensitySolver.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <mutex>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "HelmholtzEOSMixtureBackend.h"
#include "ReducingFunctions.h"
#include "CoolProp/Configuration.h"

namespace CoolProp {
namespace ChebDensity {

// the reducing function and the backend take std::vector<CoolPropDbl>; the mole fractions here are std::vector<double>
static_assert(std::is_same_v<CoolPropDbl, double>, "ChebDensity assumes CoolPropDbl is double");

namespace {

constexpr double PI = 3.14159265358979323846;

/// Chebyshev interpolant of degree M of f on [lo, hi] (values at the Chebyshev points of the first kind)
template <int M, class F>
ChebyshevBernstein::Coeffs<M> chebfit(F f, double lo, double hi) {
    ChebyshevBernstein::Coeffs<M> fv{}, c{};
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

/// (m + h u) q(u) with delta = m + h u on [lo, hi], exactly, as a series of one degree more
CoeffsG times_delta(const CoeffsQ& q, double lo, double hi) {
    const double m = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
    CoeffsG r{};
    for (int k = 0; k <= NQ; ++k) {
        r[k] += m * q[k];
        if (k == 0)
            r[1] += h * q[0];
        else {
            r[k + 1] += 0.5 * h * q[k];
            r[k - 1] += 0.5 * h * q[k];
        }
    }
    return r;
}

// Chebyshev-Lobatto points of degree NQ, u_i = cos(pi i / NQ); those of degree NQ/2, NQ/4 are every 2nd, 4th of them
const std::array<double, NQ + 1>& lobatto_points() {
    static const auto u = [] {
        std::array<double, NQ + 1> r{};
        for (int i = 0; i <= NQ; ++i)
            r[i] = std::cos(PI * i / NQ);
        r[NQ / 2] = 0.0;  // exactly
        return r;
    }();
    return u;
}
// the new points of degree 2 NQ: cos(pi (2j + 1) / (2 NQ)), between the degree-NQ points
const std::array<double, NQ>& lobatto_midpoints() {
    static const auto u = [] {
        std::array<double, NQ> r{};
        for (int j = 0; j < NQ; ++j)
            r[j] = std::cos(PI * (2 * j + 1) / (2 * NQ));
        return r;
    }();
    return u;
}
// Chebyshev coefficients of the degree-n interpolant through values at the degree-n Lobatto points:
// c_k = (2/n) sum''_j f_j cos(pi j k / n), with the j = 0, n terms halved and c_0, c_n halved
template <int n>
const std::array<std::array<double, n + 1>, n + 1>& lobatto_matrix() {
    static const auto M = [] {
        std::array<std::array<double, n + 1>, n + 1> m{};
        for (int k = 0; k <= n; ++k)
            for (int j = 0; j <= n; ++j)
                m[k][j] = 2.0 / n * ((k == 0 || k == n) ? 0.5 : 1.0) * ((j == 0 || j == n) ? 0.5 : 1.0) * std::cos(PI * j * k / n);
        return m;
    }();
    return M;
}
// degree-n interpolant from v, the values at the degree-NQ Lobatto points (every (NQ/n)-th of them is used)
template <int n>
ChebyshevBernstein::Coeffs<n> lobatto_fit(const std::array<double, NQ + 1>& v) {
    static_assert(NQ % n == 0, "the degree must divide NQ");
    const auto& M = lobatto_matrix<n>();
    ChebyshevBernstein::Coeffs<n> c{};
    for (int k = 0; k <= n; ++k) {
        double sum = 0;
        for (int j = 0; j <= n; ++j)
            sum += M[k][j] * v[static_cast<std::size_t>(j) * (NQ / n)];
        c[k] = sum;
    }
    return c;
}

void fill_powers(double D, double* pw) {
    pw[0] = 1;
    for (int k = 1; k <= Term::MAX_POW; ++k)
        pw[k] = pw[k - 1] * D;
}

bool small_int(double v, int& out) {
    if (v >= 0 && v <= Term::MAX_POW && v == std::floor(v)) {
        out = static_cast<int>(v);
        return true;
    }
    out = -1;
    return false;
}

}  // namespace

// ------------------------------------------------------------------ Term

double Term::phi(double D) const {
    double u = 0;
    if (has_cl) u -= c * std::pow(D, l);
    if (has_e1) u -= e1 * (D - eps1);
    if (has_e2) u -= e2 * (D - eps2) * (D - eps2);
    return std::pow(D, d) * std::exp(u);
}

double Term::chi(double D) const {
    double f = 0, df = 0;
    chi_d(D, nullptr, f, df);
    return f;
}

// chi = delta^d e^u (d + delta u'),  chi' = delta^(d-1) e^u [a^2 + delta u' + delta^2 u''],  a = d + delta u'
void Term::chi_d(double D, const double* pw, double& f, double& df) const {
    const bool use_pw = pw != nullptr && di >= 0 && (!has_cl || li >= 0);
    const double iD = 1 / D;
    double u = 0, du = 0, d2u = 0;
    if (has_cl) {
        const double dl = use_pw ? pw[li] : std::pow(D, l);
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
    const double base = (use_pw ? pw[di] : std::pow(D, d)) * std::exp(u), a = d + D * du;
    f = base * a;
    df = base * iD * (a * a + D * du + D * D * d2u);
}

double Term::parts(double D) const {
    double u = 0, du = 0;
    if (has_cl) {
        const double dl = std::pow(D, l);
        u -= c * dl;
        du += std::abs(c * l * dl / D);
    }
    if (has_e1) {
        u -= e1 * (D - eps1);
        du += std::abs(e1);
    }
    if (has_e2) {
        u -= e2 * (D - eps2) * (D - eps2);
        du += std::abs(2 * e2 * (D - eps2));
    }
    return std::pow(D, d) * std::exp(u) * (std::abs(d) + D * du);
}

double Term::phi_pw(double D, const double* pw) const {
    double u = 0;
    if (has_cl) u -= c * (li >= 0 ? pw[li] : std::pow(D, l));
    if (has_e1) u -= e1 * (D - eps1);
    if (has_e2) u -= e2 * (D - eps2) * (D - eps2);
    return (di >= 0 ? pw[di] : std::pow(D, d)) * std::exp(u);
}

double Term::kappa(double tau, double ltau) const {
    double u = t * ltau;
    if (has_om) u -= om * std::exp(m * ltau);
    if (has_b1) u -= b1 * (tau - g1);
    if (has_b2) u -= b2 * (tau - g2) * (tau - g2);
    return n * std::exp(u);
}

bool Term::same_delta(const Term& o) const {
    return o.d == d && o.has_cl == has_cl && (!has_cl || (o.c == c && o.l == l)) && o.has_e1 == has_e1 && (!has_e1 || (o.e1 == e1 && o.eps1 == eps1))
           && o.has_e2 == has_e2 && (!has_e2 || (o.e2 == e2 && o.eps2 == eps2));
}

// ------------------------------------------------------------------ non-analytic terms

// Delta-derivatives of order 0-2, as ResidualHelmholtzNonAnalytic::all_deltaonly computes them
NonAnalyticValues eval_nonanalytic(const std::vector<NonAnalyticTerm>& terms, double tau_in, double delta_in) {
    NonAnalyticValues v;
    const double tau = std::abs(tau_in - 1) < 10 * DBL_EPSILON ? 1.0 + 10 * DBL_EPSILON : tau_in;
    const double delta = std::abs(delta_in - 1) < 10 * DBL_EPSILON ? 1.0 + 10 * DBL_EPSILON : delta_in;
    const double dm = delta - 1, d2 = dm * dm;
    double ad = 0, add = 0;
    // The terms of one fluid share beta and often a (IAPWS-95: a = 3.5, beta = 0.3 throughout), so the powers of
    // (delta - 1)^2 are computed once per distinct value; d2^(s-1) = d2^s / d2 and Delta^(b-1) = Delta^b / Delta are
    // exact up to one rounding (d2 >= (10 eps)^2 > 0 and Delta > 0 here), unlike separate pow() calls.
    double last_beta = std::numeric_limits<double>::quiet_NaN(), pw_t = 0, pw_t1 = 0;
    double last_a = std::numeric_limits<double>::quiet_NaN(), pa = 0, pa1 = 0;
    for (const auto& el : terms) {
        if (!(el.beta == last_beta)) {
            pw_t = std::pow(d2, 1.0 / (2.0 * el.beta));  // d2^(1/(2 beta))
            pw_t1 = pw_t / d2;                           // d2^(1/(2 beta) - 1)
            last_beta = el.beta;
        }
        if (!(el.a == last_a)) {
            pa = std::pow(d2, el.a);
            pa1 = pa / d2;
            last_a = el.a;
        }
        const double theta = (1.0 - tau) + el.A * pw_t;
        const double dtheta = el.A / el.beta * pw_t1 * dm;
        const double d2theta = el.A / el.beta * (1 / el.beta - 1) * pw_t1;
        const double PSI = std::exp(-el.C * d2 - el.D * (tau - 1.0) * (tau - 1.0));
        const double dPSI = -2.0 * el.C * dm * PSI;
        const double d2PSI = (2.0 * el.C * d2 - 1.0) * 2.0 * el.C * PSI;
        const double DELTA = theta * theta + el.B * pa;
        const double dDELTA = 2 * theta * dtheta + 2 * el.B * el.a * pa1 * dm;
        const double d2DELTA = 2 * (theta * d2theta + dtheta * dtheta + el.B * (2 * el.a * el.a - el.a) * pa1);
        const double Db = std::pow(DELTA, el.b), Db1 = Db / DELTA, Db2 = Db1 / DELTA;
        const double dDb = el.b * Db1 * dDELTA;
        const double d2Db = el.b * (Db1 * d2DELTA + (el.b - 1.0) * Db2 * dDELTA * dDELTA);
        v.alphar += delta * el.n * Db * PSI;
        ad += el.n * (Db * (PSI + delta * dPSI) + dDb * delta * PSI);
        add += el.n * (Db * (2.0 * dPSI + delta * d2PSI) + 2.0 * dDb * (PSI + delta * dPSI) + d2Db * delta * PSI);
        // roundoff scale of delta * (the d(alphar)/d(delta) contribution), in units of eps: the size of its parts, with
        // Delta and dDelta/d(delta) weighted by their own cancellation (theta -> 0 and Delta -> 0 near the critical point)
        const double th_abs = std::abs(1.0 - tau) + std::abs(el.A * pw_t);
        const double eDELTA = 2 * std::abs(theta) * th_abs + std::abs(el.B * pa);  // abs. rounding of Delta / eps
        const double edDELTA = 2 * th_abs * std::abs(dtheta) + std::abs(2 * el.B * el.a * pa1 * dm);
        const double kD = eDELTA / DELTA;
        const double t1 = std::abs(Db) * (std::abs(PSI) + delta * std::abs(dPSI)) * (1 + std::abs(el.b) * kD);
        const double t2 = std::abs(el.b * Db1) * delta * std::abs(PSI) * (std::abs(dDELTA) * (1 + std::abs(el.b - 1) * kD) + edDELTA);
        v.parts += delta * std::abs(el.n) * (t1 + t2);
    }
    v.chi = delta * ad;
    v.dchi = ad + delta * add;
    return v;
}

// |chi| <= delta |n| [Delta^b (psi + delta |psi'|) + |d(Delta^b)/d(delta)| delta psi], with on [lo, hi] (d = |delta - 1|
// in [dmin, dmax]) and |1 - tau| <= dtau_max:
//   theta <= dtau_max + A dmax^(1/beta),  Delta <= theta_max^2 + B dmax^(2a),  |theta'| <= (A/beta) dmax^(1/beta - 1),
//   psi <= exp(-C dmin^2) exp(-D (tau-1)^2),  |psi'| <= 2 C dmax psi,
//   |d(Delta^b)/d(delta)| = b Delta^(b-1) |2 theta theta' + 2 a B d^(2a-1)| is bounded by
//     b >= 1:      b Delta_max^(b-1) (2 theta_max |theta'|_max + 2 a B dmax^(2a-1))
//     1/2 <= b < 1: 2 b Delta_max^(b-1/2) |theta'|_max + 2 a b B^b dmax^(2ab-1)
//       (using |theta| <= Delta^(1/2) and Delta >= B d^(2a) in Delta^(b-1)).
double NonAnalyticTerm::bound_factor(double lo, double hi, double dtau_max) const {
    const double inf = std::numeric_limits<double>::infinity();
    if (!(beta > 0 && beta <= 1 && b >= 0.5 && a >= 0.5 && 2 * a * b >= 1 && A >= 0 && B > 0 && C >= 0 && D >= 0)) return inf;
    if (!(lo >= 0 && hi >= lo && dtau_max >= 0)) return inf;
    const double dmin = (lo <= 1 && hi >= 1) ? 0.0 : std::min(std::abs(lo - 1), std::abs(hi - 1));
    const double dmax = std::max(std::abs(lo - 1), std::abs(hi - 1));
    const double th = dtau_max + A * std::pow(dmax, 1 / beta);
    const double Dl = th * th + B * std::pow(dmax, 2 * a);
    const double dth = A / beta * std::pow(dmax, 1 / beta - 1);
    const double psi = std::exp(-C * dmin * dmin), dpsi = 2 * C * dmax * psi;
    const double dDb = b >= 1 ? b * std::pow(Dl, b - 1) * (2 * th * dth + 2 * a * B * std::pow(dmax, 2 * a - 1))
                              : 2 * b * std::pow(Dl, b - 0.5) * dth + 2 * a * b * std::pow(B, b) * std::pow(dmax, 2 * a * b - 1);
    const double r = hi * std::abs(n) * (std::pow(Dl, b) * (psi + hi * dpsi) + dDb * hi * psi);
    return std::isfinite(r) ? r : inf;
}

void collect_terms(const ResidualHelmholtzGeneralizedExponential& g, int i, int j, std::vector<Term>& out) {
    for (const auto& el : g.elements) {
        Term T;
        T.i = i;
        T.j = j;
        T.n = el.n;
        T.d = el.d;
        T.t = el.t;
        // Mirrors ResidualHelmholtzGeneralizedExponential::all(); build() verifies the regrouped model against
        // the backend's own alphar, so a divergence here is caught there.
        T.has_cl = g.delta_li_in_u && std::isfinite(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON;
        T.c = el.c;
        T.l = el.l_double;
        T.has_om = g.tau_mi_in_u && std::abs(el.m_double) > 0;
        T.om = el.omega;
        T.m = el.m_double;
        T.has_e1 = g.eta1_in_u && std::isfinite(el.eta1);
        T.e1 = el.eta1;
        T.eps1 = el.epsilon1;
        T.has_e2 = g.eta2_in_u && std::isfinite(el.eta2);
        T.e2 = el.eta2;
        T.eps2 = el.epsilon2;
        T.has_b1 = g.beta1_in_u && std::isfinite(el.beta1);
        T.b1 = el.beta1;
        T.g1 = el.gamma1;
        T.has_b2 = g.beta2_in_u && std::isfinite(el.beta2);
        T.b2 = el.beta2;
        T.g2 = el.gamma2;
        small_int(T.d, T.di);
        if (T.has_cl) small_int(T.l, T.li);
        out.push_back(T);
    }
}

// ------------------------------------------------------------------ build

namespace {
constexpr std::size_t MAX_PIECES = 4096;
constexpr int NA_GRADING_LEVELS = 10;  // pieces graded toward delta = 1 down to a width of 2^-10 (see build)
// the |1 - tau| ladder of the non-analytic bounds (see build)
constexpr double NA_DTAU_RUNGS[] = {1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0};
constexpr int NA_TABLE_TAU_LEVELS = 16;  // tau-cells of the 2-D non-analytic tables graded toward tau = 1 (see NATable)
}  // namespace

// ------------------------------------------------------------------ shared 2-D tables of the non-analytic terms

// The non-analytic terms depend on (tau, delta) only -- not on the mixture -- so instead of fitting them on every piece
// at every assemble() (~10^2 us near tau = 1), each fluid's terms are tabulated once per process, on a (tau, delta) grid
// of their own: per distinct D, the smooth part g_D = chi / exp(-D (tau-1)^2) (the Gaussian in tau is steep, D ~
// 300-800, and is multiplied back exactly), as degree-NQ x NQ Lobatto interpolants on tau-cells x delta-cells.
//  - delta-cells: quarters of [0, 8], graded toward delta = 1 like the pieces (1 +- 2^-k), cut back to where the terms
//    are not negligible; a mixture's pieces include these edges, so each piece lies in one cell;
//  - tau-cells: the band where the terms are not negligible, graded toward tau = 1;
//  - a cell keeps its table only if its error, measured on the 16 x 16 grid between its nodes, is negligible against
//    the table tolerance; elsewhere -- where the valley Delta ~ 0 along tau - 1 = A |delta - 1|^(1/beta) crosses the
//    cell, for tau > 1 only -- assemble() fits at the actual tau.
struct NATable
{
    std::vector<double> Dg;                            ///< the distinct D
    std::vector<std::vector<NonAnalyticTerm>> gterms;  ///< the terms with each D, with D set to 0
    std::vector<double> tau_edges;                     ///< tau-cells (empty: the terms are negligible everywhere)
    std::vector<double> delta_edges;                   ///< delta-cells
    std::vector<int> cell;                             ///< [q * n_delta + d]: index of the stored table, -1: fit at assemble, -2: negligible
    std::vector<double> cell_err;                      ///< measured error of the table in chi (incl. the largest Gaussian factor on the cell)
    std::vector<double> cell_parts;                    ///< roundoff scale of chi on the cell
    std::vector<double> coef;                          ///< stored tables, each Dg.size() x (NQ+1) x (NQ+1): [g][m (tau)][k (delta)]
    [[nodiscard]] std::size_t n_delta() const {
        return delta_edges.empty() ? 0 : delta_edges.size() - 1;
    }
};

namespace {

// sum over the terms of their bound on [lo, hi] for |1 - tau| <= dt, with the Gaussian factor taken at |1 - tau| = dte
double na_bound(const std::vector<NonAnalyticTerm>& terms, double lo, double hi, double dt, double dte) {
    double b = 0;
    for (const auto& tm : terms)
        b += tm.bound_factor(lo, hi, dt) * std::exp(-tm.D * dte * dte);
    return b;
}

std::shared_ptr<NATable> build_nonanalytic_table(const std::vector<NonAnalyticTerm>& terms, double tol) {
    auto t = std::make_shared<NATable>();
    const double negligible = 1e-3 * tol;
    for (const auto& tm : terms) {
        auto it = std::find(t->Dg.begin(), t->Dg.end(), tm.D);
        if (it == t->Dg.end()) {
            t->Dg.push_back(tm.D);
            t->gterms.emplace_back();
            it = t->Dg.end() - 1;
        }
        NonAnalyticTerm g = tm;
        g.D = 0;
        t->gterms[static_cast<std::size_t>(it - t->Dg.begin())].push_back(g);
    }
    const std::size_t NGd = t->Dg.size();
    std::vector<double> de;
    for (int m = 0; m <= 32; ++m)
        de.push_back(0.25 * m);
    de.push_back(1.0);
    for (int k = 1; k <= NA_GRADING_LEVELS; ++k) {
        de.push_back(1.0 - std::ldexp(1.0, -k));
        de.push_back(1.0 + std::ldexp(1.0, -k));
    }
    std::sort(de.begin(), de.end());
    de.erase(std::unique(de.begin(), de.end()), de.end());
    // tau band half-width w: the largest |1 - tau| (geometric grid from 4 down) where the bound is not negligible on
    // some cell
    double w = 0;
    for (int it = 0; it < 1200 && w == 0; ++it) {
        const double dt = 4.0 * std::pow(0.97, it);
        for (std::size_t d = 0; d + 1 < de.size() && w == 0; ++d)
            if (!(na_bound(terms, de[d], de[d + 1], dt, dt) < negligible)) w = dt / 0.97;
    }
    if (w == 0) return t;  // negligible everywhere
    // delta range: up to the last cell where the bound over the band is not negligible
    std::size_t last = 1;
    for (std::size_t d = 0; d + 1 < de.size(); ++d)
        if (!(na_bound(terms, de[d], de[d + 1], w, 0.0) < negligible)) last = d + 1;
    de.resize(last + 1);
    t->delta_edges = de;
    const double tlo = std::max(1 - w, 1e-9), thi = 1 + w;
    std::vector<double> te = {tlo, thi, 1.0};
    for (int k = 1; k <= NA_TABLE_TAU_LEVELS; ++k) {
        te.push_back(1 - w * std::ldexp(1.0, -k));
        te.push_back(1 + w * std::ldexp(1.0, -k));
    }
    for (double v : te)
        if (v >= tlo && v <= thi) t->tau_edges.push_back(v);
    std::sort(t->tau_edges.begin(), t->tau_edges.end());
    t->tau_edges.erase(std::unique(t->tau_edges.begin(), t->tau_edges.end()), t->tau_edges.end());
    const int Q = static_cast<int>(t->tau_edges.size()) - 1;
    const int Pd = static_cast<int>(t->n_delta());
    t->cell.assign(static_cast<std::size_t>(Q) * Pd, -2);
    t->cell_err.assign(static_cast<std::size_t>(Q) * Pd, 0.0);
    t->cell_parts.assign(static_cast<std::size_t>(Q) * Pd, 0.0);
    const auto& ul = lobatto_points();
    const auto& um = lobatto_midpoints();
    const auto& M = lobatto_matrix<NQ>();
    constexpr int N1 = NQ + 1;
    std::vector<double> tab(NGd * N1 * N1);
    for (int qc = 0; qc < Q; ++qc) {
        const double tl = t->tau_edges[qc], th = t->tau_edges[qc + 1];
        const double dtmin = (tl <= 1 && th >= 1) ? 0.0 : std::min(std::abs(tl - 1), std::abs(th - 1));
        const double dtmax = std::max(std::abs(tl - 1), std::abs(th - 1));
        for (int d = 0; d < Pd; ++d) {
            const double dl = de[d], dh = de[d + 1];
            if (na_bound(terms, dl, dh, dtmax, dtmin) < negligible) continue;  // -2: negligible on the whole cell
            double err = 0, parts = 0;
            bool finite = true;
            for (std::size_t gi = 0; gi < NGd; ++gi) {
                const double emax = std::exp(-t->Dg[gi] * dtmin * dtmin);  // largest Gaussian factor on the cell
                double F[N1][N1], G[N1][N1];
                for (int i = 0; i < N1; ++i)
                    for (int j = 0; j < N1; ++j) {
                        const auto v = eval_nonanalytic(t->gterms[gi], tl + (th - tl) * (ul[i] + 1) / 2, dl + (dh - dl) * (ul[j] + 1) / 2);
                        F[i][j] = v.chi;
                        finite = finite && std::isfinite(v.chi) && std::isfinite(v.parts);
                        parts = std::max(parts, emax * v.parts);
                    }
                for (int i = 0; i < N1; ++i)  // delta direction
                    for (int k = 0; k < N1; ++k) {
                        double sum = 0;
                        for (int j = 0; j < N1; ++j)
                            sum += M[k][j] * F[i][j];
                        G[i][k] = sum;
                    }
                double* C = &tab[gi * N1 * N1];
                for (int m = 0; m < N1; ++m)  // tau direction
                    for (int k = 0; k < N1; ++k) {
                        double sum = 0;
                        for (int i = 0; i < N1; ++i)
                            sum += M[m][i] * G[i][k];
                        C[m * N1 + k] = sum;
                    }
                double eg = 0;
                for (double ut : um)
                    for (double ud : um) {
                        CoeffsQ dc{};
                        for (int k = 0; k < N1; ++k) {
                            CoeffsQ tc{};
                            for (int m = 0; m < N1; ++m)
                                tc[m] = C[m * N1 + k];
                            dc[k] = ChebyshevBernstein::clenshaw<NQ>(tc, ut);
                        }
                        const auto v = eval_nonanalytic(t->gterms[gi], tl + (th - tl) * (ut + 1) / 2, dl + (dh - dl) * (ud + 1) / 2);
                        const double e = std::abs(ChebyshevBernstein::clenshaw<NQ>(dc, ud) - v.chi);
                        finite = finite && std::isfinite(e) && std::isfinite(v.parts);
                        eg = std::max(eg, e);
                        parts = std::max(parts, emax * v.parts);
                    }
                err += emax * eg;
            }
            if (!finite) return nullptr;
            const std::size_t at = static_cast<std::size_t>(qc) * Pd + d;
            t->cell_err[at] = err;
            t->cell_parts[at] = parts;
            if (10 * err <= negligible) {
                t->cell[at] = static_cast<int>(t->coef.size() / (NGd * N1 * N1));
                t->coef.insert(t->coef.end(), tab.begin(), tab.end());
            } else {
                t->cell[at] = -1;
            }
        }
    }
    return t;
}

// The fluid's table, from a process-wide cache keyed on the exact term coefficients and the tolerance (compared
// exactly, so two fluids share a table only if their terms are identical, and a changed EOS gets a new one).  Entries
// are kept for the life of the process: backends are often short-lived (one per PropsSI call), and a weak entry would
// rebuild the table (~80 ms) each time; there is one entry per distinct (terms, tol), and few fluids have such terms.
// Built outside the lock; a racing duplicate build is only waste.
std::shared_ptr<const NATable> get_nonanalytic_table(const std::vector<NonAnalyticTerm>& terms, double tol) {
    std::vector<double> key = {tol, double(NQ), double(NA_GRADING_LEVELS), double(NA_TABLE_TAU_LEVELS)};
    for (const auto& tm : terms)
        key.insert(key.end(), {tm.n, tm.a, tm.b, tm.beta, tm.A, tm.B, tm.C, tm.D});
    bool cacheable = true;
    for (double v : key)
        cacheable = cacheable && std::isfinite(v);  // NaN would break the map's ordering
    if (!cacheable) return build_nonanalytic_table(terms, tol);
    static std::mutex mtx;
    static std::map<std::vector<double>, std::shared_ptr<const NATable>> cache;
    {
        std::scoped_lock lock(mtx);
        const auto it = cache.find(key);
        if (it != cache.end()) return it->second;
    }
    std::shared_ptr<const NATable> built = build_nonanalytic_table(terms, tol);
    if (!built) return nullptr;
    std::scoped_lock lock(mtx);
    return cache.emplace(key, built).first->second;  // another thread's, if it got there first
}

}  // namespace

std::shared_ptr<const Tables> Tables::build(HelmholtzEOSMixtureBackend& HEOS, const BuildOptions& opt, std::string* reason) {
    if (!(opt.delta_max > 0 && std::isfinite(opt.delta_max)) || !(opt.tau_min >= 0) || !(opt.tau_max > opt.tau_min) || !std::isfinite(opt.tau_max)
        || !(opt.tol > 0) || !(opt.min_width >= 1e-9 * opt.delta_max)) {
        throw std::invalid_argument("ChebDensity::Tables::build: invalid options");
    }
    auto decline = [&](const std::string& why) {
        if (reason) *reason = why;
        return std::shared_ptr<const Tables>();
    };
    auto& comps = HEOS.get_components();
    if (comps.empty()) return decline("no components");
    if (!HEOS.Reducing || !HEOS.residual_helmholtz) return decline("backend has no reducing function or residual model");

    std::shared_ptr<Tables> T(new Tables());
    T->m_opt = opt;
    T->m_N = static_cast<int>(comps.size());
    const int N = T->m_N;
    for (int i = 0; i < N; ++i) {
        if (comps[i].EOSVector.empty()) return decline("component " + std::to_string(i) + " has no equation of state");
        const auto& ar = comps[i].EOS().alphar;
        // Term types that do not factor as kappa(tau) phi(delta) and have no add-in here.  Checked by type, not only by
        // the numerical verification below, which samples a few points and could miss a term that is locally small.
        if (!ar.SAFT.disabled) return decline("component " + std::to_string(i) + " has SAFT association terms");
        if (ar.cubic.enabled || ar.XiangDeiters.enabled || ar.GaoB.enabled)
            return decline("component " + std::to_string(i) + " has residual terms other than generalized exponential");
        collect_terms(ar.GenExp, i, -1, T->m_terms);
        if (ar.NonAnalytic.N > 0) {  // owned copy: no pointer into the backend (COO-125)
            NAComp c;
            c.i = i;
            for (const auto& el : ar.NonAnalytic.elements)
                c.terms.push_back({el.n, el.a, el.b, el.beta, el.A, el.B, el.C, el.D});
            T->m_na.push_back(std::move(c));
        }
        T->m_Ri.push_back(comps[i].gas_constant());
    }
    const auto& Ex = HEOS.residual_helmholtz->Excess;
    T->m_F.assign(N, std::vector<double>(N, 0.0));
    if (N > 1) {
        if (Ex.F.size() != static_cast<std::size_t>(N) || Ex.DepartureFunctionMatrix.size() != static_cast<std::size_t>(N))
            return decline("excess term not sized for the components");
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j) {
                T->m_F[i][j] = Ex.F[i][j];
                if (Ex.F[i][j] == 0) continue;
                if (Ex.DepartureFunctionMatrix[i].size() != static_cast<std::size_t>(N) || !Ex.DepartureFunctionMatrix[i][j])
                    return decline("missing departure function for a pair with F_ij != 0");
                collect_terms(Ex.DepartureFunctionMatrix[i][j]->phi, i, j, T->m_terms);
            }
    }
    T->m_red.reset(HEOS.Reducing->copy());
    // The gas constant as the backend computes it (CoolProp's NORMALIZE_GAS_CONSTANTS setting, or a backend override
    // such as GERG-2008's fixed R): either one constant, or the mole-fraction average of the component values.
    {
        const auto& mf = HEOS.get_mole_fractions_ref();
        if (mf.size() != static_cast<std::size_t>(N)) return decline("mole fractions not set");
        const double Rb = HEOS.calc_gas_constant();
        double Rmix = 0;
        bool all_equal = true;
        for (int i = 0; i < N; ++i) {
            Rmix += mf[i] * T->m_Ri[i];
            all_equal = all_equal && T->m_Ri[i] == T->m_Ri[0];
        }
        if (!(Rb > 0) || !std::isfinite(Rb)) return decline("invalid gas constant");
        T->m_R_mixed = !all_equal && std::abs(Rb - Rmix) <= 1e-14 * Rb;
        T->m_R = Rb;
    }

    // group the terms by their delta-part
    for (const auto& tm : T->m_terms) {
        const auto it = std::find_if(T->m_reps.begin(), T->m_reps.end(), [&](const Term& a) { return tm.same_delta(a); });
        if (it == T->m_reps.end()) {
            T->m_group.push_back(static_cast<int>(T->m_reps.size()));
            T->m_reps.push_back(tm);
        } else {
            T->m_group.push_back(static_cast<int>(it - T->m_reps.begin()));
        }
    }

    const std::string bad = T->verify(HEOS);
    if (!bad.empty()) return decline(bad);
    for (auto& c : T->m_na) {  // the fluids' 2-D tables of their non-analytic terms: shared, built on first use
        c.table = get_nonanalytic_table(c.terms, opt.tol);
        if (!c.table) return decline("non-finite non-analytic table");
    }

    // Largest |W_g| over the tau range (x factors are at most 1, departure terms weighted by |F_ij|).  Used only to
    // steer the piece splitting; the margins in assemble() use the actual weights.
    const std::size_t NGRP = T->m_reps.size();
    std::vector<double> gmax(NGRP, 0.0);
    {
        const double tlo = opt.tau_min > 0 ? opt.tau_min : opt.tau_max / 400;
        for (std::size_t k = 0; k < T->m_terms.size(); ++k) {
            const Term& tm = T->m_terms[k];
            double km = 0;
            bool finite = true;
            for (int s = 0; s <= 400; ++s) {
                const double tau = tlo + (opt.tau_max - tlo) * s / 400.0, v = std::abs(tm.kappa(tau, std::log(tau)));
                finite = finite && std::isfinite(v);  // separately: std::max would drop a NaN
                km = std::max(km, v);
            }
            if (!finite) return decline("tau-dependence not finite on [tau_min, tau_max]");
            gmax[T->m_group[k]] += km * (tm.j < 0 ? 1.0 : std::abs(T->m_F[tm.i][tm.j]));
        }
        for (double v : gmax)
            if (!std::isfinite(v)) return decline("tau-dependence not finite on [tau_min, tau_max]");
    }

    // Roundoff scale of evaluating chi_g on [lo, hi]: the size of the parts it is summed from, at the interpolation
    // nodes.  Not the size of chi_g itself, which passes through zero (e.g. at delta = 1 for every d = l term,
    // chi = delta^d e^{-delta^l} (d - l delta^l)) while the rounding error of its evaluation does not.
    auto parts_max = [](const Term& rep, double lo, double hi) {
        double m = 0;
        bool finite = true;
        for (int j = 0; j <= NQ; ++j) {
            const double v = rep.parts(lo + (hi - lo) * (std::cos(PI * (j + 0.5) / (NQ + 1)) + 1) / 2);
            finite = finite && std::isfinite(v);
            m = std::max(m, v);
        }
        return finite ? m : std::numeric_limits<double>::quiet_NaN();  // any non-finite value -> NaN, checked by the caller
    };

    // Adaptive pieces shared by all groups.  The allowed fit error of a weighted group is tol in units of Z, relaxed
    // where the group is locally large (tol relative to its smaller end value on the piece) and to the roundoff
    // floor of its evaluation (except on the piece touching delta = 0).
    std::vector<std::pair<double, double>> todo = {{0.0, opt.delta_max}}, done;
    while (!todo.empty()) {
        const auto [lo, hi] = todo.back();
        todo.pop_back();
        bool ok = true;
        for (std::size_t g = 0; g < NGRP && ok; ++g) {
            const Term& rep = T->m_reps[g];
            const CoeffsQ c = chebfit<NQ>([&](double D) { return rep.chi(D); }, lo, hi);
            double big = 0;
            for (double v : c)
                big = std::max(big, std::abs(v));
            const double loc = std::min(std::abs(rep.chi(lo > 0 ? lo : DBL_MIN)), std::abs(rep.chi(hi)));
            const double pm = parts_max(rep, lo, hi);
            if (!std::isfinite(pm)) return decline("non-finite residual terms on [0, delta_max]");
            const double floor = lo > 0 ? 4 * DBL_EPSILON * gmax[g] * std::max(big, pm) : 0.0;
            const double allowed = std::max(opt.tol * std::max(1.0, gmax[g] * loc), floor);
            ok = (std::abs(c[NQ]) + std::abs(c[NQ - 1])) * gmax[g] <= allowed;
        }
        const double mid = 0.5 * (lo + hi);
        if (done.size() + todo.size() >= MAX_PIECES) return decline("more than " + std::to_string(MAX_PIECES) + " delta-pieces (tol unattainable?)");
        if (ok || hi - lo < opt.min_width || !(lo < mid && mid < hi)) {
            done.emplace_back(lo, hi);
        } else {
            todo.emplace_back(mid, hi);
            todo.emplace_back(lo, mid);
        }
    }
    std::sort(done.begin(), done.end());
    T->m_edges = {0.0};
    for (const auto& pr : done)
        T->m_edges.push_back(pr.second);
    if (!T->m_na.empty()) {
        // The non-analytic terms are not analytic at delta = 1 (theta carries |delta - 1|^(1/beta), and Delta^b sits on
        // top of it), so Chebyshev fits on a piece ending there converge only algebraically.  Make delta = 1 an edge and
        // grade the pieces geometrically toward it, 1 +- 2^-k: each graded piece is then at a distance from the
        // singular point comparable to its width (geometric convergence), and on the last ones the terms are small.
        // (These edges are added after the MAX_PIECES check: the grading and the tables' cell edges, a few dozen.)
        std::vector<double> extra = {1.0};
        for (int k = 1; k <= NA_GRADING_LEVELS; ++k) {
            extra.push_back(1.0 - std::ldexp(1.0, -k));
            extra.push_back(1.0 + std::ldexp(1.0, -k));
        }
        for (const auto& c : T->m_na)  // and the edges of the 2-D tables' delta-cells, so no piece straddles two cells
            extra.insert(extra.end(), c.table->delta_edges.begin(), c.table->delta_edges.end());
        for (double e : extra)
            if (e > 0 && e < opt.delta_max) T->m_edges.push_back(e);
        std::sort(T->m_edges.begin(), T->m_edges.end());
        T->m_edges.erase(std::unique(T->m_edges.begin(), T->m_edges.end()), T->m_edges.end());
    }
    const int P = static_cast<int>(T->m_edges.size()) - 1;

    T->m_C.resize(NGRP * P);
    T->m_Cn.resize(NGRP * P);
    T->m_Ct.resize(NGRP * P);
    T->m_Cp.resize(NGRP * P);
    for (std::size_t g = 0; g < NGRP; ++g) {
        const Term& rep = T->m_reps[g];
        for (int p = 0; p < P; ++p) {
            const double lo = T->m_edges[p], hi = T->m_edges[p + 1];
            const CoeffsQ c = chebfit<NQ>([&](double D) { return rep.chi(D); }, lo, hi);
            double s = 0;
            for (double v : c)
                s += std::abs(v);
            // Fit-error bound: the error measured at 64 points between the interpolation nodes, never below the
            // tail estimate (the tail alone underestimates it on pieces where the series has not settled)
            double meas = 0, pmax = parts_max(rep, lo, hi);
            bool finite = std::isfinite(pmax);
            for (int k = 0; k < 64; ++k) {
                const double u = -1 + 2 * (k + 0.5) / 64, D = lo + (hi - lo) * (u + 1) / 2;
                const double e = std::abs(ChebyshevBernstein::clenshaw<NQ>(c, u) - rep.chi(D)), pv = rep.parts(D);
                finite = finite && std::isfinite(e) && std::isfinite(pv);  // not via max(): it would drop a NaN
                meas = std::max(meas, e);
                pmax = std::max(pmax, pv);
            }
            T->m_C[p * NGRP + g] = c;
            T->m_Cn[p * NGRP + g] = s;
            T->m_Ct[p * NGRP + g] = std::max(meas, std::abs(c[NQ]) + std::abs(c[NQ - 1]));
            T->m_Cp[p * NGRP + g] = pmax;
            if (!finite || !std::isfinite(s)) return decline("non-finite fit on [0, delta_max]");
        }
    }
    // Bound on each non-analytic term per piece, without its factor exp(-D (tau-1)^2), at a ladder of |1 - tau|: the
    // bound grows monotonically with |1 - tau| (through theta_max), so assemble() takes the first rung at or above the
    // actual value -- much tighter near the critical point than the bound over the whole tau range.
    {
        const double dtau_max = std::max(std::abs(1 - opt.tau_min), std::abs(opt.tau_max - 1));
        for (double r : NA_DTAU_RUNGS)
            if (r < dtau_max) T->m_na_dtau.push_back(r);
        T->m_na_dtau.push_back(dtau_max);
        const std::size_t R = T->m_na_dtau.size();
        for (auto& c : T->m_na) {
            const std::size_t nk = c.terms.size();
            c.K.resize(static_cast<std::size_t>(P) * nk * R);
            for (int p = 0; p < P; ++p)
                for (std::size_t k = 0; k < nk; ++k)
                    for (std::size_t r = 0; r < R; ++r)
                        c.K[(p * nk + k) * R + r] = c.terms[k].bound_factor(T->m_edges[p], T->m_edges[p + 1], T->m_na_dtau[r]);
        }
    }
    // where each piece sits in its fluid's 2-D table: the cell, and the re-expansion onto the piece when it is only part
    // of the cell (the pieces' edges include the cells' edges, so a piece never straddles two cells)
    for (auto& c : T->m_na) {
        constexpr int N1 = NQ + 1;
        const auto& de = c.table->delta_edges;
        const auto& ul = lobatto_points();
        const auto& Ml = lobatto_matrix<NQ>();
        c.table_cell.assign(P, -1);
        c.table_identity.assign(P, 0);
        c.table_reexpand.assign(static_cast<std::size_t>(P) * N1 * N1, 0.0);
        for (int p = 0; p < P; ++p) {
            const double lo = T->m_edges[p], hi = T->m_edges[p + 1];
            const int pd = static_cast<int>(std::upper_bound(de.begin(), de.end(), lo) - de.begin()) - 1;
            if (pd < 0 || pd + 1 >= static_cast<int>(de.size()) || hi > de[pd + 1]) continue;  // outside the table
            c.table_cell[p] = pd;
            if (lo == de[pd] && hi == de[pd + 1]) {
                c.table_identity[p] = 1;
                continue;
            }
            // piece coefficients = M (V cell coefficients), V[j][k] = T_k(cell coordinate of the piece's Lobatto node j)
            double V[N1][N1];
            for (int j = 0; j < N1; ++j) {
                const double x = lo + (hi - lo) * (ul[j] + 1) / 2, uc = 2 * (x - de[pd]) / (de[pd + 1] - de[pd]) - 1;
                V[j][0] = 1;
                V[j][1] = uc;
                for (int k = 2; k < N1; ++k)
                    V[j][k] = 2 * uc * V[j][k - 1] - V[j][k - 2];
            }
            double* Rp = &c.table_reexpand[static_cast<std::size_t>(p) * N1 * N1];
            for (int a = 0; a < N1; ++a)
                for (int k = 0; k < N1; ++k) {
                    double sum = 0;
                    for (int j = 0; j < N1; ++j)
                        sum += Ml[a][j] * V[j][k];
                    Rp[a * N1 + k] = sum;
                }
        }
    }
    return T;
}

// The regrouped model must reproduce the backend's alphar and delta d(alphar)/d(delta): every pure component, and the
// equimolar mixture, at points spread over the rectangle.  Catches term types or backends (e.g. cubic) whose
// residual is not the sum collected here.
std::string Tables::verify(HelmholtzEOSMixtureBackend& HEOS) const {
    std::vector<std::vector<double>> xs;
    for (int i = 0; i < m_N; ++i) {
        std::vector<double> x(m_N, 0.0);
        x[i] = 1;
        xs.push_back(x);
    }
    if (m_N > 1) xs.emplace_back(m_N, 1.0 / m_N);
    const double tlo = m_opt.tau_min > 0 ? m_opt.tau_min : 0.05 * m_opt.tau_max;
    std::vector<std::pair<double, double>> pts;  // (tau, delta)
    for (double ft : {0.13, 0.5, 0.97})
        for (double fd : {0.011, 0.23, 0.61, 0.97})
            pts.emplace_back(tlo + ft * (m_opt.tau_max - tlo), fd * m_opt.delta_max);
    // the non-analytic terms are negligible except near tau = 1 = delta (exp(-D (tau-1)^2), D ~ 300-800): sample there
    // too, or a wrong evaluation of them would pass
    if (!m_na.empty())
        for (double tau : {0.98, 0.999, 1.0, 1.001, 1.02})
            for (double D : {0.9, 0.99, 1.01, 1.1})
                if (D <= m_opt.delta_max) pts.emplace_back(std::min(std::max(tau, tlo), m_opt.tau_max), D);  // clamped, not dropped
    for (const auto& x : xs)
        for (const auto& [tau, D] : pts) {
            const double lt = std::log(tau);
            // alphar, chi = delta alphar_delta and its delta-derivative (used by Newton steps on the true equation)
            double a = 0, z = 0, dz = 0, sa = 0, sz = 0, sdz = 0;
            double pw[Term::MAX_POW + 1];
            fill_powers(D, pw);
            for (const Term& tm : m_terms) {
                const double w = (tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * m_F[tm.i][tm.j]) * tm.kappa(tau, lt);
                double f = 0, df = 0;
                tm.chi_d(D, pw, f, df);
                const double va = w * tm.phi(D);
                a += va;
                z += w * f;
                dz += w * df;
                sa += std::abs(va);
                sz += std::abs(w * f);
                sdz += std::abs(w * df);
            }
            for (const auto& c : m_na) {
                const auto v = eval_nonanalytic(c.terms, tau, D);
                a += x[c.i] * v.alphar;
                z += x[c.i] * v.chi;
                dz += x[c.i] * v.dchi;
                sa += std::abs(x[c.i] * v.alphar);
                sz += std::abs(x[c.i] * v.chi);
                sdz += std::abs(x[c.i] * v.dchi);
            }
            double a_ref = 0, z_ref = 0, dz_ref = 0;
            try {
                a_ref = HEOS.calc_alphar_deriv_nocache(0, 0, x, tau, D);
                const double a1 = HEOS.calc_alphar_deriv_nocache(0, 1, x, tau, D), a2 = HEOS.calc_alphar_deriv_nocache(0, 2, x, tau, D);
                z_ref = D * a1;
                dz_ref = a1 + D * a2;
            } catch (const std::exception& e) {
                return std::string("backend alphar failed during verification: ") + e.what();
            }
            if (!(std::abs(a - a_ref) <= 1e-12 * (1 + sa)) || !(std::abs(z - z_ref) <= 1e-12 * (1 + sz))
                || !(std::abs(dz - dz_ref) <= 1e-12 * (1 + sdz)))
                return "regrouped residual does not reproduce the backend's alphar (unsupported term type)";
        }
    return "";
}

// ------------------------------------------------------------------ assemble

bool Tables::assemble(double T, const std::vector<double>& x, State& S) const {
    if (x.size() != static_cast<std::size_t>(m_N) || !(T > 0) || !std::isfinite(T)) return false;
    for (double v : x)
        if (!std::isfinite(v)) return false;
    const double Tr = m_red->Tr(x), rhor = m_red->rhormolar(x), tau = Tr / T;
    if (!(tau >= m_opt.tau_min && tau <= m_opt.tau_max && tau > 0) || !(rhor > 0) || !std::isfinite(rhor)) return false;
    const double lt = std::log(tau);
    const std::size_t NGRP = m_reps.size();
    S.W.assign(NGRP, 0.0);
    for (std::size_t k = 0; k < m_terms.size(); ++k) {
        const Term& tm = m_terms[k];
        const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * m_F[tm.i][tm.j];
        S.W[m_group[k]] += X * tm.kappa(tau, lt);
    }
    const int P = n_pieces();
    const double rnd = (45.0 + static_cast<double>(NGRP)) * DBL_EPSILON;
    S.G.resize(P);
    S.margin.resize(P);
    S.na_fits = 0;
    // the rung of the |1 - tau| ladder of the bounds: the first at or above |1 - tau|
    const std::size_t R = m_na_dtau.size();
    std::size_t rung = 0;
    while (rung + 1 < R && m_na_dtau[rung] < std::abs(1 - tau))
        ++rung;
    // tau-factor exp(-D (tau-1)^2) of each non-analytic term's bound, shared by all pieces
    std::vector<double> na_tau_factor;
    for (const NAComp& c : m_na)
        for (const auto& tm : c.terms)
            na_tau_factor.push_back(std::exp(-tm.D * (tau - 1) * (tau - 1)));
    // per non-analytic component: its tau-cell, the Chebyshev polynomials in the cell's variable, the Gaussian factors
    struct NATau
    {
        int q = -1;
        std::array<double, NQ + 1> Tm{};
        std::vector<double> ef;
    };
    std::vector<NATau> nat(m_na.size());
    for (std::size_t ci = 0; ci < m_na.size(); ++ci) {
        const NAComp& c = m_na[ci];
        if (!c.table) continue;
        const auto& te = c.table->tau_edges;
        if (te.size() < 2 || tau < te.front() || tau > te.back()) continue;
        int qc = static_cast<int>(std::upper_bound(te.begin(), te.end(), tau) - te.begin()) - 1;
        qc = std::min(std::max(qc, 0), static_cast<int>(te.size()) - 2);
        const double ut = std::min(1.0, std::max(-1.0, 2 * (tau - te[qc]) / (te[qc + 1] - te[qc]) - 1));
        auto& nt = nat[ci];
        nt.q = qc;
        nt.Tm[0] = 1;
        nt.Tm[1] = ut;
        for (int m = 2; m <= NQ; ++m)
            nt.Tm[m] = 2 * ut * nt.Tm[m - 1] - nt.Tm[m - 2];
        for (double Dv : c.table->Dg)
            nt.ef.push_back(std::exp(-Dv * (tau - 1) * (tau - 1)));
    }
    S.na_table = 0;
    for (int p = 0; p < P; ++p) {
        std::size_t ntf = 0, nci = 0;

        CoeffsQ q{};
        q[0] = 1.0;
        double scale = 1.0, fiterr = 0.0;
        // roundoff allowance per unit of scale: an empirical ~45 eps for evaluating one chi_g (the parts bound omits
        // the exp/pow rounding, (|u| + d + l) eps relative), plus one eps per group summed
        const CoeffsQ* cp = &m_C[p * NGRP];
        const double* cn = &m_Cn[p * NGRP];
        const double* ct = &m_Ct[p * NGRP];
        const double* cpart = &m_Cp[p * NGRP];
        for (std::size_t g = 0; g < NGRP; ++g) {
            const double w = S.W[g];
            for (int k = 0; k <= NQ; ++k)
                q[k] += w * cp[g][k];
            scale += std::abs(w) * (cn[g] + cpart[g]);
            fiterr += std::abs(w) * ct[g];
        }
        for (const NAComp& c : m_na) {  // non-analytic add-in at this tau
            const double xi = x[c.i];
            const std::size_t nk = c.terms.size();
            const double* tf = &na_tau_factor[ntf];
            ntf += nk;
            const NATau& nt = nat[nci++];
            if (xi == 0) continue;
            double bound = 0;
            for (std::size_t k = 0; k < nk; ++k)
                bound += c.K[(p * nk + k) * R + rung] * tf[k];
            bound *= std::abs(xi);
            if (bound < 1e-3 * m_opt.tol) {  // negligible on this piece (rigorous bound): into the margin, not the fit
                fiterr += bound;
                continue;
            }
            if (nt.q >= 0 && c.table_cell[p] >= 0) {  // from the 2-D table, if this cell has one
                const NATable& tb = *c.table;
                const std::size_t at = static_cast<std::size_t>(nt.q) * tb.n_delta() + c.table_cell[p];
                const int idx = tb.cell[at];
                if (idx >= 0) {
                    constexpr int N1 = NQ + 1;
                    const std::size_t NGd = tb.Dg.size();
                    const double* C = &tb.coef[static_cast<std::size_t>(idx) * NGd * N1 * N1];
                    CoeffsQ cc{};  // on the cell
                    for (std::size_t gi = 0; gi < NGd; ++gi) {
                        const double fac = xi * nt.ef[gi];
                        const double* Cg = C + gi * N1 * N1;
                        for (int m = 0; m < N1; ++m) {
                            const double tm = fac * nt.Tm[m];
                            for (int k = 0; k < N1; ++k)
                                cc[k] += tm * Cg[m * N1 + k];
                        }
                    }
                    CoeffsQ cf = cc;  // on the piece
                    if (c.table_identity[p] == 0) {
                        const double* Rp = &c.table_reexpand[static_cast<std::size_t>(p) * N1 * N1];
                        for (int a = 0; a < N1; ++a) {
                            double sum = 0;
                            for (int k = 0; k < N1; ++k)
                                sum += Rp[a * N1 + k] * cc[k];
                            cf[a] = sum;
                        }
                    }
                    double s1 = 0;
                    for (int k = 0; k <= NQ; ++k) {
                        q[k] += cf[k];
                        s1 += std::abs(cf[k]) + (c.table_identity[p] == 0 ? std::abs(cc[k]) : 0.0);  // re-expansion roundoff
                    }
                    // the cell's error and roundoff scale bound the piece (part of the cell)
                    scale += s1 + std::abs(xi) * tb.cell_parts[at];
                    fiterr += 10 * std::abs(xi) * tb.cell_err[at];
                    ++S.na_table;
                    continue;
                }
            }
            ++S.na_fits;
            const double lo = m_edges[p], hi = m_edges[p + 1];
            bool finite = true;
            double pmax = 0;
            auto f = [&](double D) {
                const auto v = eval_nonanalytic(c.terms, tau, D);
                finite = finite && std::isfinite(v.chi) && std::isfinite(v.parts);
                pmax = std::max(pmax, std::abs(xi) * v.parts);
                return xi * v.chi;
            };
            // Fit at this tau: degree-NQ Chebyshev-Lobatto interpolation (the ends are nodes), its error MEASURED at the
            // NQ points between the nodes (a tail estimate alone fails: on wide pieces the slowly decaying CO2 terms,
            // C = 10, beat it by ~100x).  Only the cells the valley crosses get here, and they need the full degree.
            const auto& ul = lobatto_points();
            static_assert(std::tuple_size_v<std::decay_t<decltype(ul)>> == NQ + 1, "one value per Lobatto point");
            std::array<double, NQ + 1> v{};
            for (std::size_t k = 0; k < v.size(); ++k)
                v[k] = f(lo + (hi - lo) * (ul[k] + 1) / 2);
            const CoeffsQ cf = lobatto_fit<NQ>(v);
            const double tail = std::abs(cf[NQ]) + std::abs(cf[NQ - 1]);
            double meas = 0;
            for (double u : lobatto_midpoints()) {
                const double e = std::abs(ChebyshevBernstein::clenshaw<NQ>(cf, u) - f(lo + (hi - lo) * (u + 1) / 2));
                finite = finite && std::isfinite(e);
                meas = std::max(meas, e);
            }
            for (double e : v)
                finite = finite && std::isfinite(e);
            if (!finite) return false;
            double s1 = 0;
            for (int k = 0; k <= NQ; ++k) {
                q[k] += cf[k];
                s1 += std::abs(cf[k]);
            }
            scale += s1 + pmax;
            fiterr += std::max(10 * meas, 5 * tail);
        }
        S.G[p] = times_delta(q, m_edges[p], m_edges[p + 1]);
        // |G_tables - G_true| <= delta * |Z_tables - Z_true| <= hi * (fit error + roundoff), with a factor 2 on
        // the (measured, not proven) fit error.  The roundoff term covers the evaluation of the fits and of the
        // chi_g themselves (scale includes the size of their parts), which the measured error only samples.
        S.margin[p] = m_edges[p + 1] * (2 * fiterr + rnd * scale);
        if (!std::isfinite(S.margin[p])) return false;
    }
    double Rmix = m_R;
    if (m_R_mixed) {
        Rmix = 0;
        for (int i = 0; i < m_N; ++i)
            Rmix += x[i] * m_Ri[i];
    }
    S.x = x;
    S.T = T;
    S.tau = tau;
    S.rhor = rhor;
    S.t_scale = 1.0 / (rhor * Rmix * T);
    return std::isfinite(S.t_scale) && S.t_scale > 0;
}

// ------------------------------------------------------------------ true equation

void Tables::true_G(const State& S, double D, double t, double& G, double& dG, double& scale) const {
    if (D <= 0) {  // the chi_g vanish like delta^d at 0: G(0) = -t (direct evaluation is 0/0)
        G = -t;
        dG = 1;
        scale = std::abs(t);
        return;
    }
    double pw[Term::MAX_POW + 1];
    fill_powers(D, pw);
    double Z = 1, dZ = 0, a = 0;
    for (std::size_t g = 0; g < m_reps.size(); ++g) {
        double f, df;
        m_reps[g].chi_d(D, pw, f, df);
        Z += S.W[g] * f;
        dZ += S.W[g] * df;
        a += std::abs(S.W[g] * f);
    }
    for (const NAComp& c : m_na) {
        const double xi = S.x[c.i];
        if (xi == 0) continue;
        const auto v = eval_nonanalytic(c.terms, S.tau, D);
        Z += xi * v.chi;
        dZ += xi * v.dchi;
        a += std::abs(xi) * v.parts;
    }
    G = D * Z - t;
    dG = Z + D * dZ;
    scale = D * (1 + a) + std::abs(t);
}

// ------------------------------------------------------------------ the stable root at given p

namespace {
struct Candidate
{
    double D = 0;                ///< table root in delta
    double Da = 0, Db = 0;       ///< its certified bracket
    bool negative_at_a = false;  ///< sign of G at Da
};
}  // namespace

double Tables::stable_root(const State& S, double p) const {
    if (!(p > 0) || !std::isfinite(p)) return -1;
    const double t = p * S.t_scale;
    const int P = n_pieces();
    // a root beyond delta_max: G(delta_max) <= 0 means the isotherm has not reached p yet
    {
        double G, dG, sc;
        true_G(S, m_opt.delta_max, t, G, dG, sc);
        if (!(G > 0)) return -1;
    }
    // all roots
    thread_local std::vector<ChebyshevBernstein::Root> rr;
    thread_local std::vector<Candidate> roots;
    roots.clear();
    for (int pc = 0; pc < P; ++pc) {
        CoeffsG g = S.G[pc];
        g[0] -= t;
        if (std::abs(g[0]) > ChebyshevBernstein::l1_tail<NG>(g) + S.margin[pc]) continue;  // no root on this piece
        ChebyshevBernstein::real_roots<NG>(g, S.margin[pc], rr);
        for (const auto& r : rr) {
            if (!r.certified) return -1;  // cannot vouch for the root set here
            const double D = delta_of(pc, r.u);
            if (!roots.empty() && std::abs(D - roots.back().D) < 1e-10 * roots.back().D) continue;  // shared piece edge
            roots.push_back({D, delta_of(pc, r.ua), delta_of(pc, r.ub), ChebyshevBernstein::clenshaw<NG>(g, r.ua) < 0});
        }
    }
    if (roots.empty()) return -1;
    // selection by spinodal branches
    int k = 0;
    if (roots.size() > 1) {
        // the extrema of delta Z: the first local maximum and the last local minimum, classified by the TRUE slope on
        // either side (the table's derivative root can sit well off the true one where delta Z is flat, so the probe
        // widens geometrically within the piece until the two sides differ)
        double dmax1 = 1e300, dminL = -1;
        thread_local std::vector<ChebyshevBernstein::Root> ro;
        for (int pc = 0; pc < P; ++pc) {
            const CoeffsG d = ChebyshevBernstein::derivative<NG>(S.G[pc]);
            double sa = 0;
            for (double v : d)
                sa += std::abs(v);
            if (std::abs(d[0]) > ChebyshevBernstein::l1_tail<NG>(d) + 1e-13 * sa) continue;
            ChebyshevBernstein::real_roots<NG>(d, 1e-13 * sa, ro);
            const double hmax = 0.5 * (m_edges[pc + 1] - m_edges[pc]);
            for (const auto& q : ro) {
                if (!q.certified) return -1;  // an extremum we cannot place: the branch limits are unknown
                const double D = delta_of(pc, q.u);
                bool classified = false;
                for (int s = 0; s < 40 && !classified; ++s) {  // h = 1e-7 D 8^s, within the piece
                    const double h = 1e-7 * D * std::pow(8.0, s);
                    if (h > hmax) break;
                    double G, dl, dr, sc;
                    true_G(S, D - h, 0, G, dl, sc);
                    true_G(S, D + h, 0, G, dr, sc);
                    if (dl > 0 && dr < 0) {
                        dmax1 = std::min(dmax1, D);
                        classified = true;
                    } else if (dl < 0 && dr > 0) {
                        dminL = std::max(dminL, D);
                        classified = true;
                    }
                }
                // Neither a maximum nor a minimum within the piece: an inflection where delta Z is flat, or a pair of
                // extrema closer than the probe.  Dropping it could move the branch limits past an interior root, so
                // defer rather than guess.
                if (!classified) return -1;
            }
        }
        int vap = -1, liq = -1;
        for (int j = 0; j < static_cast<int>(roots.size()); ++j) {
            double G, dG, sc;
            true_G(S, roots[j].D, t, G, dG, sc);
            if (!(dG > 0)) continue;  // mechanically unstable
            if (roots[j].D < dmax1) vap = j;
            if (roots[j].D > dminL && liq < 0) liq = j;
        }
        if (vap < 0 && liq < 0) return -1;
        if (vap < 0)
            k = liq;
        else if (liq < 0 || liq == vap)
            k = vap;
        else {
            auto g = [&](int j) { return std::log(roots[j].D) + alphar(S, roots[j].D) + t / roots[j].D; };
            k = g(vap) <= g(liq) ? vap : liq;
        }
    }
    // polish: bracketed Newton on the true equation; the certified bracket's end signs are known
    const Candidate& c = roots[k];
    double D = c.D, a = c.Da, b = c.Db, G = 0, dG = 0, sc = 0;
    for (int it = 0; it < 30; ++it) {
        true_G(S, D, t, G, dG, sc);
        if (G == 0) break;
        ((G < 0) == c.negative_at_a ? a : b) = D;
        const double step = G / dG;
        if (std::abs(step) <= 2 * DBL_EPSILON * D) {
            D -= step;
            break;
        }
        double Dn = D - step;
        if (!(Dn > std::min(a, b) && Dn < std::max(a, b))) Dn = 0.5 * (a + b);
        D = Dn;
    }
    if (!(D > 0) || !std::isfinite(D)) return -1;
    true_G(S, D, t, G, dG, sc);
    if (!(dG > 0) || !(std::abs(G) <= 1e-10 * sc)) return -1;  // mechanically stable, and actually a root
    return D * S.rhor;
}

const void* Tables::nonanalytic_table(std::size_t k) const {
    return k < m_na.size() ? m_na[k].table.get() : nullptr;
}

double Tables::alphar(const State& S, double D) const {
    double pw[Term::MAX_POW + 1];
    fill_powers(D, pw);
    double a = 0;
    for (std::size_t g = 0; g < m_reps.size(); ++g)
        a += S.W[g] * m_reps[g].phi_pw(D, pw);
    for (const NAComp& c : m_na)
        if (S.x[c.i] != 0) a += S.x[c.i] * eval_nonanalytic(c.terms, S.tau, D).alphar;
    return a;
}

}  // namespace ChebDensity
}  // namespace CoolProp
