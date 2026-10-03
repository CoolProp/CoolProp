#include "ChebDensitySolver.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
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
}

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
        // Term types that do not factor as kappa(tau) phi(delta).  The non-analytic terms are checked by type, not
        // only by the numerical verification below: they are negligible away from the critical point, so sample
        // points could miss them.
        if (ar.NonAnalytic.N > 0) return decline("component " + std::to_string(i) + " has non-analytic residual terms");
        if (!ar.SAFT.disabled) return decline("component " + std::to_string(i) + " has SAFT association terms");
        if (ar.cubic.enabled || ar.XiangDeiters.enabled || ar.GaoB.enabled)
            return decline("component " + std::to_string(i) + " has residual terms other than generalized exponential");
        collect_terms(ar.GenExp, i, -1, T->m_terms);
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
    const int P = static_cast<int>(done.size());

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
    for (const auto& x : xs)
        for (double ft : {0.13, 0.5, 0.97})
            for (double fd : {0.011, 0.23, 0.61, 0.97}) {
                const double tau = tlo + ft * (m_opt.tau_max - tlo), D = fd * m_opt.delta_max, lt = std::log(tau);
                double a = 0, z = 0, sa = 0, sz = 0;
                for (const Term& tm : m_terms) {
                    const double w = (tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * m_F[tm.i][tm.j]) * tm.kappa(tau, lt);
                    const double va = w * tm.phi(D), vz = w * tm.chi(D);
                    a += va;
                    z += vz;
                    sa += std::abs(va);
                    sz += std::abs(vz);
                }
                double a_ref = 0, z_ref = 0;
                try {
                    a_ref = HEOS.calc_alphar_deriv_nocache(0, 0, x, tau, D);
                    z_ref = D * HEOS.calc_alphar_deriv_nocache(0, 1, x, tau, D);
                } catch (const std::exception& e) {
                    return std::string("backend alphar failed during verification: ") + e.what();
                }
                if (!(std::abs(a - a_ref) <= 1e-12 * (1 + sa)) || !(std::abs(z - z_ref) <= 1e-12 * (1 + sz)))
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
    for (int p = 0; p < P; ++p) {
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
    G = D * Z - t;
    dG = Z + D * dZ;
    scale = D * (1 + a) + std::abs(t);
}

double Tables::alphar(const State& S, double D) const {
    double pw[Term::MAX_POW + 1];
    fill_powers(D, pw);
    double a = 0;
    for (std::size_t g = 0; g < m_reps.size(); ++g)
        a += S.W[g] * m_reps[g].phi_pw(D, pw);
    return a;
}

}  // namespace ChebDensity
}  // namespace CoolProp
