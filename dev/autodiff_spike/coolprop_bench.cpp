// Experiment 3: CoolProp's hand-coded multiparameter derivatives vs a value-only
// evaluation of the same terms, and vs two AD formulations of the same terms:
//   naive      : one bivariate Taylor2<4> pass through every term (generic AD)
//   separable  : per term, univariate Taylor1<4> in delta and in tau, then the outer
//                product n F_j(delta) G_i(tau) -- AD that respects the term structure
// Only the ResidualHelmholtzGeneralizedExponential block is timed (power, exponential,
// Gaussian terms); that is the whole residual for the fluids chosen.
//
// Build: see run_coolprop.py (links against a Release libCoolProp.a).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "taylor1.hpp"

using CoolProp::HelmholtzDerivatives;
using CoolProp::ResidualHelmholtzGeneralizedExponential;
constexpr int NO = 4;
constexpr int NS = (NO + 1) * (NO + 2) / 2;
constexpr int id(int i, int j) {
    return (i + j) * (i + j + 1) / 2 + j;
}  // i: tau order, j: delta order

static bool valid(double x) {
    return std::isfinite(x);
}

// -------- value only, mirroring the u-construction of GenExp::all()
double value_only(const ResidualHelmholtzGeneralizedExponential& g, double tau, double delta) {
    const double lt = std::log(tau), ld = std::log(delta);
    double s = 0;
    for (const auto& el : g.elements) {
        double u = 0;
        if (g.delta_li_in_u && valid(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON) {
            double dl;
            if (el.l_is_int) {
                dl = 1;
                for (int k = 0; k < el.l_int; ++k)
                    dl *= delta;
            } else {
                dl = std::exp(el.l_double * ld);
            }
            u -= el.c * dl;
        }
        if (g.tau_mi_in_u && std::abs(el.m_double) > 0) u -= el.omega * std::exp(el.m_double * lt);
        if (g.eta1_in_u && valid(el.eta1)) u -= el.eta1 * (delta - el.epsilon1);
        if (g.eta2_in_u && valid(el.eta2)) u -= el.eta2 * (delta - el.epsilon2) * (delta - el.epsilon2);
        if (g.beta1_in_u && valid(el.beta1)) u -= el.beta1 * (tau - el.gamma1);
        if (g.beta2_in_u && valid(el.beta2)) u -= el.beta2 * (tau - el.gamma2) * (tau - el.gamma2);
        s += el.n * std::exp(el.t * lt + el.d * ld + u);
    }
    return s;
}

// -------- hand: CoolProp all(), converted to A_ij = tau^i delta^j d^(i+j) alphar
void hand(ResidualHelmholtzGeneralizedExponential& g, double tau, double delta, double* A) {
    HelmholtzDerivatives d;
    g.all(tau, delta, d);
    const double t = tau, e = delta;
    A[id(0, 0)] = d.alphar;
    A[id(1, 0)] = t * d.dalphar_dtau;
    A[id(0, 1)] = e * d.dalphar_ddelta;
    A[id(2, 0)] = t * t * d.d2alphar_dtau2;
    A[id(1, 1)] = t * e * d.d2alphar_ddelta_dtau;
    A[id(0, 2)] = e * e * d.d2alphar_ddelta2;
    A[id(3, 0)] = t * t * t * d.d3alphar_dtau3;
    A[id(2, 1)] = t * t * e * d.d3alphar_ddelta_dtau2;
    A[id(1, 2)] = t * e * e * d.d3alphar_ddelta2_dtau;
    A[id(0, 3)] = e * e * e * d.d3alphar_ddelta3;
    A[id(4, 0)] = t * t * t * t * d.d4alphar_dtau4;
    A[id(3, 1)] = t * t * t * e * d.d4alphar_ddelta_dtau3;
    A[id(2, 2)] = t * t * e * e * d.d4alphar_ddelta2_dtau2;
    A[id(1, 3)] = t * e * e * e * d.d4alphar_ddelta3_dtau;
    A[id(0, 4)] = e * e * e * e * d.d4alphar_ddelta4;
}
void hand_deltaonly(ResidualHelmholtzGeneralizedExponential& g, double tau, double delta, double* A) {
    HelmholtzDerivatives d;
    g.all_deltaonly(tau, delta, d);
    A[0] = d.alphar;
    A[1] = d.d4alphar_ddelta4;
}

// -------- naive generic AD: the u-construction above, verbatim, in Taylor2<4>
void naive(const ResidualHelmholtzGeneralizedExponential& g, double tau, double delta, double* A) {
    using B = nd::Taylor2<NO>;
    B T(tau), D(delta);
    T.c[B::idx(1, 0)] = tau;  // tau = tau0 (1 + x), delta = delta0 (1 + y)
    D.c[B::idx(0, 1)] = delta;
    const B lt = log(T), ld = log(D);
    B s(0.0);
    for (const auto& el : g.elements) {
        B u(0.0);
        if (g.delta_li_in_u && valid(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON) u -= exp(ld * el.l_double) * el.c;
        if (g.tau_mi_in_u && std::abs(el.m_double) > 0) u -= exp(lt * el.m_double) * el.omega;
        if (g.eta1_in_u && valid(el.eta1)) u -= (D - el.epsilon1) * el.eta1;
        if (g.eta2_in_u && valid(el.eta2)) u -= (D - el.epsilon2) * (D - el.epsilon2) * el.eta2;
        if (g.beta1_in_u && valid(el.beta1)) u -= (T - el.gamma1) * el.beta1;
        if (g.beta2_in_u && valid(el.beta2)) u -= (T - el.gamma2) * (T - el.gamma2) * el.beta2;
        s += exp(lt * el.t + ld * el.d + u) * el.n;
    }
    double fi = 1;
    for (int i = 0; i <= NO; ++i, fi *= i) {
        double fj = 1;
        for (int j = 0; i + j <= NO; ++j, fj *= j)
            A[id(i, j)] = fi * fj * s(i, j);
    }
}

// -------- separable AD: univariate Taylor per factor, outer product per term
void separable(const ResidualHelmholtzGeneralizedExponential& g, double tau, double delta, double* A) {
    using U = nd::Taylor1<NO>;
    const U T = U::variable(tau, tau), D = U::variable(delta, delta);
    const U lt = log(T), ld = log(D);
    std::array<double, NS> c{};
    for (const auto& el : g.elements) {
        U ud(0.0), ut(0.0);
        if (g.delta_li_in_u && valid(el.l_double) && el.l_double > 0 && std::abs(el.c) > DBL_EPSILON) ud = ud - exp(ld * el.l_double) * el.c;
        if (g.tau_mi_in_u && std::abs(el.m_double) > 0) ut = ut - exp(lt * el.m_double) * el.omega;
        if (g.eta1_in_u && valid(el.eta1)) ud = ud - (D + (-el.epsilon1)) * el.eta1;
        if (g.eta2_in_u && valid(el.eta2)) ud = ud - (D + (-el.epsilon2)) * (D + (-el.epsilon2)) * el.eta2;
        if (g.beta1_in_u && valid(el.beta1)) ut = ut - (T + (-el.gamma1)) * el.beta1;
        if (g.beta2_in_u && valid(el.beta2)) ut = ut - (T + (-el.gamma2)) * (T + (-el.gamma2)) * el.beta2;
        const U F = exp(ld * el.d + ud), G = exp(lt * el.t + ut);
        for (int i = 0; i <= NO; ++i) {
            const double nG = el.n * G.t[i];
            for (int j = 0; i + j <= NO; ++j)
                c[id(i, j)] += nG * F.t[j];
        }
    }
    double fi = 1;
    for (int i = 0; i <= NO; ++i, fi *= i) {
        double fj = 1;
        for (int j = 0; i + j <= NO; ++j, fj *= j)
            A[id(i, j)] = fi * fj * c[id(i, j)];
    }
}

template <class F>
double time_ns(F&& f) {
    const int reps = 9, iters = 20000;
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int it = 0; it < iters; ++it)
            f(it);
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / iters);
    }
    return best;
}

int main() {
    struct Case
    {
        std::string fluid;
        double T, rho;  // K, mol/m^3
    };
    const std::vector<Case> cases = {{"Propane", 300, 11000}, {"Propane", 300, 100},   {"Nitrogen", 100, 25000},
                                     {"Nitrogen", 300, 400},  {"R1234yf", 300, 10000}, {"n-Decane", 400, 4500}};
    std::printf("fluid,T,rho,nterms,value_ns,hand_ns,hand_deltaonly_ns,naive_ns,separable_ns,maxrel_naive,maxrel_separable\n");
    volatile double sink = 0;
    for (const auto& cs : cases) {
        CoolProp::HelmholtzEOSBackend HEOS(cs.fluid);
        auto& eos = HEOS.get_components()[0].EOS();
        auto& ge = eos.alphar.GenExp;
        const double tau = eos.reduce.T / cs.T, delta = cs.rho / eos.reduce.rhomolar;
        // The GenExp block must be the whole residual for the comparison to mean anything
        const double full = eos.alphar.all(tau, delta).alphar, ge_only = value_only(ge, tau, delta);
        if (std::abs(full - ge_only) > 1e-12 * std::abs(full)) {
            std::fprintf(stderr, "%s: residual has non-GenExp terms (%.17g vs %.17g); skipping\n", cs.fluid.c_str(), full, ge_only);
            continue;
        }
        double Ah[NS], An[NS], As[NS], tmp[NS];
        hand(ge, tau, delta, Ah);
        naive(ge, tau, delta, An);
        separable(ge, tau, delta, As);
        double en = 0, es = 0;
        for (int m = 0; m < NS; ++m) {
            en = std::max(en, std::abs(An[m] - Ah[m]) / std::abs(Ah[m]));
            es = std::max(es, std::abs(As[m] - Ah[m]) / std::abs(Ah[m]));
        }
        auto jig = [&](int it) { return delta * (1.0 + 1e-12 * (it & 7)); };
        const double tv = time_ns([&](int it) { sink = sink + value_only(ge, tau, jig(it)); });
        const double th = time_ns([&](int it) {
            hand(ge, tau, jig(it), tmp);
            sink = sink + tmp[0];
        });
        const double td = time_ns([&](int it) {
            hand_deltaonly(ge, tau, jig(it), tmp);
            sink = sink + tmp[0];
        });
        const double tn = time_ns([&](int it) {
            naive(ge, tau, jig(it), tmp);
            sink = sink + tmp[0];
        });
        const double ts = time_ns([&](int it) {
            separable(ge, tau, jig(it), tmp);
            sink = sink + tmp[0];
        });
        std::printf("%s,%g,%g,%zu,%.1f,%.1f,%.1f,%.1f,%.1f,%.1e,%.1e\n", cs.fluid.c_str(), cs.T, cs.rho, ge.elements.size(), tv, th, td, tn, ts, en,
                    es);
    }
}
