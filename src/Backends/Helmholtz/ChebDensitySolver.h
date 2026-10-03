#pragma once
/**
 * Chebyshev tables of the pressure equation of a multi-fluid Helmholtz model, for finding ALL density roots
 * at given (T, p, x).
 *
 * With the reduced variables tau = Tr(x)/T and delta = rho/rhor(x), the pressure equation is
 *
 *     G(delta) = delta Z(delta) - t = 0,    t = p / (rhor(x) R T),    Z = 1 + delta d(alphar)/d(delta).
 *
 * Every term of a generalized-exponential residual (pure-fluid or departure function) factors as
 * kappa_k(tau) phi_k(delta), with phi_k = delta^d exp(-c delta^l - eta1 (delta - eps1) - eta2 (delta - eps2)^2).
 * Terms whose delta-parts are identical are grouped, so for a mixture
 *
 *     Z = 1 + sum_g W_g(T, x) chi_g(delta),     chi_g = delta d(phi_g)/d(delta),
 *     W_g = sum_{k in g} X_k kappa_k(tau),      X_k = x_i (pure-fluid term of i) or x_i x_j F_ij (departure term ij).
 *
 * The tiers:
 *  - Tables::build (once per component set): fits each chi_g on adaptive delta-pieces of [0, delta_max] by
 *    Chebyshev series of degree NQ, with a measured fit-error bound per piece;
 *  - Tables::assemble (once per (T, x)): combines the fits into G on every piece as a degree-NG series, with a
 *    bound ("margin") on |G_tables - G_true| per piece;
 *  - the caller (once per p): shifts the constant coefficient by -t and isolates the roots on every piece with
 *    ChebyshevBernstein::real_roots, using the margin as the sign tolerance.
 *
 * The non-analytic critical-region terms of IAPWS-95 and Span-Wagner type, n Delta^b delta psi, do not factor.  They
 * depend on (tau, delta) only, though, so per component they are tabulated at build as 2-D Chebyshev tables on
 * tau-cells x the delta-pieces (the steep Gaussian exp(-D (tau-1)^2) factored out and multiplied back exactly), and
 * assemble() contracts the table at its tau.  Per piece, in order:
 *  - a rigorous bound on their size (NonAnalyticTerm::bound_factor) shows them negligible: left out, the bound into
 *    the margin;
 *  - the cell's table, kept at build only where its measured error is negligible: contracted, error into the margin;
 *  - otherwise (where Delta ~ 0 along tau - 1 = A |delta - 1|^(1/beta) crosses the cell, a few pieces for tau > 1):
 *    fit at the actual tau by degree-NQ Lobatto interpolation with a measured error.
 * The pieces have an edge at delta = 1, where these terms are not analytic, and are graded geometrically toward it
 * (1 +- 2^-k), as are the tau-cells toward tau = 1, so that the interpolants converge there.
 *
 * The tables live on a hard rectangle tau in [tau_min, tau_max], delta in [0, delta_max]; nothing is
 * extrapolated, and assemble() declines a (T, x) outside it.  Models with other residual terms that do not factor
 * (SAFT association, cubic, Gao-B, ...) are declined at build.
 *
 * The margin bounds the difference from G as evaluated from the same grouped terms (true_G), not from the backend's
 * own pressure; the two agree to roundoff (build() verifies the regrouping to 1e-12 relative to the size of the terms,
 * the tests find ~1e-15).  The fit errors in the margin are measured -- the generalized-exponential fits between the
 * interpolation nodes (factor 2); the non-analytic tables on a 16 x 16 grid between their nodes (factor 10) and the
 * per-tau fits between their nodes (max of 10 x measured, 5 x tail); all doubled in the margin -- and its roundoff
 * allowance is empirical; neither is proven.
 *
 * Thread safety: a built Tables object is immutable and owns everything it reads (its own copy of the reducing
 * function, no pointers into the backend that built it), so it may be shared across threads and across backends
 * with the SAME model.  It is a snapshot: interaction parameters changed on the backend afterwards (F_ij, the
 * reducing function, departure terms) or a change of NORMALIZE_GAS_CONSTANTS / R_U_CODATA are not seen; a cache of
 * tables must be keyed on (or invalidated by) them.
 */

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "CoolProp/numerics/ChebyshevBernstein.h"

namespace CoolProp {

class HelmholtzEOSMixtureBackend;
class ReducingFunction;
class ResidualHelmholtzGeneralizedExponential;

namespace ChebDensity {

/// Degree of the per-group fits of chi_g(delta)
constexpr int NQ = 16;
/// Degree of G = delta Z - t on a piece (one more than NQ: the factor delta)
constexpr int NG = NQ + 1;
static_assert(NG <= ChebyshevBernstein::MAX_DEGREE, "no Chebyshev-to-Bernstein matrix for this degree");

using CoeffsQ = ChebyshevBernstein::Coeffs<NQ>;
using CoeffsG = ChebyshevBernstein::Coeffs<NG>;

/// One generalized-exponential term n delta^d tau^t exp(u), split as kappa(tau) phi(delta)
struct Term
{
    int i = 0, j = -1;  ///< component (j = -1 for a pure-fluid term; i < j for a departure term)
    double n = 0, d = 0, t = 0;
    bool has_cl = false, has_om = false, has_e1 = false, has_e2 = false, has_b1 = false, has_b2 = false;
    double c = 0, l = 0, om = 0, m = 0, e1 = 0, eps1 = 0, e2 = 0, eps2 = 0, b1 = 0, g1 = 0, b2 = 0, g2 = 0;
    int di = -1, li = -1;  ///< d and l as small non-negative integers (powers from a table), else -1

    /// phi = delta^d exp(u_delta)
    [[nodiscard]] double phi(double D) const;
    /// chi = delta d(phi)/d(delta)
    [[nodiscard]] double chi(double D) const;
    /// chi and d(chi)/d(delta); pw[k] = delta^k for k <= MAX_POW is used when the exponents are small integers
    void chi_d(double D, const double* pw, double& f, double& df) const;
    /// Size of the parts chi is summed from, delta^d e^u (|d| + delta |u'|, termwise): the scale of its rounding error
    [[nodiscard]] double parts(double D) const;
    /// phi, with pw as in chi_d
    [[nodiscard]] double phi_pw(double D, const double* pw) const;
    /// kappa = n tau^t exp(u_tau), with ltau = ln(tau)
    [[nodiscard]] double kappa(double tau, double ltau) const;
    /// true if phi is the same function of delta for both terms
    [[nodiscard]] bool same_delta(const Term& o) const;

    static constexpr int MAX_POW = 24;
};

/// One non-analytic critical-region term, alphar = n Delta^b delta psi with
///   Delta = theta^2 + B ((delta - 1)^2)^a,  theta = (1 - tau) + A ((delta - 1)^2)^(1/(2 beta)),
///   psi = exp(-C (delta - 1)^2 - D (tau - 1)^2)            (IAPWS-95 Table 5 form; also Span-Wagner CO2)
struct NonAnalyticTerm
{
    double n = 0, a = 0, b = 0, beta = 0, A = 0, B = 0, C = 0, D = 0;

    /// Bound on |delta d(alphar)/d(delta)| of this term over delta in [lo, hi] and |1 - tau| <= dtau_max, WITHOUT
    /// the factor exp(-D (tau - 1)^2): multiply by it for a given tau.  +inf when the parameters are outside the
    /// range the bound is derived for (beta in (0, 1], b >= 1/2, 2 a b >= 1, a >= 1/2, A, C >= 0, B > 0).
    [[nodiscard]] double bound_factor(double lo, double hi, double dtau_max) const;
};

/// Sum of non-analytic terms at (tau, delta), evaluated as ResidualHelmholtzNonAnalytic does (same expressions and the
/// same offset away from tau = 1 and delta = 1)
struct NonAnalyticValues
{
    double alphar = 0;
    double chi = 0;    ///< delta d(alphar)/d(delta)
    double dchi = 0;   ///< d(chi)/d(delta)
    double parts = 0;  ///< roundoff scale of chi: the size of what it is summed from, incl. the cancellation in Delta
};
NonAnalyticValues eval_nonanalytic(const std::vector<NonAnalyticTerm>& terms, double tau, double delta);

/// Split the terms of a generalized-exponential residual into Term (appended to out)
void collect_terms(const ResidualHelmholtzGeneralizedExponential& g, int i, int j, std::vector<Term>& out);

struct BuildOptions
{
    double delta_max = 4.0;   ///< upper end of the delta range
    double tau_min = 0.0;     ///< lower end of the tau range (0: any tau > 0)
    double tau_max = 0.0;     ///< upper end of the tau range; required (> tau_min)
    double tol = 1e-6;        ///< allowed fit error of each weighted group, in units of Z (see build())
    double min_width = 1e-3;  ///< pieces are not split below this width in delta (at least 1e-9 delta_max)
};

class Tables
{
   public:
    /// Per-(T, x) state: G on every piece and its error bound.  One per thread.
    struct State
    {
        std::vector<CoeffsG> G;      ///< Chebyshev coefficients of delta Z(delta) on piece p (t not subtracted)
        std::vector<double> margin;  ///< bound on |G_tables - G_true| on piece p
        std::vector<double> W;       ///< group weights W_g(T, x)
        std::vector<double> x;       ///< mole fractions
        double T = 0, tau = 0, rhor = 0;
        double t_scale = 0;  ///< t = p * t_scale = p / (rhor R T)
        long na_fits = 0;    ///< pieces on which a non-analytic contribution was fit (rest: bounded into the margin)
        long na_table = 0;   ///< pieces on which it came from the 2-D table instead (not counted in na_fits)
    };

    /// Build the tables for the components of HEOS (its mole fractions must be set: they are used to identify how
    /// the backend forms the mixture gas constant, not otherwise).  Returns nullptr when the model is declined; the reason is
    /// written to *reason when given.  Throws only for invalid options.
    [[nodiscard]] static std::shared_ptr<const Tables> build(HelmholtzEOSMixtureBackend& HEOS, const BuildOptions& opt,
                                                             std::string* reason = nullptr);

    /// Combine the tables at (T, x).  Returns false (S unspecified) if tau = Tr(x)/T is outside the rectangle,
    /// x has the wrong length, or anything is non-finite.
    [[nodiscard]] bool assemble(double T, const std::vector<double>& x, State& S) const;

    /// The true G(delta) = delta Z - t at the state of S (evaluated from the same grouped terms, no tables), its
    /// delta-derivative, and a roundoff scale (|G| below ~1e-15 * scale is noise)
    void true_G(const State& S, double delta, double t, double& G, double& dG, double& scale) const;
    /// alphar(tau, delta) at the state of S
    [[nodiscard]] double alphar(const State& S, double delta) const;

    /// delta at piece variable u in [-1, 1] of piece p
    [[nodiscard]] double delta_of(int p, double u) const {
        return m_edges[p] + (m_edges[p + 1] - m_edges[p]) * (u + 1) / 2;
    }
    [[nodiscard]] int n_pieces() const {
        return static_cast<int>(m_edges.size()) - 1;
    }
    [[nodiscard]] const std::vector<double>& edges() const {
        return m_edges;
    }
    [[nodiscard]] std::size_t n_terms() const {
        return m_terms.size();
    }
    [[nodiscard]] std::size_t n_groups() const {
        return m_reps.size();
    }
    [[nodiscard]] std::size_t n_components() const {
        return static_cast<std::size_t>(m_N);
    }
    [[nodiscard]] bool has_nonanalytic() const {
        return !m_na.empty();
    }
    [[nodiscard]] const BuildOptions& options() const {
        return m_opt;
    }

   private:
    Tables() = default;
    [[nodiscard]] std::string verify(HelmholtzEOSMixtureBackend& HEOS) const;
    struct NAComp;
    [[nodiscard]] bool build_na_tables(NAComp& c) const;

    BuildOptions m_opt;
    int m_N = 0;
    std::vector<Term> m_terms;
    std::vector<int> m_group;     ///< term -> group
    std::vector<Term> m_reps;     ///< one representative term per group (its delta-part)
    std::vector<double> m_edges;  ///< piece edges, m_edges[0] = 0, back() = delta_max
    std::vector<CoeffsQ> m_C;     ///< fit of chi_g on piece p at [p * n_groups + g]
    std::vector<double> m_Cn;     ///< l1 norm of each fit (roundoff scale)
    std::vector<double> m_Ct;     ///< fit-error bound of each fit
    std::vector<double> m_Cp;     ///< roundoff scale of evaluating chi_g on the piece (see Term::parts)
    struct NAComp
    {
        int i = 0;                           ///< component
        std::vector<NonAnalyticTerm> terms;  ///< owned copy of its non-analytic terms
        std::vector<double> K;               ///< bound factor of term k on piece p, rung r: [(p * terms.size() + k) * rungs + r]
        // 2-D tables of chi_NA(tau, delta) / exp(-D (tau-1)^2) per distinct D, on tau-cells x the delta-pieces (see build)
        std::vector<double> Dg;                            ///< the distinct D
        std::vector<std::vector<NonAnalyticTerm>> gterms;  ///< the terms with each D, with D set to 0
        std::vector<double> tau_edges;                     ///< tau-cells (empty: no tables)
        std::vector<int> cell;                             ///< [q * P + p]: index of the stored table, -1: fit at assemble, -2: negligible
        std::vector<double> cell_err;                      ///< [q * P + p]: measured error of the table in chi (incl. max exp factor)
        std::vector<double> cell_parts;                    ///< [q * P + p]: roundoff scale of chi on the cell
        std::vector<double> coef;                          ///< stored tables, each Dg.size() x (NQ+1) x (NQ+1): [g][m (tau)][k (delta)]
    };
    std::vector<NAComp> m_na;
    std::vector<double> m_na_dtau;            ///< the |1 - tau| ladder of the bounds (increasing; last = max over the tau range)
    std::vector<std::vector<double>> m_F;     ///< F_ij
    std::shared_ptr<ReducingFunction> m_red;  ///< owned copy
    std::vector<double> m_Ri;                 ///< component gas constants
    double m_R = 0;                           ///< gas constant, when it does not depend on x
    bool m_R_mixed = false;                   ///< true: R = sum x_i R_i
};

}  // namespace ChebDensity
}  // namespace CoolProp
