/**
 This C++ code is the implementation of the analyses presented in the paper
 I.Bell and A. Jäger, "Helmholtz energy translations for common cubic equations of state
 for use in one-fluid and multi-fluid mixture models", J. Res. NIST, 2016

 This code is in the public domain, though if used in academic work, we would appreciate
 a reference back to the paper given above.

 */

#ifndef CUBIC_H
#define CUBIC_H

#include <utility>
#include <vector>
#include <cmath>
#include <cstring>
#include <limits>
#include <array>
#include <memory>
using std::shared_ptr;
#include "CoolProp/Exceptions.h"

/// An abstract alpha function for the EOS, defining the interface for the alpha function
class AbstractCubicAlphaFunction
{
   protected:
    double a0;              ///< The constant term multiplying the alpha function
    double Tr_over_Tci,     ///< The (constant) reducing temperature divided by the critical temperature of the pure component
      sqrt_Tr_Tci;          ///< The sqrt of the (constant) reducing temperature divided by the critical temperature of the pure component
    std::vector<double> c;  ///< The vector of constants
    /// Incremented every time a parameter that affects term() changes (e.g., set_Tr_over_Tci());
    /// used by AbstractCubic's aii cache to detect when cached term() values are stale.
    unsigned long m_version = 0;

   public:
    virtual ~AbstractCubicAlphaFunction() = default;
    /**
     * \brief A "term" is the value of a_ii(tau) = a0_ii * alpha(tau), or one of its tau-derivatives,
     * for this pure component's attractive energy parameter.
     *
     * \param tau The reciprocal reduced temperature, Tr/T
     * \param itau The order of the derivative with respect to tau to return: itau=0 returns a_ii
     *             itself, itau=1 returns da_ii/dtau, itau=2 the second derivative, and so on.
     *             Implementations support itau in the range 0-4 and throw for higher orders.
     */
    virtual double term(double tau, std::size_t itau) = 0;
    /// Compute all 5 tau-derivatives at once. Default falls back to 5 separate term() calls.
    virtual void calc_all_terms(double tau, std::array<double, 5>& terms) {
        for (int k = 0; k < 5; ++k)
            terms[k] = term(tau, k);
    }
    void set_Tr_over_Tci(double Tr_over_Tci) {
        this->Tr_over_Tci = Tr_over_Tci;
        this->sqrt_Tr_Tci = sqrt(Tr_over_Tci);
        ++m_version;
    };
    /// Update the leading constant a0 that scales the alpha function (e.g., after a change to Tc or pc)
    void set_a0(double a0_new) {
        this->a0 = a0_new;
        ++m_version;
    };
    /// Monotonically increasing counter bumped whenever this alpha function's parameters change;
    /// callers that cache term() values (e.g., AbstractCubic::_ensure_aii_cache) can compare this
    /// against a previously-recorded value to detect stale cache entries.
    unsigned long version() const {
        return m_version;
    }
    AbstractCubicAlphaFunction(double a0, double Tr_over_Tci) : a0(a0), Tr_over_Tci(Tr_over_Tci), sqrt_Tr_Tci(sqrt(Tr_over_Tci)) {};
};

/// An implementation of AbstractCubicAlphaFunction for the baseline alpha function of PR or SRK
class BasicMathiasCopemanAlphaFunction : public AbstractCubicAlphaFunction
{
    double m;  ///< The term coming from the function of omega
   public:
    BasicMathiasCopemanAlphaFunction(double a0, double m_ii, double Tr_over_Tci)
      : AbstractCubicAlphaFunction(a0, Tr_over_Tci), m(m_ii) {

        };
    double term(double tau, std::size_t itau) override;
    void calc_all_terms(double tau, std::array<double, 5>& terms) override;
};

/// An implementation of AbstractCubicAlphaFunction for the Twu alpha function
class TwuAlphaFunction : public AbstractCubicAlphaFunction
{
   public:
    TwuAlphaFunction(double a0, double L, double M, double N, double Tr_over_Tci) : AbstractCubicAlphaFunction(a0, Tr_over_Tci) {
        c.resize(3);
        c[0] = L;
        c[1] = M;
        c[2] = N;
    };
    double term(double tau, std::size_t itau) override;
    void calc_all_terms(double tau, std::array<double, 5>& terms) override;
};

/// An implementation of AbstractCubicAlphaFunction for the Mathias-Copeman alpha function
class MathiasCopemanAlphaFunction : public AbstractCubicAlphaFunction
{
   public:
    MathiasCopemanAlphaFunction(double a0, double c1, double c2, double c3, double Tr_over_Tci) : AbstractCubicAlphaFunction(a0, Tr_over_Tci) {
        c.resize(3);
        c[0] = c1;
        c[1] = c2;
        c[2] = c3;
    };
    double term(double tau, std::size_t itau) override;
    void calc_all_terms(double tau, std::array<double, 5>& terms) override;
};

class AbstractCubic
{

   protected:
    double rho_r,                                               ///< The reducing density to be used [mol/m3]
      T_r;                                                      ///< The reducing temperature to be used [K]
    std::vector<double> Tc,                                     ///< Vector of critical temperatures (in K)
      pc,                                                       ///< Vector of critical pressures (in Pa)
      acentric;                                                 ///< Vector of acentric factors (unitless)
    double R_u;                                                 ///< The universal gas constant  in J/(mol*K)
    double Delta_1,                                             ///< The first cubic constant
      Delta_2;                                                  ///< The second cubic constant
    int N;                                                      ///< Number of components in the mixture
    std::vector<std::vector<double>> k;                         ///< The interaction parameters (k_ii = 0)
    double cm;                                                  ///< The volume translation parameter
    std::vector<shared_ptr<AbstractCubicAlphaFunction>> alpha;  ///< The vector of alpha functions for the pure components
    /// Cache: aii_term values for the most recent tau.  Populated lazily by _ensure_aii_cache().
    mutable double m_tau_cache;
    mutable std::vector<std::array<double, 5>> m_aii_cache;
    /// alpha[i]->version() values recorded when m_aii_cache was last populated; used together with
    /// m_tau_cache to detect when an alpha function was mutated (e.g., via set_Tr_over_Tci()) without
    /// going through one of AbstractCubic's own cache-invalidating setters.
    mutable std::vector<unsigned long> m_alpha_versions_cache;
    /// Cache: per-component covolumes b0_ii(i).  These depend only on Tc, pc and R_u (fixed at
    /// construction), so they are computed once on first use; bm_term() reads them instead of
    /// re-evaluating the R_u*Tc/pc division on every call.
    mutable std::vector<double> m_b0_ii_cache;
    mutable std::vector<char> m_b0_ii_cache_valid;
    void _ensure_aii_cache(double tau) const {
        bool stale = (std::memcmp(&tau, &m_tau_cache, sizeof(double)) != 0) || (m_alpha_versions_cache.size() != alpha.size());
        if (!stale) {
            for (int i = 0; i < N; ++i) {
                if (alpha[i]->version() != m_alpha_versions_cache[i]) {
                    stale = true;
                    break;
                }
            }
        }
        if (stale) {
            m_aii_cache.resize(N);
            m_alpha_versions_cache.resize(N);
            for (int i = 0; i < N; ++i) {
                alpha[i]->calc_all_terms(tau, m_aii_cache[i]);
                m_alpha_versions_cache[i] = alpha[i]->version();
            }
            m_tau_cache = tau;
            m_aux_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate s/aij caches built from a different aii state
        }
    }

    /// Caches built ON TOP of the aii cache, keyed on the same tau (#3192 cubic perf):
    ///   m_s_cache[i][p]    = d^p/dtau^p of sqrt(a_i(tau))  (the geometric-mean factor s_i and its tau-derivatives)
    ///   m_aij_cache[i][j]  = the 5 tau-derivatives of a_ij = (1 - k_ij) * sqrt(a_i*a_j) = (1 - k_ij) * s_i*s_j
    /// Both are composition-independent and constant at fixed tau, so they are built once per tau and reused across
    /// every phase/composition/derivative evaluation of a flash (was: re-derived ~1e6x per flash).  s_i is obtained
    /// from the cached a_i derivatives by the recurrence s^(n) = (a^(n) - sum_{p=1}^{n-1} C(n,p) s^(p) s^(n-p))/(2 s0).
    mutable double m_aux_tau_cache = std::numeric_limits<double>::quiet_NaN();
    mutable std::vector<std::array<double, 5>> m_s_cache;
    mutable std::vector<std::vector<std::array<double, 5>>> m_aij_cache;
    void _ensure_aux_cache(double tau) const {
        _ensure_aii_cache(tau);  // refreshes m_aii_cache (and NaNs m_aux_tau_cache on a rebuild)
        if (std::memcmp(&tau, &m_aux_tau_cache, sizeof(double)) == 0) return;
        // --- per-component sqrt(a_i) and its tau-derivatives, via the s^2 = a recurrence ---
        m_s_cache.resize(N);
        static const double binom[5][5] = {{1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};
        for (int i = 0; i < N; ++i) {
            const std::array<double, 5>& a = m_aii_cache[i];
            std::array<double, 5>& s = m_s_cache[i];
            s[0] = std::sqrt(a[0]);
            const double inv2s0 = 1.0 / (2.0 * s[0]);
            for (int n = 1; n < 5; ++n) {
                double conv = 0.0;
                for (int p = 1; p < n; ++p)
                    conv += binom[n][p] * s[p] * s[n - p];
                s[n] = (a[n] - conv) * inv2s0;
            }
        }
        // --- a_ij = (1-k_ij)*s_i*s_j and its tau-derivatives (Leibniz over s_i*s_j) ---
        m_aij_cache.assign(N, std::vector<std::array<double, 5>>(N));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) {
                const double one_minus_k = 1.0 - k[i][j];
                std::array<double, 5>& aij = m_aij_cache[i][j];
                for (int n = 0; n < 5; ++n) {
                    double d = 0.0;
                    for (int p = 0; p <= n; ++p)
                        d += binom[n][p] * m_s_cache[i][p] * m_s_cache[j][n - p];
                    aij[n] = one_minus_k * d;
                }
                if (j != i) m_aij_cache[j][i] = aij;  // symmetric
            }
        }
        m_aux_tau_cache = tau;
    }

    /// Composition-scoped cache of bm = sum_i x_i b_i (#3192 perf).  bm is a pure function of the
    /// composition, but the residual-Helmholtz volume helpers (c_term/A_term/psi_*) recompute it
    /// ~dozens of times per alphar-derivative evaluation.  _sync_comp(x) is called once at each
    /// alphar-level entry (alphar / d_alphar_dxi / d2_alphar_dxidxj / d3_alphar_dxidxjdxk) and rebuilds
    /// m_bm only when the composition actually changed (one O(N) memcmp -- a composition "generation"
    /// check); the nested helpers then read m_bm directly (O(1)).  Those helpers are private to
    /// GeneralizedCubic and reachable ONLY through those four entries, so m_bm is current where read.
    std::vector<double> m_comp_x_cache;
    double m_bm = std::numeric_limits<double>::quiet_NaN();
    void _sync_comp(const std::vector<double>& x) {
        const std::size_t n = static_cast<std::size_t>(N);
        if (m_comp_x_cache.size() == n && x.size() >= n && std::memcmp(x.data(), m_comp_x_cache.data(), n * sizeof(double)) == 0) {
            return;
        }
        // Virtual dispatch: linear sum_i x_i b_i for SRK/PR, quadratic sum_ij x_i x_j b_ij for VTPR.
        m_bm = bm_term(x);
        m_comp_x_cache.assign(x.begin(), x.begin() + N);
    }

   public:
    /**
     \brief The abstract base clase for the concrete implementations of the cubic equations of state

     This abstract base class describes the structure that must be implemented by concrete implementations
     of the cubic equations of state (SRK, PR, etc.).  The virtual functions must be implemented by the
     derived classes, the remaining functions are generic and are not dependent on the equation of state,
     so long as it has the formulation given in this work.

     */
    AbstractCubic(const std::vector<double>& Tc, std::vector<double> pc, std::vector<double> acentric, double R_u, double Delta_1, double Delta_2,
                  const std::vector<double>& C1 = std::vector<double>(), const std::vector<double>& C2 = std::vector<double>(),
                  const std::vector<double>& C3 = std::vector<double>());
    virtual ~AbstractCubic() = default;
    /// Set the constants for the Mathias-Copeman alpha function, or if C1,C2,C3 are all empty, set the default alpha model
    void set_alpha(const std::vector<double>& C1, const std::vector<double>& C2, const std::vector<double>& C3);
    /// Set the alpha function for the i-th component
    void set_alpha_function(std::size_t i, shared_ptr<AbstractCubicAlphaFunction>& acaf) {
        alpha[i] = acaf;
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    };
    /// Get the alpha function for the i-th component
    shared_ptr<AbstractCubicAlphaFunction> get_alpha_function(std::size_t i) {
        return alpha[i];
    };
    /// Get all the alpha functions
    const std::vector<shared_ptr<AbstractCubicAlphaFunction>>& get_all_alpha_functions() {
        return this->alpha;
    };
    /// Set all the alpha functions
    void set_all_alpha_functions(const std::vector<shared_ptr<AbstractCubicAlphaFunction>>& alpha) {
        this->alpha = alpha;
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    };

    /// Get the entire kij matrix in one shot
    const std::vector<std::vector<double>>& get_kmat() {
        return k;
    };
    /// Set the entire kij matrix in one shot
    void set_kmat(const std::vector<std::vector<double>>& k) {
        this->k = k;
        m_aux_tau_cache = std::numeric_limits<double>::quiet_NaN();  // k_ij is baked into m_aij_cache; force a rebuild
    };
    /// Set the kij factor for the ij pair
    void set_kij(std::size_t i, std::size_t j, double val) {
        k[i][j] = val;
        k[j][i] = val;
        m_aux_tau_cache = std::numeric_limits<double>::quiet_NaN();  // k_ij is baked into m_aij_cache; force a rebuild
    }
    /// Get the kij factor for the ij pair
    double get_kij(std::size_t i, std::size_t j) {
        return k[i][j];
    }
    /// Get the vector of critical temperatures (in K)
    std::vector<double>& get_Tc() {
        return Tc;
    }
    /// Get the vector of critical pressures (in Pa)
    std::vector<double>& get_pc() {
        return pc;
    }
    /// Get the vector of acentric factors
    std::vector<double>& get_acentric() {
        return acentric;
    }
    /// Read-only accessor for value of Delta_1
    double get_Delta_1() {
        return Delta_1;
    }
    /// Read-only accessor for value of Delta_2
    double get_Delta_2() {
        return Delta_2;
    }
    /// Read-only accessor for value of R_u (universal gas constant)
    double get_R_u() {
        return R_u;
    }
    /// Set the reducing temperature to be used
    void set_Tr(double Tr) {
        T_r = Tr;
        for (std::size_t i = 0; i < alpha.size(); ++i) {
            alpha[i]->set_Tr_over_Tci(T_r / Tc[i]);
        }
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    }
    /// Set the reducing density to be used
    void set_rhor(double rhor) {
        rho_r = rhor;
    }
    /// Get the reducing temperature to be used
    double get_Tr() {
        return T_r;
    }
    /// Get the reducing density to be used
    double get_rhor() {
        return rho_r;
    }

    /**
     * \brief Set the critical temperature of the i-th component [K]
     *
     * Updates Tc[i], refreshes the alpha function's a0 and Tr/Tci ratio, and
     * invalidates the b0 and aii caches so that all subsequent evaluations pick
     * up the new value consistently.
     */
    void set_Tci(std::size_t i, double Tci) {
        Tc[i] = Tci;
        // a0_ii depends on Tc[i]^2/pc[i], and Tr_over_Tci = T_r/Tc[i].
        // Refresh both in the existing alpha function so its parameters stay consistent.
        alpha[i]->set_a0(a0_ii(i));
        alpha[i]->set_Tr_over_Tci(T_r / Tc[i]);
        // b0_ii depends on Tc[i]/pc[i]; mark the affected entry as stale.
        if (i < m_b0_ii_cache.size()) {
            m_b0_ii_cache_valid[i] = 0;
        }
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate aii cache
    }
    /**
     * \brief Set the critical pressure of the i-th component [Pa]
     *
     * Updates pc[i], refreshes the alpha function's a0, and invalidates the b0
     * cache so that all subsequent evaluations pick up the new value consistently.
     */
    void set_pci(std::size_t i, double pci) {
        pc[i] = pci;
        // a0_ii depends on Tc[i]^2/pc[i]; refresh it in the existing alpha function.
        alpha[i]->set_a0(a0_ii(i));
        // b0_ii depends on Tc[i]/pc[i]; mark the affected entry as stale.
        if (i < m_b0_ii_cache.size()) {
            m_b0_ii_cache_valid[i] = 0;
        }
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate aii cache
    }
    /// Set the three Mathias-Copeman constants in one shot for the component i of a mixture
    void set_C_MC(std::size_t i, double c1, double c2, double c3) {
        alpha[i] = std::make_shared<MathiasCopemanAlphaFunction>(a0_ii(i), c1, c2, c3, T_r / Tc[i]);
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    }
    /// Set the three Twu constants in one shot for the component i of a mixture
    void set_C_Twu(std::size_t i, double L, double M, double N) {
        alpha[i] = std::make_shared<TwuAlphaFunction>(a0_ii(i), L, M, N, T_r / Tc[i]);
        m_tau_cache = std::numeric_limits<double>::quiet_NaN();  // invalidate cache
    }
    /// Get the leading constant in the expression for the pure fluid attractive energy term
    /// (must be implemented by derived classes)
    virtual double a0_ii(std::size_t i) = 0;
    /// Get the leading constant in the expression for the pure fluid covolume term
    /// (must be implemented by derived classes)
    virtual double b0_ii(std::size_t i) = 0;
    /// Get the m_ii variable in the alpha term inculuded in the attractive part
    virtual double m_ii(std::size_t i) = 0;

    /// The residual non-dimensionalized Helmholtz energy \f$\alpha^r\f$
    virtual double alphar(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta);
    /// The first composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d_alphar_dxi(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                bool xN_independent);
    /// The second composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d2_alphar_dxidxj(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                    std::size_t j, bool xN_independent);
    /// The third composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d3_alphar_dxidxjdxk(double tau, double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i,
                                       std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The n-th derivative of \f$a_m\f$ with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     */
    virtual double am_term(double tau, const std::vector<double>& x, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_am_term_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_am_term_dxidxj(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d3_am_term_dxidxjdxk(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k,
                                        bool xN_independent);

    /**
     * \brief The term \f$b_{\rm m}\f$ (mixture co-volume)
     * \param x The vector of mole fractions
     */
    virtual double bm_term(const std::vector<double>& x);
    /** \brief The first composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_bm_term_dxi(const std::vector<double>& x, std::size_t i, bool xN_independent);
    /** \brief The second composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_bm_term_dxidxj(const std::vector<double>& x, std::size_t i, std::size_t j, bool xN_independent);
    /** \brief The third composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d3_bm_term_dxidxjdxk(const std::vector<double>& x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    /**
	* \brief The term \f$c_{\rm m}\f$ (volume translation)
	*/
    virtual double cm_term();
    /// Set the volume translation parameter
    void set_cm(double val) {
        cm = val;
    }
    /// Get the volume translation parameter
    double get_cm() {
        return cm;
    }

    /// Modify the surface parameter Q_k of the sub group sgi
    virtual void set_Q_k(const size_t sgi, const double value) {
        throw CoolProp::ValueError("set_Q_k not defined for AbstractCubic");
    };
    /// Retrieve the surface parameter Q_k of the sub group sgi
    virtual double get_Q_k(const size_t sgi) const {
        throw CoolProp::ValueError("get_Q_k not defined for AbstractCubic");
    };

    /**
     * \brief The n-th \f$\tau\f$ derivative of \f$a_{ij}(\tau)\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The first index
     * \param j The second index
     * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double aij_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
    /** The n-th tau derivative of \f$u(\tau)\f$, the argument of sqrt in the cross aij term
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The first index
     * \param j The first index
     * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double u_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
    /** Take the n-th tau derivative of the \f$a_{ii}(\tau)\f$ pure fluid contribution
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The index of the component
     * \param itau The number of derivatives of \f$u\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double aii_term(double tau, std::size_t i, std::size_t itau);

    /**
     * \brief The term \f$ \psi^{(-)}\f$ and its \f$\tau\f$ and \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     */
    double psi_minus(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta);
    /**
     * \brief The third composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_psi_minus_dxi(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_psi_minus_dxidxj(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j,
                               bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_psi_minus_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j,
                                  std::size_t k, bool xN_independent);

    /**
     * \brief The term \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[ \Pi_{12} = (1+\Delta_1\bm\rhor \delta)(1+\Delta_2\bm\rhor \delta) \f]
     */
    double PI_12(double delta, const std::vector<double>& x, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_PI_12_dxi(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_PI_12_dxidxj(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_PI_12_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k,
                              bool xN_independent);

    /**
     * \brief The term \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     */
    double tau_times_a(double tau, const std::vector<double>& x, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_tau_times_a_dxi(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_tau_times_a_dxidxj(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_tau_times_a_dxidxjdxk(double tau, const std::vector<double>& x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k,
                                    bool xN_independent);

    /**
     * \brief The term \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[  \psi^{(+)} = \dfrac{\ln\left(\dfrac{\Delta_1\bm\rhor \delta+1}{\Delta_2\bm\rhor \delta+1}\right)}{\bm(\Delta_1-\Delta_2)}  \f]
     */
    double psi_plus(double delta, const std::vector<double>& x, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_psi_plus_dxi(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_psi_plus_dxidxj(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_psi_plus_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k,
                                 bool xN_independent);

    /** \brief The term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     *
     * \f$c\f$ is given by
     * \f[
     * c = \frac{1}{b_m}
     * \f]
     * \param x The vector of mole fractions
     */
    double c_term(const std::vector<double>& x) {
        return 1 / m_bm;
    };
    /**
     * \brief The first composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_c_term_dxi(const std::vector<double>& x, std::size_t i, bool xN_independent) {
        return -d_bm_term_dxi(x, i, xN_independent) / pow(m_bm, 2);
    };
    /**
     * \brief The second composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_c_term_dxidxj(const std::vector<double>& x, std::size_t i, std::size_t j, bool xN_independent) {
        double bm = m_bm;
        return (2 * d_bm_term_dxi(x, i, xN_independent) * d_bm_term_dxi(x, j, xN_independent) - bm * d2_bm_term_dxidxj(x, i, j, xN_independent))
               / pow(bm, 3);
    };
    /**
     * \brief The third composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_c_term_dxidxjdxk(const std::vector<double>& x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent) {
        double bm = m_bm;
        return 1 / pow(bm, 4)
               * (2 * bm
                    * (d_bm_term_dxi(x, i, xN_independent) * d2_bm_term_dxidxj(x, j, k, xN_independent)
                       + d_bm_term_dxi(x, j, xN_independent) * d2_bm_term_dxidxj(x, i, k, xN_independent)
                       + d_bm_term_dxi(x, k, xN_independent) * d2_bm_term_dxidxj(x, i, j, xN_independent))
                  - pow(bm, 2) * d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent)
                  - 6 * d_bm_term_dxi(x, i, xN_independent) * d_bm_term_dxi(x, j, xN_independent) * d_bm_term_dxi(x, k, xN_independent));
    };

    /**
     * \brief The term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     *
     * \f[
     * A = \log\left(\frac{\Delta_1\delta\rho_r b_m+1}{\Delta_2\delta\rho_r b+1}\right)
     * \f]
     *
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     */
    double A_term(double delta, const std::vector<double>& x) {
        double bm = m_bm;
        double cm = cm_term();
        return log((delta * rho_r * (Delta_1 * bm + cm) + 1) / (delta * rho_r * (Delta_2 * bm + cm) + 1));
    };
    /**
     * \brief The first composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_A_term_dxi(double delta, const std::vector<double>& x, std::size_t i, bool xN_independent) {
        std::size_t idelta = 0;
        return delta * rho_r * d_bm_term_dxi(x, i, xN_independent) * (Delta_1 - Delta_2) / PI_12(delta, x, idelta);
    };
    /**
     * \brief The second composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_A_term_dxidxj(double delta, const std::vector<double>& x, std::size_t i, std::size_t j, bool xN_independent) {
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, x, idelta);
        return delta * rho_r * (Delta_1 - Delta_2) / pow(PI12, 2)
               * (PI12 * d2_bm_term_dxidxj(x, i, j, xN_independent)
                  - d_PI_12_dxi(delta, x, 0, j, xN_independent) * d_bm_term_dxi(x, i, xN_independent));
    };
    /**
     * \brief The third composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_A_term_dxidxjdxk(double delta, const std::vector<double>& x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent) {
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, x, idelta);
        // The leading factor
        double lead = delta * rho_r * (Delta_1 - Delta_2) / pow(PI12, 3);
        return lead
               * (-PI12
                    * (d_PI_12_dxi(delta, x, idelta, j, xN_independent) * d2_bm_term_dxidxj(x, i, k, xN_independent)
                       + d_PI_12_dxi(delta, x, idelta, k, xN_independent) * d2_bm_term_dxidxj(x, i, j, xN_independent)
                       + d_bm_term_dxi(x, i, xN_independent) * d2_PI_12_dxidxj(delta, x, idelta, j, k, xN_independent))
                  + pow(PI12, 2) * d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent)
                  + 2 * d_PI_12_dxi(delta, x, idelta, j, xN_independent) * d_PI_12_dxi(delta, x, idelta, k, xN_independent)
                      * d_bm_term_dxi(x, i, xN_independent));
    };
    // Allows to modify the unifac interaction parameters aij, bij and cij. Only for use with VTPR backend.
    virtual void set_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter, const double value) {
        throw CoolProp::NotImplementedError("set_interaction_parameter is not implemented for this backend");
    }
    // Allows to get the unifac interaction parameters aij, bij and cij. Only for use with VTPR backend.
    virtual double get_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter) {
        throw CoolProp::NotImplementedError("get_interaction_parameter is not implemented for this backend");
    }
};

class PengRobinson : public AbstractCubic
{
   public:
    PengRobinson(const std::vector<double>& Tc, std::vector<double> pc, std::vector<double> acentric, double R_u,
                 const std::vector<double>& C1 = std::vector<double>(), const std::vector<double>& C2 = std::vector<double>(),
                 const std::vector<double>& C3 = std::vector<double>())
      : AbstractCubic(Tc, std::move(pc), std::move(acentric), R_u, 1 + sqrt(2.0), 1 - sqrt(2.0), C1, C2, C3) {
        set_alpha(C1, C2, C3);
    };

    PengRobinson(double Tc, double pc, double acentric, double R_u)
      : AbstractCubic(std::vector<double>(1, Tc), std::vector<double>(1, pc), std::vector<double>(1, acentric), R_u, 1 + sqrt(2.0), 1 - sqrt(2.0)) {
        set_alpha(std::vector<double>(), std::vector<double>(), std::vector<double>());
    };

    double a0_ii(std::size_t i) override;
    double b0_ii(std::size_t i) override;
    double m_ii(std::size_t i) override;
};

class SRK : public AbstractCubic
{
   public:
    SRK(const std::vector<double>& Tc, std::vector<double> pc, std::vector<double> acentric, double R_u,
        const std::vector<double>& C1 = std::vector<double>(), const std::vector<double>& C2 = std::vector<double>(),
        const std::vector<double>& C3 = std::vector<double>())
      : AbstractCubic(Tc, std::move(pc), std::move(acentric), R_u, 1, 0, C1, C2, C3) {
        set_alpha(C1, C2, C3);
    };
    SRK(double Tc, double pc, double acentric, double R_u)
      : AbstractCubic(std::vector<double>(1, Tc), std::vector<double>(1, pc), std::vector<double>(1, acentric), R_u, 1, 0) {
        set_alpha(std::vector<double>(), std::vector<double>(), std::vector<double>());
    };

    double a0_ii(std::size_t i) override;
    double b0_ii(std::size_t i) override;
    double m_ii(std::size_t i) override;
};

#endif
