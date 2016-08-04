/**
 This C++ code is the implementation of the analyses presented in the paper
 I.Bell and A. Jäger, "Helmholtz energy translations for common cubic equations of state
 for use in one-fluid and multi-fluid mixture models", J. Res. NIST, 2016
 
 This code is in the public domain, though if used in academic work, we would appreciate
 a reference back to the paper given above.
 
 */

#ifndef CUBIC_H
#define CUBIC_H

#include <vector>
#include <cmath>

enum AII_Model {AII_MODEL_NOT_SPECIFIED = 0, AII_MATHIAS_COPEMAN, AII_TWU};
class AbstractCubic
{
protected:
    std::vector<double> Tc, ///< Vector of critical temperatures (in K)
    pc, ///< Vector of critical pressures (in Pa)
    acentric; ///< Vector of acentric factors (unitless)
    double R_u; ///< The universal gas constant  in J/(mol*K)
    double Delta_1, ///< The first cubic constant
    Delta_2; ///< The second cubic constant
    int N; ///< Number of components in the mixture
    std::vector< std::vector<double> > k; ///< The interaction parameters (k_ii = 0)
    bool simple_aii; ///< True if the Mathias-Copeman equation for a_ii is not being used
    std::vector<double> C1, C2, C3; ///< The Mathias-Copeman coefficients for a_ii
    std::vector<double> L_Twu, M_Twu, N_Twu; ///< The Twu coefficients for a_ii
    AII_Model aii_model; ///< Enumeration for the aii model in use
    double cm; ///< The volume translation parameter
public:
    static const double rho_r, T_r;
    /**
     \brief The abstract base clase for the concrete implementations of the cubic equations of state
     
     This abstract base class describes the structure that must be implemented by concrete implementations
     of the cubic equations of state (SRK, PR, etc.).  The virtual functions must be implemented by the
     derived classes, the remaining functions are generic and are not dependent on the equation of state,
     so long as it has the formulation given in this work.
     
     */
    AbstractCubic(std::vector<double> Tc,
                  std::vector<double> pc,
                  std::vector<double> acentric,
                  double R_u,
                  double Delta_1,
                  double Delta_2,
                  std::vector<double> C1 = std::vector<double>(),
                  std::vector<double> C2 = std::vector<double>(),
                  std::vector<double> C3 = std::vector<double>()
                  )
    : Tc(Tc), pc(pc), acentric(acentric), R_u(R_u), Delta_1(Delta_1), Delta_2(Delta_2), C1(C1), C2(C2), C3(C3)
    {
        N = static_cast<int>(Tc.size());
        k.resize(N, std::vector<double>(N, 0));
        cm = 0.;
        /// If no Mathias-Copeman coefficients are passed in (all empty vectors), use the predictive scheme for m_ii
        simple_aii = (C1.empty() && C2.empty() && C3.empty() && L_Twu.empty() && M_Twu.empty() && N_Twu.empty());
    };
    /// Set the kij factor for the ij pair
    void set_kij(std::size_t i, std::size_t j, double val){ k[i][j] = val; }
    /// Get the kij factor for the ij pair
    double get_kij(std::size_t i, std::size_t j){ return k[i][j]; }
    /// Get the vector of critical temperatures (in K)
    std::vector<double> &get_Tc(){ return Tc; }
    /// Get the vector of critical pressures (in Pa)
    std::vector<double> &get_pc(){ return pc; }
    /// Get the vector of acentric factors
    std::vector<double> &get_acentric(){ return acentric; }
    /// Read-only accessor for value of Delta_1
    double get_Delta_1(){ return Delta_1; }
    /// Read-only accessor for value of Delta_2
    double get_Delta_2(){ return Delta_2; }
    /// Read-only accessor for value of R_u (universal gas constant)
    double get_R_u(){ return R_u; }
    /// Get a constant reference to a constant vector of Mathias-Copeman constants
    const std::vector<double> &get_C_ref(int i){
        switch (i){
            case 1: return C1;
            case 2: return C2;
            case 3: return C3;
            default:
                throw -1;
        }
    }
    /// Get a constant reference to a constant vector of Mathias-Copeman constants
    std::vector<double> &get_C(int i){
        switch (i){
            case 1: return C1;
            case 2: return C2;
            case 3: return C3;
            default:
                throw -1;
        }
    }
    /// Set the three Mathias-Copeman constants in one shot for a pure fluid
    void set_C_MC(double c1, double c2, double c3){
        C1.resize(1); C2.resize(1); C3.resize(1);
        C1[0] = c1, C2[0] = c2; C3[0] = c3;
        simple_aii = false; aii_model = AII_MATHIAS_COPEMAN;
    }
    /// Set the three Twu constants in one shot for a pure fluid
    void set_C_Twu(double L, double M, double N){
        L_Twu.resize(1); L_Twu[0] = L;
        M_Twu.resize(1); M_Twu[0] = M;
        N_Twu.resize(1); N_Twu[0] = N;
        simple_aii = false; aii_model = AII_TWU;
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
    virtual double alphar(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta);
    /// The first composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d_alphar_dxi(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
    /// The second composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d2_alphar_dxidxj(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /// The third composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    virtual double d3_alphar_dxidxjdxk(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /**
     * \brief The n-th derivative of \f$a_m\f$ with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     */
    virtual double am_term(double tau, const std::vector<double> &x, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_am_term_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_am_term_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
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
    virtual double d3_am_term_dxidxjdxk(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /**
     * \brief The term \f$b_{\rm m}\f$ (mixture co-volume)
     * \param x The vector of mole fractions
     */
    virtual double bm_term(const std::vector<double> &x);
    /** \brief The first composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_bm_term_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent);
    /** \brief The second composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_bm_term_dxidxj(const std::vector<double> &x, std::size_t i, std::size_t j, bool xN_independent);
    /** \brief The third composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d3_bm_term_dxidxjdxk(const std::vector<double> &x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

protected:
    /**
    * \brief The term \f$c_{\rm m}\f$ (volume translation)
    */
    virtual double cm_term();
    /// Set the volume translation parameter
    void set_cm(double val) {cm = val; }

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
    double psi_minus(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta);
    /**
     * \brief The third composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_psi_minus_dxi(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
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
    double d2_psi_minus_dxidxj(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
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
    double d3_psi_minus_dxidxjdxk(double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /**
     * \brief The term \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[ \Pi_{12} = (1+\Delta_1\bm\rhor \delta)(1+\Delta_2\bm\rhor \delta) \f]
     */
    double PI_12(double delta, const std::vector<double> &x, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_PI_12_dxi(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_PI_12_dxidxj(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
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
    double d3_PI_12_dxidxjdxk(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /**
     * \brief The term \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     */
    double tau_times_a(double tau, const std::vector<double> &x, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_tau_times_a_dxi(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_tau_times_a_dxidxj(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
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
    double d3_tau_times_a_dxidxjdxk(double tau, const std::vector<double> &x, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /**
     * \brief The term \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[  \psi^{(+)} = \dfrac{\ln\left(\dfrac{\Delta_1\bm\rhor \delta+1}{\Delta_2\bm\rhor \delta+1}\right)}{\bm(\Delta_1-\Delta_2)}  \f]
     */
    double psi_plus(double delta, const std::vector<double> &x, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_psi_plus_dxi(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_psi_plus_dxidxj(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
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
    double d3_psi_plus_dxidxjdxk(double delta, const std::vector<double> &x, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);
    
    /** \brief The term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     *
     * \f$c\f$ is given by
     * \f[
     * c = \frac{1}{b_m}
     * \f]
     * \param x The vector of mole fractions
     */
    double c_term(const std::vector<double> &x){
        return 1/bm_term(x);
    };
    /**
     * \brief The first composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_c_term_dxi(const std::vector<double> &x, std::size_t i, bool xN_independent){
        return -d_bm_term_dxi(x,i,xN_independent)/pow(bm_term(x), 2);
    };
    /**
     * \brief The second composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_c_term_dxidxj(const std::vector<double> &x, std::size_t i, std::size_t j, bool xN_independent){
        double bm = bm_term(x);
        return (2*d_bm_term_dxi(x, i, xN_independent)*d_bm_term_dxi(x, j, xN_independent) - bm*d2_bm_term_dxidxj(x, i,j,xN_independent))/pow(bm, 3);
    };
    /**
     * \brief The third composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_c_term_dxidxjdxk(const std::vector<double> &x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent){
        double bm = bm_term(x);
        return 1/pow(bm,4)*(2*bm*(d_bm_term_dxi(x, i, xN_independent)*d2_bm_term_dxidxj(x, j, k, xN_independent)
                                +d_bm_term_dxi(x, j, xN_independent)*d2_bm_term_dxidxj(x, i, k, xN_independent)
                                +d_bm_term_dxi(x, k, xN_independent)*d2_bm_term_dxidxj(x, i, j, xN_independent)
                                )
                           - pow(bm,2)*d3_bm_term_dxidxjdxk(x, i,j,k,xN_independent)
                           -6*d_bm_term_dxi(x, i, xN_independent)*d_bm_term_dxi(x, j, xN_independent)*d_bm_term_dxi(x, k, xN_independent)
                           );
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
    double A_term(double delta, const std::vector<double> &x){
        double bm = bm_term(x);
        double cm = cm_term();
        return log((delta*rho_r*(Delta_1*bm+cm)+1)/(delta*rho_r*(Delta_2*bm + cm) +1));
    };
    /**
     * \brief The first composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_A_term_dxi(double delta, const std::vector<double> &x, std::size_t i, bool xN_independent){
        std::size_t idelta = 0;
        return delta*rho_r*d_bm_term_dxi(x,i,xN_independent)*(Delta_1-Delta_2)/PI_12(delta, x, idelta);
    };
    /**
     * \brief The second composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_A_term_dxidxj(double delta, const std::vector<double> &x, std::size_t i, std::size_t j, bool xN_independent){
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, x, idelta);
        return delta*rho_r*(Delta_1-Delta_2)/pow(PI12, 2)*(PI12*d2_bm_term_dxidxj(x,i,j,xN_independent)
                                                           - d_PI_12_dxi(delta, x, 0, j,xN_independent)* d_bm_term_dxi(x,i,xN_independent));
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
    double d3_A_term_dxidxjdxk(double delta, const std::vector<double> &x, std::size_t i, std::size_t j, std::size_t k, bool xN_independent){
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, x, idelta);
        // The leading factor
        double lead = delta*rho_r*(Delta_1-Delta_2)/pow(PI12, 3);
        return lead*(-PI12*(d_PI_12_dxi(delta, x, idelta, j, xN_independent)*d2_bm_term_dxidxj(x,i,k,xN_independent)
                            +d_PI_12_dxi(delta, x, idelta, k, xN_independent)*d2_bm_term_dxidxj(x,i,j,xN_independent)
                            +d_bm_term_dxi(x,i,xN_independent)*d2_PI_12_dxidxj(delta, x, idelta, j, k, xN_independent))
                     +pow(PI12, 2)*d3_bm_term_dxidxjdxk(x, i, j, k, xN_independent)
                     +2*d_PI_12_dxi(delta, x, idelta, j, xN_independent)*d_PI_12_dxi(delta, x, idelta, k, xN_independent)*d_bm_term_dxi(x,i, xN_independent)
                     );
    };
};

class PengRobinson : public AbstractCubic
{
public:
    PengRobinson(std::vector<double> Tc,
                 std::vector<double> pc,
                 std::vector<double> acentric,
                 double R_u,
                 std::vector<double> C1 = std::vector<double>(),
                 std::vector<double> C2 = std::vector<double>(),
                 std::vector<double> C3 = std::vector<double>()
                 )
    : AbstractCubic(Tc, pc, acentric, R_u, 1+sqrt(2.0), 1-sqrt(2.0),C1,C2,C3) {};
    
    PengRobinson(double Tc,
                 double pc,
                 double acentric,
                 double R_u)
    : AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1+sqrt(2.0), 1-sqrt(2.0)) {};
    
    virtual double a0_ii(std::size_t i);
    virtual double b0_ii(std::size_t i);
    virtual double m_ii(std::size_t i);
};

class SRK : public AbstractCubic
{
public:
    SRK(std::vector<double> Tc,
        std::vector<double> pc,
        std::vector<double> acentric,
        double R_u,
        std::vector<double> C1 = std::vector<double>(),
        std::vector<double> C2 = std::vector<double>(),
        std::vector<double> C3 = std::vector<double>())
    : AbstractCubic(Tc, pc, acentric, R_u, 1, 0, C1, C2, C3) {};
    SRK(double Tc,
        double pc,
        double acentric,
        double R_u)
    : AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1, 0) {};
    
    virtual double a0_ii(std::size_t i);
    virtual double b0_ii(std::size_t i);
    virtual double m_ii(std::size_t i);
};

#endif