#ifndef SBTLBACKEND_H
#define SBTLBACKEND_H

#include "TabularBackends.h"
#include "Exceptions.h"
#include "DataStructures.h"

namespace CoolProp {

// ---------------------------------------------------------------------------
// DU-space table subclasses — forward evaluations of T(D,U) and P(D,U)
// ---------------------------------------------------------------------------

/**
 * @brief Liquid-phase DU table: D ∈ [D_crit, D_L(T_triple)], U ∈ [U_L(T_triple), U(T_max, D_crit)]
 *
 * Both axes are linearly spaced (liquid density range is only ~3×).
 * The single-phase liquid region occupies the lower-right triangle of this rectangle;
 * two-phase and out-of-range cells are left as _HUGE holes by build() and are marked
 * invalid by build_sbtl_coeffs().
 */
class DULiquidTable : public SinglePhaseGriddedTableData
{
   public:
    DULiquidTable() {
        xkey = iDmolar;
        ykey = iUmolar;
        logx = false;  // D ratio liquid: ~3–4×, linear OK
        logy = false;
    }
    void set_limits() override;
};

/**
 * @brief Gas-phase DU table: D ∈ [D_V(T_triple), D_crit], U ∈ [U_V(T_triple), U(T_max, D_V)]
 *
 * D is log-spaced because the gas density range spans ~5 orders of magnitude
 * (from ~0.27 mol/m³ at the triple point to ~17873 mol/m³ at the critical point for water).
 */
class DUGasTable : public SinglePhaseGriddedTableData
{
   public:
    DUGasTable() {
        xkey = iDmolar;
        ykey = iUmolar;
        logx = true;   // D ratio gas: ~66000× for water, log spacing essential
        logy = false;
    }
    void set_limits() override;
};

/**
 * @brief SBTLBackend implements the Spline-Based Table Look-Up (SBTL) method
 * for fast thermodynamic property evaluation.
 *
 * Based on: Kruse & Knobloch (IAPWS standard).
 * Interpolation scheme: bi-quadratic (degree 2 × 2) Lagrange polynomials on a
 * 3×3 stencil of node values.  This gives 9 coefficients per cell per property,
 * stored in the existing CellCoeffs::std::vector<double> members.
 *
 * Coefficient layout (row-major, 3×3):
 *   alpha[m*3 + n]  corresponds to  a_{mn}  in  z = Σ_{m,n=0}^{2} a_{mn} xhat^m yhat^n
 * where xhat,yhat ∈ [0,1] map to the cell interval [x_i, x_{i+1}] × [y_j, y_{j+1}].
 *
 * Inversion (e.g. find xhat given z_target and yhat) reduces to a 1-D quadratic
 * equation that is solved analytically — no Newton iteration needed.
 *
 * D,U flash: instead of inverting the pT table (requiring Newton iteration),
 * dedicated DU-space tables store T(D,U) and P(D,U) directly, so flash is a
 * single cell-lookup + polynomial evaluation — O(1) with no iteration.
 */
class SBTLBackend : public TabularBackend
{
   public:
    /// Constructor: loads or builds tables, then computes bi-quadratic coefficients
    SBTLBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend(AS), du_tables_built(false), sat_cache_built(false) {
        imposed_phase_index = iphase_not_imposed;
        if (!this->AS->get_mole_fractions().empty()) {
            check_tables();
            build_sbtl_coeffs(dataset->single_phase_logph, dataset->coeffs_ph);
            build_sbtl_coeffs(dataset->single_phase_logpT, dataset->coeffs_pT);
            is_mixture = (this->AS->get_mole_fractions().size() > 1);
            if (!is_mixture) {
                build_du_tables();
            }
        }
    }

    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        this->AS->set_mole_fractions(mole_fractions);
        is_mixture = true;
        check_tables();
        build_sbtl_coeffs(dataset->single_phase_logph, dataset->coeffs_ph);
        build_sbtl_coeffs(dataset->single_phase_logpT, dataset->coeffs_pT);
        // DU tables are not built for mixtures
    }

    std::string backend_name(void) {
        return get_backend_string(SBTL_BACKEND);
    }

    /**
     * @brief Build bi-quadratic SBTL coefficients from node values.
     *
     * For each interior cell (i, j) the stencil uses nodes (i-1..i+1, j-1..j+1).
     * Boundary cells are clamped and marked invalid (with a pointer to the nearest
     * valid interior cell), matching the CellCoeffs alt_i/alt_j pattern.
     *
     * The 9 coefficients a_{mn} (m,n ∈ {0,1,2}) stored row-major at alpha[m*3+n]
     * define the polynomial z(xhat,yhat) = Σ a_{mn} xhat^m yhat^n with
     * xhat = (x-x_i)/(x_{i+1}-x_i), yhat = (y-y_j)/(y_{j+1}-y_j).
     */
    void build_sbtl_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs);

    // -----------------------------------------------------------------------
    // Evaluation
    // -----------------------------------------------------------------------

    /// Core evaluation: bi-quadratic polynomial at (x,y) in cell (i,j)
    double evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                 const parameters output, const double x, const double y, const std::size_t i, const std::size_t j);

    double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logph, dataset->coeffs_ph, output, _hmolar, _p, i, j);
    }
    double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logpT, dataset->coeffs_pT, output, _T, _p, i, j);
    }

    // -----------------------------------------------------------------------
    // Derivatives
    // -----------------------------------------------------------------------

    double evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                                            double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);

    double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logph, dataset->coeffs_ph, output, _hmolar, _p, i, j, Nx, Ny);
    }
    double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logpT, dataset->coeffs_pT, output, _T, _p, i, j, Nx, Ny);
    }

    // -----------------------------------------------------------------------
    // Transport properties (bilinear interpolation)
    // -----------------------------------------------------------------------

    double evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j);

    double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase_transport(dataset->single_phase_logph, output, _hmolar, _p, i, j);
    }
    double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase_transport(dataset->single_phase_logpT, output, _T, _p, i, j);
    }

    // -----------------------------------------------------------------------
    // Cell finding
    // -----------------------------------------------------------------------

    void find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                          double x, double y, std::size_t& i, std::size_t& j);

    void find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                               const parameters variable1, const double value1, const parameters otherkey, const double otherval,
                               std::size_t& i, std::size_t& j);

    // -----------------------------------------------------------------------
    // Inversion (analytical quadratic formula)
    // -----------------------------------------------------------------------

    void invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                               parameters other_key, double other, double y, std::size_t i, std::size_t j);

    void invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                               parameters other_key, double other, double x, std::size_t i, std::size_t j);

    // -----------------------------------------------------------------------
    // D,U flash: direct DU-table lookup (no iteration)
    // -----------------------------------------------------------------------

    /// Override the base-class DmolarUmolar flash to use dedicated DU-space polynomial
    /// tables so that T(D,U) and P(D,U) are direct forward evaluations — no iteration.
    void flash_DmolarUmolar(CoolPropDbl D, CoolPropDbl U);

   private:
    // Dedicated (D, U)-space tables for direct T(D,U) and P(D,U) evaluation.
    // Separate tables for liquid (D > D_crit) and gas (D < D_crit) avoid the need
    // for phase-aware switching inside evaluate_single_phase.
    DULiquidTable du_liquid;
    DUGasTable    du_gas;
    std::vector<std::vector<CellCoeffs>> coeffs_du_liquid;
    std::vector<std::vector<CellCoeffs>> coeffs_du_gas;
    bool du_tables_built;

    // Saturation property cache for flash_DmolarUmolar phase detection.
    // Replaces repeated sat_eval (EOS) calls with O(1) array lookups.
    // Arrays are log-spaced in pressure from p_triple to p_critical.
    // Per SBTL PDF Appendix A4, these enable fast Newton iteration for two-phase states.
    std::vector<double> _sat_logp;  // log(p) values, uniformly spaced
    std::vector<double> _sat_DL;    // saturation liquid density  [mol/m^3]
    std::vector<double> _sat_DV;    // saturation vapor density   [mol/m^3]
    std::vector<double> _sat_UL;    // saturation liquid U        [J/mol]
    std::vector<double> _sat_UV;    // saturation vapor U         [J/mol]
    std::vector<double> _sat_T;     // saturation temperature     [K]
    double _sat_log_p_min;          // log(p_triple)
    double _sat_inv_dlogp;          // 1 / ((log(p_crit) - log(p_triple)) / (N-1))
    double _sat_DL_triple;          // cached saturation liquid density at triple point
    double _sat_DV_triple;          // cached saturation vapor density at triple point
    bool sat_cache_built;

    /// Build the liquid and gas DU tables and compute their SBTL coefficients.
    void build_du_tables();

    /// Build the saturation property cache (called from build_du_tables).
    void build_sat_cache();

    /// O(1) saturation liquid property lookup at pressure p using the precomputed cache.
    double fast_sat_lookup(parameters param, double p) const;
    /// O(1) saturation vapor property lookup at pressure p using the precomputed cache.
    double fast_sat_lookup_vapor(parameters param, double p) const;
};

}  // namespace CoolProp

#endif  // SBTLBACKEND_H
