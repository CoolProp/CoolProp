#ifndef SBTLBACKEND_H
#define SBTLBACKEND_H

#include "TabularBackends.h"
#include "Exceptions.h"
#include "DataStructures.h"

namespace CoolProp {

// ---------------------------------------------------------------------------
// DU-space table subclasses — forward evaluations of T(D,U) and P(D,U)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// DT-space table subclasses — forward evaluations of P(D,T), H(D,T), S(D,T)
// ---------------------------------------------------------------------------

/**
 * @brief Liquid-phase DT table: D ∈ [D_crit, D_L(T_triple)×1.10], T ∈ [T_triple, T_max×1.499]
 *
 * D is linearly spaced (liquid range ~3–4×).
 * P(D,T), H(D,T), S(D,T), U(D,T) are smooth functions of (D,T) in the single-phase
 * liquid region — no saturation-dome crossing, no inversion needed.
 */
class DTLiquidTable : public SinglePhaseGriddedTableData
{
   public:
    DTLiquidTable() {
        xkey = iDmolar;
        ykey = iT;
        logx = false;  // D ratio liquid: ~3–4×, linear OK
        logy = false;
    }
    void set_limits() override;
};

/**
 * @brief Gas-phase DT table: D ∈ [D_V(T_triple)×0.999, D_crit], T ∈ [T_triple, T_max×1.499]
 *
 * D is log-spaced because the gas density range spans many orders of magnitude.
 * P(D,T), H(D,T), S(D,T), U(D,T) are smooth functions in the single-phase gas region.
 */
class DTGasTable : public SinglePhaseGriddedTableData
{
   public:
    DTGasTable() {
        xkey = iDmolar;
        ykey = iT;
        logx = true;   // D ratio gas: large range, log spacing essential
        logy = false;
    }
    void set_limits() override;
};

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
 * Interpolation scheme:
 *   pT and pH tables: bi-cubic (degree-3) 4-point Lagrange on a 4×4 stencil.
 *     16 coefficients per cell, alpha[m*4+n] for m,n ∈ {0,1,2,3}.
 *     Accuracy comparable to BICUBIC.  Inversion (e.g. find xhat given z_target
 *     and yhat) solves a cubic by Newton from a quadratic seed — 2–3 iterations.
 *   DU tables: bi-quadratic (degree-2) 3-point Lagrange on a 3×3 stencil.
 *     9 coefficients per cell, alpha[m*3+n] for m,n ∈ {0,1,2}.
 *     Forward evaluation only — no inversion needed for D,U flash.
 *
 * D,U flash: instead of inverting the pT table (requiring Newton iteration),
 * dedicated DU-space tables store T(D,U) and P(D,U) directly, so flash is a
 * single cell-lookup + polynomial evaluation — O(1) with no iteration.
 */
class SBTLBackend : public TabularBackend
{
   public:
    /// Constructor: loads or builds tables, then computes bi-quadratic coefficients
    SBTLBackend(shared_ptr<CoolProp::AbstractState> AS)
      : TabularBackend(AS), du_tables_built(false), sat_cache_built(false),
        dt_tables_built(false), _dt_flash_active(false),
        _dt_tbl_ptr(nullptr), _dt_coeff_ptr(nullptr) {
        imposed_phase_index = iphase_not_imposed;
        if (!this->AS->get_mole_fractions().empty()) {
            check_tables();
            build_sbtl_coeffs(dataset->single_phase_logph, _sbtl_coeffs_ph, 4);
            build_sbtl_coeffs(dataset->single_phase_logpT, _sbtl_coeffs_pT, 4);
            is_mixture = (this->AS->get_mole_fractions().size() > 1);
            if (!is_mixture) {
                build_du_tables();
                build_dt_tables();
            }
        }
    }

    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        this->AS->set_mole_fractions(mole_fractions);
        is_mixture = true;
        check_tables();
        build_sbtl_coeffs(dataset->single_phase_logph, _sbtl_coeffs_ph, 4);
        build_sbtl_coeffs(dataset->single_phase_logpT, _sbtl_coeffs_pT, 4);
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
     * pT and pH tables use nstencil=4 (degree-3, 16 coefficients) for BICUBIC-comparable
     * accuracy.  DU tables use nstencil=3 (degree-2, 9 coefficients) since they only
     * need forward evaluation and their density axes are well-conditioned.
     */
    void build_sbtl_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, int nstencil = 3);

    // -----------------------------------------------------------------------
    // Evaluation
    // -----------------------------------------------------------------------

    /// Core evaluation: bi-quadratic polynomial at (x,y) in cell (i,j)
    double evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                 const parameters output, const double x, const double y, const std::size_t i, const std::size_t j);

    double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logph, _sbtl_coeffs_ph, output, _hmolar, _p, i, j);
    }
    double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logpT, _sbtl_coeffs_pT, output, _T, _p, i, j);
    }

    // -----------------------------------------------------------------------
    // Derivatives
    // -----------------------------------------------------------------------

    double evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                                            double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);

    double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logph, _sbtl_coeffs_ph, output, _hmolar, _p, i, j, Nx, Ny);
    }
    double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logpT, _sbtl_coeffs_pT, output, _T, _p, i, j, Nx, Ny);
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

    // -----------------------------------------------------------------------
    // D,T flash: direct DT-table forward evaluation (no inversion)
    // -----------------------------------------------------------------------

    /// Override the base-class DmolarT flash to use dedicated DT-space polynomial
    /// tables where P(D,T), H(D,T), S(D,T), U(D,T) are direct forward evaluations.
    /// This avoids the ill-conditioned D(T,P)→P inversion near the saturation boundary.
    void flash_DmolarT(CoolPropDbl D, CoolPropDbl T);

    /// Reset _dt_flash_active before every update so stale values never persist.
    void update(CoolProp::input_pairs input_pair, double val1, double val2);

    // -----------------------------------------------------------------------
    // Property accessors that respect the DT-flash cache
    // -----------------------------------------------------------------------
    /// When _dt_flash_active is set (DmolarT flash was last), return the value
    /// cached directly from the DT table instead of evaluating the pT table,
    /// which might point to an extrapolated cell after a DT→pT mapping.
    CoolPropDbl calc_hmolar(void);
    CoolPropDbl calc_smolar(void);
    CoolPropDbl calc_umolar(void);
    CoolPropDbl calc_cpmolar(void);
    CoolPropDbl calc_cvmolar(void);
    CoolPropDbl calc_rhomolar(void);

   private:
    // SBTL's own pT and pH polynomial coefficients — stored separately from
    // dataset->coeffs_ph/pT (which BICUBIC also uses) to prevent conflicts when
    // both backends share the same TabularDataSet.
    std::vector<std::vector<CellCoeffs>> _sbtl_coeffs_ph;
    std::vector<std::vector<CellCoeffs>> _sbtl_coeffs_pT;

    // Helper: return the appropriate SBTL coefficient array for the given table,
    // falling back to `default_coeffs` for non-pT/pH tables (e.g. DU tables).
    const std::vector<std::vector<CellCoeffs>>& sbtl_coeffs_for(
        const SinglePhaseGriddedTableData& table,
        const std::vector<std::vector<CellCoeffs>>& default_coeffs) const {
        if (&table == &dataset->single_phase_logph) return _sbtl_coeffs_ph;
        if (&table == &dataset->single_phase_logpT) return _sbtl_coeffs_pT;
        return default_coeffs;
    }

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

    // Dedicated (D, T)-space tables for direct P(D,T), H(D,T), S(D,T), U(D,T) evaluation.
    // Avoids the ill-conditioned inversion D(T,P)→P near the saturation boundary.
    DTLiquidTable dt_liquid;
    DTGasTable    dt_gas;
    std::vector<std::vector<CellCoeffs>> coeffs_dt_liquid;
    std::vector<std::vector<CellCoeffs>> coeffs_dt_gas;
    bool dt_tables_built;

    // State cache for flash_DmolarT: set true after a DT flash, cleared on next update().
    // When active, calc_hmolar/smolar/umolar/rhomolar return cached DT values directly.
    bool   _dt_flash_active;
    double _dt_hmolar;  // cached H from DT table [J/mol]
    double _dt_smolar;  // cached S from DT table [J/mol/K]
    double _dt_umolar;  // cached U from DT table [J/mol]
    std::size_t _dt_i;  // DT table cell row index
    std::size_t _dt_j;  // DT table cell column index
    SinglePhaseGriddedTableData*              _dt_tbl_ptr;    // pointer to active DT table
    std::vector<std::vector<CellCoeffs>>*     _dt_coeff_ptr;  // pointer to active DT coeffs

    /// Build the liquid and gas DT tables and compute their SBTL coefficients.
    void build_dt_tables();
};

}  // namespace CoolProp

#endif  // SBTLBACKEND_H
