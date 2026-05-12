#ifndef SBTLBACKEND_H
#define SBTLBACKEND_H

#include "TabularBackends.h"
#include "Exceptions.h"
#include "DataStructures.h"
#include <optional>

namespace CoolProp {

/**
 * @brief Single Chebyshev expansion in log p on one interval.
 */
class Cheb1DPiece
{
   public:
    Cheb1DPiece() = default;
    Cheb1DPiece(double log_p_lo, double log_p_hi, std::vector<double> coeffs) : log_p_lo_(log_p_lo), log_p_hi_(log_p_hi), c_(std::move(coeffs)) {}

    [[nodiscard]] bool valid() const {
        return c_.size() >= 2;
    }

    template <typename Fn>
    static Cheb1DPiece build(double p_lo, double p_hi, std::size_t N, Fn&& f) {
        const double a = std::log(p_lo), b = std::log(p_hi);
        std::vector<double> values(N + 1);
        for (std::size_t k = 0; k <= N; ++k) {
            const double theta = M_PI * static_cast<double>(k) / static_cast<double>(N);
            const double xk = std::cos(theta);
            const double log_pk = 0.5 * (a + b) + 0.5 * (b - a) * xk;
            values[k] = f(std::exp(log_pk));
        }
        std::vector<double> coeffs(N + 1, 0.0);
        for (std::size_t k = 0; k <= N; ++k) {
            double s = 0.0;
            for (std::size_t j = 0; j <= N; ++j) {
                const double wj = (j == 0 || j == N) ? 0.5 : 1.0;
                s += wj * values[j] * std::cos(M_PI * static_cast<double>(j) * static_cast<double>(k) / static_cast<double>(N));
            }
            coeffs[k] = (2.0 / static_cast<double>(N)) * s;
        }
        coeffs[0] *= 0.5;
        coeffs[N] *= 0.5;
        return Cheb1DPiece(a, b, std::move(coeffs));
    }

    [[nodiscard]] double eval(double p) const {
        const double log_p = std::log(p);
        double t = (2.0 * log_p - (log_p_lo_ + log_p_hi_)) / (log_p_hi_ - log_p_lo_);
        if (t < -1.0) t = -1.0;
        if (t > 1.0) t = 1.0;
        const std::size_t N = c_.size() - 1;
        double bk1 = 0.0, bk2 = 0.0;
        for (std::size_t k = N; k > 0; --k) {
            const double bk = c_[k] + 2.0 * t * bk1 - bk2;
            bk2 = bk1;
            bk1 = bk;
        }
        return c_[0] + t * bk1 - bk2;
    }

    [[nodiscard]] double log_p_hi() const {
        return log_p_hi_;
    }

   private:
    double log_p_lo_{0.0};
    double log_p_hi_{0.0};
    std::vector<double> c_;
};

/**
 * @brief Piecewise Chebyshev expansion in log p, used to interpolate
 * the per-isobar boundary functions h_lo(p), h_hi(p) of NormalizedPHTable.
 *
 * Subdivision lets the bulk piece converge spectrally on the smooth body
 * of the function while a near-critical piece absorbs the sqrt-cusp at
 * p_crit (h_sat,L(p) and h_sat,V(p) both terminate with vertical tangent
 * because the analytic Helmholtz EOS has β = 1/2).  Without subdivision,
 * a single Cheb covering all the way to p_crit converges only
 * algebraically and the contaminated polynomial Runge-oscillates at
 * lower p, ruining accuracy mid-table.
 *
 * Pieces are built tail-to-head: caller passes piece breakpoints in
 * ascending p; each piece's domain is contiguous with the next.
 */
class Cheb1D
{
   public:
    Cheb1D() = default;

    [[nodiscard]] bool valid() const {
        return !pieces_.empty() && pieces_.front().valid();
    }

    /// Build piecewise Cheb on the supplied breakpoints {p_0, p_1, ..., p_n}
    /// (sorted ascending), with N coeffs per piece.  Produces n pieces:
    /// [p_0, p_1], [p_1, p_2], ..., [p_{n-1}, p_n].
    template <typename Fn>
    static Cheb1D build(const std::vector<double>& breakpoints, std::size_t N, Fn&& f) {
        Cheb1D out;
        if (breakpoints.size() < 2) return out;
        out.pieces_.reserve(breakpoints.size() - 1);
        for (std::size_t k = 0; k + 1 < breakpoints.size(); ++k) {
            out.pieces_.push_back(Cheb1DPiece::build(breakpoints[k], breakpoints[k + 1], N, f));
        }
        return out;
    }

    /// Evaluate at p, dispatching to the piece whose domain contains p.
    /// Outside the global domain, the closest piece is used (with internal
    /// clamp on its parameter).
    [[nodiscard]] double eval(double p) const {
        const double log_p = std::log(p);
        // Linear scan — typically ≤ 4 pieces, faster than a bisect.
        for (std::size_t k = 0; k + 1 < pieces_.size(); ++k) {
            if (log_p <= pieces_[k].log_p_hi()) return pieces_[k].eval(p);
        }
        return pieces_.back().eval(p);
    }

   private:
    std::vector<Cheb1DPiece> pieces_;
};

/**
 * @brief Normalized-h × log-P table — coordinate-aligned PH replacement.
 *
 * Three flavours, distinguished by the `region_` tag:
 *   LIQUID:  P ∈ [P_min, P_crit), h ∈ [h_min(P), h_sat,L(P)],
 *            xnorm = (h - h_min(P)) / (h_sat,L(P) - h_min(P))
 *   VAPOR:   P ∈ [P_min, P_crit), h ∈ [h_sat,V(P), h_max(P)],
 *            xnorm = (h - h_sat,V(P)) / (h_max(P) - h_sat,V(P))
 *   SUPER:   P ∈ [P_crit, P_max],  h ∈ [h_min(P), h_max(P)],
 *            xnorm = (h - h_min(P)) / (h_max(P) - h_min(P))
 *
 * The grid is uniform in xnorm × log(P). xkey = iHmolar (still: corner values
 * are real h's, computed via h_from_xnorm at table-build time and stored).
 * h_lo_isobar and h_hi_isobar persist the per-row enthalpy bounds so xnorm_from_h
 * at lookup matches the build-time encoding. h_lo and h_hi values are populated
 * by NormalizedPHTable::populate_isobar_bounds before SinglePhaseGriddedTableData::build
 * runs.  Construction requires h_sat values for LIQUID/VAPOR — supplied by an
 * SBTLBackend* helper at populate time.
 */
class NormalizedPHTable : public SinglePhaseGriddedTableData
{
   public:
    enum Region : uint8_t
    {
        LIQUID = 0,
        VAPOR = 1,
        SUPER = 2
    };

    explicit NormalizedPHTable(Region region) : region_(region) {
        xkey = iHmolar;
        ykey = iP;
        logx = false;  // xnorm is uniform in [0, 1]
        logy = true;   // log-spaced in pressure, same as LogPHTable
    }

    void set_limits() override;

    /// Forward transform: given a query (h, P) inside this region, compute
    /// the normalized x-coordinate.  Optional `h_sat` arg is used for the
    /// region's saturation boundary (h_hi for LIQUID, h_lo for VAPOR);
    /// passing the H-superancillary value here keeps routing and lookup
    /// consistent.  If NaN, falls back to the Cheb evaluation.
    [[nodiscard]] double xnorm_from_h(double h, double P, double h_sat = std::numeric_limits<double>::quiet_NaN()) const;

    /// Inverse transform: given (xnorm ∈ [0,1], P), produce h.  Same
    /// optional h_sat argument as xnorm_from_h.
    [[nodiscard]] double h_from_xnorm(double xnorm, double P, double h_sat = std::numeric_limits<double>::quiet_NaN()) const;

    [[nodiscard]] Region region() const noexcept {
        return region_;
    }

    /// Per-Ny-row enthalpy bounds, sampled at build time:
    ///   LIQUID: h_lo = h(T_min, P), h_hi = h_sat,L(P)
    ///   VAPOR:  h_lo = h_sat,V(P), h_hi = h(T_max, P) (or pmax-isobar high-T)
    ///   SUPER:  h_lo = h(T_min, P), h_hi = h(T_max, P)
    /// Persisted via msgpack so lookup-time xnorm matches build-time encoding.
    std::vector<double> h_lo_isobar, h_hi_isobar;

    /// Spectral expansion of h_lo(p) and h_hi(p) on [log p_min, log p_max].
    /// xnorm_from_h evaluates these at lookup; h_lo_isobar/h_hi_isobar
    /// are still populated and used at table build time for cell-fill.
    Cheb1D h_lo_cheb, h_hi_cheb;

   private:
    Region region_;
};

/**
 * @brief Normalized-T × log-P table — coordinate-aligned PT replacement.
 *
 * Three flavours, distinguished by the `region_` tag:
 *   LIQUID:  P ∈ [P_min, P_crit), T ∈ [T_lo(P), T_sat,L(P)],
 *            xnorm = (T - T_lo(P)) / (T_sat,L(P) - T_lo(P))
 *   VAPOR:   P ∈ [P_min, P_crit), T ∈ [T_sat,V(P), T_hi(P)],
 *            xnorm = (T - T_sat,V(P)) / (T_hi(P) - T_sat,V(P))
 *   SUPER:   P ∈ [P_crit, P_max],  T ∈ [T_lo(P), T_hi(P)],
 *            xnorm = (T - T_lo(P)) / (T_hi(P) - T_lo(P))
 *
 * T_lo(P) = T_melt(P) when the fluid has a melting ancillary, else T_min.
 * T_hi(P) = T_max * 1.499 (extended range).
 *
 * The grid is uniform in xnorm × log(P).  xkey = iT (corner values are real
 * temperatures, computed via T_from_xnorm at table-build time and stored).
 * T_lo_isobar / T_hi_isobar persist the per-row temperature bounds so
 * xnorm_from_T at lookup matches the build-time encoding.  Same coordinate
 * trick as NormalizedPHTable: the saturation curve becomes a coordinate
 * axis (xnorm = 1 for LIQUID, xnorm = 0 for VAPOR), so cells never straddle
 * the dome and the compressed-liquid PT cell-selection 5 % residue
 * documented in dev/sbtl_pt_outstanding_work.md disappears.
 */
class NormalizedPTTable : public SinglePhaseGriddedTableData
{
   public:
    enum Region : uint8_t
    {
        LIQUID = 0,
        VAPOR = 1,
        SUPER = 2
    };

    explicit NormalizedPTTable(Region region) : region_(region) {
        xkey = iT;
        ykey = iP;
        logx = false;  // xnorm is uniform in [0, 1]
        logy = true;   // log-spaced in pressure
    }

    void set_limits() override;

    /// Forward transform: given a query (T, P) inside this region, compute
    /// the normalized x-coordinate.  Optional T_sat arg supplies the
    /// region's saturation boundary (T_hi for LIQUID, T_lo for VAPOR);
    /// falls back to the Cheb evaluation when NaN.
    [[nodiscard]] double xnorm_from_T(double T, double P, double T_sat = std::numeric_limits<double>::quiet_NaN()) const;

    /// Inverse transform: given (xnorm ∈ [0,1], P), produce T.
    [[nodiscard]] double T_from_xnorm(double xnorm, double P, double T_sat = std::numeric_limits<double>::quiet_NaN()) const;

    [[nodiscard]] Region region() const noexcept {
        return region_;
    }

    /// Per-Ny-row temperature bounds, sampled at build time.
    /// LIQUID: T_lo = T_melt(P) (or T_min), T_hi = T_sat,L(P)
    /// VAPOR:  T_lo = T_sat,V(P),           T_hi = T_max_ext
    /// SUPER:  T_lo = T_melt(P) (or T_min), T_hi = T_max_ext
    std::vector<double> T_lo_isobar, T_hi_isobar;

    /// Spectral expansion of T_lo(p) and T_hi(p) over [log p_min, log p_max].
    /// xnorm_from_T evaluates these at lookup; T_lo_isobar / T_hi_isobar
    /// remain the source of truth at table build time so cell-build and
    /// cell-lookup encode the identical xnorm.
    Cheb1D T_lo_cheb, T_hi_cheb;

   private:
    Region region_;
};

/**
 * @brief SBTLBackend implements the Spline-Based Table Look-Up (SBTL) method
 * for fast thermodynamic property evaluation.
 *
 * Based on IAPWS Guideline G13-15 (Kunick et al.).  Pure fluids only; supports
 * PT_INPUTS, HmolarP_INPUTS, and HmassP_INPUTS as the sole input pairs.
 * Mixtures (real or pseudo-pure) are rejected at construction with a clean
 * error; users wanting tabular mixture support should use the BICUBIC or
 * TTSE backends instead.
 *
 * Interpolation backbone — coordinate-aligned PH and PT tables:
 *   NormalizedPHTable (LIQUID / VAPOR / SUPER): cubic B-spline with the
 *     saturation curve as a coordinate axis at xnorm=1 (LIQUID) /
 *     xnorm=0 (VAPOR).  Routed for HmolarP_INPUTS / HmassP_INPUTS.
 *   NormalizedPTTable (LIQUID / VAPOR / SUPER): same trick with T_lo(P),
 *     T_hi(P) replacing h_lo / h_hi.  Routed for PT_INPUTS in the
 *     subcritical and supercritical regions; the cusp at (T_c, p_c) and
 *     its surrounding shoulder fall back to HEOS direct via a 2D
 *     critical-region box.
 *
 * Derivative and transport accessors (cp, cv, speed of sound, viscosity,
 * conductivity, first_partial_deriv) route through HEOS direct after
 * priming the AbstractState at the table-resolved (h, p) or (T, p) state;
 * polynomial-derivative paths are not supported and throw.
 */
/// Free function: build a 16-coefficient Hermite bicubic alpha vector
/// from 16 corner data (4 values, 4 ∂f/∂xi, 4 ∂f/∂eta, 4 ∂²f/∂xi∂eta).
/// Derivatives must already be scaled to the unit-cell coordinates
/// (multiplied by Δxi=1, Δeta=1).  Returns the alpha vector in SBTL
/// convention: alpha[4m + n] = coefficient of xi^m·eta^n in the
/// reconstructed polynomial.  Corner indexing: 00=(0,0), 10=(1,0),
/// 01=(0,1), 11=(1,1).  Exposed for unit testing of CoolProp-foi.1.
std::vector<double> hermite_bicubic_polynomial_coeffs(double f00, double f10, double f01, double f11, double fx00, double fx10, double fx01,
                                                      double fx11, double fy00, double fy10, double fy01, double fy11, double fxy00, double fxy10,
                                                      double fxy01, double fxy11);

class SBTLBackend : public TabularBackend
{
   public:
    /// Constructor: loads or builds the underlying TabularDataSet and
    /// constructs the coordinate-aligned PH and PT spline backbones.
    /// Throws ValueError for mixtures or pseudo-pure fluids.
    SBTLBackend(shared_ptr<CoolProp::AbstractState> AS)
      : TabularBackend(AS),
        normph_tables_built(false),
        normpt_tables_built(false),
        _normph_liquid(NormalizedPHTable::LIQUID),
        _normph_vapor(NormalizedPHTable::VAPOR),
        _normph_super(NormalizedPHTable::SUPER),
        _normpt_liquid(NormalizedPTTable::LIQUID),
        _normpt_vapor(NormalizedPTTable::VAPOR),
        _normpt_super(NormalizedPTTable::SUPER) {
        imposed_phase_index = iphase_not_imposed;
        // Strict pure-fluid check.  Pseudo-pure fluids (R407C, R404A,
        // R407A, R410A, R507A) have mole_fractions.size() == 1 but a
        // non-zero glide band; we reject them here so the normalized
        // routing never has to deal with T_sat,L != T_sat,V.  Some
        // single-fluid backends (e.g. IF97) don't implement
        // get_mole_fractions; we treat them as pure by construction.
        bool composition_known = true;
        std::size_t n_components = 1;
        try {
            const auto& mole_fractions = this->AS->get_mole_fractions();
            n_components = mole_fractions.size();
        } catch (...) {
            // Backend doesn't expose composition (IF97, etc.); assume
            // single-component (the backend factory only routes pure
            // fluids to those backends in the first place).
            composition_known = false;
        }
        if (composition_known && n_components == 0) {
            // No fluid set yet — defer build.  This matches the legacy
            // TTSE/BICUBIC behaviour where the ctor doesn't allocate
            // anything until update() is called for the first time.
        } else {
            auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(this->AS.get());
            if (n_components > 1 || (heos != nullptr && !heos->is_pure())) {
                throw ValueError("SBTL backend only supports pure fluids "
                                 "(no mixtures or pseudo-pure fluids).  "
                                 "Use a different tabular backend "
                                 "(BICUBIC, TTSE) for mixtures.");
            }
            build_normph_tables();
            build_normpt_tables();
        }
    }

    std::string backend_name(void) {
        return get_backend_string(SBTL_BACKEND);
    }

    // -----------------------------------------------------------------------
    // Evaluation
    // -----------------------------------------------------------------------

    /// Core evaluation: bi-cubic polynomial at (xnorm, log P) in cell (i, j)
    /// of an active normalized PH or PT table.
    double evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                 const parameters output, const double x, const double y, const std::size_t i, const std::size_t j);

    // For these accessors, the most recent update() must have routed through
    // a normalized PH or PT table; subsequent property reads evaluate against
    // the same active table using the cached xnorm in place of the physical
    // h / T.  No legacy-table fallback exists — pure-fluid PT/PH only.

    double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) {
        if (_normph_active_table == nullptr || _normph_active_coeffs == nullptr) {
            throw ValueError("SBTL: no active normalized PH table — only PH input pairs route through SBTL.");
        }
        return evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, output, _normph_xnorm, _p, i, j);
    }
    double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) {
        if (_normpt_active_table == nullptr || _normpt_active_coeffs == nullptr) {
            throw ValueError("SBTL: no active normalized PT table — only PT input pairs route through SBTL.");
        }
        return evaluate_single_phase(*_normpt_active_table, *_normpt_active_coeffs, output, _normpt_xnorm, _p, i, j);
    }

    // -----------------------------------------------------------------------
    // Derivatives — SBTL does not support polynomial derivatives.
    // Property accessors (calc_cpmolar etc.) route through HEOS direct.
    // -----------------------------------------------------------------------

    double evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                                            double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);

    double evaluate_single_phase_phmolar_derivative(parameters /*output*/, std::size_t /*i*/, std::size_t /*j*/, std::size_t /*Nx*/,
                                                    std::size_t /*Ny*/) {
        throw ValueError("SBTL: polynomial derivatives not supported; use HEOS direct via calc_first_partial_deriv or the cp/cv/c_sound accessors");
    }
    double evaluate_single_phase_pT_derivative(parameters /*output*/, std::size_t /*i*/, std::size_t /*j*/, std::size_t /*Nx*/, std::size_t /*Ny*/) {
        throw ValueError("SBTL: polynomial derivatives not supported; use HEOS direct via calc_first_partial_deriv or the cp/cv/c_sound accessors");
    }

    // -----------------------------------------------------------------------
    // Transport properties — SBTL does not interpolate transport.
    // The override of calc_viscosity / calc_conductivity routes through HEOS.
    // -----------------------------------------------------------------------

    double evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j);

    double evaluate_single_phase_phmolar_transport(parameters /*output*/, std::size_t /*i*/, std::size_t /*j*/) {
        throw ValueError("SBTL: tabular transport not supported; use HEOS direct via viscosity()/conductivity() accessors");
    }
    double evaluate_single_phase_pT_transport(parameters /*output*/, std::size_t /*i*/, std::size_t /*j*/) {
        throw ValueError("SBTL: tabular transport not supported; use HEOS direct via viscosity()/conductivity() accessors");
    }

    // -----------------------------------------------------------------------
    // Cell finding — bisect on the normalized table's xvec/yvec.  All cells
    // on a normalized table are valid by construction (cell-build skips
    // entire isobars where the bounds aren't well-defined), so the only
    // failure mode is OOB coordinates, which the routing in update() guards.
    // -----------------------------------------------------------------------

    void find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, double x, double y,
                                          std::size_t& i, std::size_t& j);

    void find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, const parameters variable1,
                               const double value1, const parameters otherkey, const double otherval, std::size_t& i, std::size_t& j);

    // -----------------------------------------------------------------------
    // Inversion — not supported (SBTL only handles forward PT / PH queries).
    // -----------------------------------------------------------------------

    void invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters other_key,
                               double other, double y, std::size_t i, std::size_t j);

    void invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters other_key,
                               double other, double x, std::size_t i, std::size_t j);

    /// Override of update() to route HmolarP/HmassP through the coordinate-aligned
    /// NormalizedPHTable fast path (and the critical-region HEOS fallback box).
    void update(CoolProp::input_pairs input_pair, double val1, double val2);

    // -----------------------------------------------------------------------
    // Normalized-PH table support
    // -----------------------------------------------------------------------

    /// Populate the per-isobar enthalpy bounds (h_lo_isobar, h_hi_isobar) for a
    /// NormalizedPHTable. Must be called after the table's xvec/yvec grids are
    /// laid down (set_limits + base build's preamble) and before the cell-fill
    /// loop. Routes to AS for h(T_min, P) / h(T_max, P), and to the SBTL sat
    /// cache / H-superancillary for h_sat,L|V(p).
    ///
    /// LIQUID:   h_lo = h(T_min, P), h_hi = h_sat,L(P)
    /// VAPOR:    h_lo = h_sat,V(P), h_hi = h(T_max@p_max, P) (cap at pmax-isobar)
    /// SUPER:    h_lo = h(T_min, P), h_hi = h(T_max_ext, P)
    void populate_normph_bounds(NormalizedPHTable& table);

    /// Build a NormalizedPHTable end-to-end: set_limits + resize + populate
    /// bounds + sample HEOS at every (xnorm[i], yvec[j]) grid point, filling
    /// the LIST_OF_MATRICES values. Cells where the EOS update throws or
    /// returns a two-phase state are left as _HUGE holes (build_bspline_coeffs
    /// will mark them invalid).
    void build_normph_table(NormalizedPHTable& table);

    /// Build all three NormalizedPHTables (LIQUID, VAPOR, SUPER) and their
    /// per-cell SBTL coefficients.  Pure fluids only.  Idempotent.
    void build_normph_tables();

    /// Populate per-isobar T bounds (T_lo_isobar, T_hi_isobar) on a
    /// NormalizedPTTable.  LIQUID: T_lo = T_melt(P) (or T_min), T_hi = T_sat,L(P).
    /// VAPOR:  T_lo = T_sat,V(P), T_hi = T_max_ext.  SUPER: T_lo = T_melt (or T_min),
    /// T_hi = T_max_ext.  Also fits Cheb1D expansions of the non-saturation
    /// boundaries (LIQUID T_lo, VAPOR T_hi, SUPER both) for fast xnorm_from_T.
    void populate_normpt_bounds(NormalizedPTTable& table);

    /// Build a NormalizedPTTable end-to-end (mirror of build_normph_table).
    /// Cell-fills via PT_INPUTS (or PQ_INPUTS at sat boundaries).
    void build_normpt_table(NormalizedPTTable& table);

    /// Build all three NormalizedPTTables and their per-cell coefficients.
    /// Pure fluids only.  Idempotent.
    void build_normpt_tables();

    /// Helper: get T_sat,L(P) and T_sat,V(P) at the given subcritical pressure
    /// via PQ_INPUTS HEOS calls.  Throws if P is out of the saturation range.
    /// Single-deep cache so adjacent PT_INPUTS queries at the same P don't
    /// repeat the PQ inverse — common in pipeline / curve-trace workloads.
    void saturation_T_LV(double p, double& T_L, double& T_V);

    /// 2D natural cubic B-spline coefficient construction with per-cell
    /// polynomial expansion (Kunick-style C^2 splines).  Replaces independent
    /// per-cell 4-point Lagrange.  Includes a self-check that verifies the
    /// per-cell polynomial reproduces the corner data values to relative 1e-8.
    void build_bspline_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs);

    // -----------------------------------------------------------------------
    // Saturation curve queries (building blocks for coordinate-aligned PH)
    // -----------------------------------------------------------------------

    /// Saturation enthalpy h_sat,L(p) [J/mol].  Prefers a precomputed
    /// H-superancillary (machine-precision); on failure, calls
    /// AS->update(PQ_INPUTS, p, 0) and returns AS->hmolar() — slightly more
    /// expensive but still correct for any pure fluid with an EOS.
    double saturation_hmolar_liquid(double p) {
        if (try_h_superanc(p, /*Q=*/0, /*out=*/_hsat_tmp)) return _hsat_tmp;
        this->AS->update(PQ_INPUTS, p, 0.0);
        return static_cast<double>(this->AS->hmolar());
    }

    /// Saturation enthalpy h_sat,V(p) [J/mol]. Same precision story as
    /// saturation_hmolar_liquid().
    double saturation_hmolar_vapor(double p) {
        if (try_h_superanc(p, /*Q=*/1, /*out=*/_hsat_tmp)) return _hsat_tmp;
        this->AS->update(PQ_INPUTS, p, 1.0);
        return static_cast<double>(this->AS->hmolar());
    }

    /// Combined liquid+vapor saturation enthalpy lookup.  Computes T_sat(p)
    /// once via the inverse-pressure superancillary, then evaluates the H
    /// expansion at both Q=0 and Q=1.  Saves ~1 superancillary eval (1/3 of
    /// the cost of two separate calls).  Per-call cost: ~0.1 µs vs 0.15 µs
    /// for two separate saturation_hmolar_*(p) calls.  Plus a 1-deep cache
    /// for repeat queries at the same P (common in pipeline workloads).
    void saturation_hmolar_LV(double p, double& h_L, double& h_V);

    /// Property accessors overridden so the coordinate-aligned PH path can
    /// route property reads through the active NormalizedPHTable cell, and
    /// so the critical-region HEOS-fallback box can return its cached values.
    CoolPropDbl calc_hmolar(void);
    CoolPropDbl calc_smolar(void);
    CoolPropDbl calc_umolar(void);
    CoolPropDbl calc_rhomolar(void);

    // Critbox-aware derivative and transport accessors.  See the matching
    // .cpp definitions for the rationale.
    CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);
    CoolPropDbl calc_cpmolar();
    CoolPropDbl calc_cvmolar();
    CoolPropDbl calc_speed_sound();
    CoolPropDbl calc_viscosity();
    CoolPropDbl calc_conductivity();

    /// Prime the underlying HEOS at the state matching the active table:
    /// HmolarP_INPUTS for normph, PT_INPUTS for normpt.  Used by the
    /// derivative / transport accessor fallbacks.
    void ensure_heos_at_active_state();

   private:
    // Try the H expansion of the fluid's superancillary at pressure p.  Returns
    // true if h was evaluated and written to `out`; returns false if the H
    // expansion is unavailable (mixture, no SUPERANCILLARY block, or H
    // expansions not pre-built — most fluid JSONs ship rhoL/rhoV/p only).
    // The caller falls back to the tabular sat cache.
    //
    // Future work: derive H expansions on the fly via SuperAncillary::
    // add_variable('H', ...).  That requires the upstream SuperAncillary
    // template body to compile when ArrayType = std::vector<double>, which
    // currently uses Eigen-only operations (.matrix() on the coefficient
    // array).  Once upstream supports std::vector<double>, this fallback
    // path tightens from ~1e-3 to machine precision.
    bool try_h_superanc(double p, short Q, double& out) const;

    // Scratch slot for try_h_superanc to avoid an extra return-by-tuple in
    // the hot path of saturation_hmolar_{liquid,vapor}().
    mutable double _hsat_tmp{};

    // Single-deep cache for the combined L/V sat lookup.  When subsequent
    // queries land at the same P (common in pipeline / curve-trace
    // workloads), we reuse the result without recomputing T_sat.
    mutable double _hsat_LV_cache_p{-1.0};
    mutable double _hsat_LV_cache_hL{0.0};
    mutable double _hsat_LV_cache_hV{0.0};

    // Same shape, T_sat,L(P) / T_sat,V(P) for the PT routing path.
    mutable double _Tsat_LV_cache_p{-1.0};
    mutable double _Tsat_LV_cache_TL{0.0};
    mutable double _Tsat_LV_cache_TV{0.0};

    // Coordinate-aligned PH tables — saturation curve becomes a coordinate axis.
    // See dev/sbtl_normalized_ph_design.md for the math.
    NormalizedPHTable _normph_liquid;
    NormalizedPHTable _normph_vapor;
    NormalizedPHTable _normph_super;
    std::vector<std::vector<CellCoeffs>> _coeffs_normph_liquid;
    std::vector<std::vector<CellCoeffs>> _coeffs_normph_vapor;
    std::vector<std::vector<CellCoeffs>> _coeffs_normph_super;
    bool normph_tables_built;

    // Set during update(HmolarP/HmassP_INPUTS) when the query routes through
    // a normalized table.  evaluate_single_phase_phmolar consults these and
    // uses _normph_xnorm in place of _hmolar so cell-eval stays consistent
    // with the cached i/j indices on the active normph table.
    NormalizedPHTable* _normph_active_table = nullptr;
    const std::vector<std::vector<CellCoeffs>>* _normph_active_coeffs = nullptr;
    double _normph_xnorm = 0.0;

    // Coordinate-aligned PT tables — same coordinate trick as the PH path.
    // See dev/sbtl_pt_outstanding_work.md for the design rationale (fix for
    // the 5 % compressed-liquid PT residue and the triple-corner hang).
    NormalizedPTTable _normpt_liquid;
    NormalizedPTTable _normpt_vapor;
    NormalizedPTTable _normpt_super;
    std::vector<std::vector<CellCoeffs>> _coeffs_normpt_liquid;
    std::vector<std::vector<CellCoeffs>> _coeffs_normpt_vapor;
    std::vector<std::vector<CellCoeffs>> _coeffs_normpt_super;
    bool normpt_tables_built;

    // Same active-table state for PT_INPUTS routing.  When a PT query
    // succeeds through one of the normpt tables, _normpt_active_* point
    // to that table + coeffs, _normpt_xnorm caches the xnorm coordinate
    // so the property accessors evaluate the polynomial at (xnorm, P)
    // (not (T, P)) — matching how cell-build encoded the data.
    NormalizedPTTable* _normpt_active_table = nullptr;
    const std::vector<std::vector<CellCoeffs>>* _normpt_active_coeffs = nullptr;
    double _normpt_xnorm = 0.0;

    // Cached molar enthalpy at the critical point — used to define the
    // tight 2D HEOS-fallback box around (h_crit, p_crit).  Computed
    // lazily on the first HmolarP query that lands near pcrit.
    mutable std::optional<double> _hmolar_crit_cached;

    // Critical-region HEOS-fallback box state cache.  When _critbox_active
    // is true, every calc_* property accessor returns the value stashed
    // here (computed by a direct AS->update(...) call) instead of routing
    // through the table machinery — the polynomial backbone has no valid
    // cell to evaluate near the critical cusp.  The full set of properties
    // is cached so derivative quantities (cp/cv/speed_sound) and transport
    // (viscosity/conductivity) all get the HEOS-direct answer too, rather
    // than falling through to a wrong-table polynomial evaluation.
    bool _critbox_active{false};
    double _critbox_smolar{0.0};
    double _critbox_umolar{0.0};
    double _critbox_cpmolar{0.0};
    double _critbox_cvmolar{0.0};
    double _critbox_speed_sound{0.0};
    double _critbox_viscosity{0.0};
    double _critbox_conductivity{0.0};
};

}  // namespace CoolProp

#endif  // SBTLBACKEND_H
