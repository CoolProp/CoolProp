#ifndef PHASE_ENVELOPE_TRACERS_H
#define PHASE_ENVELOPE_TRACERS_H

/**
 * Experimental isopleth tracers for mixture phase envelopes, selected by the
 * configuration key PHASE_ENVELOPE_ALGORITHM.  The legacy density-marching
 * tracer in PhaseEnvelopeRoutines::build stays the default; the tracers here
 * are candidates that are compared against it on the predefined-mixture
 * torture corpus before any promotion.
 *
 * All candidates share one continuation driver (Michelsen 1980 style:
 * specification equation, tangent predictor, Newton corrector, step size from
 * the corrector iteration count, critical point crossed by the simultaneous
 * sign change of all ln K).  They differ only in the unknown vector and the
 * residual/Jacobian provider, the IsoplethSystem.
 *
 * Design: docs/superpowers/specs/2026-09-11-phase-envelope-tracers-design.md
 */

#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "CoolProp/fluids/PhaseEnvelope.h"

#include <memory>
#include <string>
#include <vector>

namespace CoolProp {

class PhaseEnvelopeTracers
{
   public:
    /// One converged point on the boundary, in the storage convention of
    /// PhaseEnvelopeData: the incipient phase is stored as "liquid", the feed
    /// phase as "vapor", all the way around the envelope.
    struct TracedPoint
    {
        double T = 0, p = 0, rho_inc = 0, rho_feed = 0, h_inc = 0, h_feed = 0, s_inc = 0, s_feed = 0;
        double max_abs_lnK = 0;
        std::vector<CoolPropDbl> x_inc;
    };

    struct Options
    {
        double dS_initial = 0.05;           ///< first step, log units of the marching variable
        double dS_max_log = 0.2;            ///< step cap when the specified variable is ln T, ln p or ln rho
        double dS_max_lnK = 0.05;           ///< step cap when the specified variable is a ln K
        double dS_max_lnK_critical = 0.01;  ///< step cap for a ln K spec inside the near-critical band
        double lnK_critical = 0.01;         ///< near-critical band: max_i |ln K_i| below this
        double lnK_cap_reference = 5.0;     ///< the ln K step cap is scaled up by max(1, max|ln K| / this)
        double max_abs_lnK = 50.0;          ///< stop once max_i |ln K_i| exceeds this (or 1.2x its value at the start, whichever is larger)
        double dS_min = 1e-6;               ///< give up after halving below this
        double max_step_lnT = 0.02;         ///< point-density cap on the predicted change of ln T
        double max_step_lnmarch = 0.1;      ///< point-density cap on the predicted change of the marching variable
        int corrector_max_iter = 10;
        double corrector_tol = 1e-9;        ///< on max |F|
        double corrector_step_tol = 1e-11;  ///< a Newton step below this (log units) counts as converged if max |F| < corrector_noise_tol
        double corrector_noise_tol = 1e-6;  ///< residual noise floor accepted together with a negligible step
        std::size_t max_points = 1000;
        double p_ceiling = 1e9;  ///< [Pa] stop for open envelopes
        /// [K] stop when T drops below this after the critical crossing.  Disabled by
        /// default: HEOS.Tmin() is a mole-fraction-weighted number that real envelopes run
        /// well below (methane/ethane traces to 65 K against a weighted Tmin of 91 K), so
        /// using it as a floor truncates envelopes that would otherwise close.  A trace that
        /// runs out of EOS stops as "stalled" instead, which is reported, not silent.
        double T_floor = 0;
        double merge_rho_tol = 0.02;       ///< densities within this relative distance count as merged
        double merge_lnK_tol = 0.1;        ///< ...which is only legitimate when max|ln K| is below this (a real critical point)
        double closure_rho_ratio = 100.0;  ///< closure also needs the two densities to differ by at least this factor
        std::size_t min_points_for_built = 5;
        int start_retries = 4;  ///< decades of start pressure tried after the configured one
    };

    /** Residual and Jacobian provider for one choice of unknowns.
     *
     * Composition convention shared by all systems: z is the feed composition,
     * K_i is feed over incipient mole fraction (the PhaseEnvelopeData K = y/x
     * convention), w_i = z_i / K_i are the unnormalised incipient amounts,
     * W = sum w, and x_i = w_i / W is the incipient composition used to
     * evaluate the EOS.  The scale redundancy of K is removed by the equation
     * W - 1 = 0.
     */
    class IsoplethSystem
    {
       public:
        explicit IsoplethSystem(HelmholtzEOSMixtureBackend& HEOS);
        virtual ~IsoplethSystem() = default;
        [[nodiscard]] virtual std::string name() const = 0;
        [[nodiscard]] virtual std::size_t size() const = 0;
        [[nodiscard]] virtual std::size_t index_lnT() const = 0;
        [[nodiscard]] virtual std::size_t index_marching() const = 0;
        [[nodiscard]] std::size_t ncomp() const {
            return N;
        }
        [[nodiscard]] bool is_lnK(std::size_t i) const {
            return i < N;
        }
        /// Evaluate residuals and Jacobian at X for the specification X[ns] = S.
        /// Leaves the SatL (incipient) and SatV (feed) instances at the evaluated states.
        virtual void residual_jacobian(const Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F, Eigen::MatrixXd& J) = 0;
        /// Physical state at X; valid after residual_jacobian was evaluated at this X.
        virtual void unpack(const Eigen::VectorXd& X, TracedPoint& pt) const;
        /// Unknown vector from a converged starting saturation state (x = incipient, y = feed).
        [[nodiscard]] virtual Eigen::VectorXd pack(const SaturationSolvers::newton_raphson_saturation_options& s0) const = 0;
        /// Pressure scale for the pressure-equality residual (density system only); harmless elsewhere.
        virtual void set_p_ref(double) {}

       protected:
        HelmholtzEOSMixtureBackend* HEOS;  ///< not owned
        std::size_t N;
        std::vector<CoolPropDbl> z;
        void incipient_composition(const Eigen::VectorXd& X, std::vector<CoolPropDbl>& w, std::vector<CoolPropDbl>& x, double& W) const;
        /// Chain rule from XN_DEPENDENT composition derivatives D_k (k < N-1) to d/dlnK_j
        /// for a quantity evaluated at the normalised incipient composition.
        static double dlnK_chain(const std::vector<CoolPropDbl>& x, const std::vector<CoolPropDbl>& D, std::size_t j);
    };

    /// Unknowns [ln K_1..ln K_N, ln T, ln rho_inc, ln rho_feed]; DmolarT updates, no density root solve.
    class LnKDensitySystem : public IsoplethSystem
    {
       public:
        explicit LnKDensitySystem(HelmholtzEOSMixtureBackend& HEOS) : IsoplethSystem(HEOS) {}
        [[nodiscard]] std::string name() const override {
            return "lnK_density";
        }
        [[nodiscard]] std::size_t size() const override {
            return N + 3;
        }
        [[nodiscard]] std::size_t index_lnT() const override {
            return N;
        }
        [[nodiscard]] std::size_t index_marching() const override {
            return N + 2;
        }
        void residual_jacobian(const Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F, Eigen::MatrixXd& J) override;
        [[nodiscard]] Eigen::VectorXd pack(const SaturationSolvers::newton_raphson_saturation_options& s0) const override;
        void set_p_ref(double p) override {
            p_ref = p;
        }

       private:
        double p_ref = 1e5;
    };

    /// Unknowns [ln K_1..ln K_N, ln T, ln p]; classic Michelsen form with a density root solve per phase.
    class LnKPressureSystem : public IsoplethSystem
    {
       public:
        explicit LnKPressureSystem(HelmholtzEOSMixtureBackend& HEOS) : IsoplethSystem(HEOS) {}
        [[nodiscard]] std::string name() const override {
            return "lnK_pressure";
        }
        [[nodiscard]] std::size_t size() const override {
            return N + 2;
        }
        [[nodiscard]] std::size_t index_lnT() const override {
            return N;
        }
        [[nodiscard]] std::size_t index_marching() const override {
            return N + 1;
        }
        void residual_jacobian(const Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F, Eigen::MatrixXd& J) override;
        [[nodiscard]] Eigen::VectorXd pack(const SaturationSolvers::newton_raphson_saturation_options& s0) const override;

       private:
        double rho_inc_guess = -1, rho_feed_guess = -1;
    };

    /// Names accepted by PHASE_ENVELOPE_ALGORITHM, "legacy" first.
    static std::vector<std::string> algorithms();

    /// Construct the system for a name other than "legacy"; throws ValueError for unknown names.
    static std::unique_ptr<IsoplethSystem> make_system(HelmholtzEOSMixtureBackend& HEOS, const std::string& algorithm);

    /// Options implied by the level string of build_phase_envelope ("", "veryfine", "none", anything else).
    static Options options_for_level(const std::string& level);

    /** Legacy starting point (low-pressure dew point via preconditioner, Wilson,
     * successive substitution and one pressure-imposed Newton solve), retried at
     * ten times the pressure up to opts.start_retries times.  Returns the converged
     * state with x = incipient liquid, y = feed; s0.p is the pressure actually used.
     */
    static SaturationSolvers::newton_raphson_saturation_options starting_point(HelmholtzEOSMixtureBackend& HEOS, const Options& opts);

    /// One attempt at the legacy starting point at exactly p_start; throws on failure.
    static SaturationSolvers::newton_raphson_saturation_options starting_point_at(HelmholtzEOSMixtureBackend& HEOS, double p_start);

    /// Run the continuation with a given system and options, filling HEOS.PhaseEnvelope.
    static void run(HelmholtzEOSMixtureBackend& HEOS, IsoplethSystem& sys, const Options& opts);

    /// Why the most recent trace on this thread stopped ("closed", "floor", "ceiling", "pure", "max_points", "stalled").
    static const std::string& last_stop_reason();
    static const std::string& last_stop_detail();

    /// Entry point used by PhaseEnvelopeRoutines::build for any algorithm other than "legacy".
    static void trace(HelmholtzEOSMixtureBackend& HEOS, const std::string& algorithm, const std::string& level);
};

} /* namespace CoolProp */

#endif
