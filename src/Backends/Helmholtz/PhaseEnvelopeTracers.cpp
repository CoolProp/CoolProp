#include "PhaseEnvelopeTracers.h"
#include "MixtureDerivatives.h"
#include "CoolProp/Configuration.h"
#include "CoolProp/detail/tools.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace CoolProp {

// ---------------------------------------------------------------------------
// IsoplethSystem
// ---------------------------------------------------------------------------

PhaseEnvelopeTracers::IsoplethSystem::IsoplethSystem(HelmholtzEOSMixtureBackend& HEOS)
  : HEOS(&HEOS), N(HEOS.get_mole_fractions_ref().size()), z(HEOS.get_mole_fractions_ref()) {
    if (N < 2) {
        throw ValueError("PhaseEnvelopeTracers: a mixture with at least two components is required");
    }
    if (!HEOS.SatL || !HEOS.SatV) {
        throw ValueError("PhaseEnvelopeTracers: SatL/SatV instances are not available");
    }
}

void PhaseEnvelopeTracers::IsoplethSystem::incipient_composition(const Eigen::VectorXd& X, std::vector<CoolPropDbl>& w, std::vector<CoolPropDbl>& x,
                                                                 double& W) const {
    w.resize(N);
    x.resize(N);
    W = 0;
    for (std::size_t i = 0; i < N; ++i) {
        w[i] = z[i] * std::exp(-X[static_cast<Eigen::Index>(i)]);  // z_i / K_i
        W += w[i];
    }
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = w[i] / W;
    }
}

double PhaseEnvelopeTracers::IsoplethSystem::dlnK_chain(const std::vector<CoolPropDbl>& x, const std::vector<CoolPropDbl>& D, std::size_t j) {
    // dx_k/dlnK_j = x_j (x_k - delta_kj); with XN_DEPENDENT derivatives D_k (k < N-1)
    // the contraction over all k reduces to sum_{k<N-1} D_k dx_k/dlnK_j because the
    // normalised composition has zero-sum differentials.
    const std::size_t Nm1 = x.size() - 1;
    double s = 0;
    for (std::size_t k = 0; k < Nm1; ++k) {
        s += D[k] * x[k];
    }
    if (j < Nm1) {
        s -= D[j];
    }
    return x[j] * s;
}

void PhaseEnvelopeTracers::IsoplethSystem::unpack(const Eigen::VectorXd& X, TracedPoint& pt) const {
    HelmholtzEOSMixtureBackend &inc = *HEOS->SatL, &feed = *HEOS->SatV;
    pt.T = inc.T();
    pt.p = 0.5 * (inc.p() + feed.p());
    pt.rho_inc = inc.rhomolar();
    pt.rho_feed = feed.rhomolar();
    pt.h_inc = inc.hmolar();
    pt.h_feed = feed.hmolar();
    pt.s_inc = inc.smolar();
    pt.s_feed = feed.smolar();
    pt.x_inc = inc.get_mole_fractions();
    pt.max_abs_lnK = 0;
    for (std::size_t i = 0; i < N; ++i) {
        pt.max_abs_lnK = std::max(pt.max_abs_lnK, std::abs(X[static_cast<Eigen::Index>(i)]));
    }
}

// ---------------------------------------------------------------------------
// LnKDensitySystem
// ---------------------------------------------------------------------------

Eigen::VectorXd PhaseEnvelopeTracers::LnKDensitySystem::pack(const SaturationSolvers::newton_raphson_saturation_options& s0) const {
    Eigen::VectorXd X(static_cast<Eigen::Index>(size()));
    for (std::size_t i = 0; i < N; ++i) {
        X[static_cast<Eigen::Index>(i)] = std::log(z[i] / s0.x[i]);
    }
    X[static_cast<Eigen::Index>(N)] = std::log(s0.T);
    X[static_cast<Eigen::Index>(N + 1)] = std::log(s0.rhomolar_liq);
    X[static_cast<Eigen::Index>(N + 2)] = std::log(s0.rhomolar_vap);
    return X;
}

void PhaseEnvelopeTracers::LnKDensitySystem::residual_jacobian(const Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F,
                                                               Eigen::MatrixXd& J) {
    const auto n = static_cast<Eigen::Index>(size());
    const auto iN = static_cast<Eigen::Index>(N);
    F.setZero(n);
    J.setZero(n, n);

    std::vector<CoolPropDbl> w, x;
    double W = 0;
    incipient_composition(X, w, x, W);
    const double T = std::exp(X[iN]), rho_inc = std::exp(X[iN + 1]), rho_feed = std::exp(X[iN + 2]);

    HelmholtzEOSMixtureBackend &inc = *HEOS->SatL, &feed = *HEOS->SatV;
    inc.set_mole_fractions(x);
    inc.update(DmolarT_INPUTS, rho_inc, T);
    feed.set_mole_fractions(z);
    feed.update(DmolarT_INPUTS, rho_feed, T);

    const x_N_dependency_flag xN = XN_DEPENDENT;
    std::vector<CoolPropDbl> D(N - 1);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ii = static_cast<Eigen::Index>(i);
        F[ii] = std::log(MixtureDerivatives::fugacity_i(inc, i, xN)) - std::log(MixtureDerivatives::fugacity_i(feed, i, xN));
        for (std::size_t k = 0; k + 1 < N; ++k) {
            D[k] = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(inc, i, k, xN);
        }
        for (std::size_t j = 0; j < N; ++j) {
            J(ii, static_cast<Eigen::Index>(j)) = dlnK_chain(x, D, j);
        }
        J(ii, iN) =
          T * (MixtureDerivatives::dln_fugacity_i_dT__constrho_n(inc, i, xN) - MixtureDerivatives::dln_fugacity_i_dT__constrho_n(feed, i, xN));
        J(ii, iN + 1) = rho_inc * MixtureDerivatives::dln_fugacity_i_drho__constT_n(inc, i, xN);
        J(ii, iN + 2) = -rho_feed * MixtureDerivatives::dln_fugacity_i_drho__constT_n(feed, i, xN);
    }
    // Normalisation: W - 1 = 0
    F[iN] = W - 1;
    for (std::size_t j = 0; j < N; ++j) {
        J(iN, static_cast<Eigen::Index>(j)) = -w[j];
    }
    // Pressure equality.  p_ref (held fixed during a corrector) is floored at 1e-3 rho_inc R T:
    // the liquid pressure is a difference of terms of that order, so its round-off noise is a
    // fixed fraction of it, and a 100 Pa dew point must not demand pressure equality to 1e-9
    // of 100 Pa.
    const double p_scale = std::max(p_ref, 1e-3 * rho_inc * inc.gas_constant() * T);
    F[iN + 1] = (inc.p() - feed.p()) / p_scale;
    std::vector<CoolPropDbl> P(N - 1);
    for (std::size_t k = 0; k + 1 < N; ++k) {
        P[k] = MixtureDerivatives::dpdxj__constT_V_xi(inc, k, xN);
    }
    for (std::size_t j = 0; j < N; ++j) {
        J(iN + 1, static_cast<Eigen::Index>(j)) = dlnK_chain(x, P, j) / p_scale;
    }
    J(iN + 1, iN) = T * (MixtureDerivatives::dpdT__constV_n(inc) - MixtureDerivatives::dpdT__constV_n(feed)) / p_scale;
    J(iN + 1, iN + 1) = rho_inc * MixtureDerivatives::dpdrho__constT_n(inc) / p_scale;
    J(iN + 1, iN + 2) = -rho_feed * MixtureDerivatives::dpdrho__constT_n(feed) / p_scale;
    // Specification
    F[iN + 2] = X[static_cast<Eigen::Index>(ns)] - S;
    J(iN + 2, static_cast<Eigen::Index>(ns)) = 1;
}

// ---------------------------------------------------------------------------
// LnKPressureSystem
// ---------------------------------------------------------------------------

Eigen::VectorXd PhaseEnvelopeTracers::LnKPressureSystem::pack(const SaturationSolvers::newton_raphson_saturation_options& s0) const {
    Eigen::VectorXd X(static_cast<Eigen::Index>(size()));
    for (std::size_t i = 0; i < N; ++i) {
        X[static_cast<Eigen::Index>(i)] = std::log(z[i] / s0.x[i]);
    }
    X[static_cast<Eigen::Index>(N)] = std::log(s0.T);
    X[static_cast<Eigen::Index>(N + 1)] = std::log(s0.p);
    // The density guesses are mutable trace state; pack is const, so seed them via a cast.
    auto* self = const_cast<LnKPressureSystem*>(this);  // NOLINT(cppcoreguidelines-pro-type-const-cast)
    self->rho_inc_guess = s0.rhomolar_liq;
    self->rho_feed_guess = s0.rhomolar_vap;
    return X;
}

void PhaseEnvelopeTracers::LnKPressureSystem::residual_jacobian(const Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F,
                                                                Eigen::MatrixXd& J) {
    const auto n = static_cast<Eigen::Index>(size());
    const auto iN = static_cast<Eigen::Index>(N);
    F.setZero(n);
    J.setZero(n, n);

    std::vector<CoolPropDbl> w, x;
    double W = 0;
    incipient_composition(X, w, x, W);
    const double T = std::exp(X[iN]), p = std::exp(X[iN + 1]);

    HelmholtzEOSMixtureBackend &inc = *HEOS->SatL, &feed = *HEOS->SatV;
    inc.set_mole_fractions(x);
    inc.update_TP_guessrho(T, p, rho_inc_guess);
    rho_inc_guess = inc.rhomolar();
    feed.set_mole_fractions(z);
    feed.update_TP_guessrho(T, p, rho_feed_guess);
    rho_feed_guess = feed.rhomolar();

    const x_N_dependency_flag xN = XN_DEPENDENT;
    std::vector<CoolPropDbl> D(N - 1);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ii = static_cast<Eigen::Index>(i);
        F[ii] = std::log(MixtureDerivatives::fugacity_i(inc, i, xN)) - std::log(MixtureDerivatives::fugacity_i(feed, i, xN));
        for (std::size_t k = 0; k + 1 < N; ++k) {
            D[k] = MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(inc, i, k, xN);
        }
        for (std::size_t j = 0; j < N; ++j) {
            J(ii, static_cast<Eigen::Index>(j)) = dlnK_chain(x, D, j);
        }
        J(ii, iN) = T * (MixtureDerivatives::dln_fugacity_i_dT__constp_n(inc, i, xN) - MixtureDerivatives::dln_fugacity_i_dT__constp_n(feed, i, xN));
        J(ii, iN + 1) =
          p * (MixtureDerivatives::dln_fugacity_i_dp__constT_n(inc, i, xN) - MixtureDerivatives::dln_fugacity_i_dp__constT_n(feed, i, xN));
    }
    F[iN] = W - 1;
    for (std::size_t j = 0; j < N; ++j) {
        J(iN, static_cast<Eigen::Index>(j)) = -w[j];
    }
    F[iN + 1] = X[static_cast<Eigen::Index>(ns)] - S;
    J(iN + 1, static_cast<Eigen::Index>(ns)) = 1;
}

// ---------------------------------------------------------------------------
// Factory, options, starting point
// ---------------------------------------------------------------------------

std::vector<std::string> PhaseEnvelopeTracers::algorithms() {
    return {"legacy", "lnK_density", "lnK_pressure"};
}

std::unique_ptr<PhaseEnvelopeTracers::IsoplethSystem> PhaseEnvelopeTracers::make_system(HelmholtzEOSMixtureBackend& HEOS,
                                                                                        const std::string& algorithm) {
    if (algorithm == "lnK_density") {
        return std::unique_ptr<IsoplethSystem>(new LnKDensitySystem(HEOS));
    }
    if (algorithm == "lnK_pressure") {
        return std::unique_ptr<IsoplethSystem>(new LnKPressureSystem(HEOS));
    }
    throw ValueError(format("PHASE_ENVELOPE_ALGORITHM [%s] is not valid; options are %s", algorithm.c_str(), strjoin(algorithms(), ", ").c_str()));
}

PhaseEnvelopeTracers::Options PhaseEnvelopeTracers::options_for_level(const std::string& level) {
    Options o;
    if (level == "veryfine") {
        o.max_step_lnT = 0.01;
        o.max_step_lnmarch = 0.05;
    } else if (level == "none") {
        o.max_step_lnT = 0.05;
        o.max_step_lnmarch = 0.2;
    }
    return o;
}

SaturationSolvers::newton_raphson_saturation_options PhaseEnvelopeTracers::starting_point_at(HelmholtzEOSMixtureBackend& HEOS, double p_start) {
    SaturationSolvers::mixture_VLE_IO io;
    io.sstype = SaturationSolvers::imposed_p;
    io.Nstep_max = 20;
    io.beta = 1;
    const std::vector<CoolPropDbl>& z = HEOS.get_mole_fractions_ref();
    double Tguess = SaturationSolvers::saturation_preconditioner(HEOS, p_start, SaturationSolvers::imposed_p, z);
    Tguess = SaturationSolvers::saturation_Wilson(HEOS, 1.0, p_start, SaturationSolvers::imposed_p, z, Tguess);
    SaturationSolvers::successive_substitution(HEOS, 1.0, Tguess, p_start, z, HEOS.K, io);

    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    IO.bubble_point = false;  // dew point: liquid (x) is the incipient phase
    IO.x = io.x;
    IO.y = z;
    IO.rhomolar_liq = io.rhomolar_liq;
    IO.rhomolar_vap = io.rhomolar_vap;
    IO.T = io.T;
    IO.p = io.p;
    IO.Nstep_max = 30;
    IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;
    NR.call(HEOS, IO.y, IO.x, IO);

    bool ok = ValidNumber(IO.T) && ValidNumber(IO.p) && ValidNumber(IO.rhomolar_liq) && ValidNumber(IO.rhomolar_vap) && IO.p > 0
              && IO.rhomolar_liq > IO.rhomolar_vap;
    for (std::size_t i = 0; ok && i < IO.x.size(); ++i) {
        ok = ValidNumber(IO.x[i]) && IO.x[i] > 0 && IO.x[i] <= 1;
    }
    if (!ok) {
        throw ValueError("starting dew point is not a valid two-phase state");
    }
    return IO;
}

SaturationSolvers::newton_raphson_saturation_options PhaseEnvelopeTracers::starting_point(HelmholtzEOSMixtureBackend& HEOS, const Options& opts) {
    const bool debug = get_debug_level() > 0;
    double p_start = get_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA);
    std::string last_error;
    for (int attempt = 0; attempt <= opts.start_retries; ++attempt, p_start *= 10) {
        try {
            SaturationSolvers::newton_raphson_saturation_options IO = starting_point_at(HEOS, p_start);
            if (debug) {
                std::cout << format("PhaseEnvelopeTracers: start at p = %g Pa, T = %g K (attempt %d)\n", IO.p, IO.T, attempt);
            }
            return IO;
        } catch (std::exception& e) {
            last_error = e.what();
            if (debug) {
                std::cout << format("PhaseEnvelopeTracers: start at p = %g Pa failed: %s\n", p_start, e.what());
            }
        }
    }
    throw ValueError(format("PhaseEnvelopeTracers: unable to obtain a starting dew point after %d attempts; last error: %s", opts.start_retries + 1,
                            last_error.c_str()));
}

// ---------------------------------------------------------------------------
// Continuation driver
// ---------------------------------------------------------------------------

namespace {

thread_local std::string g_stop_reason, g_stop_detail;

Eigen::Index argmax_abs(const Eigen::VectorXd& v) {
    Eigen::Index i = 0;
    v.cwiseAbs().maxCoeff(&i);
    return i;
}

struct CorrectorResult
{
    bool ok = false;
    int iters = 0;
    std::string message;
};

/// Newton corrector on F(X; ns, S) = 0 starting from X; on success X, F, J hold the converged state.
CorrectorResult corrector(PhaseEnvelopeTracers::IsoplethSystem& sys, Eigen::VectorXd& X, std::size_t ns, double S, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, const PhaseEnvelopeTracers::Options& opts) {
    CorrectorResult res;
    for (int it = 0; it < opts.corrector_max_iter; ++it) {
        try {
            sys.residual_jacobian(X, ns, S, F, J);
        } catch (std::exception& e) {
            res.message = std::string("EOS evaluation failed: ") + e.what();
            return res;
        }
        if (!F.allFinite() || !J.allFinite()) {
            res.message = "non-finite residual or Jacobian";
            return res;
        }
        res.iters = it + 1;
        if (get_debug_level() > 1) {
            std::cout << format("    corrector it %d: max|F| = %g (row %d)\n", it, F.cwiseAbs().maxCoeff(), static_cast<int>(argmax_abs(F)));
        }
        if (F.cwiseAbs().maxCoeff() < opts.corrector_tol) {
            res.ok = true;
            return res;
        }
        Eigen::VectorXd dX = J.colPivHouseholderQr().solve(-F);
        if (!dX.allFinite()) {
            res.message = "singular Jacobian";
            return res;
        }
        // Guard against wild steps: all unknowns are logarithms, so a unit step is already huge
        const double m = dX.cwiseAbs().maxCoeff();
        if (m > 1.0) {
            dX *= 1.0 / m;
        }
        X += dX;
        if (m < opts.corrector_step_tol) {
            // The Newton step is below the resolution of the unknowns: converged to the noise
            // floor of the residuals.  Re-evaluate at the final X so F and J are consistent.
            try {
                sys.residual_jacobian(X, ns, S, F, J);
            } catch (std::exception& e) {
                res.message = std::string("EOS evaluation failed: ") + e.what();
                return res;
            }
            res.ok = F.allFinite() && J.allFinite() && F.cwiseAbs().maxCoeff() < opts.corrector_noise_tol;
            if (!res.ok) {
                res.message = format("stalled with max|F| = %g", F.cwiseAbs().maxCoeff());
            }
            return res;
        }
    }
    res.message = format("corrector did not converge in %d iterations", opts.corrector_max_iter);
    return res;
}

}  // namespace

const std::string& PhaseEnvelopeTracers::last_stop_reason() {
    return g_stop_reason;
}
const std::string& PhaseEnvelopeTracers::last_stop_detail() {
    return g_stop_detail;
}
void PhaseEnvelopeTracers::set_last_stop(const std::string& reason, const std::string& detail) {
    g_stop_reason = reason;
    g_stop_detail = detail;
}

void PhaseEnvelopeTracers::run(HelmholtzEOSMixtureBackend& HEOS, IsoplethSystem& sys, const Options& opts) {
    const bool debug = get_debug_level() > 0;
    g_stop_reason = "none";
    g_stop_detail.clear();
    const double T_floor = opts.T_floor;
    PhaseEnvelopeData& env = HEOS.PhaseEnvelope;
    const std::size_t N = sys.ncomp();
    const auto n = static_cast<Eigen::Index>(sys.size());
    const std::size_t i_lnT = sys.index_lnT(), i_march = sys.index_marching();

    env.resize(N);
    bool crossed_any = false;  // env.icrit is unsigned; -1 is its 'not set' sentinel

    // Start: the legacy low-pressure dew point, polished in the tracer's own variables.  A
    // start that cannot be polished (degenerate incipient phase, e.g. a heavy trace
    // component below its own triple point) is retried a decade higher in pressure.
    double p_start = get_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA);
    Eigen::VectorXd X, F(n);
    Eigen::MatrixXd J(n, n);
    std::size_t ns = i_march;
    CorrectorResult cr;
    bool started = false;
    std::string last_error;
    for (int attempt = 0; attempt <= opts.start_retries && !started; ++attempt, p_start *= 10) {
        try {
            SaturationSolvers::newton_raphson_saturation_options s0 = starting_point_at(HEOS, p_start);
            X = sys.pack(s0);
            sys.set_p_ref(s0.p);
            cr = corrector(sys, X, ns, X[static_cast<Eigen::Index>(ns)], F, J, opts);
            if (!cr.ok) {
                throw ValueError(format("could not polish the starting point in the %s variables: %s", sys.name().c_str(), cr.message.c_str()));
            }
            p_start = s0.p;
            started = true;
        } catch (std::exception& e) {
            last_error = e.what();
            if (debug) {
                std::cout << format("PhaseEnvelopeTracers: start at p = %g Pa failed: %s\n", p_start, e.what());
            }
        }
    }
    if (!started) {
        throw ValueError(
          format("PhaseEnvelopeTracers: unable to start after %d attempts; last error: %s", opts.start_retries + 1, last_error.c_str()));
    }
    if (debug) {
        std::cout << format("PhaseEnvelopeTracers: started at p = %g Pa\n", p_start);
    }
    auto store = [&](TracedPoint& pt) {
        env.store_variables(pt.T, pt.p, pt.rho_inc, pt.rho_feed, pt.h_inc, pt.h_feed, pt.s_inc, pt.s_feed, pt.x_inc, HEOS.get_mole_fractions_ref());
    };
    TracedPoint pt;
    sys.unpack(X, pt);
    store(pt);
    sys.set_p_ref(pt.p);
    // The degeneracy limit is relative to the spread already present at the start: a wide-boiling
    // multicomponent gas can begin at max|lnK| of 60 or more (n-octane against helium at 100 Pa),
    // and an absolute limit would reject its first point.
    const double lnK_limit = std::max(opts.max_abs_lnK, 1.2 * pt.max_abs_lnK);

    Eigen::VectorXd prev_tangent = Eigen::VectorXd::Zero(n);
    double dS = opts.dS_initial;
    int iters_prev = 3;
    std::string stop_reason;

    for (;;) {
        // Tangent from the converged Jacobian (spec row is a unit row, so dXdS[ns] = 1)
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
        rhs[n - 1] = 1;
        Eigen::VectorXd dXdS = J.colPivHouseholderQr().solve(rhs);
        if (!dXdS.allFinite()) {
            g_stop_reason = "stalled";
            stop_reason = "singular Jacobian for the tangent";
            break;
        }
        // Direction continuity: first step increases the marching variable, later steps follow the previous tangent
        const double orient = prev_tangent.isZero() ? dXdS[static_cast<Eigen::Index>(i_march)] : dXdS.dot(prev_tangent);
        const double sgn = orient >= 0 ? 1.0 : -1.0;
        dS = sgn * std::abs(dS);

        // Specification switch: largest tangent component; only ln K inside the near-critical band
        const double maxlnK = X.head(static_cast<Eigen::Index>(N)).cwiseAbs().maxCoeff();
        std::size_t ns_new = ns;
        double best = -1;
        for (std::size_t i = 0; i < sys.size(); ++i) {
            if (maxlnK < opts.lnK_critical && !sys.is_lnK(i)) {
                continue;
            }
            const double a = std::abs(dXdS[static_cast<Eigen::Index>(i)]);
            if (a > best) {
                best = a;
                ns_new = i;
            }
        }
        if (ns_new != ns) {
            const double scale = dXdS[static_cast<Eigen::Index>(ns_new)];
            dXdS /= scale;
            dS *= scale;
            ns = ns_new;
        }

        // Step size: iteration-count adaptation, caps, point-density limits
        // Grow on an easy corrector, hold on a moderate one, shrink on a hard one
        if (iters_prev <= 3) {
            dS *= 1.5;
        } else if (iters_prev >= 6) {
            dS *= 0.5;
        }
        double cap = opts.dS_max_log;
        if (sys.is_lnK(ns)) {
            // Scale the cap with the magnitude of ln K.  Deep in a low-temperature tail the
            // incipient phase goes numerically pure and ln K races to tens while T and the
            // marching density barely move, so a fixed cap spends hundreds of points there
            // for no resolution gain.  Near the critical point the scale factor is 1 and the
            // tight critical cap still applies.
            cap = (maxlnK < opts.lnK_critical) ? opts.dS_max_lnK_critical : opts.dS_max_lnK * std::max(1.0, maxlnK / opts.lnK_cap_reference);
        }
        dS = (dS >= 0 ? 1.0 : -1.0) * std::min(std::abs(dS), cap);
        {
            const double dT = std::abs(dXdS[static_cast<Eigen::Index>(i_lnT)] * dS), dM = std::abs(dXdS[static_cast<Eigen::Index>(i_march)] * dS);
            double f = 1;
            if (dT > opts.max_step_lnT) f = std::min(f, opts.max_step_lnT / dT);
            if (dM > opts.max_step_lnmarch) f = std::min(f, opts.max_step_lnmarch / dM);
            dS *= f;
        }

        // Predictor / corrector with step halving on rejection
        bool accepted = false;
        Eigen::VectorXd X_new, F_new;
        Eigen::MatrixXd J_new;
        TracedPoint pt_new;
        std::string reject;
        while (std::abs(dS) >= opts.dS_min) {
            Eigen::VectorXd X_pred = X + dXdS * dS;
            const double S_new = X_pred[static_cast<Eigen::Index>(ns)];
            X_new = X_pred;
            cr = corrector(sys, X_new, ns, S_new, F_new, J_new, opts);
            reject.clear();
            if (!cr.ok) {
                reject = cr.message;
            } else {
                const double moved = (X_new - X_pred).cwiseAbs().maxCoeff();
                const double allowed = std::max(2.0 * std::abs(dS) * dXdS.cwiseAbs().maxCoeff(), 1e-6);
                const double maxlnK_new = X_new.head(static_cast<Eigen::Index>(N)).cwiseAbs().maxCoeff();
                sys.unpack(X_new, pt_new);
                if (!ValidNumber(pt_new.p) || pt_new.p <= 0 || !ValidNumber(pt_new.T)) {
                    reject = "invalid pressure or temperature";
                } else if (moved > allowed) {
                    reject = format("corrector moved %g, predictor step allows %g", moved, allowed);
                } else if (!sys.is_lnK(ns) && maxlnK_new < 1e-8) {
                    reject = "trivial solution";
                } else if (std::abs(pt_new.rho_inc / pt_new.rho_feed - 1) < opts.merge_rho_tol && maxlnK_new > opts.merge_lnK_tol) {
                    // The two phases have the same density but different compositions.  At a
                    // genuine critical point both merge together, so density equality with a
                    // large ln K spread is a spurious root, not a point on the boundary.
                    // Without this the N2/CH4/C2/C3 bubble branch drops from 4.3 MPa to 800 Pa
                    // at constant T onto a nonsense state, which then trips the closure test.
                    reject = format("densities merged (%g vs %g) while max|lnK| = %g", pt_new.rho_inc, pt_new.rho_feed, maxlnK_new);
                } else if (!prev_tangent.isZero() && (X_new - X).dot(prev_tangent) <= 0) {
                    // The corrector landed behind the current point along the direction of
                    // travel: the predictor overshot a turning point and the Newton solve fell
                    // back onto the stretch already traced.  Without this the trace happily
                    // retraces itself for hundreds of points (seen on N2/CH4/C2/C3 below 50 K,
                    // where the EOS is extrapolating and the curve folds back).
                    reject = "corrector reversed the direction of travel";
                }
            }
            if (reject.empty()) {
                accepted = true;
                break;
            }
            if (debug) {
                std::cout << format("PhaseEnvelopeTracers: reject dS = %g (ns = %d): %s\n", dS, static_cast<int>(ns), reject.c_str());
            }
            dS *= 0.5;
        }
        if (!accepted) {
            g_stop_reason = "stalled";
            stop_reason = "step size below minimum after repeated rejections: " + reject;
            break;
        }

        // Critical crossing: every ln K changed sign
        bool crossed = true;
        for (std::size_t i = 0; i < N && crossed; ++i) {
            const auto ii = static_cast<Eigen::Index>(i);
            crossed = X[ii] * X_new[ii] < 0;
        }
        if (crossed && !crossed_any) {
            crossed_any = true;
            env.icrit = env.T.size();
        }

        prev_tangent = dXdS * (dS >= 0 ? 1.0 : -1.0);
        X = X_new;
        F = F_new;
        J = J_new;
        iters_prev = cr.iters;
        store(pt_new);
        sys.set_p_ref(pt_new.p);
        if (debug) {
            std::cout << format("PhaseEnvelopeTracers: n = %d ns = %d dS = %g iters = %d T = %g p = %g rho_inc = %g rho_feed = %g max|lnK| = %g\n",
                                static_cast<int>(env.T.size()), static_cast<int>(ns), dS, cr.iters, pt_new.T, pt_new.p, pt_new.rho_inc,
                                pt_new.rho_feed, pt_new.max_abs_lnK);
        }

        // Termination
        const std::size_t npts = env.T.size();
        const double xmax = *std::max_element(pt_new.x_inc.begin(), pt_new.x_inc.end());
        // Closure means the trace came back round to a low-pressure point on the other
        // branch, where the incipient phase is again dilute relative to the feed.  Pressure
        // alone is not enough: a collapse onto a degenerate root also drives p down, and
        // calling that "closed" reports a wrong envelope as a good one.
        const double rho_ratio = std::max(pt_new.rho_inc, pt_new.rho_feed) / std::min(pt_new.rho_inc, pt_new.rho_feed);
        if (crossed_any && npts > 5 && pt_new.p < p_start && rho_ratio > opts.closure_rho_ratio) {
            env.closed = true;
            g_stop_reason = "closed";
            stop_reason = "closed";
            break;
        }
        if (T_floor > 0 && crossed_any && pt_new.T < T_floor) {
            g_stop_reason = "floor";
            stop_reason = format("temperature %g K below floor %g K", pt_new.T, T_floor);
            break;
        }
        if (pt_new.p > opts.p_ceiling) {
            g_stop_reason = "ceiling";
            stop_reason = format("pressure above ceiling %g Pa", opts.p_ceiling);
            break;
        }
        if (xmax > 1 - 1e-9) {
            g_stop_reason = "pure";
            stop_reason = "incipient phase reached a pure component";
            break;
        }
        if (pt_new.max_abs_lnK > lnK_limit) {
            g_stop_reason = "degenerate";
            stop_reason = format("max|lnK| = %g exceeds %g; the incipient phase is numerically pure", pt_new.max_abs_lnK, lnK_limit);
            break;
        }
        if (npts >= opts.max_points) {
            g_stop_reason = "max_points";
            stop_reason = "maximum number of points";
            break;
        }
    }
    g_stop_detail = stop_reason;

    if (debug) {
        std::cout << format("PhaseEnvelopeTracers: stopped after %d points: %s\n", static_cast<int>(env.T.size()), stop_reason.c_str());
    }
    if (env.T.size() < opts.min_points_for_built) {
        throw ValueError(format("PhaseEnvelopeTracers: only %d points traced (%s)", static_cast<int>(env.T.size()), stop_reason.c_str()));
    }
    env.built = true;
}

void PhaseEnvelopeTracers::trace(HelmholtzEOSMixtureBackend& HEOS, const std::string& algorithm, const std::string& level) {
    std::unique_ptr<IsoplethSystem> sys = make_system(HEOS, algorithm);
    run(HEOS, *sys, options_for_level(level));
}

} /* namespace CoolProp */
