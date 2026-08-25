/*!
Residual Entropy Scaling (RES) transport properties.

Deliberately backend-neutral: this TU includes only AbstractState.h and RESTransport.h, NOT any
concrete backend.  That is what lets the REFPROP backend evaluate RES without pulling in the
Helmholtz backend, and it is why these routines live here rather than in TransportRoutines.cpp
(whose header includes HelmholtzEOSMixtureBackend.h).

Parameters come from AbstractState::RES_data(); the only state-dependent inputs are ordinary
AbstractState calls plus the single narrow hook calc_drhomass_dp_constT_at().
*/

#include "Backends/RES/RESTransport.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"  // format()

#include <cmath>
#include <vector>

namespace CoolProp {

// ─── Residual Entropy Scaling (RES) ──────────────────────────────────────────
// Physical constants (2018 CODATA)
namespace {
static const double RES_N_A = 6.02214076e23;  // mol^-1
static const double RES_k_B = 1.380649e-23;   // J/K
// Exponent vectors are hardcoded per the fitting scheme of the literature.
static const double vis_pow[3] = {1.8, 2.4, 2.8};
static const double tc_pow[4] = {1.0, 1.5, 2.0, 2.5};
}  // namespace

/// s_plus = -s_res/R is the model's independent variable, and every term raises it to a FRACTIONAL
/// power, so a non-positive or non-finite value yields NaN rather than an error -- which then
/// leaves the routine as a plausible-looking transport value.  Backends do produce such states:
/// the Peng-Robinson backend returns a NaN residual entropy for hydrogen forced onto its liquid
/// root at 150 K and 112 MPa.  Refuse it here rather than let a NaN escape to the caller.
static void check_s_plus(double s_plus, const char* prop) {
    if (!ValidNumber(s_plus) || s_plus <= 0) {
        throw ValueError(format("RES %s: the residual entropy is not usable at this state "
                                "(s+ = -s_res/R = %g); RES requires s+ > 0.",
                                prop, s_plus));
    }
}

/// Decide where the dilute-gas term comes from, and fetch it when the answer is "the backend".
/// Returns true and sets `value` (already the MIXTURE value) when the backend's own transport
/// model supplied it; returns false when the caller should use the fitted polynomial, plus the
/// Wilke rule for mixtures.
///
/// RES_DILUTE_AUTO reproduces the two source papers, which chose per CODE PATH and not per fluid:
///   Martinek 2025 (viscosity)    pure -> the backend's own eta0     mixture -> polynomial+Wilke
///   Li 2024      (conductivity)  pure -> polynomial                 mixture -> the backend's own lambda0
/// Both obtain that value by flashing to a near-zero pressure (1e-9 Pa and 0.1 Pa respectively)
/// and reading the ordinary transport property there, and both fall back to the polynomial if the
/// call fails.  AUTO does the same.
///
/// An EXPLICIT RES_DILUTE_BACKEND_NATIVE does NOT fall back.  Quietly substituting a different
/// model than the one the caller asked for is exactly how a confidently wrong number ships.
static bool dilute_from_backend(AbstractState& HEOS, parameters key, RESDiluteSource src, bool is_pure, double& value) {
    const bool want_native = (src == RES_DILUTE_BACKEND_NATIVE) || (src == RES_DILUTE_AUTO && ((key == iviscosity) ? is_pure : !is_pure));
    if (!want_native) {
        return false;
    }
    if (HEOS.transport_native(key, 0.0, value)) {
        return true;
    }
    if (src == RES_DILUTE_BACKEND_NATIVE) {
        throw ValueError(format("RES: the dilute-gas term for %s was set to RES_DILUTE_BACKEND_NATIVE, but %s does "
                                "not supply a native transport model to RES for this fluid.",
                                (key == iviscosity) ? "viscosity" : "conductivity", HEOS.backend_name().c_str()));
    }
    return false;
}

/// Wilke (1950) mixing rule for a vector of pure-component transport properties.
/// phi[i][j] = (1 + sqrt(eta0[i]/eta0[j]) * (M[j]/M[i])^0.25)^2 / sqrt(8*(1+M[i]/M[j]))
static double wilke_mix(const std::vector<double>& eta0, const std::vector<double>& M, const std::vector<double>& x) {
    const std::size_t n = eta0.size();
    if (n == 1) return eta0[0];
    double result = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double denom = 0;
        for (std::size_t j = 0; j < n; ++j) {
            double ratio = eta0[i] / eta0[j];
            double Mratio = M[j] / M[i];
            double phi_ij = POW2(1.0 + sqrt(ratio) * pow(Mratio, 0.25)) / sqrt(8.0 * (1.0 + M[i] / M[j]));
            denom += x[j] * phi_ij;
        }
        result += x[i] * eta0[i] / denom;
    }
    return result;
}

CoolPropDbl RESTransport::viscosity(AbstractState& HEOS) {
    const std::vector<CoolPropDbl>& z = HEOS.get_mole_fractions();
    const std::size_t N = z.size();
    const double R = HEOS.gas_constant();
    const double T = HEOS.T();
    const double rho = HEOS.rhomass();       // kg/m³
    const double M_mix = HEOS.molar_mass();  // kg/mol
    const double s_res = HEOS.smolar_residual();
    const double s_plus = -s_res / R;
    const double rhoN = rho / M_mix * RES_N_A;  // number density [m^-3]
    check_s_plus(s_plus, "viscosity");

    // Guard: alpha-function consistency
    for (std::size_t i = 0; i < N; ++i) {
        if (!HEOS.RES_data().comps[i].viscosity.provided)
            throw ValueError(format("Viscosity RES parameters are not available for component '%s'. "
                                    "Ensure RES parameters are loaded or call set_viscosity_RES_parameters().",
                                    HEOS.RES_data().comps[i].name.c_str()));
        if (!HEOS.RES_data().comps[i].viscosity.n_params_match_alpha)
            throw ValueError(format("Viscosity RES n-coefficients for component '%s' were fitted for a different "
                                    "alpha function. Call set_viscosity_RES_residual_params(%zu, n_res, xita) "
                                    "with parameters refitted for your alpha function.",
                                    HEOS.RES_data().comps[i].name.c_str(), i));
    }

    // Dilute-gas viscosity per component [Pa·s] from polynomial in T
    std::vector<double> eta0(N), M(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& d = HEOS.RES_data().comps[i].viscosity;
        const auto& c = d.n_dilute;  // [n0..n4], result in µPa·s
        double v = c[0] + T * (c[1] + T * (c[2] + T * (c[3] + T * c[4])));
        eta0[i] = v * 1e-6;  // µPa·s → Pa·s
        M[i] = d.molar_mass;
    }

    std::vector<double> xd(N);
    for (std::size_t i = 0; i < N; ++i)
        xd[i] = z[i];
    double eta0_mix = 0;
    if (!dilute_from_backend(HEOS, iviscosity, HEOS.RES_data().viscosity_dilute_source, N == 1, eta0_mix)) {
        eta0_mix = wilke_mix(eta0, M, xd);
    }

    // Mixture RES parameters: n_mix[k] = Σ x_i * n_i[k] / xita_i^pow[k]; xita_mix = 1
    double n_mix[3] = {0, 0, 0};
    if (N == 1) {
        const auto& d = HEOS.RES_data().comps[0].viscosity;
        // Pure fluid: n[k]/xita^pow[k] * s_plus^pow[k] == n[k] * (s_plus/xita)^pow[k],
        // so fold xita into s_plus rather than into the coefficients.
        double s_eff = s_plus / d.xita;
        double sum = 0;
        for (int k = 0; k < 3; ++k)
            sum += d.n_res[k] * pow(s_eff, vis_pow[k]);
        double vis_plus = exp(sum) - 1.0;
        double m = M[0] / RES_N_A;
        double vis_res = vis_plus / pow(s_plus, 2.0 / 3.0) * pow(rhoN, 2.0 / 3.0) * sqrt(m * RES_k_B * T);
        return eta0_mix + vis_res;
    }
    // Mixture path
    for (std::size_t i = 0; i < N; ++i) {
        const auto& d = HEOS.RES_data().comps[i].viscosity;
        for (int k = 0; k < 3; ++k)
            n_mix[k] += xd[i] * d.n_res[k] / pow(d.xita, vis_pow[k]);
    }
    double sum = 0;
    for (int k = 0; k < 3; ++k)
        sum += n_mix[k] * pow(s_plus, vis_pow[k]);
    double vis_plus = exp(sum) - 1.0;

    // Mixture effective molecular mass (mole-fraction weighted M for dilute ideal-gas limit)
    double M_eff = 0;
    for (std::size_t i = 0; i < N; ++i)
        M_eff += xd[i] * M[i];
    double m_mix = M_eff / RES_N_A;

    double vis_res = vis_plus / pow(s_plus, 2.0 / 3.0) * pow(rhoN, 2.0 / 3.0) * sqrt(m_mix * RES_k_B * T);
    return eta0_mix + vis_res;
}

CoolPropDbl RESTransport::conductivity(AbstractState& HEOS) {
    const std::vector<CoolPropDbl>& z = HEOS.get_mole_fractions();
    const std::size_t N = z.size();
    const double R = HEOS.gas_constant();
    const double T = HEOS.T();
    const double rho = HEOS.rhomass();
    const double M_mix = HEOS.molar_mass();
    const double s_res = HEOS.smolar_residual();
    const double s_plus = -s_res / R;
    const double rhoN = rho / M_mix * RES_N_A;
    check_s_plus(s_plus, "conductivity");

    for (std::size_t i = 0; i < N; ++i) {
        if (!HEOS.RES_data().comps[i].conductivity.provided)
            throw ValueError(format("Conductivity RES parameters are not available for component '%s'. "
                                    "Ensure RES parameters are loaded or call set_conductivity_RES_parameters().",
                                    HEOS.RES_data().comps[i].name.c_str()));
        if (!HEOS.RES_data().comps[i].conductivity.n_params_match_alpha)
            throw ValueError(format("Conductivity RES n-coefficients for component '%s' were fitted for a different "
                                    "alpha function. Call set_conductivity_RES_residual_params(%zu, n_res, xita) "
                                    "with parameters refitted for your alpha function.",
                                    HEOS.RES_data().comps[i].name.c_str(), i));
    }

    // Dilute-gas thermal conductivity per component from polynomial [W/m/K]
    std::vector<double> tc0(N), M(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& d = HEOS.RES_data().comps[i].conductivity;
        const auto& c = d.n_dilute;  // result in W/m/K
        tc0[i] = c[0] + T * (c[1] + T * (c[2] + T * (c[3] + T * c[4])));
        M[i] = d.molar_mass;
    }

    std::vector<double> xd(N);
    for (std::size_t i = 0; i < N; ++i)
        xd[i] = z[i];
    double tc0_mix = 0;
    if (!dilute_from_backend(HEOS, iconductivity, HEOS.RES_data().conductivity_dilute_source, N == 1, tc0_mix)) {
        tc0_mix = wilke_mix(tc0, M, xd);
    }

    // Mass-fraction weighted sqrt(m) for effective molecular mass (Yang 2024)
    double M_mix_actual = 0;
    for (std::size_t i = 0; i < N; ++i)
        M_mix_actual += xd[i] * M[i];
    double sqrt_m_mix = 0;
    for (std::size_t i = 0; i < N; ++i) {
        double w_i = xd[i] * M[i] / M_mix_actual;  // mass fraction
        sqrt_m_mix += w_i * sqrt(M[i] / RES_N_A);
    }
    double m_mix = sqrt_m_mix * sqrt_m_mix;  // (Σ w_i √m_i)²

    // Mixture RES parameters
    double n_mix[4] = {0, 0, 0, 0};
    double tc_plus;
    if (N == 1) {
        const auto& d = HEOS.RES_data().comps[0].conductivity;
        double s_eff = s_plus / d.xita;
        tc_plus = 0;
        for (int k = 0; k < 4; ++k)
            tc_plus += d.n_res[k] * pow(s_eff, tc_pow[k]);
    } else {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& d = HEOS.RES_data().comps[i].conductivity;
            for (int k = 0; k < 4; ++k)
                n_mix[k] += xd[i] * d.n_res[k] / pow(d.xita, tc_pow[k]);
        }
        tc_plus = 0;
        for (int k = 0; k < 4; ++k)
            tc_plus += n_mix[k] * pow(s_plus, tc_pow[k]);
    }

    double tc_res = tc_plus / pow(s_plus, 2.0 / 3.0) * RES_k_B * pow(rhoN, 2.0 / 3.0) * sqrt(RES_k_B * T / m_mix);

    // Critical enhancement (pure fluid only, requires fitted Li 2024 parameters).
    // The guards below mirror Olchowy_critical_enhancement() in the Li 2024 supporting-information
    // code (code_SI.py), which returns exactly 0 outside the near-critical region.  Without them
    // the enhancement is applied in dense-liquid / high-temperature states where the reference
    // implementation suppresses it.
    double delta_tc_crit = 0;
    {
        // Mole-fraction-mixed enhancement parameters.  For a pure fluid this collapses to that
        // one component's values.  Li 2024 mixes linearly (code_SI.py:180-186) and mixes 1/q_D
        // rather than q_D -- its parameter table stores qDinv -- so mix the reciprocal here too.
        //
        // Some fluids (D2O, HELIUM, ORTHOHYD) carry an all-zero critical-enhancement record in
        // the source tables, meaning "no parameters fitted" rather than "no enhancement".  Every
        // component must have them, otherwise a zero would silently dilute the mixed values.
        //
        // R_D is NOT taken from the parameter table.  Li's get_paramters() overwrites it with a
        // flat 1.02 for every fluid (code_SI.py:65) and never reads the table's R_D column, so
        // that column is dead data -- ConductivityRESData::R_D carries it only for completeness.
        // Using it instead moves PROPANE, whose enhancement is ~32 % of lambda, by +0.3 %.
        // gamma, by contrast, IS read per fluid (code_SI.py:66), and takes two distinct values.
        const double R_D = 1.02;
        bool crit_usable = N > 0;
        double gamma_e = 0, phi0 = 0, Gamma = 0, qD_inv = 0, t_ref_pure = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const ConductivityRESData& cd = HEOS.RES_data().comps[i].conductivity;
            if (!cd.crit_provided || !(cd.Gamma > 0) || !(cd.phi0 > 0) || !(cd.q_D > 0) || !ValidNumber(cd.gamma_uni) || !(cd.gamma_uni > 0)) {
                crit_usable = false;
                break;
            }
            gamma_e += xd[i] * cd.gamma_uni;
            phi0 += xd[i] * cd.phi0;
            Gamma += xd[i] * cd.Gamma;
            qD_inv += xd[i] / cd.q_D;
            t_ref_pure = cd.t_ref;
        }
        // A pure fluid takes its reference temperature from the table and needs it to be real;
        // for a mixture Li computes it as 1.5*Tc instead, so there is nothing to check yet.
        if (N == 1 && (!ValidNumber(t_ref_pure) || !(t_ref_pure > 0))) {
            crit_usable = false;
        }

        // Mixtures are gated by policy: the enhancement needs the MIXTURE critical point, and on
        // backends that must solve for it that is slow and unreliable.  See RESMixtureEnhancement.
        const RESMixtureEnhancement mix_pol = HEOS.RES_data().mixture_enhancement;
        if (N > 1 && !(mix_pol == RES_MIX_ENH_ON || (mix_pol == RES_MIX_ENH_AUTO && HEOS.critical_point_is_cheap()))) {
            crit_usable = false;
        }

        if (crit_usable) {
            const double nu = 0.63;  // universal exponent; gamma is per-fluid (mixed above)
            double Tc = 0, pc = 0, rhoc = 0;
            bool have_crit = true;
            if (N == 1) {
                Tc = HEOS.T_critical();
                pc = HEOS.p_critical();
                rhoc = HEOS.rhomass_critical();
            } else {
                // Under AUTO a failed critical-point lookup is not a reason to fail the whole
                // conductivity call -- Li wraps the entire mixture enhancement in try/except and
                // returns zero.  Under an explicit RES_MIX_ENH_ON the caller asked for it, so the
                // failure surfaces instead of being papered over with a silently smaller answer.
                try {
                    if (HEOS.critical_point_is_cheap()) {
                        Tc = HEOS.T_critical();
                        pc = HEOS.p_critical();
                        rhoc = HEOS.rhomass_critical();
                    } else {
                        // Where the critical point has to be SOLVED for, those three accessors run
                        // the solve three times over for one answer -- on HEOS that is ~460 ms
                        // each.  Ask once instead.  Requiring exactly one point mirrors
                        // HelmholtzEOSMixtureBackend::calc_T_critical(): with several, Tc is
                        // ambiguous and RES has no basis for choosing between them.
                        const std::vector<CriticalState> pts = HEOS.all_critical_points();
                        if (pts.size() != 1) {
                            throw ValueError(format("RES: the mixture critical enhancement needs a single critical "
                                                    "point, but the backend reported %d.",
                                                    static_cast<int>(pts.size())));
                        }
                        Tc = pts[0].T;
                        pc = pts[0].p;
                        rhoc = pts[0].rhomolar * HEOS.molar_mass();
                    }
                } catch (...) {
                    if (mix_pol == RES_MIX_ENH_ON) {
                        throw;
                    }
                    have_crit = false;
                }
            }
            if (have_crit && ValidNumber(Tc) && Tc > 0 && ValidNumber(rhoc) && rhoc > 0 && ValidNumber(pc) && pc > 0) {
                const double rho_mass = rho;
                const double delta_r = rho_mass / rhoc;
                // Li 2024 takes the mixture reference temperature as 1.5*Tc rather than from the
                // per-fluid table, which has no entry for a mixture.
                const double t_ref = (N == 1) ? t_ref_pure : 1.5 * Tc;

                // Outside the near-critical region the enhancement is identically zero (Li 2024).
                // Also skip when rho == 0, where delta_r would make Omega0 undefined.
                if (rho_mass <= 0 || delta_r >= 2.0 || T / Tc > 1.4) return tc0_mix + tc_res;

                // dp/drho at current T and at t_ref (using same delta, tau_ref)
                const double delta_st = HEOS.delta();  // rho/rhoc (molar)
                double dp_drho =
                  HEOS.gas_constant() * T * (1.0 + 2.0 * delta_st * HEOS.dalphar_dDelta() + delta_st * delta_st * HEOS.d2alphar_dDelta2());
                // dp/drho = dp_drho_molar * M  when working in mass density; use molar base
                double drhodp_t = 1.0 / dp_drho * HEOS.molar_mass();  // (kg/m³)/Pa

                double drhodp_tref = HEOS.drhomass_dp_constT_at(t_ref);

                double arg = drhodp_t - (t_ref / T) * drhodp_tref;
                // The enhancement term needs a viscosity, and which one is a real choice -- see
                // RESEnhancementViscosity.  Under AUTO we take the backend's own viscosity at the
                // current density, as Li 2024 does; only when the backend has none do we fall back
                // to the RES viscosity.  That fallback needs every component to carry viscosity
                // RES parameters valid for the current alpha function, and a fluid can legitimately
                // have conductivity parameters but not viscosity ones -- so skip the enhancement
                // rather than let viscosity() throw out of a conductivity call.
                const RESEnhancementViscosity evs = HEOS.RES_data().conductivity_enhancement_viscosity;
                double vis_native = 0;
                bool have_native = false;
                if (evs != RES_ENH_VIS_RES) {
                    have_native = HEOS.transport_native(iviscosity, HEOS.rhomolar(), vis_native) && vis_native > 0;
                    if (!have_native && evs == RES_ENH_VIS_BACKEND_NATIVE) {
                        throw ValueError(format("RES: the critical-enhancement viscosity was set to RES_ENH_VIS_BACKEND_NATIVE, "
                                                "but %s does not supply a native viscosity model to RES for this fluid.",
                                                HEOS.backend_name().c_str()));
                    }
                }
                bool vis_res_usable = true;
                for (std::size_t i = 0; i < N && vis_res_usable; ++i) {
                    const ViscosityRESData& vd = HEOS.RES_data().comps[i].viscosity;
                    vis_res_usable = vd.provided && vd.n_params_match_alpha;
                }
                if (arg > 0 && (have_native || vis_res_usable)) {
                    double phi = phi0 * pow((pc * rho_mass) / (Gamma * rhoc * rhoc), nu / gamma_e) * pow(arg, nu / gamma_e);
                    double y = phi / qD_inv;
                    double kappa_cv = HEOS.cpmass() / HEOS.cvmass();
                    double Omega = 2.0 / M_PI * ((1.0 - 1.0 / kappa_cv) * atan(y) + y / kappa_cv);
                    double Omega0 = 2.0 / M_PI * (1.0 - exp(-1.0 / (1.0 / y + POW2(y / delta_r) / 3.0)));
                    double vis = have_native ? vis_native : RESTransport::viscosity(HEOS);
                    if (phi > 0 && vis > 0) {
                        delta_tc_crit = RES_k_B * R_D * rho_mass * HEOS.cpmass() * T / (vis * phi * 6.0 * M_PI) * (Omega - Omega0);
                    }
                    // Li 2024 clamps a negative enhancement to zero rather than subtracting it.
                    if (!(delta_tc_crit > 0)) {
                        delta_tc_crit = 0;
                    }
                }
            }
        }
    }

    return tc0_mix + tc_res + delta_tc_crit;
}

// ─── End of RES section ───────────────────────────────────────────────────────

}  // namespace CoolProp
