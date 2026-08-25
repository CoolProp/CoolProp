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
static const double RES_N_A  = 6.02214076e23;   // mol^-1
static const double RES_k_B  = 1.380649e-23;    // J/K
// Exponent vectors are hardcoded per the fitting scheme of the literature.
static const double vis_pow[3]  = {1.8, 2.4, 2.8};
static const double tc_pow[4]   = {1.0, 1.5, 2.0, 2.5};
}  // namespace

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
    const double R      = HEOS.gas_constant();
    const double T      = HEOS.T();
    const double rho    = HEOS.rhomass();   // kg/m³
    const double M_mix  = HEOS.molar_mass();  // kg/mol
    const double s_res  = HEOS.smolar_residual();
    const double s_plus = -s_res / R;
    const double rhoN   = rho / M_mix * RES_N_A;  // number density [m^-3]

    // Guard: alpha-function consistency
    for (std::size_t i = 0; i < N; ++i) {
        if (!HEOS.RES_data().comps[i].viscosity.provided)
            throw ValueError(format(
                "Viscosity RES parameters are not available for component '%s'. "
                "Ensure RES parameters are loaded or call set_viscosity_RES_parameters().",
                HEOS.RES_data().comps[i].name.c_str()));
        if (!HEOS.RES_data().comps[i].viscosity.n_params_match_alpha)
            throw ValueError(format(
                "Viscosity RES n-coefficients for component '%s' were fitted for a different "
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
        M[i]    = d.molar_mass;
    }

    std::vector<double> xd(N);
    for (std::size_t i = 0; i < N; ++i) xd[i] = z[i];
    double eta0_mix = wilke_mix(eta0, M, xd);

    // Mixture RES parameters: n_mix[k] = Σ x_i * n_i[k] / xita_i^pow[k]; xita_mix = 1
    double n_mix[3] = {0, 0, 0};
    if (N == 1) {
        const auto& d = HEOS.RES_data().comps[0].viscosity;
        // Pure fluid: n[k]/xita^pow[k] * s_plus^pow[k] == n[k] * (s_plus/xita)^pow[k],
        // so fold xita into s_plus rather than into the coefficients.
        double s_eff = s_plus / d.xita;
        double sum = 0;
        for (int k = 0; k < 3; ++k) sum += d.n_res[k] * pow(s_eff, vis_pow[k]);
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
    for (int k = 0; k < 3; ++k) sum += n_mix[k] * pow(s_plus, vis_pow[k]);
    double vis_plus = exp(sum) - 1.0;

    // Mixture effective molecular mass (mole-fraction weighted M for dilute ideal-gas limit)
    double M_eff = 0;
    for (std::size_t i = 0; i < N; ++i) M_eff += xd[i] * M[i];
    double m_mix = M_eff / RES_N_A;

    double vis_res = vis_plus / pow(s_plus, 2.0 / 3.0) * pow(rhoN, 2.0 / 3.0) * sqrt(m_mix * RES_k_B * T);
    return eta0_mix + vis_res;
}

CoolPropDbl RESTransport::conductivity(AbstractState& HEOS) {
    const std::vector<CoolPropDbl>& z = HEOS.get_mole_fractions();
    const std::size_t N = z.size();
    const double R      = HEOS.gas_constant();
    const double T      = HEOS.T();
    const double rho    = HEOS.rhomass();
    const double M_mix  = HEOS.molar_mass();
    const double s_res  = HEOS.smolar_residual();
    const double s_plus = -s_res / R;
    const double rhoN   = rho / M_mix * RES_N_A;

    for (std::size_t i = 0; i < N; ++i) {
        if (!HEOS.RES_data().comps[i].conductivity.provided)
            throw ValueError(format(
                "Conductivity RES parameters are not available for component '%s'. "
                "Ensure RES parameters are loaded or call set_conductivity_RES_parameters().",
                HEOS.RES_data().comps[i].name.c_str()));
        if (!HEOS.RES_data().comps[i].conductivity.n_params_match_alpha)
            throw ValueError(format(
                "Conductivity RES n-coefficients for component '%s' were fitted for a different "
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
        M[i]   = d.molar_mass;
    }

    std::vector<double> xd(N);
    for (std::size_t i = 0; i < N; ++i) xd[i] = z[i];
    double tc0_mix = wilke_mix(tc0, M, xd);

    // Mass-fraction weighted sqrt(m) for effective molecular mass (Yang 2024)
    double M_mix_actual = 0;
    for (std::size_t i = 0; i < N; ++i) M_mix_actual += xd[i] * M[i];
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
        for (int k = 0; k < 4; ++k) tc_plus += d.n_res[k] * pow(s_eff, tc_pow[k]);
    } else {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& d = HEOS.RES_data().comps[i].conductivity;
            for (int k = 0; k < 4; ++k)
                n_mix[k] += xd[i] * d.n_res[k] / pow(d.xita, tc_pow[k]);
        }
        tc_plus = 0;
        for (int k = 0; k < 4; ++k) tc_plus += n_mix[k] * pow(s_plus, tc_pow[k]);
    }

    double tc_res = tc_plus / pow(s_plus, 2.0 / 3.0) * RES_k_B * pow(rhoN, 2.0 / 3.0) * sqrt(RES_k_B * T / m_mix);

    // Critical enhancement (pure fluid only, requires fitted Li 2024 parameters).
    // The guards below mirror Olchowy_critical_enhancement() in the Li 2024 supporting-information
    // code (code_SI.py), which returns exactly 0 outside the near-critical region.  Without them
    // the enhancement is applied in dense-liquid / high-temperature states where the reference
    // implementation suppresses it.
    double delta_tc_crit = 0;
    if (N == 1) {
        const ConductivityRESData& cd = HEOS.RES_data().comps[0].conductivity;
        // Some fluids (D2O, HELIUM, ORTHOHYD) carry an all-zero critical-enhancement record in the
        // source tables, meaning "no parameters fitted".  Li 2024 gates on t_ref > 0; do the same
        // here, since t_ref == 0 and Gamma == 0 would otherwise yield Tc/0 and a division by zero
        // that propagates NaN into the returned conductivity.
        const bool crit_usable = cd.crit_provided && ValidNumber(cd.t_ref) && cd.t_ref > 0 && cd.Gamma > 0 && cd.phi0 > 0 && cd.q_D > 0;
        if (crit_usable) {
            const double nu = 0.63, gamma_e = 1.239;  // universal exponents
            const double R_D = 1.02;  // Li 2024 uses universal values (see Python module comment)
            const double Tc   = HEOS.T_critical();
            const double pc   = HEOS.p_critical();
            const double rhoc = HEOS.rhomass_critical();
            const double rho_mass = rho;
            const double delta_r  = rho_mass / rhoc;
            const double t_ref    = cd.t_ref;
            const double q_D      = cd.q_D;
            const double Gamma    = cd.Gamma;
            const double phi0     = cd.phi0;

            // Outside the near-critical region the enhancement is identically zero (Li 2024).
            // Also skip when rho == 0, where delta_r would make Omega0 undefined.
            if (rho_mass <= 0 || delta_r >= 2.0 || T / Tc > 1.4) return tc0_mix + tc_res;

            // dp/drho at current T and at t_ref (using same delta, tau_ref)
            const double delta_st = HEOS.delta();  // rho/rhoc (molar)
            double dp_drho = HEOS.gas_constant() * T * (1.0 + 2.0 * delta_st * HEOS.dalphar_dDelta()
                                                         + delta_st * delta_st * HEOS.d2alphar_dDelta2());
            // dp/drho = dp_drho_molar * M  when working in mass density; use molar base
            double drhodp_t    = 1.0 / dp_drho * HEOS.molar_mass();  // (kg/m³)/Pa

            double drhodp_tref = HEOS.drhomass_dp_constT_at(t_ref);

            double arg = drhodp_t - (t_ref / T) * drhodp_tref;
            // The enhancement term needs a viscosity.  Li 2024 takes it from the reference EOS; we
            // use the RES viscosity, which requires this component to also carry viscosity RES
            // parameters valid for the current alpha function.  A fluid can legitimately have
            // conductivity RES parameters but not viscosity ones, so skip the enhancement rather
            // than letting viscosity_RES() throw out of a conductivity call.
            const ViscosityRESData& vd = HEOS.RES_data().comps[0].viscosity;
            const bool vis_usable = vd.provided && vd.n_params_match_alpha;
            if (arg > 0 && vis_usable) {
                double phi = phi0 * pow((pc * rho_mass) / (Gamma * rhoc * rhoc), nu / gamma_e) * pow(arg, nu / gamma_e);
                double y   = phi * q_D;
                double kappa_cv = HEOS.cpmass() / HEOS.cvmass();
                double Omega  = 2.0 / M_PI * ((1.0 - 1.0 / kappa_cv) * atan(y) + y / kappa_cv);
                double Omega0 = 2.0 / M_PI * (1.0 - exp(-1.0 / (1.0 / y + POW2(y / delta_r) / 3.0)));
                double vis = RESTransport::viscosity(HEOS);
                if (phi > 0 && vis > 0)
                    delta_tc_crit = RES_k_B * R_D * rho_mass * HEOS.cpmass() * T
                                    / (vis * phi * 6.0 * M_PI) * (Omega - Omega0);
            }
        }
    }

    return tc0_mix + tc_res + delta_tc_crit;
}

// ─── End of RES section ───────────────────────────────────────────────────────

}  // namespace CoolProp
