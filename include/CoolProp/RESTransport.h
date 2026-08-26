#ifndef COOLPROP_RESTRANSPORT_H
#define COOLPROP_RESTRANSPORT_H

/*!
Residual Entropy Scaling (RES) transport parameter storage.

This header is deliberately lightweight -- <string>, <vector>, the _HUGE definition and the
enums in DataStructures.h -- because AbstractState.h includes it in order to own an
RESTransportStore.
Do NOT add heavy includes here (CoolPropFluid.h drags in Eigen and the superancillary machinery
and must not become reachable from AbstractState.h).

The parameter records themselves live here rather than in CoolPropFluid.h so that a backend
which has no CoolPropFluid at all -- notably REFPROP -- can still carry RES parameters.
CoolPropFluid.h includes this header, so `fluid.transport.viscosity_res` keeps working
unchanged.
*/

#include "CoolProp/DataStructures.h"
#include "CoolProp/numerics/numerics.h"

#include <cstddef>
#include <string>
#include <vector>

namespace CoolProp {

/// The RES model has a FIXED number of coefficients per property.  RESTransport.cpp indexes
/// these positions unconditionally once `provided` is set, so the setters must reject anything
/// shorter -- see AbstractState::check_RES_coeff_sizes().  Keep in step with vis_pow[] / tc_pow[]
/// in src/Backends/RES/RESTransport.cpp, which carry the matching exponents.
constexpr std::size_t RES_N_DILUTE = 5;            ///< dilute-gas polynomial, n0..n4
constexpr std::size_t RES_N_RES_VISCOSITY = 3;     ///< residual viscosity coefficients
constexpr std::size_t RES_N_RES_CONDUCTIVITY = 4;  ///< residual conductivity coefficients

// Residual Entropy Scaling (RES) transport data. Martinek 2025 / Yang 2021-2025.
struct ViscosityRESData
{
    std::vector<double> n_dilute;  // 5 polynomial coefficients: eta0(T)/uPas = sum_k n_k * T^k
    std::vector<double> n_res;     // 3 coefficients for ln(vis_plus + 1)
    double xita = 1.0;             // scaling factor xi; 1.0 when individual fit was used
    int group_num = -1;
    double molar_mass = _HUGE;         // kg/mol, cached for Wilke mixing
    bool n_params_match_alpha = true;  // false after alpha-function change without new n_res
    bool provided = false;
};

// Residual Entropy Scaling (RES) thermal conductivity data. Yang 2021-2025 / Li 2024.
struct ConductivityRESData
{
    std::vector<double> n_dilute;  // 5 polynomial coefficients: lambda0(T)/Wm-1K-1
    std::vector<double> n_res;     // 4 coefficients for tc_plus
    double xita = 1.0;
    int group_num = -1;
    double molar_mass = _HUGE;
    // Olchowy-Sengers critical enhancement parameters (Li 2024 parameterization)
    // R_D is shipped for completeness but deliberately UNUSED: Li 2024's own code overwrites
    // it with a flat 1.02 for every fluid (code_SI.py:65) and never reads this column.
    double R_D = _HUGE, gamma_uni = _HUGE, Gamma = _HUGE;
    double phi0 = _HUGE, t_ref = _HUGE, q_D = 0.0;
    bool n_params_match_alpha = true;
    bool crit_provided = false;
    bool provided = false;
};

/// RES parameters for one component, independent of any equation of state.
struct RESComponentData
{
    std::string name;  ///< used only for error messages naming the offending component
    ViscosityRESData viscosity;
    ConductivityRESData conductivity;
};

/// The live, per-instance RES state owned by AbstractState.
///
/// The records in CoolPropFluid::transport are the *shipped* parameters that seed this store;
/// this is the mutable copy the setters act on and the transport routines read.
struct RESTransportStore
{
    std::vector<RESComponentData> comps;
    bool viscosity_enabled = false;
    bool conductivity_enabled = false;
    RESDiluteSource viscosity_dilute_source = RES_DILUTE_AUTO;
    RESDiluteSource conductivity_dilute_source = RES_DILUTE_AUTO;
    RESEnhancementViscosity conductivity_enhancement_viscosity = RES_ENH_VIS_AUTO;
    RESMixtureEnhancement mixture_enhancement = RES_MIX_ENH_AUTO;
};

}  // namespace CoolProp

#endif  // COOLPROP_RESTRANSPORT_H
