#ifndef COOLPROP_RESTRANSPORT_H
#define COOLPROP_RESTRANSPORT_H

/*!
Residual Entropy Scaling (RES) transport parameter storage.

This header is deliberately lightweight -- it pulls in nothing but <string>, <vector> and the
_HUGE definition -- because AbstractState.h includes it in order to own an RESTransportStore.
Do NOT add heavy includes here (CoolPropFluid.h drags in Eigen and the superancillary machinery
and must not become reachable from AbstractState.h).

The parameter records themselves live here rather than in CoolPropFluid.h so that a backend
which has no CoolPropFluid at all -- notably REFPROP -- can still carry RES parameters.
CoolPropFluid.h includes this header, so `fluid.transport.viscosity_res` keeps working
unchanged.
*/

#include "CoolProp/numerics/numerics.h"

#include <string>
#include <vector>

namespace CoolProp {

// Residual Entropy Scaling (RES) transport data. Martinek 2025 / Yang 2021-2025.
struct ViscosityRESData
{
    std::vector<double> n_dilute;      // 5 polynomial coefficients: eta0(T)/uPas = sum_k n_k * T^k
    std::vector<double> n_res;         // 3 coefficients for ln(vis_plus + 1)
    double xita = 1.0;                 // scaling factor xi; 1.0 when individual fit was used
    int group_num = -1;
    double molar_mass = _HUGE;         // kg/mol, cached for Wilke mixing
    bool n_params_match_alpha = true;  // false after alpha-function change without new n_res
    bool provided = false;
};

// Residual Entropy Scaling (RES) thermal conductivity data. Yang 2021-2025 / Li 2024.
struct ConductivityRESData
{
    std::vector<double> n_dilute;      // 5 polynomial coefficients: lambda0(T)/Wm-1K-1
    std::vector<double> n_res;         // 4 coefficients for tc_plus
    double xita = 1.0;
    int group_num = -1;
    double molar_mass = _HUGE;
    // Olchowy-Sengers critical enhancement parameters (Li 2024 parameterization)
    double R_D = _HUGE, gamma_uni = _HUGE, Gamma = _HUGE;
    double phi0 = _HUGE, t_ref = _HUGE, q_D = 0.0;
    bool n_params_match_alpha = true;
    bool crit_provided = false;
    bool provided = false;
};

/// Where the dilute-gas term of the RES model should come from.
///
/// The two source papers made OPPOSITE choices, per code path rather than per fluid:
/// Martinek 2025 uses REFPROP's native eta0 for pure fluids and the fitted polynomial with a
/// Wilke rule for mixtures; Li 2024 uses the polynomial for pure fluids and REFPROP's native
/// lambda0 for mixtures.  AUTO reproduces that, resolving to POLYNOMIAL on any backend that has
/// no native dilute-gas model -- which keeps HEOS and the cubics numerically unchanged.
enum RESDiluteSource
{
    RES_DILUTE_AUTO = 0,
    RES_DILUTE_POLYNOMIAL,
    RES_DILUTE_BACKEND_NATIVE
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
};

}  // namespace CoolProp

#endif  // COOLPROP_RESTRANSPORT_H
