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
    // R_D is shipped for completeness but deliberately UNUSED: Li 2024's own code overwrites
    // it with a flat 1.02 for every fluid (code_SI.py:65) and never reads this column.
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

/// Which viscosity the Olchowy-Sengers critical enhancement consumes.
///
/// The enhancement term contains a viscosity.  Li 2024's reference implementation feeds it
/// REFPROP's NATIVE viscosity (code_SI.py:116, `PropsSI('V', T, Dmass, 'REFPROP::'+fluid)`) --
/// NOT the RES viscosity it is in the middle of computing.  Using the RES viscosity is more
/// self-consistent, but it does not reproduce the published values: on the REFPROP backend, where
/// everything else agrees to ~0.01%, that one choice is the entire remaining gap for PROPANE
/// (-4.9%) and R143A (+1.4%).
enum RESEnhancementViscosity
{
    RES_ENH_VIS_AUTO = 0,       ///< the backend's own viscosity if it has one, else the RES viscosity
    RES_ENH_VIS_RES,            ///< always the RES viscosity (self-consistent, diverges from the paper)
    RES_ENH_VIS_BACKEND_NATIVE  ///< always the backend's own viscosity; throws if it has none
};

/// Whether the RES critical enhancement is applied to MIXTURES.
///
/// Li 2024 does apply it (code_SI.py::Olchowy_critical_enhancement_mix), and reproducing the
/// published mixture values needs it.  It is nevertheless NOT enabled by default on CoolProp's
/// own backends: the enhancement requires the MIXTURE critical point, which those backends have
/// to SOLVE for -- not robustly, and not quickly -- while the physical case for a critical
/// enhancement in mixtures is not well established in the first place.  REFPROP reports its
/// mixture critical point directly (CRITPdll), so AUTO enables it there and nowhere else.
///
/// Pure fluids are unaffected by this setting; their enhancement is always applied.
enum RESMixtureEnhancement
{
    RES_MIX_ENH_AUTO = 0,  ///< on where the backend reports a critical point directly, off elsewhere
    RES_MIX_ENH_OFF,       ///< never applied to mixtures
    RES_MIX_ENH_ON         ///< always applied; may be slow, and surfaces critical-point solver failures
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
