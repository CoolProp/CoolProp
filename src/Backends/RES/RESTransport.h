#ifndef COOLPROP_RES_RESTRANSPORT_H
#define COOLPROP_RES_RESTRANSPORT_H

/*!
Residual Entropy Scaling (RES) transport model.

Backend-neutral by construction: this header forward-declares AbstractState rather than
including any concrete backend, so REFPROP can call these without pulling in the Helmholtz
backend.  Parameters are read from AbstractState::RES_data().
*/

#include "CoolProp/detail/tools.h"  // CoolPropDbl

namespace CoolProp {

class AbstractState;

namespace RESTransport {

/// RES viscosity in Pa*s.  Works for pure fluids and mixtures.
/// Throws ValueError naming the component if its parameters are missing, or if they were fitted
/// for a different alpha function.
CoolPropDbl viscosity(AbstractState& HEOS);

/// RES thermal conductivity in W/m/K, including the Li 2024 Olchowy-Sengers critical
/// enhancement for pure fluids.  Calls viscosity() internally for the enhancement term.
CoolPropDbl conductivity(AbstractState& HEOS);

}  // namespace RESTransport
}  // namespace CoolProp

#endif  // COOLPROP_RES_RESTRANSPORT_H
