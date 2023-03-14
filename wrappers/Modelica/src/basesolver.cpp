#include "basesolver.h"
#include <math.h>
#include "CoolPropLib.h"

//! Constructor.
/*!
  The constructor is copying the medium name, library name and substance name
  to the locally defined variables.
  @param mediumName Arbitrary medium name
  @param libraryName Name of the external fluid property library
  @param substanceName Substance name
*/
BaseSolver::BaseSolver(const std::string& mediumName, const std::string& libraryName, const std::string& substanceName)
  : mediumName(mediumName), libraryName(libraryName), substanceName(substanceName) {}

//! Destructor
/*!
  The destructor for the base solver if currently not doing anything.
*/
BaseSolver::~BaseSolver() {}

//! Return molar mass (Default implementation provided)
double BaseSolver::molarMass() const {
    return _fluidConstants.MM;
}

//! Return temperature at critical point (Default implementation provided)
double BaseSolver::criticalTemperature() const {
    return _fluidConstants.Tc;
}

//! Return pressure at critical point (Default implementation provided)
double BaseSolver::criticalPressure() const {
    return _fluidConstants.pc;
}

//! Return molar volume at critical point (Default implementation provided)
double BaseSolver::criticalMolarVolume() const {
    return _fluidConstants.MM / _fluidConstants.dc;
}

//! Return density at critical point (Default implementation provided)
double BaseSolver::criticalDensity() const {
    return _fluidConstants.dc;
}

//! Return specific enthalpy at critical point (Default implementation provided)
double BaseSolver::criticalEnthalpy() const {
    return _fluidConstants.hc;
}

//! Return specific entropy at critical point (Default implementation provided)
double BaseSolver::criticalEntropy() const {
    return _fluidConstants.sc;
}

//! Set fluid constants
/*!
  This function sets the fluid constants which are defined in the
  FluidConstants record in Modelica. It should be called when a new
  solver is created.

  Must be re-implemented in the specific solver
*/
void BaseSolver::setFluidConstants() {}

//! Set state from p, h, and phase
/*!
  This function sets the thermodynamic state record for the given pressure
  p, the specific enthalpy h and the specified phase. The computed values are
  written to the ExternalThermodynamicState property struct.

  Must be re-implemented in the specific solver
  @param p Pressure
  @param h Specific enthalpy
  @param phase Phase (2 for two-phase, 1 for one-phase, 0 if not known)
  @param properties ExternalThermodynamicState property struct
*/
void BaseSolver::setState_ph(double& p, double& h, int& phase, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setState_ph() not implemented in the Solver object");
}

//! Set state from p and T
/*!
  This function sets the thermodynamic state record for the given pressure
  p and the temperature T. The computed values are
  written to the ExternalThermodynamicState property struct.

  Must be re-implemented in the specific solver
  @param p Pressure
  @param T Temperature
  @param properties ExternalThermodynamicState property struct
*/
void BaseSolver::setState_pT(double& p, double& T, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setState_pT() not implemented in the Solver object");
}

//! Set state from d, T, and phase
/*!
  This function sets the thermodynamic state record for the given density
  d, the temperature T and the specified phase. The computed values are
  written to the ExternalThermodynamicState property struct.

  Must be re-implemented in the specific solver
  @param d Density
  @param T Temperature
  @param phase Phase (2 for two-phase, 1 for one-phase, 0 if not known)
  @param properties ExternalThermodynamicState property struct
*/
void BaseSolver::setState_dT(double& d, double& T, int& phase, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setState_dT() not implemented in the Solver object");
}

//! Set state from p, s, and phase
/*!
  This function sets the thermodynamic state record for the given pressure
  p, the specific entropy s and the specified phase. The computed values are
  written to the ExternalThermodynamicState property struct.

  Must be re-implemented in the specific solver
  @param p Pressure
  @param s Specific entropy
  @param phase Phase (2 for two-phase, 1 for one-phase, 0 if not known)
  @param properties ExternalThermodynamicState property struct
*/
void BaseSolver::setState_ps(double& p, double& s, int& phase, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setState_ps() not implemented in the Solver object");
}

//! Compute Prandtl number
/*!
  This function returns the Prandtl number
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::Pr(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: Pr() not implemented in the Solver object");
    return 0;
}

//! Compute temperature
/*!
  This function returns the temperature
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::T(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: T() not implemented in the Solver object");
    return 0;
}

//! Compute velocity of sound
/*!
  This function returns the velocity of sound
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::a(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: a() not implemented in the Solver object");
    return 0;
}

//! Compute isobaric expansion coefficient
/*!
  This function returns the isobaric expansion coefficient
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::beta(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: beta() not implemented in the Solver object");
    return 0;
}

//! Compute specific heat capacity cp
/*!
  This function returns the specific heat capacity cp
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::cp(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: cp() not implemented in the Solver object");
    return 0;
}

//! Compute specific heat capacity cv
/*!
  This function returns the specific heat capacity cv
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::cv(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: cv() not implemented in the Solver object");
    return 0;
}

//! Compute density
/*!
  This function returns the density
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::d(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: d() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of density wrt enthalpy at constant pressure
/*!
  This function returns the derivative of density wrt enthalpy at constant pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::ddhp(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddhp() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of density wrt pressure at constant enthalpy
/*!
  This function returns the derivative of density wrt pressure at constant enthalpy
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::ddph(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddph() not implemented in the Solver object");
    return 0;
}

//! Compute dynamic viscosity
/*!
  This function returns the dynamic viscosity
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::eta(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: eta() not implemented in the Solver object");
    return 0;
}

//! Compute specific enthalpy
/*!
  This function returns the specific enthalpy
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::h(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: h() not implemented in the Solver object");
    return 0;
}

//! Compute compressibility
/*!
  This function returns the compressibility
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::kappa(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: kappa() not implemented in the Solver object");
    return 0;
}

//! Compute thermal conductivity
/*!
  This function returns the thermal conductivity
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::lambda(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: lambda() not implemented in the Solver object");
    return 0;
}

//! Compute pressure
/*!
  This function returns the pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::p(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: p() not implemented in the Solver object");
    return 0;
}

//! Compute phase flag
/*!
  This function returns the phase flag
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
int BaseSolver::phase(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: phase() not implemented in the Solver object");
    return 0;
}

//! Compute specific entropy
/*!
  This function returns the specific entropy
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::s(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: s() not implemented in the Solver object");
    return 0;
}

//! Compute total derivative of density ph
/*!
  This function returns the total derivative of density ph
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::d_der(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: d_der() not implemented in the Solver object");
    return 0;
}

//! Compute isentropic enthalpy
/*!
  This function returns the enthalpy at pressure p after an isentropic
  transformation from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param p New pressure
  @param properties ExternalThermodynamicState property struct corresponding to current state
*/
double BaseSolver::isentropicEnthalpy(double& p, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: isentropicEnthalpy() not implemented in the Solver object");
    return 0;
}

//! Set saturation properties from p
/*!
  This function sets the saturation properties for the given pressure p.
  The computed values are written to the ExternalSaturationProperties property struct.

  Must be re-implemented in the specific solver
  @param p Pressure
  @param properties ExternalSaturationProperties property struct
*/
void BaseSolver::setSat_p(double& p, ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setSat_p() not implemented in the Solver object");
}

//! Set saturation properties from T
/*!
  This function sets the saturation properties for the given temperature T.
  The computed values are written to the ExternalSaturationProperties property struct.

  Must be re-implemented in the specific solver
  @param T Temperature
  @param properties ExternalSaturationProperties property struct
*/
void BaseSolver::setSat_T(double& T, ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: setSat_T() not implemented in the Solver object");
}

//! Set bubble state
/*!
  This function sets the bubble state record bubbleProperties corresponding to the
  saturation data contained in the properties record.

  The default implementation of the setBubbleState function is relying on the correct
  behaviour of setState_ph with respect to the state input. Can be overridden
  in the specific solver code to get more efficient or correct handling of this situation.
  @param properties ExternalSaturationProperties record with saturation properties data
  @param phase Phase (1: one-phase, 2: two-phase)
  @param bubbleProperties ExternalThermodynamicState record where to write the bubble point properties
*/
void BaseSolver::setBubbleState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const bubbleProperties) {
    // Set the bubble state property record based on the saturation properties record
    setState_ph(properties->psat, properties->hl, phase, bubbleProperties);
}

//! Set dew state
/*!
  This function sets the dew state record dewProperties corresponding to the
  saturation data contained in the properties record.

  The default implementation of the setDewState function is relying on the correct
  behaviour of setState_ph with respect to the state input. Can be overridden
  in the specific solver code to get more efficient or correct handling of this situation.
  @param properties ExternalSaturationProperties record with saturation properties data
  @param phase Phase (1: one-phase, 2: two-phase)
  @param dewProperties ExternalThermodynamicState record where to write the dew point properties
*/
void BaseSolver::setDewState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const dewProperties) {
    // Set the dew state property record based on the saturation properties record
    setState_ph(properties->psat, properties->hv, phase, dewProperties);
}

//! Compute derivative of Ts wrt pressure
/*!
  This function returns the derivative of Ts wrt pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::dTp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dTp() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of dls wrt pressure
/*!
  This function returns the derivative of dls wrt pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::ddldp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddldp() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of dvs wrt pressure
/*!
  This function returns the derivative of dvs wrt pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::ddvdp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddvdp() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of hls wrt pressure
/*!
  This function returns the derivative of hls wrt pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::dhldp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dhldp() not implemented in the Solver object");
    return 0;
}

//! Compute derivative of hvs wrt pressure
/*!
  This function returns the derivative of hvs wrt pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::dhvdp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dhvdp() not implemented in the Solver object");
    return 0;
}

//! Compute density at bubble line
/*!
  This function returns the density at bubble line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::dl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dl() not implemented in the Solver object");
    return 0;
}

//! Compute density at dew line
/*!
  This function returns the density at dew line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::dv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dv() not implemented in the Solver object");
    return 0;
}

//! Compute enthalpy at bubble line
/*!
  This function returns the enthalpy at bubble line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::hl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: hl() not implemented in the Solver object");
    return 0;
}

//! Compute enthalpy at dew line
/*!
  This function returns the enthalpy at dew line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::hv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: hv() not implemented in the Solver object");
    return 0;
}

//! Compute surface tension
/*!
  This function returns the surface tension
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::sigma(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sigma() not implemented in the Solver object");
    return 0;
}

//! Compute entropy at bubble line
/*!
  This function returns the entropy at bubble line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::sl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sl() not implemented in the Solver object");
    return 0;
}

//! Compute entropy at dew line
/*!
  This function returns the entropy at dew line
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::sv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sv() not implemented in the Solver object");
    return 0;
}

//! Compute derivatives
/*!
  This function computes the derivatives according to the Bridgman's table.
  The computed values are written to the two phase medium property struct.
  This function can be called from within the setState_XX routines
  when implementing a new solver. Please be aware that cp, beta and
  kappa have to be provided to allow the computation of the derivatives. It
  returns false if the computation failed.

  Default implementation provided.
  @param properties ExternalThermodynamicState property record
*/
bool BaseSolver::computeDerivatives(ExternalThermodynamicState* const properties) {
    // Check whether cp is equal to zero
    if (properties->cp == 0.0) return false;
    // Check whether density is equal to zero
    if (properties->d == 0.0) return false;
    // Compute ddph
    properties->ddph =
      -(properties->T * properties->beta * properties->beta - properties->beta - properties->kappa * properties->d * properties->cp) / properties->cp;
    // Compute ddhp
    properties->ddhp = -properties->beta * properties->d / properties->cp;
    return true;
}

//! Compute saturation pressure
/*!
  This function returns the saturation pressure
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::psat(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: psat() not implemented in the Solver object");
    return 0;
}

//! Compute saturation temperature
/*!
  This function returns the saturation temperature
  from the state specified by the properties input

  Must be re-implemented in the specific solver
  @param properties ExternalSaturationProperties property struct corresponding to current state
*/
double BaseSolver::Tsat(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: Tsat() not implemented in the Solver object");
    return 0;
}
