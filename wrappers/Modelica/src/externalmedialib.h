/*!
  \file externalmedialib.h
  \brief Header file to be included in the Modelica tool, with external function interfaces

  C/C++ layer for external medium models extending from
  PartialExternalTwoPhaseMedium.

  Francesco Casella, Christoph Richter, Roberto Bonifetto
  2006-2012
  Copyright Politecnico di Milano, TU Braunschweig, Politecnico di Torino
*/

#ifndef EXTERNALMEDIALIB_H_
#define EXTERNALMEDIALIB_H_

// Constants for input choices (see ExternalMedia.Common.InputChoices)
#define CHOICE_dT 1
#define CHOICE_ph 2
#define CHOICE_ps 3
#define CHOICE_pT 4

// Define struct
//! ExternalThermodynamicState property struct
/*!
  The ExternalThermodynamicState property struct defines all the properties that
  are computed by external Modelica medium models extending from
  PartialExternalTwoPhaseMedium.
*/

typedef struct
{

    //! Temperature
    double T;
    //! Velocity of sound
    double a;
    //! Isobaric expansion coefficient
    double beta;
    //! Specific heat capacity cp
    double cp;
    //! Specific heat capacity cv
    double cv;
    //! Density
    double d;
    //! Derivative of density wrt enthalpy at constant pressure
    double ddhp;
    //! Derivative of density wrt pressure at constant enthalpy
    double ddph;
    //! Dynamic viscosity
    double eta;
    //! Specific enthalpy
    double h;
    //! Compressibility
    double kappa;
    //! Thermal conductivity
    double lambda;
    //! Pressure
    double p;
    //! Phase flag: 2 for two-phase, 1 for one-phase
    int phase;
    //! Specific entropy
    double s;

} ExternalThermodynamicState;

//! ExternalSaturationProperties property struct
/*!
  The ExternalSaturationProperties property struct defines all the saturation properties
  for the dew and the bubble line that are computed by external Modelica medium models
  extending from PartialExternalTwoPhaseMedium.
*/

typedef struct
{
    //! Saturation temperature
    double Tsat;
    //! Derivative of Ts wrt pressure
    double dTp;
    //! Derivative of dls wrt pressure
    double ddldp;
    //! Derivative of dvs wrt pressure
    double ddvdp;
    //! Derivative of hls wrt pressure
    double dhldp;
    //! Derivative of hvs wrt pressure
    double dhvdp;
    //! Density at bubble line (for pressure ps)
    double dl;
    //! Density at dew line (for pressure ps)
    double dv;
    //! Specific enthalpy at bubble line (for pressure ps)
    double hl;
    //! Specific enthalpy at dew line (for pressure ps)
    double hv;
    //! Saturation pressure
    double psat;
    //! Surface tension
    double sigma;
    //! Specific entropy at bubble line (for pressure ps)
    double sl;
    //! Specific entropy at dew line (for pressure ps)
    double sv;

} ExternalSaturationProperties;

// Define export
#ifdef __cplusplus
#    define EXPORT __declspec(dllexport)
#else
#    define EXPORT
#endif  // __cplusplus

#ifdef __cplusplus
extern "C"
{
#endif  // __cplusplus

    EXPORT double TwoPhaseMedium_getMolarMass_C_impl(const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_getCriticalTemperature_C_impl(const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_getCriticalPressure_C_impl(const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_getCriticalMolarVolume_C_impl(const char* mediumName, const char* libraryName, const char* substanceName);

    EXPORT void TwoPhaseMedium_setState_ph_C_impl(double p, double h, int phase, ExternalThermodynamicState* state, const char* mediumName,
                                                  const char* libraryName, const char* substanceName);
    EXPORT void TwoPhaseMedium_setState_pT_C_impl(double p, double T, ExternalThermodynamicState* state, const char* mediumName,
                                                  const char* libraryName, const char* substanceName);
    EXPORT void TwoPhaseMedium_setState_dT_C_impl(double d, double T, int phase, ExternalThermodynamicState* state, const char* mediumName,
                                                  const char* libraryName, const char* substanceName);
    EXPORT void TwoPhaseMedium_setState_ps_C_impl(double p, double s, int phase, ExternalThermodynamicState* state, const char* mediumName,
                                                  const char* libraryName, const char* substanceName);

    EXPORT double TwoPhaseMedium_prandtlNumber_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                      const char* substanceName);
    EXPORT double TwoPhaseMedium_temperature_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                    const char* substanceName);
    EXPORT double TwoPhaseMedium_velocityOfSound_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                        const char* substanceName);
    EXPORT double TwoPhaseMedium_isobaricExpansionCoefficient_C_impl(ExternalThermodynamicState* state, const char* mediumName,
                                                                     const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_specificHeatCapacityCp_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                               const char* substanceName);
    EXPORT double TwoPhaseMedium_specificHeatCapacityCv_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                               const char* substanceName);
    EXPORT double TwoPhaseMedium_density_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                const char* substanceName);
    EXPORT double TwoPhaseMedium_density_derh_p_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                       const char* substanceName);
    EXPORT double TwoPhaseMedium_density_derp_h_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                       const char* substanceName);
    EXPORT double TwoPhaseMedium_dynamicViscosity_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                         const char* substanceName);
    EXPORT double TwoPhaseMedium_specificEnthalpy_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                         const char* substanceName);
    EXPORT double TwoPhaseMedium_isothermalCompressibility_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                                  const char* substanceName);
    EXPORT double TwoPhaseMedium_thermalConductivity_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                            const char* substanceName);
    EXPORT double TwoPhaseMedium_pressure_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                 const char* substanceName);
    EXPORT double TwoPhaseMedium_specificEntropy_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                        const char* substanceName);
    EXPORT double TwoPhaseMedium_density_ph_der_C_impl(ExternalThermodynamicState* state, const char* mediumName, const char* libraryName,
                                                       const char* substanceName);
    EXPORT double TwoPhaseMedium_isentropicEnthalpy_C_impl(double p_downstream, ExternalThermodynamicState* refState, const char* mediumName,
                                                           const char* libraryName, const char* substanceName);

    EXPORT void TwoPhaseMedium_setSat_p_C_impl(double p, ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                               const char* substanceName);
    EXPORT void TwoPhaseMedium_setSat_T_C_impl(double T, ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                               const char* substanceName);
    EXPORT void TwoPhaseMedium_setBubbleState_C_impl(ExternalSaturationProperties* sat, int phase, ExternalThermodynamicState* state,
                                                     const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT void TwoPhaseMedium_setDewState_C_impl(ExternalSaturationProperties* sat, int phase, ExternalThermodynamicState* state,
                                                  const char* mediumName, const char* libraryName, const char* substanceName);

    EXPORT double TwoPhaseMedium_saturationTemperature_C_impl(double p, const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_saturationTemperature_derp_C_impl(double p, const char* mediumName, const char* libraryName,
                                                                   const char* substanceName);
    EXPORT double TwoPhaseMedium_saturationTemperature_derp_sat_C_impl(ExternalSaturationProperties* sat, const char* mediumName,
                                                                       const char* libraryName, const char* substanceName);

    EXPORT double TwoPhaseMedium_dBubbleDensity_dPressure_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                                 const char* substanceName);
    EXPORT double TwoPhaseMedium_dDewDensity_dPressure_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                              const char* substanceName);
    EXPORT double TwoPhaseMedium_dBubbleEnthalpy_dPressure_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                                  const char* substanceName);
    EXPORT double TwoPhaseMedium_dDewEnthalpy_dPressure_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                               const char* substanceName);
    EXPORT double TwoPhaseMedium_bubbleDensity_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                      const char* substanceName);
    EXPORT double TwoPhaseMedium_dewDensity_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                   const char* substanceName);
    EXPORT double TwoPhaseMedium_bubbleEnthalpy_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                       const char* substanceName);
    EXPORT double TwoPhaseMedium_dewEnthalpy_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                    const char* substanceName);
    EXPORT double TwoPhaseMedium_saturationPressure_C_impl(double T, const char* mediumName, const char* libraryName, const char* substanceName);
    EXPORT double TwoPhaseMedium_surfaceTension_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                       const char* substanceName);
    EXPORT double TwoPhaseMedium_bubbleEntropy_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                      const char* substanceName);
    EXPORT double TwoPhaseMedium_dewEntropy_C_impl(ExternalSaturationProperties* sat, const char* mediumName, const char* libraryName,
                                                   const char* substanceName);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif /*EXTERNALMEDIALIB_H_*/
