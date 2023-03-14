

#ifndef HUMAIR_H
#define HUMAIR_H

#include "CoolPropTools.h"

namespace HumidAir {
/* \brief Standard I/O function using base SI units exclusively
 *
 */
double HAPropsSI(const std::string& OutputName, const std::string& Input1Name, double Input1, const std::string& Input2Name, double Input2,
                 const std::string& Input3Name, double Input3);

/* \brief Standard I/O function using mixed kSI units
 *
 * \warning DEPRECATED!! Use \ref HAPropsSI
 */
double HAProps(const std::string& OutputName, const std::string& Input1Name, double Input1, const std::string& Input2Name, double Input2,
               const std::string& Input3Name, double Input3);

// -----------------------
// Extra I/O function
// -----------------------
double HAProps_Aux(const char* OutputName, double T, double p, double W, char* units);

// Properties for Ice Ih at temperatures below 273.16 K
double IceProps(const char* Name, double T, double p);

//Turn on the use of virial correlations for air and water
void UseVirialCorrelations(int flag);
void UseIsothermCompressCorrelation(int flag);
void UseIdealGasEnthalpyCorrelations(int flag);

// --------------
// Help functions
// --------------
void HAHelp(void);
int returnHumAirCode(const char* Code);

// ----------------------
// Other simple functions
// ----------------------
double cair_sat(double T);

} /* namespace HumidAir */
#endif
