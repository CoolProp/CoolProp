/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INCOMPRESSIBLEFLUID_H_
#define INCOMPRESSIBLEFLUID_H_

#include "DataStructures.h"
#include "Helmholtz.h"
#include "Solvers.h"

#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>

namespace CoolProp {

struct IncompressiblePolynomialData{
    std::vector<double> coeffs;
};

struct IncompressibleViscosityVariables{
 enum IncompressibleViscosityEnum{
     INCOMPRESSIBLE_VISCOSITY_POLYNOMIAL,
     INCOMPRESSIBLE_VISCOSITY_NOT_SET
 };
 int type;
 IncompressiblePolynomialData poly;
 IncompressibleViscosityVariables(){type = INCOMPRESSIBLE_VISCOSITY_NOT_SET;};
};

struct IncompressibleConductivityVariables{
 enum IncompressibleConductivityEnum{
     INCOMPRESSIBLE_CONDUCTIVITY_POLYNOMIAL,
     INCOMPRESSIBLE_CONDUCTIVITY_NOT_SET
 };
 int type;
 IncompressiblePolynomialData poly;
 IncompressibleConductivityVariables(){type = INCOMPRESSIBLE_CONDUCTIVITY_NOT_SET;};
};

struct IncompressibleDensityVariables{
 enum IncompressibleDensityEnum{
     INCOMPRESSIBLE_DENSITY_POLYNOMIAL,
     INCOMPRESSIBLE_DENSITY_NOT_SET
 };
 int type;
 IncompressiblePolynomialData poly;
 IncompressibleDensityVariables(){type = INCOMPRESSIBLE_DENSITY_NOT_SET;};
};

struct IncompressibleSpecificHeatVariables{
 enum IncompressibleSpecificHeatEnum{
     INCOMPRESSIBLE_SPECIFIC_HEAT_POLYNOMIAL,
     INCOMPRESSIBLE_SPECIFIC_HEAT_NOT_SET
 };
 int type;
 IncompressiblePolynomialData poly;
 IncompressibleSpecificHeatVariables(){type = INCOMPRESSIBLE_SPECIFIC_HEAT_NOT_SET;};
};

/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
struct IncompressibleFluid {

    std::string name;
    std::string reference;

    double Tmin, Tmax;

    IncompressibleViscosityVariables viscosity;
    IncompressibleConductivityVariables conductivity;
    IncompressibleSpecificHeatVariables specific_heat;
    IncompressibleDensityVariables density;

};


} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
