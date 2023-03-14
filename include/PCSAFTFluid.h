#ifndef PCSAFTFLUID_H
#define PCSAFTFLUID_H

#include <string>
#include <vector>
#include <map>

#include "rapidjson_include.h"

namespace CoolProp {

struct PCSAFTValues
{
    CoolPropDbl m;       ///< Number of segments
    CoolPropDbl sigma;   ///< Segment diameter (1/Angstrom)
    CoolPropDbl u;       ///< Dispersion energy divided by Boltzmann constant (K)
    CoolPropDbl uAB;     ///< Association energy (K)
    CoolPropDbl volA;    ///< Association volume
    std::vector<std::string> assocScheme; ///< The type of association for each associating functional group (see Huang and Radosz 1990)
    CoolPropDbl dipm;    ///< Dipole moment (Debye)
    CoolPropDbl dipnum;  ///< Number of dipole moments per molecule
    CoolPropDbl z;       ///< Charge of the compound
};

class PCSAFTFluid
{
   protected:
    std::string name;      // name of fluid
    std::string CAS;       // CAS number
    CoolPropDbl molemass;  ///< Molar mass (kg/mol)
    std::vector<std::string> aliases;
    PCSAFTValues params;

   public:
    PCSAFTFluid(){};
    PCSAFTFluid(rapidjson::Value::ValueIterator itr);
    ~PCSAFTFluid(){};

    std::string getName() const {
        return name;
    }
    std::string getCAS() const {
        return CAS;
    }
    CoolPropDbl molar_mass() const {
        return molemass;
    }
    std::vector<std::string> getAliases() const {
        return aliases;
    }
    CoolPropDbl getM() const {
        return params.m;
    }
    CoolPropDbl getSigma() const {
        return params.sigma;
    }
    CoolPropDbl getU() const {
        return params.u;
    }
    CoolPropDbl getUAB() const {
        return params.uAB;
    }
    CoolPropDbl getVolA() const {
        return params.volA;
    }
    std::vector<std::string> getAssocScheme() const {
        return params.assocScheme;
    }
    CoolPropDbl getDipm() const {
        return params.dipm;
    }
    CoolPropDbl getDipnum() const {
        return params.dipnum;
    }
    CoolPropDbl getZ() const {
        return params.z;
    }

    void calc_water_sigma(double t);
};

} /* namespace CoolProp */
#endif /* PCSAFTFLUID_H_ */
