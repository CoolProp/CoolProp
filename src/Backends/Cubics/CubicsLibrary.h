

#ifndef CUBICS_LIBRARY_H
#define CUBICS_LIBRARY_H

#include <vector>
#include <string>

namespace CoolProp{

struct CubicsValues{
    double Tc, ///< Critical temperature (K)
        pc, ///< Critical pressure (Pa)
        molemass, ///< Molar mass (kg/mol)
        acentric; ///< Acentric factor (-)
    std::string name, // name of fluid
                CAS, // CAS reference number of fluid
                BibTeX; // BibTex key(s) for the values
    std::vector<std::string> aliases;
};

/**
 * @param identifier The name or registry number of the fluid (or an alias)
 */
CubicsValues get_cubic_values(const std::string &identifier);

} /* namespace CoolProp */

#endif