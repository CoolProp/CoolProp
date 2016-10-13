

#ifndef CUBICS_LIBRARY_H
#define CUBICS_LIBRARY_H

#include <vector>
#include <string>

namespace CoolProp{
    
    namespace CubicLibrary{

        struct CubicsValues{
            double Tc, ///< Critical temperature (K)
                pc, ///< Critical pressure (Pa)
                molemass, ///< Molar mass (kg/mol)
                acentric; ///< Acentric factor (-)
            std::string name, // name of fluid
                        CAS, // CAS reference number of fluid
                        BibTeX; // BibTex key(s) for the values
            std::vector<std::string> aliases;
            std::string alpha_type; ///< The type of alpha function
            std::vector<double> alpha_coeffs; ///< The vector of coefficients for the alpha function
        };

        /**
         * @param identifier The name or registry number of the fluid (or an alias)
         */
        CubicsValues get_cubic_values(const std::string &identifier);

            
        /** \brief Add an array of fluids to the cubics library (as a JSON-formatted string)
         * @param JSON A JSON-formatted string with the fluid information
         */
        void add_fluids_as_JSON(const std::string &JSON);
        
        /// Get the schema used to validate the cubic fluids
        std::string get_cubic_fluids_schema();
        
    } /* namespace CubicLibrary */

} /* namespace CoolProp */

#endif