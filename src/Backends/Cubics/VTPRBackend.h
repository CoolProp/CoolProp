//
//  VTPRBackend.h
//  CoolProp
//
//  Created by Ian on 7/17/16.
//
//

#ifndef VTPRBackend_h
#define VTPRBackend_h

#include "CoolPropTools.h"
#include "DataStructures.h"
#include "GeneralizedCubic.h"
#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Exceptions.h"
#include <vector>
#include "CubicBackend.h"
#include "Configuration.h"
#include "UNIFAQLibrary.h"
#include "UNIFAQ.h"
#include "VTPRCubic.h"

namespace CoolProp {
    
    
class VTPRBackend : public PengRobinsonBackend  {
    
private:
    
    std::vector<double> Tc, pc, omega, m_ii;
    double R;
    std::vector<std::string> m_fluid_names;
public:
    
    VTPRBackend(const std::vector<std::string> fluid_identifiers,
                const double R_u = get_config_double(R_U_CODATA),
                bool generate_SatL_and_SatV = true)
    {
        std::vector<double> Tc, pc, acentric;
        N = fluid_identifiers.size();
        components.resize(N);
        // Extract data from the UNIFAQ parameter library
        const UNIFAQLibrary::UNIFAQParameterLibrary & lib = LoadLibrary();
        for (std::size_t i = 0; i < fluid_identifiers.size(); ++i){
            UNIFAQLibrary::Component comp = lib.get_component("name", fluid_identifiers[i]);
            // Taking the mole fraction from the cubic library. TODO: Need a better solution
            components[i] = CubicLibrary::get_cubic_values(fluid_identifiers[i]);
            Tc.push_back(comp.Tc); // [K]
            pc.push_back(comp.pc); // [Pa]
            acentric.push_back(comp.acentric); // [-]
        }
        cubic.reset(new VTPRCubic(Tc, pc, acentric, R_u, lib));
        setup(fluid_identifiers, generate_SatL_and_SatV);
    };
    
    /// Set the alpha function based on the alpha function defined in the components vector;
    void set_alpha_from_components();
    
    /// Return the fluid names
    std::vector<std::string> calc_fluid_names(void) { return m_fluid_names; }
    
    /// Set the pointer to the residual helmholtz class, etc.
    void setup(const std::vector<std::string> &names, bool generate_SatL_and_SatV = true);
    
    /// Load the UNIFAQ library if needed and get const reference to it
    const UNIFAQLibrary::UNIFAQParameterLibrary &LoadLibrary();
    
    void set_mole_fractions(const std::vector<double> &z){
        mole_fractions = z;
        mole_fractions_double = z;
        VTPRCubic * _cubic= static_cast<VTPRCubic *>(cubic.get());
        _cubic->get_unifaq().set_mole_fractions(z);
    };
    
};
    
}; /* namespace CoolProp */

#endif /* VTPRBackend_h */
