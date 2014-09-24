#ifndef MIXTURE_PARAMETERS_H
#define MIXTURE_PARAMETERS_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{
    
/** \brief Get a comma-separated list of CAS code pairs separated by commas
 * 
 * Each of the pairs will be CAS1&CAS2 ("&" delimited)
 */
std::string get_csv_mixture_binary_pairs();

std::string get_mixing_pair_info(std::string CAS1, std::string CAS2, std::string key);

class MixtureParameters
{
public:
    static void set_mixture_parameters(HelmholtzEOSMixtureBackend &HEOS);
};

} /* namespace CoolProp */
#endif