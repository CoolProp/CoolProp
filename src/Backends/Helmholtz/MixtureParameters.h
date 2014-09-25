#ifndef MIXTURE_PARAMETERS_H
#define MIXTURE_PARAMETERS_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{
    
/** \brief Get a comma-separated list of CAS code pairs
 * 
 * Each of the pairs will be CAS1&CAS2 ("&" delimited)
 */
std::string get_csv_mixture_binary_pairs();

/** \brief Get a string for the given binary pair
 * 
 * 
 */
std::string get_mixture_binary_pair_data(const std::string &CAS1, const std::string &CAS2, const std::string &param);

class MixtureParameters
{
public:
    static void set_mixture_parameters(HelmholtzEOSMixtureBackend &HEOS);
};

} /* namespace CoolProp */
#endif