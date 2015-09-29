#ifndef MIXTURE_PARAMETERS_H
#define MIXTURE_PARAMETERS_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{
    
/** \brief Get a comma-separated list of CAS code pairs
 * 
 * Each of the pairs will be CAS1&CAS2 ("&" delimited)
 */
std::string get_csv_mixture_binary_pairs();

/** \brief Get the parameters for a predefined mixture - R410A, R404A, etc. if the mixture is predefined
 * 
 */
bool is_predefined_mixture(const std::string &name, Dictionary &dict);

/** \brief Get a comma-separated list of predefined mixtures in CoolProp
 * 
 */
std::string get_csv_predefined_mixtures();
    
/** \brief Get a string for the given binary pair
 * 
 * 
 */
std::string get_mixture_binary_pair_data(const std::string &CAS1, const std::string &CAS2, const std::string &param);
    
/**
 * @brief Set a parameter for the given binary pair
 * @param CAS1 The CAS # for the first fluid (order matters!)
 * @param CAS2 The CAS # for the second fluid (order matters!)
 * @param param The parameter you want to set
 * @param val The value of the parameter
 * @return None
 */
void set_mixture_binary_pair_data(const std::string &CAS1, const std::string &CAS2, const std::string &param, const double val);

/**
 * @brief Apply a simple mixing rule for a given binary pair
 * @param CAS1 The CAS # for the first fluid (order matters!)
 * @param CAS2 The CAS # for the second fluid (order matters!)
 * @param rule The simple mixing rule to be used ("linear", "Lorentz-Berthelot")
 */
void apply_simple_mixing_rule(const std::string &CAS1, const std::string &CAS2, const std::string &rule);
    
class MixtureParameters
{
public:
    static void set_mixture_parameters(HelmholtzEOSMixtureBackend &HEOS);
};

} /* namespace CoolProp */
#endif