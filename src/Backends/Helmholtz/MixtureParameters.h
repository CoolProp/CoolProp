#ifndef MIXTURE_PARAMETERS_H
#define MIXTURE_PARAMETERS_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

/** \brief Get a comma-separated list of CAS code pairs
 *
 * Each of the pairs will be CAS1&CAS2 ("&" delimited)
 */
std::string get_csv_mixture_binary_pairs();

/** \brief Get the parameters for a predefined mixture - R410A, R404A, etc. if the mixture is predefined
 *
 */
bool is_predefined_mixture(const std::string& name, Dictionary& dict);

/** \brief Get a comma-separated list of predefined mixtures in CoolProp
 *
 */
std::string get_csv_predefined_mixtures();

/** \brief Get a string for the given binary pair
 *
 *
 */
std::string get_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& param);

/**
 * @brief Set a parameter for the given binary pair
 * @param CAS1 The CAS # for the first fluid (order matters!)
 * @param CAS2 The CAS # for the second fluid (order matters!)
 * @param param The parameter you want to set
 * @param val The value of the parameter
 * @return None
 */
void set_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& param, const double val);

/**
 * @brief Apply a simple mixing rule for a given binary pair
 * @param identifier1 The CAS # (or name) for the first fluid
 * @param identifier2 The CAS # (or name) for the second fluid
 * @param rule The simple mixing rule to be used ("linear", "Lorentz-Berthelot")
 */
void apply_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule);

class MixtureParameters
{
   public:
    static void set_mixture_parameters(HelmholtzEOSMixtureBackend& HEOS);
};

/**
 * @brief Get the allocated Departure function for a given departure function name
 * @param Name The name of the function to be used, or its alias
 * @warning The pointer points to an instance created with new, you should manage the pointer with shared_ptr or similar
 */
DepartureFunction* get_departure_function(const std::string& Name);

/// A Data structure for holding BIP coming from REFPROP
struct REFPROP_binary_element
{
    std::string CAS1, CAS2, model;
    double betaT, gammaT, betaV, gammaV, Fij;
    std::vector<std::string> comments;
};
/// A data structure for holding departure functions coming from REFPROP
struct REFPROP_departure_function
{
    short Npower, Nspecial, Nterms_power, Nterms_special;
    std::string model;
    std::vector<double> a, t, d, e, eta, epsilon, beta, gamma;
    std::vector<std::string> comments;
};

/**
 * @brief Set the departure functions in the departure function library from a string format
 * @param string_data The departure functions to be set, either provided as a JSON-formatted string
 *                    or as a string of the contents of a HMX.BNC file from REFPROP
 *
 * @note By default, if a departure function already exists in the library, this is an error,
 * unless the configuration variable OVERWRITE_DEPARTURE_FUNCTIONS is set to true
 */
void set_departure_functions(const std::string& string_data);

/**
 * @brief Set the interaction parameters from a string format
 * @param string_data The model parameters, as a JSON-formatted string
 *
 */
void set_interaction_parameters(const std::string& string_data);

} /* namespace CoolProp */
#endif
