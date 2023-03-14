#ifndef UNIFAC_H_
#define UNIFAC_H_

#include <map>

#include "UNIFACLibrary.h"
#include "CachedElement.h"
#include "Exceptions.h"

/// Structure containing data for the pure fluid in the mixture
struct ComponentData
{
    std::map<std::size_t, double> X, theta, lnGamma;
    int group_count;  ///< The total number of groups in the pure fluid
};

namespace UNIFAC {
class UNIFACMixture
{
   private:
    /// A const reference to the library of group and interaction parameters
    const UNIFACLibrary::UNIFACParameterLibrary& library;

    CoolProp::CachedElement _T;  ///< The cached temperature

    std::size_t N;  ///< Number of components

    double m_T;  ///< The temperature in K
    double T_r;  ///< Reducing temperature

    std::map<std::pair<std::size_t, std::size_t>, double> Psi_;  /// < temporary storage for Psi

    std::map<std::size_t, double> m_Xg,  ///< Map from sgi to mole fraction of group in the mixture
      m_thetag,                          ///< Map from sgi to theta for the group in the mixture
      m_lnGammag,                        ///< Map from sgi to ln(Gamma) for the group in the mixture
      m_Q;                               ///< Map from sgi to Q for the sgi

    /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
    std::map<std::pair<int, int>, UNIFACLibrary::InteractionParameters> interaction;

    /// A map from SGI to MGI
    std::map<std::size_t, std::size_t> m_sgi_to_mgi;

    /// The set of unique groups in this mixture
    std::set<std::size_t> unique_groups;

    std::vector<double> mole_fractions;

    std::vector<UNIFACLibrary::Component> components;

    std::vector<ComponentData> pure_data;

   public:
    UNIFACMixture(const UNIFACLibrary::UNIFACParameterLibrary& library, const double T_r) : library(library), T_r(T_r){};

    /**
        * \brief Set all the interaction parameters between groups
        *
        * \param subgroups A vector of the set of the unique Group forming the mixture - these
        * permutations represent the set of posisble binary interactions
        */
    void set_interaction_parameters();
    void set_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter, const double value);
    /// Get one of the mgi-mgi interaction pairs
    double get_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter);

    /// Set the mole fractions of the components in the mixtures (not the groups)
    void set_mole_fractions(const std::vector<double>& z);

    /// Get the mole fractions of the components in the mixtures (not the groups)
    const std::vector<double>& get_mole_fractions() {
        return mole_fractions;
    }

    /// Set the temperature of the components in the mixtures (not the groups)
    void set_temperature(const double T);

    /// Get the temperature
    double get_temperature() const {
        return m_T;
    }

    double Psi(std::size_t sgi1, std::size_t sgi2) const;

    double theta_pure(std::size_t i, std::size_t sgi) const;

    void activity_coefficients(double tau, const std::vector<double>& z, std::vector<double>& gamma);

    double ln_gamma_R(const double tau, std::size_t i, std::size_t itau);

    std::size_t group_count(std::size_t i, std::size_t sgi) const;

    /// Add a component with the defined groups defined by (count, sgi) pairs
    void add_component(const UNIFACLibrary::Component& comp);

    void set_components(const std::string& identifier_type, std::vector<std::string> identifiers);

    const std::vector<UNIFACLibrary::Component>& get_components() {
        return components;
    };

    void set_pure_data();

    /// Modify the surface parameter Q_k of the sub group sgi
    void set_Q_k(const size_t sgi, const double value);

    /// Get the surface parameter Q_k of the sub group sgi
    double get_Q_k(const size_t sgi) const;
};

} /* namespace UNIFAC */

#endif
