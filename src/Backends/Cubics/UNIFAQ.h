#ifndef UNIFAQ_H_
#define UNIFAQ_H_

#include <map>

#include "UNIFAQLibrary.h"
#include "CachedElement.h"

/// Structure containing data for the pure fluid in the mixture
struct ComponentData {
    std::map<std::size_t, double> X, theta, lnGamma;
    int group_count; ///< The total number of groups in the pure fluid
};

namespace UNIFAQ
{
    class UNIFAQMixture
    {
    private:
        CoolProp::CachedElement _T; ///< The cached temperature

        double m_T; ///< The temperature in K

        std::vector<double> m_r,
                            m_q,
                            m_l,
                            m_phi,
                            m_theta,
                            m_ln_Gamma_C;

        std::map<std::size_t, double> m_Xg,  ///< Map from sgi to mole fraction of group in the mixture
                                      m_thetag, ///< Map from sgi to theta for the group in the mixture
                                      m_lnGammag; ///< Map from sgi to ln(Gamma) for the group in the mixture

        /// A const reference to the library of group and interaction parameters
        const UNIFAQLibrary::UNIFAQParameterLibrary &library;

        /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
        std::map<std::pair<int, int>, UNIFAQLibrary::InteractionParameters> interaction;

        /// A map from SGI to MGI
        std::map<std::size_t, std::size_t> m_sgi_to_mgi;

        /// A vector of unique groups in this mixture
        std::vector<UNIFAQLibrary::Group> unique_groups;
    
        std::vector<double> mole_fractions;

        std::vector<UNIFAQLibrary::Component> components;

        std::vector<ComponentData> pure_data;
    
    public:
        
        UNIFAQMixture(const UNIFAQLibrary::UNIFAQParameterLibrary &library) : library(library) {};

        /** 
        * \brief Set all the interaction parameters between groups
        *
        * \param subgroups A vector of the set of the unique Group forming the mixture - these 
        * permutations represent the set of posisble binary interactions
        */
        void set_interaction_parameters();

        /// Set the mole fractions of the components in the mixtures (not the groups)
        void set_mole_fractions(const std::vector<double> &z);
        
        /// Get the mole fractions of the components in the mixtures (not the groups)
        const std::vector<double> & get_mole_fractions() { return mole_fractions; }

        /// Set the mole fractions of the components in the mixtures (not the groups) AND the mole fractions you want to use
        void set_temperature(const double T, const std::vector<double> &z);

        /// Get the temperature
        double get_temperature() const { return m_T; }

        double Psi(std::size_t sgi1, std::size_t sgi2) const;

        double theta_pure(std::size_t i, std::size_t sgi) const;

        double activity_coefficient(std::size_t i) const;

        double ln_gamma_R(std::size_t i) const;

        std::size_t group_count(std::size_t i, std::size_t sgi) const;

        /// Add a component with the defined groups defined by (count, sgi) pairs
        void add_component(const UNIFAQLibrary::Component &comp);
    
        void set_components(const std::string &identifier_type, std::vector<std::string> identifiers);
        
        const std::vector<UNIFAQLibrary::Component> & get_components() { return components; };
    };

} /* namespace UNIFAQ */

#endif