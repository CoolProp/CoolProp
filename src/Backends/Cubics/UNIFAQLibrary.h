#ifndef UNIFAQ_LIBRARY_H
#define UNIFAQ_LIBRARY_H

#include <vector>
#include <exception>

#include "rapidjson_include.h"

namespace UNIFAQLibrary{

    /// A structure containing references for a single group (its multiplicity, main group index, etc.)
    struct Group{
        int sgi, ///< Sub group index
            mgi; ///< Main group index
        double R_k, ///< R_k
            Q_k; ///< Q_k
    };

    /// A structure containing the parameters for a given mgi-mgi pair
    struct InteractionParameters{
        int mgi1, ///< The first main group index
            mgi2; ///< The second main group index
        double a_ij, ///< 
               a_ji, ///< 
               b_ij, ///< 
               b_ji, ///< 
               c_ij, ///< 
               c_ji; ///< 
        /// Swap a_ij with a_ji, b_ij with b_ji, etc.
        void swap() {
            std::swap(a_ij, a_ji);
            std::swap(b_ij, b_ji);
            std::swap(c_ij, c_ji);
        }
        /// Set all the values to 0
        void zero_out() {
            a_ij = 0; a_ji = 0;
            b_ij = 0; b_ji = 0;
            c_ij = 0; c_ji = 0;
        }
    };

    /// A structure containing a group (its count, index, etc.) for a subgroup forming a part of a component
    struct ComponentGroup {
        int count;
        UNIFAQLibrary::Group group;
        ComponentGroup(const int count, const UNIFAQLibrary::Group group) : count(count), group(group) {};
    };

    /// A structure containing the groups and additional information for a component
    struct Component{
        std::string name,            ///< A user-readable name (not guaranteed unique)
                    inchikey,        ///< The InChI key for the component
                    registry_number, ///< The registry number for the component in xxxxxxxxxx-xx-x format
                    userid;          ///< A user-specified string identifier
        double Tc,                   ///< The critical temperature in K
               pc,                   ///< The critical pressure in Pa
               acentric;             ///< The acentric factor
        std::vector<ComponentGroup> groups;
    };

    /**
     * \brief A container for the parameters for a given UNIFAQ model
     *
     * This container is intended to be sufficiently generic to allow the user to populate it with UNIFAQ parameters from
     * any of the standard UNIFAQ models
     * 
     * Input of parameters (population) is done using JSON-formatted strings, and the class can be interrogated to return
     * the desired group information and/or interaction parameters
     */
    struct UNIFAQParameterLibrary{
    private:
        bool m_populated; ///< True if the library has been populated
        std::vector<Group> groups; ///< The collection of groups forming the component from the group decomposition
        std::vector<InteractionParameters> interaction_parameters; ///< The collection of interaction parameters between main groups in the library
        std::vector<Component> components; ///< The collection of components that are included in this library

        /// Convert string to JSON document
        void jsonize(std::string &s, rapidjson::Document &doc);

        /// Populate internal data structures based on rapidjson Documents
        void populate(rapidjson::Value &group_data, rapidjson::Value &interaction_data, rapidjson::Value &decomp_data);

    public:
        UNIFAQParameterLibrary() : m_populated(false) {};
        
        /// Return true if library has been populated
        bool is_populated(){ return m_populated; };
        
        /// Populate internal data structures based on JSON-formatted strings
        void populate(std::string &group_data, std::string &interaction_data, std::string &decomp_data);
        
        /// Get the data for group with given sub group index
        Group get_group(int sgi) const;

        /// Check if the sub group index can be retrieved
        bool has_group(int sgi) const;

        /// Get the group decomposition for a given component
        Component get_component(const std::string &identifier, const std::string &value) const;
        
        /// Get the interaction parameters for given mgi-mgi pair
        InteractionParameters get_interaction_parameters(int mgi1, int mgi2) const;
    };

}; /* namespace UNIFAQLibrary*/

#endif