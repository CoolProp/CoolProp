

#ifndef PCSAFT_LIBRARY_H
#define PCSAFT_LIBRARY_H

#include <vector>
#include <string>
#include "PCSAFTFluid.h"
#include "CoolPropTools.h"
#include "rapidjson_include.h"

namespace CoolProp {

std::string get_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key);
void set_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key, const double value);

namespace PCSAFTLibrary {

class PCSAFTLibraryClass
{
   private:
    std::map<std::size_t, PCSAFTFluid> fluid_map;
    std::map<std::string, std::size_t> string_to_index_map;
    bool empty;  // Is empty
    /// Map from sorted pair of CAS numbers to interaction parameter map.  The interaction parameter map is a map from key (string) to value (double)
    std::map<std::vector<std::string>, std::vector<Dictionary>> m_binary_pair_map;

    void load_from_JSON(rapidjson::Document& doc);
    void load_from_string(const std::string& str);

   public:
    PCSAFTLibraryClass();

    bool is_empty() {
        return empty;
    };

    int add_many(rapidjson::Value& listing);

    PCSAFTFluid& get(const std::string& key);
    PCSAFTFluid& get(std::size_t key);

    std::map<std::vector<std::string>, std::vector<Dictionary>>& binary_pair_map() {
        return m_binary_pair_map;
    };

    std::string get_binary_interaction_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key);
    void set_binary_interaction_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key, const double value);
};

/** \brief Add an array of fluids to the PC-SAFT library (as a JSON-formatted string)
 * @param JSON A JSON-formatted string with the fluid information
 */
void add_fluids_as_JSON(const std::string& JSON);

/// Get the schema used to validate the PC-SAFT fluids
std::string get_pcsaft_fluids_schema();

PCSAFTLibraryClass& get_library(void);
}  // namespace PCSAFTLibrary
} /* namespace CoolProp */

#endif
