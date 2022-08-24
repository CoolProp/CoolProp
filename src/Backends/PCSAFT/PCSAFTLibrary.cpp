#include <string>
#include <map>
#include "PCSAFTLibrary.h"
#include "all_pcsaft_JSON.h"                   // Makes a std::string variable called all_pcsaft_JSON
#include "pcsaft_fluids_schema_JSON.h"         // Makes a std::string variable called pcsaft_fluids_schema_JSON
#include "mixture_binary_pairs_pcsaft_JSON.h"  // Makes a std::string variable called mixture_binary_pairs_pcsaft_JSON
#include "rapidjson_include.h"
#include "CPstrings.h"
#include "CoolProp.h"
#include "Configuration.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "CoolPropTools.h"

namespace CoolProp {

std::string get_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key) {
    return PCSAFTLibrary::get_library().get_binary_interaction_pcsaft(CAS1, CAS2, key);
}

void set_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key, const double value) {
    PCSAFTLibrary::get_library().set_binary_interaction_pcsaft(CAS1, CAS2, key, value);
}

namespace PCSAFTLibrary {

static PCSAFTLibraryClass library;

PCSAFTLibraryClass& get_library(void) {
    return library;
}

PCSAFTLibraryClass::PCSAFTLibraryClass() : empty(true) {
    // This JSON formatted string comes from the all_pcsaft_JSON.h header which is a C++-escaped version of the JSON file
    add_fluids_as_JSON(all_pcsaft_JSON);

    // Then we add the library of binary interaction parameters
    if (m_binary_pair_map.size() == 0) {
        PCSAFTLibraryClass::load_from_string(mixture_binary_pairs_pcsaft_JSON);
    }
}

// Get a PCSAFTFluid instance stored in this library
PCSAFTFluid& PCSAFTLibraryClass::get(const std::string& key) {
    // Try to find it
    std::map<std::string, std::size_t>::iterator it = string_to_index_map.find(key);
    // If it is found
    if (it != string_to_index_map.end()) {
        return get(it->second);
    } else {
        throw ValueError(format("key [%s] was not found in string_to_index_map in PCSAFTLibraryClass", key.c_str()));
    }
}

/// Get a PCSAFTFluid instance stored in this library
/**
 @param key The index of the fluid in the map
 */
PCSAFTFluid& PCSAFTLibraryClass::get(std::size_t key) {
    // Try to find it
    std::map<std::size_t, PCSAFTFluid>::iterator it = fluid_map.find(key);
    // If it is found
    if (it != fluid_map.end()) {
        return it->second;
    } else {
        throw ValueError(format("key [%d] was not found in PCSAFTLibraryClass", key));
    }
};

void add_fluids_as_JSON(const std::string& JSON) {
    // First we validate the json string against the schema;
    std::string errstr;
    cpjson::schema_validation_code val_code = cpjson::validate_schema(pcsaft_fluids_schema_JSON, JSON, errstr);
    // Then we check the validation code

    if (val_code == cpjson::SCHEMA_VALIDATION_OK) {
        rapidjson::Document dd;

        dd.Parse<0>(JSON.c_str());
        if (dd.HasParseError()) {
            throw ValueError("Unable to load all_pcsaft_JSON.json");
        } else {
            try {
                library.add_many(dd);
            } catch (std::exception& e) {
                std::cout << e.what() << std::endl;
            }
        }
    } else {
        if (get_debug_level() > 0) {
            throw ValueError(format("Unable to load PC-SAFT library with error: %s", errstr.c_str()));
        }
    }
}

int PCSAFTLibraryClass::add_many(rapidjson::Value& listing) {
    int counter = 0;
    std::string fluid_name;
    for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr) {
        try {
            PCSAFTFluid fluid(itr);
            fluid_name = fluid.getName();

            // If the fluid is ok...

            // First check that none of the identifiers are already present
            bool already_present = false;

            if (string_to_index_map.find(fluid.getCAS()) != string_to_index_map.end()
                || string_to_index_map.find(fluid_name) != string_to_index_map.end()
                || string_to_index_map.find(upper(fluid_name)) != string_to_index_map.end()) {
                already_present = true;
            } else {
                // Check the aliases
                for (std::size_t i = 0; i < fluid.getAliases().size(); ++i) {
                    if (string_to_index_map.find(fluid.getAliases()[i]) != string_to_index_map.end()) {
                        already_present = true;
                        break;
                    }
                    if (string_to_index_map.find(upper(fluid.getAliases()[i])) != string_to_index_map.end()) {
                        already_present = true;
                        break;
                    }
                }
            }

            if (already_present) {
                if (!get_config_bool(OVERWRITE_FLUIDS)) {
                    throw ValueError(format(
                      "Cannot load fluid [%s:%s] because it is already in library; consider enabling the config boolean variable OVERWRITE_FLUIDS",
                      fluid.getName().c_str(), fluid.getCAS().c_str()));
                } else {
                    // Remove the one(s) that are already there

                    // Remove the actual fluid instance
                    std::size_t index = string_to_index_map.find(fluid_name)->second;

                    if (fluid_map.find(index) != fluid_map.end()) {
                        fluid_map.erase(fluid_map.find(index));
                    }

                    if (string_to_index_map.find(fluid_name) != string_to_index_map.end()) {
                        fluid_map.erase(fluid_map.find(index));
                    }

                    // Remove the identifiers pointing to that instance
                    if (string_to_index_map.find(fluid.getCAS()) != string_to_index_map.end()) {
                        string_to_index_map.erase(string_to_index_map.find(fluid.getCAS()));
                    }
                    if (string_to_index_map.find(fluid_name) != string_to_index_map.end()) {
                        string_to_index_map.erase(string_to_index_map.find(fluid_name));
                    }
                    // Check the aliases
                    for (std::size_t i = 0; i < fluid.getAliases().size(); ++i) {
                        if (string_to_index_map.find(fluid.getAliases()[i]) != string_to_index_map.end()) {
                            string_to_index_map.erase(string_to_index_map.find(fluid.getAliases()[i]));
                        }
                        if (string_to_index_map.find(upper(fluid.getAliases()[i])) != string_to_index_map.end()) {
                            string_to_index_map.erase(string_to_index_map.find(upper(fluid.getAliases()[i])));
                        }
                    }
                }
            }

            // By now, the library has been cleared of remnants of this fluid; safe to add the fluid now.

            // Get the next index for this fluid
            std::size_t index = fluid_map.size();

            // Add index->fluid mapping
            fluid_map[index] = fluid;

            // fluid_map[index] = cpjson::json2string(fluid_json);

            // Add CAS->index mapping
            string_to_index_map[fluid.getCAS()] = index;

            // Add name->index mapping
            string_to_index_map[fluid_name] = index;

            // Add the aliases
            for (std::size_t i = 0; i < fluid.getAliases().size(); ++i) {
                string_to_index_map[fluid.getAliases()[i]] = index;

                // Add uppercase alias for EES compatibility
                string_to_index_map[upper(fluid.getAliases()[i])] = index;
            }

            counter++;
            if (get_debug_level() > 5) {
                std::cout << format("Loaded.\n");
            }
        } catch (const std::exception& e) {
            throw ValueError(format("Unable to load fluid [%s] due to error: %s", fluid_name.c_str(), e.what()));
        }
    }
    return counter;
};

std::string get_pcsaft_fluids_schema() {
    return pcsaft_fluids_schema_JSON;
}

std::string PCSAFTLibraryClass::get_binary_interaction_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key) {
    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    std::vector<std::string> CASrev;
    CASrev.push_back(CAS2);
    CASrev.push_back(CAS1);

    if (m_binary_pair_map.find(CAS) != m_binary_pair_map.end()) {
        std::vector<Dictionary>& v = m_binary_pair_map[CAS];
        try {
            if (key == "name1") {
                return v[0].get_string("name1");
            } else if (key == "name2") {
                return v[0].get_string("name2");
            } else if (key == "BibTeX") {
                return v[0].get_string("BibTeX");
            } else if (key == "kij") {
                return format("%0.16g", v[0].get_double("kij"));
            } else if (key == "kijT") {
                try {
                    return format("%0.16g", v[0].get_double("kijT"));
                } catch (ValueError) {
                    return format("%0.16g", 0.0);
                }
            } else {
            }
        } catch (...) {
        }
        throw ValueError(format("Could not match the parameter [%s] for the binary pair [%s,%s] - for now this is an error.", key.c_str(),
                                CAS1.c_str(), CAS2.c_str()));
    } else if (m_binary_pair_map.find(CASrev) != m_binary_pair_map.end()) {
        std::vector<Dictionary>& v = m_binary_pair_map[CASrev];
        try {
            if (key == "name1") {
                return v[0].get_string("name1");
            } else if (key == "name2") {
                return v[0].get_string("name2");
            } else if (key == "BibTeX") {
                return v[0].get_string("BibTeX");
            } else if (key == "kij") {
                return format("%0.16g", v[0].get_double("kij"));
            } else if (key == "kijT") {
                try {
                    return format("%0.16g", v[0].get_double("kijT"));
                } catch (ValueError) {
                    return format("%0.16g", 0.0);
                }
            } else {
            }
        } catch (...) {
        }
        throw ValueError(format("Could not match the parameter [%s] for the binary pair [%s,%s] - for now this is an error.", key.c_str(),
                                CAS1.c_str(), CAS2.c_str()));
    } else {
        // Sort, see if other order works properly
        std::sort(CAS.begin(), CAS.end());
        if (m_binary_pair_map.find(CAS) != m_binary_pair_map.end()) {
            throw ValueError(format("Could not match the binary pair [%s,%s] - order of CAS numbers is backwards; found the swapped CAS numbers.",
                                    CAS1.c_str(), CAS2.c_str()));
        } else {
            throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS1.c_str(), CAS2.c_str()));
        }
    }
}

void PCSAFTLibraryClass::set_binary_interaction_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& key, const double value) {
    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    std::vector<std::string> CASrev;
    CASrev.push_back(CAS2);
    CASrev.push_back(CAS1);

    if (m_binary_pair_map.find(CAS) != m_binary_pair_map.end()) {
        if (get_config_bool(OVERWRITE_BINARY_INTERACTION)) {
            std::vector<Dictionary>& v = m_binary_pair_map[CAS];
            if (v[0].has_number(key)) {
                v[0].add_number(key, value);
            } else {
                throw ValueError(format("Could not set the parameter [%s] for the binary pair [%s,%s] - for now this is an error", key.c_str(),
                                        CAS1.c_str(), CAS2.c_str()));
            }
        } else {
            throw ValueError(
              format("CAS pair(%s,%s) already in binary interaction map; considering enabling configuration key OVERWRITE_BINARY_INTERACTION",
                     CAS1.c_str(), CAS2.c_str()));
        }
    } else if (m_binary_pair_map.find(CASrev) != m_binary_pair_map.end()) {
        if (get_config_bool(OVERWRITE_BINARY_INTERACTION)) {
            std::vector<Dictionary>& v = m_binary_pair_map[CASrev];
            if (v[0].has_number(key)) {
                v[0].add_number(key, value);
            } else {
                throw ValueError(format("Could not set the parameter [%s] for the binary pair [%s,%s] - for now this is an error", key.c_str(),
                                        CAS1.c_str(), CAS2.c_str()));
            }
        } else {
            throw ValueError(
              format("CAS pair(%s,%s) already in binary interaction map; considering enabling configuration key OVERWRITE_BINARY_INTERACTION",
                     CAS1.c_str(), CAS2.c_str()));
        }
    } else {
        Dictionary dict;
        std::vector<std::string> CAS;
        CAS.push_back(CAS1);
        CAS.push_back(CAS2);
        dict.add_number(key, value);

        m_binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary>>(CAS, std::vector<Dictionary>(1, dict)));
    }
}

void PCSAFTLibraryClass::load_from_JSON(rapidjson::Document& doc) {
    for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr) {
        // Get the empty dictionary to be filled by the appropriate interaction parameter
        Dictionary dict;

        // Get the vector of CAS numbers
        std::vector<std::string> CAS;
        CAS.push_back(cpjson::get_string(*itr, "CAS1"));
        CAS.push_back(cpjson::get_string(*itr, "CAS2"));
        std::string name1 = cpjson::get_string(*itr, "Name1");
        std::string name2 = cpjson::get_string(*itr, "Name2");

        // Sort the CAS number vector
        std::sort(CAS.begin(), CAS.end());

        // A sort was carried out, names/CAS were swapped
        bool swapped = CAS[0].compare(cpjson::get_string(*itr, "CAS1")) != 0;

        if (swapped) {
            std::swap(name1, name2);
        }

        // Populate the dictionary with common terms
        dict.add_string("name1", name1);
        dict.add_string("name2", name2);
        dict.add_string("BibTeX", cpjson::get_string(*itr, "BibTeX"));
        if (itr->HasMember("kij")) {
            dict.add_number("kij", cpjson::get_double(*itr, "kij"));
        } else {
            std::cout << "Loading error: binary pair of " << name1 << " & " << name2 << "does not provide kij" << std::endl;
        }
        if (itr->HasMember("kijT")) {
            dict.add_number("kijT", cpjson::get_double(*itr, "kijT"));
        }

        std::map<std::vector<std::string>, std::vector<Dictionary>>::iterator it = m_binary_pair_map.find(CAS);
        if (it == m_binary_pair_map.end()) {
            // Add to binary pair map by creating one-element vector
            m_binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary>>(CAS, std::vector<Dictionary>(1, dict)));
        } else {
            if (get_config_bool(OVERWRITE_BINARY_INTERACTION)) {
                // Already there, see http://www.cplusplus.com/reference/map/map/insert/, so we are going to pop it and overwrite it
                m_binary_pair_map.erase(it);
                std::pair<std::map<std::vector<std::string>, std::vector<Dictionary>>::iterator, bool> ret;
                ret = m_binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary>>(CAS, std::vector<Dictionary>(1, dict)));
                assert(ret.second == true);
            } else {
                // Error if already in map!
                throw ValueError(
                  format("CAS pair(%s,%s) already in binary interaction map; considering enabling configuration key OVERWRITE_BINARY_INTERACTION",
                         CAS[0].c_str(), CAS[1].c_str()));
            }
        }
    }
}

void PCSAFTLibraryClass::load_from_string(const std::string& str) {
    rapidjson::Document doc;
    doc.Parse<0>(str.c_str());
    if (doc.HasParseError()) {
        throw ValueError("Unable to parse PC-SAFT binary interaction parameter string");
    }
    load_from_JSON(doc);
}

}  // namespace PCSAFTLibrary
}  // namespace CoolProp
