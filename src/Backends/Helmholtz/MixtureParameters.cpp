#include "MixtureParameters.h"
#include "CPstrings.h"
#include "mixture_departure_functions_JSON.h"  // Creates the variable mixture_departure_functions_JSON
#include "mixture_binary_pairs_JSON.h"         // Creates the variable mixture_binary_pairs_JSON
#include "predefined_mixtures_JSON.h"          // Makes a std::string variable called predefined_mixtures_JSON

namespace CoolProp {

/** \brief A library of predefined mixtures
 *
 * Each entry in the predefined mixture library contains the names and mole fractions for the binary pairs
 */
class PredefinedMixturesLibrary
{
   public:
    std::map<std::string, Dictionary> predefined_mixture_map;

    PredefinedMixturesLibrary() {
        rapidjson::Document doc;

        doc.Parse<0>(predefined_mixtures_JSON.c_str());
        if (doc.HasParseError()) {
            throw ValueError();
        }

        // Iterate over the papers in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr) {
            // Instantiate the empty dictionary to be filled
            Dictionary dict;
            // Get the name
            std::string name = cpjson::get_string(*itr, "name") + ".mix";
            // Get the fluid names
            dict.add_string_vector("fluids", cpjson::get_string_array(*itr, "fluids"));
            // Get the mole fractions
            dict.add_double_vector("mole_fractions", cpjson::get_double_array(*itr, "mole_fractions"));
            // Add to the map
            predefined_mixture_map.insert(std::pair<std::string, Dictionary>(name, dict));
            // Also add the uppercase version to the map
            predefined_mixture_map.insert(std::pair<std::string, Dictionary>(upper(name), dict));
        }
    }
};
static PredefinedMixturesLibrary predefined_mixtures_library;

std::string get_csv_predefined_mixtures() {
    std::vector<std::string> out;
    for (std::map<std::string, Dictionary>::const_iterator it = predefined_mixtures_library.predefined_mixture_map.begin();
         it != predefined_mixtures_library.predefined_mixture_map.end(); ++it) {
        out.push_back(it->first);
    }
    return strjoin(out, ",");
}

bool is_predefined_mixture(const std::string& name, Dictionary& dict) {
    std::map<std::string, Dictionary>::const_iterator iter = predefined_mixtures_library.predefined_mixture_map.find(name);
    if (iter != predefined_mixtures_library.predefined_mixture_map.end()) {
        dict = iter->second;
        return true;
    } else {
        return false;
    }
}

/** \brief A library of binary pair parameters for the mixture
 *
 * Each entry in the binary pair library includes reducing parameters as well as the name of the reducing function to be used and
 */
class MixtureBinaryPairLibrary
{
   private:
    /// Map from sorted pair of CAS numbers to reducing parameter map.  The reducing parameter map is a map from key (string) to value (double)
    std::map<std::vector<std::string>, std::vector<Dictionary>> m_binary_pair_map;

   public:
    std::map<std::vector<std::string>, std::vector<Dictionary>>& binary_pair_map() {
        // Set the default departure functions if none have been provided yet
        if (m_binary_pair_map.size() == 0) {
            load_defaults();
        }
        return m_binary_pair_map;
    };

    void load_from_string(const std::string& str) {
        rapidjson::Document doc;
        doc.Parse<0>(str.c_str());
        if (doc.HasParseError()) {
            std::cout << str << std::endl;
            throw ValueError("Unable to parse binary interaction function string");
        }
        load_from_JSON(doc);
    }

    // Load the defaults that come from the JSON-encoded string compiled into library
    // as the variable mixture_departure_functions_JSON
    void load_defaults() {
        load_from_string(mixture_binary_pairs_JSON);
    }

    /** \brief Construct the binary pair library including all the binary pairs that are possible
     *
     * The data structure also includes space for a string that gives the pointer to the departure function to be used for this binary pair.
     */
    void load_from_JSON(rapidjson::Document& doc) {

        // Iterate over the papers in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr) {
            // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
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
            dict.add_number("F", cpjson::get_double(*itr, "F"));
            if (std::abs(dict.get_number("F")) > DBL_EPSILON) {
                dict.add_string("function", cpjson::get_string(*itr, "function"));
            }

            if (itr->HasMember("xi") && itr->HasMember("zeta")) {
                dict.add_string("type", "Lemmon-xi-zeta");
                // Air and HFC mixtures from Lemmon - we could also directly do the conversion
                dict.add_number("xi", cpjson::get_double(*itr, "xi"));
                dict.add_number("zeta", cpjson::get_double(*itr, "zeta"));
            } else if (itr->HasMember("gammaT") && itr->HasMember("gammaV") && itr->HasMember("betaT") && itr->HasMember("betaV")) {
                dict.add_string("type", "GERG-2008");
                dict.add_number("gammaV", cpjson::get_double(*itr, "gammaV"));
                dict.add_number("gammaT", cpjson::get_double(*itr, "gammaT"));

                double betaV = cpjson::get_double(*itr, "betaV");
                double betaT = cpjson::get_double(*itr, "betaT");
                if (swapped) {
                    dict.add_number("betaV", 1 / betaV);
                    dict.add_number("betaT", 1 / betaT);
                } else {
                    dict.add_number("betaV", betaV);
                    dict.add_number("betaT", betaT);
                }
            } else {
                std::cout << "Loading error: binary pair of " << name1 << " & " << name2
                          << "does not provide either a) xi and zeta b) gammaT, gammaV, betaT, and betaV" << std::endl;
                continue;
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
                    ret =
                      m_binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary>>(CAS, std::vector<Dictionary>(1, dict)));
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
    /// Add a simple mixing rule
    void add_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule) {
        // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
        Dictionary dict;

        // Get the names/CAS of the compounds
        std::string CAS1, CAS2, name1 = identifier1, name2 = identifier2;
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS1, HEOS2;

        std::vector<std::string> id1split = strsplit(identifier1, '-');
        if (id1split.size() == 3) {  // Check if identifier is in CAS format
            CAS1 = identifier1;
        } else {
            std::vector<std::string> names1(1, identifier1);
            HEOS1.reset(new CoolProp::HelmholtzEOSMixtureBackend(names1));
            CAS1 = HEOS1->fluid_param_string("CAS");
        }

        std::vector<std::string> id2split = strsplit(identifier2, '-');
        if (id2split.size() == 3) {  // Check if identifier is in CAS format
            CAS2 = identifier2;
        } else {
            std::vector<std::string> names2(1, identifier2);
            HEOS2.reset(new CoolProp::HelmholtzEOSMixtureBackend(names2));
            CAS2 = HEOS2->fluid_param_string("CAS");
        }

        // Get the vector of CAS numbers
        std::vector<std::string> CAS;
        CAS.push_back(CAS1);
        CAS.push_back(CAS2);

        // Sort the CAS number vector
        std::sort(CAS.begin(), CAS.end());

        // Swap fluid names if the CAS numbers were swapped
        if (CAS[0] != CAS1) {
            std::swap(name1, name2);
        }

        // Populate the dictionary with common terms
        dict.add_string("name1", name1);
        dict.add_string("name2", name2);
        dict.add_string("BibTeX", "N/A - " + rule);
        dict.add_number("F", 0);
        dict.add_string("type", "GERG-2008");

        if (rule == "linear") {
            // Terms for linear mixing
            HEOS1.reset(new CoolProp::HelmholtzEOSMixtureBackend(std::vector<std::string>(1, name1)));
            HEOS2.reset(new CoolProp::HelmholtzEOSMixtureBackend(std::vector<std::string>(1, name2)));

            dict.add_number("gammaT", 0.5 * (HEOS1->T_critical() + HEOS2->T_critical()) / sqrt(HEOS1->T_critical() * HEOS2->T_critical()));
            double rhoc1 = HEOS1->rhomolar_critical(), rhoc2 = HEOS2->rhomolar_critical();
            dict.add_number("gammaV", 4 * (1 / rhoc1 + 1 / rhoc2) / pow(1 / pow(rhoc1, 1.0 / 3.0) + 1 / pow(rhoc2, 1.0 / 3.0), 3));
            dict.add_number("betaV", 1.0);
            dict.add_number("betaT", 1.0);
        } else if (rule == "Lorentz-Berthelot") {
            // Terms for Lorentz-Berthelot quadratic mixing

            dict.add_number("gammaT", 1.0);
            dict.add_number("gammaV", 1.0);
            dict.add_number("betaV", 1.0);
            dict.add_number("betaT", 1.0);
        } else {
            throw ValueError(format("Your simple mixing rule [%s] was not understood", rule.c_str()));
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
};
// The modifiable parameter library
static MixtureBinaryPairLibrary mixturebinarypairlibrary;
// A fixed parameter library containing the default values
static MixtureBinaryPairLibrary mixturebinarypairlibrary_default;

/// Add a simple mixing rule
void apply_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule) {
    mixturebinarypairlibrary.add_simple_mixing_rule(identifier1, identifier2, rule);
}

std::string get_csv_mixture_binary_pairs() {

    std::vector<std::string> out;
    for (std::map<std::vector<std::string>, std::vector<Dictionary>>::const_iterator it = mixturebinarypairlibrary.binary_pair_map().begin();
         it != mixturebinarypairlibrary.binary_pair_map().end(); ++it) {
        out.push_back(strjoin(it->first, "&"));
    }
    return strjoin(out, ",");
}

std::string get_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& key) {
    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    if (mixturebinarypairlibrary.binary_pair_map().find(CAS) != mixturebinarypairlibrary.binary_pair_map().end()) {
        std::vector<Dictionary>& v = mixturebinarypairlibrary.binary_pair_map()[CAS];
        try {
            if (key == "name1") {
                return v[0].get_string("name1");
            } else if (key == "name2") {
                return v[0].get_string("name2");
            } else if (key == "BibTeX") {
                return v[0].get_string("BibTeX");
            } else if (key == "function") {
                return v[0].get_string("function");
            } else if (key == "type") {
                return v[0].get_string("type");
            } else if (key == "F") {
                return format("%0.16g", v[0].get_double("F"));
            } else if (key == "xi") {
                return format("%0.16g", v[0].get_double("xi"));
            } else if (key == "zeta") {
                return format("%0.16g", v[0].get_double("zeta"));
            } else if (key == "gammaT") {
                return format("%0.16g", v[0].get_double("gammaT"));
            } else if (key == "gammaV") {
                return format("%0.16g", v[0].get_double("gammaV"));
            } else if (key == "betaT") {
                return format("%0.16g", v[0].get_double("betaT"));
            } else if (key == "betaV") {
                return format("%0.16g", v[0].get_double("betaV"));
            } else {
            }
        } catch (...) {
        }
        throw ValueError(format("Could not match the parameter [%s] for the binary pair [%s,%s] - for now this is an error.", key.c_str(),
                                CAS1.c_str(), CAS2.c_str()));
    } else {
        // Sort, see if other order works properly
        std::sort(CAS.begin(), CAS.end());
        if (mixturebinarypairlibrary.binary_pair_map().find(CAS) != mixturebinarypairlibrary.binary_pair_map().end()) {
            throw ValueError(format("Could not match the binary pair [%s,%s] - order of CAS numbers is backwards; found the swapped CAS numbers.",
                                    CAS1.c_str(), CAS2.c_str()));
        } else {
            throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS1.c_str(), CAS2.c_str()));
        }
    }
}
void set_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& key, const double value) {

    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    if (mixturebinarypairlibrary.binary_pair_map().find(CAS) != mixturebinarypairlibrary.binary_pair_map().end()) {
        std::vector<Dictionary>& v = mixturebinarypairlibrary.binary_pair_map()[CAS];
        if (v[0].has_number(key)) {
            v[0].add_number(key, value);
        } else {
            throw ValueError(format("Could not set the parameter [%s] for the binary pair [%s,%s] - for now this is an error", key.c_str(),
                                    CAS1.c_str(), CAS2.c_str()));
        }
    } else {
        // Sort, see if other order works properly
        std::sort(CAS.begin(), CAS.end());
        if (mixturebinarypairlibrary.binary_pair_map().find(CAS) != mixturebinarypairlibrary.binary_pair_map().end()) {
            throw ValueError(format("Could not match the binary pair [%s,%s] - order of CAS numbers is backwards; found the swapped CAS numbers.",
                                    CAS1.c_str(), CAS2.c_str()));
        } else {
            throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS1.c_str(), CAS2.c_str()));
        }
    }
}

std::string get_reducing_function_name(const std::string& CAS1, const std::string& CAS2) {

    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    // Sort the CAS number vector - map is based on sorted CAS codes
    std::sort(CAS.begin(), CAS.end());

    if (mixturebinarypairlibrary.binary_pair_map().find(CAS) != mixturebinarypairlibrary.binary_pair_map().end()) {
        return mixturebinarypairlibrary.binary_pair_map()[CAS][0].get_string("function");
    } else {
        throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS1.c_str(), CAS2.c_str()));
    }
}

void set_interaction_parameters(const std::string& string_data) {
    // JSON-encoded string for binary interaction parameters
    mixturebinarypairlibrary.load_from_string(string_data);
}

/** \brief A container for the departure functions for CoolProp mixtures
 */
class MixtureDepartureFunctionsLibrary
{
   private:
    /// Map from sorted pair of CAS numbers to departure term dictionary.
    std::map<std::string, Dictionary> m_departure_function_map;

   public:
    std::map<std::string, Dictionary>& departure_function_map() {
        // Set the default departure functions if none have been provided yet
        if (m_departure_function_map.size() == 0) {
            load_defaults();
        }
        return m_departure_function_map;
    };

    void load_from_string(const std::string& str) {
        rapidjson::Document doc;
        doc.Parse<0>(str.c_str());
        if (doc.HasParseError()) {
            std::cout << str << std::endl;
            throw ValueError("Unable to parse departure function string");
        }
        load_from_JSON(doc);
    }
    void load_from_JSON(rapidjson::Document& doc) {
        // Iterate over the departure functions in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr) {
            // Get the empty dictionary to be filled in
            Dictionary dict;

            // Populate the dictionary with common terms
            std::string Name = cpjson::get_string(*itr, "Name");
            std::string type = cpjson::get_string(*itr, "type");
            dict.add_string("Name", Name);
            dict.add_string("BibTeX", cpjson::get_string(*itr, "BibTeX"));
            dict.add_string_vector("aliases", cpjson::get_string_array(*itr, "aliases"));
            dict.add_string("type", type);

            // Terms for the power (common to both types)
            dict.add_double_vector("n", cpjson::get_double_array(*itr, "n"));
            dict.add_double_vector("d", cpjson::get_double_array(*itr, "d"));
            dict.add_double_vector("t", cpjson::get_double_array(*itr, "t"));

            // Now we need to load additional terms
            if (!type.compare("GERG-2008")) {
                // Number of terms that are power terms
                dict.add_number("Npower", cpjson::get_double(*itr, "Npower"));
                // Terms for the gaussian
                dict.add_double_vector("eta", cpjson::get_double_array(*itr, "eta"));
                dict.add_double_vector("epsilon", cpjson::get_double_array(*itr, "epsilon"));
                dict.add_double_vector("beta", cpjson::get_double_array(*itr, "beta"));
                dict.add_double_vector("gamma", cpjson::get_double_array(*itr, "gamma"));
            } else if (type == "Gaussian+Exponential") {
                // Number of terms that are power terms
                dict.add_number("Npower", cpjson::get_double(*itr, "Npower"));
                // The decay strength parameters
                dict.add_double_vector("l", cpjson::get_double_array(*itr, "l"));
                // Terms for the gaussian part
                dict.add_double_vector("eta", cpjson::get_double_array(*itr, "eta"));
                dict.add_double_vector("epsilon", cpjson::get_double_array(*itr, "epsilon"));
                dict.add_double_vector("beta", cpjson::get_double_array(*itr, "beta"));
                dict.add_double_vector("gamma", cpjson::get_double_array(*itr, "gamma"));
            } else if (!type.compare("Exponential")) {
                dict.add_double_vector("l", cpjson::get_double_array(*itr, "l"));
            } else {
                throw ValueError(format("It was not possible to parse departure function with type [%s]", type.c_str()));
            }
            // Add the normal name;
            add_one(Name, dict);
            std::vector<std::string> aliases = dict.get_string_vector("aliases");
            // Add the aliases too;
            for (std::vector<std::string>::const_iterator it = aliases.begin(); it != aliases.end(); ++it) {
                // Add the alias;
                add_one(*it, dict);
            }
        }
    }
    void add_one(const std::string& name, Dictionary& dict) {

        // Check if this name is already in use
        std::map<std::string, Dictionary>::iterator it = m_departure_function_map.find(name);
        if (it == m_departure_function_map.end()) {
            // Not in map, add new entry to map with dictionary as value and Name as key
            m_departure_function_map.insert(std::pair<std::string, Dictionary>(name, dict));
        } else {
            if (get_config_bool(OVERWRITE_DEPARTURE_FUNCTION)) {
                // Already there, see http://www.cplusplus.com/reference/map/map/insert/
                m_departure_function_map.erase(it);
                std::pair<std::map<std::string, Dictionary>::iterator, bool> ret;
                ret = m_departure_function_map.insert(std::pair<std::string, Dictionary>(name, dict));
                assert(ret.second == true);
            } else {
                // Error if already in map!
                //
                // Collect all the current names for departure functions for a nicer error message
                std::vector<std::string> names;
                for (std::map<std::string, Dictionary>::const_iterator it = m_departure_function_map.begin(); it != m_departure_function_map.end();
                     ++it) {
                    names.push_back(it->first);
                }
                throw ValueError(format("Name of departure function [%s] is already loaded. Current departure function names are: %s", name.c_str(),
                                        strjoin(names, ",").c_str()));
            }
        }
    }
    // Load the defaults that come from the JSON-encoded string compiled into library
    // as the variable mixture_departure_functions_JSON
    void load_defaults() {
        load_from_string(mixture_departure_functions_JSON);
    }
};
static MixtureDepartureFunctionsLibrary mixturedeparturefunctionslibrary;

DepartureFunction* get_departure_function(const std::string& Name) {
    // Get the dictionary itself
    Dictionary& dict_dep = mixturedeparturefunctionslibrary.departure_function_map()[Name];

    if (dict_dep.is_empty()) {
        throw ValueError(format("Departure function name [%s] seems to be invalid", Name.c_str()));
    }

    // These terms are common
    std::vector<double> n = dict_dep.get_double_vector("n");
    std::vector<double> d = dict_dep.get_double_vector("d");
    std::vector<double> t = dict_dep.get_double_vector("t");

    std::string type_dep = dict_dep.get_string("type");

    if (!type_dep.compare("GERG-2008")) {
        // Number of power terms needed
        int Npower = static_cast<int>(dict_dep.get_number("Npower"));
        // Terms for the gaussian
        std::vector<double> eta = dict_dep.get_double_vector("eta");
        std::vector<double> epsilon = dict_dep.get_double_vector("epsilon");
        std::vector<double> beta = dict_dep.get_double_vector("beta");
        std::vector<double> gamma = dict_dep.get_double_vector("gamma");
        return new GERG2008DepartureFunction(n, d, t, eta, epsilon, beta, gamma, Npower);
    } else if (!type_dep.compare("Exponential")) {
        // Powers of the exponents inside the exponential term
        std::vector<double> l = dict_dep.get_double_vector("l");
        return new ExponentialDepartureFunction(n, d, t, l);
    } else if (!type_dep.compare("Gaussian+Exponential")) {
        // Number of power terms needed
        int Npower = static_cast<int>(dict_dep.get_number("Npower"));
        // Powers of the exponents inside the exponential term
        std::vector<double> l = dict_dep.get_double_vector("l");
        // Terms for the gaussian
        std::vector<double> eta = dict_dep.get_double_vector("eta");
        std::vector<double> epsilon = dict_dep.get_double_vector("epsilon");
        std::vector<double> beta = dict_dep.get_double_vector("beta");
        std::vector<double> gamma = dict_dep.get_double_vector("gamma");
        return new GaussianExponentialDepartureFunction(n, d, t, l, eta, epsilon, beta, gamma, Npower);
    } else {
        throw ValueError();
    }
}
void MixtureParameters::set_mixture_parameters(HelmholtzEOSMixtureBackend& HEOS) {

    std::vector<CoolPropFluid> components = HEOS.get_components();

    std::size_t N = components.size();

    STLMatrix beta_v, gamma_v, beta_T, gamma_T;
    beta_v.resize(N, std::vector<CoolPropDbl>(N, 0));
    gamma_v.resize(N, std::vector<CoolPropDbl>(N, 0));
    beta_T.resize(N, std::vector<CoolPropDbl>(N, 0));
    gamma_T.resize(N, std::vector<CoolPropDbl>(N, 0));

    HEOS.residual_helmholtz->Excess.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            if (i == j) {
                continue;
            }

            std::string CAS1 = components[i].CAS;
            std::vector<std::string> CAS(2, "");
            CAS[0] = components[i].CAS;
            CAS[1] = components[j].CAS;
            std::sort(CAS.begin(), CAS.end());

            // The variable swapped is true if a swap occurred.
            bool swapped = (CAS1.compare(CAS[0]) != 0);

            // ***************************************************
            //         Reducing parameters for binary pair
            // ***************************************************

            if (mixturebinarypairlibrary.binary_pair_map().find(CAS) == mixturebinarypairlibrary.binary_pair_map().end()) {
                throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS[0].c_str(), CAS[1].c_str()));
            }

            // Get a reference to the first matching binary pair in the dictionary
            Dictionary& dict_red = mixturebinarypairlibrary.binary_pair_map()[CAS][0];

            // Get the name of the type being used, one of GERG-2008, Lemmon-xi-zeta, etc.
            std::string type_red = dict_red.get_string("type");

            if (!type_red.compare("GERG-2008")) {
                if (swapped) {
                    beta_v[i][j] = 1 / dict_red.get_number("betaV");
                    beta_T[i][j] = 1 / dict_red.get_number("betaT");
                } else {
                    beta_v[i][j] = dict_red.get_number("betaV");
                    beta_T[i][j] = dict_red.get_number("betaT");
                }
                gamma_v[i][j] = dict_red.get_number("gammaV");
                gamma_T[i][j] = dict_red.get_number("gammaT");
            } else if (!type_red.compare("Lemmon-xi-zeta")) {
                LemmonAirHFCReducingFunction::convert_to_GERG(components, i, j, dict_red, beta_T[i][j], beta_v[i][j], gamma_T[i][j], gamma_v[i][j]);
            } else {
                throw ValueError(format("type [%s] for reducing function for pair [%s, %s] is invalid", type_red.c_str(),
                                        dict_red.get_string("Name1").c_str(), dict_red.get_string("Name2").c_str()));
            }
            /*
            if (i == 0){
                std::cout << format("betaT %10.9Lg gammaT %10.9Lg betaV %10.9Lg gammaV %10.9Lg %s",
                                    beta_T[i][j], gamma_T[i][j], beta_v[i][j], gamma_v[i][j], get_mixture_binary_pair_data(CAS[0],CAS[1],"gammaT").c_str()) << std::endl;
            }
            */

            // ***************************************************
            //     Departure functions used in excess term
            // ***************************************************

            // Set the scaling factor F for the excess term
            HEOS.residual_helmholtz->Excess.F[i][j] = dict_red.get_number("F");

            if (std::abs(HEOS.residual_helmholtz->Excess.F[i][j]) < DBL_EPSILON) {
                // Empty departure function that will just return 0
                std::vector<double> n(1, 0), d(1, 1), t(1, 1), l(1, 0);
                HEOS.residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(new ExponentialDepartureFunction(n, d, t, l));
                continue;
            }

            // Get the name of the departure function to be used for this binary pair
            std::string Name = CoolProp::get_reducing_function_name(components[i].CAS, components[j].CAS);

            HEOS.residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(get_departure_function(Name));
        }
    }
    // We have obtained all the parameters needed for the reducing function, now set the reducing function for the mixture
    HEOS.Reducing = shared_ptr<ReducingFunction>(new GERG2008ReducingFunction(components, beta_v, gamma_v, beta_T, gamma_T));
}

void parse_HMX_BNC(const std::string& s, std::vector<REFPROP_binary_element>& BIP, std::vector<REFPROP_departure_function>& functions) {
    // Capture the betas, gammas, Fij, models
    bool block_started = false;
    std::size_t i_started = 0, i_ended = 0, i = 0;
    std::vector<std::string> lines = strsplit(s, '\n');
    for (std::vector<std::string>::iterator it = lines.begin(); it != lines.end(); ++it) {
        if (strstartswith(strstrip(*it), "#BNC")) {
            block_started = true;
            i_started = i + 1;
        }
        if (block_started && strstrip(*it).empty()) {
            i_ended = i - 1;
            break;
        }
        i++;
    }
    // Find the first line with a !
    for (i = i_started; i < i_ended; ++i) {
        if (strstrip(lines[i]) == "!") {
            i_started = i;
            break;
        }
    }
    // Find all the lines that '!' - these are delimiters
    std::vector<std::pair<std::size_t, std::size_t>> bounds;  // The left and right indices (inclusive) that form a binary pair
    std::size_t last_excalamation = i_started;
    for (i = i_started; i <= i_ended; ++i) {
        if (strstrip(lines[i]) == "!") {
            bounds.push_back(std::make_pair<std::size_t, std::size_t>(last_excalamation + 1, i - 1));
            last_excalamation = i;
        }
    }
    // Parse each chunk
    std::vector<REFPROP_binary_element> chunks;
    for (std::vector<std::pair<std::size_t, std::size_t>>::iterator it = bounds.begin(); it != bounds.end(); ++it) {
        REFPROP_binary_element bnc;
        for (std::size_t i = it->first; i <= it->second; ++i) {
            // Store comments
            if (strstartswith(lines[i], "?")) {
                bnc.comments.push_back(lines[i]);
                continue;
            }
            // Parse the line with the thermo BIP
            if (lines[i].find("/") > 0) {
                // Split at ' '
                std::vector<std::string> bits = strsplit(strstrip(lines[i]), ' ');
                // Remove empty elements
                for (std::size_t j = bits.size() - 1; j > 0; --j) {
                    if (bits[j].empty()) {
                        bits.erase(bits.begin() + j);
                    }
                }
                // Get the line that contains the thermo BIP
                if (bits[0].find("/") > 0 && bits[1].size() == 3) {
                    std::vector<std::string> theCAS = strsplit(bits[0], '/');
                    bnc.CAS1 = theCAS[0];
                    bnc.CAS2 = theCAS[1];
                    bnc.model = bits[1];
                    bnc.betaT = string2double(bits[2]);
                    bnc.gammaT = string2double(bits[3]);
                    bnc.betaV = string2double(bits[4]);
                    bnc.gammaV = string2double(bits[5]);
                    bnc.Fij = string2double(bits[6]);
                    break;
                } else if (strstrip(bits[0]) == "CAS#") {
                    break;
                } else {
                    throw CoolProp::ValueError(format("Unable to parse binary interaction line: %s", lines[i]));
                }
            }
        }
        if (!bnc.CAS1.empty()) {
            BIP.push_back(bnc);
        }
    }

    // ****************************************
    //      Parse the departure functions
    // ****************************************
    for (std::size_t i = i_ended + 1; i < lines.size(); ++i) {
        std::size_t j_end;
        // Find the end of this block
        for (j_end = i + 1; j_end < lines.size(); ++j_end) {
            if (strstrip(lines[j_end]).empty()) {
                j_end -= 1;
                break;
            }
        }

        if (strstartswith(lines[i], "#MXM")) {
            REFPROP_departure_function dep;
            dep.Npower = -1;
            dep.model = std::string(lines[i + 1].begin(), lines[i + 1].begin() + 3);
            dep.comments.push_back(lines[i + 1]);
            for (std::size_t j = i + 2; j <= j_end; ++j) {
                if (strstartswith(strstrip(lines[j]), "?")) {
                    dep.comments.push_back(lines[j]);
                    continue;
                }
                if (strstartswith(strstrip(lines[j]), "!")) {
                    j += 2;  // Skip the BIP here, not used
                    continue;
                }
                std::vector<std::string> bits = strsplit(lines[j], ' ');
                // Remove empty elements
                for (std::size_t k = bits.size() - 1; k > 0; --k) {
                    if (bits[k].empty()) {
                        bits.erase(bits.begin() + k);
                    }
                }

                if (dep.Npower < 0) {  // Not extracted yet, let's do it now
                    // Extract the number of terms
                    dep.Npower = static_cast<short>(strtol(bits[0].c_str(), NULL, 10));
                    dep.Nterms_power = static_cast<short>(strtol(bits[1].c_str(), NULL, 10));
                    dep.Nspecial = static_cast<short>(strtol(bits[3].c_str(), NULL, 10));
                    dep.Nterms_special = static_cast<short>(strtol(bits[4].c_str(), NULL, 10));
                } else {
                    dep.a.push_back(string2double(bits[0]));
                    dep.t.push_back(string2double(bits[1]));
                    dep.d.push_back(string2double(bits[2]));
                    // Extracting "polynomial" terms
                    if (dep.Nterms_power == 4) {
                        dep.e.push_back(string2double(bits[3]));
                    }
                    if (dep.Nspecial > 0) {
                        if (dep.a.size() - 1 < dep.Npower) {
                            dep.eta.push_back(0);
                            dep.epsilon.push_back(0);
                            dep.beta.push_back(0);
                            dep.gamma.push_back(0);
                        } else {
                            // Extracting "special" terms
                            dep.eta.push_back(string2double(bits[3]));
                            dep.epsilon.push_back(string2double(bits[4]));
                            dep.beta.push_back(string2double(bits[5]));
                            dep.gamma.push_back(string2double(bits[6]));
                        }
                    }
                }
            }
            functions.push_back(dep);
        }
    }
}

void set_departure_functions(const std::string& string_data) {
    if (string_data.find("#MXM") != std::string::npos) {
        // REFPROP HMX.BNC file was provided
        std::vector<REFPROP_binary_element> BIP;
        std::vector<REFPROP_departure_function> functions;
        parse_HMX_BNC(string_data, BIP, functions);

        {
            rapidjson::Document doc;
            doc.SetArray();
            for (std::vector<REFPROP_binary_element>::const_iterator it = BIP.begin(); it < BIP.end(); ++it) {
                rapidjson::Value el;
                el.SetObject();
                el.AddMember("CAS1", rapidjson::Value(it->CAS1.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                el.AddMember("CAS2", rapidjson::Value(it->CAS2.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                el.AddMember("Name1", "??", doc.GetAllocator());
                el.AddMember("Name2", "??", doc.GetAllocator());
                el.AddMember("betaT", it->betaT, doc.GetAllocator());
                el.AddMember("gammaT", it->gammaT, doc.GetAllocator());
                el.AddMember("betaV", it->betaV, doc.GetAllocator());
                el.AddMember("gammaV", it->gammaV, doc.GetAllocator());
                el.AddMember("F", it->Fij, doc.GetAllocator());
                el.AddMember("function", rapidjson::Value(it->model.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                std::string tex_string = "(from HMX.BNC format)::" + strjoin(it->comments, "\n");
                el.AddMember("BibTeX", rapidjson::Value(tex_string.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                doc.PushBack(el, doc.GetAllocator());
            }
            mixturebinarypairlibrary.load_from_JSON(doc);
        }
        {
            rapidjson::Document doc;
            doc.SetArray();
            for (std::vector<REFPROP_departure_function>::const_iterator it = functions.begin(); it < functions.end(); ++it) {
                rapidjson::Value el;
                el.SetObject();
                el.AddMember("Name", rapidjson::Value(it->model.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                std::vector<std::string> aliases;
                cpjson::set_string_array("aliases", aliases, el, doc);
                cpjson::set_double_array("n", it->a, el, doc);
                cpjson::set_double_array("d", it->d, el, doc);
                cpjson::set_double_array("t", it->t, el, doc);
                if (it->Nterms_special > 0 || it->Nterms_power == 3) {
                    el.AddMember("type", "GERG-2008", doc.GetAllocator());
                    el.AddMember("Npower", it->Npower, doc.GetAllocator());
                    if (it->Nterms_power == 3 && it->Nspecial == 0) {
                        std::vector<double> zeros(it->a.size(), 0);
                        cpjson::set_double_array("eta", zeros, el, doc);
                        cpjson::set_double_array("epsilon", zeros, el, doc);
                        cpjson::set_double_array("beta", zeros, el, doc);
                        cpjson::set_double_array("gamma", zeros, el, doc);
                    } else {
                        cpjson::set_double_array("eta", it->eta, el, doc);
                        cpjson::set_double_array("epsilon", it->epsilon, el, doc);
                        cpjson::set_double_array("beta", it->beta, el, doc);
                        cpjson::set_double_array("gamma", it->gamma, el, doc);
                    }
                } else {
                    el.AddMember("type", "Exponential", doc.GetAllocator());
                    cpjson::set_double_array("l", it->e, el, doc);
                }

                std::string tex_string = "(from HMX.BNC format)::" + strjoin(it->comments, "\n");
                el.AddMember("BibTeX", rapidjson::Value(tex_string.c_str(), doc.GetAllocator()).Move(), doc.GetAllocator());
                doc.PushBack(el, doc.GetAllocator());
            }
            mixturedeparturefunctionslibrary.load_from_JSON(doc);
        }
    } else {
        // JSON-encoded string for departure functions
        mixturedeparturefunctionslibrary.load_from_string(string_data);
    }
}

} /* namespace CoolProp */
