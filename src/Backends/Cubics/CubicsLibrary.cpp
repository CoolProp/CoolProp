#include <string>
#include <map>
#include "CubicsLibrary.h"
#include "all_cubics_JSON.h"           // Makes a std::string variable called all_cubics_JSON
#include "cubic_fluids_schema_JSON.h"  // Makes a std::string variable called cubic_fluids_schema_JSON
#include "rapidjson_include.h"
#include "CPstrings.h"
#include "CoolProp.h"
#include "Configuration.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"

namespace CoolProp {
namespace CubicLibrary {

class CubicsLibraryClass
{
   private:
    std::map<std::string, CubicsValues> fluid_map;
    std::map<std::string, std::string> aliases_map;
    bool empty;  // Is empty
   public:
    CubicsLibraryClass() : empty(true) {
        // This JSON formatted string comes from the all_cubics_JSON.h header which is a C++-escaped version of the JSON file
        add_fluids_as_JSON(all_cubics_JSON);
    };
    bool is_empty() {
        return empty;
    };
    int add_many(rapidjson::Value& listing) {
        int counter = 0;
        for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr) {
            CubicsValues val;
            val.Tc = cpjson::get_double(*itr, "Tc");
            val.pc = cpjson::get_double(*itr, "pc");
            val.acentric = cpjson::get_double(*itr, "acentric");
            val.molemass = cpjson::get_double(*itr, "molemass");
            val.name = cpjson::get_string(*itr, "name");
            val.aliases = cpjson::get_string_array(*itr, "aliases");
            val.CAS = cpjson::get_string(*itr, "CAS");
            if (itr->HasMember("rhomolarc") && (*itr)["rhomolarc"].IsNumber()) {
                val.rhomolarc = cpjson::get_double(*itr, "rhomolarc");
            }
            if (itr->HasMember("alpha") && (*itr)["alpha"].IsObject()) {
                rapidjson::Value& alpha = (*itr)["alpha"];
                val.alpha_type = cpjson::get_string(alpha, "type");
                val.alpha_coeffs = cpjson::get_double_array(alpha, "c");
            } else {
                val.alpha_type = "default";
            }
            if (itr->HasMember("alpha0") && (*itr)["alpha0"].IsArray()) {
                val.alpha0 = JSONFluidLibrary::parse_alpha0((*itr)["alpha0"]);
            }
            std::pair<std::map<std::string, CubicsValues>::iterator, bool> ret;
            ret = fluid_map.insert(std::pair<std::string, CubicsValues>(upper(val.name), val));
            if (ret.second == false && get_config_bool(OVERWRITE_FLUIDS)) {
                // Already there, see http://www.cplusplus.com/reference/map/map/insert/
                fluid_map.erase(ret.first);
                ret = fluid_map.insert(std::pair<std::string, CubicsValues>(upper(val.name), val));
                if (get_debug_level() > 0) {
                    std::cout << "added the cubic fluid: " + val.name << std::endl;
                }
                assert(ret.second == true);
            }

            for (std::vector<std::string>::const_iterator it = val.aliases.begin(); it != val.aliases.end(); ++it) {
                if (aliases_map.find(*it) == aliases_map.end()) {
                    // It's not already in aliases map
                    aliases_map.insert(std::pair<std::string, std::string>(*it, upper(val.name)));
                }
            }
            counter++;
        }
        return counter;
    };
    CubicsValues get(const std::string& identifier) {
        std::string uppercase_identifier = upper(identifier);
        // Try to find it
        std::map<std::string, CubicsValues>::iterator it = fluid_map.find(uppercase_identifier);
        // If it is found
        if (it != fluid_map.end()) {
            return it->second;
        } else {
            std::map<std::string, std::string>::iterator italias = aliases_map.find(uppercase_identifier);
            if (italias != aliases_map.end()) {
                // Alias was found, use it to get the fluid name, and then the cubic values
                return fluid_map.find(italias->second)->second;
            } else {
                throw ValueError(format("Fluid identifier [%s] was not found in CubicsLibrary", uppercase_identifier.c_str()));
            }
        }
    };
    std::string get_fluids_list() {
        std::vector<std::string> out;
        for (std::map<std::string, CubicsValues>::const_iterator it = fluid_map.begin(); it != fluid_map.end(); ++it) {
            out.push_back(it->first);
        }
        return strjoin(out, ",");
    }
};
static CubicsLibraryClass library;

void add_fluids_as_JSON(const std::string& JSON) {
    // First we validate the json string against the schema;
    std::string errstr;
    cpjson::schema_validation_code val_code = cpjson::validate_schema(cubic_fluids_schema_JSON, JSON, errstr);
    // Then we check the validation code
    if (val_code == cpjson::SCHEMA_VALIDATION_OK) {
        rapidjson::Document dd;

        dd.Parse<0>(JSON.c_str());
        if (dd.HasParseError()) {
            throw ValueError("Cubics JSON is not valid JSON");
        } else {
            try {
                library.add_many(dd);
            } catch (std::exception& e) {
                throw ValueError(format("Unable to load cubics library with error: %s", errstr.c_str()));
            }
        }
    } else {
        throw ValueError(format("Unable to validate cubics library against schema with error: %s", errstr.c_str()));
    }
}

CubicLibrary::CubicsValues get_cubic_values(const std::string& identifier) {
    return library.get(identifier);
}
std::string get_cubic_fluids_schema() {
    return cubic_fluids_schema_JSON;
}
std::string get_cubic_fluids_list() {
    return library.get_fluids_list();
}

}  // namespace CubicLibrary
}  // namespace CoolProp
