#include <string>
#include <map>
#include "CubicsLibrary.h"
#include "all_cubics_JSON.h" // Makes a std::string variable called all_cubics_JSON
#include "cubic_fluids_schema_JSON.h" // Makes a std::string variable called cubic_fluids_schema_JSON
#include "rapidjson_include.h"
#include "CPstrings.h"
#include "CoolProp.h"

namespace CoolProp{
namespace CubicLibrary{

class CubicsLibraryClass{
private:
    std::map<std::string, CubicsValues> fluid_map;
    std::map<std::string, std::string> aliases_map;
    bool empty; // Is empty
public:
    CubicsLibraryClass() : empty(true) {
        // This JSON formatted string comes from the all_cubics_JSON.h header which is a C++-escaped version of the JSON file
        add_fluids_as_JSON(all_cubics_JSON);
    };
    bool is_empty(){ return empty; };
    int add_many(rapidjson::Value &listing)
    {
        int counter = 0;
        for (rapidjson::Value:: ValueIterator itr = listing.Begin();
                itr != listing.End(); ++itr) {
            CubicsValues val;
            val.Tc = cpjson::get_double(*itr, "Tc");
            val.pc = cpjson::get_double(*itr, "pc");
            val.acentric = cpjson::get_double(*itr, "acentric");
            val.molemass = cpjson::get_double(*itr, "molemass");
            val.name = cpjson::get_string(*itr, "name");
            val.aliases = cpjson::get_string_array(*itr, "aliases");
            fluid_map.insert(std::pair<std::string, CubicsValues>(val.name, val) );

            for (std::vector<std::string>::const_iterator it = val.aliases.begin(); it != val.aliases.end(); ++it){
                if (aliases_map.find(*it) == aliases_map.end()){
                    // It's not already in aliases map
                    aliases_map.insert(std::pair<std::string, std::string>(*it, val.name) );
                }
            }
            counter ++;
        }
        return counter;
    };
    CubicsValues get(const std::string & identifier){
        std::string uppercase_identifier = upper(identifier);
        // Try to find it
        std::map<std::string, CubicsValues>::iterator it = fluid_map.find(uppercase_identifier);
        // If it is found
        if (it != fluid_map.end()) {
            return it->second;
        } else {
            std::map<std::string, std::string>::iterator italias = aliases_map.find(uppercase_identifier);
            if (italias != aliases_map.end()){
                // Alias was found, use it to get the fluid name, and then the cubic values
                return fluid_map.find(italias->second)->second;
            }
            else{
                throw ValueError(format("Fluid identifier [%s] was not found in CubicsLibrary", identifier.c_str()));
            }
        }
    };
};
static CubicsLibraryClass library;

    
void add_fluids_as_JSON(const std::string &JSON)
{
    // First we validate the json string against the schema;
    std::string errstr;
    cpjson::schema_validation_code val_code = cpjson::validate_schema(cubic_fluids_schema_JSON, JSON, errstr);
    // Then we check the validation code
    if (val_code == cpjson::SCHEMA_VALIDATION_OK){
        rapidjson::Document dd;
        
        dd.Parse<0>(JSON.c_str());
        if (dd.HasParseError()){
            throw ValueError("Unable to load all_cubics_JSON.json");
        } else{
            try{
                library.add_many(dd);
            }catch(std::exception &e){std::cout << e.what() << std::endl;}
        }
    }
    else{
        if (get_debug_level() > 0){ throw ValueError(format("Unable to load cubics library with error: %s", errstr.c_str())); }
    }
}

CubicLibrary::CubicsValues get_cubic_values(const std::string &identifier){
    return library.get(identifier);
}
std::string get_cubic_fluids_schema(){
    return cubic_fluids_schema_JSON;
}

} /* namepace CriticalLibrary */
} /* namepace CoolProp */
