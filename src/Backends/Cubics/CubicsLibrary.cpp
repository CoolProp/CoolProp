#include <string>
#include <map>
#include "CubicsLibrary.h"
#include "all_cubics_JSON.h" // Makes a std::string variable called all_cubics_JSON
#include "rapidjson_include.h"

namespace CoolProp{

class CubicsLibrary{
private:
    std::map<std::string, CubicsValues> fluid_map;
    std::map<std::string, std::string> aliases_map;
    bool empty; // Is empty
public:
    CubicsLibrary() : empty(true) {};
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
        // Try to find it
        std::map<std::string, CubicsValues>::iterator it = fluid_map.find(identifier);
        // If it is found
        if (it != fluid_map.end()) {
            return it->second;
        } else {
            std::map<std::string, std::string>::iterator italias = aliases_map.find(identifier);
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
static CubicsLibrary library;

void load_cubics_library(){
    rapidjson::Document dd;
    // This json formatted string comes from the all_cubics_JSON.h header which is a C++-escaped version of the JSON file
    dd.Parse<0>(all_cubics_JSON.c_str());
    if (dd.HasParseError()){
        throw ValueError("Unable to load all_cubics_JSON.json");
    } else{
        try{
            library.add_many(dd);
        }catch(std::exception &e){std::cout << e.what() << std::endl;}
    }
}

CubicsValues get_cubic_values(const std::string &identifier){
    if (library.is_empty()){ load_cubics_library(); }
    return library.get(identifier);
}

} /* namepace CoolProp */
