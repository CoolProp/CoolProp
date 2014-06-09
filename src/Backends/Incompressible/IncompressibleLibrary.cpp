#include "IncompressibleLibrary.h"
#include "all_incompressibles_JSON.h" // Makes a std::string variable called all_fluids_JSON

namespace CoolProp{

static JSONIncompressibleLibrary library;

void load_incompressible_library()
{
	rapidjson::Document dd;
    // This json formatted string comes from the all_fluids_JSON.h header which is a C++-escaped version of the JSON file
    dd.Parse<0>(all_incompressibles_JSON.c_str());
	if (dd.HasParseError()){
        throw ValueError("Unable to load all_fluids.json");
    } else{
        try{library.add_many(dd);}catch(std::exception &e){std::cout << e.what() << std::endl;}
    }
}

JSONIncompressibleLibrary & get_incompressible_library(void){
	if (library.is_empty()){ load_incompressible_library(); }
	return library;
}

IncompressibleFluid& get_incompressible_fluid(std::string fluid_string){
    if (library.is_empty()){ load_incompressible_library(); }
    return library.get(fluid_string);
}

std::string get_incompressible_list(void){
    if (library.is_empty()){ load_incompressible_library(); }
    return library.get_fluid_list();
};

} /* namespace CoolProp */
