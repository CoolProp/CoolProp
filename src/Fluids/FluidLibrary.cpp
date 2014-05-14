
#include "FluidLibrary.h"
#include "../Backends/ReducingFunctions.h"
#include "all_fluids_JSON.h" // Makes a std::string variable called all_fluids_JSON

namespace CoolProp{

static JSONFluidLibrary library;

void load()
{
	rapidjson::Document dd;
    // This json formatted string comes from the all_fluids_JSON.h header which is a C++-escaped version of the JSON file
    dd.Parse<0>(all_fluids_JSON.c_str());
	if (dd.HasParseError()){throw ValueError("Unable to load all_fluids.json");} else{library.add_many(dd);}
}

JSONFluidLibrary & get_library(void){
	if (library.is_empty()){ load(); }
	return library;
}

CoolPropFluid& get_fluid(std::string fluid_string){
    if (library.is_empty()){ load(); }
    return library.get(fluid_string);
}

std::string get_fluid_list(void){
    if (library.is_empty()){ load(); }
    return library.get_fluid_list();
};

} /* namespace CoolProp */