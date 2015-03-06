
#include "FluidLibrary.h"
#include "all_fluids_JSON.h" // Makes a std::string variable called all_fluids_JSON

namespace CoolProp{

static JSONFluidLibrary library;

void load()
{
    rapidjson::Document dd;
    // This json formatted string comes from the all_fluids_JSON.h header which is a C++-escaped version of the JSON file
    dd.Parse<0>(all_fluids_JSON.c_str());
    if (dd.HasParseError()){
        throw ValueError("Unable to load all_fluids.json");
    } else{
        try{library.add_many(dd);}catch(std::exception &e){std::cout << e.what() << std::endl;}
    }
}

JSONFluidLibrary & get_library(void){
    if (library.is_empty()){ load(); }
    return library;
}

CoolPropFluid get_fluid(const std::string &fluid_string){
    if (library.is_empty()){ load(); }
    return library.get(fluid_string);
}

std::string get_fluid_list(void){
    if (library.is_empty()){ load(); }
    return library.get_fluid_list();
};

void set_fluid_enthalpy_entropy_offset(const std::string &fluid, double delta_a1, double delta_a2, const std::string &ref){
    if (library.is_empty()){ load(); }
    library.set_fluid_enthalpy_entropy_offset(fluid, delta_a1, delta_a2, ref);
}

} /* namespace CoolProp */