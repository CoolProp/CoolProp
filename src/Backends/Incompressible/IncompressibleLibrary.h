
#ifndef INCOMPRESSIBLELIBRARY_H
#define INCOMPRESSIBLELIBRARY_H

#include "IncompressibleFluid.h"

#include "rapidjson/rapidjson_include.h"

#include <map>
#include <algorithm>

namespace CoolProp{

// Forward declaration of the necessary debug function to avoid including the whole header
extern int get_debug_level();

/// A container for the fluid parameters for the incompressible fluids
/**
This container holds copies of all of the fluid instances for the fluids that are loaded in incompressible.
New fluids can be added by passing in a rapidjson::Value instance to the add_one function, or
a rapidjson array of fluids to the add_many function.
*/
class JSONIncompressibleLibrary
{
    /// Map from CAS code to JSON instance.
	/** This is not practical for the incomressibles, the CAS may not be
	 *  defined for blends of heat transfer fluids and solutions.
     */
    std::map<std::size_t, IncompressibleFluid> fluid_map;
    std::vector<std::string> name_vector;
    std::map<std::string, std::size_t> string_to_index_map;
    bool _is_empty;

protected:
    /// A general function to parse the json files that hold the coefficient matrices
    IncompressibleData parse_coefficients(rapidjson::Value &obj, std::string id, bool vital);
    double parse_value(rapidjson::Value &obj, std::string id, bool vital, double def);

public:
    // Default constructor;
    JSONIncompressibleLibrary(){ _is_empty = true;};

    bool is_empty(void){ return _is_empty;};

    /// Add all the fluid entries in the rapidjson::Value instance passed in
    void add_many(rapidjson::Value &listing);
    void add_one(rapidjson::Value &fluid_json);

    /// Get an IncompressibleFluid instance stored in this library
    /**
    @param name Name of the fluid
    */
    IncompressibleFluid& get(std::string key);

    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key The index of the fluid in the map
    */
    IncompressibleFluid& get(std::size_t key);

    /// Return a comma-separated list of fluid names
    std::string get_fluid_list(void){ return strjoin(name_vector, ",");};
};

/// Get a reference to the library instance
JSONIncompressibleLibrary & get_incompressible_library(void);

/// Get a comma-separated-list of incompressible fluids that are included
std::string get_incompressible_list(void);

/// Get the fluid structure returned as a reference
IncompressibleFluid& get_incompressible_fluid(std::string fluid_string);

} /* namespace CoolProp */
#endif
