
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
    /// Map from CAS code to JSON instance.  For pseudo-pure fluids, use name in place of CAS code since no CASE number is defined for mixtures
    std::map<std::size_t, IncompressibleFluid> fluid_map;
    std::vector<std::string> name_vector;
    std::map<std::string, std::size_t> string_to_index_map;
    bool _is_empty;
protected:

    /// Parse the viscosity
    void parse_viscosity(rapidjson::Value &viscosity, IncompressibleFluid & fluid)
    {
        if (viscosity.HasMember("type")){
            std::string type = cpjson::get_string(viscosity, "type");
            if (!type.compare("polynomial")){
                fluid.viscosity.type = CoolProp::IncompressibleViscosityVariables::INCOMPRESSIBLE_VISCOSITY_POLYNOMIAL;
                fluid.viscosity.poly.coeffs = cpjson::get_double_array(viscosity["coeffs"]);
                return;
            }
            else{
                throw ValueError(format("viscosity type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
            }
        }
        else{
            throw ValueError(format("viscosity does not have \"type\" for fluid %s", fluid.name.c_str()));
        }
    };

    /// Parse the conductivity
    void parse_conductivity(rapidjson::Value &conductivity, IncompressibleFluid & fluid)
    {
        if (conductivity.HasMember("type")){
            std::string type = cpjson::get_string(conductivity, "type");
            if (!type.compare("polynomial")){
                fluid.conductivity.type = CoolProp::IncompressibleConductivityVariables::INCOMPRESSIBLE_CONDUCTIVITY_POLYNOMIAL;
                fluid.conductivity.poly.coeffs = cpjson::get_double_array(conductivity["coeffs"]);
                return;
            }
            else{
                throw ValueError(format("conductivity type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
            }
        }
        else{
            throw ValueError(format("conductivity does not have \"type\" for fluid %s", fluid.name.c_str()));
        }
    };

    /// Parse the specific_heat
    void parse_specific_heat(rapidjson::Value &specific_heat, IncompressibleFluid & fluid)
    {
        if (specific_heat.HasMember("type")){
            std::string type = cpjson::get_string(specific_heat, "type");
            if (!type.compare("polynomial")){
                fluid.specific_heat.type = CoolProp::IncompressibleSpecificHeatVariables::INCOMPRESSIBLE_SPECIFIC_HEAT_POLYNOMIAL; return;
                fluid.specific_heat.poly.coeffs = cpjson::get_double_array(specific_heat["coeffs"]);
            }
            else{
                throw ValueError(format("specific_heat type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
            }
        }
        else{
            throw ValueError(format("specific_heat does not have \"type\" for fluid %s", fluid.name.c_str()));
        }
    };

    /// Parse the density
    void parse_density(rapidjson::Value &density, IncompressibleFluid & fluid)
    {
        if (density.HasMember("type")){
            std::string type = cpjson::get_string(density, "type");
            if (!type.compare("polynomial")){
                fluid.density.type = CoolProp::IncompressibleDensityVariables::INCOMPRESSIBLE_DENSITY_POLYNOMIAL; return;
                fluid.density.poly.coeffs = cpjson::get_double_array(density["coeffs"]);
            }
            else{
                throw ValueError(format("density type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
            }
        }
        else{
            throw ValueError(format("density does not have \"type\" for fluid %s", fluid.name.c_str()));
        }
    };

    /// Validate the fluid file that was just constructed
    void validate(IncompressibleFluid & fluid)
    {
    }
public:

    // Default constructor;
    JSONIncompressibleLibrary(){
        _is_empty = true;
    };
    bool is_empty(void){ return _is_empty;};

    /// Add all the fluid entries in the rapidjson::Value instance passed in
    void add_many(rapidjson::Value &listing)
    {
        for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr)
        {
            add_one(*itr);
        }
    };
    void add_one(rapidjson::Value &fluid_json)
    {
        _is_empty = false;

        // Get the next index for this fluid
        std::size_t index = fluid_map.size();

        // Add index->fluid mapping
        fluid_map[index] = IncompressibleFluid();

        // Create an instance of the fluid
        IncompressibleFluid &fluid = fluid_map[index];

        fluid.name = cpjson::get_string(fluid_json, "name");
        fluid.Tmin = cpjson::get_double(fluid_json, "Tmin");
        fluid.Tmax = cpjson::get_double(fluid_json, "Tmax");

        parse_conductivity(fluid_json["conductivity"], fluid);
        parse_density(fluid_json["density"], fluid);
        parse_viscosity(fluid_json["viscosity"], fluid);
        parse_specific_heat(fluid_json["specific_heat"], fluid);

        // Add name->index mapping
        string_to_index_map[fluid.name] = index;

    };
    /// Get an IncompressibleFluid instance stored in this library
    /**
    @param name Name of the fluid
    */
    IncompressibleFluid& get(std::string key)
    {
        std::map<std::string, std::size_t>::iterator it;
        // Try to find it
        it = string_to_index_map.find(key);
        // If it is found
        if (it != string_to_index_map.end()){
            return get(it->second);
        }
        else{
            throw ValueError(format("key [%s] was not found in string_to_index_map in JSONIncompressibleLibrary",key.c_str()));
        }
    };
    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key The index of the fluid in the map
    */
    IncompressibleFluid& get(std::size_t key)
    {
        std::map<std::size_t, IncompressibleFluid>::iterator it;
        // Try to find it
        it = fluid_map.find(key);
        // If it is found
        if (it != fluid_map.end()){
            return it->second;
        }
        else{
            throw ValueError(format("key [%d] was not found in JSONIncompressibleLibrary",key));
        }
    };
    /// Return a comma-separated list of fluid names
    std::string get_fluid_list(void)
    {
        return strjoin(name_vector, ",");
    };
};

/// Get a reference to the library instance
JSONIncompressibleLibrary & get_incompressible_library(void);

/// Get a comma-separated-list of incompressible fluids that are included
std::string get_incompressible_list(void);

/// Get the fluid structure returned as a reference
IncompressibleFluid& get_incompressible_fluid(std::string fluid_string);

} /* namespace CoolProp */
#endif
