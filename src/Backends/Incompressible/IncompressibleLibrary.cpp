#include "IncompressibleLibrary.h"
#include "MatrixMath.h"
#include "rapidjson/rapidjson_include.h"
#include "all_incompressibles_JSON.h" // Makes a std::string variable called all_incompressibles_JSON

namespace CoolProp{

/// A general function to parse the json files that hold the coefficient matrices
IncompressibleData JSONIncompressibleLibrary::parse_coefficients(rapidjson::Value &obj, std::string id, bool vital){
	IncompressibleData fluidData;
	if (obj.HasMember(id.c_str())) {
		//rapidjson::Value value = obj[id.c_str()];
		if (obj[id.c_str()].HasMember("type")){
			if (obj[id.c_str()].HasMember("coeffs")){
				std::string type = cpjson::get_string(obj[id.c_str()], "type");
				if (!type.compare("polynomial")){
					fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
					fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(obj[id.c_str()]["coeffs"]));
					return fluidData;
				}
				else if (!type.compare("exponential")){
					fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL;
					fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(obj[id.c_str()]["coeffs"]));
					return fluidData;
				}
				else if (!type.compare("exppolynomial")){
					fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL;
					fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(obj[id.c_str()]["coeffs"]));
					return fluidData;
				}
				else if (!type.compare("expoffset")){
					fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPOFFSET;
					fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(obj[id.c_str()]["coeffs"]));
					return fluidData;
				}
				else if (!type.compare("polyoffset")){
					fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYOFFSET;
					fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(obj[id.c_str()]["coeffs"]));
					return fluidData;
				}
				else if (vital){
					throw ValueError(format("The type [%s] is not understood for [%s] of incompressible fluids. Please check your JSON file.", type.c_str(), id.c_str()));
				}
				else{
                    std::cout << format("The type [%s] is not understood for [%s] of incompressible fluids. Please check your JSON file.\n", type.c_str(), id.c_str());
				}
			}
			else{
				throw ValueError(format("Your file does not have an entry for \"coeffs\" in [%s], which is vital for this function.", id.c_str()));
			}
		}
		else{
			throw ValueError(format("Your file does not have an entry for \"type\" in [%s], which is vital for this function.", id.c_str()));
		}
	}
	else{
		if (vital) {
			throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
		}
	}
	return fluidData;
}

/// Get a double from the JSON storage if it is defined, otherwise return def
double JSONIncompressibleLibrary::parse_value(rapidjson::Value &obj, std::string id, bool vital, double def=0.0){
	if (obj.HasMember(id.c_str())) {
        return cpjson::get_double(obj, id);
    }
	else{
		if (vital) {
			throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
		}
		else{
			return def;
		}
	}
}

/// Add all the fluid entries in the rapidjson::Value instance passed in
void JSONIncompressibleLibrary::add_many(rapidjson::Value &listing) {
	for (rapidjson::Value::ValueIterator itr = listing.Begin();
			itr != listing.End(); ++itr) {
		add_one(*itr);
	}
};

void JSONIncompressibleLibrary::add_one(rapidjson::Value &fluid_json) {
	_is_empty = false;

	// Get the next index for this fluid
	std::size_t index = fluid_map.size();

	// Add index->fluid mapping
	fluid_map[index] = IncompressibleFluid();

	// Create an instance of the fluid
	IncompressibleFluid &fluid = fluid_map[index];
    try
    {

	    fluid.setName(cpjson::get_string(fluid_json, "name"));
	    fluid.setDescription(cpjson::get_string(fluid_json, "description"));
	    fluid.setReference(cpjson::get_string(fluid_json, "reference"));
	    fluid.setTmax(parse_value(fluid_json, "Tmax", true, 0.0));
	    fluid.setTmin(parse_value(fluid_json, "Tmin", true, 0.0));
	    fluid.setxmax(parse_value(fluid_json, "xmax", false, 1.0));
	    fluid.setxmin(parse_value(fluid_json, "xmin", false, 0.0));
	    fluid.setTminPsat(parse_value(fluid_json, "TminPsat", false, 0.0));

	    fluid.setTbase(parse_value(fluid_json, "Tbase", false, 0.0));
	    fluid.setxbase(parse_value(fluid_json, "xbase", false, 0.0));

	    /// Setters for the coefficients
	    fluid.setDensity(parse_coefficients(fluid_json, "density", true));
	    fluid.setSpecificHeat(parse_coefficients(fluid_json, "specific_heat", true));
	    fluid.setViscosity(parse_coefficients(fluid_json, "viscosity", false));
	    fluid.setConductivity(parse_coefficients(fluid_json, "conductivity", false));
	    fluid.setPsat(parse_coefficients(fluid_json, "saturation_pressure", false));
	    fluid.setTfreeze(parse_coefficients(fluid_json, "T_freeze", false));
	    fluid.setVolToMass(parse_coefficients(fluid_json, "volume2mass", false));
	    fluid.setMassToMole(parse_coefficients(fluid_json, "mass2mole", false));

	    fluid.set_reference_state(
			    parse_value(fluid_json, "Tref", false, 25+273.15) ,
			    parse_value(fluid_json, "pref", false, 1.01325e5) ,
			    parse_value(fluid_json, "xref", false, 0.0) ,
			    parse_value(fluid_json, "href", false, 0.0) ,
			    parse_value(fluid_json, "sref", false, 0.0)
			    );

	    /// A function to check coefficients and equation types.
	    /// \todo Implement the validation function
	    //fluid.validate();

	    // Add name->index mapping
	    string_to_index_map[fluid.getName()] = index;
    }
    catch(std::exception &e)
    {
        std::cout << format("Unable to load fluid: %s\n", fluid.getName().c_str());
        throw;
    }

};

/// Get an IncompressibleFluid instance stored in this library
/**
 @param name Name of the fluid
 */
IncompressibleFluid& JSONIncompressibleLibrary::get(std::string key) {
	std::map<std::string, std::size_t>::iterator it;
	// Try to find it
	it = string_to_index_map.find(key);
	// If it is found
	if (it != string_to_index_map.end()) {
		return get(it->second);
	} else {
		throw ValueError(
			format(
					"key [%s] was not found in string_to_index_map in JSONIncompressibleLibrary",
					key.c_str()
			)
		);
	}
};

/// Get a IncompressibleFluid instance stored in this library
/**
 @param key The index of the fluid in the map
 */
IncompressibleFluid& JSONIncompressibleLibrary::get(std::size_t key) {
	std::map<std::size_t, IncompressibleFluid>::iterator it;
	// Try to find it
	it = fluid_map.find(key);
	// If it is found
	if (it != fluid_map.end()) {
		return it->second;
	} else {
		throw ValueError(
			format("key [%d] was not found in JSONIncompressibleLibrary",key));
	}
};

























static JSONIncompressibleLibrary library;

void load_incompressible_library()
{
	rapidjson::Document dd;
    // This json formatted string comes from the all_incompressibles_JSON.h header which is a C++-escaped version of the JSON file
    dd.Parse<0>(all_incompressibles_JSON.c_str());
	if (dd.HasParseError()){
        throw ValueError("Unable to load all_incompressibles_JSON.json");
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
