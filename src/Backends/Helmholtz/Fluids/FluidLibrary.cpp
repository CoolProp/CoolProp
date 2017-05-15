
#include "FluidLibrary.h"
#include "all_fluids_JSON.h" // Makes a std::string variable called all_fluids_JSON
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"

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

void JSONFluidLibrary::set_fluid_enthalpy_entropy_offset(const std::string &fluid, double delta_a1, double delta_a2, const std::string &ref)
{
    // Try to find it
    std::map<std::string, std::size_t>::const_iterator it = string_to_index_map.find(fluid);
    if (it != string_to_index_map.end()){
        std::map<std::size_t, CoolPropFluid>::iterator it2 = fluid_map.find(it->second);
        // If it is found
        if (it2 != fluid_map.end()){
            if (!ValidNumber(delta_a1) || !ValidNumber(delta_a2) ){
                throw ValueError(format("Not possible to set reference state for fluid %s because offset values are NAN",fluid.c_str()));
            }
            it2->second.EOS().alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, ref);
            
            shared_ptr<CoolProp::HelmholtzEOSBackend> HEOS(new CoolProp::HelmholtzEOSBackend(it2->second));
            HEOS->specify_phase(iphase_gas); // Something homogeneous;
            // Calculate the new enthalpy and entropy values
            HEOS->update(DmolarT_INPUTS, it2->second.EOS().hs_anchor.rhomolar, it2->second.EOS().hs_anchor.T);
            it2->second.EOS().hs_anchor.hmolar = HEOS->hmolar();
            it2->second.EOS().hs_anchor.smolar = HEOS->smolar();
            
            double f = (HEOS->name() == "Water" || HEOS->name() == "CarbonDioxide") ? 1.00001 : 1.0;

            // Calculate the new enthalpy and entropy values at the reducing state
            HEOS->update(DmolarT_INPUTS, it2->second.EOS().reduce.rhomolar*f, it2->second.EOS().reduce.T*f);
            it2->second.EOS().reduce.hmolar = HEOS->hmolar();
            it2->second.EOS().reduce.smolar = HEOS->smolar();

            // Calculate the new enthalpy and entropy values at the critical state
            HEOS->update(DmolarT_INPUTS, it2->second.crit.rhomolar*f, it2->second.crit.T*f);
            it2->second.crit.hmolar = HEOS->hmolar();
            it2->second.crit.smolar = HEOS->smolar();

            // Calculate the new enthalpy and entropy values
            HEOS->update(DmolarT_INPUTS, it2->second.triple_liquid.rhomolar, it2->second.triple_liquid.T);
            it2->second.triple_liquid.hmolar = HEOS->hmolar();
            it2->second.triple_liquid.smolar = HEOS->smolar();

            // Calculate the new enthalpy and entropy values
            HEOS->update(DmolarT_INPUTS, it2->second.triple_vapor.rhomolar, it2->second.triple_vapor.T);
            it2->second.triple_vapor.hmolar = HEOS->hmolar();
            it2->second.triple_vapor.smolar = HEOS->smolar();

            if (!HEOS->is_pure()){
                // Calculate the new enthalpy and entropy values
                HEOS->update(DmolarT_INPUTS, it2->second.EOS().max_sat_T.rhomolar, it2->second.EOS().max_sat_T.T);
                it2->second.EOS().max_sat_T.hmolar = HEOS->hmolar();
                it2->second.EOS().max_sat_T.smolar = HEOS->smolar();
                // Calculate the new enthalpy and entropy values
                HEOS->update(DmolarT_INPUTS, it2->second.EOS().max_sat_p.rhomolar, it2->second.EOS().max_sat_p.T);
                it2->second.EOS().max_sat_p.hmolar = HEOS->hmolar();
                it2->second.EOS().max_sat_p.smolar = HEOS->smolar();
            }
        }
        else{
            throw ValueError(format("fluid [%s] was not found in JSONFluidLibrary",fluid.c_str()));
        }
    }
}
    
/// Add all the fluid entries in the JSON-encoded string passed in
void JSONFluidLibrary::add_many(const std::string &JSON_string){
    
    // First load all the baseline fluids
    if (library.is_empty()){ load(); }
    
    // Then, load the fluids we would like to add
    rapidjson::Document doc;
    cpjson::JSON_string_to_rapidjson(JSON_string, doc);
    library.add_many(doc);
};
    
void JSONFluidLibrary::add_many(rapidjson::Value &listing)
{
    for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr)
    {
        add_one(*itr);
    }
};
    
void JSONFluidLibrary::add_one(rapidjson::Value &fluid_json)
{
    _is_empty = false;
    
    std::size_t index = fluid_map.size();
    
    CoolPropFluid fluid;
    
    // Fluid name
    fluid.name = fluid_json["INFO"]["NAME"].GetString(); name_vector.push_back(fluid.name);
    
    try{
        // CAS number
        if (!fluid_json["INFO"].HasMember("CAS")){ throw ValueError(format("fluid [%s] does not have \"CAS\" member",fluid.name.c_str())); }
        fluid.CAS = fluid_json["INFO"]["CAS"].GetString();
        
        // REFPROP alias
        if (!fluid_json["INFO"].HasMember("REFPROP_NAME")){ throw ValueError(format("fluid [%s] does not have \"REFPROP_NAME\" member",fluid.name.c_str())); }
        fluid.REFPROPname = fluid_json["INFO"]["REFPROP_NAME"].GetString();
        
        if (fluid_json["INFO"].HasMember("FORMULA")){
            fluid.formula = cpjson::get_string(fluid_json["INFO"], "FORMULA");
        }
        else{ fluid.formula = "N/A"; }
        
        if (fluid_json["INFO"].HasMember("INCHI_STRING")){
            fluid.InChI = cpjson::get_string(fluid_json["INFO"], "INCHI_STRING");
        }
        else{ fluid.InChI = "N/A"; }
        
        if (fluid_json["INFO"].HasMember("INCHI_KEY")){
            fluid.InChIKey = cpjson::get_string(fluid_json["INFO"], "INCHI_KEY");
        }
        else{ fluid.InChIKey = "N/A"; }
        
        if (fluid_json["INFO"].HasMember("SMILES")){
            fluid.smiles = cpjson::get_string(fluid_json["INFO"], "SMILES");
        }
        else{ fluid.smiles = "N/A"; }
        
        if (fluid_json["INFO"].HasMember("CHEMSPIDER_ID")){
            fluid.ChemSpider_id = cpjson::get_integer(fluid_json["INFO"], "CHEMSPIDER_ID");
        }
        else{ fluid.ChemSpider_id = -1; }
        
        if (fluid_json["INFO"].HasMember("2DPNG_URL")){
            fluid.TwoDPNG_URL = cpjson::get_string(fluid_json["INFO"], "2DPNG_URL");
        }
        else{ fluid.TwoDPNG_URL = "N/A"; }
        
        // Parse the environmental parameters
        if (!(fluid_json["INFO"].HasMember("ENVIRONMENTAL"))){
            if (get_debug_level() > 0){
                std::cout << format("Environmental data are missing for fluid [%s]\n", fluid.name.c_str()) ;
            }
        }
        else{
            parse_environmental(fluid_json["INFO"]["ENVIRONMENTAL"], fluid);
        }
        
        // Aliases
        fluid.aliases = cpjson::get_string_array(fluid_json["INFO"]["ALIASES"]);
        
        // Critical state
        if (!fluid_json.HasMember("STATES")){ throw ValueError(format("fluid [%s] does not have \"STATES\" member",fluid.name.c_str())); }
        parse_states(fluid_json["STATES"], fluid);
        
        if (get_debug_level() > 5){
            std::cout << format("Loading fluid %s with CAS %s; %d fluids loaded\n", fluid.name.c_str(), fluid.CAS.c_str(), index);
        }
        
        // EOS
        parse_EOS_listing(fluid_json["EOS"], fluid);
        
        // Validate the fluid
        validate(fluid);
        
        // Ancillaries for saturation
        if (!fluid_json.HasMember("ANCILLARIES")){throw ValueError(format("Ancillary curves are missing for fluid [%s]",fluid.name.c_str()));};
        parse_ancillaries(fluid_json["ANCILLARIES"],fluid);
        
        // Surface tension
        if (!(fluid_json["ANCILLARIES"].HasMember("surface_tension"))){
            if (get_debug_level() > 0){
                std::cout << format("Surface tension curves are missing for fluid [%s]\n", fluid.name.c_str()) ;
            }
        }
        else{
            parse_surface_tension(fluid_json["ANCILLARIES"]["surface_tension"], fluid);
        }
        
        // Melting line
        if (!(fluid_json["ANCILLARIES"].HasMember("melting_line"))){
            if (get_debug_level() > 0){
                std::cout << format("Melting line curves are missing for fluid [%s]\n", fluid.name.c_str()) ;
            }
        }
        else{
            parse_melting_line(fluid_json["ANCILLARIES"]["melting_line"], fluid);
        }
        
        // Parse the environmental parameters
        if (!(fluid_json.HasMember("TRANSPORT"))){
            default_transport(fluid);
        }
        else{
            parse_transport(fluid_json["TRANSPORT"], fluid);
        }
        
        // If the fluid is ok...
        
        // First check that none of the identifiers are already present
        bool already_present = false;
        
        if (string_to_index_map.find(fluid.CAS) != string_to_index_map.end()
            || string_to_index_map.find(fluid.name) != string_to_index_map.end()
            || string_to_index_map.find(upper(fluid.name)) != string_to_index_map.end()
            ){
            already_present = true;
        }
        else{
            // Check the aliases
            for (std::size_t i = 0; i < fluid.aliases.size(); ++i)
            {
                if (string_to_index_map.find(fluid.aliases[i]) != string_to_index_map.end()){ already_present = true; break; }
                if (string_to_index_map.find(upper(fluid.aliases[i])) != string_to_index_map.end()){ already_present = true; break; }
            }
        }
        
        if (already_present){
            if (!get_config_bool(OVERWRITE_FLUIDS)){
                throw ValueError(format("Cannot load fluid [%s:%s] because it is already in library; consider enabling the config boolean variable OVERWRITE_FLUIDS", fluid.name.c_str(), fluid.CAS.c_str()));
            }
            else{
                // Remove the one(s) that are already there
                
                // Remove the actual fluid instance
                std::size_t index = string_to_index_map.find(fluid.name)->second;
                
                if (JSONstring_map.find(index) != JSONstring_map.end()){
                    JSONstring_map.erase(JSONstring_map.find(index));
                }
                
                if (string_to_index_map.find(fluid.name) != string_to_index_map.end()){
                    fluid_map.erase(fluid_map.find(index));
                }
                
                // Remove the identifiers pointing to that instance
                if(string_to_index_map.find(fluid.CAS) != string_to_index_map.end()){
                    string_to_index_map.erase(string_to_index_map.find(fluid.CAS));
                }
                if(string_to_index_map.find(fluid.name) != string_to_index_map.end()){
                    string_to_index_map.erase(string_to_index_map.find(fluid.name));
                }
                // Check the aliases
                for (std::size_t i = 0; i < fluid.aliases.size(); ++i)
                {
                    if (string_to_index_map.find(fluid.aliases[i]) != string_to_index_map.end()){
                        string_to_index_map.erase(string_to_index_map.find(fluid.aliases[i]));
                    }
                    if (string_to_index_map.find(upper(fluid.aliases[i])) != string_to_index_map.end()){
                        string_to_index_map.erase(string_to_index_map.find(upper(fluid.aliases[i])));
                    }
                }
            }
        }
        
        // By now, the library has been cleared of remnants of this fluid; safe to add the fluid now.
        
        // Get the next index for this fluid
        index = fluid_map.size();
        
        // Add index->fluid mapping
        fluid_map[index] = fluid;
        
        JSONstring_map[index] = cpjson::json2string(fluid_json);
        
        // Add CAS->index mapping
        string_to_index_map[fluid.CAS] = index;
        
        // Add name->index mapping
        string_to_index_map[fluid.name] = index;
        
        // Add the aliases
        for (std::size_t i = 0; i < fluid.aliases.size(); ++i)
        {
            string_to_index_map[fluid.aliases[i]] = index;
            
            // Add uppercase alias for EES compatibility
            string_to_index_map[upper(fluid.aliases[i])] = index;
        }
        
        if (get_debug_level() > 5){ std::cout << format("Loaded.\n"); }
    }
    catch (const std::exception &e){
        throw ValueError(format("Unable to load fluid [%s] due to error: %s",fluid.name.c_str(),e.what()));
    }
};


JSONFluidLibrary & get_library(void){
    if (library.is_empty()){ load(); }
    return library;
}

CoolPropFluid get_fluid(const std::string &fluid_string){
    if (library.is_empty()){ load(); }
    return library.get(fluid_string);
}
    
std::string get_fluid_as_JSONstring(const std::string &identifier){
    if (library.is_empty()){ load(); }
    return library.get_JSONstring(identifier);
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