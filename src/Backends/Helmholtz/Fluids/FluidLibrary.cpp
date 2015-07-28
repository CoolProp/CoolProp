
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