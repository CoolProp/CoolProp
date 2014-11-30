#include "CoolProp.h"
#include "AbstractState.h"

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#if defined(__ISWINDOWS__)
#include <windows.h>
#else
#ifndef DBL_EPSILON
    #include <limits>
    #define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif
#endif

#include <memory>

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <exception>
#include <stdio.h>
#include <string>
#include "CoolPropTools.h"
#include "Solvers.h"
#include "MatrixMath.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "Backends/Incompressible/IncompressibleLibrary.h"
#include "Backends/Incompressible/IncompressibleBackend.h"
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Helmholtz/MixtureParameters.h"
#include "DataStructures.h"

#if defined(ENABLE_CATCH)
    #include "catch.hpp"
#endif

namespace CoolProp
{

static int debug_level = 0;
static std::string error_string;
static std::string warning_string;

void set_debug_level(int level){debug_level = level;}
int get_debug_level(void){return debug_level;}

//// This is very hacky, but pull the git revision from the file
#include "gitrevision.h" // Contents are like "std::string gitrevision = "aa121435436ggregrea4t43t433";"
#include "cpversion.h" // Contents are like "char version [] = "2.5";"

void set_warning_string(std::string warning){
    warning_string = warning;
}
void set_error_string(std::string error){
    error_string = error;
}

// Return true if the string has "BACKEND::*" format where * signifies a wildcard
bool has_backend_in_string(const std::string &fluid_string, std::size_t &i)
{
    i = fluid_string.find("::");
    return i != std::string::npos;
}

void extract_backend(const std::string &fluid_string, std::string &backend, std::string &fluid)
{
    std::size_t i;
    std::string _fluid_string = fluid_string;
    // For backwards compatibility reasons, if "REFPROP-" or "REFPROP-MIX:" start 
    // the fluid_string, replace them with "REFPROP::"
    if (_fluid_string.find("REFPROP-MIX:") == 0)
    {
        _fluid_string.replace(0, 12, "REFPROP::");
    }
    if (_fluid_string.find("REFPROP-") == 0)
    {
        _fluid_string.replace(0, 8, "REFPROP::");
    }
    if (has_backend_in_string(_fluid_string, i))
    {
        // Part without the ::
        backend = _fluid_string.substr(0, i);
        // Fluid name after the ::
        fluid = _fluid_string.substr(i+2);
    }
    else
    {
        backend = "?";
        fluid = _fluid_string;
    }
    if (get_debug_level()>10) std::cout << format("%s:%d: backend extracted. backend: %s. fluid: %s\n",__FILE__,__LINE__, backend.c_str(), fluid.c_str());
}

bool has_fractions_in_string(const std::string &fluid_string)
{
    // If can find both "[" and "]", it must have mole fractions encoded as string
    return (fluid_string.find("[")!=std::string::npos && fluid_string.find("]")!=std::string::npos);
}
bool has_solution_concentration(const std::string &fluid_string)
{
    // If can find "-", expect mass fractions encoded as string
    return (fluid_string.find('-') != std::string::npos && fluid_string.find('%') != std::string::npos);
}

std::string extract_fractions(const std::string &fluid_string, std::vector<double> &fractions)
{
    if (has_fractions_in_string(fluid_string))
    {
        fractions.clear();
        std::vector<std::string> names;

        // Break up into pairs - like "Ethane[0.5]&Methane[0.5]" -> ("Ethane[0.5]","Methane[0.5]")
        std::vector<std::string> pairs = strsplit(fluid_string, '&');

        for (std::size_t i = 0; i < pairs.size(); ++i)
        {
            std::string fluid = pairs[i];

            // Must end with ']'
            if (fluid[fluid.size()-1] != ']')
                throw ValueError(format("Fluid entry [%s] must end with ']' character",pairs[i].c_str()));

            // Split at '[', but first remove the ']' from the end by taking a substring
            std::vector<std::string> name_fraction = strsplit(fluid.substr(0, fluid.size()-1), '[');

            if (name_fraction.size() != 2){throw ValueError(format("Could not break [%s] into name/fraction", fluid.substr(0, fluid.size()-1).c_str()));}

            // Convert fraction to a double
            char *pEnd;
            std::string &name = name_fraction[0], &fraction = name_fraction[1];
            double f = strtod(fraction.c_str(), &pEnd);

            // If pEnd points to the last character in the string, it wasn't able to do the conversion
            if (pEnd == &(fraction[fraction.size()-1])){throw ValueError(format("Could not convert [%s] into number", fraction.c_str()));}

            // And add to vector
            fractions.push_back(f);

            // Add name
            names.push_back(name);
        }

        if (get_debug_level()>10) std::cout << format("%s:%d: Detected fractions of %s for %s.",__FILE__,__LINE__,vec_to_string(fractions).c_str(), (strjoin(names, "&")).c_str());
        // Join fluids back together
        return strjoin(names, "&");
    }
    else if (has_solution_concentration(fluid_string))
    {
        fractions.clear();
        double x;

        std::vector<std::string> fluid_parts = strsplit(fluid_string,'-');
        // Check it worked
        if (fluid_parts.size() != 2){
            throw ValueError(format("Format of incompressible solution string [%s] is invalid, should be like \"EG-20%\" or \"EG-0.2\" ", fluid_string.c_str()) );
        }

        // Convert the concentration into a string
        char* pEnd;
        x = strtod(fluid_parts[1].c_str(), &pEnd);

        // Check if per cent or fraction syntax is used
        if (!strcmp(pEnd,"%")){    x *= 0.01;}
        fractions.push_back(x);
        if (get_debug_level()>10) std::cout << format("%s:%d: Detected incompressible concentration of %s for %s.",__FILE__,__LINE__,vec_to_string(fractions).c_str(), fluid_parts[0].c_str());
        return fluid_parts[0];
    }
    else
    { 
        return fluid_string;
    }
}

// Internal function to do the actual calculations, make this a wrapped function so
// that error bubbling can be done properly
double _PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &backend, const std::string &Ref, const std::vector<double> &z)
{
    parameters iOutput = iundefined_parameter;
    parameters iOf = iundefined_parameter, iWrt = iundefined_parameter, iConstant = iundefined_parameter, 
               iOf1 = iundefined_parameter, iWrt1 = iundefined_parameter, iConstant1 = iundefined_parameter, 
               iWrt2 = iundefined_parameter, iConstant2 = iundefined_parameter;
    double x1, x2;
    
    if (get_debug_level()>5){
        std::cout << format("%s:%d: _PropsSI(%s,%s,%g,%s,%g,%s,%s)\n",__FILE__,__LINE__,Output.c_str(),Name1.c_str(),Prop1,Name2.c_str(),Prop2,backend.c_str(),Ref.c_str(), vec_to_string(z).c_str()).c_str();
    }

    // The state we are going to use
    shared_ptr<AbstractState> State;
    
    // If the fractions of the components have been encoded in the string, extract them
    // If they have not, this function does nothing
    std::vector<double> fractions;
    if (z.empty())
    {
        // Make a one-element vector
        fractions = std::vector<double>(1, 1);
    }
    else{
        // Make a copy
        fractions = z;
    }
    
    std::string fluid_string = extract_fractions(Ref, fractions);
    
    // We are going to let the factory function load the state
    State.reset(AbstractState::factory(backend, fluid_string));
    
    // Set the fraction for the state
    if (State->using_mole_fractions()){
        State->set_mole_fractions(fractions);
    } else if (State->using_mass_fractions()){
        State->set_mass_fractions(fractions);
    } else if (State->using_volu_fractions()){
        State->set_volu_fractions(fractions);
    } else {
        if (get_debug_level()>50) std::cout << format("%s:%d: _PropsSI, could not set composition to %s, defaulting to mole fraction.\n",__FILE__,__LINE__, vec_to_string(z).c_str()).c_str();
    }
    
    // First check if it is a trivial input (critical/max parameters for instance)
    if (is_valid_parameter(Output, iOutput))
    {
        if (is_trivial_parameter(iOutput)){
            double val = State->trivial_keyed_output(iOutput);
            return val;
        }
    }
    
    parameters iName1 = get_parameter_index(Name1);
    parameters iName2 = get_parameter_index(Name2);

    // Obtain the input pair
    CoolProp::input_pairs pair = generate_update_pair(iName1, Prop1, iName2, Prop2, x1, x2);

    // Update the state
    State->update(pair, x1, x2);
    
    if (iOutput != iundefined_parameter){
        // Get the desired output
        double val = State->keyed_output(iOutput);
        
        // Return the value
        return val;
    }
    else if (is_valid_first_derivative(Output, iOf, iWrt, iConstant)){
        // Return the desired output
        double val = State->first_partial_deriv(iOf, iWrt, iConstant);
        
        // Return the value
        return val;
    }
    else if (is_valid_second_derivative(Output, iOf1, iWrt1, iConstant1, iWrt2, iConstant2)){
        // Return the desired output
        double val = State->second_partial_deriv(iOf1, iWrt1, iConstant1, iWrt2, iConstant2);
        
        // Return the value
        return val;
    }
    else{
        throw ValueError(format("Output [%s] is not a parameter or a string representation of a derivative",Output.c_str()).c_str());
    }
}
double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref, const std::vector<double> &z)
{
    std::string backend, fluid;
    // Fractions are already provided, we just need to parse the Ref string
    extract_backend(Ref, backend, fluid);
    CATCH_ALL_ERRORS_RETURN_HUGE(return _PropsSI(Output,Name1,Prop1,Name2,Prop2,backend, fluid,z);)
}
double PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *FluidName, const std::vector<double> &x)
{
    std::string _Output = Output, _Name1 = Name1, _Name2 = Name2, _FluidName = FluidName;
    return PropsSI(_Output,_Name1,Prop1,_Name2,Prop2,_FluidName, x);
}
double Props1SI(const std::string &FluidName, const std::string &Output)
{
    std::string _FluidName = FluidName, empty_string = "", _Output = Output;
    bool valid_fluid1 = is_valid_fluid_string(_FluidName);
    bool valid_fluid2 = is_valid_fluid_string(_Output);
    if (valid_fluid1 && valid_fluid2){
        set_error_string(format("Both inputs to Props1SI [%s,%s] are valid fluids", Output.c_str(), FluidName.c_str()));
        return _HUGE;
    }
    if (!valid_fluid1 && !valid_fluid2){
        set_error_string(format("Neither input to Props1SI [%s,%s] is a valid fluid", Output.c_str(), FluidName.c_str()));
        return _HUGE;
    }
    if (!valid_fluid1 && valid_fluid2){
        // They are backwards, swap
        std::swap(_Output, _FluidName);
    }
    
    // First input is the fluid, second input is the input parameter
    double val1 = PropsSI(_Output, "", 0, "", 0, _FluidName);
    if (!ValidNumber(val1)){
        set_error_string(format("Unable to use input parameter [%s] in Props1SI for fluid %s; error was %s", _Output.c_str(), _FluidName.c_str(), get_global_param_string("errstring").c_str()));
        return _HUGE;
    }
    else{
        return val1;
    }
}
double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref)
{
    std::string backend, fluid;
    #if !defined(PROPSSI_NO_ERROR_CATCH)
    // In this function the error catching happens;
    try{
    #else
    std::cout << "macro is on; error checking disabled in PropsSI" << std::endl;
    #endif
        // BEGIN OF TRY
        // Here is the real code that is inside the try block
        extract_backend(Ref, backend, fluid);
        double val = _PropsSI(Output, Name1, Prop1, Name2, Prop2, backend, fluid, std::vector<double>());
        if (get_debug_level() > 1){ std::cout << format("_PropsSI will return %g",val) << std::endl; }
        return val;
        // END OF TRY
    #if !defined(PROPSSI_NO_ERROR_CATCH)
    }
    catch(const std::exception& e){
        set_error_string(e.what() + format(" : PropsSI(\"%s\",\"%s\",%0.10g,\"%s\",%0.10g,\"%s\")",Output.c_str(),Name1.c_str(), Prop1, Name2.c_str(), Prop2, Ref.c_str())); 
        #if defined (PROPSSI_ERROR_STDOUT)
        std::cout << e.what() << std::endl; 
        #endif
        if (get_debug_level() > 1){std::cout << e.what() << std::endl;}
        return _HUGE; 
    }
    catch(...){
        return _HUGE; 
    }
    #endif
}
std::vector<double> PropsSI(const std::string &Output, const std::string &Name1, const std::vector<double> &Prop1, const std::string &Name2, const std::vector<double> Prop2, const std::string &FluidName)
{
    return PropsSI(Output, Name1, Prop1, Name2, Prop2, FluidName, std::vector<double>(1,1));
}
std::vector<double> PropsSI(const std::string &Output, const std::string &Name1, const std::vector<double> &Prop1, const std::string &Name2, const std::vector<double> Prop2, const std::string &Ref, const std::vector<double> &z)
{
    std::vector<double> out(Prop1.size(), _HUGE);
    if (Prop1.size() != Prop2.size())
    {
        throw ValueError(format("Sizes of Prop1 [%d] and Prop2 [%d] to PropsSI are not the same", Prop1.size(), Prop2.size()));
    }
    for (std::size_t i = 0; i < Prop1.size(); ++i)
    {
        out[i] = PropsSI(Output,Name1,Prop1[i],Name2,Prop2[i],Ref,z);
    }
    return out;
}
double Props1SI(std::string FluidName,std::string Output)
{
    return PropsSI(Output,"",0,"",0,FluidName);
}

bool is_valid_fluid_string(std::string &input_fluid_string)
{
    try{
        std::string backend, fluid;
        std::vector<double> fractions;
        // First try to extract backend and fractions
        extract_backend(input_fluid_string, backend, fluid);
        std::string fluid_string = extract_fractions(fluid, fractions);
        // We are going to let the factory function load the state
        shared_ptr<AbstractState> State(AbstractState::factory(backend, fluid_string));
        return true;
    }
    catch (std::exception &e){
        return false;
    }
}
double saturation_ancillary(const std::string &fluid_name, const std::string &output, int Q, const std::string &input, double value){
    
    // Generate the state instance
    std::vector<std::string> names(1, fluid_name);
    shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(names));
    
    parameters iInput = get_parameter_index(input);
    parameters iOutput = get_parameter_index(output);
    
    return HEOS->saturation_ancillary(iOutput, Q, iInput, value);
}
void set_reference_stateS(std::string Ref, std::string reference_state)
{
    shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS;
    std::vector<std::string> _comps(1, Ref);
    HEOS.reset(new CoolProp::HelmholtzEOSMixtureBackend(_comps));

    if (!reference_state.compare("IIR"))
    {
        HEOS->update(QT_INPUTS, 0, 273.15);

        // Get current values for the enthalpy and entropy
        double deltah = HEOS->hmass() - 200000; // offset from 200000 J/kg enthalpy
        double deltas = HEOS->smass() - 1000; // offset from 1000 J/kg/K entropy
        double delta_a1 = deltas/(HEOS->gas_constant()/HEOS->molar_mass());
        double delta_a2 = -deltah/(HEOS->gas_constant()/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "IIR");
        HEOS->update_states();
    }
    else if (!reference_state.compare("ASHRAE"))
    {
        HEOS->update(QT_INPUTS, 0, 243.15);

        // Get current values for the enthalpy and entropy
        double deltah = HEOS->hmass() - 0; // offset from 0 J/kg enthalpy
        double deltas = HEOS->smass() - 0; // offset from 0 J/kg/K entropy
        double delta_a1 = deltas/(HEOS->gas_constant()/HEOS->molar_mass());
        double delta_a2 = -deltah/(HEOS->gas_constant()/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "ASHRAE");
        HEOS->update_states();
    }
    else if (!reference_state.compare("NBP"))
    {
        // Saturated liquid boiling point at 1 atmosphere
        HEOS->update(PQ_INPUTS, 101325, 0);

        double deltah = HEOS->hmass() - 0; // offset from 0 kJ/kg enthalpy
        double deltas = HEOS->smass() - 0; // offset from 0 kJ/kg/K entropy
        double delta_a1 = deltas/(HEOS->gas_constant()/HEOS->molar_mass());
        double delta_a2 = -deltah/(HEOS->gas_constant()/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        if (get_debug_level() > 5){std::cout << format("[set_reference_stateD] delta_a1 %g delta_a2 %g\n",delta_a1, delta_a2);}
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "NBP");
        HEOS->update_states();
    }
    else if (!reference_state.compare("DEF"))
    {
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(0,0,"");
    }
    else if (!reference_state.compare("RESET"))
    {
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(0, 0, "");
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffsetCore.set(0, 0, "");
    }
    else
    {
        throw ValueError(format("reference state string is invalid: [%s]",reference_state.c_str()));
    }
}
void set_reference_stateD(std::string Ref, double T, double rhomolar, double h0, double s0)
{
    shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS;
    std::vector<std::string> _comps(1, Ref);
    HEOS.reset(new CoolProp::HelmholtzEOSMixtureBackend(_comps));
    
    HEOS->update(DmolarT_INPUTS, rhomolar, T);

    // Get current values for the enthalpy and entropy
    double deltah = HEOS->hmass() - h0; // offset from specified enthalpy in J/mol
    double deltas = HEOS->smass() - s0; // offset from specified entropy in J/mol/K
    double delta_a1 = deltas/(8.314472/HEOS->molar_mass());
    double delta_a2 = -deltah/(8.314472/HEOS->molar_mass()*HEOS->get_reducing_state().T);
    HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "custom");
    HEOS->update_states();
}

std::string get_BibTeXKey(std::string Ref, std::string key)
{
    std::vector<std::string> names(1, Ref);
    HelmholtzEOSMixtureBackend HEOS(names);

    if (!key.compare("EOS")){ return HEOS.get_components()[0]->pEOS->BibTeX_EOS; }
    else if (!key.compare("CP0")){ return HEOS.get_components()[0]->pEOS->BibTeX_CP0; }
    else if (!key.compare("VISCOSITY")){ return HEOS.get_components()[0]->transport.BibTeX_viscosity; }
    else if (!key.compare("CONDUCTIVITY")){ return HEOS.get_components()[0]->transport.BibTeX_conductivity; }
    else if (!key.compare("ECS_LENNARD_JONES")){ throw NotImplementedError(); }
    else if (!key.compare("ECS_VISCOSITY_FITS")){ throw NotImplementedError(); }
    else if (!key.compare("ECS_CONDUCTIVITY_FITS")){ throw NotImplementedError(); }
    else if (!key.compare("SURFACE_TENSION")){ return HEOS.get_components()[0]->ancillaries.surface_tension.BibTeX;}
    else if (!key.compare("MELTING_LINE")){ return HEOS.get_components()[0]->ancillaries.melting_line.BibTeX;}
    else{ return "Bad key";}
}
std::string get_global_param_string(std::string ParamName)
{
    if (!ParamName.compare("version")){ return version; }
    else if (!ParamName.compare("gitrevision")){
        return gitrevision;
    }
    else if (!ParamName.compare("errstring")){
        std::string temp = error_string; error_string = ""; return temp;
    }
    else if (!ParamName.compare("warnstring")){
        std::string temp = warning_string; warning_string = ""; return temp;
    }
    else if (!ParamName.compare("FluidsList") || !ParamName.compare("fluids_list") || !ParamName.compare("fluidslist")){
        return get_fluid_list();
    }
    else if (!ParamName.compare("incompressible_list_pure")){
        return get_incompressible_list_pure();
    }
    else if (!ParamName.compare("incompressible_list_solution")){
        return get_incompressible_list_solution();
    }
    else if (!ParamName.compare("mixture_binary_pairs_list")){
        return get_csv_mixture_binary_pairs();
    }
    else if (!ParamName.compare("parameter_list") ){
        return get_csv_parameter_list();
    }
    else{
        throw ValueError(format("Input value [%s] is invalid",ParamName.c_str()));
    }
};
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to get_global_param_string","[get_global_param_string]")
{
    const int num_good_inputs = 7;
    std::string good_inputs[num_good_inputs] = {"version", "gitrevision", "fluids_list", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list","parameter_list"};
    std::ostringstream ss3c;
    for (int i = 0; i<num_good_inputs; ++i){
        ss3c << "Test for" << good_inputs[i];
        SECTION(ss3c.str(), ""){          
            CHECK_NOTHROW(CoolProp::get_global_param_string(good_inputs[i]));  
        };
    }
    CHECK_THROWS(CoolProp::get_global_param_string(""));
};
#endif

std::string get_fluid_param_string(std::string FluidName, std::string ParamName)
{
    std::string backend, fluid;
    extract_backend(FluidName, backend, fluid);
    if (backend == "INCOMP"){
        try{
            shared_ptr<CoolProp::IncompressibleBackend> INCOMP(new CoolProp::IncompressibleBackend(fluid));
                    
            if (!ParamName.compare("long_name")){
                return INCOMP->calc_name();
            }
            else{
                throw ValueError(format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()));
            }
        }
        catch(std::exception &e){ throw ValueError(format("CoolProp error: %s", e.what())); }
        catch(...){ throw ValueError("CoolProp error: Indeterminate error"); }
    }
    
    try{
        std::vector<std::string> comps(1, FluidName);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(comps));
        CoolProp::CoolPropFluid *fluid = HEOS->get_components()[0];
                
        if (!ParamName.compare("aliases")){
            return strjoin(fluid->aliases, ", ");
        }
        else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number")){
            return fluid->CAS;
        }
        else if (!ParamName.compare("ASHRAE34")){
            return fluid->environment.ASHRAE34;
        }
        else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname")){
            return fluid->REFPROPname;
        }
        else if (ParamName.find("BibTeX") == 0) // Starts with "BibTeX"
        {
            std::vector<std::string> parts = strsplit(ParamName,'-');
            if (parts.size() != 2){ throw ValueError(format("Unable to parse BibTeX string %s",ParamName.c_str()));}
            return get_BibTeXKey( FluidName, parts[1]);
        }
        else{
            throw ValueError(format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()));
        }
    }
    catch(std::exception &e){ throw ValueError(format("CoolProp error: %s", e.what())); }
    catch(...){ throw ValueError("CoolProp error: Indeterminate error"); }
}
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to get_fluid_param_string", "[get_fluid_param_string]")
{
    const int num_good_inputs = 10;
    std::string good_inputs[num_good_inputs] = {"aliases", "CAS", "ASHRAE34", "REFPROPName", "BibTeX-CONDUCTIVITY", "BibTeX-EOS", "BibTeX-CP0", "BibTeX-SURFACE_TENSION","BibTeX-MELTING_LINE","BibTeX-VISCOSITY"};
    std::ostringstream ss3c;
    for (int i = 0; i < num_good_inputs; ++i){
        ss3c << "Test for" << good_inputs[i];
        SECTION(ss3c.str(), ""){          
            CHECK_NOTHROW(CoolProp::get_fluid_param_string("Water", good_inputs[i]));
        };
    }
    CHECK_THROWS(CoolProp::get_fluid_param_string("","aliases"));
    CHECK_THROWS(CoolProp::get_fluid_param_string("Water",""));
    CHECK_THROWS(CoolProp::get_fluid_param_string("Water","BibTeX-"));
};
#endif

std::string phase_lookup_string(phases Phase)
{
    switch (Phase)
    {
        case iphase_liquid: ///< Liquid
            return "liquid";
        case iphase_supercritical: ///< Supercritical (p > pc, T > Tc)
            return "supercritical";
        case iphase_supercritical_gas: ///< Supercritical gas (p < pc, T > Tc)
            return "supercritical_gas";
        case iphase_supercritical_liquid: ///< Supercritical liquid (p > pc, T < Tc)
            return "supercritical_liquid";
        case iphase_critical_point: ///< At the critical point
            return "critical_point";
        case iphase_gas: ///< Subcritical gas
            return "gas";
        case iphase_twophase: ///< Twophase
            return "twophase";
        case iphase_unknown: ///< Unknown phase
            return "unknown";
        case iphase_not_imposed:
            return "not_imposed";
    }
    throw ValueError("I should never be thrown");
}
std::string PhaseSI(const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName)
{
    double Phase_double = PropsSI("Phase",Name1,Prop1,Name2,Prop2,FluidName);
    if (!ValidNumber(Phase_double)){ return "";}
    std::size_t Phase_int = static_cast<std::size_t>(Phase_double);
    return phase_lookup_string(static_cast<phases>(Phase_int));
}
    
std::string PhaseSI(const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName, const std::vector<double> &z)
{
    double Phase_double = PropsSI("Phase",Name1,Prop1,Name2,Prop2,FluidName,z);
    if (!ValidNumber(Phase_double)){ return "";}
    std::size_t Phase_int = static_cast<std::size_t>(Phase_double);
    return phase_lookup_string(static_cast<phases>(Phase_int));
}

} /* namespace CoolProp */