#include "CoolProp.h"
#include "AbstractState.h"

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>

#endif

#if defined(__ISWINDOWS__)
#include <windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
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

void _PropsSI_initialize(const std::string &backend,
                         const std::vector<std::string> &fluid_names,
                         const std::vector<double> &z,
                         shared_ptr<AbstractState> &State){

    if (fluid_names.empty()){throw ValueError("fluid_names cannot be empty");}

    std::vector<double> fractions(1, 1.0); // Default to one component, unity fraction
    const std::vector<double> *fractions_ptr = NULL; // Pointer to the array to be used;

    if (fluid_names.size() > 1){
        // Set the pointer - we are going to use the supplied fractions; they must be provided
        fractions_ptr = &z;
        // Reset the state
        State.reset(AbstractState::factory(backend, fluid_names));
    }
    else if (fluid_names.size() == 1){
        if (has_fractions_in_string(fluid_names[0]) || has_solution_concentration(fluid_names[0])){
            // Extract fractions from the string
            std::string fluid_string = extract_fractions(fluid_names[0], fractions);
            // Set the pointer - we are going to use the extracted fractions
            fractions_ptr = &fractions;
            // Reset the state
            State.reset(AbstractState::factory(backend, fluid_string));
        }
        else{
            if (z.empty()){
                // Set the pointer - we are going to use the default fractions
                fractions_ptr = &fractions;
            }
            else{
                // Set the pointer - we are going to use the provided fractions
                fractions_ptr = &z;
            }
            // Reset the state
            State.reset(AbstractState::factory(backend, fluid_names));
        }
    }
    else { // The only path where fractions_ptr stays NULL
      throw ValueError("fractions_ptr is NULL");
    }

    // Set the fraction for the state
    if (State->using_mole_fractions()){
        // If a predefined mixture or a pure fluid, the fractions will already be set
        if (State->get_mole_fractions().empty()){
            State->set_mole_fractions(*fractions_ptr);
        }
    } else if (State->using_mass_fractions()){
        State->set_mass_fractions(*fractions_ptr);
    } else if (State->using_volu_fractions()){
        State->set_volu_fractions(*fractions_ptr);
    } else {
        if (get_debug_level()>50) std::cout << format("%s:%d: _PropsSI, could not set composition to %s, defaulting to mole fraction.\n",__FILE__,__LINE__, vec_to_string(z).c_str()).c_str();
    }
}

struct output_parameter{
	enum OutputParametersType {OUTPUT_TYPE_UNSET = 0, OUTPUT_TYPE_TRIVIAL, OUTPUT_TYPE_NORMAL, OUTPUT_TYPE_FIRST_DERIVATIVE, OUTPUT_TYPE_SECOND_DERIVATIVE};
	CoolProp::parameters Of1, Wrt1, Constant1, Wrt2, Constant2;
	OutputParametersType type;
    /// Parse a '&' separated string into a data structure with one entry per output
    /// Covers both normal and derivative outputs
    static std::vector<output_parameter> get_output_parameters(const std::vector<std::string> &Outputs){
        std::vector<output_parameter> outputs;
        for (std::vector<std::string>::const_iterator str = Outputs.begin(); str != Outputs.end(); ++str){
            output_parameter out;
            CoolProp::parameters iOutput;
            if (is_valid_parameter(*str, iOutput)){
                out.Of1 = iOutput;
                if (is_trivial_parameter(iOutput)){ out.type = OUTPUT_TYPE_TRIVIAL; }
                else{ out.type = OUTPUT_TYPE_NORMAL; }
            }
            else if (is_valid_first_derivative(*str, out.Of1, out.Wrt1, out.Constant1)){
                out.type = OUTPUT_TYPE_FIRST_DERIVATIVE;
            }
            else if (is_valid_second_derivative(*str, out.Of1, out.Wrt1, out.Constant1, out.Wrt2, out.Constant2)){
                out.type = OUTPUT_TYPE_SECOND_DERIVATIVE;
            }
            else{
                throw ValueError(format("Output string is invalid [%s]", str->c_str()));
            }
            outputs.push_back(out);
        }
        return outputs;
    };
};

void _PropsSI_outputs(shared_ptr<AbstractState> &State,
	     			 std::vector<output_parameter> output_parameters,
		    		 CoolProp::input_pairs input_pair,
			    	 const std::vector<double> &in1,
			    	 const std::vector<double> &in2,
			    	 std::vector<std::vector<double> > &IO){

	// Check the inputs
	if (in1.size() != in2.size()){ throw ValueError(format("lengths of in1 [%d] and in2 [%d] are not the same", in1.size(), in2.size()));}
	bool one_input_one_output = (in1.size() == 1 && in2.size() == 1 && output_parameters.size() == 1);
    // If all trivial outputs, never do a state update
    bool all_trivial_outputs = true;
    for (std::size_t j = 0; j < output_parameters.size(); ++j){
        if (output_parameters[j].type != output_parameter::OUTPUT_TYPE_TRIVIAL){
            all_trivial_outputs = false;
        }
    }
    parameters p1, p2;
    // If all outputs are also inputs, never do a state update
    bool all_outputs_in_inputs = true;
    if (input_pair != INPUT_PAIR_INVALID){
        // Split the input pair into parameters
        split_input_pair(input_pair, p1, p2);
        // See if each parameter is in the output vector and is a normal type input
        for (std::size_t j = 0; j < output_parameters.size(); ++j){
            if (output_parameters[j].type != output_parameter::OUTPUT_TYPE_NORMAL){
                all_outputs_in_inputs = false; break;
            }
            if (!(output_parameters[j].Of1 == p1 || output_parameters[j].Of1 == p2)){
                all_outputs_in_inputs = false; break;
            }
        }
    }
    else{
        all_outputs_in_inputs = false;
    }

	if (get_debug_level() > 100)
	{
	   std::cout << format("%s (%d): input pair = %d ",__FILE__,__LINE__, input_pair)          << std::endl;
	   std::cout << format("%s (%d): in1 = %s ",__FILE__,__LINE__, vec_to_string(in1).c_str()) << std::endl;
	   std::cout << format("%s (%d): in2 = %s ",__FILE__,__LINE__, vec_to_string(in2).c_str()) << std::endl;
	}

	// Resize the output matrix
    std::size_t N1 = std::max(static_cast<std::size_t>(1), in1.size());
    std::size_t N2 = std::max(static_cast<std::size_t>(1), output_parameters.size());
	IO.resize(N1, std::vector<double>(N2, _HUGE));

    // Throw an error if at the end, there were no successes
    bool success = false;

    if (get_debug_level() > 100)
    {
        std::cout << format("%s (%d): Iterating over %d input value pairs.",__FILE__,__LINE__,IO.size()) << std::endl;
    }
	// Iterate over the state variable inputs
	for (std::size_t i = 0; i < IO.size(); ++i){
		try{
            if (input_pair != INPUT_PAIR_INVALID && !all_trivial_outputs && !all_outputs_in_inputs){
                // Update the state since it is a valid set of inputs
                State->update(input_pair, in1[i], in2[i]);
            }
        }
        catch(std::exception &){
            if (one_input_one_output){IO.clear(); throw;} // Re-raise the exception since we want to bubble the error
            // All the outputs are filled with _HUGE; go to next input
            for (std::size_t j = 0; j < IO[i].size(); ++j){ IO[i][j] = _HUGE; }
            continue;
        }

        for (std::size_t j = 0; j < IO[i].size(); ++j){
            // If all the outputs are inputs, there is no need for a state input
            if (all_outputs_in_inputs){
                if (p1 == output_parameters[j].Of1){
                    IO[i][j] = in1[i]; success = true; continue;
                }
                else if (p2 == output_parameters[j].Of1){
                    IO[i][j] = in2[i]; success = true; continue;
                }
                else{
                    throw ValueError();
                }
            }
            try{
                output_parameter &output = output_parameters[j];
                switch (output.type){
                    case output_parameter::OUTPUT_TYPE_TRIVIAL:
                    case output_parameter::OUTPUT_TYPE_NORMAL:
                        IO[i][j] = State->keyed_output(output.Of1); break;
                    case output_parameter::OUTPUT_TYPE_FIRST_DERIVATIVE:
                        IO[i][j] = State->first_partial_deriv(output.Of1, output.Wrt1, output.Constant1); break;
                    case output_parameter::OUTPUT_TYPE_SECOND_DERIVATIVE:
                        IO[i][j] = State->second_partial_deriv(output.Of1, output.Wrt1, output.Constant1, output.Wrt2, output.Constant2); break;
                    default:
                        throw ValueError(format("")); break;
                }
                // At least one has succeeded
                success = true;
            }
            catch(std::exception &){
                if (one_input_one_output){IO.clear(); throw;} // Re-raise the exception since we want to bubble the error
                IO[i][j] = _HUGE;
            }
        }
	}
    if (success == false) { IO.clear(); throw ValueError(format("No outputs were able to be calculated"));}
}

void _PropsSImulti(const std::vector<std::string> &Outputs,
                   const std::string &Name1,
                   const std::vector<double> &Prop1,
                   const std::string &Name2,
                   const std::vector<double> &Prop2,
                   const std::string &backend,
                   const std::vector<std::string> &fluids,
                   const std::vector<double> &fractions,
                   std::vector<std::vector<double> > &IO)
{
    shared_ptr<AbstractState> State;
    CoolProp::parameters key1, key2;
    CoolProp::input_pairs input_pair;
    std::vector<output_parameter> output_parameters;
    std::vector<double> v1, v2;

    try{
        // Initialize the State class
        _PropsSI_initialize(backend, fluids, fractions, State);
    }
    catch(std::exception &e){
        // Initialization failed.  Stop.
        throw ValueError(format("Initialize failed for backend: \"%s\", fluid: \"%s\" fractions \"%s\"; error: %s",backend.c_str(), strjoin(fluids,"&").c_str(), vec_to_string(fractions, "%0.10f").c_str(), e.what()) );
    }

    try{
        // Get update pair
        is_valid_parameter(Name1, key1);
        is_valid_parameter(Name2, key2);
        input_pair = generate_update_pair(key1, Prop1, key2, Prop2, v1, v2);
    }
    catch (std::exception &e){
        // Input parameter parsing failed.  Stop
        throw ValueError(format("Input pair parsing failed for Name1: \"%s\", Name2: \"%s\"; err: %s", Name1.c_str(), Name2.c_str(), e.what()));
    }

    try{
        output_parameters = output_parameter::get_output_parameters(Outputs);
    }
    catch (std::exception &e){
        // Output parameter parsing failed.  Stop.
        throw ValueError(format("Output parameter parsing failed; error: %s", e.what()));
    }

    // Calculate the output(s).  In the case of a failure, all values will be filled with _HUGE
    _PropsSI_outputs(State, output_parameters, input_pair, v1, v2, IO);
}

std::vector<std::vector<double> > PropsSImulti(const std::vector<std::string> &Outputs,
                                               const std::string &Name1,
                                               const std::vector<double> &Prop1,
                                               const std::string &Name2,
                                               const std::vector<double> &Prop2,
                                               const std::string &backend,
                                               const std::vector<std::string> &fluids,
                                               const std::vector<double> &fractions)
{
    std::vector<std::vector<double> > IO;

    #if !defined(NO_ERROR_CATCHING)
    try{
    #endif

        // Call the subfunction that can bubble errors
        _PropsSImulti(Outputs, Name1, Prop1, Name2, Prop2, backend, fluids, fractions, IO);

        // Return the value(s)
        return IO;

    #if !defined(NO_ERROR_CATCHING)
    }
    catch(const std::exception& e){
        set_error_string(e.what());
        #if defined (PROPSSI_ERROR_STDOUT)
        std::cout << e.what() << std::endl;
        #endif
        if (get_debug_level() > 1){std::cout << e.what() << std::endl;}
    }
    catch(...){
    }
    #endif
    return std::vector<std::vector<double> >();
}
double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref)
{
    std::string backend, fluid; std::vector<double> fractions(1,1.0);
    #if !defined(NO_ERROR_CATCHING)
    try{
    #endif

        // BEGIN OF TRY
        // Here is the real code that is inside the try block
    

        extract_backend(Ref, backend, fluid);
        std::string fluid_string = fluid;
        if (has_fractions_in_string(fluid) || has_solution_concentration(fluid)){
            fluid_string = extract_fractions(fluid, fractions);
        }
        std::vector<std::vector<double> > IO;
        _PropsSImulti(strsplit(Output,'&'), Name1, std::vector<double>(1, Prop1), Name2, std::vector<double>(1, Prop2), backend, strsplit(fluid_string, '&'), fractions, IO);
        if (IO.empty()){ throw ValueError(get_global_param_string("errstring").c_str()); }
        if (IO.size()!= 1 || IO[0].size() != 1){ throw ValueError(format("output should be 1x1; error was %s", get_global_param_string("errstring").c_str())); }

        double val = IO[0][0];

        if (get_debug_level() > 1){ std::cout << format("_PropsSI will return %g",val) << std::endl; }
        return val;
        // END OF TRY
    #if !defined(NO_ERROR_CATCHING)
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
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to PropsSI","[PropsSI]")
{
    SECTION("Single state, single output"){
        CHECK(ValidNumber(CoolProp::PropsSI("T","P",101325,"Q",0,"Water")));
    };
    SECTION("Single state, single output, pure incompressible"){
        CHECK(ValidNumber(CoolProp::PropsSI("D","P",101325,"T",300,"INCOMP::DowQ")));
    };
    SECTION("Single state, trivial output, pure incompressible"){
        CHECK(ValidNumber(CoolProp::PropsSI("Tmin","P",0,"T",0,"INCOMP::DowQ")));
    };
    SECTION("Bad input pair"){
        CHECK(!ValidNumber(CoolProp::PropsSI("D","Q",0,"Q",0,"Water")));
    };
    SECTION("Single state, single output, 40% incompressible"){
        CHECK(ValidNumber(CoolProp::PropsSI("D","P",101325,"T",300,"INCOMP::MEG[0.40]")));
    };
    SECTION("Single state, single output, predefined CoolProp mixture"){
        CHECK(ValidNumber(CoolProp::PropsSI("T","Q",1,"P",3e6,"HEOS::R125[0.7]&R32[0.3]")));
    };
    SECTION("Single state, single output"){
        CHECK(ValidNumber(CoolProp::PropsSI("T","P",101325,"Q",0,"HEOS::Water")));
    };
    SECTION("Single state, single output, predefined mixture"){
        CHECK(ValidNumber(CoolProp::PropsSI("T","P",101325,"Q",0,"R410A.mix")));
    };
    SECTION("Single state, single output, predefined mixture from REFPROP"){
        CHECK(ValidNumber(CoolProp::PropsSI("T","P",101325,"Q",0,"REFPROP::R410A.mix")));
    };
    SECTION("Single state, single output, bad predefined mixture from REFPROP"){
        CHECK(!ValidNumber(CoolProp::PropsSI("T","P",101325,"Q",0,"REFPROP::RRRRRR.mix")));
    };
    SECTION("Predefined mixture"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z;
        std::vector<std::string> outputs(1,"T"); outputs.push_back("Dmolar");
        std::vector<std::vector<double> > IO;
        std::vector<std::string> fluids(1, "R410A.mix");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
    };
    SECTION("Single state, two outputs"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1,"T"); outputs.push_back("Dmolar");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(std::vector<std::vector<double> > IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
    };
    SECTION("Single state, two bad outputs"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"???????"); outputs.push_back("?????????");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
        CHECK(IO.size() == 0);
    };
    SECTION("Two states, one output"){
        std::vector<double> p(2, 101325), Q(2, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1,"T");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(std::vector<std::vector<double> > IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
    };
    SECTION("Two states, two outputs"){
        std::vector<double> p(2, 101325), Q(2, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1,"T"); outputs.push_back("Dmolar");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(std::vector<std::vector<double> > IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
    };
    SECTION("cp and its derivative representation"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"Cpmolar"); outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(!IO.empty());
        CAPTURE(IO[0][0]);
        CAPTURE(IO[0][1]);
        CHECK(std::abs(IO[0][0] - IO[0][1]) < 1e-5);
    };
    SECTION("bad fluid"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"Cpmolar"); outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "????????");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad mole fraction length"){
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"T");
        std::vector<std::string> fluids(1, "Water&Ethanol");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad input lengths"){
        std::vector<double> p(1, 101325), Q(2, 1.0), z(100, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"Cpmolar"); outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"P",p,"Q",Q,"HEOS",fluids,z););
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad input pair"){
        std::vector<double> Q(2, 1.0), z(1, 1.0);
        std::vector<std::vector<double> > IO;
        std::vector<std::string> outputs(1,"Cpmolar"); outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs,"Q",Q,"Q",Q,"HEOS",fluids,z););
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
};
#endif

/****************************************************
 *                  Props1SI                        *
 ****************************************************/

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
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to Props1SI","[Props1SI],[PropsSI]")
{
    SECTION("Good fluid, good parameter"){
        CHECK(ValidNumber(CoolProp::Props1SI("Tcrit","Water")));
    };
    SECTION("Good fluid, good parameter"){
        CHECK(ValidNumber(CoolProp::PropsSI("Tcrit","",0,"",0,"Water")));
    };
    SECTION("Good fluid, good parameter, inverted"){
        CHECK(ValidNumber(CoolProp::Props1SI("Water","Tcrit")));
    };
    SECTION("Good fluid, bad parameter"){
        CHECK(!ValidNumber(CoolProp::Props1SI("Water","????????????")));
    };
    SECTION("Bad fluid, good parameter"){
        CHECK(!ValidNumber(CoolProp::Props1SI("?????","Tcrit")));
    };
};
#endif


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
    catch (std::exception &){
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
	else if (!ParamName.compare("predefined_mixtures") ){
		return get_csv_predefined_mixtures();
    }
    else{
        throw ValueError(format("Input value [%s] is invalid",ParamName.c_str()));
    }
};
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to get_global_param_string","[get_global_param_string]")
{
    const int num_good_inputs = 8;
    std::string good_inputs[num_good_inputs] = {"version", "gitrevision", "fluids_list", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list","parameter_list","predefined_mixtures"};
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
    try {
        std::string backend, fluid;
        extract_backend(FluidName, backend, fluid);
        if (backend == "INCOMP"){
            shared_ptr<CoolProp::IncompressibleBackend> INCOMP(new CoolProp::IncompressibleBackend(fluid));

            if (!ParamName.compare("long_name")){
                return INCOMP->calc_name();
            }
            else{
                throw ValueError(format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()));
            }
        }

        std::vector<std::string> comps(1, FluidName);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(comps));
        CoolProp::CoolPropFluid *cpfluid = HEOS->get_components()[0];

        if (!ParamName.compare("aliases")){
            return strjoin(cpfluid->aliases, ", ");
        }
        else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number")){
            return cpfluid->CAS;
        }
        else if (!ParamName.compare("ASHRAE34")){
            return cpfluid->environment.ASHRAE34;
        }
        else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname")){
            return cpfluid->REFPROPname;
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

/*
std::string PhaseSI(const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName, const std::vector<double> &z)
{
    double Phase_double = PropsSI("Phase",Name1,Prop1,Name2,Prop2,FluidName,z);
    if (!ValidNumber(Phase_double)){ return "";}
    std::size_t Phase_int = static_cast<std::size_t>(Phase_double);
    return phase_lookup_string(static_cast<phases>(Phase_int));
}
*/
} /* namespace CoolProp */