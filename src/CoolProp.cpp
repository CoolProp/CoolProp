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
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Helmholtz/MixtureParameters.h"
#include "DataStructures.h"

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

//int global_Phase = -1;

void set_warning_string(std::string warning){
    warning_string = warning;
}
void set_error_string(std::string error){
	error_string = error;
}
//
//static int IsCoolPropFluid(std::string FluidName)
//{
//	// Try to get the fluid from Fluids by name
//	try
//	{
//		pFluid = Fluids.get_fluid(FluidName);
//	}
//	catch (NotImplementedError &)
//	{
//		return false;
//	}
//	// If NULL, didn't find it (or its alias)
//	if (pFluid!=NULL)
//	{
//		return true;
//	}
//	else
//		return false;
//}
//
//static int IsBrine(const char* Ref)
//{
//	// First check whether it is one of the Brines that does
//	// not have a pure-fluid equivalent in CoolProp
//    if (
//        strcmp(Ref,"HC-10")==0 ||
//        strncmp(Ref,"PG-",3)==0 ||
//		strncmp(Ref,"EG-",3)==0 ||
//		strncmp(Ref,"EA-",3)==0 ||
//		strncmp(Ref,"MA-",3)==0 ||
//		strncmp(Ref,"Glycerol-",9)==0 ||
//		strncmp(Ref,"K2CO3-",6)==0 ||
//		strncmp(Ref,"CaCl2-",6)==0 ||
//		strncmp(Ref,"MgCl2-",6)==0 ||
//		strncmp(Ref,"NaCl-",5)==0 ||
//		strncmp(Ref,"KAC-",4)==0 ||
//		strncmp(Ref,"KFO-",4)==0 ||
//		strncmp(Ref,"LiCl-",4)==0 ||
//        strncmp(Ref,"NH3/H2O-",8)==0
//       )
//    {
//        return 1;
//    }
//	// Then check for diluants that are also pure fluids in CoolProp
//	else if ( (strncmp(Ref,"Methanol-",9)==0 && Ref[8] == '-') ||
//		      (strncmp(Ref,"Ethanol-",8)==0 && Ref[7] == '-') ||
//			  (strncmp(Ref,"NH3-",4)==0 && Ref[3] == '-')
//		)
//	{
//		return 1;
//	}
//    else
//    {
//        return 0;
//    }
//}
//
//long getFluidType(std::string FluidName){
//	if (IsREFPROP(FluidName)) { return FLUID_TYPE_REFPROP;}
//	else if(IsIncompressibleLiquid(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_LIQUID;}
//	// TODO SOLUTION: Check if working
//	else if(IsIncompressibleSolution(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_SOLUTION; }
//	else {
//		// Try to get the index of the fluid
//		long iFluid = get_Fluid_index(FluidName);
//		// If iFluid is greater than -1, it is a CoolProp Fluid, otherwise not
//		if (iFluid > -1) {
//			// Get a pointer to the fluid object
//			pFluid = get_fluid(iFluid);
//			if (pFluid->pure())	{ return FLUID_TYPE_PURE;}
//			else { return FLUID_TYPE_PSEUDOPURE; }
//		} else {
//			throw ValueError(format("Bad Fluid name [%s] - not a CoolProp fluid",FluidName.c_str()));
//		}
//	}
//	return -1;
//}
//
//EXPORT_CODE int CONVENTION IsFluidType(const char *Ref, const char *Type)
//{
//	pFluid = Fluids.get_fluid(Ref);
//
//	if (IsBrine(Ref)){ // TODO Solution: Remove this part
//		if (!strcmp(Type,"Brine")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsIncompressibleSolution(Ref)){
//		if (!strcmp(Type,"Solution")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsIncompressibleLiquid(Ref)){
//		if (!strcmp(Type,"Liquid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsREFPROP(Ref)){
//		if (!strcmp(Type,"PureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (!pFluid->pure()){
//		if (!strcmp(Type,"PseudoPure") || !strcmp(Type,"PseudoPureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (pFluid->pure()){
//		if (!strcmp(Type,"PureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//    else
//    {
//        return 0;
//    }
//}
//
//
//std::string Phase_Trho(std::string Fluid, double T, double rho)
//{
//	try{
//		// Try to load the CoolProp Fluid
//		pFluid = Fluids.get_fluid(Fluid);
//		double pL,pV,rhoL,rhoV;
//		return pFluid->phase_Trho(T,rho, pL, pV, rhoL, rhoV);
//	}
//	catch(NotImplementedError &){
//		return std::string("");
//	}
//	return std::string("");
//}
//
//std::string Phase(std::string Fluid, double T, double p)
//{
//	try{
//		// Try to load the CoolProp Fluid
//		pFluid = Fluids.get_fluid(Fluid);
//		double pL,pV,rhoL,rhoV;
//		return pFluid->phase_Tp(T, p, pL, pV, rhoL, rhoV);
//	}
//	catch(NotImplementedError &){
//		return std::string("");
//	}
//	return std::string("");
//}
//
//std::string Phase_Tp(std::string Fluid, double T, double p)
//{
//	return Phase(Fluid,T,p);
//}
//

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
    // For backwards compatibility reasons, if "REFPROP-" or "REFPROP-MIX:" start the fluid_string, replace them with "REFPROP::"
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
		if (!strcmp(pEnd,"%")){	x *= 0.01;}
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
    
	// First check if it is a trivial input (critical/max parameters for instance)
	if (is_valid_parameter(Output, iOutput) && is_trivial_parameter(iOutput))
	{
		double val = State->trivial_keyed_output(iOutput);
		return val;
	}
	
	long iName1 = get_parameter_index(Name1);
	long iName2 = get_parameter_index(Name2);

	if (State->using_mole_fractions()){
		State->set_mole_fractions(fractions);
	} else if (State->using_mass_fractions()){
		State->set_mass_fractions(fractions);
	} else if (State->using_volu_fractions()){
		State->set_volu_fractions(fractions);
	} else {
		if (get_debug_level()>50) std::cout << format("%s:%d: _PropsSI, could not set composition to %s, defaulting to mole fraction.\n",__FILE__,__LINE__, vec_to_string(z).c_str()).c_str();
	}

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
    std::string  _FluidName = FluidName, empty_string = "", _Output = Output;
    double val1 = PropsSI(_FluidName, empty_string, 0, empty_string, 0, _Output);
    if (!ValidNumber(val1)){
        // flush the error
        set_error_string("");
        // Try with them flipped
        val1 = PropsSI(_Output, empty_string, 0, empty_string, 0, _FluidName);
    }
    if (!ValidNumber(val1)){
        set_error_string(format("Unable to use inputs %s,%s in Props1SI (order doesn't matter)"));
        return _HUGE;
    }
    return val1;
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
        std::cout << format("Sizes of Prop1 [%d] and Prop2 [%d] to PropsSI are not the same", Prop1.size(), Prop2.size()) << std::endl;
        return out;
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


//EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
//{
//    Prop1 = convert_from_unit_system_to_SI(iName1, Prop1, get_standard_unit_system());
//    Prop2 = convert_from_unit_system_to_SI(iName2, Prop2, get_standard_unit_system());
//	double out = IPropsSI(iOutput,iName1,Prop1,iName2,Prop2,iFluid);
//    return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
//}
//double _CoolProp_Fluid_PropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, Fluid *pFluid)
//{
//	double val = _HUGE, T = _HUGE;
//	// This private method uses the indices directly for speed
//
//	if (get_debug_level()>3){
//		std::cout << format("%s:%d: _CoolProp_Fluid_PropsSI(%d,%d,%g,%d,%g,%s)\n",__FILE__,__LINE__,iOutput,iName1, Prop1, iName2, Prop2, pFluid->get_name().c_str()).c_str();
//	}
//    if (iName1 == iT){
//        T = Prop1;}
//    else if (iName2 == iT){
//        T = Prop2;}
//
//	// Generate a State instance wrapped around the Fluid instance
//	CoolPropStateClassSI CPS(pFluid);
//
//	// Check if it is an output that doesn't require a state input
//	// Deal with it and return
//	switch (iOutput)
//	{
//        case iI:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->surface_tension_T(T);
//            }
//        case iRhosatLanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->rhosatL(T);
//            }
//        case iRhosatVanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->rhosatV(T);
//            }
//        case iPsatLanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            if (CPS.pFluid->pure()){
//                return CPS.pFluid->psat(T);
//            }
//            else{
//                return CPS.pFluid->psatL(T);
//            }
//            }
//        case iPsatVanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            if (CPS.pFluid->pure()){
//                return CPS.pFluid->psat(T);
//            }
//            else{
//                return CPS.pFluid->psatV(T);
//            }
//            }
//		case iMM:
//		case iPcrit:
//		case iTcrit:
//		case iTtriple:
//		case iPtriple:
//      case iPmax:
//      case iTmax:
//		case iRhocrit:
//		case iTmin:
//		case iAccentric:
//		case iPHASE_LIQUID:
//		case iPHASE_GAS:
//		case iPHASE_SUPERCRITICAL:
//		case iPHASE_TWOPHASE:
//		case iGWP20:
//		case iGWP100:
//		case iGWP500:
//		case iODP:
//		case iCritSplineT:
//		case iScrit:
//		case iHcrit:
//		case iTreduce:
//		case iRhoreduce:
//			return CPS.keyed_output(iOutput);
//	}
//
//	// Update the class
//	CPS.update(iName1,Prop1,iName2,Prop2);
//
//	// Debug
//	if (get_debug_level()>9){std::cout << format("%s:%d: State update successful\n",__FILE__,__LINE__).c_str();}
//
//	// Get the output
//	val = CPS.keyed_output(iOutput);
//
//	// Debug
//	if (get_debug_level()>5){std::cout << format("%s:%d: _CoolProp_Fluid_PropsSI returns: %g\n",__FILE__,__LINE__,val).c_str();}
//
//	// Return the value
//	return val;
//}
//EXPORT_CODE double CONVENTION IPropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
//{
//	pFluid = Fluids.get_fluid(iFluid);
//	// Didn't work
//	if (pFluid == NULL){
//		err_string=std::string("CoolProp error: ").append(format("%d is an invalid fluid index to IProps",iFluid));
//		return _HUGE;
//	}
//	else{
//		// In this function the error catching happens;
//		try{
//			// This is already converted to the right units since we take in SI units
//			return _CoolProp_Fluid_PropsSI(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
//		}
//		catch(std::exception &e){
//			err_string=std::string("CoolProp error: ").append(e.what());
//			return _HUGE;
//		}
//		catch(...){
//			err_string=std::string("CoolProp error: Indeterminate error");
//			return _HUGE;
//		}
//	}
//}

///// Calculate some interesting derivatives
//double _CoolProp_Deriv_Terms(long iTerm, double T, double rho, Fluid * pFluid)
//{
//	double val = _HUGE;
//	// This private method uses the indices directly for speed
//
//	if (get_debug_level()>3){
//		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
//	}
//
//	switch (iTerm) {
//		case iDERdh_dp__rho:
//		case iDERdh_dp__v:
//		case iDERZ:
//		case iDERdZ_dDelta:
//		case iDERdZ_dTau:
//		case iDERB:
//		case iDERdB_dT:
//		case iDERC:
//		case iDERdC_dT:
//		case iDERphir:
//		case iDERdphir_dTau:
//		case iDERdphir_dDelta:
//		case iDERd2phir_dTau2:
//		case iDERd2phir_dDelta2:
//		case iDERd2phir_dDelta_dTau:
//		case iDERd3phir_dDelta3:
//		case iDERd3phir_dDelta2_dTau:
//		case iDERd3phir_dDelta_dTau2:
//		case iDERd3phir_dTau3:
//		case iDERphi0:
//		case iDERdphi0_dTau:
//		case iDERd2phi0_dTau2:
//		case iDERdphi0_dDelta:
//		case iDERd2phi0_dDelta2:
//		case iDERd2phi0_dDelta_dTau:
//		case iDERd3phi0_dTau3:
//			{
//			// Generate a State instance wrapped around the Fluid instance
//			CoolPropStateClass CPS(pFluid);
//
//			// Force the update to consider the inputs as single-phase inputs
//			CPS.flag_SinglePhase =  true;
//
//			// Update the class
//			CPS.update(iT,T,iD,rho);
//
//			// Get the output value
//			val = CPS.keyed_output(iTerm);
//			break;
//			}
//
//		case iDERrho_smoothed:
//		case iDERdrho_smoothed_dh:
//		case iDERdrho_smoothed_dp:
//		case iDERdrhodh_constp_smoothed:
//		case iDERdrhodp_consth_smoothed:
//		case iDERIsothermalCompressibility:
//			{
//			// Generate a State instance wrapped around the Fluid instance
//			CoolPropStateClass CPS(pFluid);
//
//			// Update the class
//			CPS.update(iT,T,iD,rho);
//
//			// Get the output value
//			val = CPS.keyed_output(iTerm);
//			break;
//			}
//
//		default:
//			throw ValueError(format("Sorry DerivTerms is a work in progress, your derivative term [%d] is not available!",iTerm));
//	}
//
//	if (get_debug_level()>5){
//		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
//	}
//	// Return the value
//	return val;
//}
//
//// Define the functions from the header
//double DerivTerms(long iTerm, double T, double rho, Fluid * pFluid){
//	return _CoolProp_Deriv_Terms(iTerm,T,rho,pFluid);
//}
//double DerivTerms(std::string Term, double T, double rho, std::string Fluidname){
//	if (get_debug_level()>5){
//			std::cout<<__FILE__<<": "<<Term.c_str()<<",T="<<T<<",rho="<<rho<<","<<Fluidname.c_str()<<std::endl;
//		}
//		/*
//	    Derivatives are only supported for CoolProp fluids
//	    */
//	    if (IsCoolPropFluid(Fluidname))
//		{
//			pFluid = Fluids.get_fluid(Fluidname);
//			// for compatibility, replace B and C with VB and VC
//			if ((!Term.compare("B")) || (!Term.compare("C"))) {
//				Term = std::string("V").append(Term);
//			}
//			// Convert all the parameters to integers
//			long iOutput = get_param_index(Term);
//			if (iOutput<0)
//				throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Term.c_str()));
//
//			if (T<=0)
//				throw ValueError(format("Your input temperature [%f] is not valid.",T));
//
//			if (rho<=0)
//				throw ValueError(format("Your input density [%f] is not valid.",rho));
//			// Call the internal method that uses the parameters converted to longs
//			return _CoolProp_Deriv_Terms(iOutput,T,rho,pFluid);
//		}
//		else
//		{
//			throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid.",Fluidname.c_str()));
//		}
//}
//
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
        double delta_a1 = deltas/(8.314472/HEOS->molar_mass());
        double delta_a2 = -deltah/(8.314472/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "IIR");
        HEOS->update_states();
	}
	else if (!reference_state.compare("ASHRAE"))
	{
        HEOS->update(QT_INPUTS, 0, 243.15);

		// Get current values for the enthalpy and entropy
		double deltah = HEOS->hmass() - 0; // offset from 0 J/kg enthalpy
		double deltas = HEOS->smass() - 0; // offset from 0 J/kg/K entropy
		double delta_a1 = deltas/(8.314472/HEOS->molar_mass());
        double delta_a2 = -deltah/(8.314472/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "ASHRAE");
        HEOS->update_states();
	}
	else if (!reference_state.compare("NBP"))
	{
		// Saturated liquid boiling point at 1 atmosphere
        HEOS->update(PQ_INPUTS, 101325, 0);

		double deltah = HEOS->hmass() - 0; // offset from 0 kJ/kg enthalpy
		double deltas = HEOS->smass() - 0; // offset from 0 kJ/kg/K entropy
		double delta_a1 = deltas/(8.314472/HEOS->molar_mass());
        double delta_a2 = -deltah/(8.314472/HEOS->molar_mass()*HEOS->get_reducing_state().T);
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, "NBP");
        HEOS->update_states();
	}
    else if (!reference_state.compare("DEF"))
    {
        //HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(0,0);
        throw NotImplementedError("Default reference state has not been implemented yet");
    }
    else if (!reference_state.compare("RESET"))
    {
        HEOS->get_components()[0]->pEOS->alpha0.EnthalpyEntropyOffset.set(0, 0, "");
    }
	else
	{
        throw ValueError(format("reference state string is invalid: [%s]",reference_state.c_str()));
	}
}
int set_reference_stateD(std::string Ref, double T, double rhomolar, double h0, double s0)
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
	else{ return "Bad key";}
}
std::string get_global_param_string(std::string ParamName)
{
	if (!ParamName.compare("version")){
		return version;
	}
	else if (!ParamName.compare("gitrevision")){
		return gitrevision;
	}
    else if (!ParamName.compare("errstring")){
		std::string temp = error_string;
		error_string = "";
		return temp;
	}
	else if (!ParamName.compare("warnstring")){
		std::string temp = warning_string;
		warning_string = "";
		return temp;
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
	else if (!ParamName.compare("parameter_list") )
    {
        return get_csv_parameter_list();
    }
	else{
		return format("Input value [%s] is invalid",ParamName.c_str()).c_str();
	}
};
std::string get_fluid_param_string(std::string FluidName, std::string ParamName)
{
	try{
        std::vector<std::string> comps(1,FluidName);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(comps));

        CoolProp::CoolPropFluid *fluid = HEOS->get_components()[0];

		if (!ParamName.compare("aliases"))
		{
			return strjoin(fluid->aliases, ", ");
		}
		else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number"))
		{
            return fluid->CAS;
		}
		else if (!ParamName.compare("ASHRAE34"))
		{
            return fluid->environment.ASHRAE34;
		}
		else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname"))
		{
            return fluid->REFPROPname;
		}
        else if (ParamName.find("BibTeX") == 0) // Starts with "BibTeX"
        {
            std::vector<std::string> parts = strsplit(ParamName,'-');
            // 
            std::string item = parts[1];
            if (item == "EOS"){
                return fluid->pEOS->BibTeX_EOS;
            }
            else if (item == "CP0"){
                return fluid->pEOS->BibTeX_CP0;
            }
            else if (item == "SURFACE_TENSION"){
                return fluid->ancillaries.surface_tension.BibTeX;
            }
            else if (item == "MELTING_LINE"){
                return fluid->ancillaries.melting_line.BibTeX;
            }
            else if (item == "VISCOSITY"){
                return fluid->transport.BibTeX_viscosity;
            }
            else if (item == "CONDUCTIVITY"){
                return fluid->transport.BibTeX_conductivity;
            }
            else{
                return format("Could not match BibTeX item: %s", item.c_str());
            }
        }
		else
		{
			return format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str());
		}
	}
	catch(std::exception &e)
	{
		return(std::string("CoolProp error: ").append(e.what()));
	}
	catch(...){
		return(std::string("CoolProp error: Indeterminate error"));
	}
}
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
