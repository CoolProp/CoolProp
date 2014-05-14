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

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <exception>
#include <stdio.h>
#include <string>
#include "CoolPropTools.h"
#include "Solvers.h"
#include "Fluids\FluidLibrary.h"

namespace CoolProp
{

static int debug_level = 0;
static std::string error_string;
static std::string warning_string;

void set_debug_level(int level){debug_level = level;}
int get_debug_level(void){return debug_level;}

//// This is very hacky, but pull the git revision from the file
#include "gitrevision.h" // Contents are like "std::string gitrevision = "aa121435436ggregrea4t43t433";"
#include "version.h" // Contents are like "char version [] = "2.5";"

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

bool has_fractions_in_string(const std::string &fluid_string)
{
    // Start at the "::" if it is found to chop off the delimiter between backend and fluid
    std::size_t i = fluid_string.find("::");
    // If can find both "&" and "]", it must have mole fractions encoded as string
    return fluid_string.find("&",i+2) != -1 && fluid_string.find("]",i+2);
}
std::string extract_fractions(const std::string &fluid_string, std::vector<double> &fractions)
{
    fractions.clear();
    std::vector<std::string> names;
    std::string all_pairs, backend_string;
    // Start at the "::" if it is found to chop off the delimiter between backend and fluid
    int i = fluid_string.find("::");
    
    // If no backend/fluid delimiter
    if (i < 0){
        // Just use the full string
        backend_string = "";
        all_pairs = fluid_string;
    }
    else{
        // Part including the ::
        backend_string = fluid_string.substr(0, i+2);
        // Keep the part after "::"
        all_pairs = fluid_string.substr(i+2);
    }

    // Break up into pairs (like "Methane[0.5]")
    std::vector<std::string> pairs = strsplit(all_pairs, '&');

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

    // Join fluids back together
    return backend_string + strjoin(names, "&");
}

// Internal function to do the actual calculations, make this a wrapped function so
// that error bubbling can be done properly
double _PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref, const std::vector<double> &z)
{
    static std::string unknown_backend = "?";
    double x1, x2;
    if (get_debug_level()>5){
        std::cout << format("%s:%d: _PropsSI(%s,%s,%g,%s,%g,%s)\n",__FILE__,__LINE__,Output.c_str(),Name1.c_str(),Prop1,Name2.c_str(),Prop2,Ref.c_str()).c_str();
    }

    // Convert all the input and output parameters to integers
    long iOutput = get_parameter_index(Output);
    long iName1 = get_parameter_index(Name1);
    long iName2 = get_parameter_index(Name2);

    // The state we are going to use
    AbstractState *State = NULL;
    try
    {
        // We are going to let the factory function determine which backend to use
        //
        // Generate the State class pointer using the factory function with unknown backend
        State = AbstractState::factory(unknown_backend, Ref);

        if (State == NULL){ throw ValueError("unable to instantiate AbstractState*");}

        if (State->using_mole_fractions()){
            State->set_mole_fractions(z);
        }
        else{
            State->set_mass_fractions(z);
        }

        // Obtain the input pair
        long pair = generate_update_pair(iName1, Prop1, iName2, Prop2, x1, x2);

        // First check if it is a trivial input (critical/max parameters for instance)
        // TODO: check for trivial inputs that do not require the use of the eos
        /*if (State->is_trivial_output(iOutput))
        { 
            double val = State->trivial_keyed_output(iOutput);
            delete(State);
            return val;
        };*/

        // Update the state
        State->update(pair, x1, x2);

        // Return the desired output
        double val = State->keyed_output(iOutput);
        
        // Free the pointer to the State class
        delete (State);

        // Return the value
        return val;
    }
    catch(...){	
        delete(State); throw;
    }
}
double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref, const std::vector<double> &z)
{
    CATCH_ALL_ERRORS_RETURN_HUGE(return _PropsSI(Output,Name1,Prop1,Name2,Prop2,Ref,z);)
}
double PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *FluidName, const std::vector<double> &x)
{
    std::string _Output = Output, _Name1 = Name1, _Name2 = Name2, _FluidName = FluidName;
    return PropsSI(_Output,_Name1,Prop1,_Name2,Prop2,_FluidName, x);
}
double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &Ref)
{
    // In this function the error catching happens;
    try{          
        // Extract fractions if they are encoded in the fluid string
        if (has_fractions_in_string(Ref))
        {
            std::vector<double> fractions;
            // Extract the fractions and reformulate the list of fluids REFPROP::Methane[0.5]&Ethane[0.5] -> REFPROP::Methane&Ethane and [0.5,0.5]
            std::string Ref2 = extract_fractions(Ref, fractions);
            return _PropsSI(Output, Name1, Prop1, Name2, Prop2, Ref2, fractions);
        }
        else
        {
            // Here it must be a pure fluid since the fractions are not coded in the string.  If not, it will fail internally
            return _PropsSI(Output, Name1, Prop1, Name2, Prop2, Ref, std::vector<double>(1,1));
        }
    }
    catch(const std::exception& e){ set_error_string(e.what()); return _HUGE; }
    catch(...){ return _HUGE; }
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
//int set_reference_stateS(std::string Ref, std::string reference_state)
//{
//	Fluid *pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		return set_reference_stateP(pFluid, reference_state);
//	}
//	else{
//		return -1;
//	}
//}
//
//int set_reference_stateP(Fluid *pFluid, std::string reference_state)
//{
//	CoolPropStateClassSI CPS(pFluid);
//	if (!reference_state.compare("IIR"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,273.15,iQ,0);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-200000; // offset from 200 kJ/kg enthalpy
//		double deltas = s1-1000; // offset from 1 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else if (!reference_state.compare("ASHRAE"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,233.15,iQ,0);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
//		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else if (!reference_state.compare("NBP"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iP,101325.0,iQ,0); // Saturated boiling point at 1 atmosphere
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
//		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else
//	{ 
//		return -1;
//	}
//
//}
//int set_reference_stateD(std::string Ref, double T, double rho, double h0, double s0)
//{
//	pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,T,iD,rho);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-h0; // offset from given enthalpy in SI units
//		double deltas = s1-s0; // offset from given enthalpy in SI units
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else{
//		return -1;
//	}
//}

//std::string get_BibTeXKey(std::string Ref, std::string item)
//{
//	pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		
//		if (!item.compare("EOS")){ return pFluid->BibTeXKeys.EOS; }
//		else if (!item.compare("CP0")){ return pFluid->BibTeXKeys.CP0; }
//		else if (!item.compare("VISCOSITY")){ return pFluid->BibTeXKeys.VISCOSITY; }
//		else if (!item.compare("CONDUCTIVITY")){ return pFluid->BibTeXKeys.CONDUCTIVITY; }
//		else if (!item.compare("ECS_LENNARD_JONES")){ return pFluid->BibTeXKeys.ECS_LENNARD_JONES; }
//		else if (!item.compare("ECS_FITS")){ return pFluid->BibTeXKeys.ECS_FITS; }
//		else if (!item.compare("SURFACE_TENSION")){ return pFluid->BibTeXKeys.SURFACE_TENSION; }
//		else{ return "Bad key";}
//	}
//	else{
//		return std::string("");
//	}
//}
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
		error_string = std::string("");
		return temp;
	}
	else if (!ParamName.compare("warnstring")){
		std::string temp = warning_string;
		warning_string = std::string("");
		return temp;
	}
	else if (!ParamName.compare("FluidsList") || !ParamName.compare("fluids_list") || !ParamName.compare("fluidslist")){
		return get_fluid_list();
	}
	else{
		return format("Input value [%s] is invalid",ParamName.c_str()).c_str();
	}
};
//std::string get_fluid_param_string(std::string FluidName, std::string ParamName)
//{
//	try{
//		pFluid = Fluids.get_fluid(FluidName);
//		// Didn't work
//		if (pFluid == NULL){
//			err_string=std::string("CoolProp error: ").append(format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()));
//			return format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()).c_str();
//		}
//		else{
//			if (!ParamName.compare("aliases"))
//			{
//				std::vector<std::string> v = pFluid->get_aliases();
//				return strjoin(v,", ");
//			}
//			else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number"))
//			{
//				return pFluid->params.CAS;
//			}
//			else if (!ParamName.compare("ASHRAE34"))
//			{
//				return pFluid->environment.ASHRAE34;
//			}
//			else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname"))
//			{
//				return pFluid->get_REFPROPname();
//			}
//			else if (!ParamName.compare("TTSE_mode"))
//			{
//				int mode = pFluid->TTSESinglePhase.get_mode();
//				switch (mode)
//				{
//				case TTSE_MODE_TTSE:
//					return "TTSE";
//				case TTSE_MODE_BICUBIC:
//					return "BICUBIC";
//				default:
//					throw ValueError("TTSE mode is invalid");
//				}
//			}
//			else
//			{
//				return format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()).c_str();
//			}
//		}
//	}
//	catch(std::exception &e)
//	{
//		return(std::string("CoolProp error: ").append(e.what()));
//	}
//	catch(...){
//		return(std::string("CoolProp error: Indeterminate error"));
//	}
//}

} /* namespace CoolProp */