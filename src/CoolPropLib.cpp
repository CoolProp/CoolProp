#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#else
#include <fenv.h>
#endif

#include "CoolPropLib.h"
#include "CoolProp.h"
#include "HumidAirProp.h"
#include "DataStructures.h"
#include "Exceptions.h"
#include "float.h"
#include "crossplatform_shared_ptr.h"
#include "AbstractState.h"
#include "Exceptions.h"

#include <string.h>

bool str2buf(const std::string& str, char * buf, int n)
{
  if (str.size() < static_cast<unsigned int>(n)) {
    strcpy(buf, str.c_str());
    return true;
  }
  return false;
}

// In Microsoft Excel, they seem to check the FPU exception bits and error out because of it.  
// By calling the _clearfp(), we can reset these bits, and not get the error
// See also http://stackoverflow.com/questions/11685441/floating-point-error-when-calling-dll-function-from-vba/27336496#27336496
// See also http://stackoverflow.com/questions/16849009/in-linux-do-there-exist-functions-similar-to-clearfp-and-statusfp for linux and OSX
void reset_fpu()
{
    #if defined(_MSC_VER)
        _clearfp(); // For MSVC, clear the floating point error flags
    #elif defined(FE_ALL_EXCEPT)
        feclearexcept(FE_ALL_EXCEPT);
    #endif
}
double convert_from_kSI_to_SI(long iInput, double value)
{
    if (get_debug_level() > 8){
        std::cout << format("%s:%d: convert_from_kSI_to_SI(i=%d,value=%g)\n",__FILE__,__LINE__,iInput,value).c_str();
    }

    switch (iInput)
    {
    case CoolProp::iP:
    case CoolProp::iCpmass:
    case CoolProp::iCp0mass:
    case CoolProp::iSmass:
    case CoolProp::iGmass:
    case CoolProp::iCvmass:
    case CoolProp::iHmass:
    case CoolProp::iUmass:
    case CoolProp::iconductivity:
        return value*1000.0;
    case CoolProp::iDmass:
    case CoolProp::ispeed_sound:
    case CoolProp::iQ:
    case CoolProp::iviscosity:
    case CoolProp::iT:
    case CoolProp::iPrandtl:
    case CoolProp::isurface_tension:
        return value;
    default:
        throw CoolProp::ValueError(format("index [%d] is invalid in convert_from_kSI_to_SI",iInput).c_str());
    }
}

double convert_from_SI_to_kSI(long iInput, double value)
{
    if (get_debug_level() > 8){
        std::cout << format("%s:%d: convert_from_SI_to_kSI(%d,%g)\n",__FILE__,__LINE__,iInput,value).c_str();
    }

    switch (iInput)
    {
    case CoolProp::iP:
    case CoolProp::iCpmass:
    case CoolProp::iCp0mass:
    case CoolProp::iSmass:
    case CoolProp::iGmass:
    case CoolProp::iCvmass:
    case CoolProp::iHmass:
    case CoolProp::iUmass:
    case CoolProp::iconductivity:
        return value/1000.0;
    case CoolProp::iDmass:
    case CoolProp::iQ:
    case CoolProp::ispeed_sound:
    case CoolProp::iviscosity:
    case CoolProp::iT:
    case CoolProp::isurface_tension:
        return value;
    default:
        throw CoolProp::ValueError(format("index [%d] is invalid in convert_from_SI_to_kSI", iInput).c_str());
    }
}

EXPORT_CODE long CONVENTION redirect_stdout(const char* file){
    FILE *fp = freopen(file, "a+", stdout);
    reset_fpu();
    return (fp) ? 1 : 0; // 0 = failure if redirection could not be done; original stdout is already closed
}
EXPORT_CODE int CONVENTION set_reference_stateS(const char *Ref, const char *reference_state)
{
    try{
        CoolProp::set_reference_stateS(std::string(Ref), std::string(reference_state));
        reset_fpu();
        return true;
    }
    catch(...){
        reset_fpu();
        return false;
    }
}
EXPORT_CODE int CONVENTION set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0)
{
    try{
        CoolProp::set_reference_stateD(std::string(Ref), T, rho, h0, s0);
        reset_fpu();
        return true;
    }
    catch(...){
        reset_fpu();
        return false;
    }
}

// All the function interfaces that point to the single-input Props function
EXPORT_CODE double CONVENTION Props1(const char *FluidName, const char *Output){
    double val = PropsS(Output, "", 0, "", 0, FluidName);
    reset_fpu();
    return val;
}
EXPORT_CODE double CONVENTION PropsS(const char *Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char * Ref){
    double val = Props(Output,Name1[0],Prop1,Name2[0],Prop2,Ref);
    reset_fpu();
    return val;
}
EXPORT_CODE double CONVENTION Props(const char *Output, const char Name1, double Prop1, const char Name2, double Prop2, const char * Ref)
{
    try
    {
        // Get parameter indices
        std::string sName1 = std::string(1, Name1), sName2 = std::string(1, Name2);
        long iOutput = get_param_index(Output);
        long iName1 = CoolProp::get_parameter_index(sName1);
        long iName2 = CoolProp::get_parameter_index(sName2);

        // Convert inputs to SI
        Prop1 = convert_from_kSI_to_SI(iName1, Prop1);
        Prop2 = convert_from_kSI_to_SI(iName2, Prop2);

        // Call the SI function
        double val = PropsSI(Output, sName1.c_str(), Prop1, sName2.c_str(), Prop2, Ref);

        reset_fpu();

        // Convert back to unit system
        return convert_from_SI_to_kSI(iOutput, val);
    }
    catch(std::exception &e){CoolProp::set_error_string(e.what()); return _HUGE;}
    catch(...){CoolProp::set_error_string("Undefined error"); return _HUGE;}
}
EXPORT_CODE double CONVENTION saturation_ancillary(const char *fluid_name, const char *output, int Q, const char *input, double value)
{
    try
    {
        double val = CoolProp::saturation_ancillary(fluid_name, std::string(output), Q, std::string(input), value);
        reset_fpu();
        return val;
    }
    catch(std::exception &e){CoolProp::set_error_string(e.what()); return _HUGE;}
    catch(...){CoolProp::set_error_string("Undefined error"); return _HUGE;}
}
EXPORT_CODE double CONVENTION Props1SI(const char *FluidName, const char *Output)
{
    double val = CoolProp::Props1SI(std::string(FluidName), std::string(Output));
    reset_fpu();
    return val;
}
EXPORT_CODE double CONVENTION PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char * FluidName)
{
    double val = CoolProp::PropsSI(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(FluidName));
    reset_fpu();
    return val;
}
EXPORT_CODE long CONVENTION PhaseSI(const char *Name1, double Prop1, const char *Name2, double Prop2, const char * FluidName, char *phase, int n)
{
    std::string s = CoolProp::PhaseSI(std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(FluidName));
    reset_fpu();
    return str2buf(s, phase, n) ? 1 : 0;
}
/*
 * EXPORT_CODE double CONVENTION PropsSIZ(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char * FluidName, const double *z, int n)
{
    std::string _Output = Output, _Name1 = Name1, _Name2 = Name2, _FluidName = FluidName;
    double val = CoolProp::PropsSI(_Output, _Name1, Prop1, _Name2, Prop2, _FluidName, std::vector<double>(z, z+n));
    reset_fpu();
    return val;
}
 * */
EXPORT_CODE void CONVENTION propssi_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char * FluidName, double *output)
{
    *output = PropsSI(Output, Name1, *Prop1, Name2, *Prop2, FluidName);
}

EXPORT_CODE double CONVENTION K2F(double T){
    return T * 9 / 5 - 459.67;
}
EXPORT_CODE double CONVENTION F2K(double T_F){
    return (T_F + 459.67) * 5 / 9;
}
EXPORT_CODE int CONVENTION get_debug_level(){
    return CoolProp::get_debug_level();
}
EXPORT_CODE void CONVENTION set_debug_level(int level){
    CoolProp::set_debug_level(level);
}
EXPORT_CODE long CONVENTION get_param_index(const char * param){
    try{
        return CoolProp::get_parameter_index(param);
    }
    catch(...){
        return -1;
    }
}
EXPORT_CODE long CONVENTION get_input_pair_index(const char * pair){
    try{
        return CoolProp::get_input_pair_index(pair);
    }
    catch(...){
        return -1;
    }
}
EXPORT_CODE long CONVENTION get_global_param_string(const char *param, char * Output, int n)
{
    try{
        std::string s = CoolProp::get_global_param_string(param);
        return str2buf(s, Output, n) ? 1 : 0;
    }
    catch(...){
        return 0;
    }
}
EXPORT_CODE long CONVENTION get_parameter_information_string(const char *param, char * Output, int n)
{
    try{
        int key = CoolProp::get_parameter_index(param);
        std::string s = CoolProp::get_parameter_information(key, "long");
        return str2buf(s, Output, n) ? 1 : 0;
    }
    catch(...){
		str2buf(format("parameter is invalid: %s", param), Output, n);
	}
    return 0;
}
EXPORT_CODE long CONVENTION get_fluid_param_string(const char *fluid, const char *param, char * Output, int n)
{
    try{
        std::string s = CoolProp::get_fluid_param_string(std::string(fluid), std::string(param));
        return str2buf(s, Output, n) ? 1 : 0;
    }
    catch(...){
        return 0;
    }
}
EXPORT_CODE double CONVENTION HAPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char * Name3, double Prop3)
{
    double val = HumidAir::HAPropsSI(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(Name3), Prop3);
    reset_fpu();
    return val;
}
EXPORT_CODE double CONVENTION cair_sat(double T)
{
    double val = HumidAir::cair_sat(T);
    reset_fpu();
    return val;
}
EXPORT_CODE void CONVENTION hapropssi_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char * Name3, const double * Prop3, double *output)
{
    *output = HAPropsSI(Output, Name1, *Prop1, Name2, *Prop2, Name3, *Prop3);
}
EXPORT_CODE double CONVENTION HAProps(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char * Name3, double Prop3)
{
    double val = HumidAir::HAProps(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(Name3), Prop3);
    reset_fpu();
    return val;
}
EXPORT_CODE void CONVENTION haprops_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char * Name3, const double * Prop3, double *output)
{
    *output = HAProps(Output, Name1, *Prop1, Name2, *Prop2, Name3, *Prop3);
}
class AbstractStateLibrary{
private:
    std::map<std::size_t, shared_ptr<CoolProp::AbstractState> > ASlibrary;
    long next_handle;
public:
    AbstractStateLibrary(): next_handle(0){};
    long add(shared_ptr<CoolProp::AbstractState> AS){
        ASlibrary.insert(std::pair<std::size_t, shared_ptr<CoolProp::AbstractState> >(this->next_handle,  AS));
        this->next_handle++;
        return next_handle-1;
    }
    void remove(long handle){
        std::size_t count_removed = ASlibrary.erase(handle);
        if (count_removed != 1){
            throw CoolProp::HandleError("could not free handle");
        }
    }
    shared_ptr<CoolProp::AbstractState> & get(long handle){
        std::map<std::size_t, shared_ptr<CoolProp::AbstractState> >::iterator it = ASlibrary.find(handle);
        if (it != ASlibrary.end()){
            return it->second;
        }
        else{
            throw CoolProp::HandleError("could not get handle");
        }
    }
};
static AbstractStateLibrary handle_manager;

EXPORT_CODE long CONVENTION AbstractState_factory(const char* backend, const char* fluids, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, fluids));
        return handle_manager.add(AS);
    }
    catch(CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(...){
        *errcode = 3;
    }
    return -1;
}
EXPORT_CODE void CONVENTION AbstractState_free(const long handle, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        handle_manager.remove(handle);
    }
    catch(CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(...){
        *errcode = 3;
    }
}
EXPORT_CODE void CONVENTION AbstractState_set_fractions(const long handle, const double *fractions, const long N, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    std::vector<double> _fractions(fractions, fractions + N);
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);
        if (AS->using_mole_fractions()){
            AS->set_mole_fractions(_fractions);
        }
        else if (AS->using_mass_fractions()){
            AS->set_mass_fractions(_fractions);
        }
        else if (AS->using_volu_fractions()){
            AS->set_volu_fractions(_fractions);
        }
    }
    catch(CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(...){
        *errcode = 3;
    }
}
EXPORT_CODE void CONVENTION AbstractState_update(const long handle, const long input_pair, const double value1, const double value2, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);
        AS->update(static_cast<CoolProp::input_pairs>(input_pair), value1, value2);
    }
    catch(CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(...){
        *errcode = 3;
    }
}
EXPORT_CODE double CONVENTION AbstractState_keyed_output(const long handle, const long param, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);
        return AS->keyed_output(static_cast<CoolProp::parameters>(param));
    }
    catch(CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch(...){
        *errcode = 3;
    }
    return _HUGE;
}

EXPORT_CODE double CONVENTION AbstractState_first_saturation_deriv(const long handle, const long Of, const long Wrt, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);
        return AS->first_saturation_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt));
    }
    catch (CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (...){
        *errcode = 3;
    }
    return _HUGE;
}

EXPORT_CODE double CONVENTION AbstractState_first_partial_deriv(const long handle, const long Of, const long Wrt, const long Constant, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);
        return AS->first_partial_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt), static_cast<CoolProp::parameters>(Constant));
    }
    catch (CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (CoolProp::CoolPropBaseError &e){
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (...){
        *errcode = 3;
    }
    return _HUGE;
}

EXPORT_CODE void CONVENTION AbstractState_update_and_common_out(const long handle, const long input_pair, const double* value1, const double* value2, const long length, double* T, double* p, double* rhomolar, double* hmolar, double* smolar, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);

        for (int i = 0; i<length; i++){
            try{
                AS->update(static_cast<CoolProp::input_pairs>(input_pair), *(value1+i), *(value2+i));
                *(T+i) = AS->T();
                *(p+i) = AS->p();
                *(rhomolar+i) = AS->rhomolar();
                *(hmolar+i) = AS->hmolar();
                *(smolar+i) = AS->smolar();
            }
            catch (...){
                
            }
        };
    }
    catch (CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (...){
        *errcode = 3;
    }
}

EXPORT_CODE void CONVENTION AbstractState_update_and_5_out(const long handle, const long input_pair, const double* value1, const double* value2, const long length, long* outputs, double* out1, double* out2, double* out3, double* out4, double* out5, long *errcode, char *message_buffer, const long buffer_length)
{
    *errcode = 0;
    try{
        shared_ptr<CoolProp::AbstractState> &AS = handle_manager.get(handle);

        for (int i = 0; i<length; i++){
            try{
                AS->update(static_cast<CoolProp::input_pairs>(input_pair), *(value1+i), *(value2+i));
                *(out1+i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[0]));
                *(out2+i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[1]));
                *(out3+i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[2]));
                *(out4+i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[3]));
                *(out5+i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[4]));
            }
            catch (...){

            }
        };
    }
    catch (CoolProp::HandleError &e){
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)){
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        }
        else{
            *errcode = 2;
        }
    }
    catch (...){
        *errcode = 3;
    }
}


/// *********************************************************************************
/// *********************************************************************************
///                     EMSCRIPTEN (for javascript)
/// *********************************************************************************
/// *********************************************************************************

#ifdef EMSCRIPTEN

#include <emscripten/bind.h>
using namespace emscripten;

// Binding code
EMSCRIPTEN_BINDINGS(coolprop_lib_bindings) {
    function("F2K", &F2K);
}

#endif
