#if defined(_MSC_VER)
#    define _CRTDBG_MAP_ALLOC
#    define _CRT_SECURE_NO_WARNINGS
#    include <crtdbg.h>
#else
#    include <fenv.h>
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
#include "Configuration.h"
#include "Backends/Helmholtz/MixtureParameters.h"

#include <string.h>

void str2buf(const std::string& str, char* buf, int n) {
    if (str.size() < static_cast<unsigned int>(n))
        strcpy(buf, str.c_str());
    else
        throw CoolProp::ValueError("Buffer size is too small");
}
void HandleException(long* errcode, char* message_buffer, const long buffer_length) {
    try {
        throw;  // Rethrow the error, and here we handle the error
    } catch (CoolProp::HandleError& e) {
        std::string errmsg = std::string("HandleError: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)) {
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        } else {
            *errcode = 2;
        }
    } catch (CoolProp::CoolPropBaseError& e) {
        std::string errmsg = std::string("Error: ") + e.what();
        if (errmsg.size() < static_cast<std::size_t>(buffer_length)) {
            *errcode = 1;
            strcpy(message_buffer, errmsg.c_str());
        } else {
            *errcode = 2;
        }
    } catch (...) {
        *errcode = 3;
    }
}

// In Microsoft Excel, they seem to check the FPU exception bits and error out because of it.
// By calling the _clearfp(), we can reset these bits, and not get the error
// See also http://stackoverflow.com/questions/11685441/floating-point-error-when-calling-dll-function-from-vba/27336496#27336496
// See also http://stackoverflow.com/questions/16849009/in-linux-do-there-exist-functions-similar-to-clearfp-and-statusfp for linux and OSX
struct fpu_reset_guard
{
    ~fpu_reset_guard() {
#if defined(_MSC_VER)
        _clearfp();  // For MSVC, clear the floating point error flags
#elif defined(FE_ALL_EXCEPT)
        feclearexcept(FE_ALL_EXCEPT);
#endif
    }
};
double convert_from_kSI_to_SI(long iInput, double value) {
    if (get_debug_level() > 8) {
        std::cout << format("%s:%d: convert_from_kSI_to_SI(i=%d,value=%g)\n", __FILE__, __LINE__, iInput, value).c_str();
    }

    switch (iInput) {
        case CoolProp::iP:
        case CoolProp::iCpmass:
        case CoolProp::iCp0mass:
        case CoolProp::iSmass:
        case CoolProp::iGmass:
        case CoolProp::iCvmass:
        case CoolProp::iHmass:
        case CoolProp::iUmass:
        case CoolProp::iconductivity:
            return value * 1000.0;
        case CoolProp::iDmass:
        case CoolProp::ispeed_sound:
        case CoolProp::iQ:
        case CoolProp::iviscosity:
        case CoolProp::iT:
        case CoolProp::iPrandtl:
        case CoolProp::isurface_tension:
            return value;
        default:
            throw CoolProp::ValueError(format("index [%d] is invalid in convert_from_kSI_to_SI", iInput).c_str());
    }
}

double convert_from_SI_to_kSI(long iInput, double value) {
    if (get_debug_level() > 8) {
        std::cout << format("%s:%d: convert_from_SI_to_kSI(%d,%g)\n", __FILE__, __LINE__, iInput, value).c_str();
    }

    switch (iInput) {
        case CoolProp::iP:
        case CoolProp::iCpmass:
        case CoolProp::iCp0mass:
        case CoolProp::iSmass:
        case CoolProp::iGmass:
        case CoolProp::iCvmass:
        case CoolProp::iHmass:
        case CoolProp::iUmass:
        case CoolProp::iconductivity:
            return value / 1000.0;
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

EXPORT_CODE long CONVENTION redirect_stdout(const char* file) {
    FILE* fp = freopen(file, "a+", stdout);
    return (fp) ? 1 : 0;  // 0 = failure if redirection could not be done; original stdout is already closed
}
EXPORT_CODE int CONVENTION set_reference_stateS(const char* Ref, const char* reference_state) {
    fpu_reset_guard guard;
    try {
        CoolProp::set_reference_stateS(std::string(Ref), std::string(reference_state));
        return true;
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return false;
}
EXPORT_CODE int CONVENTION set_reference_stateD(const char* Ref, double T, double rhomolar, double hmolar0, double smolar0) {
    fpu_reset_guard guard;
    try {
        CoolProp::set_reference_stateD(std::string(Ref), T, rhomolar, hmolar0, smolar0);
        return true;
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return false;
}

// All the function interfaces that point to the single-input Props function
EXPORT_CODE double CONVENTION Props1(const char* FluidName, const char* Output) {
    fpu_reset_guard guard;
    double val = Props1SI(Output, FluidName);
    if (!ValidNumber(val))
        // Error code was already set in Props1SI
        return val;
    // val is valid; so, Output is already checked in Props1SI -> get_parameter_index won't throw
    CoolProp::parameters iOutput = CoolProp::get_parameter_index(Output);
    return convert_from_SI_to_kSI(iOutput, val);
}
EXPORT_CODE double CONVENTION PropsS(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Ref) {
    return Props(Output, Name1[0], Prop1, Name2[0], Prop2, Ref);
}
EXPORT_CODE double CONVENTION Props(const char* Output, const char Name1, double Prop1, const char Name2, double Prop2, const char* Ref) {
    fpu_reset_guard guard;
    try {
        // Get parameter indices
        std::string sName1 = std::string(1, Name1), sName2 = std::string(1, Name2);
        CoolProp::parameters iOutput = CoolProp::get_parameter_index(Output);
        if (!CoolProp::is_trivial_parameter(iOutput)) {
            CoolProp::parameters iName1 = CoolProp::get_parameter_index(sName1);
            CoolProp::parameters iName2 = CoolProp::get_parameter_index(sName2);

            // Convert inputs to SI
            Prop1 = convert_from_kSI_to_SI(iName1, Prop1);
            Prop2 = convert_from_kSI_to_SI(iName2, Prop2);
        }

        // Call the SI function
        double val = PropsSI(Output, sName1.c_str(), Prop1, sName2.c_str(), Prop2, Ref);

        // Convert back to unit system
        return convert_from_SI_to_kSI(iOutput, val);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return _HUGE;
}
EXPORT_CODE double CONVENTION saturation_ancillary(const char* fluid_name, const char* output, int Q, const char* input, double value) {
    fpu_reset_guard guard;
    try {
        return CoolProp::saturation_ancillary(fluid_name, std::string(output), Q, std::string(input), value);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return _HUGE;
}
EXPORT_CODE double CONVENTION Props1SI(const char* FluidName, const char* Output) {
    fpu_reset_guard guard;
    return CoolProp::Props1SI(std::string(FluidName), std::string(Output));
}
EXPORT_CODE void CONVENTION Props1SImulti(const char* Outputs, char* backend, const char* FluidNames, const double* fractions,
                                          const long length_fractions, double* result, long* resdim1) {
    fpu_reset_guard guard;
    try {
        // Outputs is a delimited string separated by LIST_STRING_DELIMITER
        std::string delim = CoolProp::get_config_string(LIST_STRING_DELIMITER);
        // strsplit only support char delimiter
        if (delim.length() > 1)
            throw CoolProp::ValueError(format("Length of string delimiter [%d] is bigger than 1 [%d]", delim.length(), delim.size()));
        std::vector<std::string> _outputs = strsplit(Outputs, delim[0]);
        // FluidNames is a delimited string separated by LIST_STRING_DELIMITER
        std::vector<std::string> _fluidNames = strsplit(FluidNames, delim[0]);
        if (_fluidNames.size() != length_fractions)
            throw CoolProp::ValueError(
              format("Length of fractions vector  [%d] is not equal to length of fluidNames vector [%d]", _fluidNames.size(), length_fractions));
        std::vector<double> _fractions(fractions, fractions + length_fractions);
        std::vector<std::vector<double>> _result = CoolProp::Props1SImulti(_outputs, backend, _fluidNames, _fractions);
        // if CoolProp::Props1SImulti fails it will return an empty vector -> set result dimensions to 0
        if (_result.size() == 0) {
            *resdim1 = 0;
        } else {
            if (_result.size() > *resdim1)
                throw CoolProp::ValueError(format("Result vector [%d] is bigger than allocated memory [%d]", _result[0].size(), *resdim1));
            *resdim1 = _result[0].size();
            for (int i = 0; i < _result[0].size(); i++) {
                result[i] = _result[0][i];
            }
        }
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
}
EXPORT_CODE double CONVENTION PropsSI(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* FluidName) {
    fpu_reset_guard guard;
    return CoolProp::PropsSI(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(FluidName));
}
EXPORT_CODE void CONVENTION PropsSImulti(const char* Outputs, const char* Name1, double* Prop1, const long size_Prop1, const char* Name2,
                                         double* Prop2, const long size_Prop2, char* backend, const char* FluidNames, const double* fractions,
                                         const long length_fractions, double* result, long* resdim1, long* resdim2) {
    fpu_reset_guard guard;
    try {
        // Outputs is a delimited string separated by LIST_STRING_DELIMITER
        std::string delim = CoolProp::get_config_string(LIST_STRING_DELIMITER);
        // strsplit only support char delimiter
        if (delim.length() > 1)
            throw CoolProp::ValueError(format("Length of string delimiter [%d] is bigger than 1 [%d]", delim.length(), delim.size()));
        std::vector<std::string> _outputs = strsplit(Outputs, delim[0]);
        if (size_Prop1 != size_Prop2)
            throw CoolProp::ValueError(
              format("Length of input parameter 1 [%d] is not equal to length of input parameter 2 [%d]", size_Prop1, size_Prop2));
        // make vectors out of double pointer
        std::vector<double> _prop1(Prop1, Prop1 + size_Prop1);
        std::vector<double> _prop2(Prop2, Prop2 + size_Prop2);
        // FluidNames is a delimited string separated by LIST_STRING_DELIMITER
        std::vector<std::string> _fluidNames = strsplit(FluidNames, delim[0]);
        if (_fluidNames.size() != length_fractions)
            throw CoolProp::ValueError(
              format("Length of fractions vector  [%d] is not equal to length of fluidNames vector [%d]", _fluidNames.size(), length_fractions));
        std::vector<double> _fractions(fractions, fractions + length_fractions);
        std::vector<std::vector<double>> _result =
          CoolProp::PropsSImulti(_outputs, std::string(Name1), _prop1, std::string(Name2), _prop2, backend, _fluidNames, _fractions);
        // if CoolProp::PropsSImulti fails it will return an empty vector -> set result dimensions to 0
        if (_result.size() == 0) {
            *resdim1 = 0;
            *resdim2 = 0;
        } else {
            if (_result.size() > *resdim1 || _result[0].size() > *resdim2)
                throw CoolProp::ValueError(
                  format("Result matrix [%d x %d] is bigger than allocated memory [%d x %d]", _result.size(), _result[0].size(), *resdim1, *resdim2));
            *resdim1 = _result.size();
            *resdim2 = _result[0].size();
            for (int i = 0; i < _result.size(); i++) {
                for (int j = 0; j < _result[i].size(); j++) {
                    result[j + _result[i].size() * i] = _result[i][j];
                }
            }
        }
    } catch (std::exception& e) {
        // set result dimensions to 0 to signalize, an error occured
        *resdim1 = 0;
        *resdim2 = 0;
        CoolProp::set_error_string(e.what());
    } catch (...) {
        // set result dimensions to 0 to signalize, an error occured
        *resdim1 = 0;
        *resdim2 = 0;
        CoolProp::set_error_string("Undefined error");
    }
}
EXPORT_CODE long CONVENTION PhaseSI(const char* Name1, double Prop1, const char* Name2, double Prop2, const char* FluidName, char* phase, int n) {
    fpu_reset_guard guard;
    try {
        std::string s = CoolProp::PhaseSI(std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(FluidName));
        str2buf(s, phase, n);
        return 1;
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return 0;
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
EXPORT_CODE void CONVENTION propssi_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                     const char* FluidName, double* output) {
    *output = PropsSI(Output, Name1, *Prop1, Name2, *Prop2, FluidName);
}

EXPORT_CODE double CONVENTION K2F(double T) {
    return T * 9 / 5 - 459.67;
}
EXPORT_CODE double CONVENTION F2K(double T_F) {
    return (T_F + 459.67) * 5 / 9;
}
EXPORT_CODE int CONVENTION get_debug_level() {
    return CoolProp::get_debug_level();
}
EXPORT_CODE void CONVENTION set_debug_level(int level) {
    CoolProp::set_debug_level(level);
}
EXPORT_CODE long CONVENTION get_param_index(const char* param) {
    try {
        return CoolProp::get_parameter_index(param);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return -1;
}
EXPORT_CODE long CONVENTION get_input_pair_index(const char* pair) {
    try {
        return CoolProp::get_input_pair_index(pair);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return -1;
}
EXPORT_CODE long CONVENTION get_global_param_string(const char* param, char* Output, int n) {
    try {
        std::string s = CoolProp::get_global_param_string(param);
        str2buf(s, Output, n);
        return 1;
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return 0;
}
EXPORT_CODE long CONVENTION get_parameter_information_string(const char* param, char* Output, int n) {
    try {
        int key = CoolProp::get_parameter_index(param);
        std::string s = CoolProp::get_parameter_information(key, Output);
        str2buf(s, Output, n);
        return 1;
    } catch (std::exception& e) {
        // if param is wrong, CoolProp::get_parameter_index throws string like
        // "Your input name [%s] is not valid in get_parameter_index (names are case sensitive)"
        // CoolProp::get_parameter_information throws string like
        // "Bad info string [%s] to get_parameter_information" (if Output is wrong)
        // or "Unable to match the key [%d] in get_parameter_information for info [%s]"
        // (see src/DataStructures.cpp)
        // if n is too small, str2buf throws string
        // "Buffer size is too small"
        CoolProp::set_error_string(format("get_parameter_information_string(\"%s\", \"%s\", %d): %s", param, Output, n, e.what()));
    } catch (...) {
        CoolProp::set_error_string(format("get_parameter_information_string(\"%s\", \"%s\", %d): Undefined error", param, Output, n));
    }
    return 0;
}
EXPORT_CODE long CONVENTION get_fluid_param_string(const char* fluid, const char* param, char* Output, int n) {
    try {
        std::string s = CoolProp::get_fluid_param_string(std::string(fluid), std::string(param));
        str2buf(s, Output, n);
        return 1;
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return 0;
}
EXPORT_CODE void CONVENTION set_config_string(const char* key, const char* val) {
    try {
        CoolProp::set_config_string(CoolProp::config_string_to_key(std::string(key)), std::string(val));
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
}
EXPORT_CODE void CONVENTION set_config_double(const char* key, const double val) {
    try {
        CoolProp::set_config_double(CoolProp::config_string_to_key(std::string(key)), val);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
}
EXPORT_CODE void CONVENTION set_config_bool(const char* key, const bool val) {
    try {
        CoolProp::set_config_bool(CoolProp::config_string_to_key(std::string(key)), val);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
}
EXPORT_CODE void CONVENTION set_departure_functions(const char* string_data, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        CoolProp::set_departure_functions(string_data);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE double CONVENTION HAPropsSI(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Name3,
                                        double Prop3) {
    fpu_reset_guard guard;
    return HumidAir::HAPropsSI(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(Name3), Prop3);
}
EXPORT_CODE double CONVENTION cair_sat(double T) {
    fpu_reset_guard guard;
    return HumidAir::cair_sat(T);
}
EXPORT_CODE void CONVENTION hapropssi_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                       const char* Name3, const double* Prop3, double* output) {
    *output = HAPropsSI(Output, Name1, *Prop1, Name2, *Prop2, Name3, *Prop3);
}
EXPORT_CODE double CONVENTION HAProps(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Name3,
                                      double Prop3) {
    fpu_reset_guard guard;
    try {
        return HumidAir::HAProps(std::string(Output), std::string(Name1), Prop1, std::string(Name2), Prop2, std::string(Name3), Prop3);
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }
    return _HUGE;
}
EXPORT_CODE void CONVENTION haprops_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                     const char* Name3, const double* Prop3, double* output) {
    *output = HAProps(Output, Name1, *Prop1, Name2, *Prop2, Name3, *Prop3);
}
class AbstractStateLibrary
{
   private:
    std::map<std::size_t, shared_ptr<CoolProp::AbstractState>> ASlibrary;
    long next_handle;

   public:
    AbstractStateLibrary() : next_handle(0){};
    long add(shared_ptr<CoolProp::AbstractState> AS) {
        ASlibrary.insert(std::pair<std::size_t, shared_ptr<CoolProp::AbstractState>>(this->next_handle, AS));
        this->next_handle++;
        return next_handle - 1;
    }
    void remove(long handle) {
        std::size_t count_removed = ASlibrary.erase(handle);
        if (count_removed != 1) {
            throw CoolProp::HandleError("could not free handle");
        }
    }
    shared_ptr<CoolProp::AbstractState>& get(long handle) {
        std::map<std::size_t, shared_ptr<CoolProp::AbstractState>>::iterator it = ASlibrary.find(handle);
        if (it != ASlibrary.end()) {
            return it->second;
        } else {
            throw CoolProp::HandleError("could not get handle");
        }
    }
};
static AbstractStateLibrary handle_manager;

EXPORT_CODE long CONVENTION AbstractState_factory(const char* backend, const char* fluids, long* errcode, char* message_buffer,
                                                  const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, fluids));
        return handle_manager.add(AS);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
    return -1;
}
EXPORT_CODE void CONVENTION AbstractState_fluid_names(const long handle, char* fluids, long* errcode, char* message_buffer,
                                                      const long buffer_length) {
    *errcode = 0;

    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        std::vector<std::string> _fluids = AS->fluid_names();
        std::string fluidsstring = strjoin(_fluids, CoolProp::get_config_string(LIST_STRING_DELIMITER));
        if (fluidsstring.size() < static_cast<std::size_t>(buffer_length)) {
            strcpy(fluids, fluidsstring.c_str());
        } else {
            throw CoolProp::ValueError(format("Length of string [%d] is greater than allocated buffer length [%d]", fluidsstring.size(),
                                              static_cast<std::size_t>(buffer_length)));
        }

    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_free(const long handle, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        handle_manager.remove(handle);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_set_fractions(const long handle, const double* fractions, const long N, long* errcode, char* message_buffer,
                                                        const long buffer_length) {
    *errcode = 0;
    std::vector<double> _fractions(fractions, fractions + N);
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        if (AS->using_mole_fractions()) {
            AS->set_mole_fractions(_fractions);
        } else if (AS->using_mass_fractions()) {
            AS->set_mass_fractions(_fractions);
        } else if (AS->using_volu_fractions()) {
            AS->set_volu_fractions(_fractions);
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_get_mole_fractions(const long handle, double* fractions, const long maxN, long* N, long* errcode,
                                                             char* message_buffer, const long buffer_length) {
    *errcode = 0;

    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        std::vector<double> _fractions = AS->get_mole_fractions();
        *N = _fractions.size();
        if (*N <= maxN) {
            for (int i = 0; i < *N; i++)
                fractions[i] = _fractions[i];
        } else {
            throw CoolProp::ValueError(format("Length of array [%d] is greater than allocated buffer length [%d]", *N, maxN));
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_get_mole_fractions_satState(const long handle, const char* saturated_state, double* fractions,
                                                                      const long maxN, long* N, long* errcode, char* message_buffer,
                                                                      const long buffer_length) {
    *errcode = 0;

    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        std::vector<double> _fractions;
        double quality = AS->Q();
        std::string string_state(saturated_state);
        if (0 <= quality && quality <= 1) {
            if (string_state == "liquid") {
                _fractions = AS->mole_fractions_liquid();
            } else if (string_state == "gas") {
                _fractions = AS->mole_fractions_vapor();
            } else {
                throw CoolProp::ValueError(
                  format("Bad info string [%s] to saturated state mole fractions, options are \"liquid\" and \"gas\"", saturated_state));
            }
        } else {
            throw CoolProp::ValueError(format("AbstractState_get_mole_fractions_satState only returns outputs for saturated states if AbstractState "
                                              "quality [%g] is within two-phase region (0 <= quality <= 1)",
                                              static_cast<double>(quality)));
        }
        *N = _fractions.size();
        if (*N <= maxN) {
            for (int i = 0; i < *N; i++) {
                fractions[i] = _fractions[i];
            }
        } else {
            throw CoolProp::ValueError(format("Length of array [%d] is greater than allocated buffer length [%d]", *N, maxN));
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_update(const long handle, const long input_pair, const double value1, const double value2, long* errcode,
                                                 char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->update(static_cast<CoolProp::input_pairs>(input_pair), value1, value2);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_specify_phase(const long handle, const char* phase, long* errcode, char* message_buffer,
                                                        const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        return AS->specify_phase(CoolProp::get_phase_index(std::string(phase)));
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE void CONVENTION AbstractState_unspecify_phase(const long handle, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        return AS->unspecify_phase();
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
EXPORT_CODE double CONVENTION AbstractState_keyed_output(const long handle, const long param, long* errcode, char* message_buffer,
                                                         const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        return AS->keyed_output(static_cast<CoolProp::parameters>(param));
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
    return _HUGE;
}

EXPORT_CODE double CONVENTION AbstractState_first_saturation_deriv(const long handle, const long Of, const long Wrt, long* errcode,
                                                                   char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        return AS->first_saturation_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt));
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
    return _HUGE;
}

EXPORT_CODE double CONVENTION AbstractState_first_partial_deriv(const long handle, const long Of, const long Wrt, const long Constant, long* errcode,
                                                                char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        return AS->first_partial_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt),
                                       static_cast<CoolProp::parameters>(Constant));
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
    return _HUGE;
}

EXPORT_CODE void CONVENTION AbstractState_update_and_common_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                                const long length, double* T, double* p, double* rhomolar, double* hmolar,
                                                                double* smolar, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);

        for (int i = 0; i < length; i++) {
            try {
                AS->update(static_cast<CoolProp::input_pairs>(input_pair), *(value1 + i), *(value2 + i));
                *(T + i) = AS->T();
                *(p + i) = AS->p();
                *(rhomolar + i) = AS->rhomolar();
                *(hmolar + i) = AS->hmolar();
                *(smolar + i) = AS->smolar();
            } catch (...) {
            }
        };
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_update_and_1_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                           const long length, const long output, double* out, long* errcode, char* message_buffer,
                                                           const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);

        for (int i = 0; i < length; i++) {
            try {
                AS->update(static_cast<CoolProp::input_pairs>(input_pair), *(value1 + i), *(value2 + i));
                *(out + i) = AS->keyed_output(static_cast<CoolProp::parameters>(output));
            } catch (...) {
            }
        };
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_update_and_5_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                           const long length, long* outputs, double* out1, double* out2, double* out3, double* out4,
                                                           double* out5, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);

        for (int i = 0; i < length; i++) {
            try {
                AS->update(static_cast<CoolProp::input_pairs>(input_pair), *(value1 + i), *(value2 + i));
                *(out1 + i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[0]));
                *(out2 + i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[1]));
                *(out3 + i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[2]));
                *(out4 + i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[3]));
                *(out5 + i) = AS->keyed_output(static_cast<CoolProp::parameters>(outputs[4]));
            } catch (...) {
            }
        };
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_set_binary_interaction_double(const long handle, const long i, const long j, const char* parameter,
                                                                        const double value, long* errcode, char* message_buffer,
                                                                        const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->set_binary_interaction_double(static_cast<std::size_t>(i), static_cast<std::size_t>(j), parameter, value);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_set_cubic_alpha_C(const long handle, const long i, const char* parameter, const double c1, const double c2,
                                                            const double c3, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->set_cubic_alpha_C(static_cast<std::size_t>(i), parameter, c1, c2, c3);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_set_fluid_parameter_double(const long handle, const long i, const char* parameter, const double value,
                                                                     long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->set_fluid_parameter_double(static_cast<std::size_t>(i), parameter, value);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_build_phase_envelope(const long handle, const char* level, long* errcode, char* message_buffer,
                                                               const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->build_phase_envelope(level);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_get_phase_envelope_data(const long handle, const long length, double* T, double* p, double* rhomolar_vap,
                                                                  double* rhomolar_liq, double* x, double* y, long* errcode, char* message_buffer,
                                                                  const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        CoolProp::PhaseEnvelopeData pe = AS->get_phase_envelope_data();
        if (pe.T.size() > static_cast<std::size_t>(length)) {
            throw CoolProp::ValueError(format("Length of phase envelope vectors [%d] is greater than allocated buffer length [%d]",
                                              static_cast<int>(pe.T.size()), static_cast<int>(length)));
        }
        std::size_t N = pe.x.size();
        for (std::size_t i = 0; i < pe.T.size(); i++) {
            *(T + i) = pe.T[i];
            *(p + i) = pe.p[i];
            *(rhomolar_vap + i) = pe.rhomolar_vap[i];
            *(rhomolar_liq + i) = pe.rhomolar_liq[i];
            for (std::size_t j = 0; j < N; ++j) {
                *(x + i * N + j) = pe.x[j][i];
                *(y + i * N + j) = pe.y[j][i];
            }
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_get_phase_envelope_data_checkedMemory(const long handle, const long length, const long maxComponents, double* T,
                                                                  double* p, double* rhomolar_vap, double* rhomolar_liq, double* x, double* y,
                                                                  long* actual_length, long* actual_components, long* errcode, char* message_buffer,
                                                                  const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        CoolProp::PhaseEnvelopeData pe = AS->get_phase_envelope_data();
        *actual_length = pe.T.size();
        if (pe.T.size() > static_cast<std::size_t>(length)) {
            throw CoolProp::ValueError(format("Length of phase envelope vectors [%d] is greater than allocated buffer length [%d]",
                                              static_cast<int>(pe.T.size()), static_cast<int>(length)));
        }
        *actual_components = pe.x.size();
        if (*actual_components > static_cast<std::size_t>(maxComponents)) {
            throw CoolProp::ValueError(format("Length of phase envelope composition vectors [%d] is greater than allocated buffer length [%d]",
                                              static_cast<int>(*actual_components), static_cast<int>(maxComponents)));
        }
        for (std::size_t i = 0; i < pe.T.size(); i++) {
            *(T + i) = pe.T[i];
            *(p + i) = pe.p[i];
            *(rhomolar_vap + i) = pe.rhomolar_vap[i];
            *(rhomolar_liq + i) = pe.rhomolar_liq[i];
            for (std::size_t j = 0; j < *actual_components; ++j) {
                *(x + i * *actual_components + j) = pe.x[j][i];
                *(y + i * *actual_components + j) = pe.y[j][i];
            }
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_build_spinodal(const long handle, long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        AS->build_spinodal();
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_get_spinodal_data(const long handle, const long length, double* tau, double* delta, double* M1,
                                                            long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        CoolProp::SpinodalData spin = AS->get_spinodal_data();
        if (spin.tau.size() > static_cast<std::size_t>(length)) {
            throw CoolProp::ValueError(format("Length of spinodal vectors [%d] is greater than allocated buffer length [%d]",
                                              static_cast<int>(spin.tau.size()), static_cast<int>(length)));
        }
        for (std::size_t i = 0; i < spin.tau.size(); ++i) {
            *(tau + i) = spin.tau[i];
            *(delta + i) = spin.delta[i];
            *(M1 + i) = spin.M1[i];
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION AbstractState_all_critical_points(const long handle, long length, double* T, double* p, double* rhomolar, long* stable,
                                                              long* errcode, char* message_buffer, const long buffer_length) {
    *errcode = 0;
    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        std::vector<CoolProp::CriticalState> pts = AS->all_critical_points();
        if (pts.size() > static_cast<std::size_t>(length)) {
            throw CoolProp::ValueError(format("Length of critical point vector [%d] is greater than allocated buffer length [%d]",
                                              static_cast<int>(pts.size()), static_cast<int>(length)));
        }
        for (std::size_t i = 0; i < pts.size(); ++i) {
            *(T + i) = pts[i].T;
            *(p + i) = pts[i].p;
            *(rhomolar + i) = pts[i].rhomolar;
            *(stable + i) = pts[i].stable;
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE double CONVENTION AbstractState_keyed_output_satState(const long handle, const char* saturated_state, const long param, long* errcode,
                                                                  char* message_buffer, const long buffer_length) {
    *errcode = 0;

    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        double quality = AS->Q();
        std::string string_state(saturated_state);
        if (0 <= quality && quality <= 1) {
            if (string_state == "liquid") {
                return AS->saturated_liquid_keyed_output(static_cast<CoolProp::parameters>(param));
            } else if (string_state == "gas") {
                return AS->saturated_vapor_keyed_output(static_cast<CoolProp::parameters>(param));
            } else {
                throw CoolProp::ValueError(
                  format("Bad info string [%s] to saturated state output, options are \"liquid\" and \"gas\"", saturated_state));
            }
        } else {
            throw CoolProp::ValueError(format("AbstractState_keyed_output_satState only returns outputs for saturated states if AbstractState "
                                              "quality [%g] is within two-phase region (0 <= quality <= 1)",
                                              static_cast<double>(quality)));
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
        return _HUGE;
    }
}

EXPORT_CODE void CONVENTION AbstractState_backend_name(const long handle, char* backend, long* errcode, char* message_buffer,
                                                       const long buffer_length) {
    *errcode = 0;

    try {
        shared_ptr<CoolProp::AbstractState>& AS = handle_manager.get(handle);
        std::string backendstring = AS->backend_name();
        if (backendstring.size() < static_cast<std::size_t>(buffer_length)) {
            strcpy(backend, backendstring.c_str());
        } else {
            throw CoolProp::ValueError(format("Length of string [%d] is greater than allocated buffer length [%d]", backendstring.size(),
                                              static_cast<std::size_t>(buffer_length)));
        }
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}

EXPORT_CODE void CONVENTION add_fluids_as_JSON(const char* backend, const char* fluidstring, long* errcode, char* message_buffer,
                                               const long buffer_length) {
    *errcode = 0;

    try {
        CoolProp::add_fluids_as_JSON(backend, fluidstring);
    } catch (...) {
        HandleException(errcode, message_buffer, buffer_length);
    }
}
