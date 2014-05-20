#ifndef COOLPROPDLL_H
#define COOLPROPDLL_H

    #include "PlatformDetermination.h"

    #if defined(COOLPROP_LIB)
    #  ifndef EXPORT_CODE
    #    if defined(__ISWINDOWS__)
    #      define EXPORT_CODE extern "C" __declspec(dllexport)
    #    else
    #      define EXPORT_CODE extern "C"
    #    endif
    #  endif
    #  ifndef CONVENTION
    #    if defined(__ISWINDOWS__)
    #      define CONVENTION __stdcall
    #    else
    #      define CONVENTION
    #    endif
    #  endif
    #else
    #  ifndef EXPORT_CODE
    #    define EXPORT_CODE
    #  endif
    #  ifndef CONVENTION
    #    define CONVENTION
    #  endif
    #endif

    // Hack for PowerPC compilation to only use extern "C"
    #if defined(__powerpc__) || defined(EXTERNC)
    #  undef EXPORT_CODE
    #  define EXPORT_CODE extern "C"
    #endif
    
    // Functions with the call type like
    // EXPORT_CODE void CONVENTION AFunction(double, double);
    // will be exported to the DLL

    /*
    ####################################################################################
    Overloads for DLL wrapping puposes

    These functions take strings which are 0-terminated.  These functions pass directly to 
    equivalently named functions in CoolProp.h that take std::string
    ####################################################################################
    */

    /**
    \overload 
    \sa Props1SI(std::string, std::string)
    */
    EXPORT_CODE double CONVENTION Props1SI(const char *FluidName, const char* Output);
    /**
    \overload 
    \sa PropsSI(std::string, std::string, double, std::string, double, std::string)
    */
    EXPORT_CODE double CONVENTION PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
    /**
    \overload 
    \sa Props(std::string, std::string, double, std::string, double, std::string)
    */
    EXPORT_CODE double CONVENTION PropsS(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
    /**
    \overload 
    \sa Props(std::string, std::string, double, std::string, double, std::string)
    */
    EXPORT_CODE double CONVENTION Props(const char *Output, const char Name1, double Prop1, const char Name2, double Prop2, const char *Ref);
    /**
    \overload 
    \sa Props1(std::string, std::string)
    */
    EXPORT_CODE double CONVENTION Props1(const char *FluidName, const char *Output);
    /**
    \overload
    \sa IsFluidType(std::string, std::string)
    */
    EXPORT_CODE int CONVENTION IsFluidType(const char *Ref, const char *Type);
    
    // When using SWIG, it is extremely difficult to deal with char* for output strings, so just use 
    // the std::string version since SWIG can handle std::string just fine
    #if defined(SWIG)
        std::string get_global_param_string(std::string ParamName);
        std::string get_fluid_param_string(std::string FluidName, std::string ParamName);
    #else
        EXPORT_CODE long CONVENTION get_global_param_string(const char *param, char *Output);
        EXPORT_CODE long CONVENTION get_fluid_param_string(const char *fluid, const char *param, char *Output);
    #endif

    /**
    \overload
    \sa set_reference_stateS(std::string, std::string)
    */
    EXPORT_CODE int CONVENTION set_reference_stateS(const char *Ref, const char *reference_state);
    /**
    \overload
    \sa set_reference_stateD(std::string, double, double, double, double)
    */
    EXPORT_CODE int CONVENTION set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
    /*
    ####################################################################################
    Implemented functions

    These functions take inputs that are compatible with DLL passing and are 
    implemented in CoolPropDLL.cpp
    ####################################################################################
    */

    /**
    \brief FORTRAN 77 style wrapper of the PropsSI function
    \sa PropsSI(std::string, std::string, double, std::string, double, std::string)
    */
    EXPORT_CODE void CONVENTION F77PropsSI(const char *Output, const char *Name1, double *Prop1, const char *Name2, double *Prop2, const char * Ref, double *output);
    /**
    
    */
    EXPORT_CODE double CONVENTION PropsSIZ(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *FluidName, const double *z, int n);

    // This version uses the indices in place of the strings for speed.  Get the parameter indices
    // from get_param_index('D') for instance and the Fluid index from get_Fluid_index('Air') for instance
    EXPORT_CODE double CONVENTION IPropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
    EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
    
    /// Convert from degrees Fahrenheit to Kelvin (useful primarily for testing)
    EXPORT_CODE double CONVENTION F2K(double T_F);
    /// Convert from Kelvin to degrees Fahrenheit (useful primarily for testing)
    EXPORT_CODE double CONVENTION K2F(double T_K);
    
    EXPORT_CODE long CONVENTION get_param_index(const char *param);
    EXPORT_CODE long CONVENTION redirect_stdout(const char *file);

    // ---------------------------------
    // Getter and setter for debug level
    // ---------------------------------

    /// Get the debug level
    /// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    EXPORT_CODE int CONVENTION get_debug_level();
    /// Set the debug level
    /// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    EXPORT_CODE void CONVENTION set_debug_level(int level);

    // ---------------------------------
    //        Humid Air Properties
    // ---------------------------------

    
    EXPORT_CODE double CONVENTION HAPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Name3, double Prop3);

    /** 
    \brief FORTRAN 77 style wrapper of the HAPropsSI function
    */
    EXPORT_CODE void CONVENTION F77HAPropsSI(const char *Output, const char *Name1, double *Prop1, const char *Name2, double *Prop2, const char *Name3, double *Prop3, double *output);

    

#endif
