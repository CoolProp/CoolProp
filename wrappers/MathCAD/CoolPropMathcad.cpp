// CoolPropMathcad.cpp : Defines the exported functions for the DLL Add-in.
//

#include <string>

#ifndef NOMINMAX // Kill windows' horrible min() and max() macros
#define NOMINMAX
#endif
#include "mcadincl.h"
#undef NOMINMAX; 

enum { MC_STRING = STRING };  // substitute enumeration variable MC_STRING for STRING, use MC_STRING below
#undef STRING                 // undefine STRING as it conflicts with STRING enum in cppformat/format.h

#include "CoolProp.h"
#include "HumidAirProp.h"

enum EC { INTERRUPTED, INSUFFICIENT_MEMORY, MUST_BE_REAL, NUMBER_OF_ERRORS = MUST_BE_REAL };   // Mathcad Error Codes

    // table of error messages
    // if user function never returns an
    // error -- you do not need to create this
    // table
    char * myErrorMessageTable[NUMBER_OF_ERRORS] =  
    {   
        "Interrupted"
        "Insufficient Memory"
        "Argument must be real"
    };
    
    // this code executes the user function CP_get_global_param_string, which is a wrapper for
    // the CoolProp.get_global_param_string() function, used to get a global string parameter from CoolProp
    LRESULT  CP_get_global_param_string(
                            LPMCSTRING ParamValue,   // output (value of parameter)
                            LPCMCSTRING ParamName )  // name of parameter (string) to retrieve
    {  
        // Invoke the std::string form of get_global_param_string() function, save result to a new string s
        std::string s = CoolProp::get_global_param_string(ParamName->str);
        char * c = new char [s.size()+1]; // creat a c-string (pointer) c with the same size as s
        // copy s into c, this process avoids the const-cast type which would result from instead
        // converting the string using s.c_str()
        std::copy(s.begin(), s.end(), c); 
        c[s.size()] = '\0';
        // assign the string to the function's output parameter
        ParamValue->str = c;

        // normal return
        return 0;
    }

    // this code executes the user function CP_get_fluid_param_string, which is a wrapper for
    // the CoolProp.get_fluid_param_string() function, used to get a fluid string parameter from CoolProp
    LRESULT  CP_get_fluid_param_string(
                            LPMCSTRING ParamValue,   // output (value of parameter)
                            LPCMCSTRING FluidName,  // name of fluidr (string) to retrieve
                            LPCMCSTRING ParamName )  // name of parameter (string) to retrieve
    {  
        // Invoke the std::string form of get_fluid_param_string() function, save result to a new string s
        std::string s = CoolProp::get_fluid_param_string(FluidName->str, ParamName->str);
        char * c = new char [s.size()+1]; // creat a c-string (pointer) c with the same size as s
        // copy s into c, this process avoids the const-cast type which would result from instead
        // converting the string using s.c_str()
        std::copy(s.begin(), s.end(), c); 
        c[s.size()] = '\0';
        // assign the string to the function's output parameter
        ParamValue->str = c;

        // normal return
        return 0;
    }

    // this code executes the user function CP_set_reference_state, which is a wrapper for
    // the CoolProp.set_reference_stateS() function, used to set the H/S reference states
    // based on a standard state string of "IIR", "ASHRAE", "NBP", or "DEF".   
    LRESULT  CP_set_reference_state(
                            LPCOMPLEXSCALAR  Conf,     // output (dummy value)
                            LPCMCSTRING FluidName,     // name of fluidr (string) to retrieve
                            LPCMCSTRING StateStr )     // name of standard state (string) to set
    {  
        // Invoke the set_reference_stateS() function, no result from this void function.
        CoolProp::set_reference_stateS(FluidName->str, StateStr->str);
        // assign the dummy return value
        Conf->real = 0;

        // normal return
        return 0;
    }

    // this code executes the user function CP_Props1SI, which is a wrapper for
    // the CoolProp.PropsSI() function, used to simply extract a
    // fluid-specific parameter that is not dependent on the state
    LRESULT  CP_Props1SI(
                            LPCOMPLEXSCALAR Prop,   // pointer to the result
                            LPCMCSTRING Fluid,      // string with a valid CoolProp fluid name
                            LPCMCSTRING PropName )  // a fluid property
    {  
        // pass the arguments to the CoolProp.Props1() function
        Prop->real = CoolProp::Props1SI(Fluid->str, PropName->str);

        // normal return
        return 0;
    }

    // this code executes the user function CP_PropsSI, which is a wrapper for
    // the CoolProp.PropsSI() function, used to extract a fluid-specific parameter that is dependent on the state
    LRESULT  CP_PropsSI(
                            LPCOMPLEXSCALAR  Prop,       // pointer to the result
                            LPCMCSTRING      OutputName, // string with a valid CoolProp OutputName
                            LPCMCSTRING      InputName1, // CoolProp InputName1
                            LPCCOMPLEXSCALAR InputProp1, // CoolProp InputProp1
                            LPCMCSTRING      InputName2, // CoolProp InputName2
                            LPCCOMPLEXSCALAR InputProp2, // CoolProp InputProp2
                            LPCMCSTRING      FluidName ) // CoolProp Fluid
    {  
        // check that the first scalar argument is real
        if (InputProp1->imag != 0.0)
            return MAKELRESULT( MUST_BE_REAL, 3);    // if not, display "Argument must be real" under scalar argument

        // check that the second scalar argument is real
        if (InputProp2->imag != 0.0)
            return MAKELRESULT( MUST_BE_REAL, 6);    // if not, display "Argument must be real" under scalar argument

        // pass the arguments to the CoolProp.Props() function
        Prop->real = CoolProp::PropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);
        
        // normal return
        return 0;
    }

    // this code executes the user function CP_HAPropsSI, which is a wrapper for
    // the CoolProp.HAPropsSI() function, used to extract humid air properties in base-SI units
    LRESULT  CP_HAPropsSI(
                            LPCOMPLEXSCALAR    Prop,       // pointer to the result
                            LPCMCSTRING        OutputName, // string with a valid CoolProp Output Name
                            LPCMCSTRING        InputName1, // CoolProp InputName1
                            LPCCOMPLEXSCALAR   InputProp1, // CoolProp InputProp1
                            LPCMCSTRING        InputName2, // CoolProp InputName2
                            LPCCOMPLEXSCALAR   InputProp2, // CoolProp InputProp2
                            LPCMCSTRING        InputName3, // CoolProp InputName3
                            LPCCOMPLEXSCALAR   InputProp3) // CoolProp InputProp3
    {  
        // check that the first scalar argument is real
        if (InputProp1->imag != 0.0)
            return MAKELRESULT( MUST_BE_REAL, 3);  // if not, display "must be real" under scalar argument

        // check that the second scalar argument is real 
        if (InputProp2->imag != 0.0)
            return MAKELRESULT( MUST_BE_REAL, 5);  // if not, display "must be real" under scalar argument

        // check that the third scalar argument is real 
        if (InputProp3->imag != 0.0)
            return MAKELRESULT( MUST_BE_REAL, 7);  // if not, display "must be real" under scalar argument

        // pass the arguments to the HumidAirProp.HAProps() function
        Prop->real = HumidAir::HAPropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, InputName3->str, InputProp3->real);
        
        // normal return
        return 0;
    }

    // fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
    FUNCTIONINFO PropsParam = 
    {
    "get_global_param_string", // Name by which MathCAD will recognize the function   
    "Name of the parameter to retrieve", // Description of input parameters
    "Returns the value of the requested CoolProps parameter", // description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CP_get_global_param_string, // Pointer to the function code. 
    MC_STRING, // Returns a MathCAD string
    1, // Number of arguments
    {MC_STRING} // Argument types 
    };	

    // fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
    FUNCTIONINFO FluidParam = 
    {
    "get_fluid_param_string", // Name by which MathCAD will recognize the function   
    "Fluid, Name of the parameter to retrieve", // Description of input parameters
    "Returns the value of the requested CoolProps parameter", // description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CP_get_fluid_param_string, // Pointer to the function code. 
    MC_STRING, // Returns a MathCAD string
    2, // Number of arguments
    {MC_STRING, MC_STRING} // Argument types 
    };

    // fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
    FUNCTIONINFO RefState = 
    {
    "set_reference_state", // Name by which MathCAD will recognize the function   
    "Fluid, Reference State String", // Description of input parameters
    "Sets the reference state to either IIR, ASHRAE, NBP, or DEF.", // description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CP_set_reference_state, // Pointer to the function code. 
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    2, // Number of arguments
    {MC_STRING, MC_STRING} // Argument types 
    };

    // fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
    FUNCTIONINFO Props1SI = 
    {
    "Props1SI", // Name by which MathCAD will recognize the function   
    "Fluid, Property Name", // Description of input parameters
    "Returns a fluid-specific parameter, where the parameter is not dependent on the fluid state", // Description of the function for the Insert Function dialog box          
    (LPCFUNCTION)CP_Props1SI, // Pointer to the function code.
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    2,    // Number of arguments
    {MC_STRING, MC_STRING}  // Argument types 
    };

    // fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
    FUNCTIONINFO PropsSI = 
    {
    "PropsSI", // Name by which MathCAD will recognize the function
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name", // Description of input parameters
    "Returns a fluid-specific parameter, where the parameter is dependent on the fluid state", // Description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CP_PropsSI,  // Pointer to the function code.  
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    6, // Number of arguments
    {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING} // Argument types 
    };

    FUNCTIONINFO HAPropsSI = 
    {
    "HAPropsSI",  // Name by which MathCAD will recognize the function
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3",  // Description of input parameters
    "Returns a parameter of humid air, where the parameter is dependent on the fluid state", // Description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CP_HAPropsSI,  // Pointer to the function code.
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    7,  // Number of arguments
    {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR} // Argument types 
    };

    // ************************************************************************************
    // DLL entry point code.  
    // ************************************************************************************
    // The _CRT_INIT function is needed if you are using Microsoft's 32 bit compiler
#ifdef _WIN32 
    extern "C" BOOL WINAPI _CRT_INIT(HINSTANCE hinstDLL, DWORD dwReason, LPVOID lpReserved);
#endif

    extern "C" BOOL WINAPI  DllEntryPoint (HINSTANCE hDLL, DWORD dwReason, LPVOID lpReserved)
    {
        switch (dwReason)
        {
                case DLL_PROCESS_ATTACH:
                    //
                    // DLL is attaching to the address space of 
                    // the current process.
                    //
                    if (!_CRT_INIT(hDLL, dwReason, lpReserved)) 
                        return FALSE;
                
                    // register the error message table
                    // Note, that if your function never returns
                    // an error -- you do not need to 
                    // register an error message table
                    if ( !CreateUserErrorMessageTable( hDLL, NUMBER_OF_ERRORS, myErrorMessageTable ) )
                        break;

                    // and if the errors register OK
                    // go ahead and
                    // register user function
                    CreateUserFunction( hDLL, &PropsParam);
                    CreateUserFunction( hDLL, &FluidParam);
                    CreateUserFunction( hDLL, &RefState);
                    CreateUserFunction( hDLL, &Props1SI );
                    CreateUserFunction( hDLL, &PropsSI );
                    CreateUserFunction( hDLL, &HAPropsSI );
                    break;

                case DLL_THREAD_ATTACH:
                case DLL_THREAD_DETACH:
                case DLL_PROCESS_DETACH:

                    if (!_CRT_INIT(hDLL, dwReason, lpReserved))
                    {
                        Sleep(1000);   // Attempt to keep CRT_INIT from detaching before all threads are closed
                        return FALSE;
                    }
                    break;                   
        }
        return TRUE;
    }


    

