// CoolPropMathcad.c : Defines the exported functions for the DLL application.
//

#include <string>

#include "mcadincl.h"
// Why does MathCAD use min and max macros???
#undef min
#undef max

#include "CoolProp.h"
#include "HumidAirProp.h"

#define  MUST_BE_REAL           1
#define  NUMBER_OF_ERRORS       1   

    // table of error messages
    // if user function never returns an
    // error -- you do not need to create this
    // table
    char * myErrorMessageTable[NUMBER_OF_ERRORS] =  
    {   
        "must be real"
    };
    
    // this code executes the user function CoolProp_get_global_param_string, which is a wrapper for
	// the CoolProp.get_global_param_string() function, used to get a string parameter from CoolProp
    LRESULT  CoolProp_get_global_param_string(
        MCSTRING * const ParamValue, // output (value of parameter)
		const MCSTRING * const ParamName ) // name of parameter (string) to retrieve
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

    // this code executes the user function CoolProp_Props1SI, which is a wrapper for
	// the PropsSI() function, used to simply extract a
	// fluid-specific parameter that is not dependent on the state
    LRESULT  CoolProp_Props1SI(
		COMPLEXSCALAR * const Prop, // pointer to the result
        const MCSTRING * const Fluid, // string with a valid CoolProp fluid name
		const MCSTRING * const PropName ) // a fluid property
    {  
        // pass the arguments to the CoolProp.Props1() function
		Prop->real = CoolProp::Props1SI(Fluid->str, PropName->str);

		// normal return
        return 0;
    }

	// this code executes the user function CoolProp_PropsSI, which is a wrapper for
	// the PropsSI() function, used to extract a fluid-specific parameter that is dependent on the state
    LRESULT  CoolProp_PropsSI(
		COMPLEXSCALAR * const Prop, // pointer to the result
        const MCSTRING * const OutputName, // string with a valid CoolProp OutputName
		const MCSTRING * const InputName1, // CoolProp InputName1
		const COMPLEXSCALAR * const InputProp1, // CoolProp InputProp1
		const MCSTRING * const InputName2, // CoolProp InputName2
		const COMPLEXSCALAR * const InputProp2, // CoolProp InputProp2
		const MCSTRING * const FluidName ) // CoolProp Fluid
    {  
		// check that the scalar arguments are real - if not, return an error and highlight the invalid argument in Mathcad
		if (InputProp1->imag != 0.0)
			// if not, display "must be real" under scalar argument
			return MAKELRESULT( MUST_BE_REAL, 3);
		// check that the scalar arguments are real - if not, return an error and highlight the invalid argument in Mathcad
		if (InputProp2->imag != 0.0)
			// if not, display "must be real" under scalar argument
			return MAKELRESULT( MUST_BE_REAL, 6);

        // pass the arguments to the CoolProp.Props() function
		Prop->real = CoolProp::PropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);
		
		// normal return
        return 0;
    }

	// this code executes the user function CoolProp_HAPropsSI, which is a wrapper for
	// the HAPropsSI() function, used to extract humid air properties in base-SI units
    LRESULT  CoolProp_HAPropsSI(
		COMPLEXSCALAR * const Prop, // pointer to the result
		const MCSTRING * const OutputName, // pointer to the result
		const MCSTRING * const InputName1, // CoolProp InputName1
		const COMPLEXSCALAR * const InputProp1, // CoolProp InputProp1
		const MCSTRING * const InputName2, // CoolProp InputName2
		const COMPLEXSCALAR * const InputProp2, // CoolProp InputProp2
		const MCSTRING * const InputName3, // CoolProp InputName3
		const COMPLEXSCALAR * const InputProp3) // CoolProp InputProp3
    {  
		// check that the scalar arguments are real - if not, return an error and highlight the invalid argument in Mathcad
		if (InputProp1->imag != 0.0)
			// if not, display "must be real" under scalar argument
			return MAKELRESULT( MUST_BE_REAL, 3);
		// check that the scalar arguments are real - if not, return an error and highlight the invalid argument in Mathcad
		if (InputProp2->imag != 0.0)
			// if not, display "must be real" under scalar argument
			return MAKELRESULT( MUST_BE_REAL, 5);
		// check that the scalar arguments are real - if not, return an error and highlight the invalid argument in Mathcad
		if (InputProp3->imag != 0.0)
			// if not, display "must be real" under scalar argument
			return MAKELRESULT( MUST_BE_REAL, 7);

        // pass the arguments to the HumidAirProp.HAProps() function
		Prop->real = HumidAir::HAProps(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, InputName3->str, InputProp3->real);
		
		// normal return
        return 0;
    }

	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO get_global_param_string = 
    {
    "get_global_param_string", // Name by which MathCAD will recognize the function   
    "string, name of the parameter to retrieve", // Description of input parameters
    "returns the value of the requested CoolProps parameter", // description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CoolProp_get_global_param_string, // Pointer to the function code. 
    STRING, // Returns a MathCAD string
    1, // Number of arguments
    {STRING} // Argument types 
    };	
	
	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO Props1SI = 
    {
    "PropsSI1", // Name by which MathCAD will recognize the function   
    "Fluid, Property Name", // Description of input parameters
    "returns a fluid-specific parameter, where the parameter is not dependent on the fluid state", // Description of the function for the Insert Function dialog box          
    (LPCFUNCTION)CoolProp_Props1SI, // Pointer to the function code.
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    2,    // Number of arguments
    {STRING, STRING}  // Argument types 
    };

	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO PropsSI = 
    {
    "PropsSI", // Name by which MathCAD will recognize the function
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name", // Description of input parameters
    "Returns a fluid-specific parameter, where the parameter is dependent on the fluid state", // Description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CoolProp_PropsSI,  // Pointer to the function code.  
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    6, // Number of arguments
    {STRING, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR, STRING} // Argument types 
    };

	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO HAPropsSI = 
    {
    "HAPropsSI",  // Name by which MathCAD will recognize the function
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3",  // Description of input parameters
    "Returns a parameter of humid air, where the parameter is dependent on the fluid state", // Description of the function for the Insert Function dialog box       
    (LPCFUNCTION)CoolProp_HAPropsSI,  // Pointer to the function code.
    COMPLEX_SCALAR, // Returns a MathCAD complex scalar
    7,  // Number of arguments
    {STRING, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR} // Argument types 
    };

    // DLL entry point code.  the _CRT_INIT function is needed
	// if you are using Microsoft's 32 bit compiler
 
    extern "C" BOOL WINAPI _CRT_INIT(HINSTANCE hinstDLL, DWORD dwReason, LPVOID lpReserved);

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
                    if ( CreateUserErrorMessageTable(
                            hDLL, NUMBER_OF_ERRORS, myErrorMessageTable ) )
                        // and if the errors register OK
                        // go ahead and
                        // register user function
                        CreateUserFunction( hDLL, &Props1SI );
						CreateUserFunction( hDLL, &PropsSI );
						CreateUserFunction( hDLL, &HAPropsSI );
						CreateUserFunction( hDLL, &get_global_param_string);
                    break;

                case DLL_THREAD_ATTACH:
                case DLL_THREAD_DETACH:
                case DLL_PROCESS_DETACH:

                    if (!_CRT_INIT(hDLL, dwReason, lpReserved)) 
                        return FALSE;
                    break;                   
        }
        return TRUE;
    }

#undef INTERRUPTED    
#undef INSUFFICIENT_MEMORY
#undef MUST_BE_REAL   
#undef NUMBER_OF_ERRORS     

    

