// CoolProp 4.0 Beta Mathcad 15 Wrapper.cpp : Defines the exported functions for the DLL application.
//

#define WIN32_LEAN_AND_MEAN


#include "CoolProp.h"
#include "FluidClass.h"
#include "HumidAirProp.h"
#include "string.h"
#include "MCADINCL.H"

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


	LRESULT  GetCoolPropsParamMathcad(LPMCSTRING ParamValue, LPCMCSTRING ParamName );
	LRESULT  CoolPropMathcad1(COMPLEXSCALAR * const Prop1, const MCSTRING * const Fluid, const MCSTRING * const PropName );
    LRESULT  CoolPropMathcad(COMPLEXSCALAR * const Prop, const MCSTRING * const OutputName, const MCSTRING * const InputName1, const COMPLEXSCALAR * const InputProp1, const MCSTRING * const InputName2, const COMPLEXSCALAR * const InputProp2, const MCSTRING * const FluidName );
    LRESULT  HumidAirPropMathcad(COMPLEXSCALAR * const Prop, const MCSTRING * const OutputName, const MCSTRING * const InputName1, const COMPLEXSCALAR * const InputProp1, const MCSTRING * const InputName2, const COMPLEXSCALAR * const InputProp2, const MCSTRING * const InputName3, const COMPLEXSCALAR * const InputProp3);
   	
	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO FluidPropsParams = 
    {
    // Name by which mathcad will recognize the function
    "FluidPropsParams",        
    
    // description of "GetCoolPropsParam" parameters to be used by the Insert Function dialog box
    "string, name of the parameter to retrieve",   
    
    // description of the function for the Insert Function dialog box       
    "returns the value of the requested CoolProps parameter",    
    
    // pointer to the executible code. i.e. code that should be executed
    // when a user types in "GetCoolPropsParam(ParamName)="
    (LPCFUNCTION)GetCoolPropsParamMathcad,  
        
    // FluidPropsParams(ParamName) returns a Mathcad string
    STRING,
        
    // FluidPropsParams(ParamName) takes on one argument
    1,   
    
    // the argument type is a Mathcad string 
    {STRING} 
    };	
	
	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO FluidProp1 = 
    {
    // Name by which mathcad will recognize the function
    "FluidProp1",        
    
    // description of "FluidProp1" parameters to be used by the Insert Function dialog box
    "Fluid, Property Name",   
    
    // description of the function for the Insert Function dialog box       
    "returns a fluid-specific parameter, where the parameter is not dependent on the fluid state",    
    
    // pointer to the executible code. i.e. code that should be executed
    // when a user types in "FluidProp1(Fluid,PropName)="
    (LPCFUNCTION)CoolPropMathcad1,  
        
    // FluidProp1(Fluid,PropName) returns a Mathcad complex scalar
    COMPLEX_SCALAR,
        
    // FluidProp1(Fluid,PropName) takes on two arguments
    2,   
    
    // the first is a Mathcad string, the second is a Mathcad string 
    {STRING, STRING} 
    };

	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO FluidProp = 
    {
    // Name by which mathcad will recognize the function
    "FluidProp",        
    
    // description of "FluidProp" parameters to be used by the Insert Function dialog box
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name",   
    
    // description of the function for the Insert Function dialog box       
    "returns a fluid-specific parameter, where the parameter is dependent on the fluid state",    
    
    // pointer to the executible code. i.e. code that should be executed
    // when a user types in "FluidProp(PropName,Input1Name,Input1Value,Input2Name,Input2Value,FluidName)="
    (LPCFUNCTION)CoolPropMathcad,  
        
    // FluidProp(Fluid,PropName) returns a Mathcad complex scalar
    COMPLEX_SCALAR,
        
    // FluidProp(Fluid,PropName) takes on two arguments
    6,   
    
    // argument types 
    {STRING, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR, STRING} 
    };

	// fill out a FUNCTIONINFO structure with the information needed for registering
	// the function with Mathcad
    FUNCTIONINFO HAProp = 
    {
    // Name by which mathcad will recognize the function
    "HAProp",        
    
    // description of "HAProp" parameters to be used by the Insert Function dialog box
    "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3",   
    
    // description of the function for the Insert Function dialog box       
    "returns a parameter of humid air, where the parameter is dependent on the fluid state and humidity",    
    
    // pointer to the executible code. i.e. code that should be executed
    // when a user types in "HumidAirProp(PropName, )="
    (LPCFUNCTION)HumidAirPropMathcad,  
        
    // FluidProp(Fluid,PropName) returns a Mathcad complex scalar
    COMPLEX_SCALAR,
        
    // FluidProp(Fluid,PropName) takes on two arguments
    7,   
    
    // argument types 
    {STRING, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR, STRING, COMPLEX_SCALAR} 
    };


    // this code executes the user function FluidLIstMathcad, which is a wrapper for
	// the CoolProp.FluidsList() function, used to get a list of all the fluids in
	// the CoolProp database
    LRESULT  GetCoolPropsParamMathcad(
        LPMCSTRING ParamValue, // output (value of parameter)
		LPCMCSTRING ParamName ) // name of parameter (string) to retrieve
    {  
        // invoke the CoolProp get_global_param_string() function, save result to a new string s
	
		std::string s = get_global_param_string(std::string(ParamName->str));
		char * c = MathcadAllocate(s.size()+1); // creat a c-string (pointer) c with the same size as s
		// copy s into c, this process avoids the const-cast type which would result from instead
		// converting the string using s.c_str()
		std::copy(s.begin(), s.end(), c); 
		c[s.size()] = '\0';
		// assign the fluids list string to the function's output parameter
		ParamValue->str = c;

		// normal return
        return 0;
    }

    // this code executes the user function CoolPropMathcad1, which is a wrapper for
	// the CoolProp.CoolProp.Props() function, "Call Type #1," used to extract a
	// fluid-specific parameter that is not dependent on the state
    LRESULT  CoolPropMathcad1(
		COMPLEXSCALAR * const Prop1, // pointer to the result
        const MCSTRING * const Fluid, // string with a valid CoolProp fluid name
		const MCSTRING * const PropName ) // a fluid property
    {  
        // pass the arguments to the CoolProp.Props1() function
		Prop1->real = Props1(Fluid->str, PropName->str);

		// normal return
        return 0;
    }

	// this code executes the user function CoolPropMathcad, which is a wrapper for
	// the CoolProp.CoolProp.Props() function, "Call Type #2," used to extract a
	// fluid-specific parameter that is dependent on the state
    LRESULT  CoolPropMathcad(
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
		Prop->real = Props(OutputName->str, *InputName1->str, InputProp1->real, *InputName2->str, InputProp2->real, FluidName->str);
		
		// normal return
        return 0;
    }



	// this code executes the user function HumidAirPropMathcad, which is a wrapper for
	// the CoolProp.HumidAirlProp.HAProps() function, used to extract a
	// fluid-specific parameter that is dependent on the state
    LRESULT  HumidAirPropMathcad(
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
		Prop->real = HAProps(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, InputName3->str, InputProp3->real);
		
		// normal return
        return 0;
    }

    // DLL entry point code.  the _CRT_INIT function is needed if you are using Microsoft's 32 bit compiler
	// Rather than initializting in the DLLMain function below, a linker option was added to the Win32 platform
	// configuration to initialize the C run-time via the following entry point option: -entry:_DllMainCRTStartup@12
	BOOL APIENTRY DllMain(HINSTANCE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved )
	{
		switch (ul_reason_for_call)
		{
		case DLL_PROCESS_ATTACH:
			{// DLL is attaching to the address space of the current process.
      
			/*if (!_CRT_INIT(hModule, ul_reason_for_call, lpReserved)) {
				return FALSE;
			}*/

			if (!CreateUserErrorMessageTable( hModule, NUMBER_OF_ERRORS, myErrorMessageTable ) )
				break;
        
			if ( CreateUserFunction( hModule, &FluidProp1 ) == NULL )
					break;
	                        
	//				CreateUserFunction( hModule, &FluidProp1 );
					CreateUserFunction( hModule, &FluidProp );
					CreateUserFunction( hModule, &HAProp );
					CreateUserFunction( hModule, &FluidPropsParams );
			break;
			}
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
		}
		return TRUE;
	}


#undef INTERRUPTED    
#undef INSUFFICIENT_MEMORY
#undef MUST_BE_REAL   
#undef NUMBER_OF_ERRORS     

    


