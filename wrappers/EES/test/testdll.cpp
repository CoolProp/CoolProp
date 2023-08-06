//============================================================================================//
//                                                                                            //
//                                  EES - CoolProp interface                                  //
//                                  -------------------------                                 //
//                                                                                            //
//  This dll is an interface between EES and CoolProp.                                        //
//  In EES, external functions need to be implemented in dynamic librairies.  The first       //
//  argument sent by EES is a 256 characters char variable. The second argument is pointer    //
//  structure conaining "double" values.  The third argument is a linked list for the         //
//  input data                                                                                //
//                                                                                            //
//  The arguments are defined as follows :                                                    //
//  - The string variable contains the the definition of the fluids and of their              //
//    concentrations with the input strings concatenated to the fluid name joined by |        //
//    (e.g. "R134a|T|P|D" or "REFPROP-R134a|O|T|P" or                                         //
//    "REFPROP-MIX:R32[0.697615]&R125[0.302385]|V|P|H" (R410A))                               //
//  - mode, which is -1 if to return a default form of the call as string, normal mode        //
//    otherwise                                                                               //
//  - The last value is a linked list of the input values                                     //
//																							  //
//  The file needs to be built in coolprop_ees.dlf, which is the standard extension           //
//  for EES external functions.  If CoolProp has been built to the static library             //
//  CoolPropStaticLibrary.lib, you can build (with visual studio) CoolProp_EES.dlf with:      //
//                                                                                            //
//     link /DEBUG /DLL main.obj CoolPropStaticLibrary.lib /OUT:COOLPROP_EES.dlf              //
//																							  //
//  Only one unit system is used (modified SI - see help). Future versions might              //
//  include a detection of EES current unit system and its definition in the dll              //
//																							  //
//  Ian Bell                                                                                  //
//  Thermodynamics Laboratory                                                                 //
//  University of Liège                                                                       //
//                                                                                            //
//  January 2013                                                                              //
//============================================================================================//

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <vector>

#include "CoolPropTools.h"
#include <algorithm>
#include <string>

#include <windows.h>

HINSTANCE CoolPropdllInstance;

typedef double(__stdcall* fp_PropsSdllTYPE)(char*, char*, double, char*, double, char*);
fp_PropsSdllTYPE PropsSdll;

typedef double(__stdcall* fp_set_debug_leveldllTYPE)(int);
fp_set_debug_leveldllTYPE set_debug_leveldll;

typedef double(__stdcall* fp_get_global_param_stringdllTYPE)(char*, char*);
fp_get_global_param_stringdllTYPE get_global_param_stringdll;

static const bool EES_DEBUG = false;

// Structure for handling ees calling syntax

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names
int main() {
    double In1 = _HUGE, In2 = _HUGE, out;  // Two inputs, one output
    int NInputs;                           // Ninputs is the number of inputs
    char NInputs_string[3], err_str[1000];

    std::string ErrorMsg, Outstr, In1str, In2str, Fluidstr;
    std::vector<std::string> fluid_split;

    Outstr = "S";
    In1str = "P";
    In1 = 2449.047245069126;
    In2str = "H";
    In2 = 306.720082386865;
    Fluidstr = "R134A";

    //This block can be used to debug the code by writing output or intermediate values to a text file

    if (EES_DEBUG) {
        FILE* fp;
        fp = fopen("log.txt", "a+");
        fprintf(fp, "%s %s %g %s %g %s\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str());
        fclose(fp);
    }

    // 32-bit code here.
    TCHAR coolpropdllstring[100] = TEXT("CoolProp.dll");
    CoolPropdllInstance = LoadLibrary(coolpropdllstring);

    if (CoolPropdllInstance == NULL) {
        printf("CoolProp.dll could not be loaded");
        return 0;
    }
    PropsSdll = (fp_PropsSdllTYPE)GetProcAddress(CoolPropdllInstance, "_PropsS@32");
    set_debug_leveldll = (fp_set_debug_leveldllTYPE)GetProcAddress(CoolPropdllInstance, "_set_debug_level@4");
    get_global_param_stringdll = (fp_get_global_param_stringdllTYPE)GetProcAddress(CoolPropdllInstance, "_get_global_param_string@8");

    out = PropsSdll((char*)Outstr.c_str(), (char*)In1str.c_str(), In1, (char*)In2str.c_str(), In2, (char*)Fluidstr.c_str());

    printf(format("output: %g\n", out).c_str());

    if (fabs(out) > 1e90) {
        char err_chars[10000];
        get_global_param_stringdll("errstring", err_chars);
        std::string err_str = err_chars;
        // There was an error
        if (EES_DEBUG) {
            FILE* fp;
            fp = fopen("log.txt", "a+");
            fprintf(fp, "Error: %s \n", err_str.c_str());
            fclose(fp);
        }
        printf(err_str.c_str());
        return 0;
    } else {
        // Check if there was a warning
        char warn_chars[10000];
        get_global_param_stringdll("warnstring", warn_chars);
        std::string warn_string = warn_chars;
        if (!warn_string.empty()) {
            if (EES_DEBUG) {
                FILE* fp;
                fp = fopen("log.txt", "a+");
                fprintf(fp, "Warning: %s \n", warn_string.c_str());
                fclose(fp);
            }
            // There was a warning, write it back
            printf(warn_string.c_str());
        }
        return out;
    }
    return 0;
}
