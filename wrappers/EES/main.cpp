//============================================================================================//
//                                                                                            //
//                                  EES - CoolProp interface                                  //
//                                  -------------------------                                 //
//                                                                                            //
//  This dll is an interface between EES and CoolProp.                                        //
//  In EES, external functions need to be implemented in dynamic libraries.  The first       //
//  argument sent by EES is a 256 characters char variable. The second argument is pointer    //
//  structure containing "double" values.  The third argument is a linked list for the         //
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
#include <algorithm>
#include <string>
#include <stdio.h>
#include <string.h>
#include <vector>
#include "CoolProp.h"
#include "CoolPropLib.h"
#include "CoolPropTools.h"

static bool EES_DEBUG = false;

// Structure for handling ees calling syntax
struct EesParamRec
{
    double value;
    struct EesParamRec* next;
};

using namespace CoolProp;

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names
extern "C"
{
    __declspec(dllexport) double COOLPROP_EES(char fluid[256], int mode, struct EesParamRec* input_rec) {
        double In1 = _HUGE, In2 = _HUGE, out;  // Two inputs, one output
        int NInputs;                           // Ninputs is the number of inputs
        char NInputs_string[3], err_str[1000];
        std::string fluid_string = fluid;

        std::vector<double> z;

        std::string ErrorMsg, Outstr, In1str, In2str, Fluidstr, Units;
        std::vector<std::string> fluid_split;

        if (mode == -1) {
            strcpy(fluid, "T = PropsSI('T','P',101325,'Q',0,'Water')");
            return 0;
        }

        // Split the string that is passed in at the '~' delimiter that was used to join it
        fluid_split = strsplit(fluid_string, '~');
        if (fluid_split.size() != 5) {
            sprintf(err_str, "fluid[%s] length[%d] not 5 elements long", fluid_string.c_str(), fluid_split.size());
            strcpy(fluid, err_str);
            if (EES_DEBUG) {
                FILE* fp;
                fp = fopen("log.txt", "a+");
                fprintf(fp, "%s %s %g %s %g %s\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str());
                fprintf(fp, "%s\n", err_str);
                fclose(fp);
            }
            return 0;
        } else {
            Fluidstr = upper(fluid_split[0]);
            Outstr = upper(fluid_split[1]);
            In1str = upper(fluid_split[2]);
            In2str = upper(fluid_split[3]);
            Units = upper(fluid_split[4]);
        }

        if (Fluidstr.find("$DEBUG") != std::string::npos) {
            EES_DEBUG = true;
            Fluidstr = Fluidstr.substr(0, Fluidstr.find("$DEBUG"));
        } else {
            EES_DEBUG = false;
        }

        // Check the number of inputs
        NInputs = 0;
        EesParamRec* aninput_rec = input_rec;
        while (aninput_rec != 0) {
            if (NInputs >= 2) {
                z.push_back(aninput_rec->value);
            }
            aninput_rec = aninput_rec->next;
            NInputs++;
        };

        if (NInputs < 2) {
            sprintf(NInputs_string, "Number of inputs [%d] < 2", NInputs);
            strcpy(fluid, NInputs_string);
            return 0;
        }

        // TODO: check that the number of components agrees with the length of array

        // Get the inputs from the pointer structure sent by EES:
        In1 = input_rec->value;
        input_rec = input_rec->next;
        In2 = input_rec->value;

        //This block can be used to debug the code by writing output or intermediate values to a text file

        if (EES_DEBUG) {
            FILE* fp;
            fp = fopen("log.txt", "a+");
            fprintf(fp, "Inputs: %s %s %g %s %g %s %s\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str(), Units.c_str());
            fclose(fp);
        }

        if (EES_DEBUG) {
            // This redirects standard output to log_stdout.txt
            freopen("log_stdout.txt", "w", stdout);
            ::set_debug_level(100000);  // Maximum debugging
        }

        try {
            if (!Units.compare("SI")) {
                if (z.size() > 0) {
                    std::string backend, fluid;
                    extract_backend(Fluidstr, backend, fluid);
                    // Vectorize the inputs
                    std::vector<std::string> fluids = strsplit(fluid, '&');
                    std::vector<std::string> outputs(1, Outstr);
                    std::vector<double> val1(1, In1);
                    std::vector<double> val2(1, In2);
                    // Mole fractions are given, we use the advanced PropsSImulti function
                    std::vector<std::vector<double>> IO = PropsSImulti(outputs, In1str, val1, In2str, val2, backend, fluids, z);
                    if (IO.size() != 1 || IO[0].size() != 1) {
                        out = _HUGE;
                    } else {
                        out = IO[0][0];
                    }
                } else {
                    // Mole fractions are not given
                    out = PropsSI(Outstr, In1str, In1, In2str, In2, Fluidstr);
                }
            } else {
                if (In1str.size() != 0) {
                    strcpy(fluid, format("Input #1 [%s] can only be 1 character long for coolprop()", In1str.c_str()).c_str());
                }
                if (In2str.size() != 0) {
                    strcpy(fluid, format("Input #2 [%s] can only be 1 character long for coolprop()", In2str.c_str()).c_str());
                }
                // Mole fractions are not given
                out = Props(Outstr.c_str(), In1str[0], In1, In2str[0], In2, Fluidstr.c_str());
            }
        } catch (...) {
            std::string err_str = format("Uncaught error: \"%s\",\"%s\",%g,\"%s\",%g,\"%s\"\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(),
                                         In2, Fluidstr.c_str());
            // There was an error
            if (EES_DEBUG) {
                FILE* fp;
                fp = fopen("log.txt", "a+");
                fprintf(fp, "Error: %s \n", err_str.c_str());
                fclose(fp);
            }
            strcpy(fluid, err_str.c_str());

            return 0.0;
        }

        if (!ValidNumber(out)) {
            std::string err_str = CoolProp::get_global_param_string("errstring");
            // There was an error
            if (EES_DEBUG) {
                FILE* fp;
                fp = fopen("log.txt", "a+");
                fprintf(fp, "Error: %s \n", err_str.c_str());
                fclose(fp);
            }
            strcpy(fluid, err_str.c_str());
            return 0.0;
        } else {
            // Check if there was a warning
            std::string warn_string = CoolProp::get_global_param_string("warnstring");
            if (!warn_string.empty()) {
                if (EES_DEBUG) {
                    FILE* fp;
                    fp = fopen("log.txt", "a+");
                    fprintf(fp, "Warning: %s \n", warn_string.c_str());
                    fclose(fp);
                }
                // There was a warning, write it back
                strcpy(fluid, warn_string.c_str());
            }
            if (EES_DEBUG) {
                FILE* fp;
                fp = fopen("log.txt", "a+");
                fprintf(fp, "Output: %g\n", out);
                fclose(fp);
            }
            return out;
        }
    }
};
