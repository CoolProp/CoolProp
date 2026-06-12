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
#include <vector>
#include "CoolProp/CoolProp.h"
#include "CoolProp/CoolPropLib.h"
#include "CoolProp/detail/tools.h"

static bool EES_DEBUG = false;

// Structure for handling ees calling syntax
struct EesParamRec
{
    double value;
    struct EesParamRec* next;
};

using namespace CoolProp;

// EES always passes `fluid` as a fixed 256-char buffer (see the header comment
// at the top of this file).  Route every write to it through this helper so a
// long fluid/input string or error/warning message is truncated to fit instead
// of overrunning the buffer; the result is always NUL-terminated.
static void set_fluid(char* fluid, const std::string& message) {
    const std::size_t n = std::min(message.size(), static_cast<std::size_t>(255));
    message.copy(fluid, n);
    fluid[n] = '\0';
}

// Append a line to the EES debug log, skipping silently if the file can't be
// opened.  Debug logging must never crash or abort the EES call, so the fopen
// result is always checked before use.
static void log_debug(const std::string& line) {
    FILE* fp = fopen("log.txt", "a+");
    if (fp != nullptr) {
        fputs(line.c_str(), fp);
        fclose(fp);
    }
}

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names
extern "C"
{
    __declspec(dllexport) double COOLPROP_EES(char fluid[256], int mode, struct EesParamRec* input_rec) {
        double In1 = _HUGE, In2 = _HUGE, out = _HUGE;  // Two inputs, one output
        int NInputs = 0;                               // Ninputs is the number of inputs
        std::string fluid_string = fluid;

        std::vector<double> z;

        std::string Outstr, In1str, In2str, Fluidstr, Units;
        std::vector<std::string> fluid_split;

        if (mode == -1) {
            set_fluid(fluid, "T = PropsSI('T','P',101325,'Q',0,'Water')");
            return 0;
        }

        // Split the string that is passed in at the '~' delimiter that was used to join it
        fluid_split = strsplit(fluid_string, '~');
        if (fluid_split.size() != 5) {
            const std::string msg = format("fluid[%s] length[%d] not 5 elements long", fluid_string.c_str(), static_cast<int>(fluid_split.size()));
            set_fluid(fluid, msg);
            if (EES_DEBUG) {
                log_debug(format("%s %s %g %s %g %s\n%s\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str(), msg.c_str()));
            }
            return 0;
        } else {
            Fluidstr = upper(fluid_split[0]);
            Outstr = upper(fluid_split[1]);
            In1str = upper(fluid_split[2]);
            In2str = upper(fluid_split[3]);
            Units = upper(fluid_split[4]);
        }

        const std::size_t debug_pos = Fluidstr.find("$DEBUG");
        if (debug_pos != std::string::npos) {
            EES_DEBUG = true;
            Fluidstr.resize(debug_pos);
        } else {
            EES_DEBUG = false;
        }

        // Check the number of inputs
        EesParamRec* aninput_rec = input_rec;
        while (aninput_rec != nullptr) {
            if (NInputs >= 2) {
                z.push_back(aninput_rec->value);
            }
            aninput_rec = aninput_rec->next;
            NInputs++;
        };

        if (NInputs < 2) {
            set_fluid(fluid, format("Number of inputs [%d] < 2", NInputs));
            return 0;
        }

        // TODO: check that the number of components agrees with the length of array

        // Get the inputs from the pointer structure sent by EES:
        In1 = input_rec->value;
        input_rec = input_rec->next;
        In2 = input_rec->value;

        //This block can be used to debug the code by writing output or intermediate values to a text file

        if (EES_DEBUG) {
            log_debug(
              format("Inputs: %s %s %g %s %g %s %s\n", Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str(), Units.c_str()));
        }

        if (EES_DEBUG) {
            // This redirects standard output to log_stdout.txt; only crank up
            // the debug output if the redirect succeeded, so we don't spew to a
            // broken stdout when the log file can't be opened.
            if (freopen("log_stdout.txt", "w", stdout) != nullptr) {
                ::set_debug_level(100000);  // Maximum debugging
            }
        }

        try {
            if (Units == "SI") {
                if (!z.empty()) {
                    std::string backend, fluid_only;
                    extract_backend(Fluidstr, backend, fluid_only);
                    // Vectorize the inputs
                    std::vector<std::string> fluids = strsplit(fluid_only, '&');
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
                if (In1str.size() != 1) {
                    set_fluid(fluid, format("Input #1 [%s] can only be 1 character long for coolprop()", In1str.c_str()));
                    return 0;
                }
                if (In2str.size() != 1) {
                    set_fluid(fluid, format("Input #2 [%s] can only be 1 character long for coolprop()", In2str.c_str()));
                    return 0;
                }
                // Mole fractions are not given
                out = Props(Outstr.c_str(), In1str[0], In1, In2str[0], In2, Fluidstr.c_str());
            }
        } catch (...) {
            std::string error_message = format("Uncaught error: \"%s\",\"%s\",%g,\"%s\",%g,\"%s\"\n", Outstr.c_str(), In1str.c_str(), In1,
                                               In2str.c_str(), In2, Fluidstr.c_str());
            // There was an error
            if (EES_DEBUG) {
                log_debug(format("Error: %s \n", error_message.c_str()));
            }
            set_fluid(fluid, error_message);

            return 0.0;
        }

        if (!ValidNumber(out)) {
            std::string error_message = CoolProp::get_global_param_string("errstring");
            // There was an error
            if (EES_DEBUG) {
                log_debug(format("Error: %s \n", error_message.c_str()));
            }
            set_fluid(fluid, error_message);
            return 0.0;
        } else {
            // Check if there was a warning
            std::string warn_string = CoolProp::get_global_param_string("warnstring");
            if (!warn_string.empty()) {
                if (EES_DEBUG) {
                    log_debug(format("Warning: %s \n", warn_string.c_str()));
                }
                // There was a warning, write it back
                set_fluid(fluid, warn_string);
            }
            if (EES_DEBUG) {
                log_debug(format("Output: %g\n", out));
            }
            return out;
        }
    }
};
