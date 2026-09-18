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
//  - The string variable carries five fields joined by ~, in the order fluid,                //
//    output key, first input key, second input key and unit system, and the                  //
//    fluid field holds the mixture composition when there is one                             //
//    (e.g. "R134a~D~T~P~SI" or "REFPROP-R134a~O~T~P~SI" or                                   //
//    "REFPROP-MIX:R32[0.697615]&R125[0.302385]~V~P~H~SI" (R410A))                            //
//  - mode, which EES passes BY REFERENCE (see the F-Chart help, "External                    //
//    Functions" and the Visual C++ skeleton).  EES asks for a description of                 //
//    the call with mode = -1, for the units of the inputs with mode = -2 and                 //
//    for the units of the output with mode = -3.  Any other value means a                    //
//    normal call.  On the way back the mode says what happened: 0 and an                     //
//    empty string for a normal result, a positive value and a message for an                 //
//    error, which stops the calculation, and a negative value and a message                  //
//    for a warning, which does not.                                                          //
//  - The last value is a linked list of the input values                                     //
//																							  //
//  The file needs to be built in coolprop_ees.dlf, which is the standard extension           //
//  for EES external functions.  If CoolProp has been built to the static library             //
//  CoolPropStaticLibrary.lib, you can build (with visual studio) CoolProp_EES.dlf with:      //
//                                                                                            //
//     link /DEBUG /DLL main.obj CoolPropStaticLibrary.lib /OUT:COOLPROP_EES.dlf              //
//																							  //
//  Base SI units are the only ones this library accepts (K, Pa, J, mass).  The               //
//  unit system is named in the last field of the string and CoolProp.LIB checks              //
//  that the EES unit system matches before it calls in here.                                 //
//																							  //
//  Ian Bell                                                                                  //
//  Thermodynamics Laboratory                                                                 //
//  University of Liege                                                                       //
//                                                                                            //
//  January 2013                                                                              //
//============================================================================================//

#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
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

// Append a line to the EES debug log, skipping silently if the file cannot be
// opened.  Debug logging must never crash or abort the EES call, so the stream
// state is always checked before use.  The log holds the fluid strings of the
// calling model and lands in the folder EES runs from, where it inherits that
// folder's access rights.  Windows ignores std::filesystem::permissions apart
// from the read-only flag, so there is nothing portable to tighten here.
static void log_debug(const std::string& line) {
    std::ofstream log_file("log.txt", std::ios::app);
    if (!log_file) {
        return;
    }
    log_file << line;
}

// Report an error back to EES: the message goes into the string and the mode
// has to be positive, which makes EES stop the calculation and show it.  A
// warning uses a negative mode instead, see the header comment.  An empty
// message is replaced: stopping the calculation without saying why would be
// worse than the wrong-but-silent zero this used to return, and CoolProp does
// hand out an empty error string now and then because reading it clears it.
static void set_error(char* fluid, int& mode, const std::string& message) {
    set_fluid(fluid, message.empty() ? std::string("CoolProp failed without reporting a reason") : message);
    mode = 1;
}

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names
extern "C"
{
    // The mode argument is a reference because that is how EES passes it, see
    // the file header.  Taking it by value (as this wrapper did until 2026)
    // reads the low bits of the pointer instead of the mode, so the requests
    // below were never served.
    __declspec(dllexport) double COOLPROP_EES(char fluid[256], int& mode, struct EesParamRec* input_rec) {
        double In1 = _HUGE, In2 = _HUGE, out = _HUGE;  // Two inputs, one output
        int NInputs = 0;                               // Ninputs is the number of inputs

        std::vector<double> z;

        std::string Outstr, In1str, In2str, Fluidstr, Units;
        std::vector<std::string> fluid_split;

        // The three requests below answer in the string and leave the mode as
        // EES set it, which is what the F-Chart example does.  The 0 / positive
        // / negative convention described in the header applies to a normal
        // call, where the mode says what came of the calculation.
        if (mode == -1) {
            // EES asks for an example of the call format
            set_fluid(fluid, "T = PropsSI('T','P',101325,'Q',0,'Water')");
            return 0;
        }

        if (mode == -2 || mode == -3) {
            // EES asks for the units of the inputs (-2) or of the output (-3).
            // Both depend on the property keys that are encoded in the fluid
            // string of the actual call, so there is no fixed answer here.  An
            // empty string tells EES that no units are declared.
            set_fluid(fluid, "");
            return 0;
        }

        // Only a normal call fills the buffer, so read it here rather than
        // above.  The length is bounded because EES does not promise a
        // terminating NUL for the requests handled above.
        const std::string fluid_string(fluid, ::strnlen(fluid, 255));

        // Split the string that is passed in at the '~' delimiter that was used to join it
        fluid_split = strsplit(fluid_string, '~');
        if (fluid_split.size() != 5) {
            const std::string msg = format("fluid[%s] length[%d] not 5 elements long", fluid_string.c_str(), static_cast<int>(fluid_split.size()));
            set_error(fluid, mode, msg);
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
            set_error(fluid, mode, format("Number of inputs [%d] < 2", NInputs));
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
                // The deprecated coolprop() and coolpropsi() functions were the
                // only ones that put anything but SI here, and they were removed
                // from CoolProp.LIB in September 2026.  Reaching this point means
                // the library file and this DLL come from different releases, so
                // say that rather than guess which units the numbers are in.
                // The remedy comes first: a long unit string pushes the end of
                // the message past the 255 characters the buffer holds.
                set_error(fluid, mode,
                          format("Use PropsSI, and install CoolProp.LIB and COOLPROP_EES from the same CoolProp release. "
                                 "The deprecated coolprop() and coolpropsi() functions were removed, so unit system [%s] is no "
                                 "longer supported.",
                                 Units.c_str()));
                return 0;
            }
        } catch (...) {
            std::string error_message = format("Uncaught error: \"%s\",\"%s\",%g,\"%s\",%g,\"%s\"\n", Outstr.c_str(), In1str.c_str(), In1,
                                               In2str.c_str(), In2, Fluidstr.c_str());
            // There was an error
            if (EES_DEBUG) {
                log_debug(format("Error: %s \n", error_message.c_str()));
            }
            set_error(fluid, mode, error_message);

            return 0.0;
        }

        if (!ValidNumber(out)) {
            std::string error_message = CoolProp::get_global_param_string("errstring");
            // There was an error
            if (EES_DEBUG) {
                log_debug(format("Error: %s \n", error_message.c_str()));
            }
            set_error(fluid, mode, error_message);
            return 0.0;
        } else {
            // A normal call returns the null string and sets the mode to 0.
            set_fluid(fluid, "");
            mode = 0;
            // Check if there was a warning
            std::string warn_string = CoolProp::get_global_param_string("warnstring");
            if (!warn_string.empty()) {
                if (EES_DEBUG) {
                    log_debug(format("Warning: %s \n", warn_string.c_str()));
                }
                // There was a warning, write it back.  A negative mode makes
                // EES show the message without stopping the calculation.
                set_fluid(fluid, warn_string);
                mode = -1;
            }
            if (EES_DEBUG) {
                log_debug(format("Output: %g\n", out));
            }
            return out;
        }
    }
};
