#if defined(_MSC_VER)
#    ifndef _CRTDBG_MAP_ALLOC
#        define _CRTDBG_MAP_ALLOC
#    endif
#    ifndef _CRT_SECURE_NO_WARNINGS
#        define _CRT_SECURE_NO_WARNINGS
#    endif
#    include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "AbstractState.h"

#if defined(__ISWINDOWS__)
#    include <windows.h>
#    ifdef min
#        undef min
#    endif
#    ifdef max
#        undef max
#    endif
#else
#    ifndef DBL_EPSILON
#        include <limits>
#        define DBL_EPSILON std::numeric_limits<double>::epsilon()
#    endif
#endif

#include <memory>

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <exception>
#include <stdio.h>
#include <string>
#include <locale>
#include "CoolPropTools.h"
#include "Solvers.h"
#include "MatrixMath.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "Backends/Incompressible/IncompressibleLibrary.h"
#include "Backends/Incompressible/IncompressibleBackend.h"
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Helmholtz/MixtureParameters.h"
#include "DataStructures.h"
#include "Backends/REFPROP/REFPROPMixtureBackend.h"
#include "Backends/Cubics/CubicsLibrary.h"
#include "Backends/PCSAFT/PCSAFTLibrary.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>
#endif

namespace CoolProp {

static int debug_level = 0;
static std::string error_string;
static std::string warning_string;

void set_debug_level(int level) {
    debug_level = level;
}
int get_debug_level(void) {
    return debug_level;
}

//// This is very hacky, but pull the git revision from the file
#include "gitrevision.h"  // Contents are like "std::string gitrevision = "aa121435436ggregrea4t43t433";"
#include "cpversion.h"    // Contents are like "char version [] = "2.5";"

void set_warning_string(const std::string& warning) {
    warning_string = warning;
}
void set_error_string(const std::string& error) {
    error_string = error;
}

// Return true if the string has "BACKEND::*" format where * signifies a wildcard
bool has_backend_in_string(const std::string& fluid_string, std::size_t& i) {
    i = fluid_string.find("::");
    return i != std::string::npos;
}

void extract_backend(std::string fluid_string, std::string& backend, std::string& fluid) {
    std::size_t i;
    // For backwards compatibility reasons, if "REFPROP-" or "REFPROP-MIX:" start
    // the fluid_string, replace them with "REFPROP::"
    if (fluid_string.find("REFPROP-MIX:") == 0) {
        fluid_string.replace(0, 12, "REFPROP::");
    }
    if (fluid_string.find("REFPROP-") == 0) {
        fluid_string.replace(0, 8, "REFPROP::");
    }
    if (has_backend_in_string(fluid_string, i)) {
        // Part without the ::
        backend = fluid_string.substr(0, i);
        // Fluid name after the ::
        fluid = fluid_string.substr(i + 2);
    } else {
        backend = "?";
        fluid = fluid_string;
    }
    if (get_debug_level() > 10)
        std::cout << format("%s:%d: backend extracted. backend: %s. fluid: %s\n", __FILE__, __LINE__, backend.c_str(), fluid.c_str());
}

bool has_fractions_in_string(const std::string& fluid_string) {
    // If can find both "[" and "]", it must have mole fractions encoded as string
    return (fluid_string.find("[") != std::string::npos && fluid_string.find("]") != std::string::npos);
}
bool has_solution_concentration(const std::string& fluid_string) {
    // If can find "-", expect mass fractions encoded as string
    return (fluid_string.find('-') != std::string::npos && fluid_string.find('%') != std::string::npos);
}

struct delim : std::numpunct<char>
{
    char m_c;
    delim(char c) : m_c(c){};
    char do_decimal_point() const {
        return m_c;
    }
};

std::string extract_fractions(const std::string& fluid_string, std::vector<double>& fractions) {

    if (has_fractions_in_string(fluid_string)) {
        fractions.clear();
        std::vector<std::string> names;

        // Break up into pairs - like "Ethane[0.5]&Methane[0.5]" -> ("Ethane[0.5]","Methane[0.5]")
        std::vector<std::string> pairs = strsplit(fluid_string, '&');

        for (std::size_t i = 0; i < pairs.size(); ++i) {
            const std::string& fluid = pairs[i];

            // Must end with ']'
            if (fluid[fluid.size() - 1] != ']') throw ValueError(format("Fluid entry [%s] must end with ']' character", pairs[i].c_str()));

            // Split at '[', but first remove the ']' from the end by taking a substring
            std::vector<std::string> name_fraction = strsplit(fluid.substr(0, fluid.size() - 1), '[');

            if (name_fraction.size() != 2) {
                throw ValueError(format("Could not break [%s] into name/fraction", fluid.substr(0, fluid.size() - 1).c_str()));
            }

            // Convert fraction to a double
            const std::string &name = name_fraction[0], &fraction = name_fraction[1];
            // The default locale for conversion from string to double is en_US with . as the deliminter
            // Good: 0.1234 Bad: 0,1234
            // But you can change the punctuation character for fraction parsing
            // with the configuration variable FLOAT_PUNCTUATION to change the locale to something more convenient for you (e.g., a ',')
            // See also http://en.cppreference.com/w/cpp/locale/numpunct/decimal_point
            std::stringstream ssfraction(fraction);
            char c = get_config_string(FLOAT_PUNCTUATION)[0];
            ssfraction.imbue(std::locale(ssfraction.getloc(), new delim(c)));
            double f;
            ssfraction >> f;
            if (ssfraction.rdbuf()->in_avail() != 0) {
                throw ValueError(format("fraction [%s] was not converted fully", fraction.c_str()));
            }

            if (f > 1 || f < 0) {
                throw ValueError(format("fraction [%s] was not converted to a value between 0 and 1 inclusive", fraction.c_str()));
            }

            if ((f > 10 * DBL_EPSILON) ||  // Only push component if fraction is positive and non-zero
                (pairs.size() == 1))       // ..or if there is only one fluid (i.e. INCOMP backend )
            {
                // And add to vector
                fractions.push_back(f);

                // Add name
                names.push_back(name);
            }
        }

        if (get_debug_level() > 10)
            std::cout << format("%s:%d: Detected fractions of %s for %s.", __FILE__, __LINE__, vec_to_string(fractions).c_str(),
                                (strjoin(names, "&")).c_str());
        // Join fluids back together
        return strjoin(names, "&");
    } else if (has_solution_concentration(fluid_string)) {
        fractions.clear();
        double x;

        std::vector<std::string> fluid_parts = strsplit(fluid_string, '-');
        // Check it worked
        if (fluid_parts.size() != 2) {
            throw ValueError(
              format("Format of incompressible solution string [%s] is invalid, should be like \"EG-20%\" or \"EG-0.2\" ", fluid_string.c_str()));
        }

        // Convert the concentration into a string
        char* pEnd;
        x = strtod(fluid_parts[1].c_str(), &pEnd);

        // Check if per cent or fraction syntax is used
        if (!strcmp(pEnd, "%")) {
            x *= 0.01;
        }
        fractions.push_back(x);
        if (get_debug_level() > 10)
            std::cout << format("%s:%d: Detected incompressible concentration of %s for %s.", __FILE__, __LINE__, vec_to_string(fractions).c_str(),
                                fluid_parts[0].c_str());
        return fluid_parts[0];
    } else {
        return fluid_string;
    }
}

void _PropsSI_initialize(const std::string& backend, const std::vector<std::string>& fluid_names, const std::vector<double>& z,
                         shared_ptr<AbstractState>& State) {

    if (fluid_names.empty()) {
        throw ValueError("fluid_names cannot be empty");
    }

    std::vector<double> fractions(1, 1.0);            // Default to one component, unity fraction
    const std::vector<double>* fractions_ptr = NULL;  // Pointer to the array to be used;

    if (fluid_names.size() > 1) {
        // Set the pointer - we are going to use the supplied fractions; they must be provided
        fractions_ptr = &z;
        // Reset the state
        State.reset(AbstractState::factory(backend, fluid_names));
    } else if (fluid_names.size() == 1) {
        if (has_fractions_in_string(fluid_names[0]) || has_solution_concentration(fluid_names[0])) {
            // Extract fractions from the string
            std::string fluid_string = extract_fractions(fluid_names[0], fractions);
            // Set the pointer - we are going to use the extracted fractions
            fractions_ptr = &fractions;
            // Reset the state
            State.reset(AbstractState::factory(backend, fluid_string));
        } else {
            if (z.empty()) {
                // Set the pointer - we are going to use the default fractions
                fractions_ptr = &fractions;
            } else {
                // Set the pointer - we are going to use the provided fractions
                fractions_ptr = &z;
            }
            // Reset the state
            State.reset(AbstractState::factory(backend, fluid_names));
        }
    } else {  // The only path where fractions_ptr stays NULL
        throw ValueError("fractions_ptr is NULL");
    }
    if (!State->available_in_high_level()) {
        throw ValueError(
          "This AbstractState derived class cannot be used in the high-level interface; see www.coolprop.org/dev/coolprop/LowLevelAPI.html");
    }

    // Set the fraction for the state
    if (State->using_mole_fractions()) {
        // If a predefined mixture or a pure fluid, the fractions will already be set
        if (State->get_mole_fractions().empty()) {
            State->set_mole_fractions(*fractions_ptr);
        }
    } else if (State->using_mass_fractions()) {
        State->set_mass_fractions(*fractions_ptr);
    } else if (State->using_volu_fractions()) {
        State->set_volu_fractions(*fractions_ptr);
    } else {
        if (get_debug_level() > 50)
            std::cout << format("%s:%d: _PropsSI, could not set composition to %s, defaulting to mole fraction.\n", __FILE__, __LINE__,
                                vec_to_string(z).c_str())
                           .c_str();
    }
}

struct output_parameter
{
    enum OutputParametersType
    {
        OUTPUT_TYPE_UNSET = 0,
        OUTPUT_TYPE_TRIVIAL,
        OUTPUT_TYPE_NORMAL,
        OUTPUT_TYPE_FIRST_DERIVATIVE,
        OUTPUT_TYPE_FIRST_SATURATION_DERIVATIVE,
        OUTPUT_TYPE_SECOND_DERIVATIVE
    };
    CoolProp::parameters Of1, Wrt1, Constant1, Wrt2, Constant2;
    OutputParametersType type;
    /// Parse a '&' separated string into a data structure with one entry per output
    /// Covers both normal and derivative outputs
    static std::vector<output_parameter> get_output_parameters(const std::vector<std::string>& Outputs) {
        std::vector<output_parameter> outputs;
        for (std::vector<std::string>::const_iterator str = Outputs.begin(); str != Outputs.end(); ++str) {
            output_parameter out;
            CoolProp::parameters iOutput;
            if (is_valid_parameter(*str, iOutput)) {
                out.Of1 = iOutput;
                if (is_trivial_parameter(iOutput)) {
                    out.type = OUTPUT_TYPE_TRIVIAL;
                } else {
                    out.type = OUTPUT_TYPE_NORMAL;
                }
            } else if (is_valid_first_saturation_derivative(*str, out.Of1, out.Wrt1)) {
                out.type = OUTPUT_TYPE_FIRST_SATURATION_DERIVATIVE;
            } else if (is_valid_first_derivative(*str, out.Of1, out.Wrt1, out.Constant1)) {
                out.type = OUTPUT_TYPE_FIRST_DERIVATIVE;
            } else if (is_valid_second_derivative(*str, out.Of1, out.Wrt1, out.Constant1, out.Wrt2, out.Constant2)) {
                out.type = OUTPUT_TYPE_SECOND_DERIVATIVE;
            } else {
                throw ValueError(format("Output string is invalid [%s]", str->c_str()));
            }
            outputs.push_back(out);
        }
        return outputs;
    };
};

void _PropsSI_outputs(shared_ptr<AbstractState>& State, const std::vector<output_parameter>& output_parameters, CoolProp::input_pairs input_pair,
                      const std::vector<double>& in1, const std::vector<double>& in2, std::vector<std::vector<double>>& IO) {

    // Check the inputs
    if (in1.size() != in2.size()) {
        throw ValueError(format("lengths of in1 [%d] and in2 [%d] are not the same", in1.size(), in2.size()));
    }
    bool one_input_one_output = (in1.size() == 1 && in2.size() == 1 && output_parameters.size() == 1);
    // If all trivial outputs, never do a state update
    bool all_trivial_outputs = true;
    for (std::size_t j = 0; j < output_parameters.size(); ++j) {
        if (output_parameters[j].type != output_parameter::OUTPUT_TYPE_TRIVIAL) {
            all_trivial_outputs = false;
        }
    }
    parameters p1, p2;
    // If all outputs are also inputs, never do a state update
    bool all_outputs_in_inputs = true;
    if (input_pair != INPUT_PAIR_INVALID) {
        // Split the input pair into parameters
        split_input_pair(input_pair, p1, p2);
        // See if each parameter is in the output vector and is a normal type input
        for (std::size_t j = 0; j < output_parameters.size(); ++j) {
            if (output_parameters[j].type != output_parameter::OUTPUT_TYPE_NORMAL) {
                all_outputs_in_inputs = false;
                break;
            }
            if (!(output_parameters[j].Of1 == p1 || output_parameters[j].Of1 == p2)) {
                all_outputs_in_inputs = false;
                break;
            }
        }
    } else {
        if (!all_trivial_outputs) {
            throw ValueError(format("Input pair variable is invalid and output(s) are non-trivial; cannot do state update"));
        }
        all_outputs_in_inputs = false;
    }

    if (get_debug_level() > 100) {
        std::cout << format("%s (%d): input pair = %d ", __FILE__, __LINE__, input_pair) << std::endl;
        std::cout << format("%s (%d): in1 = %s ", __FILE__, __LINE__, vec_to_string(in1).c_str()) << std::endl;
        std::cout << format("%s (%d): in2 = %s ", __FILE__, __LINE__, vec_to_string(in2).c_str()) << std::endl;
    }

    // Get configuration variable for line tracing, see #1443
    const bool use_guesses = get_config_bool(USE_GUESSES_IN_PROPSSI);
    GuessesStructure guesses;

    // Resize the output matrix
    std::size_t N1 = std::max(static_cast<std::size_t>(1), in1.size());
    std::size_t N2 = std::max(static_cast<std::size_t>(1), output_parameters.size());
    IO.resize(N1, std::vector<double>(N2, _HUGE));

    // Throw an error if at the end, there were no successes
    bool success = false;
    bool success_inner = false;

    if (get_debug_level() > 100) {
        std::cout << format("%s (%d): Iterating over %d input value pairs.", __FILE__, __LINE__, IO.size()) << std::endl;
    }

    // Iterate over the state variable inputs
    for (std::size_t i = 0; i < IO.size(); ++i) {
        // Reset the success indicator for the current state point
        success_inner = false;
        try {
            if (input_pair != INPUT_PAIR_INVALID && !all_trivial_outputs && !all_outputs_in_inputs) {
                // Update the state since it is a valid set of inputs
                if (!use_guesses || i == 0) {
                    State->update(input_pair, in1[i], in2[i]);
                } else {
                    State->update_with_guesses(input_pair, in1[i], in2[i], guesses);
                    guesses.clear();
                }
            }
        } catch (...) {
            if (one_input_one_output) {
                IO.clear();
                throw;
            }  // Re-raise the exception since we want to bubble the error
            // All the outputs are filled with _HUGE; go to next input
            for (std::size_t j = 0; j < IO[i].size(); ++j) {
                IO[i][j] = _HUGE;
            }
            continue;
        }

        for (std::size_t j = 0; j < IO[i].size(); ++j) {
            // If all the outputs are inputs, there is no need for a state input
            if (all_outputs_in_inputs) {
                if (p1 == output_parameters[j].Of1) {
                    IO[i][j] = in1[i];
                    success_inner = true;
                    continue;
                } else if (p2 == output_parameters[j].Of1) {
                    IO[i][j] = in2[i];
                    success_inner = true;
                    continue;
                } else {
                    throw ValueError();
                }
            }
            try {
                const output_parameter& output = output_parameters[j];
                switch (output.type) {
                    case output_parameter::OUTPUT_TYPE_TRIVIAL:
                    case output_parameter::OUTPUT_TYPE_NORMAL:
                        IO[i][j] = State->keyed_output(output.Of1);
                        if (use_guesses) {
                            switch (output.Of1) {
                                case iDmolar:
                                    guesses.rhomolar = IO[i][j];
                                    break;
                                case iT:
                                    guesses.T = IO[i][j];
                                    break;
                                case iP:
                                    guesses.p = IO[i][j];
                                    break;
                                case iHmolar:
                                    guesses.hmolar = IO[i][j];
                                    break;
                                case iSmolar:
                                    guesses.smolar = IO[i][j];
                                    break;
                                default:
                                    throw ValueError("Don't understand this parameter");
                            }
                        }
                        break;
                    case output_parameter::OUTPUT_TYPE_FIRST_DERIVATIVE:
                        IO[i][j] = State->first_partial_deriv(output.Of1, output.Wrt1, output.Constant1);
                        break;
                    case output_parameter::OUTPUT_TYPE_FIRST_SATURATION_DERIVATIVE:
                        IO[i][j] = State->first_saturation_deriv(output.Of1, output.Wrt1);
                        break;
                    case output_parameter::OUTPUT_TYPE_SECOND_DERIVATIVE:
                        IO[i][j] = State->second_partial_deriv(output.Of1, output.Wrt1, output.Constant1, output.Wrt2, output.Constant2);
                        break;
                    default:
                        throw ValueError(format(""));
                        break;
                }
                // At least one has succeeded
                success_inner = true;
            } catch (...) {
                if (one_input_one_output) {
                    IO.clear();
                    throw;
                }  // Re-raise the exception since we want to bubble the error
                IO[i][j] = _HUGE;
            }
        }
        // We want to have at least rhomolar and T, but we do not raise errors here
        if (use_guesses && success_inner) {
            if (!ValidNumber(guesses.rhomolar)) {
                try {
                    guesses.rhomolar = State->rhomolar();
                } catch (...) {
                    guesses.rhomolar = _HUGE;
                }
            }
            if (!ValidNumber(guesses.T)) {
                try {
                    guesses.T = State->T();
                } catch (...) {
                    guesses.T = _HUGE;
                }
            }
        }
        // Save the success indicator, just a single valid output is enough
        success |= success_inner;
    }
    if (success == false) {
        IO.clear();
        throw ValueError(format("No outputs were able to be calculated"));
    }
}

bool StripPhase(std::string& Name, shared_ptr<AbstractState>& State)
// Parses an imposed phase out of the Input Name string using the "|" delimiter
{
    std::vector<std::string> strVec = strsplit(Name, '|');  // Split input key string in to vector containing input key [0] and phase string [1]
    if (strVec.size() > 1) {                                // If there is a phase string (contains "|" character)
        // Check for invalid backends for setting phase in PropsSI
        std::string strBackend = State->backend_name();
        if (strBackend == get_backend_string(INCOMP_BACKEND))
            throw ValueError("Cannot set phase on Incompressible Fluid; always liquid phase");  // incompressible fluids are always "liquid".
        if (strBackend == get_backend_string(IF97_BACKEND))
            throw ValueError("Can't set phase on IF97 Backend");  // IF97 has to calculate it's own phase region
        if (strBackend == get_backend_string(TTSE_BACKEND))
            throw ValueError("Can't set phase on TTSE Backend in PropsSI");  // Shouldn't be calling from High-Level anyway
        if (strBackend == get_backend_string(BICUBIC_BACKEND))
            throw ValueError("Can't set phase on BICUBIC Backend in PropsSI");  // Shouldn't be calling from High-Level anyway
        if (strBackend == get_backend_string(VTPR_BACKEND))
            throw ValueError("Can't set phase on VTPR Backend in PropsSI");  // VTPR has no phase functions to call

        phases imposed = iphase_not_imposed;  // Initialize imposed phase
        if (strVec.size() > 2)                // If there's more than on phase separator, throw error
        {
            throw ValueError(format("Invalid phase format: \"%s\"", Name));
        }
        // Handle prefixes of iphase_, phase_, or <none>
        std::basic_string<char>::iterator str_Iter;
        std::string strPhase = strVec[1];  //    Create a temp string so we can modify the prefix
        if (strPhase.find("iphase_") != strPhase.npos) {
            str_Iter = strPhase.erase(strPhase.begin());
        }  // Change "iphase_" to "phase_"
        if (strPhase.find("phase_") == strPhase.npos) {
            strPhase.insert(0, "phase_");
        }  // Prefix with "phase_" if missing
        // See if phase is a valid phase string, updating imposed while we're at it...
        if (!is_valid_phase(strPhase, imposed)) {
            throw ValueError(format("Phase string \"%s\" is not a valid phase", strVec[1]));  // throw error with original string if not valid
        }
        // Parsed phase string was valid
        Name = strVec[0];               //     Update input name to just the key string part
        State->specify_phase(imposed);  //     Update the specified phase on the backend State
        return true;                    //     Return true because a valid phase string was found
    }
    return false;  // Return false if there was no phase string on this key.
}

void _PropsSImulti(const std::vector<std::string>& Outputs, const std::string& Name1, const std::vector<double>& Prop1, const std::string& Name2,
                   const std::vector<double>& Prop2, const std::string& backend, const std::vector<std::string>& fluids,
                   const std::vector<double>& fractions, std::vector<std::vector<double>>& IO) {
    shared_ptr<AbstractState> State;
    CoolProp::parameters key1 = INVALID_PARAMETER, key2 = INVALID_PARAMETER;  // Initialize to invalid parameter values
    CoolProp::input_pairs input_pair = INPUT_PAIR_INVALID;                    // Initialize to invalid input pair
    std::vector<output_parameter> output_parameters;
    std::vector<double> v1, v2;

    try {
        // Initialize the State class
        _PropsSI_initialize(backend, fluids, fractions, State);
    } catch (std::exception& e) {
        // Initialization failed.  Stop.
        throw ValueError(format("Initialize failed for backend: \"%s\", fluid: \"%s\" fractions \"%s\"; error: %s", backend.c_str(),
                                strjoin(fluids, "&").c_str(), vec_to_string(fractions, "%0.10f").c_str(), e.what()));
    }

    //strip any imposed phase from input key strings here
    std::string N1 = Name1;                  // Make Non-constant copy of Name1 that we can modify
    std::string N2 = Name2;                  // Make Non-constant copy of Name2 that we can modify
    bool HasPhase1 = StripPhase(N1, State);  // strip phase string from first name if needed
    bool HasPhase2 = StripPhase(N2, State);  // strip phase string from second name if needed
    if (HasPhase1 && HasPhase2)              // if both Names have a phase string, don't allow it.
        throw ValueError("Phase can only be specified on one of the input key strings");

    try {
        // Get update pair
        if (is_valid_parameter(N1, key1) && is_valid_parameter(N2, key2)) input_pair = generate_update_pair(key1, Prop1, key2, Prop2, v1, v2);
    } catch (std::exception& e) {
        // Input parameter parsing failed.  Stop
        throw ValueError(format("Input pair parsing failed for Name1: \"%s\", Name2: \"%s\"; err: %s", Name1.c_str(), Name2.c_str(), e.what()));
    }

    try {
        output_parameters = output_parameter::get_output_parameters(Outputs);
    } catch (std::exception& e) {
        // Output parameter parsing failed.  Stop.
        throw ValueError(format("Output parameter parsing failed; error: %s", e.what()));
    }

    // Calculate the output(s).  In the case of a failure, all values will be filled with _HUGE
    _PropsSI_outputs(State, output_parameters, input_pair, v1, v2, IO);
}

std::vector<std::vector<double>> PropsSImulti(const std::vector<std::string>& Outputs, const std::string& Name1, const std::vector<double>& Prop1,
                                              const std::string& Name2, const std::vector<double>& Prop2, const std::string& backend,
                                              const std::vector<std::string>& fluids, const std::vector<double>& fractions) {
    std::vector<std::vector<double>> IO;

#if !defined(NO_ERROR_CATCHING)
    try {
#endif

        // Call the subfunction that can bubble errors
        _PropsSImulti(Outputs, Name1, Prop1, Name2, Prop2, backend, fluids, fractions, IO);

        // Return the value(s)
        return IO;

#if !defined(NO_ERROR_CATCHING)
    } catch (const std::exception& e) {
        set_error_string(e.what());
#    if defined(PROPSSI_ERROR_STDOUT)
        std::cout << e.what() << std::endl;
#    endif
        if (get_debug_level() > 1) {
            std::cout << e.what() << std::endl;
        }
    } catch (...) {
    }
#endif
    return std::vector<std::vector<double>>();
}
double PropsSI(const std::string& Output, const std::string& Name1, double Prop1, const std::string& Name2, double Prop2, const std::string& Ref) {
#if !defined(NO_ERROR_CATCHING)
    try {
#endif

        // BEGIN OF TRY
        // Here is the real code that is inside the try block

        std::string backend, fluid;
        extract_backend(Ref, backend, fluid);
        std::vector<double> fractions(1, 1.0);
        // extract_fractions checks for has_fractions_in_string / has_solution_concentration; no need to double check
        std::string fluid_string = extract_fractions(fluid, fractions);
        std::vector<std::vector<double>> IO;
        _PropsSImulti(strsplit(Output, '&'), Name1, std::vector<double>(1, Prop1), Name2, std::vector<double>(1, Prop2), backend,
                      strsplit(fluid_string, '&'), fractions, IO);
        if (IO.empty()) {
            throw ValueError(get_global_param_string("errstring").c_str());
        }
        if (IO.size() != 1 || IO[0].size() != 1) {
            throw ValueError(format("output should be 1x1; error was %s", get_global_param_string("errstring").c_str()));
        }

        double val = IO[0][0];

        if (get_debug_level() > 1) {
            std::cout << format("_PropsSI will return %g", val) << std::endl;
        }
        return val;
        // END OF TRY
#if !defined(NO_ERROR_CATCHING)
    } catch (const std::exception& e) {
        set_error_string(
          e.what()
          + format(" : PropsSI(\"%s\",\"%s\",%0.10g,\"%s\",%0.10g,\"%s\")", Output.c_str(), Name1.c_str(), Prop1, Name2.c_str(), Prop2, Ref.c_str()));
#    if defined(PROPSSI_ERROR_STDOUT)
        std::cout << e.what() << std::endl;
#    endif
        if (get_debug_level() > 1) {
            std::cout << e.what() << std::endl;
        }
        return _HUGE;
    } catch (...) {
        return _HUGE;
    }
#endif
}

bool add_fluids_as_JSON(const std::string& backend, const std::string& fluidstring) {
    if (backend == "SRK" || backend == "PR") {
        CubicLibrary::add_fluids_as_JSON(fluidstring);
        return true;
    } else if (backend == "HEOS") {
        JSONFluidLibrary::add_many(fluidstring);
        return true;
    } else if (backend == "PCSAFT") {
        PCSAFTLibrary::add_fluids_as_JSON(fluidstring);
        return true;
    } else {
        throw ValueError(format("You have provided an invalid backend [%s] to add_fluids_as_JSON; valid options are SRK, PR, HEOS", backend.c_str()));
    }
}
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to PropsSI", "[PropsSI]") {
    SECTION("Single state, single output") {
        CHECK(ValidNumber(CoolProp::PropsSI("T", "P", 101325, "Q", 0, "Water")));
    };
    SECTION("Single state, single output, saturation derivative") {
        CHECK(ValidNumber(CoolProp::PropsSI("d(P)/d(T)|sigma", "P", 101325, "Q", 0, "Water")));
    };
    SECTION("Single state, single output, pure incompressible") {
        CHECK(ValidNumber(CoolProp::PropsSI("D", "P", 101325, "T", 300, "INCOMP::DowQ")));
    };
    SECTION("Single state, trivial output, pure incompressible") {
        CHECK(ValidNumber(CoolProp::PropsSI("Tmin", "P", 0, "T", 0, "INCOMP::DowQ")));
    };
    SECTION("Bad input pair") {
        CHECK(!ValidNumber(CoolProp::PropsSI("D", "Q", 0, "Q", 0, "Water")));
    };
    SECTION("Single state, single output, 40% incompressible") {
        CHECK(ValidNumber(CoolProp::PropsSI("D", "P", 101325, "T", 300, "INCOMP::MEG[0.40]")));
    };
    SECTION("Single state, single output, predefined CoolProp mixture") {
        CHECK(ValidNumber(CoolProp::PropsSI("T", "Q", 1, "P", 3e6, "HEOS::R125[0.7]&R32[0.3]")));
    };
    SECTION("Single state, single output") {
        CHECK(ValidNumber(CoolProp::PropsSI("T", "P", 101325, "Q", 0, "HEOS::Water")));
    };
    SECTION("Single state, single output, predefined mixture") {
        CHECK(ValidNumber(CoolProp::PropsSI("T", "P", 101325, "Q", 0, "R410A.mix")));
    };
    SECTION("Single state, single output, predefined mixture from REFPROP") {
        CHECK(ValidNumber(CoolProp::PropsSI("T", "P", 101325, "Q", 0, "REFPROP::R410A.MIX")));
    };
    SECTION("Single state, single output, bad predefined mixture from REFPROP") {
        CHECK(!ValidNumber(CoolProp::PropsSI("T", "P", 101325, "Q", 0, "REFPROP::RRRRRR.mix")));
    };
    SECTION("Predefined mixture") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z;
        std::vector<std::string> outputs(1, "T");
        outputs.push_back("Dmolar");
        std::vector<std::vector<double>> IO;
        std::vector<std::string> fluids(1, "R410A.mix");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
    };
    SECTION("Single state, two outputs") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1, "T");
        outputs.push_back("Dmolar");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
    };
    SECTION("Single state, two bad outputs") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "???????");
        outputs.push_back("?????????");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
        CHECK(IO.size() == 0);
    };
    SECTION("Two states, one output") {
        std::vector<double> p(2, 101325), Q(2, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1, "T");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
    };
    SECTION("Two states, two outputs") {
        std::vector<double> p(2, 101325), Q(2, 1.0), z(1, 1.0);
        std::vector<std::string> outputs(1, "T");
        outputs.push_back("Dmolar");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
    };
    SECTION("cp and its derivative representation") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "Cpmolar");
        outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(!IO.empty());
        CAPTURE(IO[0][0]);
        CAPTURE(IO[0][1]);
        CHECK(std::abs(IO[0][0] - IO[0][1]) < 1e-5);
    };
    SECTION("bad fluid") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "Cpmolar");
        outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "????????");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad mole fraction length") {
        std::vector<double> p(1, 101325), Q(1, 1.0), z(1, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "T");
        std::vector<std::string> fluids(1, "Water&Ethanol");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad input lengths") {
        std::vector<double> p(1, 101325), Q(2, 1.0), z(100, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "Cpmolar");
        outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "P", p, "Q", Q, "HEOS", fluids, z));
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
    SECTION("bad input pair") {
        std::vector<double> Q(2, 1.0), z(1, 1.0);
        std::vector<std::vector<double>> IO;
        std::vector<std::string> outputs(1, "Cpmolar");
        outputs.push_back("d(Hmolar)/d(T)|P");
        std::vector<std::string> fluids(1, "Water");
        CHECK_NOTHROW(IO = CoolProp::PropsSImulti(outputs, "Q", Q, "Q", Q, "HEOS", fluids, z));
        std::string errstring = get_global_param_string("errstring");
        CAPTURE(errstring);
        REQUIRE(IO.empty());
    };
};
#endif

/****************************************************
 *                  Props1SI                        *
 ****************************************************/

double Props1SI(std::string FluidName, std::string Output) {
    bool valid_fluid1 = is_valid_fluid_string(FluidName);
    bool valid_fluid2 = is_valid_fluid_string(Output);
    if (valid_fluid1 && valid_fluid2) {
        set_error_string(format("Both inputs to Props1SI [%s,%s] are valid fluids", Output.c_str(), FluidName.c_str()));
        return _HUGE;
    }
    if (!valid_fluid1 && !valid_fluid2) {
        set_error_string(format("Neither input to Props1SI [%s,%s] is a valid fluid", Output.c_str(), FluidName.c_str()));
        return _HUGE;
    }
    if (!valid_fluid1 && valid_fluid2) {
        // They are backwards, swap
        std::swap(Output, FluidName);
    }

    // First input is the fluid, second input is the input parameter
    double val1 = PropsSI(Output, "", 0, "", 0, FluidName);
    if (!ValidNumber(val1)) {
        set_error_string(format("Unable to use input parameter [%s] in Props1SI for fluid %s; error was %s", Output.c_str(), FluidName.c_str(),
                                get_global_param_string("errstring").c_str()));
        return _HUGE;
    } else {
        return val1;
    }
}

std::vector<std::vector<double>> Props1SImulti(const std::vector<std::string>& Outputs, const std::string& backend, const std::vector<std::string>& fluids, const std::vector<double>& fractions) {
    std::vector<double> zero_vector(1, 0.);
    std::vector<std::vector<double>> val1 = PropsSImulti(Outputs, "", zero_vector, "", zero_vector, backend, fluids, fractions);
    // error handling is done in PropsSImulti, val1 will be an empty vector if an error occured
    return val1;
}


#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to Props1SI", "[Props1SI],[PropsSI]") {
    SECTION("Good fluid, good parameter") {
        CHECK(ValidNumber(CoolProp::Props1SI("Tcrit", "Water")));
    };
    SECTION("Good fluid, good parameter") {
        CHECK(ValidNumber(CoolProp::PropsSI("Tcrit", "", 0, "", 0, "Water")));
    };
    SECTION("Good fluid, good parameter, inverted") {
        CHECK(ValidNumber(CoolProp::Props1SI("Water", "Tcrit")));
    };
    SECTION("Good fluid, bad parameter") {
        CHECK(!ValidNumber(CoolProp::Props1SI("Water", "????????????")));
    };
    SECTION("Bad fluid, good parameter") {
        CHECK(!ValidNumber(CoolProp::Props1SI("?????", "Tcrit")));
    };
};
#endif

bool is_valid_fluid_string(const std::string& input_fluid_string) {
    try {
        std::string backend, fluid;
        std::vector<double> fractions;
        // First try to extract backend and fractions
        extract_backend(input_fluid_string, backend, fluid);
        std::string fluid_string = extract_fractions(fluid, fractions);
        // We are going to let the factory function load the state
        shared_ptr<AbstractState> State(AbstractState::factory(backend, fluid_string));
        return true;
    } catch (...) {
        return false;
    }
}
double saturation_ancillary(const std::string& fluid_name, const std::string& output, int Q, const std::string& input, double value) {

    // Generate the state instance
    std::vector<std::string> names(1, fluid_name);
    CoolProp::HelmholtzEOSMixtureBackend HEOS(names);

    parameters iInput = get_parameter_index(input);
    parameters iOutput = get_parameter_index(output);

    return HEOS.saturation_ancillary(iOutput, Q, iInput, value);
}
void set_reference_stateS(const std::string& fluid_string, const std::string& reference_state) {
    std::string backend, fluid;
    extract_backend(fluid_string, backend, fluid);
    if (backend == "REFPROP") {

        int ierr = 0, ixflag = 1;
        double h0 = 0, s0 = 0, t0 = 0, p0 = 0;
        char herr[255], hrf[4];
        double x0[1] = {1};
        const char* refstate = reference_state.c_str();
        if (strlen(refstate) > 3) {
            if (reference_state == "ASHRAE") {
                strcpy(hrf, "ASH");
            } else {
                throw ValueError(format("Reference state string [%s] is more than 3 characters long", reference_state.c_str()));
            }
        } else {
            strcpy(hrf, refstate);
        }
        REFPROP_SETREF(hrf, ixflag, x0, h0, s0, t0, p0, ierr, herr, 3, 255);
    } else if (backend == "HEOS" || backend == "?") {
        CoolProp::HelmholtzEOSMixtureBackend HEOS(std::vector<std::string>(1, fluid));
        if (!reference_state.compare("IIR")) {
            if (HEOS.Ttriple() > 273.15) {
                throw ValueError(format("Cannot use IIR reference state; Ttriple [%Lg] is greater than 273.15 K", HEOS.Ttriple()));
            }
            HEOS.update(QT_INPUTS, 0, 273.15);

            // Get current values for the enthalpy and entropy
            double deltah = HEOS.hmass() - 200000;  // offset from 200000 J/kg enthalpy
            double deltas = HEOS.smass() - 1000;    // offset from 1000 J/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(fluid, delta_a1, delta_a2, "IIR");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("ASHRAE")) {
            if (HEOS.Ttriple() > 233.15) {
                throw ValueError(format("Cannot use ASHRAE reference state; Ttriple [%Lg] is greater than than 233.15 K", HEOS.Ttriple()));
            }
            HEOS.update(QT_INPUTS, 0, 233.15);

            // Get current values for the enthalpy and entropy
            double deltah = HEOS.hmass() - 0;  // offset from 0 J/kg enthalpy
            double deltas = HEOS.smass() - 0;  // offset from 0 J/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(fluid, delta_a1, delta_a2, "ASHRAE");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("NBP")) {
            if (HEOS.p_triple() > 101325) {
                throw ValueError(format("Cannot use NBP reference state; p_triple [%Lg Pa] is greater than than 101325 Pa", HEOS.p_triple()));
            }
            // Saturated liquid boiling point at 1 atmosphere
            HEOS.update(PQ_INPUTS, 101325, 0);

            double deltah = HEOS.hmass() - 0;  // offset from 0 kJ/kg enthalpy
            double deltas = HEOS.smass() - 0;  // offset from 0 kJ/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(fluid, delta_a1, delta_a2, "NBP");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("DEF")) {
            set_fluid_enthalpy_entropy_offset(fluid, 0, 0, "DEF");
        } else if (!reference_state.compare("RESET")) {
            set_fluid_enthalpy_entropy_offset(fluid, 0, 0, "RESET");
        } else {
            throw ValueError(format("Reference state string is invalid: [%s]", reference_state.c_str()));
        }
    }
}
void set_reference_stateD(const std::string& Ref, double T, double rhomolar, double hmolar0, double smolar0) {
    std::vector<std::string> _comps(1, Ref);
    CoolProp::HelmholtzEOSMixtureBackend HEOS(_comps);

    HEOS.update(DmolarT_INPUTS, rhomolar, T);

    // Get current values for the enthalpy and entropy
    double deltah = HEOS.hmolar() - hmolar0;  // offset from specified enthalpy in J/mol
    double deltas = HEOS.smolar() - smolar0;  // offset from specified entropy in J/mol/K
    double delta_a1 = deltas / (HEOS.gas_constant());
    double delta_a2 = -deltah / (HEOS.gas_constant() * HEOS.get_reducing_state().T);
    set_fluid_enthalpy_entropy_offset(Ref, delta_a1, delta_a2, "custom");
}

std::string get_global_param_string(const std::string& ParamName) {
    if (!ParamName.compare("version")) {
        return version;
    } else if (!ParamName.compare("gitrevision")) {
        return gitrevision;
    } else if (!ParamName.compare("errstring")) {
        std::string temp = error_string;
        error_string = "";
        return temp;
    } else if (!ParamName.compare("warnstring")) {
        std::string temp = warning_string;
        warning_string = "";
        return temp;
    } else if (!ParamName.compare("FluidsList") || !ParamName.compare("fluids_list") || !ParamName.compare("fluidslist")) {
        return get_fluid_list();
    } else if (!ParamName.compare("incompressible_list_pure")) {
        return get_incompressible_list_pure();
    } else if (!ParamName.compare("incompressible_list_solution")) {
        return get_incompressible_list_solution();
    } else if (!ParamName.compare("mixture_binary_pairs_list")) {
        return get_csv_mixture_binary_pairs();
    } else if (!ParamName.compare("parameter_list")) {
        return get_csv_parameter_list();
    } else if (!ParamName.compare("predefined_mixtures")) {
        return get_csv_predefined_mixtures();
    } else if (!ParamName.compare("HOME")) {
        return get_home_dir();
    } else if (ParamName == "REFPROP_version") {
        return REFPROPMixtureBackend::version();
    } else if (ParamName == "cubic_fluids_schema") {
        return CoolProp::CubicLibrary::get_cubic_fluids_schema();
    } else if (ParamName == "cubic_fluids_list") {
        return CoolProp::CubicLibrary::get_cubic_fluids_list();
    } else if (ParamName == "pcsaft_fluids_schema") {
        return CoolProp::PCSAFTLibrary::get_pcsaft_fluids_schema();
    } else {
        throw ValueError(format("Input parameter [%s] is invalid", ParamName.c_str()));
    }
};
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to get_global_param_string", "[get_global_param_string]") {
    const int num_good_inputs = 8;
    std::string good_inputs[num_good_inputs] = {
      "version",        "gitrevision",        "fluids_list", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list",
      "parameter_list", "predefined_mixtures"};
    std::ostringstream ss3c;
    for (int i = 0; i < num_good_inputs; ++i) {
        ss3c << "Test for" << good_inputs[i];
        SECTION(ss3c.str(), "") {
            CHECK_NOTHROW(CoolProp::get_global_param_string(good_inputs[i]));
        };
    }
    CHECK_THROWS(CoolProp::get_global_param_string(""));
};
#endif
std::string get_fluid_param_string(const std::string& FluidName, const std::string& ParamName) {
    std::string backend, fluid;
    extract_backend(FluidName, backend, fluid);
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, fluid));
    return AS->fluid_param_string(ParamName);
}
#if defined(ENABLE_CATCH)
TEST_CASE("Check inputs to get_fluid_param_string", "[get_fluid_param_string]") {
    const int num_good_inputs = 10;
    std::string good_inputs[num_good_inputs] = {"aliases",
                                                "CAS",
                                                "ASHRAE34",
                                                "REFPROPName",
                                                "BibTeX-CONDUCTIVITY",
                                                "BibTeX-EOS",
                                                "BibTeX-CP0",
                                                "BibTeX-SURFACE_TENSION",
                                                "BibTeX-MELTING_LINE",
                                                "BibTeX-VISCOSITY"};
    std::ostringstream ss3c;
    for (int i = 0; i < num_good_inputs; ++i) {
        ss3c << "Test for" << good_inputs[i];
        SECTION(ss3c.str(), "") {
            CHECK_NOTHROW(CoolProp::get_fluid_param_string("Water", good_inputs[i]));
        };
    }
    CHECK_THROWS(CoolProp::get_fluid_param_string("", "aliases"));
    CHECK_THROWS(CoolProp::get_fluid_param_string("Water", ""));
    CHECK_THROWS(CoolProp::get_fluid_param_string("Water", "BibTeX-"));
    CHECK(CoolProp::get_fluid_param_string("Water", "pure") == "true");
    CHECK(CoolProp::get_fluid_param_string("R410A", "pure") == "false");
};
#endif

std::string phase_lookup_string(phases Phase) {
    switch (Phase) {
        case iphase_liquid:  ///< Liquid
            return "liquid";
        case iphase_supercritical:  ///< Supercritical (p > pc, T > Tc)
            return "supercritical";
        case iphase_supercritical_gas:  ///< Supercritical gas (p < pc, T > Tc)
            return "supercritical_gas";
        case iphase_supercritical_liquid:  ///< Supercritical liquid (p > pc, T < Tc)
            return "supercritical_liquid";
        case iphase_critical_point:  ///< At the critical point
            return "critical_point";
        case iphase_gas:  ///< Subcritical gas
            return "gas";
        case iphase_twophase:  ///< Twophase (between saturation curves - inclusive)
            return "twophase";
        case iphase_unknown:  ///< Unknown phase
            return "unknown";
        case iphase_not_imposed:
            return "not_imposed";
    }
    throw ValueError("I should never be thrown");
}
std::string PhaseSI(const std::string& Name1, double Prop1, const std::string& Name2, double Prop2, const std::string& FluidName) {
    double Phase_double = PropsSI("Phase", Name1, Prop1, Name2, Prop2, FluidName);  // Attempt to get "Phase" from PropsSI()
    if (!ValidNumber(Phase_double)) {                                               // if the returned phase is invalid...
        std::string strPhase = phase_lookup_string(iphase_unknown);                 //     phase is unknown.
        std::string strError = get_global_param_string("errstring").c_str();        //     fetch waiting error string
        if (strError != "") {                                                       //     if error string is not empty,
            strPhase.append(": " + strError);                                       //        append it to the phase string.
        }
        return strPhase;                                             //     return the "unknown" phase string
    }                                                                // else
    std::size_t Phase_int = static_cast<std::size_t>(Phase_double);  //     convert returned phase to int
    return phase_lookup_string(static_cast<phases>(Phase_int));      //     return phase as a string
}

/*
std::string PhaseSI(const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName, const std::vector<double> &z)
{
    double Phase_double = PropsSI("Phase",Name1,Prop1,Name2,Prop2,FluidName,z);
    if (!ValidNumber(Phase_double)){ return "";}
    std::size_t Phase_int = static_cast<std::size_t>(Phase_double);
    return phase_lookup_string(static_cast<phases>(Phase_int));
}
*/
} /* namespace CoolProp */
