
/*!
From REFPROP:
temperature                     K
pressure, fugacity              kPa
density                         mol/L
composition                     mole fraction
quality                         mole basis (moles vapor/total moles)
enthalpy, internal energy       J/mol
Gibbs, Helmholtz free energy    J/mol
entropy, heat capacity          J/(mol.K)
speed of sound                  m/s
Joule-Thomson coefficient       K/kPa
d(p)/d(rho)                     kPa.L/mol
d2(p)/d(rho)2                   kPa.(L/mol)^2
viscosity                       microPa.s (10^-6 Pa.s)
thermal conductivity            W/(m.K)
dipole moment                   debye
surface tension                 N/m
*/

#define REFPROP_IMPLEMENTATION
#define REFPROP_CSTYLE_REFERENCES
#include "externals/REFPROP-headers/REFPROP_lib.h"
#undef REFPROP_IMPLEMENTATION
#undef REFPROP_CSTYLE_REFERENCES

#include "CoolPropTools.h"
#include "REFPROPMixtureBackend.h"
#include "REFPROPBackend.h"
#include "Exceptions.h"
#include "Configuration.h"
#include "CPfilepaths.h"
#include "CoolProp.h"
#include "Solvers.h"
#include "IdealCurves.h"
#include "DataStructures.h"
#include "AbstractState.h"

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include "crossplatform_shared_ptr.h"
#include <stdlib.h>

#if defined(_MSC_VER)
#    define _CRTDBG_MAP_ALLOC
#    ifndef _CRT_SECURE_NO_WARNINGS
#        define _CRT_SECURE_NO_WARNINGS
#    endif
#    include <crtdbg.h>
#    include <sys/stat.h>
#else
#    include <sys/stat.h>
#endif

std::string LoadedREFPROPRef;

static bool dbg_refprop = false;
static const unsigned int number_of_endings = 5;
std::string endings[number_of_endings] = {"", ".FLD", ".fld", ".PPF", ".ppf"};

static char default_reference_state[] = "DEF";

// Default location, can be over-ridden by configuration variable
#if defined(__powerpc__) || defined(__ISLINUX__) || defined(__ISAPPLE__)
char refpropPath[] = "/opt/refprop";
#elif defined(__ISWINDOWS__)
char refpropPath[] = "";
#else
#    pragma error
#endif

/// Find either FLUIDS or fluids folder relative to the root path provided; return the path
std::string get_casesensitive_fluids(const std::string& root) {
    std::string joined = join_path(root, "fluids");
    if (path_exists(joined)) {
        return joined;
    } else {
        std::string ucase_joined = join_path(root, "FLUIDS");
        if (path_exists(ucase_joined)) {
            return ucase_joined;
        } else {
            throw CoolProp::ValueError(format("fluid directories \"FLUIDS\" or \"fluids\" could not be found in the directory [%s]", root));
        }
    }
}
std::string get_REFPROP_fluid_path_prefix() {
    std::string rpPath = refpropPath;
    // Allow the user to specify an alternative REFPROP path by configuration value
    std::string alt_refprop_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_PATH);
    if (!alt_refprop_path.empty()) {
        // The alternative path has been set, so we give all fluid paths as relative to this directory
        //if (!endswith(alt_refprop_path, separator)){
        //    throw CoolProp::ValueError(format("ALTERNATIVE_REFPROP_PATH [%s] must end with a path sparator, typically a slash character", alt_refprop_path.c_str()));
        //}
        if (!path_exists(alt_refprop_path)) {
            throw CoolProp::ValueError(format("ALTERNATIVE_REFPROP_PATH [%s] could not be found", alt_refprop_path.c_str()));
        }
        return get_casesensitive_fluids(alt_refprop_path);
    }
#if defined(__ISWINDOWS__)
    return rpPath;
#elif defined(__ISLINUX__) || defined(__ISAPPLE__)
    return get_casesensitive_fluids(rpPath);
#else
    throw CoolProp::NotImplementedError("This function should not be called.");
    return rpPath;
#endif
}
std::string get_REFPROP_mixtures_path_prefix() {
    std::string rpPath = refpropPath;
    // Allow the user to specify an alternative REFPROP path by configuration value
    std::string alt_refprop_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_PATH);
    std::string separator = get_separator();
    if (!alt_refprop_path.empty()) {
        //if (!endswith(alt_refprop_path, separator)) {
        //    throw CoolProp::ValueError(format("ALTERNATIVE_REFPROP_PATH [%s] must end with a path sparator, typically a slash character", alt_refprop_path.c_str()));
        //}
        if (!path_exists(alt_refprop_path)) {
            throw CoolProp::ValueError(format("ALTERNATIVE_REFPROP_PATH [%s] could not be found", alt_refprop_path.c_str()));
        }
        // The alternative path has been set
        return join_path(alt_refprop_path, "mixtures");
    }
#if defined(__ISWINDOWS__)
    return rpPath;
#elif defined(__ISLINUX__) || defined(__ISAPPLE__)
    return join_path(rpPath, "mixtures");
#else
    throw CoolProp::NotImplementedError("This function should not be called.");
    return rpPath;
#endif
}

/// Construct the path to the HMX.BNC file
std::string get_REFPROP_HMX_BNC_path() {
    std::string alt_hmx_bnc_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_HMX_BNC_PATH);
    // Use the alternative HMX.BNC path if provided - replace all the path to HMX.BNC with provided path
    if (!alt_hmx_bnc_path.empty()) {
        return alt_hmx_bnc_path;
    } else {
        // Otherwise fall back to default paths; get_REFPROP_fluid_path_prefix will query ALTERNATIVE_REFPROP_PATH
        return join_path(get_REFPROP_fluid_path_prefix(), "HMX.BNC");
    }
}

namespace CoolProp {

class REFPROPGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        REFPROPMixtureBackend::REFPROP_supported();
        if (fluid_names.size() == 1) {
            return new REFPROPBackend(fluid_names[0]);
        } else {
            return new REFPROPMixtureBackend(fluid_names);
        }
    };
};
// This static initialization will cause the generator to register
static GeneratorInitializer<REFPROPGenerator> refprop_gen(REFPROP_BACKEND_FAMILY);

void REFPROPMixtureBackend::construct(const std::vector<std::string>& fluid_names) {
    // Do the REFPROP instantiation for this fluid
    _mole_fractions_set = false;

    // Force loading of REFPROP
    REFPROP_supported();

    // Try to add this fluid to REFPROP - might want to think about making array of
    // components and setting mole fractions if they change a lot.
    this->set_REFPROP_fluids(fluid_names);

    // Bump the number of REFPROP backends that are in existence;
    REFPROPMixtureBackend::instance_counter++;

    // Set the imposed phase index to default value
    imposed_phase_index = iphase_not_imposed;
}

REFPROPMixtureBackend::~REFPROPMixtureBackend() {
    // Decrement the counter for the number of instances
    REFPROPMixtureBackend::instance_counter--;
    // Unload the shared library when the last instance is about to be destroyed
    if (REFPROPMixtureBackend::instance_counter == 0) {
        force_unload_REFPROP();
    }
}
void REFPROPMixtureBackend::check_loaded_fluid() {
    this->set_REFPROP_fluids(this->fluid_names);
}

std::size_t REFPROPMixtureBackend::instance_counter = 0;  // initialise with 0
bool REFPROPMixtureBackend::_REFPROP_supported = true;    // initialise with true
bool REFPROPMixtureBackend::REFPROP_supported() {
    /*
     * Here we build the bridge from macro definitions
     * into the actual code. This is also going to be
     * the central place to handle error messages on
     * unsupported platforms.
     */

    // Abort check if Refprop has been loaded.
    if (RefpropdllInstance != NULL) return true;

    // Store result of previous check.
    if (_REFPROP_supported) {
        // Either Refprop is supported or it is the first check.
        std::string rpv(STRINGIFY(RPVersion));
        if (rpv.compare("NOTAVAILABLE") != 0) {
            // Function names were defined in "REFPROP_lib.h",
            // This platform theoretically supports Refprop.
            std::string err;
            const std::string alt_rp_path = get_config_string(ALTERNATIVE_REFPROP_PATH);
            const std::string alt_rp_name = get_config_string(ALTERNATIVE_REFPROP_LIBRARY_PATH);
            bool loaded_REFPROP = false;
            if (!alt_rp_name.empty()) {
                loaded_REFPROP = ::load_REFPROP(err, "", alt_rp_name);
            } else {
                if (alt_rp_path.empty()) {
                    loaded_REFPROP = ::load_REFPROP(err, refpropPath, "");
                } else {
                    loaded_REFPROP = ::load_REFPROP(err, alt_rp_path, "");
                }
            }

            if (loaded_REFPROP) {
                return true;
            } else {
                printf("Good news: It is possible to use REFPROP on your system! However, the library \n");
                printf("could not be loaded. Please make sure that REFPROP is available on your system.\n\n");
                printf("Neither found in current location nor found in system PATH.\n");
                printf("If you already obtained a copy of REFPROP from http://www.nist.gov/srd/, \n");
                printf("add location of REFPROP to the PATH environment variable or your library path.\n\n");
                printf("In case you do not use Windows, have a look at https://github.com/jowr/librefprop.so \n");
                printf("to find instructions on how to compile your own version of the REFPROP library.\n\n");
                printf("ALTERNATIVE_REFPROP_PATH: %s\n", alt_rp_path.c_str());
                printf("ERROR: %s\n", err.c_str());
                _REFPROP_supported = false;
                return false;
            }
        } else {
            // No definition of function names, we do not expect
            // the Refprop library to be available.
            _REFPROP_supported = false;
            return false;
        }
    } else {
        return false;
    }
    return false;
}
std::string REFPROPMixtureBackend::version() {
    int N = -1;
    int ierr = 0;
    char fluids[10000] = "", hmx[] = "HMX.BNC", default_reference_state[] = "DEF", herr[255] = "";
    if (!REFPROP_supported()) {
        return "n/a";
    };
    // Pad the version string with NULL characters
    for (int i = 0; i < 255; ++i) {
        herr[i] = '\0';
    }
    SETUPdll(&N, fluids, hmx, default_reference_state, &ierr, herr,
             10000,              // Length of component_string (see PASS_FTN.for from REFPROP)
             refpropcharlength,  // Length of path_HMX_BNC
             lengthofreference,  // Length of reference
             errormessagelength  // Length of error message
    );
    if (strlen(herr) == 0) {
        return format("%g", ((double)ierr) / 10000.0);
    } else {
        std::string s(herr, herr + 254);
        return strstrip(s);
    }
}

void REFPROPMixtureBackend::set_REFPROP_fluids(const std::vector<std::string>& fluid_names) {
    // If the name of the refrigerant doesn't match
    // that of the currently loaded refrigerant, fluids must be loaded
    if (!cached_component_string.empty() && LoadedREFPROPRef == cached_component_string) {
        if (CoolProp::get_debug_level() > 5) {
            std::cout << format("%s:%d: The current fluid can be reused; %s and %s match \n", __FILE__, __LINE__, cached_component_string.c_str(),
                                LoadedREFPROPRef.c_str());
        }
        if (dbg_refprop)
            std::cout << format("%s:%d: The current fluid can be reused; %s and %s match \n", __FILE__, __LINE__, cached_component_string.c_str(),
                                LoadedREFPROPRef.c_str());
        int N = static_cast<int>(this->fluid_names.size());
        if (N > ncmax) {
            throw ValueError(format("Size of fluid vector [%d] is larger than the maximum defined by REFPROP [%d]", fluid_names.size(), ncmax));
        }
        // this->Ncomp = N; ( this should not get set because it is already set and is always 1 for predefined mixtures )
        mole_fractions.resize(ncmax);
        mole_fractions_liq.resize(ncmax);
        mole_fractions_vap.resize(ncmax);
        return;
    } else {
        int ierr = 0;
        this->fluid_names = fluid_names;
        char component_string[10000], herr[errormessagelength];
        std::string components_joined = strjoin(fluid_names, "|");
        std::string components_joined_raw = strjoin(fluid_names, "|");
        std::string fdPath = get_REFPROP_fluid_path_prefix();
        int N = static_cast<int>(fluid_names.size());

        // Get path to HMX.BNC file
        char hmx_bnc[255];
        const std::string HMX_path = get_REFPROP_HMX_BNC_path();
        const char* _HMX_path = HMX_path.c_str();
        if (strlen(_HMX_path) > refpropcharlength) {
            throw ValueError(format("Full HMX path (%s) is too long; max length is 255 characters", _HMX_path));
        }
        strcpy(hmx_bnc, _HMX_path);

        if (get_config_bool(REFPROP_USE_GERG)) {
            int iflag = 1,  // Tell REFPROP to use GERG04; 0 unsets GERG usage
              ierr = 0;
            char herr[255];
            GERG04dll(&N, &iflag, &ierr, herr, 255);
        }

        // Check platform support
        if (!REFPROP_supported()) {
            throw NotImplementedError("You cannot use the REFPROPMixtureBackend.");
        }

        if (N == 1 && upper(components_joined_raw).find(".MIX") != std::string::npos) {
            // It's a predefined mixture
            ierr = 0;
            std::vector<double> x(ncmax);
            char mix[255], reference_state[4] = "DEF";

            std::string path_to_MIX_file = join_path(get_REFPROP_mixtures_path_prefix(), components_joined_raw);
            const char* _components_joined_raw = path_to_MIX_file.c_str();
            if (strlen(_components_joined_raw) > 255) {
                throw ValueError(format("components (%s) is too long", components_joined_raw.c_str()));
            }
            strcpy(mix, _components_joined_raw);

            SETMIXdll(mix, hmx_bnc, reference_state, &N, component_string, &(x[0]), &ierr, herr, 255, 255, 3, 10000, 255);
            if (static_cast<int>(ierr) <= 0) {
                this->Ncomp = N;
                mole_fractions.resize(ncmax);
                mole_fractions_liq.resize(ncmax);
                mole_fractions_vap.resize(ncmax);
                LoadedREFPROPRef = mix;
                cached_component_string = mix;
                this->fluid_names.clear();
                this->fluid_names.push_back(components_joined_raw);
                if (CoolProp::get_debug_level() > 5) {
                    std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n", __FILE__, __LINE__, components_joined.c_str());
                }
                if (dbg_refprop) std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n", __FILE__, __LINE__, components_joined.c_str());
                if (get_config_bool(REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS) && ierr == -117) {
                    throw ValueError(format("Interaction parameter estimation has been disabled: %s", herr));
                }
                set_mole_fractions(std::vector<CoolPropDbl>(x.begin(), x.begin() + N));
                if (get_config_bool(REFPROP_USE_PENGROBINSON)) {
                    int iflag = 2;  // Tell REFPROP to use Peng-Robinson;
                    PREOSdll(&iflag);
                } else {
                    int iflag = 0;  // Tell REFPROP to use normal Helmholtz models
                    PREOSdll(&iflag);
                }
                return;
            } else {
                if (get_debug_level() > 0) {
                    std::cout << format("%s:%d Unable to load predefined mixture [%s] with ierr: [%d] and herr: [%s]\n", __FILE__, __LINE__, mix,
                                        ierr, herr);
                }
                throw ValueError(format("Unable to load mixture: %s", components_joined_raw.c_str()));
            }
        }

        // Loop over the file names - first we try with nothing, then .fld, then .FLD, then .ppf - means you can't mix and match
        for (unsigned int k = 0; k < number_of_endings; k++) {
            // Build the mixture string
            for (unsigned int j = 0; j < (unsigned int)N; j++) {
                if (j == 0) {
                    components_joined = join_path(fdPath, upper(fluid_names[j]) + endings[k]);
                } else {
                    components_joined += "|" + join_path(fdPath, upper(fluid_names[j]) + endings[k]);
                }
            }

            if (dbg_refprop)
                std::cout << format("%s:%d: The fluid %s has not been loaded before, current value is %s \n", __FILE__, __LINE__,
                                    components_joined_raw.c_str(), LoadedREFPROPRef.c_str());

            // Copy over the list of components
            const char* _components_joined = components_joined.c_str();
            if (strlen(_components_joined) > 10000) {
                throw ValueError(format("components_joined (%s) is too long", _components_joined));
            }
            strcpy(component_string, _components_joined);
            // Pad the fluid string all the way to 10k characters with spaces to deal with string parsing bug in REFPROP in SETUPdll
            for (int i = static_cast<int>(components_joined.size()); i < 10000; ++i) {
                component_string[i] = ' ';
            }

            ierr = 0;
            //...Call SETUP to initialize the program
            SETUPdll(&N, component_string, hmx_bnc, default_reference_state, &ierr, herr,
                     10000,              // Length of component_string (see PASS_FTN.for from REFPROP)
                     refpropcharlength,  // Length of path_HMX_BNC
                     lengthofreference,  // Length of reference
                     errormessagelength  // Length of error message
            );
            if (get_config_bool(REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS) && ierr == -117) {
                throw ValueError(format("Interaction parameter estimation has been disabled: %s", herr));
            }
            if (get_config_bool(REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS) && ierr == 117) {
                ierr = 0;
            }
            if (static_cast<int>(ierr) <= 0)  // Success (or a warning, which is silently squelched for now)
            {
                this->Ncomp = N;
                mole_fractions.resize(ncmax);
                mole_fractions_liq.resize(ncmax);
                mole_fractions_vap.resize(ncmax);
                LoadedREFPROPRef = _components_joined;
                cached_component_string = _components_joined;
                if (CoolProp::get_debug_level() > 5) {
                    std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n", __FILE__, __LINE__, components_joined.c_str());
                }
                if (dbg_refprop) std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n", __FILE__, __LINE__, components_joined.c_str());

                if (get_config_bool(REFPROP_USE_PENGROBINSON)) {
                    int iflag = 2;  // Tell REFPROP to use Peng-Robinson;
                    PREOSdll(&iflag);
                } else {
                    int iflag = 0;  // Tell REFPROP to use normal Helmholtz models
                    PREOSdll(&iflag);
                }
                return;
            } else if (k < number_of_endings - 1) {  // Keep going
                if (CoolProp::get_debug_level() > 5) {
                    std::cout << format("REFPROP error/warning [ierr: %d]: %s", ierr, herr) << std::endl;
                }
                continue;
            } else {
                if (CoolProp::get_debug_level() > 5) {
                    std::cout << format("k: %d #endings: %d", k, number_of_endings) << std::endl;
                }
                throw ValueError(format("Could not load these fluids: %s", components_joined_raw.c_str()));
            }
        }
    }
}
std::string REFPROPMixtureBackend::fluid_param_string(const std::string& ParamName) {
    if (ParamName == "CAS") {
        //        subroutine NAME (icomp,hnam,hn80,hcasn)
        //        c
        //        c  provides name information for specified component
        //            c
        //            c  input:
        //            c    icomp--component number in mixture; 1 for pure fluid
        //                c  outputs:
        //                c     hnam--component name [character*12]
        //                c     hn80--component name--long form [character*80]
        //                c    hcasn--CAS (Chemical Abstracts Service) number [character*12]
        std::vector<std::string> CASvec;
        for (int icomp = 1L; icomp <= static_cast<int>(fluid_names.size()); ++icomp) {
            char hnam[13], hn80[81], hcasn[13];
            NAMEdll(&icomp, hnam, hn80, hcasn, 12, 80, 12);
            hcasn[12] = '\0';
            std::string casn = hcasn;
            strstrip(casn);
            CASvec.push_back(casn);
        }
        return strjoin(CASvec, "&");
    } else if (ParamName == "name") {
        int icomp = 1L;
        char hnam[13], hn80[81], hcasn[13];
        NAMEdll(&icomp, hnam, hn80, hcasn, 12, 80, 12);
        hnam[12] = '\0';
        std::string name = hnam;
        strstrip(name);
        return name;
    } else if (ParamName == "long_name") {
        int icomp = 1L;
        char hnam[13], hn80[81], hcasn[13];
        NAMEdll(&icomp, hnam, hn80, hcasn, 12, 80, 12);
        hn80[80] = '\0';
        std::string n80 = hn80;
        strstrip(n80);
        return n80;
    } else {
        throw ValueError(format("parameter to fluid_param_string is invalid: %s", ParamName.c_str()));
    }
};
int REFPROPMixtureBackend::match_CAS(const std::string& CAS) {
    for (int icomp = 1L; icomp <= static_cast<int>(fluid_names.size()); ++icomp) {
        char hnam[13], hn80[81], hcasn[13];
        NAMEdll(&icomp, hnam, hn80, hcasn, 12, 80, 12);
        hcasn[12] = '\0';
        std::string casn = hcasn;
        strstrip(casn);
        if (casn == CAS) {
            return icomp;
        }
    }
    throw ValueError(format("Unable to match CAS number [%s]", CAS.c_str()));
}
/// Set binary mixture floating point parameter
void REFPROPMixtureBackend::set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter,
                                                          const double value) {
    std::size_t i = match_CAS(CAS1) - 1, j = match_CAS(CAS2) - 1;
    return set_binary_interaction_double(i, j, parameter, value);
};
/// Get binary mixture double value
double REFPROPMixtureBackend::get_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) {
    std::size_t i = match_CAS(CAS1) - 1, j = match_CAS(CAS2) - 1;
    return get_binary_interaction_double(i, j, parameter);
}
/// Get binary mixture string value
std::string REFPROPMixtureBackend::get_binary_interaction_string(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) {

    int icomp, jcomp;
    char hmodij[4], hfmix[255], hbinp[255], hfij[255], hmxrul[255];
    double fij[6];

    icomp = match_CAS(CAS1);
    jcomp = match_CAS(CAS2);

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, hfij, hbinp, hmxrul, 3, 255, 255, 255, 255);

    std::string shmodij(hmodij);
    if (shmodij.find("KW") == 0 || shmodij.find("GE") == 0)  // Starts with KW or GE
    {
        if (parameter == "model") {
            return shmodij;
        } else {
            throw ValueError(format(" I don't know what to do with your parameter [%s]", parameter.c_str()));
            return "";
        }
    } else {
        //throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij));
        return "";
    }
}
/// Set binary mixture string value
void REFPROPMixtureBackend::set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                          const std::string& value) {
    // bound-check indices
    if (i < 0 || i >= Ncomp) {
        if (j < 0 || j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp-1));
        }
    } else if (j < 0 || j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp-1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1, ierr = 0L;
    char hmodij[4], hfmix[255], hbinp[255], hfij[255], hmxrul[255];
    double fij[6];
    char herr[255];

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, hfij, hbinp, hmxrul, 3, 255, 255, 255, 255);

    if (parameter == "model") {
        if (value.length() > 4) {
            throw ValueError(format("Model parameter (%s) is longer than 4 characters.", value));
        } else {
            strcpy(hmodij, value.c_str());
        }
    } else {
        throw ValueError(format("I don't know what to do with your parameter [%s]", parameter.c_str()));
    }
    SETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, &ierr, herr, 3, 255, 255);
    if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("Unable to set parameter[%s] to value[%s]: %s", parameter.c_str(), value.c_str(), herr));
    }
}
/// Set binary mixture string parameter (EXPERT USE ONLY!!!)
void REFPROPMixtureBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                          const double value) {
    // bound-check indices
    if (i < 0 || i >= Ncomp) {
        if (j < 0 || j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp-1));
        }
    } else if (j < 0 || j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp-1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1, ierr = 0L;
    char hmodij[4], hfmix[255], hbinp[255], hfij[255], hmxrul[255];
    double fij[6];
    char herr[255];

    // Get the prior state
    GETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, hfij, hbinp, hmxrul, 3, 255, 255, 255, 255);

    std::string shmodij(hmodij);
    if (shmodij.find("KW") == 0 || shmodij.find("GE") == 0)  // Starts with KW or GE
    {
        if (parameter == "betaT") {
            fij[0] = value;
        } else if (parameter == "gammaT") {
            fij[1] = value;
        } else if (parameter == "betaV") {
            fij[2] = value;
        } else if (parameter == "gammaV") {
            fij[3] = value;
        } else if (parameter == "Fij") {
            fij[4] = value;
        } else {
            throw ValueError(format("I don't know what to do with your parameter [%s]", parameter.c_str()));
        }
        SETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, &ierr, herr, 3, 255, 255);
        if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("Unable to set parameter[%s] to value[%g]: %s", parameter.c_str(), value, herr));
        }
    } else {
        throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij));
    }
}

/// Get binary mixture double value (EXPERT USE ONLY!!!)
double REFPROPMixtureBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) {
    // bound-check indices
    if (i < 0 || i >= Ncomp) {
        if (j < 0 || j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp-1));
        }
    } else if (j < 0 || j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp-1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1;
    char hmodij[4], hfmix[255], hbinp[255], hfij[255], hmxrul[255];
    double fij[6];

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij, fij, hfmix, hfij, hbinp, hmxrul, 3, 255, 255, 255, 255);

    std::string shmodij(hmodij);
    if (shmodij.find("KW") == 0 || shmodij.find("GE") == 0)  // Starts with KW or GE
    {
        double val;
        if (parameter == "betaT") {
            val = fij[0];
        } else if (parameter == "gammaT") {
            val = fij[1];
        } else if (parameter == "betaV") {
            val = fij[2];
        } else if (parameter == "gammaV") {
            val = fij[3];
        } else if (parameter == "Fij") {
            val = fij[4];
        } else {
            throw ValueError(format(" I don't know what to do with your parameter [%s]", parameter.c_str()));
            return _HUGE;
        }
        return val;
    } else {
        //throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij));
        return _HUGE;
    }
}
void REFPROPMixtureBackend::set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
    if (mole_fractions.size() != this->Ncomp) {
        throw ValueError(
          format("Size of mole fraction vector [%d] does not equal that of component vector [%d]", mole_fractions.size(), this->Ncomp));
    }
    this->mole_fractions = std::vector<CoolPropDbl>(ncmax, 0.0);
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
        this->mole_fractions[i] = static_cast<double>(mole_fractions[i]);
    }
    this->mole_fractions_long_double = mole_fractions;  // same size as Ncomp
    _mole_fractions_set = true;
    clear_comp_change();
}
void REFPROPMixtureBackend::set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
    if (mass_fractions.size() != this->Ncomp) {
        throw ValueError(
          format("size of mass fraction vector [%d] does not equal that of component vector [%d]", mass_fractions.size(), this->Ncomp));
    }
    std::vector<double> moles(this->Ncomp);
    double sum_moles = 0.0;
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    for (int i = 1L; i <= static_cast<int>(this->Ncomp); ++i) {
        INFOdll(&i, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        moles[i - 1] = static_cast<double>(mass_fractions[i - 1]) / (wmm / 1000.0);
        sum_moles += moles[i - 1];
    }
    for (std::size_t i = 0; i < this->Ncomp; ++i) {
        moles[i] = moles[i] / sum_moles;
    }
    this->set_mole_fractions(moles);
};
void REFPROPMixtureBackend::check_status(void) {
    if (!_mole_fractions_set) {
        throw ValueError("Mole fractions not yet set");
    }
}

void REFPROPMixtureBackend::limits(double& Tmin, double& Tmax, double& rhomolarmax, double& pmax) {
    /*
     *
          subroutine LIMITS (htyp,x,tmin,tmax,Dmax,pmax)
    c
    c  returns limits of a property model as a function of composition
    c
    c  Pure fluid limits are read in from the .fld files; for mixtures, a
    c  simple mole fraction weighting in reduced variables is used.
    c
    c  inputs:
    c     htyp--flag indicating which models are to be checked [character*3]
    c           'EOS':  equation of state for thermodynamic properties
    c           'ETA':  viscosity
    c           'TCX':  thermal conductivity
    c           'STN':  surface tension
    c        x--composition array [mol frac]
    c  outputs:
    c     tmin--minimum temperature for model specified by htyp [K]
    c     tmax--maximum temperature [K]
    c     Dmax--maximum density [mol/L]
    c     pmax--maximum pressure [kPa]
     *
     */
    this->check_loaded_fluid();
    double Dmax_mol_L, pmax_kPa;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &Tmin, &Tmax, &Dmax_mol_L, &pmax_kPa, 3);
    pmax = pmax_kPa * 1000;
    rhomolarmax = Dmax_mol_L * 1000;
}
CoolPropDbl REFPROPMixtureBackend::calc_pmax(void) {
    double Tmin, Tmax, rhomolarmax, pmax;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(pmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_Tmax(void) {
    double Tmin, Tmax, rhomolarmax, pmax;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(Tmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_Tmin(void) {
    double Tmin, Tmax, rhomolarmax, pmax;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(Tmin);
};
CoolPropDbl REFPROPMixtureBackend::calc_T_critical() {
    this->check_loaded_fluid();
    int ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr, 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(Tcrit);
};
CoolPropDbl REFPROPMixtureBackend::calc_p_critical() {
    this->check_loaded_fluid();
    int ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr, 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(pcrit_kPa * 1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_rhomolar_critical() {
    int ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr, 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(dcrit_mol_L * 1000);
};
void REFPROPMixtureBackend::calc_reducing_state() {
    this->check_loaded_fluid();
    double rhored_mol_L = 0, Tr = 0;
    REDXdll(&(mole_fractions[0]), &Tr, &rhored_mol_L);
    _reducing.T = Tr;
    _reducing.rhomolar = rhored_mol_L * 1000;
}
CoolPropDbl REFPROPMixtureBackend::calc_T_reducing() {
    this->check_loaded_fluid();
    double rhored_mol_L = 0, Tr = 0;
    REDXdll(&(mole_fractions[0]), &Tr, &rhored_mol_L);
    return static_cast<CoolPropDbl>(Tr);
};
CoolPropDbl REFPROPMixtureBackend::calc_rhomolar_reducing() {
    this->check_loaded_fluid();
    double rhored_mol_L = 0, Tr = 0;
    REDXdll(&(mole_fractions[0]), &Tr, &rhored_mol_L);
    return static_cast<CoolPropDbl>(rhored_mol_L * 1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_acentric_factor() {
    //     subroutine INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas) (see calc_Ttriple())
    this->check_loaded_fluid();
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    int icomp = 1L;
    // Check if more than one
    if (Ncomp == 1) {
        // Get value for first component
        INFOdll(&icomp, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        return static_cast<CoolPropDbl>(acf);
    } else {
        throw CoolProp::ValueError("acentric factor only available for pure components in REFPROP backend");
    }
};
CoolPropDbl REFPROPMixtureBackend::calc_Ttriple() {
    //     subroutine INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
    //     c
    //     c  provides fluid constants for specified component
    //     c
    //     c  input:
    //     c    icomp--component number in mixture; 1 for pure fluid
    //     c  outputs:
    //     c      wmm--molecular weight [g/mol]
    //     c     ttrp--triple point temperature [K]
    //     c    tnbpt--normal boiling point temperature [K]
    //     c       tc--critical temperature [K]
    //     c       pc--critical pressure [kPa]
    //     c       Dc--critical density [mol/L]
    //     c       Zc--compressibility at critical point [pc/(Rgas*Tc*Dc)]
    //     c      acf--acentric factor [-]
    //     c      dip--dipole moment [debye]
    //     c     Rgas--gas constant [J/mol-K]
    this->check_loaded_fluid();
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    int icomp = 1L;
    // Check if more than one
    if (Ncomp == 1) {
        // Get value for first component
        INFOdll(&icomp, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        return static_cast<CoolPropDbl>(ttrp);
    } else {
        double Tmin, Tmax, rhomolarmax, pmax;
        limits(Tmin, Tmax, rhomolarmax, pmax);
        return static_cast<CoolPropDbl>(Tmin);
    }
};
CoolPropDbl REFPROPMixtureBackend::calc_p_triple() {
    this->check_loaded_fluid();
    double p_kPa = _HUGE;
    double rho_mol_L = _HUGE, rhoLmol_L = _HUGE, rhoVmol_L = _HUGE, hmol = _HUGE, emol = _HUGE, smol = _HUGE, cvmol = _HUGE, cpmol = _HUGE, w = _HUGE;
    int ierr = 0;
    char herr[errormessagelength + 1];
    int kq = 1;
    double __T = Ttriple(), __Q = 0;
    TQFLSHdll(&__T, &__Q, &(mole_fractions[0]), &kq, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
              &(mole_fractions_vap[0]),                 // Saturation terms
              &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
              &ierr, herr, errormessagelength);         // Error terms
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    return p_kPa * 1000;
};
CoolPropDbl REFPROPMixtureBackend::calc_dipole_moment() {
    //     subroutine INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
    //     c
    //     c  provides fluid constants for specified component
    //     c
    //     c  input:
    //     c    icomp--component number in mixture; 1 for pure fluid
    //     c  outputs:
    //     c      wmm--molecular weight [g/mol]
    //     c     ttrp--triple point temperature [K]
    //     c    tnbpt--normal boiling point temperature [K]
    //     c       tc--critical temperature [K]
    //     c       pc--critical pressure [kPa]
    //     c       Dc--critical density [mol/L]
    //     c       Zc--compressibility at critical point [pc/(Rgas*Tc*Dc)]
    //     c      acf--acentric factor [-]
    //     c      dip--dipole moment [debye]
    //     c     Rgas--gas constant [J/mol-K]
    this->check_loaded_fluid();
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    int icomp = 1L;
    // Check if more than one
    if (Ncomp == 1) {
        // Get value for first component
        INFOdll(&icomp, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        return static_cast<CoolPropDbl>(dip * 3.33564e-30);
    } else {
        throw ValueError(format("dipole moment is only available for pure fluids"));
    }
};
CoolPropDbl REFPROPMixtureBackend::calc_gas_constant() {
    this->check_loaded_fluid();
    double Rmix = 0;
    RMIX2dll(&(mole_fractions[0]), &Rmix);
    return static_cast<CoolPropDbl>(Rmix);
};
CoolPropDbl REFPROPMixtureBackend::calc_molar_mass(void) {
    this->check_loaded_fluid();
    double wmm_kg_kmol;
    WMOLdll(&(mole_fractions[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
    _molar_mass = wmm_kg_kmol / 1000;             // kg/mol
    return static_cast<CoolPropDbl>(_molar_mass.pt());
};
CoolPropDbl REFPROPMixtureBackend::calc_Bvirial(void) {
    double b;
    VIRBdll(&_T, &(mole_fractions[0]), &b);
    return b * 0.001;  // 0.001 to convert from l/mol to m^3/mol
}
CoolPropDbl REFPROPMixtureBackend::calc_dBvirial_dT(void) {
    double b;
    DBDTdll(&_T, &(mole_fractions[0]), &b);
    return b * 0.001;  // 0.001 to convert from l/mol to m^3/mol
}
CoolPropDbl REFPROPMixtureBackend::calc_Cvirial(void) {
    double c;
    VIRCdll(&_T, &(mole_fractions[0]), &c);
    return c * 1e-6;  // 1e-6 to convert from (l/mol)^2 to (m^3/mol)^2
}
double REFPROPMixtureBackend::calc_melt_Tmax() {
    this->check_loaded_fluid();
    int ierr = 0;
    char herr[255];
    double tmin, tmax, Dmax_mol_L, pmax_kPa, Tmax_melt;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &tmin, &tmax, &Dmax_mol_L, &pmax_kPa, 3);
    // Get the maximum temperature for the melting curve by using the maximum pressure
    MELTPdll(&pmax_kPa, &(mole_fractions[0]), &Tmax_melt, &ierr, herr, errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return Tmax_melt;
}
CoolPropDbl REFPROPMixtureBackend::calc_melting_line(int param, int given, CoolPropDbl value) {
    this->check_loaded_fluid();
    int ierr = 0;
    char herr[255];

    if (param == iP && given == iT) {
        double _T = static_cast<double>(value), p_kPa;
        MELTTdll(&_T, &(mole_fractions[0]), &p_kPa, &ierr, herr, errormessagelength);  // Error message
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("%s", herr).c_str());
        }  //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
        return p_kPa * 1000;
    } else if (param == iT && given == iP) {
        double p_kPa = static_cast<double>(value) / 1000.0, _T;
        MELTPdll(&p_kPa, &(mole_fractions[0]), &_T, &ierr, herr, errormessagelength);  // Error message
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("%s", herr).c_str());
        }  //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
        return _T;
    } else {
        throw ValueError(format("calc_melting_line(%s,%s,%Lg) is an invalid set of inputs ", get_parameter_information(param, "short").c_str(),
                                get_parameter_information(given, "short").c_str(), value));
    }
}
bool REFPROPMixtureBackend::has_melting_line() {
    this->check_loaded_fluid();

    int ierr = 0;
    char herr[255];
    double _T = 300, p_kPa;
    MELTTdll(&_T, &(mole_fractions[0]), &p_kPa, &ierr, herr, errormessagelength);  // Error message
    if (static_cast<int>(ierr) == 1) {
        return false;
    } else {
        return true;
    }
}

const std::vector<CoolPropDbl> REFPROPMixtureBackend::calc_mass_fractions() {
    // mass fraction is mass_i/total_mass;
    // REFPROP yields mm in kg/kmol, CP uses base SI units of kg/mol;
    CoolPropDbl mm = molar_mass();
    std::vector<CoolPropDbl> mass_fractions(mole_fractions_long_double.size());
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    // FORTRAN is 1-based indexing!
    for (int i = 1L; i <= static_cast<int>(mole_fractions_long_double.size()); ++i) {
        // Get value for first component
        INFOdll(&i, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        mass_fractions[i - 1] = (wmm / 1000.0) * mole_fractions_long_double[i - 1] / mm;
    }
    return mass_fractions;
}

CoolPropDbl REFPROPMixtureBackend::calc_PIP(void) {
    // Calculate the PIP factor of Venkatharathnam and Oellrich, "Identification of the phase of a fluid using
    // partial derivatives of pressure, volume,and temperature without reference to saturation properties:
    // Applications in phase equilibria calculations"
    double t = _T, rho = _rhomolar / 1000.0,  // mol/dm^3
      p = 0, e = 0, h = 0, s = 0, cv = 0, cp = 0, w = 0, Z = 0, hjt = 0, A = 0, G = 0, xkappa = 0, beta = 0, dPdrho = 0, d2PdD2 = 0, dPT = 0,
           drhodT = 0, drhodP = 0, d2PT2 = 0, d2PdTD = 0, spare3 = 0, spare4 = 0;
    //subroutine THERM2 (t,rho,x,p,e,h,s,cv,cp,w,Z,hjt,A,G,
    // &                   xkappa,beta,dPdrho,d2PdD2,dPT,drhodT,drhodP,
    // &                   d2PT2,d2PdTD,spare3,spare4);
    THERM2dll(&t, &rho, &(mole_fractions[0]), &p, &e, &h, &s, &cv, &cp, &w, &Z, &hjt, &A, &G, &xkappa, &beta, &dPdrho, &d2PdD2, &dPT, &drhodT,
              &drhodP, &d2PT2, &d2PdTD, &spare3, &spare4);
    return 2 - rho * (d2PdTD / dPT - d2PdD2 / dPdrho);
};

CoolPropDbl REFPROPMixtureBackend::calc_viscosity(void) {
    this->check_loaded_fluid();
    double eta, tcx, rhomol_L = 0.001 * _rhomolar;
    int ierr = 0;
    char herr[255];
    TRNPRPdll(&_T, &rhomol_L, &(mole_fractions[0]),  // Inputs
              &eta, &tcx,                            // Outputs
              &ierr, herr, errormessagelength);      // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _viscosity = 1e-6 * eta;
    _conductivity = tcx;
    return static_cast<double>(_viscosity);
}
CoolPropDbl REFPROPMixtureBackend::calc_conductivity(void) {
    // Calling viscosity also caches conductivity, use that to save calls
    calc_viscosity();
    return static_cast<double>(_conductivity);
}
CoolPropDbl REFPROPMixtureBackend::calc_surface_tension(void) {
    this->check_loaded_fluid();
    double sigma, rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    char herr[255];
    SURFTdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
             &sigma,                                 // Outputs
             &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _surface_tension = sigma;
    return static_cast<double>(_surface_tension);
}
CoolPropDbl REFPROPMixtureBackend::calc_fugacity_coefficient(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> fug_cof;
    fug_cof.resize(mole_fractions.size());
    char herr[255];
    FUGCOFdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
              &(fug_cof[0]),                          // Outputs
              &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(fug_cof[i]);
}
CoolPropDbl REFPROPMixtureBackend::calc_fugacity(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> f(mole_fractions.size());
    char herr[255];
    FGCTY2dll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
              &(f[0]),                                // Outputs
              &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(f[i] * 1000);
}
CoolPropDbl REFPROPMixtureBackend::calc_chemical_potential(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> chem_pot(mole_fractions.size());
    char herr[255];
    CHEMPOTdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
               &(chem_pot[0]),                         // Outputs
               &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(chem_pot[i]);
}

void REFPROPMixtureBackend::calc_phase_envelope(const std::string& type) {
    this->check_loaded_fluid();
    double rhoymin, rhoymax, c = 0;
    int ierr = 0;
    char herr[255];
    SATSPLNdll(&(mole_fractions[0]),              // Inputs
               &ierr, herr, errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr).c_str());
    }

    // Clear the phase envelope data
    PhaseEnvelope = PhaseEnvelopeData();
    /*
    subroutine SPLNVAL (isp,iderv,a,f,ierr,herr)
    c
    c  calculates the function value of a spline
    c
    c  inputs:
    c      isp--indicator for which spline to use (1-nc: composition,
    c           nc+1: temperature, nc+2: pressure, nc+3: density,
    c           nc+4: enthalpy, nc+5: entropy)
    c    iderv--values of -1 and -2 return lower and upper root values,
    c           value of 0 returns spline function, value of 1 returns
    c           derivative of spline function with respect to density
    c        a--root value
    c  outputs:
    c        f--value of spline function at input root
    c     ierr--error flag:   0 = successful
    c     herr--error string (character*255)
    */
    int N = 500;
    int isp = 0, iderv = -1;
    if (SPLNVALdll == NULL) {
        std::string rpv = get_global_param_string("REFPROP_version");
        throw ValueError(
          format("Your version of REFFPROP(%s) does not have the SPLNVALdll function; cannot extract phase envelope values", rpv.c_str()));
    };
    SPLNVALdll(&isp, &iderv, &c, &rhoymin, &ierr, herr, errormessagelength);
    iderv = -2;
    SPLNVALdll(&isp, &iderv, &c, &rhoymax, &ierr, herr, errormessagelength);
    int nc = this->Ncomp;
    double ratio = pow(rhoymax / rhoymin, 1 / double(N));
    for (double rho_molL = rhoymin; rho_molL < rhoymax; rho_molL *= ratio) {
        double y;
        iderv = 0;

        PhaseEnvelope.x.resize(nc);
        PhaseEnvelope.y.resize(nc);
        for (isp = 1; isp <= nc; ++isp) {
            SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
            PhaseEnvelope.x[isp - 1].push_back(y);
            PhaseEnvelope.y[isp - 1].push_back(get_mole_fractions()[isp - 1]);
        }

        PhaseEnvelope.rhomolar_vap.push_back(rho_molL * 1000);
        PhaseEnvelope.lnrhomolar_vap.push_back(log(rho_molL * 1000));
        isp = nc + 1;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
        double T = y;
        PhaseEnvelope.T.push_back(y);
        PhaseEnvelope.lnT.push_back(log(y));
        isp = nc + 2;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
        PhaseEnvelope.p.push_back(y * 1000);
        PhaseEnvelope.lnp.push_back(log(y * 1000));
        isp = nc + 3;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
        PhaseEnvelope.rhomolar_liq.push_back(y * 1000);
        PhaseEnvelope.lnrhomolar_liq.push_back(log(y * 1000));
        PhaseEnvelope.Q.push_back(static_cast<double>(y > rho_molL));
        isp = nc + 4;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
        PhaseEnvelope.hmolar_vap.push_back(y);
        isp = nc + 5;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr, errormessagelength);
        PhaseEnvelope.smolar_vap.push_back(y);

        // Other outputs that could be useful
        int ierr = 0;
        char herr[255];
        double p_kPa, emol, hmol, smol, cvmol, cpmol, w, hjt, eta, tcx;
        // "Vapor"
        THERMdll(&T, &rho_molL, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
        PhaseEnvelope.cpmolar_vap.push_back(cpmol);
        PhaseEnvelope.cvmolar_vap.push_back(cvmol);
        PhaseEnvelope.speed_sound_vap.push_back(w);
        TRNPRPdll(&T, &rho_molL, &(mole_fractions[0]),  // Inputs
                  &eta, &tcx,                           // Outputs
                  &ierr, herr, errormessagelength);     // Error message
        PhaseEnvelope.viscosity_vap.push_back(eta / 1e6);
        PhaseEnvelope.conductivity_vap.push_back(tcx);
    }
}
CoolPropDbl REFPROPMixtureBackend::calc_cpmolar_idealgas(void) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    double p0, e0, h0, s0, cv0, cp0, w0, A0, G0;
    THERM0dll(&_T, &rho_mol_L, &(mole_fractions[0]), &p0, &e0, &h0, &s0, &cv0, &cp0, &w0, &A0, &G0);
    return static_cast<CoolPropDbl>(cp0);
}

phases REFPROPMixtureBackend::GetRPphase() {
    phases RPphase = iphase_unknown;
    if (ValidNumber(_Q)) {
        if ((_Q >= 0.00) && (_Q <= 1.00)) {  // CoolProp includes Q = 1 or 0 in the two phase region,
            RPphase = iphase_twophase;       // whereas RefProp designates saturated liquid and saturated vapor.
        } else if (_Q > 1.00) {              // Above saturation curve
            RPphase = iphase_gas;
            if (_T >= calc_T_critical()) {  //     ....AND T >= Tcrit
                RPphase = iphase_supercritical_gas;
            }
        } else if (_Q < 0.00) {  // Below saturation curve
            RPphase = iphase_liquid;
            if (_p >= calc_p_critical()) {  //     ....AND P >= Pcrit
                RPphase = iphase_supercritical_liquid;
            }
        } else {                       // RefProp might return Q = 920 for Metastable
            RPphase = iphase_unknown;  // but CoolProp doesn't have an enumerator for this state,
        }                              // so it's unknown as well.

        if ((_Q == 999) || (_Q == -997)) {                                // One last check for _Q == 999||-997 (Supercritical)
            RPphase = iphase_supercritical;                               // T >= Tcrit AND P >= Pcrit
            if ((std::abs(_T - calc_T_critical()) < 10 * DBL_EPSILON) &&  // IF (T == Tcrit) AND
                (std::abs(_p - calc_p_critical()) < 10 * DBL_EPSILON)) {  //    (P == Pcrit) THEN
                RPphase = iphase_critical_point;                          //    at critical point.
            };
        }
    } else {
        RPphase = iphase_unknown;
    }
    return RPphase;
}

void REFPROPMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    this->check_loaded_fluid();
    double rho_mol_L = _HUGE, rhoLmol_L = _HUGE, rhoVmol_L = _HUGE, hmol = _HUGE, emol = _HUGE, smol = _HUGE, cvmol = _HUGE, cpmol = _HUGE, w = _HUGE,
           q = _HUGE, mm = _HUGE, p_kPa = _HUGE, hjt = _HUGE;
    int ierr = 0;
    char herr[errormessagelength + 1] = " ";

    clear();

    // Check that mole fractions have been set, etc.
    check_status();

    // Get the molar mass of the fluid for the given composition
    WMOLdll(&(mole_fractions[0]), &mm);  // returns mole mass in kg/kmol
    _molar_mass = 0.001 * mm;            // [kg/mol]

    switch (input_pair) {
        case PT_INPUTS: {
            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            _T = value2;  // Want p in [kPa] in REFPROP

            if (imposed_phase_index == iphase_not_imposed) {
                // Use flash routine to find properties
                TPFLSHdll(&_T, &p_kPa, &(mole_fractions[0]), &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                          &(mole_fractions_vap[0]),                                                       // Saturation terms
                          &q, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);  //
                if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    throw ValueError(format("PT: %s", herr).c_str());
                }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            } else {
                //c  inputs:
                //c        t--temperature [K]
                //c        p--pressure [kPa]
                //c        x--composition [array of mol frac]
                //c      kph--phase flag:  1 = liquid
                //c                        2 = vapor
                //c                 N.B.:  0 = stable phase--NOT ALLOWED (use TPFLSH)
                //c                            (unless an initial guess is supplied for rho)
                //c                       -1 = force the search in the liquid phase (for metastable points)
                //c                       -2 = force the search in the vapor phase (for metastable points)
                //c   kguess--input flag:  1 = first guess for rho provided
                //c                        0 = no first guess provided
                //c      rho--first guess for molar density [mol/L], only if kguess = 1
                //c
                //c  outputs:
                //c      rho--molar density [mol/L]
                //c     ierr--error flag:  0 = successful
                //c                      200 = CRITP did not converge
                //c                      201 = illegal input (kph <= 0)
                //c                      202 = liquid-phase iteration did not converge
                //c                      203 = vapor-phase iteration did not converge
                //c     herr--error string (character*255 variable if ierr<>0)
                int kph = -10, kguess = 0;
                if (imposed_phase_index == iphase_liquid || imposed_phase_index == iphase_supercritical_liquid) {
                    kph = 1;
                } else if (imposed_phase_index == iphase_gas || imposed_phase_index == iphase_supercritical_gas) {
                    kph = 2;
                } else {
                    throw ValueError(format("PT: cannot use this imposed phase for PT inputs"));
                }
                // Calculate rho from TP
                TPRHOdll(&_T, &p_kPa, &(mole_fractions[0]), &kph, &kguess, &rho_mol_L, &ierr, herr, errormessagelength);

                // Calculate everything else
                THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
            }

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmolarT_INPUTS: {
            // Unit conversion for REFPROP
            _rhomolar = value1;
            rho_mol_L = 0.001 * value1;
            _T = value2;  // Want rho in [mol/L] in REFPROP

            if (imposed_phase_index == iphase_not_imposed || imposed_phase_index == iphase_twophase) {
                // Use flash routine to find properties
                TDFLSHdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                          &(mole_fractions_vap[0]),  // Saturation terms
                          &q, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
                if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    throw ValueError(format("DmolarT: %s", herr).c_str());
                }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            } else {
                // phase is imposed
                // Calculate everything else
                THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
            }

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassT_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmolarT_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
        case DmolarP_INPUTS: {
            // Unit conversion for REFPROP
            rho_mol_L = 0.001 * value1;
            p_kPa = 0.001 * value2;  // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
            PDFLSHdll(&p_kPa, &rho_mol_L, &(mole_fractions[0]), &_T, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),                     // Saturation terms
                      &q, &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);             // Error terms
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarP: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _rhomolar = value1;
            _p = value2;
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassP_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmolarP_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
        case DmolarHmolar_INPUTS: {
            // Unit conversion for REFPROP
            _rhomolar = value1;
            rho_mol_L = 0.001 * value1;
            hmol = value2;  // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DHFLSH (D,h,z,t,p,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
            DHFLSHdll(&rho_mol_L, &hmol, &(mole_fractions[0]), &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),  // Saturation terms
                      &q, &emol, &smol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarHmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassHmass_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(DmolarHmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case DmolarSmolar_INPUTS: {
            // Unit conversion for REFPROP
            _rhomolar = value1;
            rho_mol_L = 0.001 * value1;
            smol = value2;  // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
            DSFLSHdll(&rho_mol_L, &smol, &(mole_fractions[0]), &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),  // Saturation terms
                      &q, &emol, &hmol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarSmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassSmass_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(DmolarSmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case DmolarUmolar_INPUTS: {
            // Unit conversion for REFPROP
            _rhomolar = value1;
            rho_mol_L = 0.001 * value1;
            emol = value2;  // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            DEFLSHdll(&rho_mol_L, &emol, &(mole_fractions[0]), &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),  // Saturation terms
                      &q, &hmol, &hmol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarUmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;
            if (0) _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;         // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassUmass_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // U: [J/mol] * [kg/mol] -> [J/mol]
            update(DmolarUmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case HmolarP_INPUTS: {
            // Unit conversion for REFPROP
            hmol = value1;
            p_kPa = 0.001 * value2;  // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            PHFLSHdll(&p_kPa, &hmol, &(mole_fractions[0]), &_T, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &emol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("HmolarPmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value2;
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case HmassP_INPUTS: {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(HmolarP_INPUTS, value1 * (double)_molar_mass, value2);
            return;
        }
        case PSmolar_INPUTS: {
            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            smol = value2;  // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            PSFLSHdll(&p_kPa, &smol, &(mole_fractions[0]), &_T, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &emol, &hmol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("PSmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case PSmass_INPUTS: {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(PSmolar_INPUTS, value1, value2 * (double)_molar_mass);
            return;
        }
        case PUmolar_INPUTS: {
            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            emol = value2;  // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine PEFLSH (p,e,z,t,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            PEFLSHdll(&p_kPa, &emol, &(mole_fractions[0]), &_T, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("PUmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case PUmass_INPUTS: {
            // Call again, but this time with molar units
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(PUmolar_INPUTS, value1, value2 * (double)_molar_mass);
            return;
        }
        case HmolarSmolar_INPUTS: {
            // Unit conversion for REFPROP
            hmol = value1;
            smol = value2;

            HSFLSHdll(&hmol, &smol, &(mole_fractions[0]), &_T, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),          // Saturation terms
                      &q, &emol, &cvmol, &cpmol, &w,     // Other thermodynamic terms
                      &ierr, herr, errormessagelength);  // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("HmolarSmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case HmassSmass_INPUTS: {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol/K]
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(HmolarSmolar_INPUTS, value1 * (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case SmolarUmolar_INPUTS: {
            // Unit conversion for REFPROP
            smol = value1;
            emol = value2;

            // from REFPROP: subroutine ESFLSH (e,s,z,t,p,D,Dl,Dv,x,y,q,h,cv,cp,w,ierr,herr)
            ESFLSHdll(&emol, &smol, &(mole_fractions[0]), &_T, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),          // Saturation terms
                      &q, &smol, &cvmol, &cpmol, &w,     // Other thermodynamic terms
                      &ierr, herr, errormessagelength);  // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("SmolarUmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case SmassUmass_INPUTS: {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K],
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(SmolarUmolar_INPUTS, value1 * (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case SmolarT_INPUTS: {
            // Unit conversion for REFPROP
            smol = value1;
            _T = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            int kr = 1;

            // from REFPROP: subroutine TSFLSH (t,s,z,kr,p,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
            TSFLSHdll(&_T, &smol, &(mole_fractions[0]), &kr, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &emol, &hmol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("SmolarT: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case SmassT_INPUTS: {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(SmolarT_INPUTS, value1 * (double)_molar_mass, value2);
            return;
        }
        case HmolarT_INPUTS: {
            // Unit conversion for REFPROP
            hmol = value1;
            _T = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            int kr = 1;

            // from REFPROP: subroutine THFLSH (t,h,z,kr,p,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
            THFLSHdll(&_T, &hmol, &(mole_fractions[0]), &kr, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &emol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("HmolarT: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case HmassT_INPUTS: {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(HmolarT_INPUTS, value1 * (double)_molar_mass, value2);
            return;
        }
        case TUmolar_INPUTS: {
            // Unit conversion for REFPROP
            _T = value1;
            emol = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            int kr = 1;

            // from REFPROP: subroutine TEFLSH (t,e,z,kr,p,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            TEFLSHdll(&_T, &emol, &(mole_fractions[0]), &kr, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),              // Saturation terms
                      &q, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);      // Error terms

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("TUmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;                     // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;          // 1000 for conversion from mol/L to mol/m3
            if (0) _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;         // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case TUmass_INPUTS: {
            // Call again, but this time with molar units
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(TUmolar_INPUTS, value1, value2 * (double)_molar_mass);
            return;
        }
        case PQ_INPUTS: {

            // c  Estimate temperature, pressure, and compositions to be used
            // c  as initial guesses to SATTP
            // c
            // c  inputs:
            // c   iFlsh--Phase flag:    0 - Flash calculation (T and P known)
            // c                         1 - T and xliq known, P and xvap returned
            // c                         2 - T and xvap known, P and xliq returned
            // c                         3 - P and xliq known, T and xvap returned
            // c                         4 - P and xvap known, T and xliq returned
            // c                         if this value is negative, the retrograde point will be returned
            // c        t--temperature [K] (input or output)
            // c        p--pressure [MPa] (input or output)
            // c        x--composition [array of mol frac]
            // c   iGuess--if set to 1, all inputs are used as initial guesses for the calculation
            // c  outputs:
            // c        d--overall molar density [mol/L]
            // c       Dl--molar density [mol/L] of saturated liquid
            // c       Dv--molar density [mol/L] of saturated vapor
            // c     xliq--liquid phase composition [array of mol frac]
            // c     xvap--vapor phase composition [array of mol frac]
            // c        q--quality
            // c     ierr--error flag:   0 = successful
            // c                         1 = unsuccessful
            // c     herr--error string (character*255 variable if ierr<>0)

            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            q = value2;  // Want p in [kPa] in REFPROP

            int iFlsh = 0, iGuess = 0, ierr = 0;
            if (std::abs(q) < 1e-10) {
                iFlsh = 3;  // bubble point
            } else if (std::abs(q - 1) < 1e-10) {
                iFlsh = 4;  // dew point
            }
            if (iFlsh != 0) {
                // SATTP (t,p,x,iFlsh,iGuess,d,Dl,Dv,xliq,xvap,q,ierr,herr)
                SATTPdll(&_T, &p_kPa, &(mole_fractions[0]), &iFlsh, &iGuess, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                         &(mole_fractions_vap[0]), &q, &ierr, herr, errormessagelength);
                if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    ierr = 0;
                    // SATPdll(p, z, kph, T, Dl, Dv, x, y, ierr, herr)
                    //
                    //kph--Phase flag : 1 - Input x is liquid composition(bubble point)
                    //                - 1 - Force calculation in the liquid phase even if T<Ttrp
                    //                  2 - Input x is vapor composition(dew point)
                    //                - 2 - Force calculation in the vapor phase even if T<Ttrp
                    //                  3 - Input x is liquid composition along the freezing line(melting line)
                    //                  4 - Input x is vapor composition along the sublimation line
                    SATPdll(&_p, &(mole_fractions[0]), &iFlsh, &_T, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]), &(mole_fractions_vap[0]), &ierr,
                            herr, errormessagelength);
                    rho_mol_L = (iFlsh == 1) ? rhoLmol_L : rhoVmol_L;
                }
                if (ierr <= 0L) {
                    // Calculate everything else
                    THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
                }
            }
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD) || iFlsh == 0) {
                // From REFPROP:
                //additional input--only for TQFLSH and PQFLSH
                //     kq--flag specifying units for input quality
                //         kq = 1 quality on MOLAR basis [moles vapor/total moles]
                //         kq = 2 quality on MASS basis [mass vapor/total mass]
                int kq = 1;
                ierr = 0;
                // Use flash routine to find properties
                PQFLSHdll(&p_kPa, &q, &(mole_fractions[0]), &kq, &_T, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                          &(mole_fractions_vap[0]),                 // Saturation terms
                          &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                          &ierr, herr, errormessagelength);         // Error terms
            }

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("PQ: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case QT_INPUTS: {
            // Unit conversion for REFPROP
            q = value1;
            _T = value2;

            // Use flash routine to find properties
            int iFlsh = 0, iGuess = 0;
            if (std::abs(q) < 1e-10) {
                iFlsh = 1;  // bubble point with T given
            } else if (std::abs(q - 1) < 1e-10) {
                iFlsh = 2;  // dew point with T given
            }
            if (iFlsh != 0) {
                // SATTP (t,p,x,iFlsh,iGuess,d,Dl,Dv,xliq,xvap,q,ierr,herr)
                SATTPdll(&_T, &p_kPa, &(mole_fractions[0]), &iFlsh, &iGuess, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                         &(mole_fractions_vap[0]), &q, &ierr, herr, errormessagelength);

                if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    ierr = 0;
                    // SATTdll(T, z, kph, P, Dl, Dv, x, y, ierr, herr)
                    //
                    //kph--Phase flag : 1 - Input x is liquid composition(bubble point)
                    //                - 1 - Force calculation in the liquid phase even if T<Ttrp
                    //                  2 - Input x is vapor composition(dew point)
                    //                - 2 - Force calculation in the vapor phase even if T<Ttrp
                    //                  3 - Input x is liquid composition along the freezing line(melting line)
                    //                  4 - Input x is vapor composition along the sublimation line
                    SATTdll(&_T, &(mole_fractions[0]), &iFlsh, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]), &(mole_fractions_vap[0]),
                            &ierr, herr, errormessagelength);
                    rho_mol_L = (iFlsh == 1) ? rhoLmol_L : rhoVmol_L;
                }
                if (ierr <= 0L) {
                    // Calculate everything else since we were able to carry out a flash call
                    THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
                }
            }
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD) || iFlsh == 0) {
                ierr = 0;
                /* From REFPROP:
                additional input--only for TQFLSH and PQFLSH
                kq--flag specifying units for input quality
                kq = 1 quality on MOLAR basis [moles vapor/total moles]
                kq = 2 quality on MASS basis [mass vapor/total mass]
                */
                int kq = 1;
                TQFLSHdll(&_T, &q, &(mole_fractions[0]), &kq, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                          &(mole_fractions_vap[0]),                 // Saturation terms
                          &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                          &ierr, herr, errormessagelength);         // Error terms
            }

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("TQ(%s): %s", LoadedREFPROPRef.c_str(), herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        default: {
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
        }
    };
    // Set these common variables that are used in every flash calculation
    _hmolar = hmol;
    _smolar = smol;
    _umolar = emol;
    _cvmolar = cvmol;
    _cpmolar = cpmol;
    _speed_sound = w;
    _tau = calc_T_reducing() / _T;
    _delta = _rhomolar / calc_rhomolar_reducing();
    _gibbsmolar = hmol - _T * smol;
    _Q = q;
    if (imposed_phase_index == iphase_not_imposed) {  // If phase is imposed, _phase will already be set.
        if (Ncomp == 1) {                             // Only set _phase for pure fluids
            _phase = GetRPphase();                    // Set the CoolProp _phase variable based on RefProp's quality value (q)
        }
    }
}

void REFPROPMixtureBackend::update_with_guesses(CoolProp::input_pairs input_pair, double value1, double value2, const GuessesStructure& guesses) {
    this->check_loaded_fluid();
    double rho_mol_L = _HUGE, hmol = _HUGE, emol = _HUGE, smol = _HUGE, cvmol = _HUGE, cpmol = _HUGE, w = _HUGE, q = _HUGE, p_kPa = _HUGE,
           hjt = _HUGE;
    int ierr = 0;
    char herr[errormessagelength + 1];

    clear();

    // Check that mole fractions have been set, etc.
    check_status();

    switch (input_pair) {
        case PT_INPUTS: {
            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            _T = value2;  // Want p in [kPa] in REFPROP

            //THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
            int kguess = 1;  // guess provided
            if (!ValidNumber(guesses.rhomolar)) {
                throw ValueError(format("rhomolar must be provided in guesses"));
            }

            int kph = (guesses.rhomolar > calc_rhomolar_critical())
                        ? 1
                        : 2;  // liquid  if density > rhoc, vapor otherwise - though we are passing the guess density
            rho_mol_L = guesses.rhomolar / 1000.0;
            TPRHOdll(&_T, &p_kPa, &(mole_fractions[0]), &kph, &kguess, &rho_mol_L, &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("PT: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case PQ_INPUTS: {
            // Unit conversion for REFPROP
            p_kPa = 0.001 * value1;
            q = value2;  // Want p in [kPa] in REFPROP
            double rhoLmol_L = -1, rhoVmol_L = -1;
            int iFlsh = 0,
                iGuess = 1,  // Use guesses
              ierr = 0;

            if (std::abs(value2) < 1e-10) {
                iFlsh = 3;  // bubble point
                if (!guesses.x.empty()) {
                    mole_fractions = guesses.x;
                    while (mole_fractions.size() < ncmax) {
                        mole_fractions.push_back(0.0);
                    }
                } else {
                    throw ValueError(format("x must be provided in guesses"));
                }
            } else if (std::abs(value2 - 1) < 1e-10) {
                iFlsh = 4;  // dew point
                if (!guesses.y.empty()) {
                    mole_fractions = guesses.y;
                    while (mole_fractions.size() < ncmax) {
                        mole_fractions.push_back(0.0);
                    }
                } else {
                    throw ValueError(format("y must be provided in guesses"));
                }
            } else {
                throw ValueError(format("For PQ w/ guesses, Q must be either 0 or 1"));
            }
            if (get_debug_level() > 9) {
                std::cout << format("guesses.T: %g\n", guesses.T);
            }
            if (!ValidNumber(guesses.T)) {
                throw ValueError(format("T must be provided in guesses"));
            } else {
                _T = guesses.T;
            }

            // SATTP (t,p,x,iFlsh,iGuess,d,Dl,Dv,xliq,xvap,q,ierr,herr)
            SATTPdll(&_T, &p_kPa, &(mole_fractions[0]), &iFlsh, &iGuess, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                     &(mole_fractions_vap[0]), &q, &ierr, herr, errormessagelength);

            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("PQ: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        default: {
            throw CoolProp::ValueError(format("Unable to match given input_pair in update_with_guesses"));
        }
    }

    // Calculate everything else
    THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);

    // Set these common variables that are used in every flash calculation
    _hmolar = hmol;
    _smolar = smol;
    _umolar = emol;
    _cvmolar = cvmol;
    _cpmolar = cpmol;
    _speed_sound = w;
    _tau = calc_T_reducing() / _T;
    _delta = _rhomolar / calc_rhomolar_reducing();
    _Q = q;
}
CoolPropDbl REFPROPMixtureBackend::call_phixdll(int itau, int idel) {
    this->check_loaded_fluid();
    double val = 0, tau = _tau, delta = _delta;
    if (PHIXdll == NULL) {
        throw ValueError("PHIXdll function is not available in your version of REFPROP. Please upgrade");
    }
    PHIXdll(&itau, &idel, &tau, &delta, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val) / pow(static_cast<CoolPropDbl>(_delta), idel) / pow(static_cast<CoolPropDbl>(_tau), itau);
}
CoolPropDbl REFPROPMixtureBackend::call_phi0dll(int itau, int idel) {
    this->check_loaded_fluid();
    double val = 0, tau = _tau, __T = T(), __rho = rhomolar() / 1000;
    if (PHI0dll == NULL) {
        throw ValueError("PHI0dll function is not available in your version of REFPROP. Please upgrade");
    }
    PHI0dll(&itau, &idel, &__T, &__rho, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val) / pow(static_cast<CoolPropDbl>(_delta), idel) / pow(tau, itau);
}
/// Calculate excess properties
void REFPROPMixtureBackend::calc_excess_properties() {
    this->check_loaded_fluid();
    int ierr = 0;
    char herr[255];
    double T_K = _T, p_kPa = _p / 1000.0, rho = 1, vE = -1, eE = -1, hE = -1, sE = -1, aE = -1, gE = -1;
    int kph = 1;

    //    subroutine EXCESS(t, p, x, kph, rho, vE, eE, hE, sE, aE, gE, ierr, herr)
    //    c
    //    c  compute excess properties as a function of temperature, pressure,
    //    c  and composition.
    //    c
    //    c  inputs :
    //c        t--temperature[K]
    //    c        p--pressure[kPa]
    //    c        x--composition[array of mol frac]
    //    c      kph--phase flag : 1 = liquid
    //    c                        2 = vapor
    //    c                        0 = stable phase
    //    c  outputs :
    //c       rho--molar density[mol / L](if input less than 0, used as initial guess)
    //    c        vE--excess volume[L / mol]
    //    c        eE--excess energy[J / mol]
    //    c        hE--excess enthalpy[J / mol]
    //    c        sE--excess entropy[J / mol - K]
    //    c        aE--excess Helmholtz energy[J / mol]
    //    c        gE--excess Gibbs energy[J / mol]
    //    c      ierr--error flag : 0 = successful
    //    c                        55 = T, p inputs in different phase for the pure fluids
    //    c      herr--error string(character * 255 variable if ierr<>0)
    EXCESSdll(&T_K, &p_kPa, &(mole_fractions[0]), &kph, &rho, &vE, &eE, &hE, &sE, &aE, &gE, &ierr, herr, errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("EXCESSdll: %s", herr).c_str());
    }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _volumemolar_excess = vE;
    _umolar_excess = eE;
    _hmolar_excess = hE;
    _smolar_excess = sE;
    _helmholtzmolar_excess = aE;
    _gibbsmolar_excess = gE;
}

void REFPROPMixtureBackend::calc_true_critical_point(double& T, double& rho) {

    class wrapper : public FuncWrapperND
    {
       public:
        const std::vector<double> z;
        wrapper(const std::vector<double>& z) : z(z){};
        std::vector<double> call(const std::vector<double>& x) {
            std::vector<double> r(2);
            double dpdrho__constT = _HUGE, d2pdrho2__constT = _HUGE;
            DPDDdll(const_cast<double*>(&(x[0])), const_cast<double*>(&(x[1])), const_cast<double*>(&(z[0])), &dpdrho__constT);
            DPDD2dll(const_cast<double*>(&(x[0])), const_cast<double*>(&(x[1])), const_cast<double*>(&(z[0])), &d2pdrho2__constT);
            r[0] = dpdrho__constT;
            r[1] = d2pdrho2__constT;
            return r;
        };
    };
    wrapper resid(mole_fractions);

    T = calc_T_critical();
    double rho_moldm3 = calc_rhomolar_critical() / 1000.0;
    std::vector<double> x(2, T);
    x[1] = rho_moldm3;
    std::vector<double> xfinal = NDNewtonRaphson_Jacobian(&resid, x, 1e-9, 30);
    T = xfinal[0];
    rho = xfinal[1] * 1000.0;
}

CoolPropDbl REFPROPMixtureBackend::calc_saturated_liquid_keyed_output(parameters key) {
    if (_rhoLmolar) {
        if (key == iDmolar) {
            return _rhoLmolar;
        } else if (key == iDmass) {
            return static_cast<double>(_rhoLmolar) * calc_saturated_liquid_keyed_output(imolar_mass);
        } else if (key == imolar_mass) {
            double wmm_kg_kmol = 0;
            WMOLdll(&(mole_fractions_liq[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
            return wmm_kg_kmol / 1000;                        // kg/mol
        } else {
            throw ValueError("Invalid parameter. Only mass and molar density are available with RefProp");
            return _HUGE;
        }
    }
    throw ValueError("The saturated liquid state has not been set.");
    return _HUGE;
}
CoolPropDbl REFPROPMixtureBackend::calc_saturated_vapor_keyed_output(parameters key) {
    if (_rhoVmolar) {
        if (key == iDmolar) {
            return _rhoVmolar;
        } else if (key == iDmass) {
            return static_cast<double>(_rhoVmolar) * calc_saturated_vapor_keyed_output(imolar_mass);
        } else if (key == imolar_mass) {
            double wmm_kg_kmol = 0;
            WMOLdll(&(mole_fractions_vap[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
            return wmm_kg_kmol / 1000;                        // kg/mol
        } else {
            throw ValueError("Invalid key.");
            return _HUGE;
        }
    }
    throw ValueError("The saturated vapor state has not been set.");
    return _HUGE;
}

void REFPROPMixtureBackend::calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) {
    if (type == "Joule-Thomson") {
        JouleThomsonCurveTracer JTCT(this, 1e5, 800);
        JTCT.trace(T, p);
    } else if (type == "Joule-Inversion") {
        JouleInversionCurveTracer JICT(this, 1e5, 800);
        JICT.trace(T, p);
    } else if (type == "Ideal") {
        IdealCurveTracer ICT(this, 1e5, 800);
        ICT.trace(T, p);
    } else if (type == "Boyle") {
        BoyleCurveTracer BCT(this, 1e5, 800);
        BCT.trace(T, p);
    } else {
        throw ValueError(format("Invalid ideal curve type: %s", type.c_str()));
    }
};

bool force_load_REFPROP() {
    std::string err;
    if (!load_REFPROP(err)) {
        if (CoolProp::get_debug_level() > 5) {
            std::cout << format("Error while loading REFPROP: %s", err) << std::endl;
        }
        LoadedREFPROPRef = "";
        return false;
    } else {
        LoadedREFPROPRef = "";
        return true;
    }
}
bool force_unload_REFPROP() {
    std::string err;
    if (!unload_REFPROP(err)) {
        if (CoolProp::get_debug_level() > 5) {
            std::cout << format("Error while unloading REFPROP: %s", err) << std::endl;
        }
        LoadedREFPROPRef = "";
        return false;
    } else {
        LoadedREFPROPRef = "";
        return true;
    }
}

void REFPROP_SETREF(char hrf[3], int ixflag, double x0[1], double& h0, double& s0, double& T0, double& p0, int& ierr, char herr[255], int l1,
                    int l2) {
    std::string err;
    bool loaded_REFPROP = ::load_REFPROP(err);
    if (!loaded_REFPROP) {
        throw ValueError(format("Not able to load REFPROP; err is: %s", err.c_str()));
    }
    SETREFdll(hrf, &ixflag, x0, &h0, &s0, &T0, &p0, &ierr, herr, l1, l2);
}

} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include "CoolProp.h"
#    include <catch2/catch_all.hpp>

TEST_CASE("Check REFPROP and CoolProp values agree", "[REFPROP]") {
    SECTION("Saturation densities agree within 0.5% at T/Tc = 0.9") {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"), ',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it) {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it), "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            CHECK_NOTHROW(S1->update(CoolProp::QT_INPUTS, 0, Tr * 0.9));
            double rho_CP = S1->rhomolar();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            CHECK_NOTHROW(S2->update(CoolProp::QT_INPUTS, 0, Tr * 0.9));
            double rho_RP = S2->rhomolar();

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(rho_CP);
            CAPTURE(rho_RP);

            double DH = (rho_RP - rho_CP) / rho_RP;
            CHECK(std::abs(DH) < 0.05);
        }
    }
    SECTION("Saturation specific heats agree within 0.5% at T/Tc = 0.9") {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"), ',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it) {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it), "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            S1->update(CoolProp::QT_INPUTS, 0, Tr * 0.9);
            double cp_CP = S1->cpmolar();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            S2->update(CoolProp::QT_INPUTS, 0, Tr * 0.9);
            double cp_RP = S2->cpmolar();

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(cp_CP);
            CAPTURE(cp_RP);
            CAPTURE(0.9 * Tr);

            double Dcp = (cp_RP - cp_CP) / cp_RP;
            CHECK(std::abs(Dcp) < 0.05);
        }
    }
    SECTION("Enthalpy and entropy reference state") {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"), ',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it) {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it), "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            CHECK_NOTHROW(S1->update(CoolProp::QT_INPUTS, 0, 0.9 * Tr));
            double h_CP = S1->hmass();
            double s_CP = S1->smass();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            CHECK_NOTHROW(S2->update(CoolProp::QT_INPUTS, 0, 0.9 * Tr));
            double h_RP = S2->hmass();
            double s_RP = S2->smass();

            double delta_a1 = (s_CP - s_RP) / (S1->gas_constant() / S1->molar_mass());
            double delta_a2 = -(h_CP - h_RP) / (S1->gas_constant() / S1->molar_mass() * S1->get_reducing_state().T);
            CAPTURE(format("%0.16f", delta_a1));
            CAPTURE(format("%0.16f", delta_a2));

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(h_CP);
            CAPTURE(h_RP);
            CAPTURE(s_CP);
            CAPTURE(s_RP);
            double DH = (S1->hmass() - S2->hmass());
            double DS = (S1->smass() - S2->smass());

            CHECK(std::abs(DH / h_RP) < 0.01);
            CHECK(std::abs(DS / s_RP) < 0.01);
        }
    }
}

TEST_CASE("Check trivial inputs for REFPROP work", "[REFPROP_trivial]") {
    const int num_inputs = 6;
    std::string inputs[num_inputs] = {"T_triple", "T_critical", "p_critical", "molar_mass", "rhomolar_critical", "rhomass_critical"};
    for (int i = 0; i < num_inputs; ++i) {
        std::ostringstream ss;
        ss << "Check " << inputs[i];
        SECTION(ss.str(), "") {
            double cp_val = CoolProp::PropsSI(inputs[i], "P", 0, "T", 0, "HEOS::Water");
            double rp_val = CoolProp::PropsSI(inputs[i], "P", 0, "T", 0, "REFPROP::Water");

            std::string errstr = CoolProp::get_global_param_string("errstring");
            CAPTURE(errstr);
            double err = (cp_val - rp_val) / cp_val;
            CHECK(err < 1e-3);
        }
    }
}

TEST_CASE("Check PHI0 derivatives", "[REFPROP_PHI0]") {

    const int num_inputs = 3;
    std::string inputs[num_inputs] = {"DALPHA0_DDELTA_CONSTTAU", "D2ALPHA0_DDELTA2_CONSTTAU", "D3ALPHA0_DDELTA3_CONSTTAU"};
    for (int i = 0; i < num_inputs; ++i) {
        std::ostringstream ss;
        ss << "Check " << inputs[i];
        SECTION(ss.str(), "") {
            double cp_val = CoolProp::PropsSI(inputs[i], "P", 101325, "T", 298, "HEOS::Water");
            double rp_val = CoolProp::PropsSI(inputs[i], "P", 101325, "T", 298, "REFPROP::Water");

            std::string errstr = CoolProp::get_global_param_string("errstring");
            CAPTURE(errstr);
            double err = std::abs((cp_val - rp_val) / cp_val);
            CHECK(err < 1e-12);
        }
    }
}

#endif
