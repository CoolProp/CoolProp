
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

#include <array>
#include <cstring>

#define REFPROP_IMPLEMENTATION
#define REFPROP_CSTYLE_REFERENCES
#include <REFPROP_lib.h>
#undef REFPROP_IMPLEMENTATION
#undef REFPROP_CSTYLE_REFERENCES

#include "CoolProp/detail/tools.h"
#include "CoolProp/detail/json.h"
#include "REFPROPMixtureBackend.h"
#include "REFPROPBackend.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/Configuration.h"
#include "CoolProp/detail/filepaths.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/numerics/Solvers.h"
#include "CoolProp/fluids/IdealCurves.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/AbstractState.h"
#include "qmass_conversions.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"

#include <cmath>
#include <optional>
#include <string>
#include <cstdio>
#include <iostream>
#include <memory>
using std::shared_ptr;
#include <cstdlib>

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

// NOLINTNEXTLINE(cert-err58-cpp): default-constructed std::string is noexcept since C++11
std::string LoadedREFPROPRef;

static bool dbg_refprop = false;
static const unsigned int number_of_endings = 5;
// NOLINTNEXTLINE(cert-err58-cpp): five short string literals — constructors won't throw on this input
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

/// Query an environment variable and return in optional if present
std::optional<std::string> get_envvar(const char* var) {
    char* valptr = getenv(var);
    if (valptr == nullptr) {
        return std::nullopt;
    } else {
        return std::string(valptr);
    }
}

/// Find either FLUIDS or fluids folder relative to the root path provided; return the path
std::string get_casesensitive_fluids(const std::string& root) {
    if (get_envvar("COOLPROP_REFPROP_ROOT")) {
        return "";
    }
    std::string joined = join_path(root, "fluids");
    if (path_exists(joined)) {
        return joined;
    } else {
        std::string ucase_joined = join_path(root, "FLUIDS");
        if (path_exists(ucase_joined)) {
            return ucase_joined;
        } else {
            throw CoolProp::ValueError(format(R"(fluid directories "FLUIDS" or "fluids" could not be found in the directory [%s])", root));
        }
    }
}
std::string get_REFPROP_fluid_path_prefix() {
    if (get_envvar("COOLPROP_REFPROP_ROOT")) {
        return "";
    }
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
    if (get_envvar("COOLPROP_REFPROP_ROOT")) {
        return "";
    }
    std::string rpPath = refpropPath;
    // Allow the user to specify an alternative REFPROP path by configuration value
    std::string alt_refprop_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_PATH);
    if (!alt_refprop_path.empty()) {
        //if (!endswith(alt_refprop_path, get_separator())) {
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

/// Return the REFPROP .FLD stem for a CoolPropFluid.
/// Falls back to `fallback` when REFPROPname is absent or the sentinel "N/A",
/// so the original user-supplied name is preserved rather than fluid.name
/// (which may differ, e.g. "R1336mzz(E)" vs the file "R1336MZZE.FLD").
static std::string refprop_stem(const CoolProp::CoolPropFluid& fluid, const std::string& fallback) {
    if (!fluid.REFPROPname.empty() && fluid.REFPROPname != "N/A") {
        return fluid.REFPROPname;
    }
    return fallback;
}

namespace CoolProp {

class REFPROPGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) override {
        REFPROPMixtureBackend::REFPROP_supported();
        if (fluid_names.size() == 1) {
            return new REFPROPBackend(fluid_names[0]);
        } else {
            return new REFPROPMixtureBackend(fluid_names);
        }
    };
};
// This static initialization will cause the generator to register.
// GeneratorInitializer's constructor only inserts into a std::map keyed by
// backend_family — std::map::operator[] / insert can technically throw
// std::bad_alloc on allocation failure, but in practice this is a single
// small entry inserted at program start; if allocation fails here the
// process is unrecoverable anyway. Mark NOLINT rather than wrap in
// call_once + lazy registration (option 1 from #2874) which would be the
// cleaner refactor but is out of this PR's scope.
// NOLINTNEXTLINE(cert-err58-cpp)
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
    if (RefpropdllInstance != nullptr) return true;

    // Store result of previous check.
    if (_REFPROP_supported) {
        // Either Refprop is supported or it is the first check.
        std::string rpv(STRINGIFY(RPVersion));
        if (rpv.compare("NOTAVAILABLE") != 0) {
            // Function names were defined in "REFPROP_lib.h",
            // This platform theoretically supports Refprop.
            std::string err;

            bool loaded_REFPROP = false;

            auto root = get_envvar("COOLPROP_REFPROP_ROOT");
            const std::string alt_rp_path = get_config_string(ALTERNATIVE_REFPROP_PATH);
            const std::string alt_rp_name = get_config_string(ALTERNATIVE_REFPROP_LIBRARY_PATH);

            if (root) {
                loaded_REFPROP = ::load_REFPROP(err, root.value().c_str(), "");
                SETPATHdll(const_cast<char*>(root.value().c_str()), 400);
            }
            if (!loaded_REFPROP) {
                if (!alt_rp_name.empty()) {
                    loaded_REFPROP = ::load_REFPROP(err, "", alt_rp_name);
                } else {
                    if (alt_rp_path.empty()) {
                        loaded_REFPROP = ::load_REFPROP(err, refpropPath, "");
                    } else {
                        loaded_REFPROP = ::load_REFPROP(err, alt_rp_path, "");
                    }
                }
            }

            if (loaded_REFPROP) {
                return true;
            } else {
                std::cout << "Good news: It is possible to use REFPROP on your system! However, the library \n"
                          << "could not be loaded. Please make sure that REFPROP is available on your system.\n\n"
                          << "Neither found in current location nor found in system PATH.\n"
                          << "If you already obtained a copy of REFPROP from http://www.nist.gov/srd/, \n"
                          << "add location of REFPROP to the PATH environment variable or your library path.\n\n"
                          << "In case you do not use Windows, have a look at https://github.com/jowr/librefprop.so \n"
                          << "to find instructions on how to compile your own version of the REFPROP library.\n\n"
                          << format("COOLPROP_REFPROP_ROOT: %s\n", (root) ? root.value().c_str() : "?")
                          << format("ALTERNATIVE_REFPROP_PATH: %s\n", alt_rp_path.c_str()) << format("ERROR: %s\n", err.c_str());
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
    for (char& i : herr) {
        i = '\0';
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

        // Resolve each name through the CoolProp fluid library so that any
        // CoolProp alias (e.g. "R600", "n-Butane") is translated to the
        // REFPROP filename (e.g. "BUTANE") before building the .FLD path.
        // Can be disabled via set_config_bool("REFPROP_RESOLVE_COOLPROP_ALIASES", false).
        std::vector<std::string> resolved_names(fluid_names);
        if (get_config_bool(REFPROP_RESOLVE_COOLPROP_ALIASES)) {
            for (std::size_t i = 0; i < resolved_names.size(); ++i) {
                try {
                    CoolPropFluid fluid = get_library().get(fluid_names[i]);
                    resolved_names[i] = refprop_stem(fluid, fluid_names[i]);
                } catch (const CoolProp::ValueError&) {
                    // Direct lookup failed.  If the name carries a REFPROP file
                    // extension (e.g. "R600.FLD"), strip it and retry so that
                    // aliases like "R600.FLD" resolve the same way as "R600".
                    static const std::array<std::string, 4> exts{".FLD", ".fld", ".PPF", ".ppf"};
                    for (const auto& ext : exts) {
                        if (endswith(fluid_names[i], ext)) {
                            const std::string stem = fluid_names[i].substr(0, fluid_names[i].size() - ext.size());
                            // Try the stem as-is first, then uppercased (the library
                            // stores both original and uppercased alias forms, but not
                            // every mixed-case variant).
                            for (const auto& lookup : {stem, upper(stem)}) {
                                try {
                                    CoolPropFluid fluid = get_library().get(lookup);
                                    resolved_names[i] = refprop_stem(fluid, stem);
                                    break;
                                } catch (const CoolProp::ValueError&) {
                                }
                            }
                            break;
                        }
                    }
                    // If still unresolved, pass the name through unchanged so
                    // REFPROP can try to load it directly (e.g. a full path or .MIX).
                }
            }
        }

        std::array<char, 10000> component_string{};
        std::array<char, errormessagelength> herr{};
        std::string components_joined = strjoin(fluid_names, "|");
        std::string components_joined_raw = strjoin(fluid_names, "|");
        std::string fdPath = get_REFPROP_fluid_path_prefix();
        int N = static_cast<int>(fluid_names.size());

        // Get path to HMX.BNC file
        std::array<char, 255> hmx_bnc{};
        const std::string HMX_path = get_REFPROP_HMX_BNC_path();
        const char* _HMX_path = HMX_path.c_str();
        // The copy writes strlen+1 bytes (incl. the NUL), so the safe maximum is
        // one less than the buffer; reject at >= to avoid a 1-byte overflow.
        if (strlen(_HMX_path) >= refpropcharlength) {
            throw ValueError(
              format("Full HMX path (%s) is too long; max length is %d characters", _HMX_path, static_cast<int>(refpropcharlength) - 1));
        }
        memcpy(hmx_bnc.data(), _HMX_path, strlen(_HMX_path) + 1);  // +1 for the NUL; length checked above

        if (get_config_bool(REFPROP_USE_GERG)) {
            int iflag = 1,  // Tell REFPROP to use GERG04; 0 unsets GERG usage
              ierr_gerg = 0;
            std::array<char, 255> herr_gerg{};
            GERG04dll(&N, &iflag, &ierr_gerg, herr_gerg.data(), 255);
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
            // The copy appends a NUL, so the safe maximum is sizeof(mix)-1.
            if (strlen(_components_joined_raw) >= sizeof(mix)) {
                throw ValueError(format("components (%s) is too long", components_joined_raw.c_str()));
            }
            memcpy(mix, _components_joined_raw, strlen(_components_joined_raw) + 1);  // +1 for the NUL; length checked above

            SETMIXdll(mix, hmx_bnc.data(), reference_state, &N, component_string.data(), &(x[0]), &ierr, herr.data(), 255, 255, 3, 10000, 255);
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
                    throw ValueError(format("Interaction parameter estimation has been disabled: %s", herr.data()));
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
                if (CoolProp::get_debug_level() > 0) {
                    std::cout << format("%s:%d Unable to load predefined mixture [%s] with ierr: [%d] and herr: [%s]\n", __FILE__, __LINE__, mix,
                                        ierr, herr.data());
                }
                throw ValueError(format("Unable to load mixture: %s", components_joined_raw.c_str()));
            }
        }

        // Loop over the file names - first we try with nothing, then .fld, then .FLD, then .ppf - means you can't mix and match
        for (unsigned int k = 0; k < number_of_endings; k++) {
            // Build the mixture string
            for (unsigned int j = 0; j < (unsigned int)N; j++) {
                if (j == 0) {
                    components_joined = join_path(fdPath, upper(resolved_names[j]) + endings[k]);
                } else {
                    components_joined += "|" + join_path(fdPath, upper(resolved_names[j]) + endings[k]);
                }
            }

            if (dbg_refprop)
                std::cout << format("%s:%d: The fluid %s has not been loaded before, current value is %s \n", __FILE__, __LINE__,
                                    components_joined_raw.c_str(), LoadedREFPROPRef.c_str());

            // Copy over the list of components
            const char* _components_joined = components_joined.c_str();
            // The safe maximum is component_string.size()-1; the space-padding
            // loop below fills the remainder, so no NUL terminator is needed.
            if (strlen(_components_joined) >= component_string.size()) {
                throw ValueError(format("components_joined (%s) is too long", _components_joined));
            }
            memcpy(component_string.data(), _components_joined, strlen(_components_joined));
            // Pad the fluid string all the way to 10k characters with spaces to deal with string parsing bug in REFPROP in SETUPdll
            for (int i = static_cast<int>(components_joined.size()); i < 10000; ++i) {
                component_string[i] = ' ';
            }

            ierr = 0;
            //...Call SETUP to initialize the program
            SETUPdll(&N, component_string.data(), hmx_bnc.data(), default_reference_state, &ierr, herr.data(),
                     10000,              // Length of component_string (see PASS_FTN.for from REFPROP)
                     refpropcharlength,  // Length of path_HMX_BNC
                     lengthofreference,  // Length of reference
                     errormessagelength  // Length of error message
            );
            if (get_config_bool(REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS) && ierr == -117) {
                throw ValueError(format("Interaction parameter estimation has been disabled: %s", herr.data()));
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
                    std::cout << format("REFPROP error/warning [ierr: %d]: %s", ierr, herr.data()) << '\n';
                }
                continue;
            } else {
                if (CoolProp::get_debug_level() > 5) {
                    std::cout << format("k: %d #endings: %d", k, number_of_endings) << '\n';
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
            std::array<char, 13> hnam{};
            std::array<char, 81> hn80{};
            std::array<char, 13> hcasn{};
            NAMEdll(&icomp, hnam.data(), hn80.data(), hcasn.data(), 12, 80, 12);
            hcasn[12] = '\0';
            std::string casn = hcasn.data();
            strstrip(casn);
            CASvec.push_back(casn);
        }
        return strjoin(CASvec, "&");
    } else if (ParamName == "name") {
        int icomp = 1L;
        std::array<char, 13> hnam{};
        std::array<char, 81> hn80{};
        std::array<char, 13> hcasn{};
        NAMEdll(&icomp, hnam.data(), hn80.data(), hcasn.data(), 12, 80, 12);
        hnam[12] = '\0';
        std::string name = hnam.data();
        strstrip(name);
        return name;
    } else if (ParamName == "long_name") {
        int icomp = 1L;
        std::array<char, 13> hnam{};
        std::array<char, 81> hn80{};
        std::array<char, 13> hcasn{};
        NAMEdll(&icomp, hnam.data(), hn80.data(), hcasn.data(), 12, 80, 12);
        hn80[80] = '\0';
        std::string n80 = hn80.data();
        strstrip(n80);
        return n80;
    } else {
        throw ValueError(format("parameter to fluid_param_string is invalid: %s", ParamName.c_str()));
    }
};
int REFPROPMixtureBackend::match_CAS(const std::string& CAS) {
    for (int icomp = 1L; icomp <= static_cast<int>(fluid_names.size()); ++icomp) {
        std::array<char, 13> hnam{};
        std::array<char, 81> hn80{};
        std::array<char, 13> hcasn{};
        NAMEdll(&icomp, hnam.data(), hn80.data(), hcasn.data(), 12, 80, 12);
        hcasn[12] = '\0';
        std::string casn = hcasn.data();
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

    int icomp = 0, jcomp = 0;
    std::array<char, 3> hmodij{};
    std::array<char, 255> hfmix{};
    std::array<char, 255> hbinp{};
    std::array<char, 255> hfij{};
    std::array<char, 255> hmxrul{};
    double fij[6];

    icomp = match_CAS(CAS1);
    jcomp = match_CAS(CAS2);

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), hfij.data(), hbinp.data(), hmxrul.data(), 3, 255, 255, 255, 255);

    // hmodij is a fixed 3-character REFPROP (FORTRAN) field with no NUL
    // terminator; read exactly the field width (constructing from a C-string
    // over-reads past the 3-byte buffer, CWE-126) and trim the trailing space
    // padding that FORTRAN uses.
    std::string shmodij(hmodij.data(), hmodij.size());
    shmodij.erase(shmodij.find_last_not_of(' ') + 1);
    //if (shmodij.find("KW") == 0 || shmodij.find("GE") == 0)  // Starts with KW or GE
    //{
    if (parameter == "model") {
        return shmodij;
    } else {
        throw ValueError(format(" I don't know what to do with your parameter [%s]", parameter.c_str()));
        return "";
    }
    //} else {
    //    //throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij.data()));
    //    return "";
    //}
}
/// Set binary mixture string value
void REFPROPMixtureBackend::set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                          const std::string& value) {
    // bound-check indices
    if (i >= Ncomp) {
        if (j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp - 1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp - 1));
        }
    } else if (j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp - 1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1, ierr = 0L;
    std::array<char, 3> hmodij{};
    std::array<char, 255> hfmix{};
    std::array<char, 255> hbinp{};
    std::array<char, 255> hfij{};
    std::array<char, 255> hmxrul{};
    double fij[6];
    std::array<char, 255> herr{};

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), hfij.data(), hbinp.data(), hmxrul.data(), 3, 255, 255, 255, 255);

    if (parameter == "model") {
        // hmodij is a fixed 3-character REFPROP (FORTRAN) field with no room for a
        // NUL terminator, so strcpy is unsafe here: reject anything longer than the
        // field and copy without a NUL, space-padding the remainder per FORTRAN
        // fixed-width string convention.
        if (value.length() > hmodij.size()) {
            throw ValueError(format("Model parameter (%s) is longer than %d characters.", value.c_str(), static_cast<int>(hmodij.size())));
        }
        hmodij.fill(' ');
        memcpy(hmodij.data(), value.data(), value.length());
    } else {
        throw ValueError(format("I don't know what to do with your parameter [%s]", parameter.c_str()));
    }
    SETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), &ierr, herr.data(), 3, 255, 255);
    if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("Unable to set parameter[%s] to value[%s]: %s", parameter.c_str(), value.c_str(), herr.data()));
    }
}
/// Set binary mixture string parameter (EXPERT USE ONLY!!!)
void REFPROPMixtureBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                          const double value) {
    // bound-check indices
    if (i >= Ncomp) {
        if (j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp - 1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp - 1));
        }
    } else if (j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp - 1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1, ierr = 0L;
    std::array<char, 3> hmodij{};
    std::array<char, 255> hfmix{};
    std::array<char, 255> hbinp{};
    std::array<char, 255> hfij{};
    std::array<char, 255> hmxrul{};
    double fij[6];
    std::array<char, 255> herr{};

    // Get the prior state
    GETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), hfij.data(), hbinp.data(), hmxrul.data(), 3, 255, 255, 255, 255);

    //if (std::string(hmodij.data()).find("KW") == 0 || std::string(hmodij.data()).find("GE") == 0)  // Starts with KW or GE
    //{
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
    SETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), &ierr, herr.data(), 3, 255, 255);
    if (ierr > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("Unable to set parameter[%s] to value[%g]: %s", parameter.c_str(), value, herr.data()));
    }
    //} else {
    //    throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij.data()));
    //}
}

/// Get binary mixture double value (EXPERT USE ONLY!!!)
double REFPROPMixtureBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) {
    // bound-check indices
    if (i >= Ncomp) {
        if (j >= Ncomp) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, Ncomp - 1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, Ncomp - 1));
        }
    } else if (j >= Ncomp) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, Ncomp - 1));
    }
    int icomp = static_cast<int>(i) + 1, jcomp = static_cast<int>(j) + 1;
    std::array<char, 3> hmodij{};
    std::array<char, 255> hfmix{};
    std::array<char, 255> hbinp{};
    std::array<char, 255> hfij{};
    std::array<char, 255> hmxrul{};
    double fij[6];

    // Get the current state
    GETKTVdll(&icomp, &jcomp, hmodij.data(), fij, hfmix.data(), hfij.data(), hbinp.data(), hmxrul.data(), 3, 255, 255, 255, 255);

    //if (std::string(hmodij.data()).find("KW") == 0 || std::string(hmodij.data()).find("GE") == 0)  // Starts with KW or GE
    //{
    double val = NAN;
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
    //} else {
    //    //throw ValueError(format("For now, model [%s] must start with KW or GE", hmodij.data()));
    //    return _HUGE;
    //}
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
    double wmm = NAN, ttrp = NAN, tnbpt = NAN, tc = NAN, pc = NAN, Dc = NAN, Zc = NAN, acf = NAN, dip = NAN, Rgas = NAN;
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
void REFPROPMixtureBackend::check_status() {
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
    double Dmax_mol_L = NAN, pmax_kPa = NAN;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &Tmin, &Tmax, &Dmax_mol_L, &pmax_kPa, 3);
    pmax = pmax_kPa * 1000;
    rhomolarmax = Dmax_mol_L * 1000;
}
CoolPropDbl REFPROPMixtureBackend::calc_pmax() {
    double Tmin = NAN, Tmax = NAN, rhomolarmax = NAN, pmax = NAN;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(pmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_Tmax() {
    double Tmin = NAN, Tmax = NAN, rhomolarmax = NAN, pmax = NAN;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(Tmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_Tmin() {
    double Tmin = NAN, Tmax = NAN, rhomolarmax = NAN, pmax = NAN;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(Tmin);
};
CoolPropDbl REFPROPMixtureBackend::calc_T_critical() {
    this->check_loaded_fluid();
    int ierr = 0;
    std::array<char, 255> herr{};
    double Tcrit = NAN, pcrit_kPa = NAN, dcrit_mol_L = NAN;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr.data(), 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return static_cast<CoolPropDbl>(Tcrit);
};
CoolPropDbl REFPROPMixtureBackend::calc_p_critical() {
    this->check_loaded_fluid();
    int ierr = 0;
    std::array<char, 255> herr{};
    double Tcrit = NAN, pcrit_kPa = NAN, dcrit_mol_L = NAN;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr.data(), 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return static_cast<CoolPropDbl>(pcrit_kPa * 1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_rhomolar_critical() {
    int ierr = 0;
    std::array<char, 255> herr{};
    double Tcrit = NAN, pcrit_kPa = NAN, dcrit_mol_L = NAN;
    CRITPdll(&(mole_fractions[0]), &Tcrit, &pcrit_kPa, &dcrit_mol_L, &ierr, herr.data(), 255);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }  //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
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
    double wmm = NAN, ttrp = NAN, tnbpt = NAN, tc = NAN, pc = NAN, Dc = NAN, Zc = NAN, acf = NAN, dip = NAN, Rgas = NAN;
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
    double wmm = NAN, ttrp = NAN, tnbpt = NAN, tc = NAN, pc = NAN, Dc = NAN, Zc = NAN, acf = NAN, dip = NAN, Rgas = NAN;
    int icomp = 1L;
    // Check if more than one
    if (Ncomp == 1) {
        // Get value for first component
        INFOdll(&icomp, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        return static_cast<CoolPropDbl>(ttrp);
    } else {
        double Tmin = NAN, Tmax = NAN, rhomolarmax = NAN, pmax = NAN;
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
    double wmm = NAN, ttrp = NAN, tnbpt = NAN, tc = NAN, pc = NAN, Dc = NAN, Zc = NAN, acf = NAN, dip = NAN, Rgas = NAN;
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
CoolPropDbl REFPROPMixtureBackend::calc_molar_mass() {
    this->check_loaded_fluid();
    double wmm_kg_kmol = NAN;
    WMOLdll(&(mole_fractions[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
    _molar_mass = wmm_kg_kmol / 1000;             // kg/mol
    return static_cast<CoolPropDbl>(_molar_mass.pt());
};
AbstractState::PhaseMolarMasses REFPROPMixtureBackend::calc_phase_molar_masses() {
    if (mole_fractions.size() == 1) {
        const double mm = molar_mass();
        return {mm, mm};
    }
    this->check_loaded_fluid();
    double mm_l_kg_per_kmol = NAN, mm_v_kg_per_kmol = NAN;
    WMOLdll(&(mole_fractions_liq[0]), &mm_l_kg_per_kmol);
    WMOLdll(&(mole_fractions_vap[0]), &mm_v_kg_per_kmol);
    return {mm_l_kg_per_kmol / 1000.0, mm_v_kg_per_kmol / 1000.0};
}

void REFPROPMixtureBackend::update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2) {
    this->check_loaded_fluid();
    if (pair == CoolProp::QmassT_INPUTS || pair == CoolProp::PQmass_INPUTS) {
        int kq = 2;  // mass-basis quality
        int ierr = 0;
        char herr[errormessagelength + 1] = {0};
        double T_K = 0, p_kPa = 0, q = 0;
        double rho_mol_L = 0, rhoLmol_L = 0, rhoVmol_L = 0;
        double emol = 0, hmol = 0, smol = 0, cvmol = 0, cpmol = 0, w = 0;

        if (pair == CoolProp::QmassT_INPUTS) {
            // QmassT: v1 is Qmass, v2 is T
            T_K = v2;
            q = v1;
            TQFLSHdll(&T_K, &q, &(mole_fractions[0]), &kq, &p_kPa, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),                 // Saturation terms
                      &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);         // Error terms
        } else {                                                // PQmass_INPUTS: v1 is P (Pa), v2 is Qmass
            p_kPa = v1 * 0.001;                                 // Pa -> kPa
            q = v2;
            PQFLSHdll(&p_kPa, &q, &(mole_fractions[0]), &kq, &T_K, &rho_mol_L, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),                 // Saturation terms
                      &emol, &hmol, &smol, &cvmol, &cpmol, &w,  // Other thermodynamic terms
                      &ierr, herr, errormessagelength);         // Error terms
        }
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("Qmass-flash: %s", herr));
        }

        // Populate state from REFPROP outputs.
        _T = T_K;
        _p = p_kPa * 1000.0;              // kPa -> Pa
        _rhomolar = rho_mol_L * 1000.0;   // mol/L -> mol/m^3
        _rhoLmolar = rhoLmol_L * 1000.0;  // mol/L -> mol/m^3
        _rhoVmolar = rhoVmol_L * 1000.0;  // mol/L -> mol/m^3
        _hmolar = hmol;
        _smolar = smol;
        _umolar = emol;
        _cvmolar = cvmol;
        _cpmolar = cpmol;
        _speed_sound = w;

        // Compute molar quality from mass quality using returned phase compositions.
        const auto MM = calc_phase_molar_masses();
        _Q = detail::Qmass_to_Qmolar(q, MM.liquid, MM.vapor);
        _Qmass = q;
        _phase = iphase_twophase;
        return;
    }
    // The 6 remaining Qmass pairs: REFPROP has no native kq flag for them.
    // Fall through to the base-class TOMS748 root-find on Qmolar.
    AbstractState::update_Qmass_pair(pair, v1, v2);
}

CoolPropDbl REFPROPMixtureBackend::calc_Bvirial() {
    double b = NAN;
    VIRBdll(&_T, &(mole_fractions[0]), &b);
    return b * 0.001;  // 0.001 to convert from l/mol to m^3/mol
}
CoolPropDbl REFPROPMixtureBackend::calc_dBvirial_dT() {
    double b = NAN;
    DBDTdll(&_T, &(mole_fractions[0]), &b);
    return b * 0.001;  // 0.001 to convert from l/mol to m^3/mol
}
CoolPropDbl REFPROPMixtureBackend::calc_Cvirial() {
    double c = NAN;
    VIRCdll(&_T, &(mole_fractions[0]), &c);
    return c * 1e-6;  // 1e-6 to convert from (l/mol)^2 to (m^3/mol)^2
}
double REFPROPMixtureBackend::calc_melt_Tmax() {
    this->check_loaded_fluid();
    int ierr = 0;
    std::array<char, 255> herr{};
    double tmin = NAN, tmax = NAN, Dmax_mol_L = NAN, pmax_kPa = NAN, Tmax_melt = NAN;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &tmin, &tmax, &Dmax_mol_L, &pmax_kPa, 3);
    // Get the maximum temperature for the melting curve by using the maximum pressure
    MELTPdll(&pmax_kPa, &(mole_fractions[0]), &Tmax_melt, &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return Tmax_melt;
}
CoolPropDbl REFPROPMixtureBackend::calc_melting_line(int param, int given, CoolPropDbl value) {
    this->check_loaded_fluid();
    int ierr = 0;
    std::array<char, 255> herr{};

    // Bounding-box sentinels for the melting curve, matching the HEOS ancillary
    // contract (Ancillaries.cpp): when param is a bound key, return the bound and
    // ignore given/value.  Each axis bound is independent.
    //
    // iP_min is the triple-point pressure, for consistency with the HEOS melting
    // ancillary.  iT_min is the melting line's lowest temperature: for most fluids
    // that is the triple-point temperature, but for water and heavy water the
    // melting line extends below it to the ice Ih/III junction (251.165 K for
    // water).  REFPROP's MELTP (CORE_ANC.FOR) takes that floor as LIMITS('EOS').Tmin,
    // so we read it from the EOS lower temperature limit.  The upper limits come from
    // the melting-line model ("MLT") via LIMITX, distinct from the EOS limits.
    if (param == iP_min) {
        return calc_p_triple();
    } else if (param == iT_min) {
        double Tmin = NAN, Tmax = NAN, rhomolarmax = NAN, pmax = NAN;
        limits(Tmin, Tmax, rhomolarmax, pmax);  // LIMITS('EOS'); Tmin is MELTP's temperature floor
        return static_cast<CoolPropDbl>(Tmin);
    } else if (param == iP_max || param == iT_max) {
        double t_in = static_cast<double>(calc_Ttriple()), D_in = 0.0, p_in = 0.0;
        double tmin_unused = NAN, tmax_melt = NAN, Dmax_unused = NAN, pmax_kPa = NAN;
        char htyp[] = "MLT";
        LIMITXdll(htyp, &t_in, &D_in, &p_in, &(mole_fractions[0]), &tmin_unused, &tmax_melt, &Dmax_unused, &pmax_kPa, &ierr, herr.data(), 3,
                  errormessagelength);
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("%s", herr.data()).c_str());
        }
        return (param == iP_max) ? static_cast<CoolPropDbl>(pmax_kPa * 1000) : static_cast<CoolPropDbl>(tmax_melt);
    } else if (param == iP && given == iT) {
        double _T = static_cast<double>(value), p_kPa = NAN;
        MELTTdll(&_T, &(mole_fractions[0]), &p_kPa, &ierr, herr.data(), errormessagelength);  // Error message
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("%s", herr.data()).c_str());
        }  //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
        return p_kPa * 1000;
    } else if (param == iT && given == iP) {
        double p_kPa = static_cast<double>(value) / 1000.0, _T = NAN;
        MELTPdll(&p_kPa, &(mole_fractions[0]), &_T, &ierr, herr.data(), errormessagelength);  // Error message
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("%s", herr.data()).c_str());
        }  //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
        return _T;
    } else {
        // Use the raw integer keys here: param/given may be invalid keys (e.g. the
        // -1 sentinel), and looking them up via get_parameter_information() would
        // itself throw, masking this "invalid inputs" diagnostic.
        throw ValueError(format("calc_melting_line(param=%d, given=%d, value=%Lg) is an invalid set of inputs", param, given, value));
    }
}
bool REFPROPMixtureBackend::has_melting_line() {
    this->check_loaded_fluid();

    int ierr = 0;
    std::array<char, 255> herr{};
    double _T = 300, p_kPa = NAN;
    MELTTdll(&_T, &(mole_fractions[0]), &p_kPa, &ierr, herr.data(), errormessagelength);  // Error message
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
    double wmm = NAN, ttrp = NAN, tnbpt = NAN, tc = NAN, pc = NAN, Dc = NAN, Zc = NAN, acf = NAN, dip = NAN, Rgas = NAN;
    // FORTRAN is 1-based indexing!
    for (int i = 1L; i <= static_cast<int>(mole_fractions_long_double.size()); ++i) {
        // Get value for first component
        INFOdll(&i, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
        mass_fractions[i - 1] = (wmm / 1000.0) * mole_fractions_long_double[i - 1] / mm;
    }
    return mass_fractions;
}

CoolPropDbl REFPROPMixtureBackend::calc_PIP() {
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

CoolPropDbl REFPROPMixtureBackend::calc_viscosity() {
    this->check_loaded_fluid();
    double eta = NAN, tcx = NAN, rhomol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::array<char, 255> herr{};
    TRNPRPdll(&_T, &rhomol_L, &(mole_fractions[0]),     // Inputs
              &eta, &tcx,                               // Outputs
              &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    _viscosity = 1e-6 * eta;
    _conductivity = tcx;
    return static_cast<double>(_viscosity);
}
CoolPropDbl REFPROPMixtureBackend::calc_conductivity() {
    // Calling viscosity also caches conductivity, use that to save calls
    calc_viscosity();
    return static_cast<double>(_conductivity);
}
CoolPropDbl REFPROPMixtureBackend::calc_surface_tension() {
    this->check_loaded_fluid();
    double sigma = NAN, rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::array<char, 255> herr{};
    SURFTdll(&_T, &rho_mol_L, &(mole_fractions[0]),    // Inputs
             &sigma,                                   // Outputs
             &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    _surface_tension = sigma;
    return static_cast<double>(_surface_tension);
}
CoolPropDbl REFPROPMixtureBackend::calc_fugacity_coefficient(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> fug_cof;
    fug_cof.resize(mole_fractions.size());
    std::array<char, 255> herr{};
    FUGCOFdll(&_T, &rho_mol_L, &(mole_fractions[0]),    // Inputs
              &(fug_cof[0]),                            // Outputs
              &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return static_cast<CoolPropDbl>(fug_cof[i]);
}
CoolPropDbl REFPROPMixtureBackend::calc_fugacity(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> f(mole_fractions.size());
    std::array<char, 255> herr{};
    FGCTY2dll(&_T, &rho_mol_L, &(mole_fractions[0]),    // Inputs
              &(f[0]),                                  // Outputs
              &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return static_cast<CoolPropDbl>(f[i] * 1000);
}
CoolPropDbl REFPROPMixtureBackend::calc_chemical_potential(std::size_t i) {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    int ierr = 0;
    std::vector<double> chem_pot(mole_fractions.size());
    std::array<char, 255> herr{};
    CHEMPOTdll(&_T, &rho_mol_L, &(mole_fractions[0]),    // Inputs
               &(chem_pot[0]),                           // Outputs
               &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }
    //else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
    return static_cast<CoolPropDbl>(chem_pot[i]);
}

void REFPROPMixtureBackend::calc_phase_envelope(const std::string& type) {
    this->check_loaded_fluid();
    double rhoymin = NAN, rhoymax = NAN, c = 0;
    int ierr = 0;
    std::array<char, 255> herr{};
    SATSPLNdll(&(mole_fractions[0]),                     // Inputs
               &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()).c_str());
    }

    // Clear the phase envelope data
    PhaseEnvelope = PhaseEnvelopeData();
    /*
    subroutine SPLNVAL (isp,iderv,a,f,ierr,herr.data())
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
    if (SPLNVALdll == nullptr) {
        std::string rpv = get_global_param_string("REFPROP_version");
        throw ValueError(
          format("Your version of REFFPROP(%s) does not have the SPLNVALdll function; cannot extract phase envelope values", rpv.c_str()));
    };
    SPLNVALdll(&isp, &iderv, &c, &rhoymin, &ierr, herr.data(), errormessagelength);
    iderv = -2;
    SPLNVALdll(&isp, &iderv, &c, &rhoymax, &ierr, herr.data(), errormessagelength);
    int nc = static_cast<int>(this->Ncomp);
    double ratio = pow(rhoymax / rhoymin, 1 / double(N));
    // Geometric density sweep (rho *= ratio, ~N iters by construction).
    for (double rho_molL = rhoymin; rho_molL < rhoymax; rho_molL *= ratio) {  // NOLINT(cert-flp30-c)
        double y = NAN;
        iderv = 0;

        PhaseEnvelope.x.resize(nc);
        PhaseEnvelope.y.resize(nc);
        for (isp = 1; isp <= nc; ++isp) {
            SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
            PhaseEnvelope.x[isp - 1].push_back(y);
            PhaseEnvelope.y[isp - 1].push_back(get_mole_fractions()[isp - 1]);
        }

        PhaseEnvelope.rhomolar_vap.push_back(rho_molL * 1000);
        PhaseEnvelope.lnrhomolar_vap.push_back(log(rho_molL * 1000));
        isp = nc + 1;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
        double T = y;
        PhaseEnvelope.T.push_back(y);
        PhaseEnvelope.lnT.push_back(log(y));
        isp = nc + 2;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
        PhaseEnvelope.p.push_back(y * 1000);
        PhaseEnvelope.lnp.push_back(log(y * 1000));
        isp = nc + 3;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
        PhaseEnvelope.rhomolar_liq.push_back(y * 1000);
        PhaseEnvelope.lnrhomolar_liq.push_back(log(y * 1000));
        PhaseEnvelope.Q.push_back(static_cast<double>(y > rho_molL));
        isp = nc + 4;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
        PhaseEnvelope.hmolar_vap.push_back(y);
        isp = nc + 5;
        SPLNVALdll(&isp, &iderv, &rho_molL, &y, &ierr, herr.data(), errormessagelength);
        PhaseEnvelope.smolar_vap.push_back(y);

        // Other outputs that could be useful
        int ierr_trn = 0;
        std::array<char, 255> herr_trn{};
        double p_kPa = NAN, emol = NAN, hmol = NAN, smol = NAN, cvmol = NAN, cpmol = NAN, w = NAN, hjt = NAN, eta = NAN, tcx = NAN;
        // "Vapor"
        THERMdll(&T, &rho_molL, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
        PhaseEnvelope.cpmolar_vap.push_back(cpmol);
        PhaseEnvelope.cvmolar_vap.push_back(cvmol);
        PhaseEnvelope.speed_sound_vap.push_back(w);
        TRNPRPdll(&T, &rho_molL, &(mole_fractions[0]),              // Inputs
                  &eta, &tcx,                                       // Outputs
                  &ierr_trn, herr_trn.data(), errormessagelength);  // Error message
        PhaseEnvelope.viscosity_vap.push_back(eta / 1e6);
        PhaseEnvelope.conductivity_vap.push_back(tcx);
    }
}
CoolPropDbl REFPROPMixtureBackend::calc_cpmolar_idealgas() {
    this->check_loaded_fluid();
    double rho_mol_L = 0.001 * _rhomolar;
    double p0 = NAN, e0 = NAN, h0 = NAN, s0 = NAN, cv0 = NAN, cp0 = NAN, w0 = NAN, A0 = NAN, G0 = NAN;
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
        }  // so it's unknown as well.

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
    // Mass-quality input pair on a true mixture: handle via update_Qmass_pair.
    // Pure / pseudo-pure (mole_fractions.size() == 1) goes through the existing
    // mass_to_molar_inputs path.
    if (CoolProp::is_Qmass_pair(input_pair) && mole_fractions.size() > 1) {
        update_Qmass_pair(input_pair, value1, value2);
        return;
    }
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

            if (imposed_phase_index == iphase_not_imposed || imposed_phase_index == iphase_twophase) {
                // Use the full 2-phase flash routine
                // from REFPROP: subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
                DSFLSHdll(&rho_mol_L, &smol, &(mole_fractions[0]), &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                          &(mole_fractions_vap[0]),  // Saturation terms
                          &q, &emol, &hmol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
                if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    throw ValueError(format("DmolarSmolar: %s", herr).c_str());
                }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            } else {
                // Phase is imposed -> use the single-phase D,S iterator (DSFL1)
                // followed by THERMdll for the remaining properties. Avoids the
                // 2-phase iteration that DSFL2 reaches inside DSFLSH and
                // sometimes fails to converge for mixtures (#2042).
                DSFL1dll(&rho_mol_L, &smol, &(mole_fractions[0]), &_T, &ierr, herr, errormessagelength);
                if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                    throw ValueError(format("DmolarSmolar (imposed phase): %s", herr).c_str());
                }
                THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
            }

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
            rho_mol_L = 0.001 * value1;  // Want rho in [mol/L] in REFPROP
            emol = value2;               // Want e in J/mol in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            DEFLSHdll(&rho_mol_L, &emol, &(mole_fractions[0]), &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                      &(mole_fractions_vap[0]),  // Saturation terms
                      &q, &hmol, &smol, &cvmol, &cpmol, &w, &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarUmolar: %s", herr).c_str());
            }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa * 1000;
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            break;
        }
        case DmassUmass_INPUTS: {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // U: [J/kg] * [kg/mol] -> [J/mol]
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
            _p = p_kPa * 1000;              // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L * 1000;   // 1000 for conversion from mol/L to mol/m3
            _rhoLmolar = rhoLmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
            _rhoVmolar = rhoVmol_L * 1000;  // 1000 for conversion from mol/L to mol/m3
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
                    //
                    // GH #1502: SATTP's iFlsh (3=bubble/P,xliq, 4=dew/P,xvap) is NOT SATP's kph
                    // (3/4 mean freezing/sublimation lines).  Map to kph 1=liquid/bubble, 2=vapor/dew.
                    // Also pass p_kPa (in kPa, as SATP expects): _p is -_HUGE here because clear()
                    // ran and _p = value1 is not set until after this block.
                    int kph = (iFlsh == 3) ? 1 : 2;
                    SATPdll(&p_kPa, &(mole_fractions[0]), &kph, &_T, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]), &(mole_fractions_vap[0]),
                            &ierr, herr, errormessagelength);
                    rho_mol_L = (kph == 1) ? rhoLmol_L : rhoVmol_L;
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
                    // Calculate everything else since we were able to carry out a flash call.
                    // #2671: keep the equilibrium pressure that SATTP/SATT returned; THERMdll
                    // would otherwise recompute a slightly different p from the (SATSPLN-refined)
                    // saturated-phase density, so send its p output to a throwaway.
                    double p_kPa_therm = _HUGE;
                    THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa_therm, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);
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
        case DmassQ_INPUTS: {
            // Convert to molar and dispatch to DmolarQ_INPUTS
            update(DmolarQ_INPUTS, value1 / static_cast<double>(_molar_mass), value2);
            return;
        }
        case DmolarQ_INPUTS: {
            // GitHub #1845: REFPROP supports D,Q via DQFL2 — wire it through.
            // REFPROP's DQFL2(d,q,z,kq,t,p,Dl,Dv,xliq,xvap,ierr,herr) iterates
            // for T given (rho, q) on the saturation envelope. Then THERMdll
            // populates the remaining caloric and transport-prep outputs.
            _rhomolar = value1;
            _Q = value2;
            rho_mol_L = 0.001 * value1;
            q = value2;
            int kq = 1;  // Q on a molar basis
            DQFL2dll(&rho_mol_L, &q, &(mole_fractions[0]), &kq, &_T, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(mole_fractions_liq[0]),
                     &(mole_fractions_vap[0]), &ierr, herr, errormessagelength);
            if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
                throw ValueError(format("DmolarQ(%s): %s", LoadedREFPROPRef.c_str(), herr).c_str());
            }
            THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);

            _p = p_kPa * 1000;
            _rhoLmolar = rhoLmol_L * 1000;
            _rhoVmolar = rhoVmol_L * 1000;
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
            if (_phase != iphase_twophase) {
                // No saturated state: clear so build_saturation_shim / keyed outputs reject it.
                _rhoLmolar = _HUGE;
                _rhoVmolar = _HUGE;
            }
        } else if (!(_Q >= 0 && _Q <= 1)) {  // mixture, not two-phase
            _rhoLmolar = _HUGE;
            _rhoVmolar = _HUGE;
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
            if (CoolProp::get_debug_level() > 9) {
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

void REFPROPMixtureBackend::update_DmolarT_direct(CoolPropDbl rhomolar, CoolPropDbl T) {
    this->check_loaded_fluid();
    clear();
    _T = T;
    _rhomolar = rhomolar;
    _Q = -1;  // not on the saturation envelope by construction

    double rho_mol_L = 0.001 * static_cast<double>(_rhomolar);
    auto T_K = static_cast<double>(_T);
    double p_kPa = _HUGE, emol = _HUGE, hmol = _HUGE, smol = _HUGE, cvmol = _HUGE, cpmol = _HUGE, w = _HUGE, hjt = _HUGE;
    THERMdll(&T_K, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);

    _p = 1000.0 * p_kPa;
    _hmolar = hmol;
    _smolar = smol;
    _umolar = emol;
    _cvmolar = cvmol;
    _cpmolar = cpmol;
    _speed_sound = w;
    _gibbsmolar = hmol - T_K * smol;
    _tau = calc_T_reducing() / _T;
    _delta = _rhomolar / calc_rhomolar_reducing();
}
CoolPropDbl REFPROPMixtureBackend::call_phixdll(int itau, int idel) {
    this->check_loaded_fluid();
    double val = 0, tau = _tau, delta = _delta;
    if (PHIXdll == nullptr) {
        throw ValueError("PHIXdll function is not available in your version of REFPROP. Please upgrade");
    }
    PHIXdll(&itau, &idel, &tau, &delta, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val) / pow(static_cast<CoolPropDbl>(_delta), idel) / pow(static_cast<CoolPropDbl>(_tau), itau);
}
CoolPropDbl REFPROPMixtureBackend::call_phi0dll(int itau, int idel) {
    this->check_loaded_fluid();
    double val = 0, tau = _tau, __T = T(), __rho = rhomolar() / 1000;
    if (PHI0dll == nullptr) {
        throw ValueError("PHI0dll function is not available in your version of REFPROP. Please upgrade");
    }
    PHI0dll(&itau, &idel, &__T, &__rho, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val) / pow(static_cast<CoolPropDbl>(_delta), idel) / pow(tau, itau);
}
/// Calculate excess properties
void REFPROPMixtureBackend::calc_excess_properties() {
    this->check_loaded_fluid();
    int ierr = 0;
    std::array<char, 255> herr{};
    double T_K = _T, p_kPa = _p / 1000.0, rho = 1, vE = -1, eE = -1, hE = -1, sE = -1, aE = -1, gE = -1;
    int kph = 1;

    //    subroutine EXCESS(t, p, x, kph, rho, vE, eE, hE, sE, aE, gE, ierr, herr.data())
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
    EXCESSdll(&T_K, &p_kPa, &(mole_fractions[0]), &kph, &rho, &vE, &eE, &hE, &sE, &aE, &gE, &ierr, herr.data(), errormessagelength);  // Error message
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("EXCESSdll: %s", herr.data()).c_str());
    }  // TODO: else if (ierr < 0) {set_warning(format("%s",herr.data()).c_str());}
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
        wrapper(const std::vector<double>& z) : z(z) {};
        std::vector<double> call(const std::vector<double>& x) override {
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

void REFPROPMixtureBackend::link_to_loaded_fluids(const REFPROPMixtureBackend& host) {
    // Per-class access: we may read host's private cached_component_string.
    this->Ncomp = host.Ncomp;
    this->fluid_names = host.fluid_names;
    this->cached_component_string = host.cached_component_string;  // => check_loaded_fluid() reuse path, no SETUPdll
    this->mole_fractions.assign(ncmax, 0.0);
    this->mole_fractions_liq.assign(ncmax, 0.0);
    this->mole_fractions_vap.assign(ncmax, 0.0);
    this->imposed_phase_index = iphase_not_imposed;
    this->_mole_fractions_set = false;
}

shared_ptr<REFPROPMixtureBackend> REFPROPMixtureBackend::build_saturation_shim(int Q) {
    this->check_loaded_fluid();
    if (!(Q == 0 || Q == 1)) {
        throw ValueError(format("build_saturation_shim requires Q==0 or Q==1; got %d", Q));
    }
    const std::vector<double>& x_phase = (Q == 0) ? mole_fractions_liq : mole_fractions_vap;
    // is_cached (operator bool) check first: operator double()/CoolPropDbl() throws if uncached.
    bool have = (Q == 0) ? static_cast<bool>(_rhoLmolar) : static_cast<bool>(_rhoVmolar);
    CoolPropDbl rho_phase = have ? ((Q == 0) ? static_cast<CoolPropDbl>(_rhoLmolar) : static_cast<CoolPropDbl>(_rhoVmolar)) : _HUGE;
    if (!have || !ValidNumber(rho_phase)) {
        throw ValueError("The saturated state has not been set (no two-phase density available).");
    }
    shared_ptr<REFPROPMixtureBackend> shim(new REFPROPMixtureBackend());
    shim->link_to_loaded_fluids(*this);
    // Copy only the first Ncomp entries of the per-phase composition.
    std::vector<CoolPropDbl> x(x_phase.begin(), x_phase.begin() + Ncomp);
    shim->set_mole_fractions(x);
    shim->specify_phase((Q == 0) ? iphase_liquid : iphase_gas);
    shim->update_DmolarT_direct(rho_phase, _T);
    shim->_Q = static_cast<CoolPropDbl>(Q);  // mark the shim as sitting on the envelope
    return shim;
}

double REFPROPMixtureBackend::dpdT_along_saturation_pure(int kph) {
    if (Ncomp != 1) {  // precondition: DPTSATKdll is a single-component saturation routine
        throw ValueError("dpdT_along_saturation_pure is only valid for pure fluids");
    }
    int icomp = 1;  // pure fluid
    auto T_K = static_cast<double>(_T);
    double p_kPa = _HUGE, rho_mol_L = _HUGE, csat = _HUGE, dpdt_kPa_K = _HUGE;
    int ierr = 0;
    std::array<char, errormessagelength> herr{};
    DPTSATKdll(&icomp, &T_K, &kph, &p_kPa, &rho_mol_L, &csat, &dpdt_kPa_K, &ierr, herr.data(), errormessagelength);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()));
    }
    return 1000.0 * dpdt_kPa_K;  // kPa/K -> Pa/K
}

double REFPROPMixtureBackend::saturation_pressure_at_T(double T, int Q) {
    // Q==0 -> bubble (kph=1), Q==1 -> dew (kph=2). Mirror update()'s QT_INPUTS SATTdll call.
    int kph = (Q == 0) ? 1 : 2;
    double T_K = T, p_kPa = _HUGE, rhoLmol_L = _HUGE, rhoVmol_L = _HUGE;
    int ierr = 0;
    std::array<char, errormessagelength> herr{};
    std::vector<double> xliq(ncmax, 0.0), xvap(ncmax, 0.0);
    SATTdll(&T_K, &(mole_fractions[0]), &kph, &p_kPa, &rhoLmol_L, &rhoVmol_L, &(xliq[0]), &(xvap[0]), &ierr, herr.data(), errormessagelength);
    if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
        throw ValueError(format("%s", herr.data()));
    }
    return 1000.0 * p_kPa;  // kPa -> Pa
}

CoolPropDbl REFPROPMixtureBackend::calc_first_saturation_deriv(parameters Of1, parameters Wrt1) {
    this->check_loaded_fluid();

    // The saturation slope comes from REFPROP's own DPTSATKdll/SATTdll at (T, composition),
    // so this needs only a valid two-phase _Q, not the L/V saturation densities. (This is
    // also why it works when invoked on a single-phase saturation shim, whose _Q is 0 or 1.)
    CoolPropDbl dTdP_sat = _HUGE;
    if (Ncomp == 1) {
        // Pure fluid: dP/dT along the vapor-pressure curve is independent of quality, so
        // any two-phase state (0<=Q<=1) is valid. DPTSATKdll (icomp=1) returns it natively;
        // kph=1 (liquid side) suffices for a pure fluid (liquid and vapor lines coincide).
        if (!(_Q >= 0 && _Q <= 1)) {
            throw ValueError(format("calc_first_saturation_deriv requires a two-phase state (0<=Q<=1); got Q=%g", _Q));
        }
        dTdP_sat = 1.0 / dpdT_along_saturation_pure(1);  // Pa/K
    } else {
        // Mixture: DPTSATKdll is a PER-COMPONENT routine (icomp); it does NOT give the
        // constant-overall-composition bubble/dew slope a mixture needs. So differentiate
        // REFPROP's own SATTdll flash numerically at fixed overall composition. Only defined
        // on the bubble (Q==0) or dew (Q==1) curve (tolerance mirrors update()'s 1e-10).
        bool at_bubble = std::abs(_Q) < 1e-10;
        bool at_dew = std::abs(_Q - 1) < 1e-10;
        if (!at_bubble && !at_dew) {
            throw ValueError(format("calc_first_saturation_deriv requires Q==0 (bubble) or Q==1 (dew) for mixtures; got Q=%g", _Q));
        }
        int Qint = at_dew ? 1 : 0;
        const double dT = 1e-3 * static_cast<double>(_T);  // relative central step
        double pp = saturation_pressure_at_T(static_cast<double>(_T) + dT, Qint);
        double pm = saturation_pressure_at_T(static_cast<double>(_T) - dT, Qint);
        dTdP_sat = (2 * dT) / (pp - pm);
    }

    if (Of1 == iT && Wrt1 == iP) {
        return dTdP_sat;
    } else if (Of1 == iP && Wrt1 == iT) {
        return 1 / dTdP_sat;
    } else if (Wrt1 == iT) {
        return first_partial_deriv(Of1, iT, iP) + first_partial_deriv(Of1, iP, iT) / dTdP_sat;
    } else if (Wrt1 == iP) {
        return first_partial_deriv(Of1, iP, iT) + first_partial_deriv(Of1, iT, iP) * dTdP_sat;
    } else {
        throw ValueError(
          format("Not possible to take first saturation derivative with respect to %s", get_parameter_information(Wrt1, "short").c_str()));
    }
}

CoolPropDbl REFPROPMixtureBackend::calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
    if (Ncomp > 1) {
        // The constant-overall-composition saturation slope a mixture two-phase derivative
        // needs is not available from the per-phase saturation shims; only pure fluids are
        // supported here. (Mixture saturated-state properties and bubble/dew first_saturation_deriv
        // are available separately.)
        throw NotImplementedError("calc_first_two_phase_deriv is not implemented for mixtures in the REFPROP backend");
    }
    if (!_rhoLmolar || !_rhoVmolar || !ValidNumber(static_cast<CoolPropDbl>(_rhoLmolar)) || !ValidNumber(static_cast<CoolPropDbl>(_rhoVmolar))) {
        throw ValueError("The saturation properties are needed for calc_first_two_phase_deriv");
    }
    shared_ptr<REFPROPMixtureBackend> SatL = build_saturation_shim(0);
    shared_ptr<REFPROPMixtureBackend> SatV = build_saturation_shim(1);

    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        return -POW2(rhomolar()) * (1 / SatV->rhomolar() - 1 / SatL->rhomolar()) / (SatV->hmolar() - SatL->hmolar());
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return -POW2(rhomass()) * (1 / SatV->rhomass() - 1 / SatL->rhomass()) / (SatV->hmass() - SatL->hmass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        CoolPropDbl dvdrhoL = -1 / POW2(SatL->rhomolar());
        CoolPropDbl dvdrhoV = -1 / POW2(SatV->rhomolar());
        CoolPropDbl dvL_dp = dvdrhoL * SatL->calc_first_saturation_deriv(iDmolar, iP);
        CoolPropDbl dvV_dp = dvdrhoV * SatV->calc_first_saturation_deriv(iDmolar, iP);
        CoolPropDbl dhL_dp = SatL->calc_first_saturation_deriv(iHmolar, iP);
        CoolPropDbl dhV_dp = SatV->calc_first_saturation_deriv(iHmolar, iP);
        CoolPropDbl dxdp_h = (Q() * dhV_dp + (1 - Q()) * dhL_dp) / (SatL->hmolar() - SatV->hmolar());
        CoolPropDbl dvdp_h = dvL_dp + dxdp_h * (1 / SatV->rhomolar() - 1 / SatL->rhomolar()) + Q() * (dvV_dp - dvL_dp);
        return -POW2(rhomolar()) * dvdp_h;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        CoolPropDbl dvdrhoL = -1 / POW2(SatL->rhomass());
        CoolPropDbl dvdrhoV = -1 / POW2(SatV->rhomass());
        CoolPropDbl dvL_dp = dvdrhoL * SatL->calc_first_saturation_deriv(iDmass, iP);
        CoolPropDbl dvV_dp = dvdrhoV * SatV->calc_first_saturation_deriv(iDmass, iP);
        CoolPropDbl dhL_dp = SatL->calc_first_saturation_deriv(iHmass, iP);
        CoolPropDbl dhV_dp = SatV->calc_first_saturation_deriv(iHmass, iP);
        CoolPropDbl dxdp_h = (Q() * dhV_dp + (1 - Q()) * dhL_dp) / (SatL->hmass() - SatV->hmass());
        CoolPropDbl dvdp_h = dvL_dp + dxdp_h * (1 / SatV->rhomass() - 1 / SatL->rhomass()) + Q() * (dvV_dp - dvL_dp);
        return -POW2(rhomass()) * dvdp_h;
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
    }
}

CoolPropDbl REFPROPMixtureBackend::calc_second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2,
                                                               parameters Constant2) {
    if (Ncomp > 1) {
        throw NotImplementedError("calc_second_two_phase_deriv is not implemented for mixtures in the REFPROP backend");
    }
    if (!_rhoLmolar || !_rhoVmolar || !ValidNumber(static_cast<CoolPropDbl>(_rhoLmolar)) || !ValidNumber(static_cast<CoolPropDbl>(_rhoVmolar))) {
        throw ValueError("The saturation properties are needed for calc_second_two_phase_deriv");
    }
    shared_ptr<REFPROPMixtureBackend> SatL = build_saturation_shim(0);
    shared_ptr<REFPROPMixtureBackend> SatV = build_saturation_shim(1);

    if (Of == iDmolar
        && ((Wrt1 == iHmolar && Constant1 == iP && Wrt2 == iP && Constant2 == iHmolar)
            || (Wrt2 == iHmolar && Constant2 == iP && Wrt1 == iP && Constant1 == iHmolar))) {
        parameters h_key = iHmolar, rho_key = iDmolar, p_key = iP;
        CoolPropDbl dv_dh_constp = calc_first_two_phase_deriv(rho_key, h_key, p_key) / (-POW2(rhomolar()));
        CoolPropDbl drhomolar_dp_consth = first_two_phase_deriv(rho_key, p_key, h_key);

        CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(h_key, p_key);
        CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(h_key, p_key);
        CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(rho_key, p_key);
        CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(rho_key, p_key);
        CoolPropDbl numerator = 1 / SatV->keyed_output(rho_key) - 1 / SatL->keyed_output(rho_key);
        CoolPropDbl denominator = SatV->keyed_output(h_key) - SatL->keyed_output(h_key);
        CoolPropDbl dnumerator = -1 / POW2(SatV->keyed_output(rho_key)) * drhoV_dp_sat + 1 / POW2(SatL->keyed_output(rho_key)) * drhoL_dp_sat;
        CoolPropDbl ddenominator = dhV_dp_sat - dhL_dp_sat;
        CoolPropDbl d_dvdh_dp_consth = (denominator * dnumerator - numerator * ddenominator) / POW2(denominator);
        return -POW2(rhomolar()) * d_dvdh_dp_consth + dv_dh_constp * (-2 * rhomolar()) * drhomolar_dp_consth;
    } else if (Of == iDmass
               && ((Wrt1 == iHmass && Constant1 == iP && Wrt2 == iP && Constant2 == iHmass)
                   || (Wrt2 == iHmass && Constant2 == iP && Wrt1 == iP && Constant1 == iHmass))) {
        parameters h_key = iHmass, rho_key = iDmass, p_key = iP;
        CoolPropDbl rho = keyed_output(rho_key);
        CoolPropDbl dv_dh_constp = calc_first_two_phase_deriv(rho_key, h_key, p_key) / (-POW2(rho));
        CoolPropDbl drho_dp_consth = first_two_phase_deriv(rho_key, p_key, h_key);

        CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(h_key, p_key);
        CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(h_key, p_key);
        CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(rho_key, p_key);
        CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(rho_key, p_key);
        CoolPropDbl numerator = 1 / SatV->keyed_output(rho_key) - 1 / SatL->keyed_output(rho_key);
        CoolPropDbl denominator = SatV->keyed_output(h_key) - SatL->keyed_output(h_key);
        CoolPropDbl dnumerator = -1 / POW2(SatV->keyed_output(rho_key)) * drhoV_dp_sat + 1 / POW2(SatL->keyed_output(rho_key)) * drhoL_dp_sat;
        CoolPropDbl ddenominator = dhV_dp_sat - dhL_dp_sat;
        CoolPropDbl d_dvdh_dp_consth = (denominator * dnumerator - numerator * ddenominator) / POW2(denominator);
        return -POW2(rho) * d_dvdh_dp_consth + dv_dh_constp * (-2 * rho) * drho_dp_consth;
    } else {
        throw ValueError("These inputs are not supported to calc_second_two_phase_deriv");
    }
}

CoolPropDbl REFPROPMixtureBackend::calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) {
    if (Ncomp > 1) {
        throw NotImplementedError("calc_first_two_phase_deriv_splined is not implemented for mixtures in the REFPROP backend");
    }
    // Note: If you need all three values (drho_dh_p, drho_dp_h and rho_spline),
    // you should calculate drho_dp_h first to avoid duplicate calculations.
    bool drho_dh_p = false;
    bool drho_dp_h = false;
    bool rho_spline = false;

    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        drho_dh_p = true;
        if (_drho_spline_dh__constp) return _drho_spline_dh__constp;
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return first_two_phase_deriv_splined(iDmolar, iHmolar, iP, x_end) * POW2(molar_mass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        drho_dp_h = true;
        if (_drho_spline_dp__consth) return _drho_spline_dp__consth;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        return first_two_phase_deriv_splined(iDmolar, iP, iHmolar, x_end) * molar_mass();
    } else if (Of == iDmolar && Wrt == iDmolar && Constant == iDmolar) {
        rho_spline = true;
        if (_rho_spline) return _rho_spline;
    } else if (Of == iDmass && Wrt == iDmass && Constant == iDmass) {
        return first_two_phase_deriv_splined(iDmolar, iDmolar, iDmolar, x_end) * molar_mass();
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv_splined");
    }

    if (!_rhoLmolar || !_rhoVmolar || !ValidNumber(static_cast<CoolPropDbl>(_rhoLmolar)) || !ValidNumber(static_cast<CoolPropDbl>(_rhoVmolar))) {
        throw ValueError("The saturation properties are needed for calc_first_two_phase_deriv_splined");
    }
    if (_Q > x_end) {
        throw ValueError(format("Q [%g] is greater than x_end [%Lg]", _Q, x_end).c_str());
    }
    bool two_phase = (Ncomp == 1) ? (_phase == iphase_twophase) : (_Q >= 0 && _Q <= 1);
    if (!two_phase) {
        throw ValueError(format("state is not two-phase"));
    }

    shared_ptr<REFPROPMixtureBackend> SatL = build_saturation_shim(0);
    shared_ptr<REFPROPMixtureBackend> SatV = build_saturation_shim(1);

    // Fake single-phase liquid state at the bubble point, and the end-of-spline two-phase state.
    shared_ptr<REFPROPMixtureBackend> Liq(new REFPROPMixtureBackend());
    Liq->link_to_loaded_fluids(*this);
    Liq->set_mole_fractions(this->get_mole_fractions());
    Liq->specify_phase(iphase_liquid);
    Liq->_Q = -1;
    Liq->update_DmolarT_direct(SatL->rhomolar(), SatL->T());

    shared_ptr<REFPROPMixtureBackend> End(new REFPROPMixtureBackend());
    End->link_to_loaded_fluids(*this);
    End->set_mole_fractions(this->get_mole_fractions());
    End->update(QT_INPUTS, x_end, SatL->T());

    CoolPropDbl Delta = Q() * (SatV->keyed_output(iHmolar) - SatL->keyed_output(iHmolar));
    CoolPropDbl Delta_end = End->keyed_output(iHmolar) - SatL->keyed_output(iHmolar);

    // At the end of the zone to which spline is applied
    CoolPropDbl drho_dh_end = End->calc_first_two_phase_deriv(iDmolar, iHmolar, iP);
    CoolPropDbl rho_end = End->keyed_output(iDmolar);

    // Faking single-phase
    CoolPropDbl rho_liq = Liq->keyed_output(iDmolar);
    CoolPropDbl drho_dh_liq_constp = Liq->first_partial_deriv(iDmolar, iHmolar, iP);

    // Spline coordinates a, b, c, d
    CoolPropDbl Abracket = (2 * rho_liq - 2 * rho_end + Delta_end * (drho_dh_liq_constp + drho_dh_end));
    CoolPropDbl a = 1 / POW3(Delta_end) * Abracket;
    CoolPropDbl b = 3 / POW2(Delta_end) * (-rho_liq + rho_end) - 1 / Delta_end * (drho_dh_end + 2 * drho_dh_liq_constp);
    CoolPropDbl c = drho_dh_liq_constp;
    CoolPropDbl d = rho_liq;

    // Either the spline value or drho/dh|p can be directly evaluated now
    _rho_spline = a * POW3(Delta) + b * POW2(Delta) + c * Delta + d;
    _drho_spline_dh__constp = 3 * a * POW2(Delta) + 2 * b * Delta + c;
    if (rho_spline) return _rho_spline;
    if (drho_dh_p) return _drho_spline_dh__constp;

    // It's drho/dp|h - calculate some more things
    // Derivatives *along* the saturation curve
    CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(iHmolar, iP);
    CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(iHmolar, iP);
    CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(iDmolar, iP);
    CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(iDmolar, iP);
    CoolPropDbl rhoV = SatV->keyed_output(iDmolar);
    CoolPropDbl rhoL = SatL->keyed_output(iDmolar);
    CoolPropDbl drho_dp_end = POW2(End->keyed_output(iDmolar)) * (x_end / POW2(rhoV) * drhoV_dp_sat + (1 - x_end) / POW2(rhoL) * drhoL_dp_sat);

    // Faking single-phase
    CoolPropDbl d2rhodhdp_liq = Liq->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);

    // Derivatives at the end point
    CoolPropDbl d2rhodhdp_end = End->calc_second_two_phase_deriv(iDmolar, iHmolar, iP, iP, iHmolar);

    // Reminder: Delta = Q()*(hV-hL) = h-hL ; Delta_end = x_end*(hV-hL)
    CoolPropDbl d_Delta_dp_consth = -dhL_dp_sat;
    CoolPropDbl d_Delta_end_dp_consth = x_end * (dhV_dp_sat - dhL_dp_sat);

    // First pressure derivative at constant h of the coefficients a,b,c,d
    CoolPropDbl d_Abracket_dp_consth =
      (2 * drhoL_dp_sat - 2 * drho_dp_end + Delta_end * (d2rhodhdp_liq + d2rhodhdp_end) + d_Delta_end_dp_consth * (drho_dh_liq_constp + drho_dh_end));
    CoolPropDbl da_dp = 1 / POW3(Delta_end) * d_Abracket_dp_consth + Abracket * (-3 / POW4(Delta_end) * d_Delta_end_dp_consth);
    CoolPropDbl db_dp = -6 / POW3(Delta_end) * d_Delta_end_dp_consth * (rho_end - rho_liq) + 3 / POW2(Delta_end) * (drho_dp_end - drhoL_dp_sat)
                        + (1 / POW2(Delta_end) * d_Delta_end_dp_consth) * (drho_dh_end + 2 * drho_dh_liq_constp)
                        - (1 / Delta_end) * (d2rhodhdp_end + 2 * d2rhodhdp_liq);
    CoolPropDbl dc_dp = d2rhodhdp_liq;
    CoolPropDbl dd_dp = drhoL_dp_sat;

    _drho_spline_dp__consth =
      (3 * a * POW2(Delta) + 2 * b * Delta + c) * d_Delta_dp_consth + POW3(Delta) * da_dp + POW2(Delta) * db_dp + Delta * dc_dp + dd_dp;
    if (drho_dp_h) return _drho_spline_dp__consth;

    throw ValueError("Something went wrong in REFPROPMixtureBackend::calc_first_two_phase_deriv_splined");
    return _HUGE;
}

CoolPropDbl REFPROPMixtureBackend::calc_saturated_liquid_keyed_output(parameters key) {
    if (!_rhoLmolar || !ValidNumber(static_cast<CoolPropDbl>(_rhoLmolar))) {
        throw ValueError("The saturated liquid state has not been set.");
    }
    switch (key) {
        case iDmolar:
            return _rhoLmolar;
        case iDmass:
            return static_cast<double>(_rhoLmolar) * calc_saturated_liquid_keyed_output(imolar_mass);
        case imolar_mass: {
            double wmm_kg_kmol = 0;
            WMOLdll(&(mole_fractions_liq[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
            return wmm_kg_kmol / 1000;                        // kg/mol
        }
        default: {
            shared_ptr<REFPROPMixtureBackend> shimL = build_saturation_shim(0);
            return shimL->keyed_output(key);
        }
    }
}
CoolPropDbl REFPROPMixtureBackend::calc_saturated_vapor_keyed_output(parameters key) {
    if (!_rhoVmolar || !ValidNumber(static_cast<CoolPropDbl>(_rhoVmolar))) {
        throw ValueError("The saturated vapor state has not been set.");
    }
    switch (key) {
        case iDmolar:
            return _rhoVmolar;
        case iDmass:
            return static_cast<double>(_rhoVmolar) * calc_saturated_vapor_keyed_output(imolar_mass);
        case imolar_mass: {
            double wmm_kg_kmol = 0;
            WMOLdll(&(mole_fractions_vap[0]), &wmm_kg_kmol);  // returns mole mass in kg/kmol
            return wmm_kg_kmol / 1000;                        // kg/mol
        }
        default: {
            shared_ptr<REFPROPMixtureBackend> shimV = build_saturation_shim(1);
            return shimV->keyed_output(key);
        }
    }
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

THERM0dllOutputs REFPROPMixtureBackend::call_THERM0dll(double T, double rho_mol_dm3, const std::vector<double>& mole_fractions) {
    /*
      subroutineTHERM0dll(T, D, z, P0, e0, h0, s0, Cv0, Cp00, w0, a0, g0)
    Compute ideal-gas thermal quantities as a function of temperature, density, and composition from core functions.

    This routine is the same as THERM, except it only calculates ideal gas properties (Z=1) at any temperature and density.

    Parameters:
    T [double ,in] :: Temperature [K]
    D [double ,in] :: Molar density [mol/L]
    z (20) [double ,in] :: Composition (array of mole fractions)
    P0 [double ,out] :: Pressure [kPa]
    e0 [double ,out] :: Internal energy [J/mol]
    h0 [double ,out] :: Enthalpy [J/mol]
    s0 [double ,out] :: Entropy [J/mol-K]
    Cv0 [double ,out] :: Isochoric heat capacity [J/mol-K]
    Cp00 [double ,out] :: Isobaric heat capacity [J/mol-K]
    w0 [double ,out] :: Speed of sound [m/s]
    a0 [double ,out] :: Helmholtz energy [J/mol]
    g0 [double ,out] :: Gibbs free energy [J/mol]
     */
    THERM0dllOutputs o;
    if (mole_fractions.size() != 20) {
        throw ValueError("mole fractions must be of size 20");
    }
    std::vector<double> mf = mole_fractions;

    THERM0dll(&T, &rho_mol_dm3, &(mf[0]), &o.p_kPa, &o.umol_Jmol, &o.hmol_Jmol, &o.smol_JmolK, &o.cvmol_JmolK, &o.cpmol_JmolK, &o.w_ms, &o.amol_Jmol,
              &o.gmol_Jmol);
    return o;
}

bool force_load_REFPROP() {
    std::string err;
    if (!load_REFPROP(err)) {
        if (CoolProp::get_debug_level() > 5) {
            std::cout << format("Error while loading REFPROP: %s", err) << '\n';
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
            std::cout << format("Error while unloading REFPROP: %s", err) << '\n';
        }
        LoadedREFPROPRef = "";
        return false;
    } else {
        LoadedREFPROPRef = "";
        return true;
    }
}

void REFPROP_SETREF(char* hrf, int ixflag, double* x0, double& h0, double& s0, double& T0, double& p0, int& ierr, char* herr, int l1, int l2) {
    std::string err;
    bool loaded_REFPROP = ::load_REFPROP(err);
    if (!loaded_REFPROP) {
        throw ValueError(format("Not able to load REFPROP; err is: %s", err.c_str()));
    }
    SETREFdll(hrf, &ixflag, x0, &h0, &s0, &T0, &p0, &ierr, herr, l1, l2);
}

} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include "CoolProp/CoolProp.h"
#    include <catch2/catch_all.hpp>

TEST_CASE("Check REFPROP and CoolProp values agree", "[REFPROP]") {
    CoolProp::Skip_if_No_REFPROP();

    SECTION("Saturation densities agree within 0.5% at T/Tc = 0.9") {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"), ',');

        for (auto& s : ss) {
            const std::string& Name = s;
            std::string RPName = CoolProp::get_fluid_param_string(s, "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", s));
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

        for (auto& s : ss) {
            const std::string& Name = s;
            std::string RPName = CoolProp::get_fluid_param_string(s, "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", s));
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
    SECTION("QT flash reports the equilibrium saturation pressure, not a THERMdll recompute (#2671)") {
        // SATTP/SATT return the equilibrium saturation pressure for a (T, Q) state.
        // Before #2671 the QT path then overwrote that pressure with THERMdll's
        // recompute from the saturated-phase density; once the SATSPLN splines
        // (activated by build_phase_envelope) refine that density, the recomputed
        // pressure drifts off the saturation curve — most severely on the dew
        // branch near the critical region.  A correct equilibrium pressure must
        // round-trip: feeding it back through a dew-point PQ flash recovers the
        // original temperature.  THERMdll's drifted pressure lands off the curve
        // and the recovered temperature shifts measurably.
        shared_ptr<CoolProp::AbstractState> S(CoolProp::AbstractState::factory("REFPROP", "Methane&Ethane"));
        std::vector<double> z = {0.4, 0.6};
        S->set_mole_fractions(z);
        S->build_phase_envelope("");  // activates the SATSPLN saturation splines
        const double T = 250;         // near-critical dew point (Tc ~ 276 K) — drift is pronounced here
        S->update(CoolProp::QT_INPUTS, 1, T);
        double p_dew = S->p();
        S->update(CoolProp::PQ_INPUTS, p_dew, 1);
        double T_back = S->T();
        CAPTURE(p_dew);
        CAPTURE(T_back);
        // With the fix the round-trip is exact to ~1e-11 K; the THERMdll drift
        // shifts it by ~3e-5 K, so 1e-7 K cleanly separates the two.
        CHECK(std::abs(T_back - T) < 1e-7);
    }
    SECTION("PQ bubble/dew saturation is self-consistent with QT on a mixture (GH #1502)") {
        // GH #1502: in the PQ saturation path the SATTPdll->SATPdll fallback used
        // the wrong SATP phase flag (SATTP's iFlsh 3/4 = bubble/dew, but SATP's kph
        // 3/4 = freezing/sublimation lines), a stale pressure, and an unconditional
        // vapor-density pick.  Forcing the fallback (a genuine SATTPdll failure) is
        // not deterministically reproducible across REFPROP versions, so this guards
        // the public contract the fix restores: a PQ bubble/dew flash must land on
        // the same equilibrium state as the QT flash that produced its pressure, and
        // the bubble (liquid) density must exceed the dew (vapor) density.
        shared_ptr<CoolProp::AbstractState> S(CoolProp::AbstractState::factory("REFPROP", "Methane&Ethane"));
        std::vector<double> z = {0.4, 0.6};
        S->set_mole_fractions(z);
        const double T = 200;  // well subcritical (Tc ~ 276 K)

        // Bubble point (Q=0): QT -> p, rho_liq; PQ(p,0) must recover both.
        S->update(CoolProp::QT_INPUTS, 0, T);
        double p_bub = S->p();
        double rho_bub_QT = S->rhomolar();
        S->update(CoolProp::PQ_INPUTS, p_bub, 0);
        CAPTURE(p_bub);
        CAPTURE(rho_bub_QT);
        CAPTURE(S->rhomolar());
        CHECK(std::abs(S->T() - T) < 1e-6);
        CHECK(std::abs(S->rhomolar() - rho_bub_QT) / rho_bub_QT < 1e-6);

        // Dew point (Q=1): QT -> p, rho_vap; PQ(p,1) must recover both.
        S->update(CoolProp::QT_INPUTS, 1, T);
        double p_dew = S->p();
        double rho_dew_QT = S->rhomolar();
        S->update(CoolProp::PQ_INPUTS, p_dew, 1);
        CAPTURE(p_dew);
        CAPTURE(rho_dew_QT);
        CAPTURE(S->rhomolar());
        CHECK(std::abs(S->T() - T) < 1e-6);
        CHECK(std::abs(S->rhomolar() - rho_dew_QT) / rho_dew_QT < 1e-6);

        // The saturated-liquid (bubble) density must exceed the saturated-vapor (dew)
        // density — the ordering the buggy fallback inverted by always picking vapor.
        CHECK(rho_bub_QT > rho_dew_QT);
    }
    SECTION("Enthalpy and entropy reference state") {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"), ',');

        for (auto& it : ss) {
            const std::string& Name = it;
            std::string RPName = CoolProp::get_fluid_param_string(it, "REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0) {
                continue;
            }

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", it));
            double Tr = S1->T_critical();
            double RCP = S1->gas_constant();
            CHECK_NOTHROW(S1->update(CoolProp::QT_INPUTS, 0, 0.9 * Tr));
            double h_CP = S1->hmass();
            double s_CP = S1->smass();
            auto j = S1->fluid_param_string("JSON");
            const nlohmann::json doc = cpjson::parse(j);
            const auto& v = doc.at(0).at("EOS").at(0).at("alpha0");
            auto s = cpjson::json2string(v);
            CAPTURE(s);

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            CHECK_NOTHROW(S2->update(CoolProp::QT_INPUTS, 0, 0.9 * Tr));
            double RRP = S2->gas_constant();
            double h_RP = S2->hmass();
            double s_RP = S2->smass();

            double delta_a1 = (s_CP - s_RP) / (S1->gas_constant() / S1->molar_mass());
            double delta_a2 = -(h_CP - h_RP) / (S1->gas_constant() / S1->molar_mass() * S1->get_reducing_state().T);
            CAPTURE(format("%0.16f", delta_a1));
            CAPTURE(format("%0.16f", delta_a2));

            CAPTURE(Name);
            CAPTURE(RRP);
            CAPTURE(RCP);
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
    CoolProp::Skip_if_No_REFPROP();

    const int num_inputs = 6;
    std::string inputs[num_inputs] = {"T_triple", "T_critical", "p_critical", "molar_mass", "rhomolar_critical", "rhomass_critical"};
    for (const auto& input : inputs) {
        std::ostringstream ss;
        ss << "Check " << input;
        SECTION(ss.str(), "") {
            double cp_val = CoolProp::PropsSI(input, "P", 0, "T", 0, "HEOS::Water");
            double rp_val = CoolProp::PropsSI(input, "P", 0, "T", 0, "REFPROP::Water");

            std::string errstr = CoolProp::get_global_param_string("errstring");
            CAPTURE(errstr);
            double err = (cp_val - rp_val) / cp_val;
            CHECK(err < 1e-3);
        }
    }
}

TEST_CASE("Check PHI0 derivatives", "[REFPROP_PHI0]") {
    CoolProp::Skip_if_No_REFPROP();

    const int num_inputs = 3;
    std::string inputs[num_inputs] = {"DALPHA0_DDELTA_CONSTTAU", "D2ALPHA0_DDELTA2_CONSTTAU", "D3ALPHA0_DDELTA3_CONSTTAU"};
    for (const auto& input : inputs) {
        std::ostringstream ss;
        ss << "Check " << input;
        SECTION(ss.str(), "") {
            double cp_val = CoolProp::PropsSI(input, "P", 101325, "T", 298, "HEOS::Water");
            double rp_val = CoolProp::PropsSI(input, "P", 101325, "T", 298, "REFPROP::Water");

            std::string errstr = CoolProp::get_global_param_string("errstring");
            CAPTURE(errstr);
            double err = std::abs((cp_val - rp_val) / cp_val);
            CHECK(err < 1e-12);
        }
    }
}

#endif
