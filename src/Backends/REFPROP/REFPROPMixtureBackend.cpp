/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

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

#include "CoolPropTools.h"
#if defined(__powerpc__)
    void *RefpropdllInstance=NULL;
    char refpropPath[] = "/opt/refprop";
#elif defined(__ISLINUX__)
    #include <dlfcn.h>
    void *RefpropdllInstance=NULL;
    char refpropPath[] = "/opt/refprop";
#elif defined(__ISWINDOWS__)
    #include <windows.h>
    HINSTANCE RefpropdllInstance=NULL;
    char refpropPath[] = "";
#elif defined(__ISAPPLE__)
    #include <dlfcn.h>
    void *RefpropdllInstance=NULL;
    char refpropPath[] = "/opt/refprop";
#else
    #pragma error
#endif

enum DLLNameManglingStyle{ NO_NAME_MANGLING = 0, LOWERCASE_NAME_MANGLING, LOWERCASE_AND_UNDERSCORE_NAME_MANGLING };

#include "REFPROP_lib.h"
#include "REFPROPMixtureBackend.h"
#include "Exceptions.h"
#include "Configuration.h"
#include "CoolProp.h"

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include "crossplatform_shared_ptr.h"

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

std::string LoadedREFPROPRef;

// Check windows
#if _WIN32 || _WIN64
   #if _WIN64
     #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

// Check GCC
#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

static bool dbg_refprop = false;
static std::string RPVersion_loaded = "";
static const unsigned int number_of_endings = 5;
std::string endings[number_of_endings] = {"", ".FLD", ".fld", ".PPF", ".ppf"};

static char rel_path_HMC_BNC[] = "HMX.BNC";
static char default_reference_state[] = "DEF";

/* Define functions as pointers and initialise them to NULL
* Declare the functions for direct access
*
* Example: SETPATHdll_POINTER SETPATHdll;
*
* ***MAGIC WARNING**!! X Macros in use
* See http://stackoverflow.com/a/148610
* See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
*/
#define X(name)  name ## _POINTER name;
 LIST_OF_REFPROP_FUNCTION_NAMES
#undef X

void *getFunctionPointer(const char * name, DLLNameManglingStyle mangling_style = NO_NAME_MANGLING)
{
    std::string function_name;
    switch(mangling_style){
        case NO_NAME_MANGLING:
            function_name = name; break;
        case LOWERCASE_NAME_MANGLING:
            function_name = lower(name); break;
        case LOWERCASE_AND_UNDERSCORE_NAME_MANGLING:
            function_name = lower(name) + "_"; break;
    }
    #if defined(__ISWINDOWS__)
        return (void *) GetProcAddress(RefpropdllInstance, function_name.c_str());
    #elif defined(__ISLINUX__)
        return dlsym(RefpropdllInstance, function_name.c_str());
    #elif defined(__ISAPPLE__)
        return dlsym(RefpropdllInstance, function_name.c_str());
    #else
        throw CoolProp::NotImplementedError("This function should not be called.");
        return NULL;
    #endif
}

//Moved pointer handling to a function, helps to maintain
//an overview and structures OS dependent parts
double setFunctionPointers()
{
    if (RefpropdllInstance==NULL)
    {
        printf("REFPROP is not loaded, make sure you call this function after loading the library.\n");
        return -_HUGE;
    }
    /* First determine the type of name mangling in use.
     * A) RPVersion -> RPVersion
     * B) RPVersion -> rpversion
     * C) RPVersion -> rpversion_
     */
     DLLNameManglingStyle mangling_style = NO_NAME_MANGLING; // defaults to no mangling

     SETUPdll = (SETUPdll_POINTER) getFunctionPointer("SETUPdll");
     if (SETUPdll == NULL){ // some mangling in use
         SETUPdll = (SETUPdll_POINTER) getFunctionPointer("setupdll");
         if (SETUPdll != NULL){
            mangling_style = LOWERCASE_NAME_MANGLING;
         }
         else{
             SETUPdll = (SETUPdll_POINTER) getFunctionPointer("setupdll_");
             if (SETUPdll != NULL){
                 mangling_style = LOWERCASE_AND_UNDERSCORE_NAME_MANGLING;
             }
             else{
                 throw CoolProp::ValueError("Could not load the symbol SETUPdll or any of its mangled forms; REFPROP shared library broken");
             }
         }
     }

    /* Set the pointers, platform independent
     *
     * Example: RPVersion = (RPVersion_POINTER) getFunctionPointer(STRINGIFY(RPVersion));
     *
     * ***MAGIC WARNING**!! X Macros in use
     * See http://stackoverflow.com/a/148610
     * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
     */
    #define X(name)  name = (name ## _POINTER) getFunctionPointer(STRINGIFY(name), mangling_style);
       LIST_OF_REFPROP_FUNCTION_NAMES
    #undef X

    return COOLPROP_OK;
}

std::string get_REFPROP_fluid_path()
{
    std::string rpPath = refpropPath;
    // Allow the user to specify an alternative REFPROP path by configuration value
    std::string alt_refprop_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_PATH);
    if (!alt_refprop_path.empty()){ rpPath = alt_refprop_path; }
    #if defined(__ISWINDOWS__)
        return rpPath;
    #elif defined(__ISLINUX__)
        return rpPath + std::string("/fluids/");
    #elif defined(__ISAPPLE__)
        return rpPath + std::string("/fluids/");
    #else
        throw CoolProp::NotImplementedError("This function should not be called.");
        return rpPath;
    #endif
}
bool load_REFPROP()
{
    // If REFPROP is not loaded
    if (RefpropdllInstance==NULL)
    {

        // Load it
        #if defined(__ISWINDOWS__)
            /* We need this logic on windows because if you use the bitness
             * macros it requires that the build bitness and the target bitness
             * are the same which is in general not the case.  Therefore, checking
             * both is safe
             */
            // First try to load the 64-bit version
            // 64-bit code here.
            TCHAR refpropdllstring[100] = TEXT("refprp64.dll");
            RefpropdllInstance = LoadLibrary(refpropdllstring);

            if (RefpropdllInstance==NULL){
                // That didn't work, let's try the 32-bit version
                // 32-bit code here.
                TCHAR refpropdllstring32[100] = TEXT("refprop.dll");
                RefpropdllInstance = LoadLibrary(refpropdllstring32);
            }

        #elif defined(__ISLINUX__)
            RefpropdllInstance = dlopen ("librefprop.so", RTLD_NOW);
        #elif defined(__ISAPPLE__)
            RefpropdllInstance = dlopen ("librefprop.dylib", RTLD_NOW);
        #else
            throw CoolProp::NotImplementedError("We should not reach this point.");
            RefpropdllInstance = NULL;
        #endif

        if (RefpropdllInstance==NULL)
        {
            #if defined(__ISWINDOWS__)
                              printf("Could not load refprop.dll \n\n");
                throw CoolProp::AttributeError("Could not load refprop.dll, make sure it is in your system search path. In case you run 64bit and you have a REFPROP license, try installing the 64bit DLL from NIST.");
            #elif defined(__ISLINUX__)
                fputs (dlerror(), stderr);
                              printf("Could not load librefprop.so \n\n");
                throw CoolProp::AttributeError("Could not load librefprop.so, make sure it is in your system search path.");
            #elif defined(__ISAPPLE__)
                fputs (dlerror(), stderr);
                              printf("Could not load librefprop.dylib \n\n");
                throw CoolProp::AttributeError("Could not load librefprop.dylib, make sure it is in your system search path.");
            #else
                throw CoolProp::NotImplementedError("Something is wrong with the platform definition, you should not end up here.");
            #endif
            return false;
        }

        if (setFunctionPointers()!=COOLPROP_OK)
        {
                          printf("There was an error setting the REFPROP function pointers, check types and names in header file.\n");
            throw CoolProp::AttributeError("There was an error setting the REFPROP function pointers, check types and names in header file.");
            return false;
        }
        char rpv[1000];
        RPVersion(rpv, 1000);
        RPVersion_loaded = rpv;
        return true;
    }
    return true;
}

namespace CoolProp {

REFPROPMixtureBackend::REFPROPMixtureBackend(const std::vector<std::string>& fluid_names) {
    // Do the REFPROP instantiation for this fluid
    _mole_fractions_set = false;

    // Try to add this fluid to REFPROP - might want to think about making array of
    // components and setting mole fractions if they change a lot.
    this->set_REFPROP_fluids(fluid_names);
    
    // Bump the number of REFPROP backends that are in existence;
    REFPROPMixtureBackend::instance_counter++;
}

REFPROPMixtureBackend::~REFPROPMixtureBackend() {
    // Decrement the counter for the number of instances
    REFPROPMixtureBackend::instance_counter--;
    // Unload the shared library when the last instance is about to be destroyed
    if (REFPROPMixtureBackend::instance_counter == 0){
        if (RefpropdllInstance!=NULL) {
            
            // Unload it
            #if defined(__ISWINDOWS__)
                FreeLibrary(RefpropdllInstance);
                //delete RefpropdllInstance;
                RefpropdllInstance = NULL;
            #elif defined(__ISLINUX__)
                dlclose (RefpropdllInstance);
                //delete RefpropdllInstance;
                RefpropdllInstance = NULL;
            #elif defined(__ISAPPLE__)
                dlclose (RefpropdllInstance);
                //delete RefpropdllInstance;
                RefpropdllInstance = NULL;
            #else
                throw CoolProp::NotImplementedError("We should not reach this point.");
                //delete RefpropdllInstance;
                RefpropdllInstance = NULL;
            #endif
            LoadedREFPROPRef = "";
        }
    }
}
void REFPROPMixtureBackend::check_loaded_fluid()
{
    this->set_REFPROP_fluids(this->fluid_names);
}

std::size_t REFPROPMixtureBackend::instance_counter = 0; // initialise with 0
bool REFPROPMixtureBackend::_REFPROP_supported = true; // initialise with true
bool REFPROPMixtureBackend::REFPROP_supported () {
    /*
     * Here we build the bridge from macro definitions
     * into the actual code. This is also going to be
     * the central place to handle error messages on
     * unsupported platforms.
     */

    // Abort check if Refprop has been loaded.
    if (RefpropdllInstance!=NULL) return true;


    // Store result of previous check.
    if (_REFPROP_supported) {
        // Either Refprop is supported or it is the first check.
        std::string rpv(STRINGIFY(RPVersion));
        if (rpv.compare("NOTAVAILABLE")!=0) {
            // Function names were defined in "REFPROP_lib.h",
            // This platform theoretically supports Refprop.
            if (load_REFPROP()) {
                return true;
            }
            else {
                printf("Good news: It is possible to use REFPROP on your system! However, the library \n");
                printf("could not be loaded. Please make sure that REFPROP is available on your system.\n\n");
                printf("Neither found in current location nor found in system PATH.\n");
                printf("If you already obtained a copy of REFPROP from http://www.nist.gov/srd/, \n");
                printf("add location of REFPROP to the PATH environment variable or your library path.\n\n");
                printf("In case you do not use Windows, have a look at https://github.com/jowr/librefprop.so \n");
                printf("to find instructions on how to compile your own version of the REFPROP library.\n\n");
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

void REFPROPMixtureBackend::set_REFPROP_fluids(const std::vector<std::string> &fluid_names)
{
    // If the name of the refrigerant doesn't match
    // that of the currently loaded refrigerant, fluids must be loaded
    if (!cached_component_string.empty() && LoadedREFPROPRef == cached_component_string)
    {
        if (CoolProp::get_debug_level() > 5){ std::cout << format("%s:%d: The current fluid can be reused; %s and %s match \n",__FILE__,__LINE__,cached_component_string.c_str(),LoadedREFPROPRef.c_str()); }
        if (dbg_refprop) std::cout << format("%s:%d: The current fluid can be reused; %s and %s match \n",__FILE__,__LINE__,cached_component_string.c_str(),LoadedREFPROPRef.c_str());
        long N = static_cast<long>(fluid_names.size());
        this->Ncomp = N;
        mole_fractions.resize(N);
        mole_fractions_liq.resize(N);
        mole_fractions_vap.resize(N);
        return;
    }
    else
    {
        long ierr=0;
        this->fluid_names = fluid_names;
        char component_string[10000], herr[errormessagelength];
        std::string components_joined = strjoin(fluid_names,"|");
        std::string components_joined_raw = strjoin(fluid_names,"|");
        std::string fdPath = get_REFPROP_fluid_path();
        long N = static_cast<long>(fluid_names.size());

        // Check platform support
        if(!REFPROP_supported()){ throw NotImplementedError("You cannot use the REFPROPMixtureBackend."); }
    
        if (N == 1 && upper(components_joined_raw).find(".MIX") != std::string::npos){
            // It's a predefined mixture
            ierr = 0;
            std::vector<double> x(ncmax);
            char mix[255];
            strcpy(mix, components_joined_raw.c_str());
            char hmx_bnc[255] = "HMX.BNC", reference_state[4] = "DEF";
            std::string alt_hmx_bnc_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_HMX_BNC_PATH);
            if (!alt_hmx_bnc_path.empty()){
                strcpy(hmx_bnc, alt_hmx_bnc_path.c_str());
            }
            SETMIXdll(mix, 
                      hmx_bnc, 
                      reference_state, 
                      &N, 
                      component_string, 
                      &(x[0]), 
                      &ierr, 
                      herr,
                      255,
                      255,
                      3,
                      10000,
                      255);
            if (static_cast<int>(ierr) <= 0){
                this->Ncomp = N;
                mole_fractions.resize(N);
                mole_fractions_liq.resize(N);
                mole_fractions_vap.resize(N);
                LoadedREFPROPRef = mix;
                cached_component_string = mix;
                if (CoolProp::get_debug_level() > 5){ std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n",__FILE__,__LINE__, components_joined.c_str()); }
                if (dbg_refprop) std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n",__FILE__,__LINE__, components_joined.c_str());
                set_mole_fractions(std::vector<CoolPropDbl>(x.begin(), x.begin()+N));
                return;
            }
            else{
                throw ValueError(format("Unable to load mixture: %s",components_joined_raw.c_str()));
            }
        }
        // Loop over the file names - first we try with nothing, then .fld, then .FLD, then .ppf - means you can't mix and match
        for (unsigned int k = 0; k < number_of_endings; k++)
        {
            // Build the mixture string
            for (unsigned int j = 0; j < (unsigned int)N; j++)
            {
                if (j == 0){
                    components_joined = fdPath + upper(fluid_names[j])+endings[k];
                }
                else{
                    components_joined += "|" + fdPath + upper(fluid_names[j])+endings[k];
                }
            }

            if (dbg_refprop) std::cout << format("%s:%d: The fluid %s has not been loaded before, current value is %s \n",__FILE__,__LINE__,components_joined_raw.c_str(),LoadedREFPROPRef.c_str());
            char path_HMX_BNC[refpropcharlength+1];
            strcpy(path_HMX_BNC, fdPath.c_str());
            strcat(path_HMX_BNC, rel_path_HMC_BNC);
            std::string alt_hmx_bnc_path = CoolProp::get_config_string(ALTERNATIVE_REFPROP_HMX_BNC_PATH);
            if (!alt_hmx_bnc_path.empty()){
                strcpy(path_HMX_BNC, alt_hmx_bnc_path.c_str());
            }
            strcpy(component_string, components_joined.c_str());

            ierr = 0;
            //...Call SETUP to initialize the program
            SETUPdll(&N, component_string, path_HMX_BNC, default_reference_state,
                     &ierr, herr,
                     10000, // Length of component_string (see PASS_FTN.for from REFPROP)
                     refpropcharlength, // Length of path_HMX_BNC
                     lengthofreference, // Length of reference
                     errormessagelength // Length of error message
                     );

            if (static_cast<int>(ierr) <= 0) // Success (or a warning, which is silently squelched for now)
            {
                this->Ncomp = N;
                mole_fractions.resize(N);
                mole_fractions_liq.resize(N);
                mole_fractions_vap.resize(N);
                LoadedREFPROPRef = component_string;
                cached_component_string = component_string;
                if (CoolProp::get_debug_level() > 5){ std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n",__FILE__,__LINE__, components_joined.c_str()); }
                if (dbg_refprop) std::cout << format("%s:%d: Successfully loaded REFPROP fluid: %s\n",__FILE__,__LINE__, components_joined.c_str());
                return;
            }
            else if (k < number_of_endings-1){ // Keep going
                if (CoolProp::get_debug_level() > 5){std::cout << format("REFPROP error/warning [ierr: %d]: %s",ierr, herr) << std::endl;}
                continue;
            }
            else
            {
                if (CoolProp::get_debug_level() > 5){std::cout << format("k: %d #endings: %d", k, number_of_endings) << std::endl;}
                throw ValueError(format("Could not load these fluids: %s", components_joined_raw.c_str()));
            }
        }
    }
}
void REFPROPMixtureBackend::set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions)
{
    if (mole_fractions.size() != this->Ncomp)
    {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]",mole_fractions.size(), this->Ncomp));
    }
    this->mole_fractions.resize(mole_fractions.size());
    for (std::size_t i = 0; i < mole_fractions.size(); ++i)
    {
        this->mole_fractions[i] = static_cast<double>(mole_fractions[i]);
    }
    this->mole_fractions_long_double = mole_fractions;
    _mole_fractions_set = true;
}
void REFPROPMixtureBackend::set_mass_fractions(const std::vector<CoolPropDbl> &mole_fractions)
{
    throw NotImplementedError("Mass fractions not currently supported");
}
void REFPROPMixtureBackend::check_status(void)
{
    if (!_mole_fractions_set){ throw ValueError("Mole fractions not yet set");}
}

void REFPROPMixtureBackend::limits(double &Tmin, double &Tmax, double &rhomolarmax, double &pmax)
{
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
    double Dmax_mol_L,pmax_kPa;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &Tmin, &Tmax, &Dmax_mol_L, &pmax_kPa, 3);
    pmax = pmax_kPa*1000;
    rhomolarmax = Dmax_mol_L*1000;
}
CoolPropDbl REFPROPMixtureBackend::calc_pmax(void){
    double Tmin, Tmax, rhomolarmax, pmax;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(pmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_Tmax(void){
    double Tmin, Tmax, rhomolarmax, pmax;
    limits(Tmin, Tmax, rhomolarmax, pmax);
    return static_cast<CoolPropDbl>(Tmax);
};
CoolPropDbl REFPROPMixtureBackend::calc_T_critical(){
    this->check_loaded_fluid();
    long ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]),&Tcrit,&pcrit_kPa,&dcrit_mol_L,&ierr,herr,255);
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); } //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(Tcrit);
};
CoolPropDbl REFPROPMixtureBackend::calc_p_critical(){
    this->check_loaded_fluid();
    long ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]),&Tcrit,&pcrit_kPa,&dcrit_mol_L,&ierr,herr,255); if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); } //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(pcrit_kPa*1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_rhomolar_critical(){
    long ierr = 0;
    char herr[255];
    double Tcrit, pcrit_kPa, dcrit_mol_L;
    CRITPdll(&(mole_fractions[0]),&Tcrit,&pcrit_kPa,&dcrit_mol_L,&ierr,herr,255); if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); } //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(dcrit_mol_L*1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_T_reducing(){
    this->check_loaded_fluid();
    double rhored_mol_L = 0, Tr = 0;
    REDXdll(&(mole_fractions[0]), &Tr, &rhored_mol_L);
    return static_cast<CoolPropDbl>(Tr);
};
CoolPropDbl REFPROPMixtureBackend::calc_rhomolar_reducing(){
    this->check_loaded_fluid();
    double rhored_mol_L = 0, Tr = 0;
    REDXdll(&(mole_fractions[0]), &Tr, &rhored_mol_L);
    return static_cast<CoolPropDbl>(rhored_mol_L*1000);
};
CoolPropDbl REFPROPMixtureBackend::calc_Ttriple(){
    this->check_loaded_fluid();
    if (mole_fractions.size() != 1){throw ValueError("calc_Ttriple cannot be evaluated for mixtures");}
    long icomp = 1;
    double wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas;
    INFOdll(&icomp, &wmm, &ttrp, &tnbpt, &tc, &pc, &Dc, &Zc, &acf, &dip, &Rgas);
    return static_cast<CoolPropDbl>(ttrp);
};
CoolPropDbl REFPROPMixtureBackend::calc_gas_constant(){
    this->check_loaded_fluid();
    double Rmix = 0;
    RMIX2dll(&(mole_fractions[0]), &Rmix);
    return static_cast<CoolPropDbl>(Rmix);
};
CoolPropDbl REFPROPMixtureBackend::calc_molar_mass(void)
{
    this->check_loaded_fluid();
    double wmm_kg_kmol;
    WMOLdll(&(mole_fractions[0]), &wmm_kg_kmol); // returns mole mass in kg/kmol
    _molar_mass = wmm_kg_kmol/1000; // kg/mol
    return static_cast<CoolPropDbl>(_molar_mass.pt());
};

double REFPROPMixtureBackend::calc_melt_Tmax()
{
    this->check_loaded_fluid();
    long ierr = 0;
    char herr[255];
    double tmin,tmax,Dmax_mol_L,pmax_kPa, Tmax_melt;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &tmin, &tmax, &Dmax_mol_L, &pmax_kPa, 3);
    // Get the maximum temperature for the melting curve by using the maximum pressure
    MELTPdll(&pmax_kPa, &(mole_fractions[0]),
             &Tmax_melt,
             &ierr,herr,errormessagelength);      // Error message
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return Tmax_melt;
}
CoolPropDbl REFPROPMixtureBackend::calc_melting_line(int param, int given, CoolPropDbl value)
{
    this->check_loaded_fluid();
    long ierr = 0;
    char herr[255];

    if (param == iP && given == iT){
        double _T = static_cast<double>(value), p_kPa;
        MELTTdll(&_T, &(mole_fractions[0]),
             &p_kPa,
             &ierr,herr,errormessagelength);      // Error message
        if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); } //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
        return p_kPa*1000;
    }
    else if (param == iT && given == iP){
        double p_kPa = static_cast<double>(value), _T;
        MELTPdll(&p_kPa, &(mole_fractions[0]),
             &_T,
             &ierr,herr,errormessagelength);      // Error message
        if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); } //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
        return p_kPa*1000;
    }
    else{
        throw ValueError(format("calc_melting_line(%s,%s,%Lg) is an invalid set of inputs ",
                                get_parameter_information(param,"short").c_str(),
                                get_parameter_information(given,"short").c_str(),
                                value
                                )
                        );
    }
}


CoolPropDbl REFPROPMixtureBackend::calc_viscosity(void)
{
    this->check_loaded_fluid();
    double eta, tcx, rhomol_L = 0.001*_rhomolar;
    long ierr = 0;
    char herr[255];
    TRNPRPdll(&_T,&rhomol_L,&(mole_fractions[0]),  // Inputs
              &eta,&tcx,                           // Outputs
              &ierr,herr,errormessagelength);      // Error message
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _viscosity = 1e-6*eta;
    _conductivity = tcx;
    return static_cast<double>(_viscosity);
}
CoolPropDbl REFPROPMixtureBackend::calc_conductivity(void)
{
    // Calling viscosity also caches conductivity, use that to save calls
    calc_viscosity();
    return static_cast<double>(_conductivity);
}
CoolPropDbl REFPROPMixtureBackend::calc_surface_tension(void)
{
    this->check_loaded_fluid();
    double sigma, rho_mol_L = 0.001*_rhomolar;
    long ierr = 0;
    char herr[255];
    SURFTdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
             &sigma,                                 // Outputs
             &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _surface_tension = sigma;
    return static_cast<double>(_surface_tension);
}
CoolPropDbl REFPROPMixtureBackend::calc_fugacity_coefficient(int i)
{
    this->check_loaded_fluid();
    double rho_mol_L = 0.001*_rhomolar;
    long ierr = 0;
    std::vector<double> fug_cof;
    fug_cof.resize(mole_fractions.size());
    char herr[255];
    FUGCOFdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
             &(fug_cof[0]),                   // Outputs
             &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); }
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<CoolPropDbl>(fug_cof[i]);
}

void REFPROPMixtureBackend::calc_phase_envelope(const std::string &type)
{
    this->check_loaded_fluid();
    long ierr = 0;
    char herr[255];
    SATSPLNdll(&(mole_fractions[0]),  // Inputs
               &ierr, herr, errormessagelength);       // Error message
    if (static_cast<int>(ierr) > 0) { throw ValueError(format("%s",herr).c_str()); }
}
CoolPropDbl REFPROPMixtureBackend::calc_cpmolar_idealgas(void)
{
    this->check_loaded_fluid();
    double rho_mol_L = 0.001*_rhomolar;
    double p0, e0, h0, s0, cv0, cp0, w0, A0, G0;
    THERM0dll(&_T,&rho_mol_L,&(mole_fractions[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
    return static_cast<CoolPropDbl>(cp0);
}

void REFPROPMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2)
{
    this->check_loaded_fluid();
    double rho_mol_L=_HUGE, rhoLmol_L=_HUGE, rhoVmol_L=_HUGE,
        hmol=_HUGE,emol=_HUGE,smol=_HUGE,cvmol=_HUGE,cpmol=_HUGE,
        w=_HUGE,q=_HUGE, mm=_HUGE, p_kPa = _HUGE;
    long ierr = 0;
    char herr[errormessagelength+1];

    clear();

    // Check that mole fractions have been set, etc.
    check_status();

    // Get the molar mass of the fluid for the given composition
    WMOLdll(&(mole_fractions[0]), &mm); // returns mole mass in kg/kmol
    _molar_mass = 0.001*mm; // [kg/mol]

    switch(input_pair)
    {
        case PT_INPUTS:
        {
            // Unit conversion for REFPROP
            p_kPa = 0.001*value1; _T = value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            TPFLSHdll(&_T,&p_kPa,&(mole_fractions[0]),&rho_mol_L,
                      &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions[0]), // Saturation terms
                      &q,&emol,&hmol,&smol,&cvmol,&cpmol,&w,
                      &ierr,herr,errormessagelength); //
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("PT: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmolarT_INPUTS:
        {
            // Unit conversion for REFPROP
            _rhomolar = value1; rho_mol_L = 0.001*value1; _T = value2; // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            TDFLSHdll(&_T,&rho_mol_L,&(mole_fractions[0]),&p_kPa,
                      &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                      &q,&emol,&hmol,&smol,&cvmol,&cpmol,&w,
                      &ierr,herr,errormessagelength);
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("DmolarT: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000;
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmassT_INPUTS:
        {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmolarT_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
        case DmolarP_INPUTS:
        {
            // Unit conversion for REFPROP
            rho_mol_L = 0.001*value1; p_kPa = 0.001*value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
            PDFLSHdll(&p_kPa,&rho_mol_L,&(mole_fractions[0]),&_T,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&hmol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength);  // Error terms
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("DmolarP: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _rhomolar = value1;
            _p = value2;
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmassP_INPUTS:
        {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmolarP_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
        case DmolarHmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            _rhomolar = value1; rho_mol_L = 0.001*value1; hmol = value2; // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DHFLSH (D,h,z,t,p,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
            DHFLSHdll(&rho_mol_L,&hmol,&(mole_fractions[0]),&_T,&p_kPa,
                      &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                      &q,&emol,&smol,&cvmol,&cpmol,&w,
                      &ierr,herr,errormessagelength);
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("DmolarHmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000;
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmassHmass_INPUTS:
        {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(DmolarHmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case DmolarSmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            _rhomolar = value1; rho_mol_L = 0.001*value1; smol = value2; // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
            DSFLSHdll(&rho_mol_L,&smol,&(mole_fractions[0]),&_T,&p_kPa,
                      &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                      &q,&emol,&hmol,&cvmol,&cpmol,&w,
                      &ierr,herr,errormessagelength);
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("DmolarSmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000;
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmassSmass_INPUTS:
        {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(DmolarSmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass );
            return;
        }
        case DmolarUmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            _rhomolar = value1; rho_mol_L = 0.001*value1; emol = value2; // Want rho in [mol/L] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            DEFLSHdll(&rho_mol_L,&emol,&(mole_fractions[0]),&_T,&p_kPa,
                      &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                      &q,&hmol,&hmol,&cvmol,&cpmol,&w,
                      &ierr,herr,errormessagelength);
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("DmolarUmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000;
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case DmassUmass_INPUTS:
        {
            // Call again, but this time with molar units
            // D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            // U: [J/mol] * [kg/mol] -> [J/mol]
            update(DmolarUmolar_INPUTS, value1 / (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case HmolarP_INPUTS:
        {
            // Unit conversion for REFPROP
            hmol = value1; p_kPa = 0.001*value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            PHFLSHdll(&p_kPa,&hmol,&(mole_fractions[0]),&_T,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms
            if (static_cast<int>(ierr) > 0) { throw ValueError(format("HmolarPmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value2;
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case HmassP_INPUTS:
        {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(HmolarP_INPUTS, value1 * (double)_molar_mass, value2);
            return;
        }
        case PSmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            p_kPa = 0.001*value1; smol = value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            PSFLSHdll(&p_kPa,&smol,&(mole_fractions[0]),&_T,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&hmol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("PSmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case PSmass_INPUTS:
        {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(PSmolar_INPUTS, value1, value2*(double)_molar_mass);
            return;
        }
        case PUmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            p_kPa = 0.001*value1; emol = value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            // from REFPROP: subroutine PEFLSH (p,e,z,t,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            PEFLSHdll(&p_kPa,&emol,&(mole_fractions[0]),&_T,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&hmol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("PUmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case PUmass_INPUTS:
        {
            // Call again, but this time with molar units
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(PUmolar_INPUTS, value1, value2*(double)_molar_mass);
            return;
        }
        case HmolarSmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            hmol = value1; smol = value2;

            HSFLSHdll(&hmol,&smol,&(mole_fractions[0]),&_T,&p_kPa,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("HmolarSmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case HmassSmass_INPUTS:
        {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol/K]
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(HmolarSmolar_INPUTS, value1 * (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case SmolarUmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            smol = value1; emol = value2;

            // from REFPROP: subroutine ESFLSH (e,s,z,t,p,D,Dl,Dv,x,y,q,h,cv,cp,w,ierr,herr)
            ESFLSHdll(&emol,&smol,&(mole_fractions[0]),&_T,&p_kPa,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("SmolarUmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case SmassUmass_INPUTS:
        {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K],
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(SmolarUmolar_INPUTS, value1 * (double)_molar_mass, value2 * (double)_molar_mass);
            return;
        }
        case SmolarT_INPUTS:
        {
            // Unit conversion for REFPROP
            smol = value1; _T = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            long kr = 1;

            // from REFPROP: subroutine TSFLSH (t,s,z,kr,p,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
            TSFLSHdll(&_T,&smol,&(mole_fractions[0]),&kr,&p_kPa,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&hmol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("SmolarT: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case SmassT_INPUTS:
        {
            // Call again, but this time with molar units
            // S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(SmolarT_INPUTS, value1 * (double)_molar_mass, value2 );
            return;
        }
        case HmolarT_INPUTS:
        {
            // Unit conversion for REFPROP
            hmol = value1; _T = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            long kr = 1;

            // from REFPROP: subroutine THFLSH (t,h,z,kr,p,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
            THFLSHdll(&_T,&hmol,&(mole_fractions[0]),&kr,&p_kPa,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&emol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("HmolarT: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case HmassT_INPUTS:
        {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(HmolarT_INPUTS, value1 * (double)_molar_mass, value2 );
            return;
        }
        case TUmolar_INPUTS:
        {
            // Unit conversion for REFPROP
            _T = value1; emol = value2;

            /*
            c  additional input--only for THFLSH, TSFLSH, and TEFLSH
            c       kr--flag specifying desired root for multi-valued inputs:
            c           1 = return lower density root
            c           2 = return higher density root
            */
            long kr = 1;

            // from REFPROP: subroutine TEFLSH (t,e,z,kr,p,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
            TEFLSHdll(&_T,&emol,&(mole_fractions[0]),&kr,&p_kPa,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &q,&hmol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("TUmolar: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case TUmass_INPUTS:
        {
            // Call again, but this time with molar units
            // U: [J/kg] * [kg/mol] -> [J/mol]
            update(TUmolar_INPUTS, value1, value2 * (double)_molar_mass);
            return;
        }
        case PQ_INPUTS:
        {
            /* From REFPROP:
            additional input--only for TQFLSH and PQFLSH
                 kq--flag specifying units for input quality
                     kq = 1 quality on MOLAR basis [moles vapor/total moles]
                     kq = 2 quality on MASS basis [mass vapor/total mass]
            */
            long kq = 1;

            // Unit conversion for REFPROP
            p_kPa = 0.001*value1; _Q = value2; // Want p in [kPa] in REFPROP

            // Use flash routine to find properties
            PQFLSHdll(&p_kPa,&_Q,&(mole_fractions[0]),&kq,&_T,&rho_mol_L,
                &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &emol,&hmol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) { throw ValueError(format("PQ: %s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

            // Set all cache values that can be set with unit conversion to SI
            _p = value1;
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        case QT_INPUTS:
        {
            /* From REFPROP:
            additional input--only for TQFLSH and PQFLSH
                 kq--flag specifying units for input quality
                     kq = 1 quality on MOLAR basis [moles vapor/total moles]
                     kq = 2 quality on MASS basis [mass vapor/total mass]
            */
            long kq = 1;

            // Unit conversion for REFPROP
            _Q = value1; _T = value2;

            // Use flash routine to find properties
            TQFLSHdll(&_T,&_Q,&(mole_fractions[0]),&kq,&p_kPa,&rho_mol_L,
                 &rhoLmol_L,&rhoVmol_L,&(mole_fractions_liq[0]),&(mole_fractions_vap[0]), // Saturation terms
                &emol,&hmol,&smol,&cvmol,&cpmol,&w, // Other thermodynamic terms
                &ierr,herr,errormessagelength); // Error terms

            if (static_cast<int>(ierr) > 0) {
                throw ValueError(format("TQ(%s): %s",LoadedREFPROPRef.c_str(), herr).c_str());
                }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());

            // Set all cache values that can be set with unit conversion to SI
            _p = p_kPa*1000; // 1000 for conversion from kPa to Pa
            _rhomolar = rho_mol_L*1000; // 1000 for conversion from mol/L to mol/m3
            if (0)
            {
                _rhoLmolar = rhoLmol_L*1000; // 1000 for conversion from mol/L to mol/m3
                _rhoVmolar = rhoVmol_L*1000; // 1000 for conversion from mol/L to mol/m3
            }
            break;
        }
        default:
        {
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
    _tau = calc_T_critical()/_T;
    _delta = _rhomolar/calc_rhomolar_critical();
    _Q = q;
}
CoolPropDbl REFPROPMixtureBackend::call_phixdll(long itau, long idel)
{
    this->check_loaded_fluid();
    double val = 0, tau = _tau, delta = _delta;
    if (PHIXdll == NULL){throw ValueError("PHIXdll function is not available in your version of REFPROP. Please upgrade");}
    PHIXdll(&itau, &idel, &tau, &delta, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val)/pow(static_cast<CoolPropDbl>(_delta),idel)/pow(static_cast<CoolPropDbl>(_tau),itau);
}
CoolPropDbl REFPROPMixtureBackend::call_phi0dll(long itau, long idel)
{
    this->check_loaded_fluid();
    throw ValueError("Temporarily the PHI0dll function is not available for REFPROP");
    double val = 0, tau = _tau, delta = _delta, __T = T(), __rho = rhomolar()/1000;
    if (PHI0dll == NULL){throw ValueError("PHI0dll function is not available in your version of REFPROP. Please upgrade");}
    PHI0dll(&itau, &idel, &__T, &__rho, &(mole_fractions[0]), &val);
    return static_cast<CoolPropDbl>(val)/pow(delta,idel)/pow(tau,itau);
}

} /* namespace CoolProp */


#ifdef ENABLE_CATCH
#include "CoolProp.h"
#include "catch.hpp"

TEST_CASE("Check REFPROP and CoolProp values agree","[REFPROP]")
{
    SECTION("Saturation densities agree within 0.5% at T/Tc = 0.9")
    {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"),',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it)
        {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it),"REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0){continue;}

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            CHECK_NOTHROW(S1->update(CoolProp::QT_INPUTS, 0, Tr*0.9););
            double rho_CP = S1->rhomolar();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            CHECK_NOTHROW(S2->update(CoolProp::QT_INPUTS, 0, Tr*0.9););
            double rho_RP = S2->rhomolar();

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(rho_CP);
            CAPTURE(rho_RP);

            double DH = (rho_RP-rho_CP)/rho_RP;
            CHECK(std::abs(DH) < 0.05);
        }
    }
    SECTION("Saturation specific heats agree within 0.5% at T/Tc = 0.9")
    {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"),',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it)
        {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it),"REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0){continue;}

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            S1->update(CoolProp::QT_INPUTS, 0, Tr*0.9);
            double cp_CP = S1->cpmolar();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            S2->update(CoolProp::QT_INPUTS, 0, Tr*0.9);
            double cp_RP = S2->cpmolar();

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(cp_CP);
            CAPTURE(cp_RP);
            CAPTURE(0.9*Tr);

            double Dcp = (cp_RP-cp_CP)/cp_RP;
            CHECK(std::abs(Dcp) < 0.05);
        }
    }
    SECTION("Enthalpy and entropy reference state")
    {
        std::vector<std::string> ss = strsplit(CoolProp::get_global_param_string("FluidsList"),',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it)
        {
            std::string Name = (*it);
            std::string RPName = CoolProp::get_fluid_param_string((*it),"REFPROP_name");

            // Skip fluids not in REFPROP
            if (RPName.find("N/A") == 0){continue;}

            shared_ptr<CoolProp::AbstractState> S1(CoolProp::AbstractState::factory("HEOS", (*it)));
            double Tr = S1->T_critical();
            CHECK_NOTHROW(S1->update(CoolProp::QT_INPUTS, 0, 0.9*Tr););
            double h_CP = S1->hmass();
            double s_CP = S1->smass();

            shared_ptr<CoolProp::AbstractState> S2(CoolProp::AbstractState::factory("REFPROP", RPName));
            CHECK_NOTHROW(S2->update(CoolProp::QT_INPUTS, 0, 0.9*Tr););
            double h_RP = S2->hmass();
            double s_RP = S2->smass();

            CAPTURE(Name);
            CAPTURE(RPName);
            CAPTURE(h_CP);
            CAPTURE(h_RP);
            CAPTURE(s_CP);
            CAPTURE(s_RP);
            double DH = (S1->hmass()-S2->hmass());
            double DS = (S1->smass()-S2->smass());

            CHECK(std::abs(DH/h_RP) < 0.01);
            CHECK(std::abs(DS/s_RP) < 0.01);
        }
    }
}

TEST_CASE("Check trivial inputs for REFPROP work", "[REFPROP_trivial]")
{
    const int num_inputs = 6;
    std::string inputs[num_inputs] = {"T_triple", "T_critical", "p_critical", "molar_mass", "rhomolar_critical", "rhomass_critical"};
    for (int i = 0; i < num_inputs; ++i){
        std::ostringstream ss;
        ss << "Check " << inputs[i];
        SECTION(ss.str(),"")
        {
            double cp_val = CoolProp::PropsSI(inputs[i],"P",0,"T",0,"HEOS::Water");
            double rp_val = CoolProp::PropsSI(inputs[i],"P",0,"T",0,"REFPROP::Water");
            
            std::string errstr = CoolProp::get_global_param_string("errstring");
            CAPTURE(errstr);
            double err = (cp_val - rp_val)/cp_val;
            CHECK(err < 1e-3);
        }
    }
}

#endif
