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
#if defined(__ISWINDOWS__)
    #include <windows.h>
    HINSTANCE RefpropdllInstance=NULL;
    char refpropPath[] = "";
#elif defined(__ISLINUX__)
    #include <dlfcn.h>
    void *RefpropdllInstance=NULL;
    char refpropPath[] = "/opt/refprop";
#elif defined(__ISAPPLE__)
    #include <dlfcn.h>
    void *RefpropdllInstance=NULL;
    char refpropPath[] = "/opt/refprop";
#else
    #pragma error
#endif

#include "REFPROP_lib.h"
#include "REFPROPMixtureBackend.h"
#include "Exceptions.h"

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <assert.h>



#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

std::string LoadedREFPROPRef;

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

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

static char rel_path_HMC_BNC[] = "HMX.BNC";
static char default_reference_state[] = "DEF";

// Define functions as pointers and initialise them to NULL
// Declare the functions for direct access
 RPVersion_POINTER RPVersion;
 SETPATHdll_POINTER SETPATHdll;
 ABFL1dll_POINTER ABFL1dll;
 ABFL2dll_POINTER ABFL2dll;
 ACTVYdll_POINTER ACTVYdll;
 AGdll_POINTER AGdll;
 CCRITdll_POINTER CCRITdll;
 CP0dll_POINTER CP0dll;
 CRITPdll_POINTER CRITPdll;
 CSATKdll_POINTER CSATKdll;
 CV2PKdll_POINTER CV2PKdll;
 CVCPKdll_POINTER CVCPKdll;
 CVCPdll_POINTER CVCPdll;
 DBDTdll_POINTER DBDTdll;
 DBFL1dll_POINTER DBFL1dll;
 DBFL2dll_POINTER DBFL2dll;
 DDDPdll_POINTER DDDPdll;
 DDDTdll_POINTER DDDTdll;
 DEFLSHdll_POINTER DEFLSHdll;
 DHD1dll_POINTER DHD1dll;
 DHFLSHdll_POINTER DHFLSHdll;
 DHFL1dll_POINTER DHFL1dll;
 DHFL2dll_POINTER DHFL2dll;
 DIELECdll_POINTER DIELECdll;
 DOTFILLdll_POINTER DOTFILLdll;
 DPDD2dll_POINTER DPDD2dll;
 DPDDKdll_POINTER DPDDKdll;
 DPDDdll_POINTER DPDDdll;
 DPDTKdll_POINTER DPDTKdll;
 DPDTdll_POINTER DPDTdll;
 DPTSATKdll_POINTER DPTSATKdll;
 DSFLSHdll_POINTER DSFLSHdll;
 DSFL1dll_POINTER DSFL1dll;
 DSFL2dll_POINTER DSFL2dll;
 ENTHALdll_POINTER ENTHALdll;
 ENTROdll_POINTER ENTROdll;
 ESFLSHdll_POINTER ESFLSHdll;
 FGCTYdll_POINTER FGCTYdll;
 FPVdll_POINTER FPVdll;
 FUGCOFdll_POINTER FUGCOFdll;
 GERG04dll_POINTER GERG04dll;
 GETFIJdll_POINTER GETFIJdll;
 GETKTVdll_POINTER GETKTVdll;
 GIBBSdll_POINTER GIBBSdll;
 HSFLSHdll_POINTER HSFLSHdll;
 INFOdll_POINTER INFOdll;
 LIMITKdll_POINTER LIMITKdll;
 LIMITSdll_POINTER LIMITSdll;
 LIMITXdll_POINTER LIMITXdll;
 MELTPdll_POINTER MELTPdll;
 MELTTdll_POINTER MELTTdll;
 MLTH2Odll_POINTER MLTH2Odll;
 NAMEdll_POINTER NAMEdll;
 PDFL1dll_POINTER PDFL1dll;
 PDFLSHdll_POINTER PDFLSHdll;
 PEFLSHdll_POINTER PEFLSHdll;
 PHFL1dll_POINTER PHFL1dll;
 PHFLSHdll_POINTER PHFLSHdll;
 PQFLSHdll_POINTER PQFLSHdll;
 PREOSdll_POINTER PREOSdll;
 PRESSdll_POINTER PRESSdll;
 PSFL1dll_POINTER PSFL1dll;
 PSFLSHdll_POINTER PSFLSHdll;
 PUREFLDdll_POINTER PUREFLDdll;
 QMASSdll_POINTER QMASSdll;
 QMOLEdll_POINTER QMOLEdll;
 RESIDUALdll_POINTER RESIDUALdll;
 SATDdll_POINTER SATDdll;
 SATEdll_POINTER SATEdll;
 SATHdll_POINTER SATHdll;
 SATPdll_POINTER SATPdll;
 SATSdll_POINTER SATSdll;
 SATTdll_POINTER SATTdll;
 SETAGAdll_POINTER SETAGAdll;
 SETKTVdll_POINTER SETKTVdll;
 SETMIXdll_POINTER SETMIXdll;
 SETMODdll_POINTER SETMODdll;
 SETREFdll_POINTER SETREFdll;
 SETUPdll_POINTER SETUPdll;
//  SPECGRdll_POINTER SPECGRdll; // not found in library
 SUBLPdll_POINTER SUBLPdll;
 SUBLTdll_POINTER SUBLTdll;
 SURFTdll_POINTER SURFTdll;
 SURTENdll_POINTER SURTENdll;
 TDFLSHdll_POINTER TDFLSHdll;
 TEFLSHdll_POINTER TEFLSHdll;
 THERM0dll_POINTER THERM0dll;
 THERM2dll_POINTER THERM2dll;
 THERM3dll_POINTER THERM3dll;
 THERMdll_POINTER THERMdll;
 THFLSHdll_POINTER THFLSHdll;
 TPFLSHdll_POINTER TPFLSHdll;
 TPFL2dll_POINTER TPFL2dll;
 TPRHOdll_POINTER TPRHOdll;
 TQFLSHdll_POINTER TQFLSHdll;
 TRNPRPdll_POINTER TRNPRPdll;
 TSFLSHdll_POINTER TSFLSHdll;
 VIRBdll_POINTER VIRBdll;
 VIRCdll_POINTER VIRCdll;
 WMOLdll_POINTER WMOLdll;
 XMASSdll_POINTER XMASSdll;
 XMOLEdll_POINTER XMOLEdll;

void *getFunctionPointer(char * name)
{
    #if defined(__ISWINDOWS__)
        return (void *) GetProcAddress(RefpropdllInstance,name);
    #elif defined(__ISLINUX__)
        return dlsym(RefpropdllInstance,name);
    #elif defined(__ISAPPLE__)
        return dlsym(RefpropdllInstance,name);
    #else
        throw CoolProp::NotImplementedError("This function should not be called.");
        return NULL;
    #endif
}

//#include <dlfcn.h>
//void *RefpropdllInstance=NULL;
//char refpropPath[] = "/opt/refprop";

//Moved pointer handling to a function, helps to maintain
//an overview and structures OS dependent parts
double setFunctionPointers()
{
    if (RefpropdllInstance==NULL)
    {
        printf("REFPROP is not loaded, make sure you call this function after loading the library.\n");
        return -_HUGE;
    }
    // set the pointers, platform independent
    RPVersion = (RPVersion_POINTER) getFunctionPointer((char *)RPVersion_NAME);
    ABFL1dll = (ABFL1dll_POINTER) getFunctionPointer((char *)ABFL1dll_NAME);
    ABFL2dll = (ABFL2dll_POINTER) getFunctionPointer((char *)ABFL2dll_NAME);
    ACTVYdll = (ACTVYdll_POINTER) getFunctionPointer((char *)ACTVYdll_NAME);
    AGdll = (AGdll_POINTER) getFunctionPointer((char *)AGdll_NAME);
    CCRITdll = (CCRITdll_POINTER) getFunctionPointer((char *)CCRITdll_NAME);
    CP0dll = (CP0dll_POINTER) getFunctionPointer((char *)CP0dll_NAME);
    CRITPdll = (CRITPdll_POINTER) getFunctionPointer((char *)CRITPdll_NAME);
    CSATKdll = (CSATKdll_POINTER) getFunctionPointer((char *)CSATKdll_NAME);
    CV2PKdll = (CV2PKdll_POINTER) getFunctionPointer((char *)CV2PKdll_NAME);
    CVCPKdll = (CVCPKdll_POINTER) getFunctionPointer((char *)CVCPKdll_NAME);
    CVCPdll = (CVCPdll_POINTER) getFunctionPointer((char *)CVCPdll_NAME);
    DBDTdll = (DBDTdll_POINTER) getFunctionPointer((char *)DBDTdll_NAME);
    DBFL1dll = (DBFL1dll_POINTER) getFunctionPointer((char *)DBFL1dll_NAME);
    DBFL2dll = (DBFL2dll_POINTER) getFunctionPointer((char *)DBFL2dll_NAME);
    DDDPdll = (DDDPdll_POINTER) getFunctionPointer((char *)DDDPdll_NAME);
    DDDTdll = (DDDTdll_POINTER) getFunctionPointer((char *)DDDTdll_NAME);
    DEFLSHdll = (DEFLSHdll_POINTER) getFunctionPointer((char *)DEFLSHdll_NAME);
    DHD1dll = (DHD1dll_POINTER) getFunctionPointer((char *)DHD1dll_NAME);
    DHFLSHdll = (DHFLSHdll_POINTER) getFunctionPointer((char *)DHFLSHdll_NAME);
    DIELECdll = (DIELECdll_POINTER) getFunctionPointer((char *)DIELECdll_NAME);
    DOTFILLdll = (DOTFILLdll_POINTER) getFunctionPointer((char *)DOTFILLdll_NAME);
    DPDD2dll = (DPDD2dll_POINTER) getFunctionPointer((char *)DPDD2dll_NAME);
    DPDDKdll = (DPDDKdll_POINTER) getFunctionPointer((char *)DPDDKdll_NAME);
    DPDDdll = (DPDDdll_POINTER) getFunctionPointer((char *)DPDDdll_NAME);
    DPDTKdll = (DPDTKdll_POINTER) getFunctionPointer((char *)DPDTKdll_NAME);
    DPDTdll = (DPDTdll_POINTER) getFunctionPointer((char *)DPDTdll_NAME);
    DPTSATKdll = (DPTSATKdll_POINTER) getFunctionPointer((char *)DPTSATKdll_NAME);
    DSFLSHdll = (DSFLSHdll_POINTER) getFunctionPointer((char *)DSFLSHdll_NAME);
    ENTHALdll = (ENTHALdll_POINTER) getFunctionPointer((char *)ENTHALdll_NAME);
    ENTROdll = (ENTROdll_POINTER) getFunctionPointer((char *)ENTROdll_NAME);
    ESFLSHdll = (ESFLSHdll_POINTER) getFunctionPointer((char *)ESFLSHdll_NAME);
    FGCTYdll = (FGCTYdll_POINTER) getFunctionPointer((char *)FGCTYdll_NAME);
    FPVdll = (FPVdll_POINTER) getFunctionPointer((char *)FPVdll_NAME);
    FUGCOFdll = (FUGCOFdll_POINTER) getFunctionPointer((char *)FUGCOFdll_NAME);
    GERG04dll = (GERG04dll_POINTER) getFunctionPointer((char *)GERG04dll_NAME);
    GETFIJdll = (GETFIJdll_POINTER) getFunctionPointer((char *)GETFIJdll_NAME);
    GETKTVdll = (GETKTVdll_POINTER) getFunctionPointer((char *)GETKTVdll_NAME);
    GIBBSdll = (GIBBSdll_POINTER) getFunctionPointer((char *)GIBBSdll_NAME);
    HSFLSHdll = (HSFLSHdll_POINTER) getFunctionPointer((char *)HSFLSHdll_NAME);
    INFOdll = (INFOdll_POINTER) getFunctionPointer((char *)INFOdll_NAME);
    LIMITKdll = (LIMITKdll_POINTER) getFunctionPointer((char *)LIMITKdll_NAME);
    LIMITSdll = (LIMITSdll_POINTER) getFunctionPointer((char *)LIMITSdll_NAME);
    LIMITXdll = (LIMITXdll_POINTER) getFunctionPointer((char *)LIMITXdll_NAME);
    MELTPdll = (MELTPdll_POINTER) getFunctionPointer((char *)MELTPdll_NAME);
    MELTTdll = (MELTTdll_POINTER) getFunctionPointer((char *)MELTTdll_NAME);
    MLTH2Odll = (MLTH2Odll_POINTER) getFunctionPointer((char *)MLTH2Odll_NAME);
    NAMEdll = (NAMEdll_POINTER) getFunctionPointer((char *)NAMEdll_NAME);
    PDFL1dll = (PDFL1dll_POINTER) getFunctionPointer((char *)PDFL1dll_NAME);
    PDFLSHdll = (PDFLSHdll_POINTER) getFunctionPointer((char *)PDFLSHdll_NAME);
    PEFLSHdll = (PEFLSHdll_POINTER) getFunctionPointer((char *)PEFLSHdll_NAME);
    PHFL1dll = (PHFL1dll_POINTER) getFunctionPointer((char *)PHFL1dll_NAME);
    PHFLSHdll = (PHFLSHdll_POINTER) getFunctionPointer((char *)PHFLSHdll_NAME);
    PQFLSHdll = (PQFLSHdll_POINTER) getFunctionPointer((char *)PQFLSHdll_NAME);
    PREOSdll = (PREOSdll_POINTER) getFunctionPointer((char *)PREOSdll_NAME);
    PRESSdll = (PRESSdll_POINTER) getFunctionPointer((char *)PRESSdll_NAME);
    PSFL1dll = (PSFL1dll_POINTER) getFunctionPointer((char *)PSFL1dll_NAME);
    PSFLSHdll = (PSFLSHdll_POINTER) getFunctionPointer((char *)PSFLSHdll_NAME);
    PUREFLDdll = (PUREFLDdll_POINTER) getFunctionPointer((char *)PUREFLDdll_NAME);
    RESIDUALdll = (RESIDUALdll_POINTER) getFunctionPointer((char *)RESIDUALdll_NAME);
    QMASSdll = (QMASSdll_POINTER) getFunctionPointer((char *)QMASSdll_NAME);
    QMOLEdll = (QMOLEdll_POINTER) getFunctionPointer((char *)QMOLEdll_NAME);
    SATDdll = (SATDdll_POINTER) getFunctionPointer((char *)SATDdll_NAME);
    SATEdll = (SATEdll_POINTER) getFunctionPointer((char *)SATEdll_NAME);
    SATHdll = (SATHdll_POINTER) getFunctionPointer((char *)SATHdll_NAME);
    SATPdll = (SATPdll_POINTER) getFunctionPointer((char *)SATPdll_NAME);
    SATSdll = (SATSdll_POINTER) getFunctionPointer((char *)SATSdll_NAME);
    SATTdll = (SATTdll_POINTER) getFunctionPointer((char *)SATTdll_NAME);
    SETAGAdll = (SETAGAdll_POINTER) getFunctionPointer((char *)SETAGAdll_NAME);
    SETKTVdll = (SETKTVdll_POINTER) getFunctionPointer((char *)SETKTVdll_NAME);
    SETMIXdll = (SETMIXdll_POINTER) getFunctionPointer((char *)SETMIXdll_NAME);
    SETMODdll = (SETMODdll_POINTER) getFunctionPointer((char *)SETMODdll_NAME);
    SETREFdll = (SETREFdll_POINTER) getFunctionPointer((char *)SETREFdll_NAME);
    SETUPdll = (SETUPdll_POINTER) getFunctionPointer((char *)SETUPdll_NAME);
//		SPECGRdll = (SPECGRdll_POINTER) getFunctionPointer((char *)SPECGRdll_NAME); // not in library
    SUBLPdll = (SUBLPdll_POINTER) getFunctionPointer((char *)SUBLPdll_NAME);
    SUBLTdll = (SUBLTdll_POINTER) getFunctionPointer((char *)SUBLTdll_NAME);
    SURFTdll = (SURFTdll_POINTER) getFunctionPointer((char *)SURFTdll_NAME);
    SURTENdll = (SURTENdll_POINTER) getFunctionPointer((char *)SURTENdll_NAME);
    TDFLSHdll = (TDFLSHdll_POINTER) getFunctionPointer((char *)TDFLSHdll_NAME);
    TEFLSHdll = (TEFLSHdll_POINTER) getFunctionPointer((char *)TEFLSHdll_NAME);
    THERM0dll = (THERM0dll_POINTER) getFunctionPointer((char *)THERM0dll_NAME);
    THERM2dll = (THERM2dll_POINTER) getFunctionPointer((char *)THERM2dll_NAME);
    THERM3dll = (THERM3dll_POINTER) getFunctionPointer((char *)THERM3dll_NAME);
    THERMdll = (THERMdll_POINTER) getFunctionPointer((char *)THERMdll_NAME);
    THFLSHdll = (THFLSHdll_POINTER) getFunctionPointer((char *)THFLSHdll_NAME);
    TPFLSHdll = (TPFLSHdll_POINTER) getFunctionPointer((char *)TPFLSHdll_NAME);
    TPRHOdll = (TPRHOdll_POINTER) getFunctionPointer((char *)TPRHOdll_NAME);
    TQFLSHdll = (TQFLSHdll_POINTER) getFunctionPointer((char *)TQFLSHdll_NAME);
    TRNPRPdll = (TRNPRPdll_POINTER) getFunctionPointer((char *)TRNPRPdll_NAME);
    TSFLSHdll = (TSFLSHdll_POINTER) getFunctionPointer((char *)TSFLSHdll_NAME);
    VIRBdll = (VIRBdll_POINTER) getFunctionPointer((char *)VIRBdll_NAME);
    VIRCdll = (VIRCdll_POINTER) getFunctionPointer((char *)VIRCdll_NAME);
    WMOLdll = (WMOLdll_POINTER) getFunctionPointer((char *)WMOLdll_NAME);
    XMASSdll = (XMASSdll_POINTER) getFunctionPointer((char *)XMASSdll_NAME);
    XMOLEdll = (XMOLEdll_POINTER) getFunctionPointer((char *)XMOLEdll_NAME);
    return COOLPROP_OK;
}

std::string get_REFPROP_fluid_path()
{
    std::string rpPath = refpropPath;
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
            #if defined(ENV64BIT)
                // 64-bit code here.
                TCHAR refpropdllstring[100] = TEXT("refprp64.dll");
                RefpropdllInstance = LoadLibrary(refpropdllstring);
            #elif defined (ENV32BIT)
                // 32-bit code here.
                TCHAR refpropdllstring[100] = TEXT("refprop.dll");
                RefpropdllInstance = LoadLibrary(refpropdllstring);
            #else
                // INCREASE ROBUSTNESS. ALWAYS THROW AN ERROR ON THE ELSE.
                #error "Must define either ENV32BIT or ENV64BIT"
            #endif
            
        #elif defined(__ISLINUX__)
            RefpropdllInstance = dlopen ("librefprop.so", RTLD_LAZY);
        #elif defined(__ISAPPLE__)
            RefpropdllInstance = dlopen ("librefprop.dylib", RTLD_LAZY);
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

        #if defined(__ISWINDOWS__)
        
        // Get data associated with path using the windows libraries, 
        // and if you can (result == 0), the path exists
        #ifdef __MINGW32__
            struct stat buf;
            if ( stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
                throw CoolProp::ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
            }
        #else
            struct _stat buf;
            if ( _stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
                throw CoolProp::ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
            }
        #endif
        #endif

        if (setFunctionPointers()!=COOLPROP_OK)
        {
                          printf("There was an error setting the REFPROP function pointers, check types and names in header file.\n");
            throw CoolProp::AttributeError("There was an error setting the REFPROP function pointers, check types and names in header file.");
            return false;
        }
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

}

REFPROPMixtureBackend::~REFPROPMixtureBackend() {
    // TODO Auto-generated destructor stub
}

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
        std::string rpv(RPVersion_NAME);
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
    long ierr=0;
    char component_string[10000], herr[errormessagelength];
    std::string components_joined = strjoin(fluid_names,"|");
    std::string fdPath = get_REFPROP_fluid_path();
    long N = static_cast<long>(fluid_names.size());

    assert(N < 20);

    // Check platform support
    if(!REFPROP_supported()){
        throw NotImplementedError("You cannot use the REFPROPMixtureBackend.");
    }

    // Build the mixture string
    for (unsigned int j = 0; j < (unsigned int)N; j++)
    {
        if (j == 0){
            components_joined = fdPath + fluid_names[j]+".fld";
        }
        else{
            components_joined += "|" + fdPath + fluid_names[j]+".fld";
        }
    }

    // Load REFPROP if it isn't loaded yet
    load_REFPROP();
    
    // If the name of the refrigerant doesn't match 
    // that of the currently loaded refrigerant
    if (LoadedREFPROPRef.compare(components_joined))
    {
        char path_HMX_BNC[refpropcharlength];
        //strcpy(path_HMX_BNC,fdPath.c_str());
        strcpy(path_HMX_BNC, rel_path_HMC_BNC);
        strcpy(component_string, components_joined.c_str());

        //...Call SETUP to initialize the program
        SETUPdll(&N, component_string, path_HMX_BNC, default_reference_state,
                 &ierr, herr,
                 10000, // Length of component_string (see PASS_FTN.for from REFPROP)
                 refpropcharlength, // Length of path_HMX_BNC
                 lengthofreference, // Length of reference
                 errormessagelength // Length of error message
                 );

        if (ierr == 0) // Success
        {
            Ncomp = N;
            mole_fractions.resize(N);
            mole_fractions_liq.resize(N);
            mole_fractions_vap.resize(N);
        }
        else if (ierr > 0) // Error
        {
            throw ValueError(format("%s",herr));
        }
        else // Warning
        {
            throw ValueError(format("%s",herr));
        }
    }
}
void REFPROPMixtureBackend::set_mole_fractions(const std::vector<long double> &mole_fractions)
{
    if (mole_fractions.size() != Ncomp)
    {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]",mole_fractions.size(), Ncomp));
    }
    this->mole_fractions.resize(mole_fractions.size());
    for (std::size_t i = 0; i < mole_fractions.size(); ++i)
    {
        this->mole_fractions[i] = static_cast<double>(mole_fractions[i]);
    }
    _mole_fractions_set = true;
}
void REFPROPMixtureBackend::set_mass_fractions(const std::vector<long double> &mole_fractions)
{
    throw NotImplementedError("Mass fractions not currently supported");
}
void REFPROPMixtureBackend::check_status(void)
{
    if (!_mole_fractions_set){ throw ValueError("Mole fractions not yet set");}
}
double REFPROPMixtureBackend::calc_melt_Tmax()
{
    long ierr;
    char herr[255];
    double tmin,tmax,Dmax_mol_L,pmax_kPa, Tmax_melt;
    char htyp[] = "EOS";
    LIMITSdll(htyp, &(mole_fractions[0]), &tmin, &tmax, &Dmax_mol_L, &pmax_kPa, 3);
    // Get the maximum temperature for the melting curve by using the maximum pressure
    MELTPdll(&pmax_kPa, &(mole_fractions[0]),
             &Tmax_melt,
             &ierr,herr,errormessagelength);      // Error message
    if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } 
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return Tmax_melt;
}
long double REFPROPMixtureBackend::calc_melt_p_T(long double T)
{
    double _T = static_cast<double>(T), p_kPa;
    long ierr;
    char herr[255];

    if (T > calc_melt_Tmax())
    {
        throw ValueError(format("Melting temperature [%g] is out of range",T));
    }

    MELTTdll(&_T, &(mole_fractions[0]),
             &p_kPa,
             &ierr,herr,errormessagelength);      // Error message
    if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } 
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return p_kPa*1000;
}
long double REFPROPMixtureBackend::calc_melt_T_p(long double p)
{
    throw NotImplementedError();
}
long double REFPROPMixtureBackend::calc_melt_rho_T(long double T)
{
    throw NotImplementedError();
}

long double REFPROPMixtureBackend::calc_viscosity(void)
{
    double eta, tcx, rhomol_L = 0.001*_rhomolar;
    long ierr;
    char herr[255];
    TRNPRPdll(&_T,&rhomol_L,&(mole_fractions[0]),  // Inputs
              &eta,&tcx,                           // Outputs
              &ierr,herr,errormessagelength);      // Error message
    if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } 
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _viscosity = 1e-6*eta;
    _conductivity = tcx;
    return static_cast<double>(_viscosity);
}
long double REFPROPMixtureBackend::calc_conductivity(void)
{
    // Calling viscosity also caches conductivity, use that to save calls
    calc_viscosity();
    return static_cast<double>(_conductivity);
}
long double REFPROPMixtureBackend::calc_surface_tension(void)
{
    double sigma, rho_mol_L = 0.001*_rhomolar;
    long ierr;
    char herr[255];
    SURFTdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
             &sigma,                                 // Outputs
             &ierr, herr, errormessagelength);       // Error message
    if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } 
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    _surface_tension = sigma;
    return static_cast<double>(_surface_tension);
}
long double REFPROPMixtureBackend::calc_fugacity_coefficient(int i)
{
    double rho_mol_L = 0.001*_rhomolar;
    long ierr;
    std::vector<double> fug_cof;
    fug_cof.resize(mole_fractions.size());
    char herr[255];
    FUGCOFdll(&_T, &rho_mol_L, &(mole_fractions[0]),  // Inputs
             &(fug_cof[0]),                   // Outputs
             &ierr, herr, errormessagelength);       // Error message
    if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } 
    //else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
    return static_cast<long double>(fug_cof[i]);
}
    
void REFPROPMixtureBackend::update(long input_pair, double value1, double value2)
{
    double rho_mol_L=_HUGE, rhoLmol_L=_HUGE, rhoVmol_L=_HUGE,
        hmol=_HUGE,emol=_HUGE,smol=_HUGE,cvmol=_HUGE,cpmol=_HUGE,
        w=_HUGE,q=_HUGE, mm=_HUGE, p_kPa = _HUGE;
    long ierr;
    char herr[255];

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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            
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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            
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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            
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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
            
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
            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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

            if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); }// TODO: else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

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
}

} /* namespace CoolProp */
