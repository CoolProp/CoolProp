
#ifndef REFPROP_LIB_H
#define REFPROP_LIB_H

/*
// The idea here is to have a common header for Windows
// and gcc-like systems. The Windows branch should cover the
// functions provided by the .dll and the gcc part covers
// the compiled .so/.dym file. Name changes caused by gfortran
// are respected and should be accounted for.
*/
// Get the platform identifiers
#include "CoolPropTools.h"

// Do some manual changes to the function names
// if needed, uses CoolProp platform detection.
#if defined(__ISWINDOWS__)
// Define compiler specific calling conventions
// for the shared library.
#  define CALLCONV __stdcall // __declspec(dllexport)
// Do not redefine function names for the shared library,
// in this case it is the REFPROP.dll and no special
// names are needed. Macros still need a value for the
// name function used below.
#  define RPVersion RPVersion
#  define SETPATHdll SETPATHdll
#  define ABFL1dll ABFL1dll
#  define ABFL2dll ABFL2dll
#  define ACTVYdll ACTVYdll
#  define AGdll AGdll
#  define CCRITdll CCRITdll
#  define CP0dll CP0dll
#  define CRITPdll CRITPdll
#  define CSATKdll CSATKdll
#  define CV2PKdll CV2PKdll
#  define CVCPKdll CVCPKdll
#  define CVCPdll CVCPdll
#  define DBDTdll DBDTdll
#  define DBFL1dll DBFL1dll
#  define DBFL2dll DBFL2dll
#  define DDDPdll DDDPdll
#  define DDDTdll DDDTdll
#  define DEFLSHdll DEFLSHdll
#  define DHD1dll DHD1dll
#  define DHFL1dll DHFL1dll
#  define DHFL2dll DHFL2dll
#  define DHFLSHdll DHFLSHdll
#  define DIELECdll DIELECdll
#  define DOTFILLdll DOTFILLdll
#  define DPDD2dll DPDD2dll
#  define DPDDKdll DPDDKdll
#  define DPDDdll DPDDdll
#  define DPDTKdll DPDTKdll
#  define DPDTdll DPDTdll
#  define DPTSATKdll DPTSATKdll
#  define DSFLSHdll DSFLSHdll
#  define DSFL1dll DSFL1dll
#  define DSFL2dll DSFL2dll
#  define ENTHALdll ENTHALdll
#  define ENTROdll ENTROdll
#  define ESFLSHdll ESFLSHdll
#  define FGCTYdll FGCTYdll
#  define FUGCOFdll FUGCOFdll
#  define FPVdll FPVdll
#  define GERG04dll GERG04dll
#  define GETFIJdll GETFIJdll
#  define GETKTVdll GETKTVdll
#  define GIBBSdll GIBBSdll
#  define HSFLSHdll HSFLSHdll
#  define INFOdll INFOdll
#  define LIMITKdll LIMITKdll
#  define LIMITSdll LIMITSdll
#  define LIMITXdll LIMITXdll
#  define MELTPdll MELTPdll
#  define MELTTdll MELTTdll
#  define MLTH2Odll MLTH2Odll
#  define NAMEdll NAMEdll
#  define PASSCMNdll PASSCMN
#  define PDFL1dll PDFL1dll
#  define PDFLSHdll PDFLSHdll
#  define PEFLSHdll PEFLSHdll
#  define PHFL1dll PHFL1dll
#  define PHFLSHdll PHFLSHdll
#  define PHIXdll PHIXdll
#  define PHI0dll PHI0dll
#  define PQFLSHdll PQFLSHdll
#  define PREOSdll PREOSdll
#  define PRESSdll PRESSdll
#  define PSFL1dll PSFL1dll
#  define PSFLSHdll PSFLSHdll
#  define PUREFLDdll PUREFLDdll
#  define QMASSdll QMASSdll
#  define QMOLEdll QMOLEdll
#  define RESIDUALdll RESIDUALdll
#  define REDXdll REDXdll
#  define RMIX2dll RMIX2dll
#  define SATDdll SATDdll
#  define SATEdll SATEdll
#  define SATHdll SATHdll
#  define SATPdll SATPdll
#  define SATSdll SATSdll
#  define SATTdll SATTdll
#  define SATSPLNdll SATSPLNdll
#  define SETAGAdll SETAGAdll
#  define SETKTVdll SETKTVdll
#  define SETMIXdll SETMIXdll
#  define SETMODdll SETMODdll
#  define SETREFdll SETREFdll
#  define SETUPdll SETUPdll
//#  define SPECGRdll SPECGRdll // not found in library
#  define SUBLPdll SUBLPdll
#  define SUBLTdll SUBLTdll
#  define SURFTdll SURFTdll
#  define SURTENdll SURTENdll
#  define TDFLSHdll TDFLSHdll
#  define TEFLSHdll TEFLSHdll
#  define THERM0dll THERM0dll
#  define THERM2dll THERM2dll
#  define THERM3dll THERM3dll
#  define THERMdll THERMdll
#  define THFLSHdll THFLSHdll
#  define TPFLSHdll TPFLSHdll
#  define TPFL2dll TPFL2dll
#  define TPRHOdll TPRHOdll
#  define TQFLSHdll TQFLSHdll
#  define TRNPRPdll TRNPRPdll
#  define TSFLSHdll TSFLSHdll
#  define VIRBdll VIRBdll
#  define VIRCdll VIRCdll
#  define WMOLdll WMOLdll
#  define XMASSdll XMASSdll
#  define XMOLEdll XMOLEdll
#elif defined(__ISLINUX__) // defined(__ISWINDOWS__)
// Define compiler specific calling conventions
// for the shared library.
#  define CALLCONV
// Define function names for the shared library,
// in this case it is the librefprop.so and the
// names might change on some systems during
// the compilation of the Fortran files.
// Possible other branches for this code could be:
// #    if !defined(_AIX)
// #    if !defined(__hpux)
// #    ifdef _CRAY
// However, I cannot test that and therefore do not include it.
#  define RPVersion rpversion_
#  define SETPATHdll setpathdll_
#  define ABFL1dll abfl1dll_
#  define ABFL2dll abfl2dll_
#  define ACTVYdll actvydll_
#  define AGdll agdll_
#  define CCRITdll ccritdll_
#  define CP0dll cp0dll_
#  define CRITPdll critpdll_
#  define CSATKdll csatkdll_
#  define CV2PKdll cv2pkdll_
#  define CVCPKdll cvcpkdll_
#  define CVCPdll cvcpdll_
#  define DBDTdll dbdtdll_
#  define DBFL1dll dbfl1dll_
#  define DBFL2dll dbfl2dll_
#  define DDDPdll dddpdll_
#  define DDDTdll dddtdll_
#  define DEFLSHdll deflshdll_
#  define DHD1dll dhd1dll_
#  define DHFL1dll dhfl1dll_
#  define DHFL2dll dhfl2dll_
#  define DHFLSHdll dhflshdll_
#  define DIELECdll dielecdll_
#  define DOTFILLdll dotfilldll_
#  define DPDD2dll dpdd2dll_
#  define DPDDKdll dpddkdll_
#  define DPDDdll dpdddll_
#  define DPDTKdll dpdtkdll_
#  define DPDTdll dpdtdll_
#  define DPTSATKdll dptsatkdll_
#  define DSFLSHdll dsflshdll_
#  define DSFL1dll dsfl1dll_
#  define DSFL2dll dsfl2dll_
#  define ENTHALdll enthaldll_
#  define ENTROdll entrodll_
#  define ESFLSHdll esflshdll_
#  define FGCTYdll fgctydll_
#  define FUGCOFdll fugcofdll_
#  define FPVdll fpvdll_
#  define GERG04dll gerg04dll_
#  define GETFIJdll getfijdll_
#  define GETKTVdll getktvdll_
#  define GIBBSdll gibbsdll_
#  define HSFLSHdll hsflshdll_
#  define INFOdll infodll_
#  define LIMITKdll limitkdll_
#  define LIMITSdll limitsdll_
#  define LIMITXdll limitxdll_
#  define MELTPdll meltpdll_
#  define MELTTdll melttdll_
#  define MLTH2Odll mlth2odll_
#  define NAMEdll namedll_
#  define PASSCMNdll passcmn_
#  define PDFL1dll pdfl1dll_
#  define PDFLSHdll pdflshdll_
#  define PEFLSHdll peflshdll_
#  define PHFL1dll phfl1dll_
#  define PHFLSHdll phflshdll_
#  define PHIXdll phixdll_
#  define PHI0dll phi0dll_
#  define PQFLSHdll pqflshdll_
#  define PREOSdll preosdll_
#  define PRESSdll pressdll_
#  define PSFL1dll psfl1dll_
#  define PSFLSHdll psflshdll_
#  define PUREFLDdll pureflddll_
#  define QMASSdll qmassdll_
#  define QMOLEdll qmoledll_
#  define RESIDUALdll residualdll_
#  define RMIX2dll rmix2dll_
#  define REDXdll redxdll_
#  define SATDdll satddll_
#  define SATEdll satedll_
#  define SATHdll sathdll_
#  define SATPdll satpdll_
#  define SATSdll satsdll_
#  define SATTdll sattdll_
#  define SATSPLNdll satsplndll_
#  define SETAGAdll setagadll_
#  define SETKTVdll setktvdll_
#  define SETMIXdll setmixdll_
#  define SETMODdll setmoddll_
#  define SETREFdll setrefdll_
#  define SETUPdll setupdll_
//#  define SPECGRdll specgrdll_ // not found in library
#  define SUBLPdll sublpdll_
#  define SUBLTdll subltdll_
#  define SURFTdll surftdll_
#  define SURTENdll surtendll_
#  define TDFLSHdll tdflshdll_
#  define TEFLSHdll teflshdll_
#  define THERM0dll therm0dll_
#  define THERM2dll therm2dll_
#  define THERM3dll therm3dll_
#  define THERMdll thermdll_
#  define THFLSHdll thflshdll_
#  define TPFLSHdll tpflshdll_
#  define TPFL2dll tpfl2dll_
#  define TPRHOdll tprhodll_
#  define TQFLSHdll tqflshdll_
#  define TRNPRPdll trnprpdll_
#  define TSFLSHdll tsflshdll_
#  define VIRBdll virbdll_
#  define VIRCdll vircdll_
#  define WMOLdll wmoldll_
#  define XMASSdll xmassdll_
#  define XMOLEdll xmoledll_
#elif defined(__ISAPPLE__) // defined(__ISLINUX__)
// Define compiler specific calling conventions
// for the shared library.
#  define CALLCONV
#  define RPVersion rpversion_
#  define SETPATHdll setpathdll_
#  define ABFL1dll abfl1dll_
#  define ABFL2dll abfl2dll_
#  define ACTVYdll actvydll_
#  define AGdll agdll_
#  define CCRITdll ccritdll_
#  define CP0dll cp0dll_
#  define CRITPdll critpdll_
#  define CSATKdll csatkdll_
#  define CV2PKdll cv2pkdll_
#  define CVCPKdll cvcpkdll_
#  define CVCPdll cvcpdll_
#  define DBDTdll dbdtdll_
#  define DBFL1dll dbfl1dll_
#  define DBFL2dll dbfl2dll_
#  define DDDPdll dddpdll_
#  define DDDTdll dddtdll_
#  define DEFLSHdll deflshdll_
#  define DHD1dll dhd1dll_
#  define DHFL1dll dhfl1dll_
#  define DHFL2dll dhfl2dll_
#  define DHFLSHdll dhflshdll_
#  define DIELECdll dielecdll_
#  define DOTFILLdll dotfilldll_
#  define DPDD2dll dpdd2dll_
#  define DPDDKdll dpddkdll_
#  define DPDDdll dpdddll_
#  define DPDTKdll dpdtkdll_
#  define DPDTdll dpdtdll_
#  define DPTSATKdll dptsatkdll_
#  define DSFLSHdll dsflshdll_
#  define DSFL1dll dsfl1dll_
#  define DSFL2dll dsfl2dll_
#  define ENTHALdll enthaldll_
#  define ENTROdll entrodll_
#  define ESFLSHdll esflshdll_
#  define FGCTYdll fgctydll_
#  define FUGCOFdll fugcofdll_
#  define FPVdll fpvdll_
#  define GERG04dll gerg04dll_
#  define GETFIJdll getfijdll_
#  define GETKTVdll getktvdll_
#  define GIBBSdll gibbsdll_
#  define HSFLSHdll hsflshdll_
#  define INFOdll infodll_
#  define LIMITKdll limitkdll_
#  define LIMITSdll limitsdll_
#  define LIMITXdll limitxdll_
#  define MELTPdll meltpdll_
#  define MELTTdll melttdll_
#  define MLTH2Odll mlth2odll_
#  define NAMEdll namedll_
#  define PASSCMNdll passcmn_
#  define PDFL1dll pdfl1dll_
#  define PDFLSHdll pdflshdll_
#  define PEFLSHdll peflshdll_
#  define PHFL1dll phfl1dll_
#  define PHFLSHdll phflshdll_
#  define PHIXdll phixdll_
#  define PHI0dll phi0dll_
#  define PQFLSHdll pqflshdll_
#  define PREOSdll preosdll_
#  define PRESSdll pressdll_
#  define PSFL1dll psfl1dll_
#  define PSFLSHdll psflshdll_
#  define PUREFLDdll pureflddll_
#  define QMASSdll qmassdll_
#  define QMOLEdll qmoledll_
#  define RESIDUALdll residualdll_
#  define REDXdll redxdll_
#  define RMIX2dll rmix2dll_
#  define SATDdll satddll_
#  define SATEdll satedll_
#  define SATHdll sathdll_
#  define SATPdll satpdll_
#  define SATSdll satsdll_
#  define SATTdll sattdll_
#  define SATSPLNdll satsplndll_
#  define SETAGAdll setagadll_
#  define SETKTVdll setktvdll_
#  define SETMIXdll setmixdll_
#  define SETMODdll setmoddll_
#  define SETREFdll setrefdll_
#  define SETUPdll setupdll_
//#  define SPECGRdll specgrdll_ // not found in library
#  define SUBLPdll sublpdll_
#  define SUBLTdll subltdll_
#  define SURFTdll surftdll_
#  define SURTENdll surtendll_
#  define TDFLSHdll tdflshdll_
#  define TEFLSHdll teflshdll_
#  define THERM0dll therm0dll_
#  define THERM2dll therm2dll_
#  define THERM3dll therm3dll_
#  define THERMdll thermdll_
#  define THFLSHdll thflshdll_
#  define TPFLSHdll tpflshdll_
#  define TPFL2dll tpfl2dll_
#  define TPRHOdll tprhodll_
#  define TQFLSHdll tqflshdll_
#  define TRNPRPdll trnprpdll_
#  define TSFLSHdll tsflshdll_
#  define VIRBdll virbdll_
#  define VIRCdll vircdll_
#  define WMOLdll wmoldll_
#  define XMASSdll xmassdll_
#  define XMOLEdll xmoledll_
#else // #elif defined(__ISAPPLE__)
// Set some dummy names for the compiler
#  define CALLCONV
#  define RPVersion NOTAVAILABLE
#  define SETPATHdll setpathdll
#  define ABFL1dll abfl1dll
#  define ABFL2dll abfl2dll
#  define ACTVYdll actvydll
#  define AGdll agdll
#  define CCRITdll ccritdll
#  define CP0dll cp0dll
#  define CRITPdll critpdll
#  define CSATKdll csatkdll
#  define CV2PKdll cv2pkdll
#  define CVCPKdll cvcpkdll
#  define CVCPdll cvcpdll
#  define DBDTdll dbdtdll
#  define DBFL1dll dbfl1dll
#  define DBFL2dll dbfl2dll
#  define DDDPdll dddpdll
#  define DDDTdll dddtdll
#  define DEFLSHdll deflshdll
#  define DHD1dll dhd1dll
#  define DHFL1dll dhfl1dll
#  define DHFL2dll dhfl2dll
#  define DHFLSHdll dhflshdll
#  define DIELECdll dielecdll
#  define DOTFILLdll dotfilldll
#  define DPDD2dll dpdd2dll
#  define DPDDKdll dpddkdll
#  define DPDDdll dpdddll
#  define DPDTKdll dpdtkdll
#  define DPDTdll dpdtdll
#  define DPTSATKdll dptsatkdll
#  define DSFLSHdll dsflshdll
#  define DSFL1dll dsfl1dll
#  define DSFL2dll dsfl2dll
#  define ENTHALdll enthaldll
#  define ENTROdll entrodll
#  define ESFLSHdll esflshdll
#  define FGCTYdll fgctydll
#  define FPVdll fpvdll
#  define GERG04dll gerg04dll
#  define GETFIJdll getfijdll
#  define GETKTVdll getktvdll
#  define GIBBSdll gibbsdll
#  define HSFLSHdll hsflshdll
#  define INFOdll infodll
#  define LIMITKdll limitkdll
#  define LIMITSdll limitsdll
#  define LIMITXdll limitxdll
#  define MELTPdll meltpdll
#  define MELTTdll melttdll
#  define MLTH2Odll mlth2odll
#  define NAMEdll namedll
#  define PASSCMNdll passcmn
#  define PDFL1dll pdfl1dll
#  define PDFLSHdll pdflshdll
#  define PEFLSHdll peflshdll
#  define PHFL1dll phfl1dll
#  define PHFLSHdll phflshdll
#  define PHIXdll phixdll
#  define PHI0dll phi0dll
#  define PQFLSHdll pqflshdll
#  define PREOSdll preosdll
#  define PRESSdll pressdll
#  define PSFL1dll psfl1dll
#  define PSFLSHdll psflshdll
#  define PUREFLDdll pureflddll
#  define QMASSdll qmassdll
#  define QMOLEdll qmoledll
#  define RESIDUALdll residualdll
#  define REDXdll redxdll
#  define RMIX2dll rmix2dll
#  define SATDdll satddll
#  define SATEdll satedll
#  define SATHdll sathdll
#  define SATPdll satpdll
#  define SATSdll satsdll
#  define SATTdll sattdll
#  define SATSPLNdll satsplndll
#  define SETAGAdll setagadll
#  define SETKTVdll setktvdll
#  define SETMIXdll setmixdll
#  define SETMODdll setmoddll
#  define SETREFdll setrefdll
#  define SETUPdll setupdll
//#  define SPECGRdll specgrdll // not found in library
#  define SUBLPdll sublpdll
#  define SUBLTdll subltdll
#  define SURFTdll surftdll
#  define SURTENdll surtendll
#  define TDFLSHdll tdflshdll
#  define TEFLSHdll teflshdll
#  define THERM0dll therm0dll
#  define THERM2dll therm2dll
#  define THERM3dll therm3dll
#  define THERMdll thermdll
#  define THFLSHdll thflshdll
#  define TPFLSHdll tpflshdll
#  define TPFL2dll tpfl2dll
#  define TPRHOdll tprhodll
#  define TQFLSHdll tqflshdll
#  define TRNPRPdll trnprpdll
#  define TSFLSHdll tsflshdll
#  define VIRBdll virbdll
#  define VIRCdll vircdll
#  define WMOLdll wmoldll
#  define XMASSdll xmassdll
#  define XMOLEdll xmoledll
#endif // else branch
//
//
// Only continue if function names have been defined.
// We might want to include some more tests here...
#if defined(RPVersion)
// define new macros for function names
// http://stackoverflow.com/questions/195975/how-to-make-a-char-string-from-a-c-macros-value
#include <string.h>
#define STR_VALUE(arg)      #arg
#define FUNCTION_NAME(name) STR_VALUE(name)
#define STRINGIFY(name) STR_VALUE(name)
//
// I'll try to follow this example from:
// http://www.gershnik.com/tips/cpp.asp
// function type: typedef void [compiler stuff]  func_t(int, float);
// function declaration: func_t func;
// pointer type: typedef func_t * func_ptr;
#ifdef __cplusplus
extern "C" {
#endif
  // For C calling conventions, replaced all "double &" with "double *", and "long &" with "long *"
  typedef void (CALLCONV RPVersion_TYPE)( char* );
  typedef void (CALLCONV SETPATHdll_TYPE)( const char* );
  typedef void (CALLCONV ABFL1dll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV ABFL2dll_TYPE)(double *,double *,double *,long *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV ACTVYdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV AGdll_TYPE)(double *,double *,double *,double *,double *);
  typedef void (CALLCONV CCRITdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV CP0dll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV CRITPdll_TYPE)(double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV CSATKdll_TYPE)(long *,double *,long *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV CV2PKdll_TYPE)(long *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV CVCPKdll_TYPE)(long *,double *,double *,double *,double *);
  typedef void (CALLCONV CVCPdll_TYPE)(double *,double *,double *,double *,double *);
  typedef void (CALLCONV DBDTdll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV DBFL1dll_TYPE)(double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DBFL2dll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DDDPdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DDDTdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DEFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DHD1dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV DHFL1dll_TYPE)(double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DHFL2dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DHFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DIELECdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DOTFILLdll_TYPE)(long *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DPDD2dll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DPDDKdll_TYPE)(long *,double *,double *,double *);
  typedef void (CALLCONV DPDDdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DPDTKdll_TYPE)(long *,double *,double *,double *);
  typedef void (CALLCONV DPDTdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV DPTSATKdll_TYPE)(long *,double *,long *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DSFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DSFL1dll_TYPE)(double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV DSFL2dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV ENTHALdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV ENTROdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV ESFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV FGCTYdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV FPVdll_TYPE)(double *,double *,double *,double *,double *);
  typedef void (CALLCONV FUGCOFdll_TYPE)(double *,double *,double *,double*, long *,char*,long );
  typedef void (CALLCONV GERG04dll_TYPE)(long *,long *,long *,char*,long );
  typedef void (CALLCONV GETFIJdll_TYPE)(char*,double *,char*,char*,long ,long ,long );
  typedef void (CALLCONV GETKTVdll_TYPE)(long *,long *,char*,double *,char*,char*,char*,char*,long ,long ,long ,long ,long );
  typedef void (CALLCONV GIBBSdll_TYPE)(double *,double *,double *,double *,double *);
  typedef void (CALLCONV HSFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV INFOdll_TYPE)(long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV LIMITKdll_TYPE)(char*,long *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long ,long );
  typedef void (CALLCONV LIMITSdll_TYPE)(char*,double *,double *,double *,double *,double *,long );
  typedef void (CALLCONV LIMITXdll_TYPE)(char*,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long ,long );
  typedef void (CALLCONV MELTPdll_TYPE)(double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV MELTTdll_TYPE)(double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV MLTH2Odll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV NAMEdll_TYPE)(long *,char*,char*,char*,long ,long ,long );
  typedef void (CALLCONV PASSCMNdll_TYPE)(char *,long *,long *,long *,char *,long*,double *, double *, long*, char*, long, long, long);
  typedef void (CALLCONV PDFL1dll_TYPE)(double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PDFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PEFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PHFL1dll_TYPE)(double *,double *,double *,long *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PHFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PHIXdll_TYPE)(long *,long *,double *,double *,double *, double *);
  typedef void (CALLCONV PHI0dll_TYPE)(long *,long *,double *,double *,double *, double *);
  typedef void (CALLCONV PQFLSHdll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PREOSdll_TYPE)(long *);
  typedef void (CALLCONV PRESSdll_TYPE)(double *,double *,double *,double *);
  typedef void (CALLCONV PSFL1dll_TYPE)(double *,double *,double *,long *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PSFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV PUREFLDdll_TYPE)(long *);
  typedef void (CALLCONV QMASSdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV QMOLEdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV RESIDUALdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV REDXdll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV RMIX2dll_TYPE)(double *,double *);
  typedef void (CALLCONV SATDdll_TYPE)(double *,double *,long *,long *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATEdll_TYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATHdll_TYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATPdll_TYPE)(double *,double *,long *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATSdll_TYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATTdll_TYPE)(double *,double *,long *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SATSPLNdll_TYPE)(double *,long *,char*,long );
  typedef void (CALLCONV SETAGAdll_TYPE)(long *,char*,long );
  typedef void (CALLCONV SETKTVdll_TYPE)(long *,long *,char*,double *,char*,long *,char*,long ,long ,long );
  typedef void (CALLCONV SETMIXdll_TYPE)(char*,char*,char*,long *,char*,double *,long *,char*,long ,long ,long ,long ,long );
  typedef void (CALLCONV SETMODdll_TYPE)(long *,char*,char*,char*,long *,char*,long ,long ,long ,long );
  typedef void (CALLCONV SETREFdll_TYPE)(char*,long *,double *,double *,double *,double *,double *,long *,char*,long ,long );
  typedef void (CALLCONV SETUPdll_TYPE)(long *,char*,char*,char*,long *,char*,long ,long ,long ,long );
  //typedef void (CALLCONV SETUPdll_TYPE)(long *,char*,char*,char*,long *,char*);
//  typedef void (CALLCONV SPECGRdll_TYPE)(double *,double *,double *,double *); // not found in library
  typedef void (CALLCONV SUBLPdll_TYPE)(double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SUBLTdll_TYPE)(double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SURFTdll_TYPE)(double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV SURTENdll_TYPE)(double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TDFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TEFLSHdll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV THERM0dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV THERM2dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV THERM3dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV THERMdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
  typedef void (CALLCONV THFLSHdll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TPFLSHdll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TPFL2dll_TYPE)(double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TPRHOdll_TYPE)(double *,double *,double *,long *,long *,double *,long *,char*,long );
  typedef void (CALLCONV TQFLSHdll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TRNPRPdll_TYPE)(double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV TSFLSHdll_TYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
  typedef void (CALLCONV VIRBdll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV VIRCdll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV WMOLdll_TYPE)(double *,double *);
  typedef void (CALLCONV XMASSdll_TYPE)(double *,double *,double *);
  typedef void (CALLCONV XMOLEdll_TYPE)(double *,double *,double *);
  //
  // Define explicit function pointers
  typedef RPVersion_TYPE * RPVersion_POINTER;
  typedef SETPATHdll_TYPE * SETPATHdll_POINTER;
  typedef ABFL1dll_TYPE * ABFL1dll_POINTER;
  typedef ABFL2dll_TYPE * ABFL2dll_POINTER;
  typedef ACTVYdll_TYPE * ACTVYdll_POINTER;
  typedef AGdll_TYPE * AGdll_POINTER;
  typedef CCRITdll_TYPE * CCRITdll_POINTER;
  typedef CP0dll_TYPE * CP0dll_POINTER;
  typedef CRITPdll_TYPE * CRITPdll_POINTER;
  typedef CSATKdll_TYPE * CSATKdll_POINTER;
  typedef CV2PKdll_TYPE * CV2PKdll_POINTER;
  typedef CVCPKdll_TYPE * CVCPKdll_POINTER;
  typedef CVCPdll_TYPE * CVCPdll_POINTER;
  typedef DBDTdll_TYPE * DBDTdll_POINTER;
  typedef DBFL1dll_TYPE * DBFL1dll_POINTER;
  typedef DBFL2dll_TYPE * DBFL2dll_POINTER;
  typedef DDDPdll_TYPE * DDDPdll_POINTER;
  typedef DDDTdll_TYPE * DDDTdll_POINTER;
  typedef DEFLSHdll_TYPE * DEFLSHdll_POINTER;
  typedef DHD1dll_TYPE * DHD1dll_POINTER;
  typedef DHFLSHdll_TYPE * DHFLSHdll_POINTER;
  typedef DHFL1dll_TYPE * DHFL1dll_POINTER;
  typedef DHFL2dll_TYPE * DHFL2dll_POINTER;
  typedef DIELECdll_TYPE * DIELECdll_POINTER;
  typedef DOTFILLdll_TYPE * DOTFILLdll_POINTER;
  typedef DPDD2dll_TYPE * DPDD2dll_POINTER;
  typedef DPDDKdll_TYPE * DPDDKdll_POINTER;
  typedef DPDDdll_TYPE * DPDDdll_POINTER;
  typedef DPDTKdll_TYPE * DPDTKdll_POINTER;
  typedef DPDTdll_TYPE * DPDTdll_POINTER;
  typedef DPTSATKdll_TYPE * DPTSATKdll_POINTER;
  typedef DSFLSHdll_TYPE * DSFLSHdll_POINTER;
  typedef DSFL1dll_TYPE * DSFL1dll_POINTER;
  typedef DSFL2dll_TYPE * DSFL2dll_POINTER;
  typedef ENTHALdll_TYPE * ENTHALdll_POINTER;
  typedef ENTROdll_TYPE * ENTROdll_POINTER;
  typedef ESFLSHdll_TYPE * ESFLSHdll_POINTER;
  typedef FGCTYdll_TYPE * FGCTYdll_POINTER;
  typedef FPVdll_TYPE * FPVdll_POINTER;
  typedef FUGCOFdll_TYPE * FUGCOFdll_POINTER;
  typedef GERG04dll_TYPE * GERG04dll_POINTER;
  typedef GETFIJdll_TYPE * GETFIJdll_POINTER;
  typedef GETKTVdll_TYPE * GETKTVdll_POINTER;
  typedef GIBBSdll_TYPE * GIBBSdll_POINTER;
  typedef HSFLSHdll_TYPE * HSFLSHdll_POINTER;
  typedef INFOdll_TYPE * INFOdll_POINTER;
  typedef LIMITKdll_TYPE * LIMITKdll_POINTER;
  typedef LIMITSdll_TYPE * LIMITSdll_POINTER;
  typedef LIMITXdll_TYPE * LIMITXdll_POINTER;
  typedef MELTPdll_TYPE * MELTPdll_POINTER;
  typedef MELTTdll_TYPE * MELTTdll_POINTER;
  typedef MLTH2Odll_TYPE * MLTH2Odll_POINTER;
  typedef NAMEdll_TYPE * NAMEdll_POINTER;
  typedef PASSCMNdll_TYPE * PASSCMNdll_POINTER;
  typedef PDFL1dll_TYPE * PDFL1dll_POINTER;
  typedef PDFLSHdll_TYPE * PDFLSHdll_POINTER;
  typedef PEFLSHdll_TYPE * PEFLSHdll_POINTER;
  typedef PHFL1dll_TYPE * PHFL1dll_POINTER;
  typedef PHFLSHdll_TYPE * PHFLSHdll_POINTER;
  typedef PHIXdll_TYPE * PHIXdll_POINTER;
  typedef PHI0dll_TYPE * PHI0dll_POINTER;
  typedef PQFLSHdll_TYPE * PQFLSHdll_POINTER;
  typedef PREOSdll_TYPE * PREOSdll_POINTER;
  typedef PRESSdll_TYPE * PRESSdll_POINTER;
  typedef PSFL1dll_TYPE * PSFL1dll_POINTER;
  typedef PSFLSHdll_TYPE * PSFLSHdll_POINTER;
  typedef PUREFLDdll_TYPE * PUREFLDdll_POINTER;
  typedef QMASSdll_TYPE * QMASSdll_POINTER;
  typedef QMOLEdll_TYPE * QMOLEdll_POINTER;
  typedef RESIDUALdll_TYPE * RESIDUALdll_POINTER;
  typedef REDXdll_TYPE * REDXdll_POINTER;
  typedef RMIX2dll_TYPE * RMIX2dll_POINTER;
  typedef SATDdll_TYPE * SATDdll_POINTER;
  typedef SATEdll_TYPE * SATEdll_POINTER;
  typedef SATHdll_TYPE * SATHdll_POINTER;
  typedef SATPdll_TYPE * SATPdll_POINTER;
  typedef SATSdll_TYPE * SATSdll_POINTER;
  typedef SATTdll_TYPE * SATTdll_POINTER;
  typedef SATSPLNdll_TYPE * SATSPLNdll_POINTER;
  typedef SETAGAdll_TYPE * SETAGAdll_POINTER;
  typedef SETKTVdll_TYPE * SETKTVdll_POINTER;
  typedef SETMIXdll_TYPE * SETMIXdll_POINTER;
  typedef SETMODdll_TYPE * SETMODdll_POINTER;
  typedef SETREFdll_TYPE * SETREFdll_POINTER;
  typedef SETUPdll_TYPE * SETUPdll_POINTER;
//  typedef SPECGRdll_TYPE * SPECGRdll_POINTER; // not found in library
  typedef SUBLPdll_TYPE * SUBLPdll_POINTER;
  typedef SUBLTdll_TYPE * SUBLTdll_POINTER;
  typedef SURFTdll_TYPE * SURFTdll_POINTER;
  typedef SURTENdll_TYPE * SURTENdll_POINTER;
  typedef TDFLSHdll_TYPE * TDFLSHdll_POINTER;
  typedef TEFLSHdll_TYPE * TEFLSHdll_POINTER;
  typedef THERM0dll_TYPE * THERM0dll_POINTER;
  typedef THERM2dll_TYPE * THERM2dll_POINTER;
  typedef THERM3dll_TYPE * THERM3dll_POINTER;
  typedef THERMdll_TYPE * THERMdll_POINTER;
  typedef THFLSHdll_TYPE * THFLSHdll_POINTER;
  typedef TPFLSHdll_TYPE * TPFLSHdll_POINTER;
  typedef TPFL2dll_TYPE * TPFL2dll_POINTER;
  typedef TPRHOdll_TYPE * TPRHOdll_POINTER;
  typedef TQFLSHdll_TYPE * TQFLSHdll_POINTER;
  typedef TRNPRPdll_TYPE * TRNPRPdll_POINTER;
  typedef TSFLSHdll_TYPE * TSFLSHdll_POINTER;
  typedef VIRBdll_TYPE * VIRBdll_POINTER;
  typedef VIRCdll_TYPE * VIRCdll_POINTER;
  typedef WMOLdll_TYPE * WMOLdll_POINTER;
  typedef XMASSdll_TYPE * XMASSdll_POINTER;
  typedef XMOLEdll_TYPE * XMOLEdll_POINTER;
#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus
#endif // defined(RPversion)
#endif // REFPROP_LIB_H
