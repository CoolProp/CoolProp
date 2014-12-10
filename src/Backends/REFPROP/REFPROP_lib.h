
#ifndef REFPROP_LIB_H
#define REFPROP_LIB_H

/* See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 * This will be used to generate function names, pointers, etc. below
 */
#define LIST_OF_REFPROP_FUNCTION_NAMES \
	X(RPVersion) \
	X(SETPATHdll) \
	X(ABFL1dll) \
	X(ABFL2dll) \
	X(ACTVYdll) \
	X(AGdll) \
	X(CCRITdll) \
	X(CP0dll) \
	X(CRITPdll) \
	X(CSATKdll) \
	X(CV2PKdll) \
	X(CVCPKdll) \
	X(CVCPdll) \
	X(DBDTdll) \
	X(DBFL1dll) \
	X(DBFL2dll) \
	X(DDDPdll) \
	X(DDDTdll) \
	X(DEFLSHdll) \
	X(DHD1dll) \
	X(DHFL1dll) \
	X(DHFL2dll) \
	X(DHFLSHdll) \
	X(DIELECdll) \
	X(DOTFILLdll) \
	X(DPDD2dll) \
	X(DPDDKdll) \
	X(DPDDdll) \
	X(DPDTKdll) \
	X(DPDTdll) \
	X(DPTSATKdll) \
	X(DSFLSHdll) \
	X(DSFL1dll) \
	X(DSFL2dll) \
	X(ENTHALdll) \
	X(ENTROdll) \
	X(ESFLSHdll) \
	X(FGCTYdll) \
	X(FUGCOFdll) \
	X(FPVdll) \
	X(GERG04dll) \
	X(GETFIJdll) \
	X(GETKTVdll) \
	X(GIBBSdll) \
	X(HSFLSHdll) \
	X(INFOdll) \
	X(LIMITKdll) \
	X(LIMITSdll) \
	X(LIMITXdll) \
	X(MELTPdll) \
	X(MELTTdll) \
	X(MLTH2Odll) \
	X(NAMEdll) \
	X(PASSCMNdll) \
	X(PDFL1dll) \
	X(PDFLSHdll) \
	X(PEFLSHdll) \
	X(PHFL1dll) \
	X(PHFLSHdll) \
	X(PHIXdll) \
	X(PHI0dll) \
	X(PQFLSHdll) \
	X(PREOSdll) \
	X(PRESSdll) \
	X(PSFL1dll) \
	X(PSFLSHdll) \
	X(PUREFLDdll) \
	X(QMASSdll) \
	X(QMOLEdll) \
	X(RESIDUALdll) \
	X(REDXdll) \
	X(RMIX2dll) \
	X(SATDdll) \
	X(SATEdll) \
	X(SATHdll) \
	X(SATPdll) \
	X(SATSdll) \
	X(SATTdll) \
	X(SATSPLNdll) \
	X(SETAGAdll) \
	X(SETKTVdll) \
	X(SETMIXdll) \
	X(SETMODdll) \
	X(SETREFdll) \
	X(SETUPdll) \
	X(SPECGRdll) \
	X(SUBLPdll) \
	X(SUBLTdll) \
	X(SURFTdll) \
	X(SURTENdll) \
	X(TDFLSHdll) \
	X(TEFLSHdll) \
	X(THERM0dll) \
	X(THERM2dll) \
	X(THERM3dll) \
	X(THERMdll) \
	X(THFLSHdll) \
	X(TPFLSHdll) \
	X(TPFL2dll) \
	X(TPRHOdll) \
	X(TQFLSHdll) \
	X(TRNPRPdll) \
	X(TSFLSHdll) \
	X(VIRBdll) \
	X(VIRCdll) \
	X(WMOLdll) \
	X(XMASSdll) \
	X(XMOLEdll)

// Get the platform identifiers
#include "CoolPropTools.h"

// Define compiler specific calling conventions
// for the shared library.
#if defined(__ISWINDOWS__)
	#define CALLCONV __stdcall
#else
	#define CALLCONV
#endif

// define new macros for function names
// http://stackoverflow.com/questions/195975/how-to-make-a-char-string-from-a-c-macros-value
#include <string.h>
#define STR_VALUE(arg)      #arg
#define FUNCTION_NAME(name) STR_VALUE(name)
#define STRINGIFY(name) STR_VALUE(name)

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
  typedef void (CALLCONV SPECGRdll_TYPE)(double *,double *,double *,double *);
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
  
  /* Define explicit function pointers
   * Each will look something like: typedef RPVersion_TYPE * RPVersion_POINTER;
   * 
   * The ## are needed to escape the _ character in the variable names
   * 
   * ***MAGIC WARNING**!! X Macros in use
   * See http://stackoverflow.com/a/148610
   * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
   */
  #define X(name) typedef name ## _TYPE * name ## _POINTER;
     LIST_OF_REFPROP_FUNCTION_NAMES
  #undef X
  
#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus

#endif // REFPROP_LIB_H