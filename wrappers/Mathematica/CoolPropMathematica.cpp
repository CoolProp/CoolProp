
/* Include required header */
#include "CoolProp.h"
#include "HumidAirProp.h"
#include "WolframLibrary.h"

/* Return the version of Library Link */
extern "C" DLLEXPORT mint WolframLibrary_getVersion( ) {
	return WolframLibraryVersion;
}

/* Initialize Library */
extern "C" DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {
	return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
extern "C" DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {
	return;
}

/* Adds one to the input, returning the result  */
extern "C" DLLEXPORT int plus_one( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
	mreal I1;
	I1 = MArgument_getReal(Args[0]);
	MArgument_setReal(Res, I1+1.0);
	return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int PropsSI( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    if (Argc != 6){ return LIBRARY_FUNCTION_ERROR; }
    char *outstring = MArgument_getUTF8String(Args[0]);
    char *In1string = MArgument_getUTF8String(Args[1]);
    mreal In1val = MArgument_getReal(Args[2]);
    char *In2string = MArgument_getUTF8String(Args[3]);
    mreal In2val = MArgument_getReal(Args[4]);
    char *Fluidstring = MArgument_getUTF8String(Args[5]);
    
    // PropsS version takes all strings, not single-character inputs
	double val = CoolProp::PropsSI(outstring,In1string,(double)In1val,In2string,(double)In2val,Fluidstring);
	MArgument_setReal(Res, val);
	return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int HAPropsSI( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    if (Argc != 7){ return LIBRARY_FUNCTION_ERROR; }
    char *outstring = MArgument_getUTF8String(Args[0]);
    char *In1string = MArgument_getUTF8String(Args[1]);
    mreal In1val = MArgument_getReal(Args[2]);
    char *In2string = MArgument_getUTF8String(Args[3]);
    mreal In2val = MArgument_getReal(Args[4]);
    char *In3string = MArgument_getUTF8String(Args[5]);
    mreal In3val = MArgument_getReal(Args[6]);
    
	double val = HumidAir::HAPropsSI(outstring,In1string,(double)In1val,In2string,(double)In2val,In3string,(double)In3val);
	MArgument_setReal(Res, val);
	return LIBRARY_NO_ERROR;
}

	