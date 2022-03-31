
/* Include required header */
#include "CoolProp.h"
#include "HumidAirProp.h"
#include "WolframLibrary.h"

/* Return the version of Library Link */
extern "C" DLLEXPORT mint WolframLibrary_getVersion() {
    return WolframLibraryVersion;
}

/* Initialize Library */
extern "C" DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
    return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
extern "C" DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
    return;
}

/* Adds one to the input, returning the result */
extern "C" DLLEXPORT int plus_one(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;
    double x = MArgument_getReal(Args[0]);
    MArgument_setReal(Res, x + 1.0);
    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int PropsSI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 6) return LIBRARY_FUNCTION_ERROR;

    char* Output = MArgument_getUTF8String(Args[0]);
    char* Name1 = MArgument_getUTF8String(Args[1]);
    double Prop1 = MArgument_getReal(Args[2]);
    char* Name2 = MArgument_getUTF8String(Args[3]);
    double Prop2 = MArgument_getReal(Args[4]);
    char* FluidName = MArgument_getUTF8String(Args[5]);

    double val = CoolProp::PropsSI(Output, Name1, Prop1, Name2, Prop2, FluidName);

    libData->UTF8String_disown(Output);
    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(FluidName);

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int HAPropsSI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 7) return LIBRARY_FUNCTION_ERROR;

    char* Output = MArgument_getUTF8String(Args[0]);
    char* Name1 = MArgument_getUTF8String(Args[1]);
    double Prop1 = MArgument_getReal(Args[2]);
    char* Name2 = MArgument_getUTF8String(Args[3]);
    double Prop2 = MArgument_getReal(Args[4]);
    char* Name3 = MArgument_getUTF8String(Args[5]);
    double Prop3 = MArgument_getReal(Args[6]);

    double val = HumidAir::HAPropsSI(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3);

    libData->UTF8String_disown(Output);
    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(Name3);

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}
