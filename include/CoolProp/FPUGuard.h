#ifndef COOLPROP_FPUGUARD_H
#define COOLPROP_FPUGUARD_H

// RAII guard that masks IEEE-754 floating-point exceptions on construction and
// restores the caller's FP environment on destruction.
//
// CoolProp uses NaN and _HUGE (HUGE_VAL) as in-band sentinels.  Arithmetic and
// comparisons on those values set the FP exception *status* flags (FE_INVALID,
// FE_OVERFLOW, FE_DIVBYZERO).  Host environments react two ways:
//   * Polling hosts (Excel/VBA) read the status word after a call returns and
//     report a spurious FP error.
//   * Trapping hosts (Delphi, some TRNSYS/Fortran builds) run with FP
//     exceptions *unmasked*, so the CPU raises a hardware trap / SIGFPE at the
//     offending instruction -- mid-calculation, before any cleanup can run.
//
// Clearing the status flags on exit (the historical fpu_reset_guard behavior)
// only helps polling hosts.  This guard additionally *masks* the exceptions on
// entry so trapping hosts do not fault, then restores the caller's environment
// on exit (which also wipes any flags we raised, covering the polling case).
//
// Windows is the platform that hurts most here: Delphi's RTL *unmasks* FP
// exceptions by default, and Excel/VBA poll the status word.  So the Windows
// masking path is keyed on _WIN32 (any compiler), not just _MSC_VER, so that
// MinGW-built DLLs are protected too -- _controlfp_s lives in <float.h> for
// both MSVC and MinGW-w64.
//
// See GitHub issue #3012 and bd CoolProp-i3o.

#if defined(_WIN32)
#    include <float.h>  // _controlfp_s, _clearfp, _MCW_EM (MS CRT extension)
#    define COOLPROP_FPUGUARD_WIN
#elif defined(__GLIBC__)
#    include <cfenv>  // fegetexcept/fedisableexcept/feenableexcept are glibc-only
#    if defined(FE_ALL_EXCEPT)
#        define COOLPROP_FPUGUARD_GLIBC
#    endif
#else
#    include <cfenv>  // standard feclearexcept only (no portable masking API)
#endif

namespace CoolProp {

class fpu_guard
{
   public:
    fpu_guard() {
#if defined(COOLPROP_FPUGUARD_WIN)
        // Read the current control word, then set every exception-mask bit
        // (_MCW_EM): a set bit *disables* (masks) that exception's trap.
        unsigned int dummy = 0;
        _controlfp_s(&saved_cw, 0, 0);           // snapshot current control word
        _controlfp_s(&dummy, _MCW_EM, _MCW_EM);  // mask all FP exceptions
#elif defined(COOLPROP_FPUGUARD_GLIBC)
        // fegetexcept() returns the set of currently *enabled* (trapping)
        // exceptions; remember it so we can restore the caller's choice.
        saved_excepts = fegetexcept();
        fedisableexcept(FE_ALL_EXCEPT);  // mask all
#endif
        // Platforms without a portable masking API (macOS/BSD, PowerPC) fall
        // through with no masking; they do not enable FP traps by default, so
        // clearing status flags on exit (see ~fpu_guard) is the sufficient,
        // historical mitigation there.
    }
    fpu_guard(const fpu_guard&) = delete;
    fpu_guard(fpu_guard&&) = delete;
    fpu_guard& operator=(const fpu_guard&) = delete;
    fpu_guard& operator=(fpu_guard&&) = delete;
    ~fpu_guard() {
#if defined(COOLPROP_FPUGUARD_WIN)
        _clearfp();  // clear flags we raised
        unsigned int dummy = 0;
        _controlfp_s(&dummy, saved_cw, _MCW_EM);  // restore caller's masking
#elif defined(COOLPROP_FPUGUARD_GLIBC)
        feclearexcept(FE_ALL_EXCEPT);   // clear flags we raised
        feenableexcept(saved_excepts);  // restore caller's trapping set
#elif defined(FE_ALL_EXCEPT)
        feclearexcept(FE_ALL_EXCEPT);  // fallback: clear status flags only
#endif
    }

   private:
#if defined(COOLPROP_FPUGUARD_WIN)
    unsigned int saved_cw{0};
#elif defined(COOLPROP_FPUGUARD_GLIBC)
    int saved_excepts{0};
#endif
};

}  // namespace CoolProp

#endif  // COOLPROP_FPUGUARD_H
