// Regression tests for CoolProp::fpu_guard (GitHub #3012 / bd CoolProp-i3o).
//
// CoolProp uses NaN and _HUGE as in-band sentinels.  Host environments that run
// with FP exceptions *unmasked* (Delphi, some TRNSYS/Fortran builds) take a
// hardware trap / SIGFPE the instant CoolProp touches one of those values,
// crashing mid-calculation.  fpu_guard masks FP exceptions on construction and
// restores the caller's environment on destruction.
//
// These tests assert the guard's contract directly (it is a header-only type,
// so it is reachable from the Catch2 runner even though CoolPropLib.cpp -- the
// C-interface that wires the guard into every entry point -- is excluded from
// that binary by CMakeLists.txt).  The decisive check is the divide-and-invalid
// inside the guarded scope: with FE_INVALID *enabled* and no masking, the
// volatile 0.0/0.0 below would raise SIGFPE and abort the process.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>
#    include "CoolProp/FPUGuard.h"
#    include <cmath>  // std::isnan
#    include <cfenv>  // std::fetestexcept / std::feclearexcept

// fegetexcept/feenableexcept/fedisableexcept are GNU glibc extensions; the
// "trapping host" half of the contract can only be exercised where they exist.
#    if defined(__GLIBC__) && defined(FE_ALL_EXCEPT)

TEST_CASE("fpu_guard masks FP traps inside its scope and survives a NaN op", "[fpu_guard][fpu][3012]") {
    std::feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID);  // emulate a trapping host (Delphi/TRNSYS)

    {
        CoolProp::fpu_guard guard;
        // Inside the guard FE_INVALID must be masked, so this 0.0/0.0 (which a
        // trapping host would fault on) produces a quiet NaN instead of SIGFPE.
        CHECK((fegetexcept() & FE_INVALID) == 0);
        volatile double zero = 0.0;
        volatile double nan = zero / zero;
        CHECK(std::isnan(nan));
    }

    // After the guard the caller's trapping configuration is restored ...
    CHECK((fegetexcept() & FE_INVALID) != 0);
    // ... and the status flags the guard raised are cleared, so a polling host
    // (Excel/VBA) sees a clean status word.
    CHECK(std::fetestexcept(FE_INVALID) == 0);

    fedisableexcept(FE_ALL_EXCEPT);
    std::feclearexcept(FE_ALL_EXCEPT);
}

TEST_CASE("fpu_guard leaves an already-masked environment unchanged", "[fpu_guard][fpu][3012]") {
    fedisableexcept(FE_ALL_EXCEPT);  // nothing trapping (the common host default)
    std::feclearexcept(FE_ALL_EXCEPT);

    {
        CoolProp::fpu_guard guard;
        CHECK((fegetexcept() & FE_ALL_EXCEPT) == 0);
    }

    // Restoring an empty trapping set must not enable anything spuriously.
    CHECK((fegetexcept() & FE_ALL_EXCEPT) == 0);
    std::feclearexcept(FE_ALL_EXCEPT);
}
#    endif  // __GLIBC__ && FE_ALL_EXCEPT

TEST_CASE("fpu_guard is constructible and clears status flags on every platform", "[fpu_guard][fpu][3012]") {
    // Platform-independent smoke test: the guard must compile, run a NaN op
    // without aborting even where the GNU trapping API is unavailable, and -- on
    // every platform path (glibc/MSVC/macOS) -- leave a clean status word once
    // destroyed, since each ~fpu_guard clears the flags it raised.
    std::feclearexcept(FE_ALL_EXCEPT);
    {
        CoolProp::fpu_guard guard;
        volatile double zero = 0.0;
        volatile double nan = zero / zero;  // raises FE_INVALID inside the guard
        CHECK(std::isnan(nan));
    }
    // The guard is now destroyed; the flag it raised must have been cleared.
    CHECK(std::fetestexcept(FE_INVALID) == 0);
}

#endif  // ENABLE_CATCH
