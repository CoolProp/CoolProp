#ifndef COOLPROP_DETAIL_CONFIGURATION_KEYS_H
#define COOLPROP_DETAIL_CONFIGURATION_KEYS_H

// Lean home of the configuration-key X-macro and the `configuration_keys`
// enum.  Split out of Configuration.h so the enum can be consumed -- e.g. by
// the generated Cython constants module -- WITHOUT pulling in Configuration.h's
// JSON-(de)serialization includes.  This header has no includes by design;
// keep it that way.

#define CONFIGURATION_KEYS_ENUM                                                                                                                      \
    X(NORMALIZE_GAS_CONSTANTS, "NORMALIZE_GAS_CONSTANTS", true, "If true, for mixtures, the molar gas constant (R) will be set to the CODATA value") \
    X(CRITICAL_WITHIN_1UK, "CRITICAL_WITHIN_1UK", true,                                                                                              \
      "If true, any temperature within 1 uK of the critical temperature will be considered to be AT the critical point")                             \
    X(CRITICAL_SPLINES_ENABLED, "CRITICAL_SPLINES_ENABLED", true,                                                                                    \
      "If true, the critical splines will be used in the near-vicinity of the critical point")                                                       \
    X(SAVE_RAW_TABLES, "SAVE_RAW_TABLES", false, "If true, the raw, uncompressed tables will also be written to file")                               \
    X(ALTERNATIVE_TABLES_DIRECTORY, "ALTERNATIVE_TABLES_DIRECTORY", "",                                                                              \
      "If provided, this path will be the root directory for the tabular data.  Otherwise, ${HOME}/.CoolProp/Tables is used")                        \
    X(ALTERNATIVE_SVDTABLES_DIRECTORY, "ALTERNATIVE_SVDTABLES_DIRECTORY", "",                                                                        \
      "If provided, this path will be the root directory for the SVDSBTL on-disk cache (both .svd.bin.z surfaces and .critpatch.bin sidecars).  "    \
      "Otherwise, ${HOME}/.CoolProp/SVDTables is used.  Useful when ${HOME} is read-only (CI containers), on shared workstations, or for "           \
      "centrally-managed cache directories.")                                                                                                        \
    X(ALTERNATIVE_REFPROP_PATH, "ALTERNATIVE_REFPROP_PATH", "",                                                                                      \
      "An alternative path to be provided to the directory that contains REFPROP's fluids and mixtures directories.  If provided, the SETPATH "      \
      "function will be called with this directory prior to calling any REFPROP functions.")                                                         \
    X(ALTERNATIVE_REFPROP_HMX_BNC_PATH, "ALTERNATIVE_REFPROP_HMX_BNC_PATH", "",                                                                      \
      "An alternative path to the HMX.BNC file.  If provided, it will be passed into REFPROP's SETUP or SETMIX routines")                            \
    X(ALTERNATIVE_REFPROP_LIBRARY_PATH, "ALTERNATIVE_REFPROP_LIBRARY_PATH", "",                                                                      \
      "An alternative path to the shared library file.  If provided, it will be used to load REFPROP")                                               \
    X(REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS, "REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS", false,                                           \
      "If true, if the binary interaction parameters in REFPROP are estimated, throw an error rather than silently continuing")                      \
    X(REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS, "REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS", false,                         \
      "If true, if the binary interaction parameters in REFPROP are unable to be estimated, silently continue rather than failing")                  \
    X(REFPROP_USE_GERG, "REFPROP_USE_GERG", false,                                                                                                   \
      "If true, rather than using the highly-accurate pure fluid equations of state, use the pure-fluid EOS from GERG-2008")                         \
    X(REFPROP_ERROR_THRESHOLD, "REFPROP_ERROR_THRESHOLD", static_cast<int>(0), "The highest acceptable error code without throwing an exception")    \
    X(REFPROP_USE_PENGROBINSON, "REFPROP_USE_PENGROBINSON", false,                                                                                   \
      "If true, rather than using the highly-accurate pure fluid equations of state, use the Peng-Robinson EOS")                                     \
    X(REFPROP_RESOLVE_COOLPROP_ALIASES, "REFPROP_RESOLVE_COOLPROP_ALIASES", true,                                                                    \
      "If true (default), fluid names passed to the REFPROP backend are resolved through the CoolProp alias table before building the .FLD path.  "  \
      "Set to false to disable alias resolution and pass names directly to REFPROP, which may be needed for corner cases.")                          \
    X(MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB, "MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB", 1.0,                                                                 \
      "The maximum allowed size of the directory that is used to store tabular data")                                                                \
    X(DONT_CHECK_PROPERTY_LIMITS, "DONT_CHECK_PROPERTY_LIMITS", false,                                                                               \
      "If true, when possible, CoolProp will skip checking whether values are inside the property limits")                                           \
    X(HENRYS_LAW_TO_GENERATE_VLE_GUESSES, "HENRYS_LAW_TO_GENERATE_VLE_GUESSES", false,                                                               \
      "If true, when doing water-based mixture dewpoint calculations, use Henry's Law to generate guesses for liquid-phase composition")             \
    X(PHASE_ENVELOPE_STARTING_PRESSURE_PA, "PHASE_ENVELOPE_STARTING_PRESSURE_PA", 100.0, "Starting pressure [Pa] for phase envelope construction")   \
    X(R_U_CODATA, "R_U_CODATA", 8.31446261815324,                                                                                                    \
      "The value for the ideal gas constant in J/mol/K according to CODATA 2022.  This value is used to harmonize all the ideal gas constants. "     \
      "This is especially important in the critical region.")                                                                                        \
    X(VTPR_UNIFAC_PATH, "VTPR_UNIFAC_PATH", "", "The path to the directory containing the UNIFAC JSON files.  Should be slash terminated")           \
    X(SPINODAL_MINIMUM_DELTA, "SPINODAL_MINIMUM_DELTA", 0.5,                                                                                         \
      "The minimal delta to be used in tracing out the spinodal; make sure that the EOS has a spinodal at this value of delta=rho/rho_r")            \
    X(OVERWRITE_FLUIDS, "OVERWRITE_FLUIDS", false,                                                                                                   \
      "If true, and a fluid is added to the fluids library that is already there, rather than not adding the fluid (and probably throwing an "       \
      "exception), overwrite it")                                                                                                                    \
    X(OVERWRITE_DEPARTURE_FUNCTION, "OVERWRITE_DEPARTURE_FUNCTION", false,                                                                           \
      "If true, and a departure function to be added is already there, rather than not adding the departure function (and probably throwing an "     \
      "exception), overwrite it")                                                                                                                    \
    X(OVERWRITE_BINARY_INTERACTION, "OVERWRITE_BINARY_INTERACTION", false,                                                                           \
      "If true, and a pair of binary interaction pairs to be added is already there, rather than not adding the binary interaction pair (and "       \
      "probably throwing an exception), overwrite it")                                                                                               \
    X(USE_GUESSES_IN_PROPSSI, "USE_GUESSES_IN_PROPSSI", false,                                                                                       \
      "If true, calls to the vectorized versions of PropsSI use the previous state as guess value while looping over the input vectors, only makes " \
      "sense when working with a single fluid and with points that are not too far from each other.")                                                \
    X(ASSUME_CRITICAL_POINT_STABLE, "ASSUME_CRITICAL_POINT_STABLE", false,                                                                           \
      "If true, evaluation of the stability of critical point will be skipped and point will be assumed to be stable")                               \
    X(VTPR_ALWAYS_RELOAD_LIBRARY, "VTPR_ALWAYS_RELOAD_LIBRARY", false,                                                                               \
      "If true, the library will always be reloaded, no matter what is currently loaded")                                                            \
    X(FLOAT_PUNCTUATION, "FLOAT_PUNCTUATION", ".", "The first character of this string will be used as the separator between the number fraction.")  \
    X(ENABLE_SUPERANCILLARIES, "ENABLE_SUPERANCILLARIES", true, "If true, the superancillary functions will be used for VLE of pure fluids")         \
    X(ENABLE_MELTING_CALORIC_HS, "ENABLE_MELTING_CALORIC_HS", true,                                                                                  \
      "If true, HS_flash may seed the cold-compressed-liquid corner from the melting-line caloric Chebyshev (cascade leg 4). "                       \
      "Set false to force the legacy fallback for that region. Default: true")                                                                       \
    X(HSU_D_TWOPHASE_EOS_POLISH, "HSU_D_TWOPHASE_EOS_POLISH", true,                                                                                  \
      "If true (default), the superancillary D+{H,S,U} two-phase flash refines its fast superancillary-based saturation temperature with a short "   \
      "full-EOS secant polish, giving an EOS-exact result. If false, the (slightly faster) superancillary solution is returned directly with a "     \
      "~1e-8 deviation; still far inside typical tolerances. No effect on single-phase states.")                                                     \
    X(LIST_STRING_DELIMITER, "LIST_STRING_DELIMITER", ",", "The delimiter to be used when converting a list of strings to a string")                 \
    X(ALLOW_SVDSBTL_IN_PROPSSI, "ALLOW_SVDSBTL_IN_PROPSSI", false,                                                                                   \
      "If true, the SVDSBTL backend is usable through the high-level PropsSI interface.  By default it is rejected (mirroring the BICUBIC/TTSE "     \
      "tabular backends) because PropsSI rebuilds the AbstractState on every call, and loading an SVDSurface from cache costs ~80 ms -- over four "  \
      "orders of magnitude slower per call than the actual SVD eval.  For batched workloads use AbstractState directly and call update() in a "      \
      "loop, or use fast_evaluate for a vectorized batch.")                                                                                          \
    X(SVDSBTL_SAMPLING_THREADS, "SVDSBTL_SAMPLING_THREADS", static_cast<int>(1),                                                                     \
      "Number of worker threads for SVDSBTL table-build sampling.  1 (default) = serial (no extra threads).  N > 1 = use N worker threads, each "    \
      "with its own source AbstractState.  0 = auto (use std::thread::hardware_concurrency()).  Typical 4-8x build-time speedup at N >= 4 on a "     \
      "multi-core machine.  Default is 1 because: (a) REFPROP is process-global and not thread-safe under SETUPdll (the parallel path falls back "   \
      "to serial when the source is REFPROP regardless of this setting); (b) per-worker source instantiation roughly multiplies peak build-time "    \
      "memory footprint by N, which can matter on capped CI containers.  Set > 1 when build cost on a fresh ~/.CoolProp/SVDTables/ cache or "        \
      "post-rev-bump rebuild becomes the bottleneck.")                                                                                               \
    X(TABULAR_NX, "TABULAR_NX", static_cast<int>(200),                                                                                               \
      "Number of x-axis grid points (T for PT table, h for PH table) for the BICUBIC and TTSE tabular backends. Increase for higher accuracy in "    \
      "regions with steep gradients (e.g. near the critical point). Memory and build cost scale as O(Nx*Ny). Tables auto-rebuild when changed.")     \
    X(MIXTURE_STABILITY_ALGORITHM, "MIXTURE_STABILITY_ALGORITHM", 1, "0: legacy, 1: Michelsen (default)")                                            \
    X(TABULAR_NY, "TABULAR_NY", static_cast<int>(200),                                                                                               \
      "Number of y-axis grid points (log P) for the BICUBIC and TTSE tabular backends. Increase for higher accuracy. Memory and build cost scale "   \
      "as O(Nx*Ny). Tables auto-rebuild when changed.")

// Use preprocessor to create the Enum
enum configuration_keys
{
#define X(Enum, String, Default, Desc) Enum,
    CONFIGURATION_KEYS_ENUM
#undef X
};

#endif  // COOLPROP_DETAIL_CONFIGURATION_KEYS_H
