#######################################
#      COMMON WRAPPER UTILITIES       #
#-------------------------------------#
# Shared configuration and utilities  #
# for SWIG-based language wrappers    #
#######################################

# Common SWIG interface file
set(COOLPROP_SWIG_INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

# Common SWIG dependencies (headers that trigger SWIG regeneration)
# Already set in Sources.cmake, but repeated here for clarity
set(SWIG_DEPENDENCIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/AbstractState.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Configuration.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PhaseEnvelope.h)

# Disable internal error catching for SWIG wrappers (let language handle errors)
macro(disable_error_catching)
  add_definitions(-DNO_ERROR_CATCHING)
endmacro()

# Helper macro to set up SWIG module with common configuration
macro(setup_swig_module MODULE_NAME LANGUAGE EXTRA_FLAGS)
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  set(I_FILE "${COOLPROP_SWIG_INTERFACE}")
  set(SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_DEPENDENCIES})

  # Combine COOLPROP_SWIG_OPTIONS (from command line) with wrapper-specific flags
  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}" "${EXTRA_FLAGS}")
  string(REPLACE " " ";" SWIG_OPTIONS "${SWIG_OPTIONS}")

  set_source_files_properties(${I_FILE} PROPERTIES
                               SWIG_FLAGS "${SWIG_OPTIONS}"
                               CPLUSPLUS ON)
endmacro()
