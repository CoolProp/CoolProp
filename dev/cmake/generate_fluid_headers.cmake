
#if(NOT COOLPROP_DEBUG_CMAKE)
#  message(FATAL_ERROR "Debug option missing, aborting.")
#else
if(NOT COOLPROP_SOURCES_ROOT)
  message(FATAL_ERROR "Root directory missing, aborting.")
endif()

find_package(PythonInterp 2.7 REQUIRED QUIET)
set(COOLPROP_HEADER_SOURCES_PREFIX "${COOLPROP_SOURCES_ROOT}/dev")
set(COOLPROP_HEADER_TARGETS_PREFIX "${COOLPROP_SOURCES_ROOT}/include")
if(NOT COOLPROP_DEBUG_CMAKE)
  set(COOLPROP_VARIABLE_OPTIONS_TMP "QUIET")
else()
  set(COOLPROP_VARIABLE_OPTIONS_TMP "DEBUG")
endif()
macro(add_header_generator COOLPROP_SOURCE_FILE_TMP COOLPROP_TARGET_FILE_TMP COOLPROP_VARIABLE_NAME_TMP)
  #add_custom_command(
  #  OUTPUT "${COOLPROP_HEADER_TARGETS_PREFIX}/${COOLPROP_TARGET_FILE_TMP}"
  #  COMMAND "${PYTHON_EXECUTABLE}" "${COOLPROP_HEADER_SOURCES_PREFIX}/generate_headers.py" "${COOLPROP_HEADER_SOURCES_PREFIX}/${COOLPROP_SOURCE_FILE_TMP}" "${COOLPROP_HEADER_TARGETS_PREFIX}/${COOLPROP_TARGET_FILE_TMP}" "${COOLPROP_VARIABLE_NAME_TMP}" "${COOLPROP_VARIABLE_OPTIONS_TMP}"
  #  DEPENDS "${COOLPROP_HEADER_SOURCES_PREFIX}/generate_headers.py" "${COOLPROP_HEADER_SOURCES_PREFIX}/${COOLPROP_SOURCE_FILE_TMP}"
  #)
  #add_dependencies(generate_fluid_headers "${COOLPROP_HEADER_TARGETS_PREFIX}/${COOLPROP_TARGET_FILE_TMP}")
  execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" "${COOLPROP_HEADER_SOURCES_PREFIX}/generate_headers.py" "${COOLPROP_HEADER_SOURCES_PREFIX}/${COOLPROP_SOURCE_FILE_TMP}" "${COOLPROP_HEADER_TARGETS_PREFIX}/${COOLPROP_TARGET_FILE_TMP}" "${COOLPROP_VARIABLE_NAME_TMP}" "${COOLPROP_VARIABLE_OPTIONS_TMP}"
  )
endmacro()
add_header_generator("all_fluids.json" "all_fluids_JSON.h" "all_fluids_JSON")
add_header_generator("all_incompressibles.json" "all_incompressibles_JSON.h" "all_incompressibles_JSON")
add_header_generator("mixtures/mixture_departure_functions.json" "mixture_departure_functions_JSON.h" "mixture_departure_functions_JSON")
add_header_generator("mixtures/mixture_binary_pairs.json" "mixture_binary_pairs_JSON.h" "mixture_binary_pairs_JSON")
add_header_generator("mixtures/predefined_mixtures.json" "predefined_mixtures_JSON.h" "predefined_mixtures_JSON")
add_header_generator("cubics/all_cubic_fluids.json" "all_cubics_JSON.h" "all_cubics_JSON")
add_header_generator("cubics/cubic_fluids_schema.json" "cubic_fluids_schema_JSON.h" "cubic_fluids_schema_JSON")