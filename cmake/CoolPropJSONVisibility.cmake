# Hide nlohmann/valijson symbols from a shared product's dynamic export table at
# LINK time (replaces the compile-time visibility pragma; see CoolProp-xa8w.6).
# ELF: --version-script hide-list; Mach-O: -unexported_symbols_list. MSVC: no-op
# (exports are opt-in via src/CoolPropLib.def).
#
# Paths are resolved relative to THIS module file so the helper works whether
# included from the repo-root CMakeLists or the Python wrapper's own cmake
# project (where CMAKE_SOURCE_DIR is the wrapper dir, not the repo root).
# CMAKE_CURRENT_LIST_DIR evaluates at parse time here (outside the function
# body), so it correctly captures <repo>/cmake regardless of where the caller
# lives.  The function captures it via _COOLPROP_JSON_LINKER_DIR.
get_filename_component(_COOLPROP_JSON_LINKER_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../dev/linker" ABSOLUTE)

function(coolprop_hide_json_symbols target)
    if(NOT TARGET ${target})
        return()
    endif()
    if(APPLE)
        set(_exp "${_COOLPROP_JSON_LINKER_DIR}/coolprop_hide_json.exp")
        if(NOT EXISTS "${_exp}")
            message(FATAL_ERROR "coolprop_hide_json_symbols: hide-list not found: ${_exp}")
        endif()
        target_link_options(${target} PRIVATE "LINKER:-unexported_symbols_list,${_exp}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_exp}")
    elseif(UNIX)
        set(_map "${_COOLPROP_JSON_LINKER_DIR}/coolprop_hide_json.map")
        if(NOT EXISTS "${_map}")
            message(FATAL_ERROR "coolprop_hide_json_symbols: hide-list not found: ${_map}")
        endif()
        target_link_options(${target} PRIVATE "LINKER:--version-script=${_map}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_map}")
    endif()
    # MSVC/WIN32: nothing.
endfunction()
