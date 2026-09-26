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
# lives.
#
# The result is stored as an INTERNAL cache entry, not a normal variable.  A
# normal variable is visible only in the directory scope that include()d this
# module and its children, but the function below is global and may be CALLED
# from a parent scope: the standalone wrappers/Python/nanobind project does
# add_subdirectory(<repo root>), so this module is included in the child scope
# while coolprop_hide_json_symbols(CoolProp) runs in the parent, where a normal
# variable would be empty and the hide-list path would resolve to
# "/coolprop_hide_json.exp".  INTERNAL implies FORCE, so the entry is
# re-derived on every configure.
get_filename_component(_coolprop_json_linker_dir "${CMAKE_CURRENT_LIST_DIR}/../dev/linker" ABSOLUTE)
set(_COOLPROP_JSON_LINKER_DIR "${_coolprop_json_linker_dir}"
    CACHE INTERNAL "Directory holding the coolprop_hide_json link-time hide-lists")
unset(_coolprop_json_linker_dir)

function(coolprop_hide_json_symbols target)
    # Fail closed.  Every call site names a target it has just created, so a
    # missing target means a rename (e.g. UseSWIG real-name drift) has
    # detached the call from its product, and that product would ship with
    # nlohmann/valijson symbols exported.  Most wrappers have no nm gate in CI
    # to catch that, so it must not be skipped silently.
    if(NOT TARGET ${target})
        message(FATAL_ERROR "coolprop_hide_json_symbols: no target named '${target}'; "
                            "call it after the target is created, with the target's real name")
    endif()
    if(NOT _COOLPROP_JSON_LINKER_DIR)
        message(FATAL_ERROR "coolprop_hide_json_symbols: _COOLPROP_JSON_LINKER_DIR is empty; "
                            "include(CoolPropJSONVisibility) before calling it")
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
