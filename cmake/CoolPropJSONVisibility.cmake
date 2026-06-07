# Hide nlohmann/valijson symbols from a shared product's dynamic export table at
# LINK time (replaces the compile-time visibility pragma; see CoolProp-xa8w.6).
# ELF: --version-script hide-list; Mach-O: -unexported_symbols_list. MSVC: no-op
# (exports are opt-in via src/CoolPropLib.def).
function(coolprop_hide_json_symbols target)
    if(NOT TARGET ${target})
        return()
    endif()
    if(APPLE)
        set(_exp "${CMAKE_SOURCE_DIR}/dev/linker/coolprop_hide_json.exp")
        target_link_options(${target} PRIVATE "LINKER:-unexported_symbols_list,${_exp}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_exp}")
    elseif(UNIX)
        set(_map "${CMAKE_SOURCE_DIR}/dev/linker/coolprop_hide_json.map")
        target_link_options(${target} PRIVATE "LINKER:--version-script=${_map}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_map}")
    endif()
    # MSVC/WIN32: nothing.
endfunction()
