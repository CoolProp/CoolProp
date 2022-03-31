#include <stddef.h>
#include <stdlib.h>
#include <lua.h>
#include <lauxlib.h>
#include "../../../include/CoolPropLib.h"
#include "compat.h"

static int lua_coolprop_Props1SI(lua_State* L) {
    lua_pushnumber(L, Props1SI(luaL_checkstring(L, 1), luaL_checkstring(L, 2)));
    return 1;
}

static int lua_coolprop_PropsSI(lua_State* L) {
    lua_pushnumber(L, PropsSI(luaL_checkstring(L, 1), luaL_checkstring(L, 2), luaL_checknumber(L, 3), luaL_checkstring(L, 4), luaL_checknumber(L, 5),
                              luaL_checkstring(L, 6)));
    return 1;
}

static int lua_coolprop_PhaseSI(lua_State* L) {
    char phase[256];
    int ok = PhaseSI(luaL_checkstring(L, 1), luaL_checknumber(L, 2), luaL_checkstring(L, 3), luaL_checknumber(L, 4), luaL_checkstring(L, 5), phase,
                     sizeof(phase));
    if (ok == 1) {
        lua_pushstring(L, phase);
        return 1;
    } else {
        return 0;
    }
}

static int lua_coolprop_get_global_param_string(lua_State* L) {
    char output[4096];
    int ok = get_global_param_string(luaL_checkstring(L, 1), output, sizeof(output));
    if (ok == 1) {
        lua_pushstring(L, output);
        return 1;
    }
    return 0;
}

static int lua_coolprop_get_parameter_information_string(lua_State* L) {
    char output[4096];
    int ok = get_parameter_information_string(luaL_checkstring(L, 1), output, sizeof(output));
    if (ok == 1) {
        lua_pushstring(L, output);
        return 1;
    }
    return 0;
}

/*
static int lua_coolprop_get_mixture_binary_pair_data(lua_State *L) {
    lua_pushinteger(L, get_mixture_binary_pair_data(
        luaL_checkstring(L, 1),
        luaL_checkstring(L, 2),
        luaL_checkstring(L, 3)
    ));
    return 1;
}
*/

static int lua_coolprop_get_fluid_param_string(lua_State* L) {
    char output[4096];
    int ok = get_fluid_param_string(luaL_checkstring(L, 1), luaL_checkstring(L, 2), output, sizeof(output));
    if (ok == 1) {
        lua_pushstring(L, output);
        return 1;
    }
    return 0;
}

static int lua_coolprop_set_reference_stateS(lua_State* L) {
    lua_pushboolean(L, set_reference_stateS(luaL_checkstring(L, 1), luaL_checkstring(L, 2)));
    return 1;
}

static int lua_coolprop_set_reference_stateD(lua_State* L) {
    lua_pushboolean(L, set_reference_stateD(luaL_checkstring(L, 1), luaL_checknumber(L, 2), luaL_checknumber(L, 3), luaL_checknumber(L, 4),
                                            luaL_checknumber(L, 5)));
    return 1;
}

static int lua_coolprop_F2K(lua_State* L) {
    lua_pushnumber(L, F2K(luaL_checknumber(L, 1)));
    return 1;
}

static int lua_coolprop_K2F(lua_State* L) {
    lua_pushnumber(L, K2F(luaL_checknumber(L, 1)));
    return 1;
}

static int lua_coolprop_get_param_index(lua_State* L) {
    lua_pushnumber(L, get_param_index(luaL_checkstring(L, 1)));
    return 1;
}

static int lua_coolprop_redirect_stdout(lua_State* L) {
    lua_pushboolean(L, redirect_stdout(luaL_checkstring(L, 1)));
    return 1;
}

static int lua_coolprop_get_debug_level(lua_State* L) {
    lua_pushinteger(L, get_debug_level());
    return 1;
}

static int lua_coolprop_set_debug_level(lua_State* L) {
    set_debug_level(luaL_checknumber(L, 1));
    return 0;
}

static int lua_coolprop_saturation_ancillary(lua_State* L) {
    lua_pushnumber(L, saturation_ancillary(luaL_checkstring(L, 1), luaL_checkstring(L, 2), luaL_checkinteger(L, 3), luaL_checkstring(L, 4),
                                           luaL_checknumber(L, 5)));
    return 1;
}

static int lua_coolprop_HAPropsSI(lua_State* L) {
    lua_pushnumber(L, HAPropsSI(luaL_checkstring(L, 1), luaL_checkstring(L, 2), luaL_checknumber(L, 3), luaL_checkstring(L, 4),
                                luaL_checknumber(L, 5), luaL_checkstring(L, 6), luaL_checknumber(L, 7)));
    return 1;
}

static struct luaL_Reg funcs[] = {{"Props1SI", lua_coolprop_Props1SI},
                                  {"PropsSI", lua_coolprop_PropsSI},
                                  {"PhaseSI", lua_coolprop_PhaseSI},
                                  {"get_global_param_string", lua_coolprop_get_global_param_string},
                                  {"get_parameter_information_string", lua_coolprop_get_parameter_information_string},
                                  //{"get_mixture_binary_pair_data", lua_coolprop_get_mixture_binary_pair_data},
                                  {"get_fluid_param_string", lua_coolprop_get_fluid_param_string},
                                  {"set_reference_stateS", lua_coolprop_set_reference_stateS},
                                  {"set_reference_stateD", lua_coolprop_set_reference_stateD},
                                  {"F2K", lua_coolprop_F2K},
                                  {"K2F", lua_coolprop_K2F},
                                  {"get_param_index", lua_coolprop_get_param_index},
                                  {"redirect_stdout", lua_coolprop_redirect_stdout},
                                  {"get_debug_level", lua_coolprop_get_debug_level},
                                  {"set_debug_level", lua_coolprop_set_debug_level},
                                  {"saturation_ancillary", lua_coolprop_saturation_ancillary},
                                  {"HAPropsSI", lua_coolprop_HAPropsSI},
                                  {NULL, NULL}};

int luaopen_coolprop_capi(lua_State* L) {
    luaL_newlib(L, funcs);
    return 1;
}
