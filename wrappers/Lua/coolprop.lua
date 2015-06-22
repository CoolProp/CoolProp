
-- Try the FFI-enabled binding if running LuaJIT
if jit and jit.status() then
    local ok, mod = pcall(require, "coolprop.ffi")
    if ok then
        return mod
    end
    -- If loading via FFI fails, try Lua/C API,
    -- which might contain a statically-linked CoolProp.
end

local ok, capi = pcall(require, "coolprop.capi")
assert(ok, "Unable to load Lua CoolProp C API Wrapper.")

local huge = math.huge

local coolprop = {}

function coolprop.Props1SI(fluidname, output)
    local v = capi.Props1SI(fluidname, output)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end

function coolprop.PropsSI(output, name1, prop1, name2, prop2, ref)
    local v = capi.PropsSI(output, name1, prop1, name2, prop2, ref)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end

function coolprop.PhaseSI(name1, prop1, name2, prop2, ref)
    local phase, err = capi.PhaseSI(name1, prop1, name2, prop2, ref or "")
    if not phase then
        return nil, err or coolprop.error()
    end
    return phase
end

function coolprop.get_global_param_string(param)
    local v, err = capi.get_global_param_string(param)
    if not v then
        return nil, err or coolprop.error()
    end
    return v
end

function coolprop.get_parameter_information_string(key)
    local val, err = capi.get_parameter_information_string(key)
    if not val then
        return nil, err or coolprop.error()
    end
    return val
end

function coolprop.get_mixture_binary_pair_data(cas1, cas2, key)
    return tonumber(capi.get_mixture_binary_pair_data(cas1, cas2, key))
end

function coolprop.get_fluid_param_string(fluid, param)
    local val, err = capi.get_fluid_param_string(fluid, param)
    if not val then
        return nil, err or coolprop.error()
    end
    return val
end

function coolprop.set_reference_stateS(ref, state)
    return capi.set_reference_stateS(ref, state)
end

function coolprop.set_reference_stateD(ref, t, rho, h0, s0)
    return capi.set_reference_stateD(ref, t, rho, h0, s0)
end

function coolprop.F2K(f)
    return capi.F2K(f)
end

function coolprop.K2F(k)
    return capi.K2F(k)
end

function coolprop.get_param_index(param)
    local value = capi.get_param_index(param)
    if value == -1 then return nil end
    return tonumber(value)
end

function coolprop.saturation_ancillary(fluid, output, q, input, value)
    local v = capi.saturation_ancillary(fluid, output, q, input, value)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end

function coolprop.redirect_stdout(file)
    return capi.redirect_stdout(file)
end

function coolprop.get_debug_level()
    return capi.get_debug_level()
end

function coolprop.set_debug_level(level)
    capi.set_debug_level(level)
end

function coolprop.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)
    local v = capi.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end

function coolprop.error()
    return coolprop.get_global_param_string("errstring")
end

function coolprop.FluidsList()
    return coolprop.get_global_param_string("FluidsList")
end

function coolprop.version()
    return coolprop.get_global_param_string("version")
end

function coolprop.gitrevision()
    return coolprop.get_global_param_string("gitrevision")
end

return coolprop
