local ffi        = require "ffi"
local ffi_new    = ffi.new
local ffi_str    = ffi.string
local ffi_load   = ffi.load
local ffi_cdef   = ffi.cdef
local ffi_typeof = ffi.typeof
local pcall      = pcall
local huge       = math.huge
ffi_cdef[[
double Props1SI(const char *FluidName, const char* Output);
double PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
long   PhaseSI(const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref, char *phase, int n);
long   get_global_param_string(const char *param, char *Output, int n);
long   get_parameter_information_string(const char *key, char *Output, int n);
long   get_mixture_binary_pair_data(const char *CAS1, const char *CAS2, const char *key);
long   get_fluid_param_string(const char *fluid, const char *param, char *Output, int n);
int    set_reference_stateS(const char *Ref, const char *reference_state);
int    set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
double F2K(double T_F);
double K2F(double T_K);
long   get_param_index(const char *param);
long   redirect_stdout(const char *file);
int    get_debug_level();
void   set_debug_level(int level);
double saturation_ancillary(const char *fluid_name, const char *output, int Q, const char *input, double value);
double HAPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Name3, double Prop3);
]]
local ok, newtab = pcall(require, "table.new")
if not ok then newtab = function() return {} end end
local but = ffi_typeof "char[?]"
local buf = ffi_new(but, 256)
local bub = ffi_new(but, 4096)
local ok, lib = pcall(ffi_load, "CoolProp")
if not ok then ok, lib = pcall(ffi_load, "coolprop") end
assert(ok, "Unable to load CoolProp. Please check that the CoolProp shared library is in a default search path for dynamic libraries of your operating system.")
local coolprop = newtab(0, 21)
function coolprop.Props1SI(fluidname, output)
    local v = lib.Props1SI(fluidname, output)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end
function coolprop.PropsSI(output, name1, prop1, name2, prop2, ref)
    local v = lib.PropsSI(output, name1, prop1, name2, prop2, ref)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end
function coolprop.PhaseSI(name1, prop1, name2, prop2, ref)
    if lib.PhaseSI(name1, prop1, name2, prop2, ref or "", buf, 256) == 1 then
        return ffi_str(buf)
    end
    return nil
end
function coolprop.get_global_param_string(param)
    if lib.get_global_param_string(param, bub, 4096) == 1 then
        return ffi_str(bub)
    end
    return nil
end
function coolprop.get_parameter_information_string(key)
    local value = lib.get_parameter_information_string(key, bub, 4096)
    if value == 0 then return nil end
    return ffi_str(bub)
end
function coolprop.get_mixture_binary_pair_data(cas1, cas2, key)
    return tonumber(lib.get_mixture_binary_pair_data(cas1, cas2, key))
end
function coolprop.get_fluid_param_string(fluid, param)
    if lib.get_fluid_param_string(fluid, param, bub, 4096) == 1 then
        return ffi_str(bub)
    end
    return nil
end
function coolprop.set_reference_stateS(ref, state)
    return lib.set_reference_stateS(ref, state) == 1
end
function coolprop.set_reference_stateD(ref, t, rho, h0, s0)
    return lib.set_reference_stateD(ref, t, rho, h0, s0) == 1
end
function coolprop.F2K(f)
    return lib.F2K(f)
end
function coolprop.K2F(k)
    return lib.K2F(k)
end
function coolprop.get_param_index(param)
    local value = lib.get_param_index(param)
    if value == -1 then return nil end
    tonumber(value)
end
function coolprop.saturation_ancillary(fluid, output, q, input, value)
    local v = lib.saturation_ancillary(fluid, output, q, input, value)
    if v == huge then
        return nil, coolprop.error()
    end
    return v
end
function coolprop.redirect_stdout(file)
    return lib.redirect_stdout(file) == 1
end
function coolprop.get_debug_level()
    return lib.get_debug_level()
end
function coolprop.set_debug_level(level)
    lib.set_debug_level(level)
end
function coolprop.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)
    local v = lib.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)
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
