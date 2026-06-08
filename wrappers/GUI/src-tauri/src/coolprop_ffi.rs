use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_long};

const BUF: c_long = 8192;

extern "C" {
    fn AbstractState_factory(
        backend: *const c_char,
        fluid_names: *const c_char,
        errcode: *mut c_long,
        message_buffer: *mut c_char,
        buffer_size: c_long,
    ) -> c_long;

    fn AbstractState_free(
        handle: c_long,
        errcode: *mut c_long,
        message_buffer: *mut c_char,
        buffer_size: c_long,
    );

    fn AbstractState_update(
        handle: c_long,
        input_pair: c_long,
        value1: c_double,
        value2: c_double,
        errcode: *mut c_long,
        message_buffer: *mut c_char,
        buffer_size: c_long,
    );

    fn AbstractState_keyed_output(
        handle: c_long,
        param: c_long,
        errcode: *mut c_long,
        message_buffer: *mut c_char,
        buffer_size: c_long,
    ) -> c_double;

    fn get_param_index(param: *const c_char) -> c_long;
    fn get_input_pair_index(input_pair: *const c_char) -> c_long;

    fn Props1SI(fluid_name: *const c_char, output: *const c_char) -> c_double;

    fn get_global_param_string(
        param: *const c_char,
        result: *mut c_char,
        result_size: c_long,
        errcode: *mut c_long,
        message_buffer: *mut c_char,
        buffer_size: c_long,
    );

    fn HAPropsSI(
        output: *const c_char,
        name1: *const c_char,
        val1: c_double,
        name2: *const c_char,
        val2: c_double,
        name3: *const c_char,
        val3: c_double,
    ) -> c_double;
}

fn cstring(s: &str) -> Result<CString, String> {
    CString::new(s).map_err(|e| e.to_string())
}

fn check(errcode: c_long, buf: &[c_char]) -> Result<(), String> {
    if errcode != 0 {
        let msg = unsafe { CStr::from_ptr(buf.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        Err(msg)
    } else {
        Ok(())
    }
}

pub fn cp_create_state(backend: &str, fluid: &str) -> Result<c_long, String> {
    let b = cstring(backend)?;
    let f = cstring(fluid)?;
    let mut ec: c_long = 0;
    let mut buf = vec![0i8; BUF as usize];
    let handle = unsafe {
        AbstractState_factory(b.as_ptr(), f.as_ptr(), &mut ec, buf.as_mut_ptr(), BUF)
    };
    check(ec, &buf)?;
    Ok(handle)
}

pub fn cp_free_state(handle: c_long) -> Result<(), String> {
    let mut ec: c_long = 0;
    let mut buf = vec![0i8; BUF as usize];
    unsafe { AbstractState_free(handle, &mut ec, buf.as_mut_ptr(), BUF) };
    check(ec, &buf)
}

pub fn cp_update_state(handle: c_long, input_pair: &str, v1: f64, v2: f64) -> Result<(), String> {
    let ip_c = cstring(input_pair)?;
    let ip_idx = unsafe { get_input_pair_index(ip_c.as_ptr()) };
    if ip_idx == -1 {
        return Err(format!("Unknown input pair: {input_pair}"));
    }
    let mut ec: c_long = 0;
    let mut buf = vec![0i8; BUF as usize];
    unsafe { AbstractState_update(handle, ip_idx, v1, v2, &mut ec, buf.as_mut_ptr(), BUF) };
    check(ec, &buf)
}

pub fn cp_keyed_output(handle: c_long, param: &str) -> Result<f64, String> {
    let p_c = cstring(param)?;
    let p_idx = unsafe { get_param_index(p_c.as_ptr()) };
    if p_idx == -1 {
        return Err(format!("Unknown parameter: {param}"));
    }
    let mut ec: c_long = 0;
    let mut buf = vec![0i8; BUF as usize];
    let v = unsafe { AbstractState_keyed_output(handle, p_idx, &mut ec, buf.as_mut_ptr(), BUF) };
    check(ec, &buf)?;
    Ok(v)
}

/// Simple single-output lookup (no error code — returns 1e300 on failure).
pub fn cp_props1si(fluid: &str, output: &str) -> Result<f64, String> {
    let f = cstring(fluid)?;
    let o = cstring(output)?;
    let v = unsafe { Props1SI(f.as_ptr(), o.as_ptr()) };
    if v.is_finite() && v.abs() < 1e20 {
        Ok(v)
    } else {
        Err(format!("Props1SI({fluid}, {output}) returned no-data value"))
    }
}

/// Humid air property lookup via HAPropsSI (stateless, 3-input).
/// On failure, retrieves the actual error message from CoolProp's global errstring.
pub fn cp_ha_props(
    output: &str,
    n1: &str, v1: f64,
    n2: &str, v2: f64,
    n3: &str, v3: f64,
) -> Result<f64, String> {
    let o  = cstring(output)?;
    let c1 = cstring(n1)?;
    let c2 = cstring(n2)?;
    let c3 = cstring(n3)?;
    let v = unsafe { HAPropsSI(o.as_ptr(), c1.as_ptr(), v1, c2.as_ptr(), v2, c3.as_ptr(), v3) };
    if v.is_finite() && v.abs() < 1e20 {
        Ok(v)
    } else {
        let msg = cp_get_global_param("errstring")
            .unwrap_or_else(|_| format!("HAPropsSI({output}) returned no-data value"));
        Err(msg)
    }
}

pub fn cp_get_global_param(param: &str) -> Result<String, String> {
    let p_c = cstring(param)?;
    let result_size: c_long = 65536;
    let mut result = vec![0i8; result_size as usize];
    let mut ec: c_long = 0;
    let mut buf = vec![0i8; BUF as usize];
    unsafe {
        get_global_param_string(
            p_c.as_ptr(),
            result.as_mut_ptr(),
            result_size,
            &mut ec,
            buf.as_mut_ptr(),
            BUF,
        )
    };
    check(ec, &buf)?;
    Ok(unsafe { CStr::from_ptr(result.as_ptr()) }
        .to_string_lossy()
        .into_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fluids_list_is_non_empty() {
        let s = cp_get_global_param("FluidsList").expect("FluidsList lookup failed");
        assert!(s.contains("Water"), "FluidsList missing Water: {s:?}");
        assert!(s.split(',').count() > 50, "FluidsList suspiciously short");
    }

    #[test]
    fn props1si_water_critical_point() {
        let tcrit = cp_props1si("Water", "Tcrit").expect("Tcrit lookup failed");
        // IAPWS-95: 647.096 K
        assert!((tcrit - 647.096).abs() < 0.01, "Tcrit = {tcrit}, expected ~647.096");
    }

    #[test]
    fn props1si_invalid_fluid_returns_error() {
        let r = cp_props1si("ThisFluidDoesNotExist", "Tcrit");
        assert!(r.is_err(), "Expected error for unknown fluid");
    }

    #[test]
    fn create_state_invalid_backend_errors() {
        let r = cp_create_state("NOT_A_BACKEND", "Water");
        assert!(r.is_err(), "Expected error for invalid backend, got {r:?}");
    }

    #[test]
    fn create_state_invalid_fluid_errors() {
        let r = cp_create_state("HEOS", "ThisFluidDoesNotExist");
        assert!(r.is_err(), "Expected error for unknown fluid, got {r:?}");
    }

    #[test]
    fn full_state_lifecycle_water_pt() {
        let h = cp_create_state("HEOS", "Water").expect("create_state");
        cp_update_state(h, "PT_INPUTS", 101325.0, 300.0).expect("update_state");
        let rho = cp_keyed_output(h, "Dmass").expect("keyed_output Dmass");
        // Density of liquid water at 1 atm, 300 K — 996 kg/m³.
        assert!((rho - 996.0).abs() < 2.0, "Dmass = {rho}, expected ~996");
        cp_free_state(h).expect("free_state");
    }

    #[test]
    fn unknown_input_pair_errors() {
        let h = cp_create_state("HEOS", "Water").expect("create_state");
        let r = cp_update_state(h, "NOT_A_PAIR", 1.0, 1.0);
        cp_free_state(h).ok();
        assert!(r.is_err(), "Expected error for unknown input pair");
    }

    #[test]
    fn unknown_output_param_errors() {
        let h = cp_create_state("HEOS", "Water").expect("create_state");
        cp_update_state(h, "PT_INPUTS", 101325.0, 300.0).expect("update_state");
        let r = cp_keyed_output(h, "NOT_A_PARAM");
        cp_free_state(h).ok();
        assert!(r.is_err(), "Expected error for unknown parameter");
    }

    #[test]
    fn happrops_known_values() {
        // Dry-bulb 20°C, RH 50%, P=101325 Pa → W ≈ 0.0073 kg/kg.
        let w = cp_ha_props("W", "T", 293.15, "R", 0.5, "P", 101325.0)
            .expect("HAPropsSI W");
        assert!((w - 0.00729).abs() < 1e-3, "W = {w}, expected ~0.0073");
    }
}
