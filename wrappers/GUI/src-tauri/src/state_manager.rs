use std::collections::HashMap;
use std::os::raw::c_long;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use once_cell::sync::Lazy;
use serde::Serialize;

use crate::coolprop_ffi;

// ---------- global handle table ------------------------------------------

static STATES: Lazy<Mutex<HashMap<u64, c_long>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

static NEXT_ID: AtomicU64 = AtomicU64::new(1);

// ---------- basic AbstractState commands ---------------------------------

#[tauri::command]
pub fn create_state(backend: String, fluid: String) -> Result<u64, String> {
    let handle = coolprop_ffi::cp_create_state(&backend, &fluid)?;
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    STATES.lock().unwrap().insert(id, handle);
    Ok(id)
}

#[tauri::command]
pub fn update_state(id: u64, input_pair: String, v1: f64, v2: f64) -> Result<(), String> {
    let handle = *STATES
        .lock()
        .unwrap()
        .get(&id)
        .ok_or_else(|| format!("Invalid state id: {id}"))?;
    coolprop_ffi::cp_update_state(handle, &input_pair, v1, v2)
}

#[tauri::command]
pub fn get_property(id: u64, param: String) -> Result<f64, String> {
    let handle = *STATES
        .lock()
        .unwrap()
        .get(&id)
        .ok_or_else(|| format!("Invalid state id: {id}"))?;
    coolprop_ffi::cp_keyed_output(handle, &param)
}

#[tauri::command]
pub fn free_state(id: u64) -> Result<(), String> {
    let handle = STATES
        .lock()
        .unwrap()
        .remove(&id)
        .ok_or_else(|| format!("Invalid state id: {id}"))?;
    coolprop_ffi::cp_free_state(handle)
}

// ---------- fluid list ---------------------------------------------------

#[tauri::command]
pub fn get_fluids_list() -> Result<Vec<String>, String> {
    let raw = coolprop_ffi::cp_get_global_param("FluidsList")?;
    let mut list: Vec<String> = raw
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    list.sort_unstable();
    Ok(list)
}

// ---------- fluid EOS limits ---------------------------------------------

#[derive(Serialize)]
pub struct FluidLimits {
    pub t_min: f64, // triple-point temperature (K)  — lower bound for T sweep
    pub t_max: f64, // critical temperature (K)       — upper bound for T sweep
    pub p_min: f64, // triple-point pressure (Pa)     — lower bound for P sweep
    pub p_max: f64, // critical pressure (Pa)         — upper bound for P sweep
}

#[tauri::command]
pub fn get_fluid_limits(fluid: String) -> Result<FluidLimits, String> {
    Ok(FluidLimits {
        t_min: coolprop_ffi::cp_props1si(&fluid, "T_triple")?,
        t_max: coolprop_ffi::cp_props1si(&fluid, "T_critical")?,
        p_min: coolprop_ffi::cp_props1si(&fluid, "p_triple")?,
        p_max: coolprop_ffi::cp_props1si(&fluid, "p_critical")?,
    })
}

// ---------- saturation table (batch, single IPC round-trip) --------------

#[derive(Serialize)]
pub struct SatPoint {
    pub independent: f64,
    pub liquid: HashMap<String, f64>,
    pub vapor: HashMap<String, f64>,
}

fn collect_props(
    handle: c_long,
    pair: &str,
    v1: f64,
    v2: f64,
    params: &[String],
) -> HashMap<String, f64> {
    let mut m = HashMap::new();
    if coolprop_ffi::cp_update_state(handle, pair, v1, v2).is_ok() {
        for p in params {
            if let Ok(v) = coolprop_ffi::cp_keyed_output(handle, p) {
                m.insert(p.clone(), v);
            }
        }
    }
    m
}

// ---------- humid air (HAPropsSI — stateless) ----------------------------

/// Compute multiple humid air outputs for a single state point.
/// Always supply P + two independent inputs (e.g. T + R).
/// `outputs`: list of HAPropsSI output names (e.g. ["W","H","S","V","C","M","K"]).
#[tauri::command]
pub fn compute_humid_air(
    n1: String, v1: f64,
    n2: String, v2: f64,
    n3: String, v3: f64,
    outputs: Vec<String>,
) -> Result<HashMap<String, f64>, String> {
    // Validate inputs with a probe output that cannot be trivially shortcut
    // (HAPropsSI returns an input value unchanged if the requested output matches
    // one of the inputs, bypassing all validation). Pick the first requested output
    // that is not one of the three input names.
    let probe = outputs.iter()
        .find(|o| *o != &n1 && *o != &n2 && *o != &n3)
        .cloned()
        .unwrap_or_else(|| "M".to_string()); // M (viscosity) is never a valid HA input
    if let Err(e) = coolprop_ffi::cp_ha_props(&probe, &n1, v1, &n2, v2, &n3, v3) {
        return Err(e);
    }

    let mut result = HashMap::new();
    for out in &outputs {
        if let Ok(v) = coolprop_ffi::cp_ha_props(out, &n1, v1, &n2, v2, &n3, v3) {
            result.insert(out.clone(), v);
        }
    }
    Ok(result)
}

// ---------- saturation table (batch, single IPC round-trip) --------------

/// Compute a saturation table using an existing AbstractState (by id).
///
/// `by_temperature`: if true, sweep T in [min_val, max_val] using QT_INPUTS;
///                   if false, sweep P in [min_val, max_val] using PQ_INPUTS.
/// `params`: list of CoolProp parameter names to return (T and P are always added).
#[tauri::command]
pub fn compute_saturation_table(
    id: u64,
    by_temperature: bool,
    min_val: f64,
    max_val: f64,
    n_points: usize,
    params: Vec<String>,
) -> Result<Vec<SatPoint>, String> {
    if n_points < 2 {
        return Err("n_points must be >= 2".to_string());
    }

    let handle = *STATES
        .lock()
        .unwrap()
        .get(&id)
        .ok_or_else(|| format!("Invalid state id: {id}"))?;

    // Always include T and P
    let mut all_params = vec!["T".to_string(), "P".to_string()];
    for p in &params {
        if p != "T" && p != "P" {
            all_params.push(p.clone());
        }
    }

    let mut points = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;
        let val = min_val + t * (max_val - min_val);

        let (pair, liq_v1, liq_v2, vap_v1, vap_v2) = if by_temperature {
            ("QT_INPUTS", 0.0f64, val, 1.0f64, val)
        } else {
            ("PQ_INPUTS", val, 0.0f64, val, 1.0f64)
        };

        let liquid = collect_props(handle, pair, liq_v1, liq_v2, &all_params);
        let vapor  = collect_props(handle, pair, vap_v1, vap_v2, &all_params);

        if !liquid.is_empty() || !vapor.is_empty() {
            points.push(SatPoint { independent: val, liquid, vapor });
        }
    }

    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fluids_list_command_returns_water() {
        let list = get_fluids_list().expect("FluidsList");
        assert!(list.iter().any(|f| f == "Water"), "FluidsList missing Water");
        assert!(list.windows(2).all(|w| w[0] <= w[1]), "FluidsList not sorted");
    }

    #[test]
    fn fluid_limits_water() {
        let lim = get_fluid_limits("Water".to_string()).expect("Water limits");
        // IAPWS-95 references: T_triple=273.16, T_crit=647.096, P_crit=22.064 MPa.
        assert!((lim.t_min - 273.16).abs() < 0.05, "t_min = {}", lim.t_min);
        assert!((lim.t_max - 647.096).abs() < 0.05, "t_max = {}", lim.t_max);
        assert!(lim.p_min > 0.0 && lim.p_min < 1000.0, "p_min = {}", lim.p_min);
        assert!(lim.p_max > 1e6 && lim.p_max < 1e9, "p_max = {}", lim.p_max);
    }

    #[test]
    fn create_update_get_free_lifecycle() {
        let id = create_state("HEOS".into(), "Water".into()).expect("create_state");
        update_state(id, "PT_INPUTS".into(), 101325.0, 300.0).expect("update_state");
        let rho = get_property(id, "Dmass".into()).expect("get_property");
        assert!((rho - 996.0).abs() < 2.0, "rho = {rho}");
        free_state(id).expect("free_state");
        // After free, the id should be invalid.
        let r = get_property(id, "Dmass".into());
        assert!(r.is_err(), "Expected error after free, got {r:?}");
    }

    #[test]
    fn invalid_state_id_errors() {
        let r = get_property(99_999_999, "Dmass".into());
        assert!(r.is_err(), "Expected error for invalid id");
    }

    #[test]
    fn saturation_table_too_few_points_errors() {
        let id = create_state("HEOS".into(), "Water".into()).unwrap();
        let r = compute_saturation_table(id, true, 280.0, 600.0, 1, vec![]);
        free_state(id).ok();
        assert!(r.is_err(), "Expected error for n_points < 2");
    }

    #[test]
    fn saturation_table_water_t_sweep() {
        let id = create_state("HEOS".into(), "Water".into()).unwrap();
        let pts = compute_saturation_table(
            id,
            true, // by temperature
            300.0,
            500.0,
            10,
            vec!["Dmass".into(), "Hmass".into()],
        )
        .expect("sat table");
        free_state(id).ok();
        assert!(!pts.is_empty(), "no points returned");
        // Each point should have liquid + vapor with T, P, Dmass, Hmass.
        for p in &pts {
            for side in [&p.liquid, &p.vapor] {
                assert!(side.contains_key("T"));
                assert!(side.contains_key("P"));
                assert!(side.contains_key("Dmass"));
                assert!(side.contains_key("Hmass"));
            }
            // Liquid density >> vapor density on the saturation curve.
            assert!(p.liquid["Dmass"] > p.vapor["Dmass"]);
        }
    }

    #[test]
    fn humid_air_command_basic_outputs() {
        let res = compute_humid_air(
            "P".into(), 101325.0,
            "T".into(), 293.15,
            "R".into(), 0.5,
            vec!["W".into(), "H".into(), "B".into()],
        ).expect("compute_humid_air");
        let w = res["W"];
        assert!((w - 0.00729).abs() < 1e-3, "W = {w}");
        assert!(res.contains_key("H"));
        assert!(res.contains_key("B"));
    }

    #[test]
    fn humid_air_command_invalid_input_errors() {
        let r = compute_humid_air(
            "P".into(), 101325.0,
            "T".into(), 293.15,
            "NOT_AN_INPUT".into(), 0.5,
            vec!["W".into()],
        );
        assert!(r.is_err(), "Expected error for invalid input name");
    }
}
