mod coolprop_ffi;
mod state_manager;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            state_manager::create_state,
            state_manager::update_state,
            state_manager::get_property,
            state_manager::free_state,
            state_manager::get_fluids_list,
            state_manager::get_fluid_limits,
            state_manager::compute_saturation_table,
            state_manager::compute_humid_air,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
