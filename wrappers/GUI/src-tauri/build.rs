fn main() {
    tauri_build::build();

    // Build CoolProp as a static library from the monorepo root.
    // src-tauri/ is three levels below the CoolProp root.
    let coolprop_root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..");

    let dst = cmake::Config::new(&coolprop_root)
        .define("COOLPROP_STATIC_LIBRARY", "ON")
        .define("COOLPROP_SHARED_LIBRARY", "OFF")
        .define("BUILD_TESTING", "OFF")
        .define("COOLPROP_CATCH_MODULE", "OFF")
        .define("CMAKE_BUILD_TYPE", "Release")
        // COOLPROP_LIB makes CoolPropLib.h wrap all exported functions with
        // extern "C", giving C linkage symbols that the Rust FFI can call.
        .define("CMAKE_CXX_FLAGS", "-DCOOLPROP_LIB")
        // std::filesystem (used in CPfilepaths.cpp) requires macOS 10.15+
        .define("CMAKE_OSX_DEPLOYMENT_TARGET", "10.15")
        .build();

    // CoolProp doesn't use cmake's install step; the lib stays in the build dir.
    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=CoolProp");

    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap();
    match target_os.as_str() {
        "macos" => println!("cargo:rustc-link-lib=c++"),
        "linux" => {
            println!("cargo:rustc-link-lib=stdc++");
            println!("cargo:rustc-link-lib=m");
        }
        "windows" => {} // MSVC links C++ stdlib automatically
        _ => {}
    }
}
