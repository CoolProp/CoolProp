fn main() {
    tauri_build::build();

    // Build CoolProp as a static library from the monorepo root.
    // src-tauri/ is three levels below the CoolProp root.
    let coolprop_root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..");

    let mut cfg = cmake::Config::new(&coolprop_root);

    // Pass through CPM dependency-source overrides (CPM_<Name>_SOURCE) so the
    // CoolProp build works in network-restricted environments where the
    // FetchContent tarball downloads are blocked and local clones are used
    // instead. No effect when the variables are unset.
    for (key, value) in std::env::vars() {
        if key.starts_with("CPM_") && key.ends_with("_SOURCE") {
            cfg.define(&key, &value);
        }
    }

    let dst = cfg
        .define("COOLPROP_STATIC_LIBRARY", "ON")
        .define("COOLPROP_SHARED_LIBRARY", "OFF")
        .define("BUILD_TESTING", "OFF")
        .define("COOLPROP_CATCH_MODULE", "OFF")
        .define("CMAKE_BUILD_TYPE", "Release")
        // COOLPROP_LIB makes CoolPropLib.h wrap all exported functions with
        // extern "C", giving C linkage symbols that the Rust FFI can call.
        .define("CMAKE_CXX_FLAGS", "-DCOOLPROP_LIB")
        // The Rust crate-type includes "cdylib" (for mobile targets), so every
        // object linked into it must be position-independent. Without -fPIC
        // here, Linux rust-lld rejects refs to std::bad_alloc, stderr, etc.
        .define("CMAKE_POSITION_INDEPENDENT_CODE", "ON")
        // std::filesystem (used in CPfilepaths.cpp) requires macOS 10.15+
        .define("CMAKE_OSX_DEPLOYMENT_TARGET", "10.15")
        .build();

    // CoolProp doesn't use cmake's install step; the lib stays in the build dir.
    // Single-config generators (Make/Ninja, used on Linux & macOS) put the lib
    // directly in build/. MSVC's multi-config generator nests it under Release/
    // (or Debug/), so add both search paths.
    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-search=native={}/build/Release", dst.display());
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
