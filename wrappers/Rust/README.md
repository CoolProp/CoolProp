# Welcome to coolprop-rs

**EXPERIMENTAL** CoolProp Wrapper for Rust

The wrapper uses rust-bindgen to create the bindings to the C++ shared library of CoolProp and adds some wrapper functions. The wrapper should work on all platforms that CoolProp and Rust work on.

## Prerequisites

The CoolProp shared library (_libCoolProp.so_) and header file (_CoolPropLib.h_) must be installed on the computer in the systems folder. On Linux e. g. _/usr/lib64_ and _/usr/include_. Instructions to compile and install CoolProp for your system can be found on the project page of [CoolProp](https://github.com/CoolProp/CoolProp).

## Installation

The wrapper gets published on [crates.io](https://crates.io/crates/coolprop-rs) as `coolplot-rs` and you can add the library in your project.

```toml
[dependencies]
coolprop-rs = "0.2"
```

## Examples

At the moment the wrapper provides access to either the full C++ bindings directly or a small subset of methods which utilizes Rust types, error handling and unit testing.

The C++ bindings can be used with:

Rust:

```Rust
use coolprop-rs::bindings::*;
```

The subset of Rust methods are at the moment **PropsSI()** and **HAPropsSI()**:

Rust:

```Rust
use coolprop-rs;
println!("{:?}", coolprop-rs::PropsSI("H", "T", 300.0, "Q", 1.0, "R134a").unwrap());
println!("{:?}", coolprop-rs::HAPropsSI("H", "T", 300.0, "P", 100000.0, "R", 0.0).unwrap());
```

Output:

```bash
413265.6843372975
27013.112479771713
```

## License

This Rust package is released under the terms of the MIT license.
