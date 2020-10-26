use std::error;
use std::fmt;
use std::panic;

mod bindings;

// macro_rules! c_str {
//     ($s:expr) => {{
//         concat!($s, "\0").as_ptr() as *const i8
//     }};
// }

#[derive(Debug)]
pub struct CoolPropError;

// pub enum CoolPropError {
//     Error
// }

impl error::Error for CoolPropError {
    // fn description(&self) -> &str {
    //     // Both underlying errors already impl `Error`, so we defer to their
    //     // implementations.
    //     match *self {
    //         CliError::Io(ref err) => err.description(),
    //         // Normally we can just write `err.description()`, but the error
    //         // type has a concrete method called `description`, which conflicts
    //         // with the trait method. For now, we must explicitly call
    //         // `description` through the `Error` trait.
    //         CliError::Parse(ref err) => error::Error::description(err),
    //     }
    // }

    // fn cause(&self) -> Option<&error::Error> {
    //     match *self {
    //         // N.B. Both of these implicitly cast `err` from their concrete
    //         // types (either `&io::Error` or `&num::ParseIntError`)
    //         // to a trait object `&Error`. This works because both error types
    //         // implement `Error`.
    //         CliError::Io(ref err) => Some(err),
    //         CliError::Parse(ref err) => Some(err),
    //     }
    // }
}

impl fmt::Display for CoolPropError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            // Underlying errors already impl `Display`, so we defer to
            // their implementations.
            CoolPropError => write!(f, "CoolProp error: {}", CoolPropError),
        }
    }
}

pub fn props_si(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    refr: &str,
) -> Result<f64, CoolPropError> {
    let result = panic::catch_unwind(|| unsafe {
        bindings::PropsSI(
            format!("{}{}", output, "\0").as_ptr() as *const i8,
            format!("{}{}", name1, "\0").as_ptr() as *const i8,
            prop1,
            format!("{}{}", name2, "\0").as_ptr() as *const i8,
            prop2,
            format!("{}{}", refr, "\0").as_ptr() as *const i8,
        )
    });
    match result {
        Ok(result) => Ok(result),
        Err(_) => Err(CoolPropError),
    }
}

pub fn ha_props_si(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    name3: &str,
    prop3: f64,
) -> Result<f64, CoolPropError> {
    let result = panic::catch_unwind(|| unsafe {
        bindings::HAPropsSI(
            format!("{}{}", output, "\0").as_ptr() as *const i8,
            format!("{}{}", name1, "\0").as_ptr() as *const i8,
            prop1,
            format!("{}{}", name2, "\0").as_ptr() as *const i8,
            prop2,
            format!("{}{}", name3, "\0").as_ptr() as *const i8,
            prop3,
        )
    });
    match result {
        Ok(result) => Ok(result),
        Err(_) => Err(CoolPropError),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_props_functions() {
        assert_eq!(
            props_si("H", "T", 300.0, "Q", 1.0, "R134a").unwrap(),
            413265.6843372975
        );
        assert_eq!(
            ha_props_si("H", "T", 300.0, "P", 100000.0, "R", 0.0).unwrap(),
            27013.112479771713
        );
    }
}
