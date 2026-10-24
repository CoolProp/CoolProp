#ifndef COOLPROP_CUBIC_CONSTANTS_H
#define COOLPROP_CUBIC_CONSTANTS_H

// The critical-point constants of the SRK and Peng-Robinson equations of state, computed from their
// closed forms (Bell and Deiters, Ind. Eng. Chem. Res., 2021) rather than written out as digits: #2745
// transcribed them with doubled and dropped digits, which put SRK's own critical point ~1e-9 below the
// input Tc and pc.  std::cbrt and std::sqrt are not constexpr before C++26, so each value is computed
// once, on first use (a function-local static, which also avoids static-initialisation-order issues).
// The test "Cubic Omega constants reproduce the critical point" checks the defining property.

#include <cmath>

namespace CoolProp {
namespace CubicConstants {

/// SRK: Omega_a = 1/(9 (2^(1/3) - 1)) = 0.42748023354034140...
inline double SRK_Omega_a() {
    static const double value = 1.0 / (9.0 * (std::cbrt(2.0) - 1.0));
    return value;
}

/// SRK: Omega_b = (2^(1/3) - 1)/3 = 0.08664034996495772...
inline double SRK_Omega_b() {
    static const double value = (std::cbrt(2.0) - 1.0) / 3.0;
    return value;
}

/// SRK: Z_c = 1/3 exactly.
inline double SRK_Zc() {
    return 1.0 / 3.0;
}

namespace detail {
struct PRCritical
{
    double Omega_a, Omega_b, Zc;
};
/// Peng-Robinson critical constants.  With eta = b/v_c = 1/(1 + (4 - sqrt 8)^(1/3) + (4 + sqrt 8)^(1/3)),
/// and in units where b = 1 and RT = 1 (so v_c = 1/eta): dp/dv = 0 gives a = g^2/((v - 1)^2 (2v + 2)) with
/// g = v^2 + 2v - 1, then p_c = 1/(v - 1) - a/g, Omega_b = p_c, Omega_a = a p_c and Z_c = p_c v_c.
inline const PRCritical& PR_critical() {
    static const PRCritical value = [] {
        const double s8 = std::sqrt(8.0);
        const double eta = 1.0 / (1.0 + std::cbrt(4.0 - s8) + std::cbrt(4.0 + s8));
        const double v = 1.0 / eta;
        const double g = v * v + 2.0 * v - 1.0;
        const double a = g * g / ((v - 1.0) * (v - 1.0) * (2.0 * v + 2.0));
        const double p = 1.0 / (v - 1.0) - a / g;
        return PRCritical{a * p, p, p * v};
    }();
    return value;
}
}  // namespace detail

/// Peng-Robinson: Omega_a = 0.45723552892138218...
inline double PR_Omega_a() {
    return detail::PR_critical().Omega_a;
}

/// Peng-Robinson: Omega_b = 0.07779607390388845...
inline double PR_Omega_b() {
    return detail::PR_critical().Omega_b;
}

/// Peng-Robinson: Z_c = 0.30740130869870384...
inline double PR_Zc() {
    return detail::PR_critical().Zc;
}

}  // namespace CubicConstants
}  // namespace CoolProp

#endif
