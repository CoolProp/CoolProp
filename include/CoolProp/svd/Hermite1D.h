#ifndef COOLPROP_SVD_HERMITE1D_H
#define COOLPROP_SVD_HERMITE1D_H

namespace CoolProp {
namespace svd {

// Cubic-Hermite interpolation on a single interval.
//
// Given values y0,y1 and slopes m0,m1 at the two ends of [x0, x1] with
// width h = x1 - x0 and t = (x - x0)/h in [0, 1], the interpolant is
//
//   p(t) = h00(t) y0 + h*h10(t) m0 + h01(t) y1 + h*h11(t) m1
//
// with the four Hermite basis polynomials
//
//   h00(t) =  2 t^3 - 3 t^2 + 1
//   h10(t) =      t^3 - 2 t^2 + t
//   h01(t) = -2 t^3 + 3 t^2
//   h11(t) =      t^3 -   t^2
//
// Slopes are an *input*: this header has no opinion on where they come
// from. Use natural-cubic-spline slopes, FD slopes, PCHIP slopes,
// analytic derivatives — anything that gives one number per knot.

struct HermiteBasis
{
    double h00;
    double h10;
    double h01;
    double h11;
};

inline HermiteBasis hermite_basis(double t) noexcept {
    const double t2 = t * t;
    const double t3 = t2 * t;
    return {2.0 * t3 - 3.0 * t2 + 1.0, t3 - 2.0 * t2 + t, -2.0 * t3 + 3.0 * t2, t3 - t2};
}

inline double hermite_eval(double y0, double y1, double m0, double m1, double h, double t) noexcept {
    const HermiteBasis b = hermite_basis(t);
    return b.h00 * y0 + b.h10 * h * m0 + b.h01 * y1 + b.h11 * h * m1;
}

// Derivative dp/dx at parameter t (note dt/dx = 1/h, so we divide by h).
inline double hermite_eval_deriv(double y0, double y1, double m0, double m1, double h, double t) noexcept {
    const double t2 = t * t;
    // d/dt of (h00, h10, h01, h11)
    const double dh00 = 6.0 * t2 - 6.0 * t;
    const double dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    const double dh01 = -6.0 * t2 + 6.0 * t;
    const double dh11 = 3.0 * t2 - 2.0 * t;
    const double dp_dt = dh00 * y0 + dh10 * h * m0 + dh01 * y1 + dh11 * h * m1;
    return dp_dt / h;
}

}  // namespace svd
}  // namespace CoolProp

#endif  // COOLPROP_SVD_HERMITE1D_H
