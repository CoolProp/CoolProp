#ifndef COOLPROP_REGION_CONSTANT_CURVE_H
#define COOLPROP_REGION_CONSTANT_CURVE_H

#include "CoolProp/region/BoundaryCurve.h"

namespace CoolProp {
namespace region {

// Trivial BoundaryCurve b(a) = c, used for fixed ceilings/floors such as
// a T_max or p_max plane.  All three queries are O(1).
class ConstantCurve final : public BoundaryCurve
{
   public:
    ConstantCurve(double a_lo, double a_hi, double b_value) noexcept : a_lo_(a_lo), a_hi_(a_hi), b_(b_value) {}

    [[nodiscard]] double eval(double /*a*/) const noexcept override {
        return b_;
    }

    // db/da is exactly zero — analytic derivative of the constant function.
    [[nodiscard]] double eval_da(double /*a*/) const noexcept override {
        return 0.0;
    }

    [[nodiscard]] std::pair<double, double> bounds() const noexcept override {
        return {b_, b_};
    }

    [[nodiscard]] std::pair<double, double> a_range() const noexcept override {
        return {a_lo_, a_hi_};
    }

   private:
    double a_lo_;
    double a_hi_;
    double b_;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_CONSTANT_CURVE_H
