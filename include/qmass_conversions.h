#ifndef QMASS_CONVERSIONS_H
#define QMASS_CONVERSIONS_H

namespace CoolProp {
namespace detail {

/// Mass-basis vapor quality from molar quality and phase molar masses.
/// All molar masses in the same units (kg/mol).
inline double Qmolar_to_Qmass(double Qmolar, double MM_liquid, double MM_vapor) {
    const double mv = Qmolar * MM_vapor;
    const double ml = (1.0 - Qmolar) * MM_liquid;
    return mv / (mv + ml);
}

/// Inverse of Qmolar_to_Qmass. Same units convention.
inline double Qmass_to_Qmolar(double Qmass, double MM_liquid, double MM_vapor) {
    const double nv = Qmass / MM_vapor;
    const double nl = (1.0 - Qmass) / MM_liquid;
    return nv / (nv + nl);
}

}  // namespace detail
}  // namespace CoolProp

#endif  // QMASS_CONVERSIONS_H
