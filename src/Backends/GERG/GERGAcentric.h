// GENERATED FILE -- do not edit by hand.
//
// Acentric factors DERIVED FROM THE GERG PURE EQUATIONS OF STATE, produced by
// dev/gerg/compute_acentric.py from teqp 0.23.2
// (https://github.com/usnistgov/teqp).  CoolProp's build has no dependency on
// teqp; this header is committed, exactly like GERGAncillaries.h and
// GERGReferenceValues.h.
//
// GERG-2004/GERG-2008 publish no acentric factor.  These are NOT CoolProp's
// tabulated acentric factors for the same fluids -- those belong to each
// fluid's REFERENCE equation of state (Setzmann-Wagner for methane, IAPWS-95
// for water, ...), and a backend named GERG-2004/GERG-2008 must not carry
// data traceable to a different equation.  Each value here is Pitzer's
// defining relation evaluated on GERG's own shortened technical form:
//
//     omega = -1 - log10( p_sat(0.7 * T_c) / p_c )
//
//   * p_sat(0.7*T_c) is a CONVERGED saturation solve on the GERG pure EOS
//     (worst |p' - p''|/p'' across all 23: 5.55e-12), NOT the fitted pV
//     ancillary -- an ancillary is a seed, never an answer.
//   * p_c is the pressure the EOS produces at GERG's TABULATED reducing state
//     (Table A3.5), i.e. the same number make_gerg_fluid stores in
//     EOS.reduce.p and the backend reports from p_critical(); T_c likewise.
//     That consistency is load-bearing -- the Wilson K-factor seed reads
//     p_c/T_c and omega off the same fluid.
//
// See dev/gerg/compute_acentric.py for both judgement calls in full, and for
// why helium and hydrogen -- whose enforced Tmin equals their T_c, so 0.7*T_c
// is below the range a caller can reach -- are computed the same way as the
// other 21 rather than left out.
//
// RECOMPUTE WHENEVER A GERG COEFFICIENT TABLE CHANGES.
// dev/gerg/verify_transcription.py re-derives this whole table and diffs it
// against these literals, so a stale table is caught mechanically.
//
// 23 rows, not 21: GERG-2008 moves the reducing parameters of carbon monoxide
// and isopentane, which moves their saturation curve, so each model has its
// own value.
#ifndef GERGACENTRIC_H_
#define GERGACENTRIC_H_

#include <map>
#include <string>

#include "GERGData.h"

namespace CoolProp {
namespace GERG {
namespace detail {

/// Keyed on the GERG component name.  GERG-2004 is the base table and
/// GERG-2008 overrides it, mirroring pure_info_2004/pure_info_2008_overrides
/// in GERGData.h.
inline const std::map<std::string, double>& acentric_2004() {
    static const std::map<std::string, double> data = {
      {"methane", 0.011384874033444126},  {"nitrogen", 0.037320074493879085}, {"carbondioxide", 0.22495338333932602},
      {"ethane", 0.0994943063462157},     {"propane", 0.15292673801041246},   {"n-butane", 0.19923650256828696},
      {"isobutane", 0.18459379857066072}, {"n-pentane", 0.2515823575488525},  {"isopentane", 0.2270240996550561},
      {"n-hexane", 0.3001762417308751},   {"n-heptane", 0.34864119722310827}, {"n-octane", 0.3948641294704356},
      {"hydrogen", -0.21865183144536715}, {"oxygen", 0.02168668809639951},    {"carbonmonoxide", 0.050916570539008665},
      {"water", 0.34496944017437237},     {"helium", -0.3859393712713193},    {"argon", -0.0024086270235262885}};
    return data;
}

/// Entries GERG-2008 changes or adds relative to GERG-2004.
inline const std::map<std::string, double>& acentric_2008_overrides() {
    static const std::map<std::string, double> data = {{"carbonmonoxide", 0.05026121956623508},
                                                       {"isopentane", 0.22745573457718304},
                                                       {"hydrogensulfide", 0.1004219317938635},
                                                       {"n-nonane", 0.443445888629316},
                                                       {"n-decane", 0.48801513532779306}};
    return data;
}

}  // namespace detail
}  // namespace GERG
}  // namespace CoolProp
#endif /* GERGACENTRIC_H_ */
