#ifndef GERGDATA_H_
#define GERGDATA_H_

// Coefficient tables for the GERG-2004 and GERG-2008 equations of state.
//
// Transcribed from teqp (https://github.com/usnistgov/teqp),
// include/teqp/models/GERG/GERG.hpp, which is the reference implementation
// this backend is validated against.  Original data: Kunz, Klimeck, Wagner,
// Jaeschke, "The GERG-2004 Wide-Range Equation of State for Natural Gases and
// Other Mixtures", GERG TM15 (2007), and Kunz & Wagner, J. Chem. Eng. Data 57
// (2012) 3032-3091.
//
// This header is PRIVATE to the GERG backend.  It is not installed and must
// not be included from include/CoolProp/.

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"

namespace CoolProp {
namespace GERG {

/// Which GERG model year this backend instance represents.
enum class GERGModel : std::uint8_t
{
    GERG_2004,
    GERG_2008
};

struct PureInfo
{
    double rhoc_molm3;  ///< Reducing density, mol/m^3
    double Tc_K;        ///< Reducing temperature, K
    double M_kgmol;     ///< Molar mass, kg/mol
};

/// Pure-fluid residual Helmholtz coefficients: alpha^r = sum_i n_i
/// delta^{d_i} tau^{t_i} exp(-c_i delta^{l_i}).  c_i is 1.0 exactly when
/// l_i > 0 and 0.0 otherwise; CoolProp's ResidualHelmholtzGeneralizedExponential
/// (add_Power, include/CoolProp/fluids/Helmholtz.h) derives c from l the same
/// way, so c is kept here only for one-to-one comparison against teqp.
struct PureCoeffs
{
    std::vector<double> n, t, d, c, l;
};

/// The two gas constants that appear in the GERG ideal-gas part.  R* is the
/// value used when the ideal-gas c_p correlations were fitted; R is the value
/// the equation of state itself uses.  teqp GERG.hpp:470-471.
constexpr double R_GERG = 8.314472;      ///< J/mol/K
constexpr double RSTAR_GERG = 8.314510;  ///< J/mol/K

/// Ideal-gas Helmholtz coefficients, GERG-2004 monograph Table A3.1.
///
/// Both vectors have length 8 so that the indices match the monograph:
/// n0[1..7] and theta0[4..7] are used, the rest are zero padding.
///
///   alpha^o_i = ln(rho/rho_c,i)
///             + (R*/R) * [  n0[1] + n0[2]*Tc/T + n0[3]*ln(Tc/T)
///                         + n0[4]*ln(|sinh(theta0[4]*Tc/T)|)
///                         + n0[6]*ln(|sinh(theta0[6]*Tc/T)|)
///                         - n0[5]*ln(|cosh(theta0[5]*Tc/T)|)
///                         - n0[7]*ln(|cosh(theta0[7]*Tc/T)|) ]
///
/// SIGN CONVENTION (load-bearing for Task 8): every n0 is stored exactly as
/// published, UN-NEGATED.  Some published n0 are themselves negative (methane's
/// n0[7] is -4.46921); "stored as published" is the rule, not "stored
/// positive".  The minus in front of the two cosh terms belongs to the
/// *expression above*, not to the coefficient.  CoolProp's
/// IdealHelmholtzGERG2004Cosh::all (src/Helmholtz.cpp) accumulates
/// `-n[i]*log(|cosh(...)|)`, i.e. it applies that minus itself, so it must be
/// handed n0[5] and n0[7] unmodified.  Negating them on the way in would
/// silently flip those contributions to h and s while leaving p, c_v and w
/// untouched.  IdealHelmholtzGERG2004Sinh accumulates `+n[i]*log(|sinh(...)|)`,
/// so n0[4] and n0[6] likewise go in unmodified.
///
/// UNITS: n0 and theta0 are dimensionless; theta0 multiplies Tc/T.  The
/// (R*/R) prefactor scales the WHOLE bracket and is NOT folded into the
/// stored coefficients, so Task 8 must apply it itself (CoolProp's
/// GERG2004Cosh/GERG2004Sinh terms do not know about it).
///
/// MIXTURES: Tc here is the COMPONENT's own PureInfo::Tc_K, not the mixture
/// reducing temperature T_red -- and NO Tc/T_red folding is needed or wanted.
/// An earlier draft of this comment claimed the opposite; it was wrong, and
/// folding would double-apply a rescale CoolProp already performs.  BOTH
/// branches of HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache hand a
/// component's alpha0 container that component's own Tc,i/T: the pure branch
/// evaluates at taustar = Tc/Tr*tau, the mixture branch at
/// tau_i = T_ci*tau/Tr, and both equal Tc,i/T.  So the coefficients go in
/// exactly as tabulated.  See the long comment on make_gerg_fluid in
/// GERGBackend.cpp, which also explains why the sinh/cosh terms are
/// re-expressed as Lead + PlanckEinsteinGeneralized instead of being handed to
/// IdealHelmholtzGERG2004Sinh/Cosh (those DO consult T_red, and the mixture
/// branch sets it to Tr, which would leave a spurious factor Tc,i/Tr inside
/// every sinh and cosh -- invisible for a pure fluid, where Tr == Tc).
struct AlphaigCoeffs
{
    std::vector<double> n0, theta0;
};

/// Binary reducing-parameter coefficients for one pair of components.
/// Field order matches teqp's struct EXACTLY (betaV, gammaV, betaT, gammaT --
/// note V before T) so that brace-initialised rows copied from teqp GERG.hpp
/// land in the right fields; reordering these compiles fine and silently
/// corrupts every mixture.
struct BetasGammas
{
    double betaV, gammaV, betaT, gammaT;
};

namespace detail {

/// Table A3.5, GERG-2004 monograph.  Tabulated in mol/dm^3, K, kg/kmol;
/// converted to mol/m^3, K, kg/mol on the way out of get_pure_info.
/// Transcribed from teqp GERG.hpp:438-457 (GERG2004::get_pure_info data_map).
inline const std::map<std::string, PureInfo>& pure_info_2004() {
    static const std::map<std::string, PureInfo> data = {{"methane", {10.139342719, 190.564000000, 16.042460}},
                                                         {"nitrogen", {11.183900000, 126.192000000, 28.013400}},
                                                         {"carbondioxide", {10.624978698, 304.128200000, 44.009500}},
                                                         {"ethane", {6.870854540, 305.322000000, 30.069040}},
                                                         {"propane", {5.000043088, 369.825000000, 44.095620}},
                                                         {"n-butane", {3.920016792, 425.125000000, 58.122200}},
                                                         {"isobutane", {3.860142940, 407.817000000, 58.122200}},
                                                         {"n-pentane", {3.215577588, 469.700000000, 72.148780}},
                                                         {"isopentane", {3.271018581, 460.350000000, 72.148780}},
                                                         {"n-hexane", {2.705877875, 507.820000000, 86.175360}},
                                                         {"n-heptane", {2.315324434, 540.130000000, 100.201940}},
                                                         {"n-octane", {2.056404127, 569.320000000, 114.228520}},
                                                         {"hydrogen", {14.940000000, 33.190000000, 2.015880}},
                                                         {"oxygen", {13.630000000, 154.595000000, 31.998800}},
                                                         {"carbonmonoxide", {10.850000000, 132.800000000, 28.010100}},
                                                         {"water", {17.873716090, 647.096000000, 18.015280}},
                                                         {"helium", {17.399000000, 5.195300000, 4.002602}},
                                                         {"argon", {13.407429659, 150.687, 39.948000}}};
    return data;
}

/// Entries that GERG-2008 changes or adds relative to GERG-2004.  Everything
/// else falls through to pure_info_2004().
/// Transcribed from teqp GERG.hpp:980-985 (GERG2008::get_pure_info data_map).
inline const std::map<std::string, PureInfo>& pure_info_2008_overrides() {
    static const std::map<std::string, PureInfo> data = {
      {"carbonmonoxide", {10.85, 132.86, 28.010100}},  // changed from GERG-2004
      {"isopentane", {3.271, 460.35, 72.148780}},      // changed from GERG-2004
      {"n-nonane", {1.81, 594.55, 128.2551}},          // new in GERG-2008
      {"n-decane", {1.64, 617.7, 142.28168}},          // new in GERG-2008
      {"hydrogensulfide", {10.19, 373.1, 34.08088}}    // new in GERG-2008
    };
    return data;
}

/// CAS number -> GERG component name.  Used by resolve_component so that
/// CoolProp aliases (CO2, R744, 124-38-9, ...) reach the right component.
/// This table does not exist in teqp; CAS numbers verified against this
/// CoolProp build's fluid library (see task-2-report.md for the verification
/// transcript).
inline const std::map<std::string, std::string>& cas_to_gerg() {
    static const std::map<std::string, std::string> data = {
      {"74-82-8", "methane"},    {"7727-37-9", "nitrogen"}, {"124-38-9", "carbondioxide"},    {"74-84-0", "ethane"},
      {"74-98-6", "propane"},    {"106-97-8", "n-butane"},  {"75-28-5", "isobutane"},         {"109-66-0", "n-pentane"},
      {"78-78-4", "isopentane"}, {"110-54-3", "n-hexane"},  {"142-82-5", "n-heptane"},        {"111-65-9", "n-octane"},
      {"1333-74-0", "hydrogen"}, {"7782-44-7", "oxygen"},   {"630-08-0", "carbonmonoxide"},   {"7732-18-5", "water"},
      {"7440-59-7", "helium"},   {"7440-37-1", "argon"},    {"7783-06-4", "hydrogensulfide"}, {"111-84-2", "n-nonane"},
      {"124-18-5", "n-decane"}};
    return data;
}

}  // namespace detail

/// GERG2004::component_names, teqp GERG.hpp:431.
inline const std::vector<std::string>& component_names(GERGModel model) {
    static const std::vector<std::string> names_2004 = {"methane",   "nitrogen",  "carbondioxide",  "ethane",   "propane",   "n-butane",
                                                        "isobutane", "n-pentane", "isopentane",     "n-hexane", "n-heptane", "n-octane",
                                                        "hydrogen",  "oxygen",    "carbonmonoxide", "water",    "helium",    "argon"};
    // GERG2008::component_names, teqp GERG.hpp:973, is names_2004 plus these three.
    static const std::vector<std::string> names_2008 = [] {
        std::vector<std::string> v = names_2004;
        v.insert(v.end(), {"hydrogensulfide", "n-nonane", "n-decane"});
        return v;
    }();
    return (model == GERGModel::GERG_2004) ? names_2004 : names_2008;
}

inline PureInfo get_pure_info(GERGModel model, const std::string& gerg_name) {
    PureInfo data{};
    bool found = false;
    if (model == GERGModel::GERG_2008) {
        const auto& ov = detail::pure_info_2008_overrides();
        auto it = ov.find(gerg_name);
        if (it != ov.end()) {
            data = it->second;
            found = true;
        }
    }
    if (!found) {
        // GERG-2004 does not contain the three fluids added in GERG-2008.
        const auto& names = component_names(model);
        if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
            throw ValueError(format("[%s] is not a component of this GERG model", gerg_name.c_str()));
        }
        const auto& base = detail::pure_info_2004();
        auto it = base.find(gerg_name);
        if (it == base.end()) {
            throw ValueError(format("Unable to load GERG pure info for [%s]", gerg_name.c_str()));
        }
        data = it->second;
    }
    data.rhoc_molm3 *= 1000;  // mol/dm^3 -> mol/m^3
    data.M_kgmol /= 1000;     // kg/kmol -> kg/mol
    return data;
}

/// Resolve a CoolProp fluid name/alias/CAS to the corresponding GERG
/// component name, throwing ValueError if it cannot be resolved or if it
/// falls outside the given model's published component set.  Defined in
/// GERGBackend.cpp because it needs get_fluid_param_string(), which would
/// otherwise pull all of CoolProp.h into this data-only header.
std::string resolve_component(GERGModel model, const std::string& user_name);

/// Pure-fluid residual coefficients (GERG-2004 monograph Table A3.2 and its
/// GERG-2008 overrides/additions).  Throws ValueError if gerg_name is not a
/// component of the given model.  Defined in GERGBackend.cpp: the tables are
/// the largest data block in this backend, and keeping them out of line
/// keeps this header's compile time down.
PureCoeffs get_pure_coeffs(GERGModel model, const std::string& gerg_name);

/// Solve the 2x2 linear system that makes the IDEAL-GAS enthalpy and entropy
/// both vanish at the reference state (T0, rho0), returning {n0[1], n0[2]}.
///
/// Direct port of teqp's AlphaigCoeffs::recalc_integration_constants
/// (GERG.hpp:49-74), with the 2x2 Eigen solve written out by hand.  Only
/// c.n0[3..7] and c.theta0[4..7] are read; c.n0[1] and c.n0[2] are ignored.
///
/// @param c      Coefficients as published (n0[1], n0[2] are not used)
/// @param T0     Reference temperature, K
/// @param Tci    Reducing temperature of the component, K
/// @param rho0   Molar density at the reference state, mol/m^3
/// @param rhoci  Reducing molar density of the component, mol/m^3
/// @param Rstar_R  the ratio R*/R
std::pair<double, double> recalc_integration_constants(const AlphaigCoeffs& c, double T0, double Tci, double rho0, double rhoci, double Rstar_R);

/// Ideal-gas coefficients for one component, with n0[1] and n0[2] REPLACED by
/// the values that give h = s = 0 for the ideal gas at 298.15 K and 101325 Pa.
/// The published integration constants are deliberately discarded, exactly as
/// teqp does in GERG200XAlphaig::get_coeffs (GERG.hpp:370-382); they refer to a
/// different reference state and using them would leave p, c_v and w correct
/// while h and s were quietly wrong.  Throws ValueError if gerg_name is not a
/// component of the given model.
AlphaigCoeffs get_alphaig_coeffs(GERGModel model, const std::string& gerg_name);

/// Binary reducing-parameter beta/gamma values (GERG-2004 monograph Table
/// A3.8; GERG-2008 manuscript Table A8 overrides 15 of those pairs and adds
/// 57 new ones for the three fluids introduced in GERG-2008).  Defined in
/// GERGBackend.cpp: this is the largest data block in the backend, and
/// keeping it out of line keeps this header's compile time down.
///
/// ORIENTATION (load-bearing for Task 9): the underlying table stores each
/// pair once, in one order.  This accessor returns values oriented for the
/// (f1, f2) order GIVEN, not the order stored.  If the stored row is for
/// (f2, f1), betaV and betaT are returned RECIPROCATED (1/betaV, 1/betaT);
/// gammaV and gammaT are symmetric between the two orderings and are
/// returned unchanged.  This mirrors teqp GERG.hpp:757-763 exactly.  Throws
/// ValueError if the pair is not found in either order.
BetasGammas get_betasgammas(GERGModel model, const std::string& f1, const std::string& f2);

/// Departure-function coefficients for one pair of components (GERG-2004
/// monograph Table A3.6 and the "generalized" departure function it defines).
/// Field names match CoolProp's GERG2008DepartureFunction constructor
/// (src/Backends/Helmholtz/ExcessHEFunction.h:108), which takes
/// (n, d, t, eta, epsilon, beta, gamma, Npower) -- Task 9 passes these fields
/// by name in that order, so the struct's OWN field order below does not
/// matter, only the names do.
struct DepartureCoeffs
{
    std::vector<double> n, t, d, eta, beta, gamma, epsilon;
};

/// GERG-2004 monograph Table A3.6.  Fewer than half of all component pairs
/// have a departure function at all; most pairs use the corresponding-states
/// term alone.  Returns false (and leaves F untouched) when the pair has no
/// departure function, mirroring teqp's `ok_missing = true` path
/// (GERG.hpp:768-802) without requiring <optional> in this header.
///
/// GERG-2008 reuses GERG-2004's F_ij unchanged (teqp GERG.hpp:1000,
/// `using GERG2004::get_Fij;`), so there is one shared table for both models
/// and the model argument is accepted only for interface symmetry with the
/// other accessors in this header.
///
/// F_ij is SYMMETRIC: F_ij == F_ji exactly, unlike the betas/gammas above --
/// do not reciprocate anything on a reversed-order lookup.
bool get_Fij(GERGModel model, const std::string& f1, const std::string& f2, double& F);

/// GERG-2004 monograph Table A3.6 departure-function coefficients: the 7
/// pairs with fluid-pair-specific coefficients, plus the "generalized"
/// departure function (shared, scaled per-pair by get_Fij) used by 8 more
/// pairs.  15 pairs total carry a departure function; every other pair
/// throws ValueError.
///
/// GERG-2008 reuses GERG-2004's departure coefficients unchanged (teqp
/// GERG.hpp:1001, `using GERG2004::get_departurecoeffs;`), so there is one
/// shared table for both models; the model argument exists only for
/// interface symmetry.
DepartureCoeffs get_departurecoeffs(GERGModel model, const std::string& f1, const std::string& f2);

/// One saturation ancillary, in exactly the schema dev/fluids/*.json's
/// ANCILLARIES block uses, so that make_gerg_fluid can hand it straight to
/// CoolProp's own SaturationAncillaryFunction and the shipped evaluation path
/// (Ancillaries.cpp) is reused rather than reimplemented.  `type` names the
/// functional form with the same strings cpjson::make_saturation_ancillary
/// accepts (FluidLibraryFactories.h:41-73):
///
///   "rhoLnoexp"  (TYPE_NOT_EXPONENTIAL, using_tau_r = false)
///        y = reducing_value * (1 + sum_i n_i theta^{t_i})
///   "rhoV", "pV" (TYPE_EXPONENTIAL, using_tau_r = true)
///        y = reducing_value * exp(T_r/T * sum_i n_i theta^{t_i})
///
///   theta = 1 - T/T_r
///
/// UNITS: reducing_value carries the units of the output -- mol/m^3 for
/// "rhoLnoexp"/"rhoV", Pa for "pV".  n and t are dimensionless; T_r, Tmin and
/// Tmax are in K.
///
/// AN ANCILLARY IS A SEED, NOT AN ANSWER.  QT_flash hands these to
/// SaturationSolvers::saturation_T_pure_Maxwell, which iterates to the true
/// GERG saturation state; no ancillary value is ever returned to a caller.
/// That is also why no superancillary may be attached to a GERG fluid:
/// FlashRoutines::sat_superanc_path_applies (FlashRoutines.cpp:558) routes
/// pure-fluid saturation straight to the Chebyshev expansion and RETURNS THAT
/// as the answer, which for a GERG fluid would be a reference-EOS saturation
/// density labelled GERG-2008.
///
/// max_rel_dev is the worst relative deviation of this fit from the traced
/// teqp VLE curve over [Tmin, T_r], as measured by dev/gerg/fit_ancillaries.py
/// and carried along so a test can assert against the fit's own claim rather
/// than a hard-coded number.
///
/// A t_i == 0 term appears in every fit deliberately: GERG's Table A3.5
/// reducing parameters are fitted quantities rather than the true critical
/// point of the shortened form (they differ by up to +1.10 K and -5.4% in
/// density), so the usual "the ancillary is exact at the critical point"
/// property does not hold and the value at T = T_r has to stay free.  See the
/// module docstring of dev/gerg/fit_ancillaries.py.
struct AncillaryCoeffs
{
    std::vector<double> n, t;
    double reducing_value, T_r, Tmin, Tmax, max_rel_dev;
    std::string type;
};

/// The saturation state at the low-temperature end of the fitted ancillary
/// range, traced with teqp rather than guessed.  make_gerg_fluid writes it
/// into EOS.sat_min_liquid / EOS.sat_min_vapor (which is what
/// HelmholtzEOSMixtureBackend::calc_Tmin_sat and calc_pmin_sat return, and
/// hence what bounds QT_flash from below) and into
/// CoolPropFluid::triple_liquid / triple_vapor.
///
/// THIS IS NOT A TRIPLE POINT.  GERG publishes none, and EOS.Ttriple stays 0
/// (see make_gerg_fluid).  CoolProp's `triple_liquid`/`triple_vapor` slots are
/// read by saturation_T_pure_Maxwell (VLERoutines.cpp:965-980) purely as the
/// low-temperature end of the saturation curve, to sanity-band the ancillary
/// seed and to build a linear fallback seed; the lowest temperature at which
/// this backend has saturation data is exactly the right value for that use.
struct SatEndState
{
    double T_K, p_Pa, rhoL_molm3, rhoV_molm3;
};

/// Saturation ancillary for one component, fitted against the GERG pure EOS
/// (dev/gerg/fit_ancillaries.py, table in GERGAncillaries.h).  `which` is
/// "rhoL", "rhoV" or "pV"; anything else throws ValueError, as does a
/// gerg_name that is not a component of the given model.  The returned
/// `type` is filled in from `which` by this accessor -- the generated table
/// leaves it empty -- so the table cannot disagree with the accessor about
/// which functional form a row is.
AncillaryCoeffs get_ancillary(GERGModel model, const std::string& gerg_name, const std::string& which);

/// Evaluate an ancillary at temperature T, returning the same units as
/// AncillaryCoeffs::reducing_value (Pa for "pV", mol/m^3 for "rhoL"/"rhoV").
/// Mirrors SaturationAncillaryFunction::evaluate exactly, including its
/// NaN-for-T-above-T_r behaviour; it exists so a test can evaluate a fit
/// without building a whole fluid.
double evaluate_ancillary(const AncillaryCoeffs& anc, double T);

/// Saturation state at the low-temperature end of the fitted range.  Throws
/// ValueError if gerg_name is not a component of the given model.
SatEndState get_sat_min_state(GERGModel model, const std::string& gerg_name);

/// Acentric factor of one GERG pure EOS, derived from that EOS by Pitzer's
/// defining relation omega = -1 - log10(p_sat(0.7 T_c) / p_c) with a CONVERGED
/// saturation solve and the pressure at GERG's tabulated reducing state
/// (dev/gerg/compute_acentric.py, table in GERGAcentric.h).  GERG publishes no
/// acentric factor, and CoolProp's tabulated one belongs to a different
/// equation of state, so it is computed rather than borrowed.
///
/// GERG-2008 has its own value for carbon monoxide and isopentane (it moved
/// their reducing parameters), hence 23 rows for 21 components.  Throws
/// ValueError if gerg_name is not a component of the given model, or if the
/// generated table holds a non-finite value.
double get_acentric_factor(GERGModel model, const std::string& gerg_name);

/// Number of leading terms in dc.n (etc.) with eta == 0 and beta == 0 -- the
/// polynomial block that CoolProp's GERG2008DepartureFunction constructor
/// (ExcessHEFunction.h:108) expects as a contiguous prefix, with every
/// remaining term forming the trailing Gaussian block.  Throws ValueError if
/// a term AFTER that leading run is also polynomial (eta == 0 && beta == 0):
/// CoolProp's constructor cannot represent a non-contiguous split, and doing
/// so silently would push Gaussian terms into the polynomial block with no
/// error anywhere.
std::size_t departure_Npower(const DepartureCoeffs& dc);

}  // namespace GERG
}  // namespace CoolProp
#endif /* GERGDATA_H_ */
