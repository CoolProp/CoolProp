#include "GERGBackend.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <memory>

#include "CoolProp/CoolProp.h"
#include "CoolProp/Configuration.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/fluids/Ancillaries.h"
#include "GERGAcentric.h"
#include "GERGAncillaries.h"
#include "GERGData.h"
#include "Backends/Helmholtz/ExcessHEFunction.h"
#include "Backends/Helmholtz/ReducingFunctions.h"

namespace CoolProp {

GERGMixtureBackend::GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names) : m_model(model) {
    if (names.empty()) {
        throw ValueError("GERG backend requires at least one component");
    }
    std::vector<CoolPropFluid> fluids;
    fluids.reserve(names.size());
    for (const auto& user_name : names) {
        // resolve_component throws if the fluid is outside this model's
        // published component set, so an unsupported component can never
        // reach make_gerg_fluid.
        fluids.push_back(GERG::make_gerg_fluid(model, GERG::resolve_component(model, user_name)));
    }
    // m_model is initialised in the member-init list above, so it is already
    // valid here.  The call is EXPLICITLY QUALIFIED because that is what it
    // means: inside a GERGMixtureBackend constructor body the dynamic type is
    // GERGMixtureBackend, so an unqualified virtual call resolves to this
    // class's override and NOT to any further-derived override -- which is
    // exactly the behaviour wanted, but reads as an accident and is what
    // cppcheck's virtualCallInConstructor flags.  Writing the qualification
    // out makes the static binding intentional rather than incidental, and
    // keeps it correct if this class is ever subclassed.
    GERGMixtureBackend::set_components(fluids);
    if (names.size() == 1) {
        set_mole_fractions(std::vector<CoolPropDbl>(1, 1.0));
    }
}

GERGMixtureBackend::GERGMixtureBackend(GERGModel model, const std::vector<CoolPropFluid>& fluids, bool generate_SatL_and_SatV) : m_model(model) {
    // Qualified for the same reason as the constructor above.
    GERGMixtureBackend::set_components(fluids, generate_SatL_and_SatV);
}

void GERGMixtureBackend::set_components(const std::vector<CoolPropFluid>& comps, bool generate_SatL_and_SatV) {
    // Build everything EXCEPT the linked saturation states.  The base
    // implementation calls set_mixture_parameters() unqualified for N > 1, so
    // that part reaches this class's override.
    HelmholtzEOSMixtureBackend::set_components(comps, false);

    if (!generate_SatL_and_SatV) {
        return;
    }
    // Mirror HelmholtzEOSMixtureBackend::set_components' SatL/SatV block, but
    // with GERG-typed states.  Constructing them with
    // generate_SatL_and_SatV = false is what stops the recursion.
    //
    // The PREVIOUS SatL/SatV (if any -- set_components is only called from
    // constructors today, so on a freshly-constructed object there is
    // nothing to remove) are erased from linked_states before being replaced,
    // so calling set_components twice on one object does not leave the
    // superseded SatL/SatV in linked_states, where sync_linked_states would
    // keep writing to them forever.  This erases ONLY the two SatL/SatV
    // entries -- by pointer identity, not by clearing the whole vector -- so
    // TPD_state/critical_state/transient_pure_state (HelmholtzEOSMixtureBackend.h:85,93,101)
    // are left untouched if they happen to already be linked.  An earlier
    // version called linked_states.clear() unconditionally: harmless today
    // because those three are only ever created lazily, long after
    // construction, but a future caller of set_components on a live object
    // would silently strand them (non-null members, never reachable through
    // linked_states again) -- trading one latent hazard for another instead
    // of fixing it.
    linked_states.erase(std::remove_if(linked_states.begin(), linked_states.end(),
                                       [this](const shared_ptr<HelmholtzEOSMixtureBackend>& s) { return s == SatL || s == SatV; }),
                        linked_states.end());

    SatL.reset(new GERGMixtureBackend(m_model, comps, false));
    SatL->specify_phase(iphase_liquid);
    linked_states.push_back(SatL);
    SatL->clear();

    SatV.reset(new GERGMixtureBackend(m_model, comps, false));
    SatV->specify_phase(iphase_gas);
    SatV->clear();
    linked_states.push_back(SatV);
}

void GERGMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    HelmholtzEOSMixtureBackend::update(input_pair, value1, value2);
    check_gerg_range_of_validity();
}

void GERGMixtureBackend::check_gerg_range_of_validity() {
    // T only, deliberately NOT p: EOS.limits.pmax (70 MPa, make_gerg_fluid)
    // is the mixture MODEL's published operating envelope, but "GERG pure
    // fluids/mixtures reproduce teqp" (this file) legitimately evaluates
    // fixed (T, rho) reference points whose pressure is far outside it --
    // including deep into the two-phase dome, where teqp's single-phase EOS
    // (and this backend, matching it) is expected to return a large or
    // negative p.  A pmax check here does not distinguish "unphysical
    // operating point" from "single-phase EOS evaluated where a real fluid
    // would be two-phase", and adding one made exactly those reference tests
    // fail with p up to 3.1e10 Pa at a perfectly ordinary (250 K, 5000
    // mol/m^3) grid point.  Matches every other property-limit guard in
    // HelmholtzEOSMixtureBackend.cpp (melting-line, Tmax_sat, ...): a caller
    // who has explicitly opted out gets no check, rather than this backend
    // enforcing a range CoolProp's own configuration says to ignore.
    if (CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS)) {
        return;
    }
    // ####################################################################
    // COMPONENT-WISE bounds, NOT Tmin()/Tmax().
    //
    // AbstractState::Tmin()/Tmax() resolve, for a mixture, to
    // HelmholtzEOSMixtureBackend::calc_Tmin/calc_Tmax
    // (HelmholtzEOSMixtureBackend.cpp:1305-1319), which are
    // `sum_i x_i * limits.T{min,max}` -- a mole-fraction weighting that is
    // NOT normalised by sum(x).  Nothing in set_mole_fractions or PropsSI's
    // `Fluid[x]&Fluid[x]` syntax requires the composition to sum to 1 (the
    // string parser rejects an individual fraction > 1, not the sum), so a
    // range guard built on them SCALES ITS OWN BOUNDS WITH sum(x) and fails
    // open: measured on this backend before this change, a 50/50 methane +
    // ethane mixture entered as [0.6]&[0.6] (sum = 1.2) reported
    // [72, 840] K, and sum(x) = 10 reported [600, 7000] K -- 900 K accepted,
    // rho = 309.922 mol/m^3 returned, no error anywhere.  The benign-input
    // version of the same coupling is sharper still: a gas analysis rounded
    // to sum(x) = 0.999999 moves Tmax to 699.9993 K and throws a spurious
    // OutOfRangeError at exactly 700 K.
    //
    // The bounds below are read straight off the per-component EOS limits
    // that make_gerg_fluid set, so they are independent of composition
    // entirely.  Choice of aggregation, made deliberately: the INTERSECTION
    // of the components' ranges -- max of the Tmin values, min of the Tmax
    // values.  For a strict backend the defensible reading of "the mixture
    // is valid where the model is valid" is "where EVERY constituent is
    // valid"; the union (or a weighted average) would let a mixture be
    // evaluated at a temperature at which one of its own components has been
    // declared out of range.  Concretely, with every GERG component carrying
    // Tmax = 700 K and Tmin = min(60 K, Tc), this makes the mixture range
    // exactly [60, 700] K unless EVERY component has Tc < 60 K (helium
    // alone: [5.1953, 700]; helium + hydrogen: [33.19, 700]), which is the
    // published mixture-model range of Kunz & Wagner 2012 section 4.1 rather
    // than an artefact of the composition vector.  NOTE this is a
    // deliberate behaviour change from the mole-fraction-weighted average:
    // a helium-rich mixture that used to be accepted at 45 K is now
    // rejected, because methane's own 60 K limit is part of the mixture.
    //
    // Deliberately NOT fixed by changing calc_Tmin/calc_Tmax in the shared
    // Helmholtz backend: those are used by every other backend and their
    // normalisation is a separate (ABI-affecting) question.  Tmin()/Tmax()
    // therefore still report the weighted values for a GERG mixture; only
    // this guard is composition-independent.
    // ####################################################################
    double Tlo = -std::numeric_limits<double>::infinity();
    double Thi = std::numeric_limits<double>::infinity();
    for (const CoolPropFluid& c : get_components()) {
        Tlo = std::max(Tlo, static_cast<double>(c.EOS().limits.Tmin));
        Thi = std::min(Thi, static_cast<double>(c.EOS().limits.Tmax));
    }
    // !ValidNumber first: a NaN _T makes BOTH comparisons false, so the bare
    // range test would wave it through.  Unreachable today (the base class
    // rejects a non-finite T before this runs), but a range check that passes
    // on NaN is a fail-open waiting for the day that stops being true.
    if (!ValidNumber(_T) || _T < Tlo || _T > Thi) {
        throw CoolProp::OutOfRangeError(format("Temperature [%g K] is outside the GERG range of validity [%g, %g] K", _T, Tlo, Thi));
    }
}

CoolPropDbl GERGMixtureBackend::calc_gas_constant() {
    // See the declaration in GERGBackend.h: the inherited implementation
    // returns the CODATA universal gas constant for any MIXTURE while
    // NORMALIZE_GAS_CONSTANTS is set (CoolProp's default), which is not the R
    // the GERG equation of state is written against.
    return GERG::R_GERG;
}

HelmholtzEOSMixtureBackend* GERGMixtureBackend::get_copy(bool generate_SatL_and_SatV) {
    auto* ptr = new GERGMixtureBackend(m_model, components, generate_SatL_and_SatV);
    ptr->sync_linked_states(this);
    return ptr;
}

void GERGMixtureBackend::set_mixture_parameters() {
    // ####################################################################
    // The reason this backend exists.
    //
    // HelmholtzEOSMixtureBackend::set_mixture_parameters delegates to
    // MixtureParameters::set_mixture_parameters, which resolves every binary
    // pair through CoolProp's GLOBAL binary-interaction-parameter library by
    // CAS number.  Measured against dev/mixtures/mixture_binary_pairs.json in
    // this tree (888 rows in total): all 210 GERG-2008 binary pairs ARE in
    // that library, 194 of them as the Kunz-JCED-2012 row GERG publishes --
    // but 16 as a LATER refit, 15 Gernert-Thesis-2013 and 1
    // Tkaczuk-JPCRD-2020.  Twelve of those 16 shift the mixture reducing
    // temperature by more than 0.03 K at z = (0.35, 0.65); four leave
    // beta/gamma numerically unchanged but attach a DIFFERENT departure
    // function.  Reaching this path from a backend named GERG2004/GERG2008
    // would return entirely plausible numbers that are not GERG, with no
    // error anywhere.  (Mutation-verified: forcing the inherited path, with a
    // CAS written into each fluid so the lookup succeeds, fails 294 assertions
    // -- 78 of them on alphar alone.)
    //
    // Everything below therefore comes from the GERG tables in GERGData.h /
    // this file, keyed on the GERG component name that make_gerg_fluid stored
    // in CoolPropFluid::name.  Note that make_gerg_fluid deliberately leaves
    // CoolPropFluid::CAS empty, so if this override were ever bypassed the
    // inherited path would throw rather than answer -- but that is a backstop,
    // not the mechanism: see set_components above, which is what keeps the
    // linked SatL/SatV states from taking the inherited path.
    // ####################################################################
    const std::vector<CoolPropFluid> comps = get_components();
    const std::size_t N = comps.size();

    STLMatrix beta_v(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix gamma_v(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix beta_T(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix gamma_T(N, std::vector<CoolPropDbl>(N, 0));

    residual_helmholtz->Excess.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            const std::string& f1 = comps[i].name;
            const std::string& f2 = comps[j].name;

            // get_betasgammas returns values already oriented for the (f1, f2)
            // order asked for, so beta_*[i][j] takes them as-is.  The
            // TRANSPOSED entry needs the reciprocal of the betas and the SAME
            // gammas -- gamma_v/gamma_T are symmetric under exchange, the
            // betas are not (GERGData.h, get_betasgammas' orientation note).
            const GERG::BetasGammas bg = GERG::get_betasgammas(m_model, f1, f2);
            beta_v[i][j] = bg.betaV;
            beta_v[j][i] = 1.0 / bg.betaV;
            gamma_v[i][j] = bg.gammaV;
            gamma_v[j][i] = bg.gammaV;
            beta_T[i][j] = bg.betaT;
            beta_T[j][i] = 1.0 / bg.betaT;
            gamma_T[i][j] = bg.gammaT;
            gamma_T[j][i] = bg.gammaT;

            double F = 0.0;
            const bool has_departure = GERG::get_Fij(m_model, f1, f2, F);
            residual_helmholtz->Excess.F[i][j] = has_departure ? F : 0.0;
            residual_helmholtz->Excess.F[j][i] = residual_helmholtz->Excess.F[i][j];

            DepartureFunctionPointer dep;
            if (has_departure) {
                const GERG::DepartureCoeffs dc = GERG::get_departurecoeffs(m_model, f1, f2);
                // GERG2008DepartureFunction's parameter order is
                // (n, d, t, eta, epsilon, beta, gamma, Npower) -- epsilon
                // BEFORE beta/gamma, which is NOT the order DepartureCoeffs
                // declares its members.  Every argument is therefore named
                // explicitly here; a transposed pair compiles cleanly and is
                // silently wrong.
                dep =
                  std::make_shared<GERG2008DepartureFunction>(dc.n, dc.d, dc.t, dc.eta, dc.epsilon, dc.beta, dc.gamma, GERG::departure_Npower(dc));
            } else {
                // Most GERG pairs have no departure function at all.  A null
                // pointer here would be a latent segfault rather than a zero:
                // ExcessTerm::update() dereferences EVERY off-diagonal entry
                // (ExcessHEFunction.h:239-248) and ExcessTerm::copy() does too
                // (:215-227), regardless of F being zero.  Install the same
                // zero-valued placeholder MixtureParameters uses
                // (MixtureParameters.cpp:629-634).
                const std::vector<double> n(1, 0), d(1, 1), t(1, 1), l(1, 0);
                dep = std::make_shared<ExponentialDepartureFunction>(n, d, t, l);
            }
            // alphar_ij == alphar_ji, and ExcessTerm::update calls
            // update(tau, delta) on each entry with identical arguments, so
            // sharing one instance across both slots is safe.
            residual_helmholtz->Excess.DepartureFunctionMatrix[i][j] = dep;
            residual_helmholtz->Excess.DepartureFunctionMatrix[j][i] = dep;
        }
    }

    // GERG2008ReducingFunction reads pFluids[i].EOS().reduce.T / .rhomolar,
    // which make_gerg_fluid sets from the GERG reducing state (Table A3.5).
    Reducing = std::make_shared<GERG2008ReducingFunction>(comps, beta_v, gamma_v, beta_T, gamma_T);
}

namespace GERG {

namespace {

// Pure-fluid residual coefficient tables.
//
// Transcribed from teqp (https://github.com/usnistgov/teqp),
// include/teqp/models/GERG/GERG.hpp: GERG2004::get_pure_coeffs (lines
// 511-591) and GERG2008::get_pure_coeffs (lines 1105-1131).  Preserving
// teqp's own split -- one shared 12-term t/d/c/l set for most fluids, one
// shared 24-term set for methane/nitrogen/ethane, and fully independent
// tables for carbon dioxide, hydrogen, water, and helium -- keeps this file
// diffable by eye against the source.

/// Shared exponent set: t, d, c, l (n differs per fluid).
struct SharedExponents
{
    std::vector<double> t, d, c, l;
};

/// 12-term set shared (GERG-2004) by propane, n-butane, isobutane,
/// n-pentane, isopentane, n-hexane, n-heptane, n-octane, oxygen, carbon
/// monoxide, and argon; also used by every fluid GERG-2008 overrides or adds
/// within this same family (carbon monoxide, isopentane, hydrogen sulfide,
/// n-nonane, n-decane).  teqp GERG.hpp:537-540 and :1122-1125.
const SharedExponents& main12_exponents() {
    static const SharedExponents e = {{0.250, 1.125, 1.500, 1.375, 0.250, 0.875, 0.625, 1.750, 3.625, 3.625, 14.500, 12.000},
                                      {1, 1, 1, 2, 3, 7, 2, 5, 1, 4, 3, 4},
                                      {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1},
                                      {0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3}};
    return e;
}

/// 24-term set shared by methane, nitrogen, and ethane. teqp GERG.hpp:546-549.
const SharedExponents& mne24_exponents() {
    static const SharedExponents e = {{0.125, 1.125, 0.375, 1.125, 0.625, 1.500, 0.625,  2.625,  2.750,  2.125,  2.000,  1.750,
                                       4.500, 4.750, 5.000, 4.000, 4.500, 7.500, 14.000, 11.500, 26.000, 28.000, 30.000, 16.000},
                                      {1, 1, 2, 2, 4, 4, 1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7},
                                      {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                                      {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6}};
    return e;
}

/// n-only table for the mne24 family, GERG-2004. teqp GERG.hpp:515-517.
const std::map<std::string, std::vector<double>>& n_mne_2004() {
    static const std::map<std::string, std::vector<double>> data = {
      {"methane",
       {0.57335704239162,     -0.16760687523730e1, 0.23405291834916,     -0.21947376343441,    0.16369201404128e-1,  0.15004406389280e-1,
        0.98990489492918e-1,  0.58382770929055,    -0.74786867560390,    0.30033302857974,     0.20985543806568,     -0.18590151133061e-1,
        -0.15782558339049,    0.12716735220791,    -0.32019743894346e-1, -0.68049729364536e-1, 0.24291412853736e-1,  0.51440451639444e-2,
        -0.19084949733532e-1, 0.55229677241291e-2, -0.44197392976085e-2, 0.40061416708429e-1,  -0.33752085907575e-1, -0.25127658213357e-2}},
      {"nitrogen",
       {0.59889711801201,     -0.16941557480731e1, 0.24579736191718,     -0.23722456755175,    0.17954918715141e-1,  0.14592875720215e-1,
        0.10008065936206,     0.73157115385532,    -0.88372272336366,    0.31887660246708,     0.20766491728799,     -0.19379315454158e-1,
        -0.16936641554983,    0.13546846041701,    -0.33066712095307e-1, -0.60690817018557e-1, 0.12797548292871e-1,  0.58743664107299e-2,
        -0.18451951971969e-1, 0.47226622042472e-2, -0.52024079680599e-2, 0.43563505956635e-1,  -0.36251690750939e-1, -0.28974026866543e-2}},
      {"ethane", {0.63596780450714,     -0.17377981785459e1, 0.28914060926272,     -0.33714276845694,   0.22405964699561e-1,  0.15715424886913e-1,
                  0.11450634253745,     0.10612049379745e1,  -0.12855224439423e1,  0.39414630777652,    0.31390924682041,     -0.21592277117247e-1,
                  -0.21723666564905,    -0.28999574439489,   0.42321173025732,     0.46434100259260e-1, -0.13138398329741,    0.11492850364368e-1,
                  -0.33387688429909e-1, 0.15183171583644e-1, -0.47610805647657e-2, 0.46917166277885e-1, -0.39401755804649e-1, -0.32569956247611e-2}}};
    return data;
}

/// n-only table for the main12 family, GERG-2004. teqp GERG.hpp:521-531.
const std::map<std::string, std::vector<double>>& n_main_2004() {
    static const std::map<std::string, std::vector<double>> data = {
      {"propane",
       {0.10403973107358e1, -0.28318404081403e1, 0.84393809606294, -0.76559591850023e-1, 0.94697373057280e-1, 0.24796475497006e-3, 0.27743760422870,
        -0.43846000648377e-1, -0.26991064784350, -0.69313413089860e-1, -0.29632145981653e-1, 0.14040126751380e-1}},
      {"n-butane",
       {0.10626277411455e1, -0.28620951828350e1, 0.88738233403777, -0.12570581155345, 0.10286308708106, 0.25358040602654e-3, 0.32325200233982,
        -0.37950761057432e-1, -0.32534802014452, -0.79050969051011e-1, -0.20636720547775e-1, 0.57053809334750e-2}},
      {"isobutane",
       {0.10429331589100e1, -0.28184272548892e1, 0.86176232397850, -0.10613619452487, 0.98615749302134e-1, 0.23948208682322e-3, 0.30330004856950,
        -0.41598156135099e-1, -0.29991937470058, -0.80369342764109e-1, -0.29761373251151e-1, 0.13059630303140e-1}},
      {"n-pentane",
       {0.10968643098001e1, -0.29988888298061e1, 0.99516886799212, -0.16170708558539, 0.11334460072775, 0.26760595150748e-3, 0.40979881986931,
        -0.40876423083075e-1, -0.38169482469447, -0.10931956843993, -0.32073223327990e-1, 0.16877016216975e-1}},
      {"isopentane",
       {0.11017531966644e1, -0.30082368531980e1, 0.99411904271336, -0.14008636562629, 0.11193995351286, 0.29548042541230e-3, 0.36370108598133,
        -0.48236083488293e-1, -0.35100280270615, -0.10185043812047, -0.35242601785454e-1, 0.19756797599888e-1}},
      {"n-hexane",
       {0.10553238013661e1, -0.26120615890629e1, 0.76613882967260, -0.29770320622459, 0.11879907733358, 0.27922861062617e-3, 0.46347589844105,
        0.11433196980297e-1, -0.48256968738131, -0.93750558924659e-1, -0.67273247155994e-2, -0.51141583585428e-2}},
      {"n-heptane",
       {0.10543747645262e1, -0.26500681506144e1, 0.81730047827543, -0.30451391253428, 0.12253868710800, 0.27266472743928e-3, 0.49865825681670,
        -0.71432815084176e-3, -0.54236895525450, -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3}},
      {"n-octane",
       {0.10722544875633e1, -0.24632951172003e1, 0.65386674054928, -0.36324974085628, 0.12713269626764, 0.30713572777930e-3, 0.52656856987540,
        0.19362862857653e-1, -0.58939426849155, -0.14069963991934, -0.78966330500036e-2, 0.33036597968109e-2}},
      {"oxygen",
       {0.88878286369701, -0.24879433312148e1, 0.59750190775886, 0.96501817061881e-2, 0.71970428712770e-1, 0.22337443000195e-3, 0.18558686391474,
        -0.38129368035760e-1, -0.15352245383006, -0.26726814910919e-1, -0.25675298677127e-1, 0.95714302123668e-2}},
      {"carbonmonoxide",
       {0.92310041400851, -0.24885845205800e1, 0.58095213783396, 0.28859164394654e-1, 0.70256257276544e-1, 0.21687043269488e-3, 0.13758331015182,
        -0.51501116343466e-1, -0.14865357483379, -0.38857100886810e-1, -0.29100433948943e-1, 0.14155684466279e-1}},
      {"argon",
       {0.85095714803969, -0.24003222943480e1, 0.54127841476466, 0.16919770692538e-1, 0.68825965019035e-1, 0.21428032815338e-3, 0.17429895321992,
        -0.33654495604194e-1, -0.13526799857691, -0.16387350791552e-1, -0.24987666851475e-1, 0.88769204815709e-2}}};
    return data;
}

/// Carbon dioxide, GERG-2004: own 22-term set. teqp GERG.hpp:552-559.
PureCoeffs carbondioxide_2004() {
    return {{0.52646564804653,    -0.14995725042592e1,  0.27329786733782,     0.12949500022786,     0.15404088341841,    -0.58186950946814,
             -0.18022494838296,   -0.95389904072812e-1, -0.80486819317679e-2, -0.35547751273090e-1, -0.28079014882405,   -0.82435890081677e-1,
             0.10832427979006e-1, -0.67073993161097e-2, -0.46827907600524e-2, -0.28359911832177e-1, 0.19500174744098e-1, -0.21609137507166,
             0.43772794926972,    -0.22130790113593,    0.15190189957331e-1,  -0.15380948953300e-1},
            {0.000, 1.250, 1.625, 0.375, 0.375,  1.375,  1.125,  1.375,  0.125,  1.625,  3.750,
             3.500, 7.500, 8.000, 6.000, 16.000, 11.000, 24.000, 26.000, 28.000, 24.000, 26.000},
            {1, 1, 2, 3, 3, 3, 4, 5, 6, 6, 1, 4, 1, 1, 3, 3, 4, 5, 5, 5, 5, 5},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 6, 6}};
}

/// Hydrogen, GERG-2004: own 14-term set. teqp GERG.hpp:562-568.
PureCoeffs hydrogen_2004() {
    return {{0.53579928451252e1, -0.62050252530595e1, 0.13830241327086, -0.71397954896129e-1, 0.15474053959733e-1, -0.14976806405771,
             -0.26368723988451e-1, 0.56681303156066e-1, -0.60063958030436e-1, -0.45043942027132, 0.42478840244500, -0.21997640827139e-1,
             -0.10499521374530e-1, -0.28955902866816e-2},
            {0.500, 0.625, 0.375, 0.625, 1.125, 2.625, 0.000, 0.250, 1.375, 4.000, 4.250, 5.000, 8.000, 8.000},
            {1, 1, 2, 2, 4, 1, 5, 5, 5, 1, 1, 2, 5, 1},
            {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 5}};
}

/// Water, GERG-2004: own 16-term set. teqp GERG.hpp:571-577.
PureCoeffs water_2004() {
    return {{0.82728408749586, -0.18602220416584e1, -0.11199009613744e1, 0.15635753976056, 0.87375844859025, -0.36674403715731, 0.53987893432436e-1,
             0.10957690214499e1, 0.53213037828563e-1, 0.13050533930825e-1, -0.41079520434476, 0.14637443344120, -0.55726838623719e-1,
             -0.11201774143800e-1, -0.66062758068099e-2, 0.46918522004538e-2},
            {0.500, 1.250, 1.875, 0.125, 1.500, 1.000, 0.750, 1.500, 0.625, 2.625, 5.000, 4.000, 4.500, 3.000, 4.000, 6.000},
            {1, 1, 1, 2, 2, 3, 4, 1, 5, 5, 1, 2, 4, 4, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 5, 5}};
}

/// Helium, GERG-2004: own 12-term set (distinct from main12_exponents).
/// teqp GERG.hpp:580-586.
PureCoeffs helium_2004() {
    return {{-0.45579024006737, 0.12516390754925e1, -0.15438231650621e1, 0.20467489707221e-1, -0.34476212380781, -0.20858459512787e-1,
             0.16227414711778e-1, -0.57471818200892e-1, 0.19462416430715e-1, -0.33295680123020e-1, -0.10863577372367e-1, -0.22173365245954e-1},
            {0.000, 0.125, 0.750, 1.000, 0.750, 2.625, 0.125, 1.250, 2.000, 1.000, 4.500, 5.000},
            {1, 1, 1, 4, 1, 3, 5, 5, 5, 2, 1, 2},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 3, 3}};
}

/// n-only table for members of the main12 family that GERG-2008 overrides
/// (carbon monoxide, isopentane) or adds (hydrogen sulfide, n-nonane,
/// n-decane).  Fluids not in this table fall through to GERG-2004.
/// teqp GERG.hpp:1110-1116.
const std::map<std::string, std::vector<double>>& n_main_2008_overrides() {
    static const std::map<std::string, std::vector<double>> data = {
      {"carbonmonoxide",
       {0.90554, -0.24515e1, 0.53149, 0.24173e-1, 0.72156e-1, 0.18818e-3, 0.19405, -0.43268e-1, -0.12778, -0.27896e-1, -0.34154e-1, 0.16329e-1}},
      {"isopentane",
       {0.10963e1, -0.30402e1, 0.10317e1, -0.15410, 0.11535, 0.29809e-3, 0.39571, -0.45881e-1, -0.35804, -0.10107, -0.35484e-1, 0.18156e-1}},
      {"hydrogensulfide",
       {0.87641, -0.20367e1, 0.21634, -0.50199e-1, 0.66994e-1, 0.19076e-3, 0.20227, -0.45348e-2, -0.22230, -0.34714e-1, -0.14885e-1, 0.74154e-2}},
      {"n-decane", {0.10461e1, -0.24807e1, 0.74372, -0.52579, 0.15315, 0.32865e-3, 0.84178, 0.55424e-1, -0.73555, -0.18507, -0.20775e-1, 0.12335e-1}},
      {"n-nonane",
       {0.11151e1, -0.27020e1, 0.83416, -0.38828, 0.13760, 0.28185e-3, 0.62037, 0.15847e-1, -0.61726, -0.15043, -0.12982e-1, 0.44325e-2}}};
    return data;
}

// Binary reducing-parameter tables (betaV, gammaV, betaT, gammaT).
//
// Transcribed from teqp (https://github.com/usnistgov/teqp),
// include/teqp/models/GERG/GERG.hpp: GERG2004::get_betasgammas (lines
// 595-766) and GERG2008::get_betasgammas (lines 1003-1103).  teqp keys these
// tables on std::pair<std::string,std::string> hashed with boost::hash; this
// backend uses std::map on the same pair type instead so as not to take a
// Boost dependency for a lookup table -- std::map orders pairs natively.
//
// Each row is stored EXACTLY ONCE, in the order teqp lists it.  The public
// accessor get_betasgammas() (below, outside this anonymous namespace) looks
// up (f1,f2) and, failing that, (f2,f1) with betaV/betaT reciprocated.

using BIPKey = std::pair<std::string, std::string>;

/// Table A3.8, GERG-2004 monograph.  153 pairs: every pair of the 18
/// GERG-2004 components, none missing, none extra.  teqp GERG.hpp:598-752.
const std::map<BIPKey, BetasGammas>& betasgammas_2004() {
    static const std::map<BIPKey, BetasGammas> data = {
      {{"methane", "nitrogen"}, {0.998721377, 1.013950311, 0.998098830, 0.979273013}},
      {{"methane", "carbondioxide"}, {0.999518072, 1.002806594, 1.022624490, 0.975665369}},
      {{"methane", "ethane"}, {0.997547866, 1.006617867, 0.996336508, 1.049707697}},
      {{"methane", "propane"}, {1.004827070, 1.038470657, 0.989680305, 1.098655531}},
      {{"methane", "n-butane"}, {0.979105972, 1.045375122, 0.994174910, 1.171607691}},
      {{"methane", "isobutane"}, {1.011240388, 1.054319053, 0.980315756, 1.161117729}},
      {{"methane", "n-pentane"}, {0.948330120, 1.124508039, 0.992127525, 1.249173968}},
      {{"methane", "isopentane"}, {1.000000000, 1.343685343, 1.000000000, 1.188899743}},
      {{"methane", "n-hexane"}, {0.958015294, 1.052643846, 0.981844797, 1.330570181}},
      {{"methane", "n-heptane"}, {0.962050831, 1.156655935, 0.977431529, 1.379850328}},
      {{"methane", "n-octane"}, {0.994740603, 1.116549372, 0.957473785, 1.449245409}},
      {{"methane", "hydrogen"}, {1.000000000, 1.018702573, 1.000000000, 1.352643115}},
      {{"methane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 0.950000000}},
      {{"methane", "carbonmonoxide"}, {0.997340772, 1.006102927, 0.987411732, 0.987473033}},
      {{"methane", "water"}, {1.012783169, 1.585018334, 1.063333913, 0.775810513}},
      {{"methane", "helium"}, {1.000000000, 0.881405683, 1.000000000, 3.159776855}},
      {{"methane", "argon"}, {1.034630259, 1.014678542, 0.990954281, 0.989843388}},
      {{"nitrogen", "carbondioxide"}, {0.977794634, 1.047578256, 1.005894529, 1.107654104}},
      {{"nitrogen", "ethane"}, {0.978880168, 1.042352891, 1.007671428, 1.098650964}},
      {{"nitrogen", "propane"}, {0.974424681, 1.081025408, 1.002677329, 1.201264026}},
      {{"nitrogen", "n-butane"}, {0.996082610, 1.146949309, 0.994515234, 1.304886838}},
      {{"nitrogen", "isobutane"}, {0.986415830, 1.100576129, 0.992868130, 1.284462634}},
      {{"nitrogen", "n-pentane"}, {1.000000000, 1.078877166, 1.000000000, 1.419029041}},
      {{"nitrogen", "isopentane"}, {1.000000000, 1.154135439, 1.000000000, 1.381770770}},
      {{"nitrogen", "n-hexane"}, {1.000000000, 1.195952177, 1.000000000, 1.472607971}},
      {{"nitrogen", "n-heptane"}, {1.000000000, 1.404554090, 1.000000000, 1.520975334}},
      {{"nitrogen", "n-octane"}, {1.000000000, 1.186067025, 1.000000000, 1.733280051}},
      {{"nitrogen", "hydrogen"}, {0.972532065, 0.970115357, 0.946134337, 1.175696583}},
      {{"nitrogen", "oxygen"}, {0.999521770, 0.997082328, 0.997190589, 0.995157044}},
      {{"nitrogen", "carbonmonoxide"}, {1.000000000, 1.008690943, 1.000000000, 0.993425388}},
      {{"nitrogen", "water"}, {1.000000000, 1.094749685, 1.000000000, 0.968808467}},
      {{"nitrogen", "helium"}, {0.969501055, 0.932629867, 0.692868765, 1.471831580}},
      {{"nitrogen", "argon"}, {1.004166412, 1.002212182, 0.999069843, 0.990034831}},
      {{"carbondioxide", "ethane"}, {1.002525718, 1.032876701, 1.013871147, 0.900949530}},
      {{"carbondioxide", "propane"}, {0.996898004, 1.047596298, 1.033620538, 0.908772477}},
      {{"carbondioxide", "n-butane"}, {1.174760923, 1.222437324, 1.018171004, 0.911498231}},
      {{"carbondioxide", "isobutane"}, {1.076551882, 1.081909003, 1.023339824, 0.929982936}},
      {{"carbondioxide", "n-pentane"}, {1.024311498, 1.068406078, 1.027000795, 0.979217302}},
      {{"carbondioxide", "isopentane"}, {1.060793104, 1.116793198, 1.019180957, 0.961218039}},
      {{"carbondioxide", "n-hexane"}, {1.000000000, 0.851343711, 1.000000000, 1.038675574}},
      {{"carbondioxide", "n-heptane"}, {1.205469976, 1.164585914, 1.011806317, 1.046169823}},
      {{"carbondioxide", "n-octane"}, {1.026169373, 1.104043935, 1.029690780, 1.074455386}},
      {{"carbondioxide", "hydrogen"}, {0.904142159, 1.152792550, 0.942320195, 1.782924792}},
      {{"carbondioxide", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"carbondioxide", "carbonmonoxide"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"carbondioxide", "water"}, {0.949055959, 1.542328793, 0.997372205, 0.775453996}},
      {{"carbondioxide", "helium"}, {0.846647561, 0.864141549, 0.768377630, 3.207456948}},
      {{"carbondioxide", "argon"}, {1.008392428, 1.029205465, 0.996512863, 1.050971635}},
      {{"ethane", "propane"}, {0.997607277, 1.003034720, 0.996199694, 1.014730190}},
      {{"ethane", "n-butane"}, {0.999157205, 1.006179146, 0.999130554, 1.034832749}},
      {{"ethane", "isobutane"}, {1.000000000, 1.006616886, 1.000000000, 1.033283811}},
      {{"ethane", "n-pentane"}, {0.993851009, 1.026085655, 0.998688946, 1.066665676}},
      {{"ethane", "isopentane"}, {1.000000000, 1.045439246, 1.000000000, 1.021150247}},
      {{"ethane", "n-hexane"}, {1.000000000, 1.169701102, 1.000000000, 1.092177796}},
      {{"ethane", "n-heptane"}, {1.000000000, 1.057666085, 1.000000000, 1.134532014}},
      {{"ethane", "n-octane"}, {1.007469726, 1.071917985, 0.984068272, 1.168636194}},
      {{"ethane", "hydrogen"}, {0.925367171, 1.106072040, 0.932969831, 1.902008495}},
      {{"ethane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"ethane", "carbonmonoxide"}, {1.000000000, 1.201417898, 1.000000000, 1.069224728}},
      {{"ethane", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"ethane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"ethane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"propane", "n-butane"}, {0.999795868, 1.003264179, 1.000310289, 1.007392782}},
      {{"propane", "isobutane"}, {0.999243146, 1.001156119, 0.998012298, 1.005250774}},
      {{"propane", "n-pentane"}, {1.044919431, 1.019921513, 0.996484021, 1.008344412}},
      {{"propane", "isopentane"}, {1.040459289, 0.999432118, 0.994364425, 1.003269500}},
      {{"propane", "n-hexane"}, {1.000000000, 1.057872566, 1.000000000, 1.025657518}},
      {{"propane", "n-heptane"}, {1.000000000, 1.079648053, 1.000000000, 1.050044169}},
      {{"propane", "n-octane"}, {1.000000000, 1.102764612, 1.000000000, 1.063694129}},
      {{"propane", "hydrogen"}, {1.000000000, 1.074006110, 1.000000000, 2.308215191}},
      {{"propane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"propane", "carbonmonoxide"}, {1.000000000, 1.108143673, 1.000000000, 1.197564208}},
      {{"propane", "water"}, {1.000000000, 1.011759763, 1.000000000, 0.600340961}},
      {{"propane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"propane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-butane", "isobutane"}, {1.000880464, 1.000414440, 1.000077547, 1.001432824}},
      {{"n-butane", "n-pentane"}, {1.000000000, 1.018159650, 1.000000000, 1.002143640}},
      {{"n-butane", "isopentane"}, {1.000000000, 1.002728262, 1.000000000, 1.000792201}},
      {{"n-butane", "n-hexane"}, {1.000000000, 1.034995284, 1.000000000, 1.009157060}},
      {{"n-butane", "n-heptane"}, {1.000000000, 1.019174227, 1.000000000, 1.021283378}},
      {{"n-butane", "n-octane"}, {1.000000000, 1.046905515, 1.000000000, 1.033180106}},
      {{"n-butane", "hydrogen"}, {1.000000000, 1.232939523, 1.000000000, 2.509259945}},
      {{"n-butane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-butane", "carbonmonoxide"}, {1.000000000, 1.084740904, 1.000000000, 1.174055065}},
      {{"n-butane", "water"}, {1.000000000, 1.223638763, 1.000000000, 0.615512682}},
      {{"n-butane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-butane", "argon"}, {1.000000000, 1.214638734, 1.000000000, 1.245039498}},
      {{"isobutane", "n-pentane"}, {1.000000000, 1.002779804, 1.000000000, 1.002495889}},
      {{"isobutane", "isopentane"}, {1.000000000, 1.002284197, 1.000000000, 1.001835788}},
      {{"isobutane", "n-hexane"}, {1.000000000, 1.010493989, 1.000000000, 1.006018054}},
      {{"isobutane", "n-heptane"}, {1.000000000, 1.021668316, 1.000000000, 1.009885760}},
      {{"isobutane", "n-octane"}, {1.000000000, 1.032807063, 1.000000000, 1.013945424}},
      {{"isobutane", "hydrogen"}, {1.000000000, 1.147595688, 1.000000000, 1.895305393}},
      {{"isobutane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isobutane", "carbonmonoxide"}, {1.000000000, 1.087272232, 1.000000000, 1.161523504}},
      {{"isobutane", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isobutane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isobutane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-pentane", "isopentane"}, {1.000000000, 1.000024352, 1.000000000, 1.000050537}},
      {{"n-pentane", "n-hexane"}, {1.000000000, 1.002480637, 1.000000000, 1.000761237}},
      {{"n-pentane", "n-heptane"}, {1.000000000, 1.008972412, 1.000000000, 1.002441051}},
      {{"n-pentane", "n-octane"}, {1.000000000, 1.069223964, 1.000000000, 1.016422347}},
      {{"n-pentane", "hydrogen"}, {1.000000000, 1.188334783, 1.000000000, 2.013859174}},
      {{"n-pentane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-pentane", "carbonmonoxide"}, {1.000000000, 1.119954454, 1.000000000, 1.206195595}},
      {{"n-pentane", "water"}, {1.000000000, 0.956677310, 1.000000000, 0.447666011}},
      {{"n-pentane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-pentane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isopentane", "n-hexane"}, {1.000000000, 1.002996055, 1.000000000, 1.001204174}},
      {{"isopentane", "n-heptane"}, {1.000000000, 1.009928531, 1.000000000, 1.003194615}},
      {{"isopentane", "n-octane"}, {1.000000000, 1.017880981, 1.000000000, 1.005647480}},
      {{"isopentane", "hydrogen"}, {1.000000000, 1.184339122, 1.000000000, 1.996386669}},
      {{"isopentane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isopentane", "carbonmonoxide"}, {1.000000000, 1.116693501, 1.000000000, 1.199475627}},
      {{"isopentane", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isopentane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"isopentane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-hexane", "n-heptane"}, {1.000000000, 1.001508227, 1.000000000, 0.999762786}},
      {{"n-hexane", "n-octane"}, {1.000000000, 1.006268954, 1.000000000, 1.001633952}},
      {{"n-hexane", "hydrogen"}, {1.000000000, 1.243461678, 1.000000000, 3.021197546}},
      {{"n-hexane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-hexane", "carbonmonoxide"}, {1.000000000, 1.155145836, 1.000000000, 1.233435828}},
      {{"n-hexane", "water"}, {1.000000000, 1.170217596, 1.000000000, 0.569681333}},
      {{"n-hexane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-hexane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-heptane", "n-octane"}, {1.000000000, 1.006767176, 1.000000000, 0.998793111}},
      {{"n-heptane", "hydrogen"}, {1.000000000, 1.159131722, 1.000000000, 3.169143057}},
      {{"n-heptane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-heptane", "carbonmonoxide"}, {1.000000000, 1.190354273, 1.000000000, 1.256295219}},
      {{"n-heptane", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-heptane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-heptane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-octane", "hydrogen"}, {1.000000000, 1.305249405, 1.000000000, 2.191555216}},
      {{"n-octane", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-octane", "carbonmonoxide"}, {1.000000000, 1.219206702, 1.000000000, 1.276744779}},
      {{"n-octane", "water"}, {1.000000000, 0.599484191, 1.000000000, 0.662072469}},
      {{"n-octane", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"n-octane", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"hydrogen", "oxygen"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"hydrogen", "carbonmonoxide"}, {1.000000000, 1.121416201, 1.000000000, 1.377504607}},
      {{"hydrogen", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"hydrogen", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"hydrogen", "argon"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"oxygen", "carbonmonoxide"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"oxygen", "water"}, {1.000000000, 1.143174289, 1.000000000, 0.964767932}},
      {{"oxygen", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"oxygen", "argon"}, {0.999746847, 0.993907223, 1.000023103, 0.990430423}},
      {{"carbonmonoxide", "water"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"carbonmonoxide", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"carbonmonoxide", "argon"}, {1.000000000, 1.159720623, 1.000000000, 0.954215746}},
      {{"water", "helium"}, {1.000000000, 1.000000000, 1.000000000, 1.000000000}},
      {{"water", "argon"}, {1.000000000, 1.038993495, 1.000000000, 1.070941866}},
      {{"helium", "argon"}, {1.00000000, 1.00000000, 1.00000000, 1.00000000}},
    };
    return data;
}

/// Table A8, GERG-2008 manuscript.  Pairs unchanged from GERG-2004 are left
/// out and fall through to betasgammas_2004().  72 pairs: 15 that GERG-2008
/// revises (the EOS changed for isopentane and CO, so their estimated
/// interaction parameters were recalculated with the "new" Tc) plus 57 new
/// pairs involving the three fluids GERG-2008 adds (hydrogen sulfide,
/// n-nonane, n-decane: 3 new x 18 old + 3 pairs among the new three).  teqp
/// GERG.hpp:1008-1088.
const std::map<BIPKey, BetasGammas>& betasgammas_2008_overrides() {
    static const std::map<BIPKey, BetasGammas> data = {
      // This set has different values than in GERG-2004. The change
      // occurs because the EOS has changed for isopentane and CO and thus
      // the estimated interaction parameters need to be calculated with
      // the "new" Tc.
      {{"ethane", "isopentane"}, {1.0, 1.045439935, 1.0, 1.021150247}},
      {{"n-butane", "carbonmonoxide"}, {1.0, 1.084740904, 1.0, 1.173916162}},
      {{"n-butane", "isopentane"}, {1.0, 1.002728434, 1.0, 1.000792201}},
      {{"isobutane", "carbonmonoxide"}, {1.0, 1.087272232, 1.0, 1.161390082}},
      {{"isobutane", "isopentane"}, {1.0, 1.002284353, 1.0, 1.001835788}},
      {{"n-pentane", "carbonmonoxide"}, {1.0, 1.119954454, 1.0, 1.206043295}},
      {{"isopentane", "carbonmonoxide"}, {1.0, 1.116694577, 1.0, 1.199326059}},
      {{"n-hexane", "carbonmonoxide"}, {1.0, 1.155145836, 1.0, 1.233272781}},
      {{"n-heptane", "carbonmonoxide"}, {1.0, 1.190354273, 1.0, 1.256123503}},
      {{"n-octane", "carbonmonoxide"}, {1.0, 1.219206702, 1.0, 1.276565536}},
      {{"n-pentane", "isopentane"}, {1.0, 1.000024335, 1.0, 1.000050537}},
      {{"isopentane", "n-hexane"}, {1.0, 1.002995876, 1.0, 1.001204174}},
      {{"isopentane", "n-heptane"}, {1.0, 1.009928206, 1.0, 1.003194615}},
      {{"isopentane", "n-octane"}, {1.0, 1.017880545, 1.0, 1.00564748}},
      {{"isopentane", "hydrogen"}, {1.0, 1.184340443, 1.0, 1.996386669}},

      // The 57 new pairs added in GERG-2008.
      {{"methane", "n-nonane"}, {1.002852287, 1.141895355, 0.947716769, 1.528532478}},
      {{"methane", "n-decane"}, {1.033086292, 1.146089637, 0.937777823, 1.568231489}},
      {{"methane", "hydrogensulfide"}, {1.012599087, 1.040161207, 1.011090031, 0.961155729}},
      {{"nitrogen", "n-nonane"}, {1.0, 1.100405929, 0.95637945, 1.749119996}},
      {{"nitrogen", "n-decane"}, {1.0, 1.0, 0.957934447, 1.822157123}},
      {{"nitrogen", "hydrogensulfide"}, {0.910394249, 1.256844157, 1.004692366, 0.9601742}},
      {{"carbondioxide", "n-nonane"}, {1.0, 0.973386152, 1.00768862, 1.140671202}},
      {{"carbondioxide", "n-decane"}, {1.000151132, 1.183394668, 1.02002879, 1.145512213}},
      {{"carbondioxide", "hydrogensulfide"}, {0.906630564, 1.024085837, 1.016034583, 0.92601888}},
      {{"ethane", "n-nonane"}, {1.0, 1.14353473, 1.0, 1.05603303}},
      {{"ethane", "n-decane"}, {0.995676258, 1.098361281, 0.970918061, 1.237191558}},
      {{"ethane", "hydrogensulfide"}, {1.010817909, 1.030988277, 0.990197354, 0.90273666}},
      {{"propane", "n-nonane"}, {1.0, 1.199769134, 1.0, 1.109973833}},
      {{"propane", "n-decane"}, {0.984104227, 1.053040574, 0.985331233, 1.140905252}},
      {{"propane", "hydrogensulfide"}, {0.936811219, 1.010593999, 0.992573556, 0.905829247}},
      {{"n-butane", "n-nonane"}, {1.0, 1.049219137, 1.0, 1.014096448}},
      {{"n-butane", "n-decane"}, {0.976951968, 1.027845529, 0.993688386, 1.076466918}},
      {{"n-butane", "hydrogensulfide"}, {0.908113163, 1.033366041, 0.985962886, 0.926156602}},
      {{"isobutane", "n-nonane"}, {1.0, 1.047298475, 1.0, 1.017817492}},
      {{"isobutane", "n-decane"}, {1.0, 1.060243344, 1.0, 1.021624748}},
      {{"isobutane", "hydrogensulfide"}, {1.012994431, 0.988591117, 0.974550548, 0.937130844}},
      {{"n-pentane", "n-nonane"}, {1.0, 1.034910633, 1.0, 1.103421755}},
      {{"n-pentane", "n-decane"}, {1.0, 1.016370338, 1.0, 1.049035838}},
      {{"n-pentane", "hydrogensulfide"}, {0.984613203, 1.076539234, 0.962006651, 0.959065662}},
      {{"isopentane", "n-nonane"}, {1.0, 1.028994325, 1.0, 1.008191499}},
      {{"isopentane", "n-decane"}, {1.0, 1.039372957, 1.0, 1.010825138}},
      {{"isopentane", "hydrogensulfide"}, {1.0, 0.835763343, 1.0, 0.982651529}},
      {{"n-hexane", "n-nonane"}, {1.0, 1.02076168, 1.0, 1.055369591}},
      {{"n-hexane", "n-decane"}, {1.001516371, 1.013511439, 0.99764101, 1.028939539}},
      {{"n-hexane", "hydrogensulfide"}, {0.754473958, 1.339283552, 0.985891113, 0.956075596}},
      {{"n-heptane", "n-nonane"}, {1.0, 1.001370076, 1.0, 1.001150096}},
      {{"n-heptane", "n-decane"}, {1.0, 1.002972346, 1.0, 1.002229938}},
      {{"n-heptane", "hydrogensulfide"}, {0.828967164, 1.087956749, 0.988937417, 1.013453092}},
      {{"n-octane", "n-nonane"}, {1.0, 1.001357085, 1.0, 1.000235044}},
      {{"n-octane", "n-decane"}, {1.0, 1.002553544, 1.0, 1.007186267}},
      {{"n-octane", "hydrogensulfide"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-nonane", "n-decane"}, {1.0, 1.00081052, 1.0, 1.000182392}},
      {{"n-nonane", "hydrogen"}, {1.0, 1.342647661, 1.0, 2.23435404}},
      {{"n-nonane", "oxygen"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-nonane", "carbonmonoxide"}, {1.0, 1.252151449, 1.0, 1.294070556}},
      {{"n-nonane", "water"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-nonane", "hydrogensulfide"}, {1.0, 1.082905109, 1.0, 1.086557826}},
      {{"n-nonane", "helium"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-nonane", "argon"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-decane", "hydrogen"}, {1.695358382, 1.120233729, 1.064818089, 3.786003724}},
      {{"n-decane", "oxygen"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-decane", "carbonmonoxide"}, {1.0, 0.87018496, 1.049594632, 1.803567587}},
      {{"n-decane", "water"}, {1.0, 0.551405318, 0.897162268, 0.740416402}},
      {{"n-decane", "hydrogensulfide"}, {0.975187766, 1.171714677, 0.973091413, 1.103693489}},
      {{"n-decane", "helium"}, {1.0, 1.0, 1.0, 1.0}},
      {{"n-decane", "argon"}, {1.0, 1.0, 1.0, 1.0}},
      {{"hydrogen", "hydrogensulfide"}, {1.0, 1.0, 1.0, 1.0}},
      {{"oxygen", "hydrogensulfide"}, {1.0, 1.0, 1.0, 1.0}},
      {{"carbonmonoxide", "hydrogensulfide"}, {0.795660392, 1.101731308, 1.025536736, 1.022749748}},
      {{"water", "hydrogensulfide"}, {1.0, 1.014832832, 1.0, 0.940587083}},
      {{"hydrogensulfide", "helium"}, {1.0, 1.0, 1.0, 1.0}},
      {{"hydrogensulfide", "argon"}, {1.0, 1.0, 1.0, 1.0}},
    };
    return data;
}

// Ideal-gas coefficient tables, GERG-2004 monograph Table A3.1.
//
// Stored exactly as teqp stores them: {n0[1..7], theta0[4..7]}, i.e. 7 n and
// 4 theta with no padding, so the literals line up one-for-one with teqp
// GERG.hpp:470-487 (GERG-2004) and :1135-1141 (GERG-2008).  The padding to
// length 8 that makes the monograph's 1-based indices work is applied in
// pad_alphaig() below, just as teqp does on the way out of its accessor.
//
// n0[1] and n0[2] as tabulated here are the PUBLISHED integration constants.
// get_alphaig_coeffs discards them; see recalc_integration_constants.

using RawAlphaig = std::pair<std::vector<double>, std::vector<double>>;

/// teqp GERG.hpp:470-487 (GERG2004::get_alphaig_coeffs dict).
const std::map<std::string, RawAlphaig>& alphaig_2004() {
    static const std::map<std::string, RawAlphaig> data = {
      {"methane",
       {{19.597538587, -83.959667892, 3.000880, 0.763150, 0.00460, 8.744320, -4.469210000}, {4.306474465, 0.936220902, 5.577233895, 5.722644361}}},
      {"nitrogen", {{11.083437707, -22.202102428, 2.500310, 0.137320, -0.14660, 0.900660, 0}, {5.251822620, -5.393067706, 13.788988208, 0}}},
      {"carbondioxide",
       {{11.925182741, -16.118762264, 2.500020, 2.044520, -1.060440, 2.033660, 0.013930000}, {3.022758166, -2.844425476, 1.589964364, 1.121596090}}},
      {"ethane",
       {{24.675465518, -77.425313760, 3.002630, 4.339390, 1.237220, 13.19740, -6.019890000}, {1.831882406, 0.731306621, 3.378007481, 3.508721939}}},
      {"propane",
       {{31.602934734, -84.463284382, 3.029390, 6.605690, 3.1970, 19.19210, -8.372670000}, {1.297521801, 0.543210978, 2.583146083, 2.777773271}}},
      {"n-butane",
       {{20.884168790, -91.638478026, 3.339440, 9.448930, 6.894060, 24.46180, 14.782400000}, {1.101487798, 0.431957660, 4.502440459, 2.124516319}}},
      {"isobutane",
       {{20.413751434, -94.467620036, 3.067140, 8.975750, 5.251560, 25.14230, 16.138800000}, {1.074673199, 0.485556021, 4.671261865, 2.191583480}}},
      {"n-pentane", {{14.536635738, -89.919548319, 3.0, 8.950430, 21.8360, 33.40320, 0}, {0.380391739, 1.789520971, 3.777411113, 0}}},
      {"isopentane", {{15.449937973, -101.298172792, 3.0, 11.76180, 20.11010, 33.16880, 0}, {0.635392636, 1.977271641, 4.169371131, 0}}},
      {"n-hexane", {{14.345993081, -96.165722367, 3.0, 11.69770, 26.81420, 38.61640, 0}, {0.359036667, 1.691951873, 3.596924107, 0}}},
      {"n-heptane", {{15.063809621, -97.345252349, 3.0, 13.72660, 30.47070, 43.55610, 0}, {0.314348398, 1.548136560, 3.259326458, 0}}},
      {"n-octane", {{15.864709639, -97.370667555, 3.0, 15.68650, 33.80290, 48.17310, 0}, {0.279143540, 1.431644769, 2.973845992, 0}}},
      {"hydrogen",
       {{13.796474934, -175.864487294, 1.479060, 0.958060, 0.454440, 1.560390, -1.375600000},
        {6.891654113, 9.847634830, 49.765290750, 50.367279301}}},
      {"oxygen", {{10.001874708, -14.996095135, 2.501460, 1.075580, 1.013340, 0, 0}, {14.461722565, 7.223325463, 0, 0}}},
      {"carbonmonoxide", {{10.814500335, -19.843695435, 2.500550, 1.028650, 0.004930, 0, 0}, {11.675075301, 5.305158133, 0, 0}}},
      {"water", {{8.203553050, -11.996306443, 3.003920, 0.010590, 0.987630, 3.069040, 0}, {0.415386589, 1.763895929, 3.874803739, 0}}},
      {"helium", {{13.628441975, -143.470759602, 1.5, 0, 0, 0, 0}, {0, 0, 0, 0}}},
      {"argon", {{8.316662546, -4.946502600, 1.50, 0, 0, 0, 0}, {0, 0, 0, 0}}}};
    return data;
}

/// Entries GERG-2008 changes or adds relative to GERG-2004; everything else
/// falls through to alphaig_2004().  teqp GERG.hpp:1135-1141.
const std::map<std::string, RawAlphaig>& alphaig_2008_overrides() {
    static const std::map<std::string, RawAlphaig> data = {
      {"carbonmonoxide",
       {{10.813340744, -19.834733959, 2.50055, 1.02865, 0.00493, 0, 0}, {11.669802800, 5.302762306, 0, 0}}},  // changed in GERG-2008
      {"isopentane",
       {{15.449907693, -101.298172792, 3.0, 11.76180, 20.11010, 33.16880, 0}, {0.635392636, 1.977271641, 4.169371131, 0}}},  // changed in GERG-2008
      {"n-nonane", {{16.313913248, -102.160247463, 3.0, 18.02410, 38.12350, 53.34150, 0}, {0.263819696, 1.370586158, 2.848860483, 0}}},
      {"n-decane", {{15.870791919, -108.858547525, 3.0, 21.00690, 43.49310, 58.36570, 0}, {0.267034159, 1.353835195, 2.833479035, 0}}},
      {"hydrogensulfide", {{9.336197742, -16.266508995, 3.0, 3.11942, 1.00243, 0, 0}, {4.914580541, 2.270653980, 0, 0}}}};
    return data;
}

/// Zero-pad the raw {7 n, 4 theta} tables up to the monograph's 1-based
/// indexing: n0[1..7], theta0[4..7].  teqp GERG.hpp:497-506.
AlphaigCoeffs pad_alphaig(const std::string& gerg_name, const RawAlphaig& raw) {
    if (raw.first.size() != 7) {
        throw ValueError(format("[%s] does not have 7 n coefficients in ideal gas", gerg_name.c_str()));
    }
    if (raw.second.size() != 4) {
        throw ValueError(format("[%s] does not have 4 theta coefficients in ideal gas", gerg_name.c_str()));
    }
    AlphaigCoeffs c;
    c.n0 = raw.first;
    c.n0.insert(c.n0.begin(), 0.0);  // 0-pad so that indexing matches GERG-2004
    c.theta0 = raw.second;
    c.theta0.insert(c.theta0.begin(), 4, 0.0);  // 0-pad so that indexing matches GERG-2004
    return c;
}

// Departure-function tables: F_ij scaling factors and departure coefficients.
//
// Transcribed from teqp (https://github.com/usnistgov/teqp),
// include/teqp/models/GERG/GERG.hpp: get_Fij (lines 768-802, Table A3.6) and
// get_departurecoeffs (lines 805-900), both in namespace GERG2004.
// GERG-2008 reuses both unchanged (teqp GERG.hpp:999-1001, `using
// GERG2004::get_Fij; using GERG2004::get_departurecoeffs;`), so there is a
// single shared table here rather than a 2004/2008 override pair like the
// tables above.
//
// Exactly 15 of the 153 possible GERG-2004 pairs carry a departure function.
// 7 pairs (methane/nitrogen, methane/carbondioxide, methane/ethane,
// methane/propane, nitrogen/carbondioxide, nitrogen/ethane, methane/hydrogen)
// have fluid-pair-specific coefficients and F_ij == 1 exactly. The other 8
// pairs share one "generalized" departure function, scaled per-pair by
// F_ij; 7 of those 8 have F_ij != 1, and the 8th (methane/n-butane) also
// uses the generalized form but happens to have F_ij == 1.

/// Table A3.6, GERG-2004 monograph.  15 pairs.  teqp GERG.hpp:771-787.
/// SYMMETRIC: F_ij == F_ji exactly.  Unlike BetasGammas above, nothing is
/// reciprocated on a reversed-order lookup.
const std::map<BIPKey, double>& fij_table() {
    static const std::map<BIPKey, double> data = {
      {{"methane", "nitrogen"}, 1.0},
      {{"methane", "carbondioxide"}, 1.0},
      {{"methane", "ethane"}, 1.0},
      {{"methane", "propane"}, 1.0},
      {{"methane", "n-butane"}, 1.0},
      {{"methane", "isobutane"}, 0.771035405688},
      {{"methane", "hydrogen"}, 1.0},
      {{"nitrogen", "carbondioxide"}, 1.0},
      {{"nitrogen", "ethane"}, 1.0},
      {{"ethane", "propane"}, 0.130424765150},
      {{"ethane", "n-butane"}, 0.281570073085},
      {{"ethane", "isobutane"}, 0.260632376098},
      {{"propane", "n-butane"}, 0.0312572600489},
      {{"propane", "isobutane"}, -0.0551609771024},
      {{"n-butane", "isobutane"}, -0.0551240293009},
    };
    return data;
}

/// teqp GERG.hpp:818-827.
DepartureCoeffs departure_methane_nitrogen() {
    DepartureCoeffs dc;
    dc.n = {-0.98038985517335e-2, 0.42487270143005e-3, -0.34800214576142e-1, -0.13333813013896, -0.11993694974627e-1,
            0.69243379775168e-1,  -0.31022508148249,   0.24495491753226,     0.22369816716981};
    dc.d = {1, 4, 1, 2, 2, 2, 2, 2, 3};
    dc.t = {0.000, 1.850, 7.850, 5.400, 0.000, 0.750, 2.800, 4.450, 4.250};
    dc.eta = {0, 0, 1.000, 1.000, 0.250, 0.000, 0.000, 0.000, 0.000};
    dc.epsilon = {0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    dc.beta = {0, 0, 1.0, 1.0, 2.5, 3.0, 3.0, 3.0, 3.0};
    dc.gamma = {0, 0, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:828-837.
DepartureCoeffs departure_methane_carbondioxide() {
    DepartureCoeffs dc;
    dc.n = {-0.10859387354942, 0.80228576727389e-1, -0.93303985115717e-2, 0.40989274005848e-1, -0.24338019772494, 0.23855347281124};
    dc.d = {1, 2, 3, 1, 2, 3};
    dc.t = {2.600, 1.950, 0.000, 3.950, 7.950, 8.000};
    dc.eta = {0, 0, 0, 1.000, 0.500, 0.000};
    dc.epsilon = {0, 0, 0, 0.5, 0.5, 0.5};
    dc.beta = {0, 0, 0, 1.0, 2.0, 3.0};
    dc.gamma = {0, 0, 0, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:838-847 (there stored as sortpair("ethane","methane")).
DepartureCoeffs departure_methane_ethane() {
    DepartureCoeffs dc;
    dc.n = {-0.80926050298746e-3, -0.75381925080059e-3, -0.41618768891219e-1, -0.23452173681569,  0.14003840584586,    0.63281744807738e-1,
            -0.34660425848809e-1, -0.23918747334251,    0.19855255066891e-2,  0.61777746171555e1, -0.69575358271105e1, 0.10630185306388e1};
    dc.t = {0.650, 1.550, 3.100, 5.900, 7.050, 3.350, 1.200, 5.800, 2.700, 0.450, 0.550, 1.950};
    dc.d = {3, 4, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3};
    dc.eta = {0, 0, 1.000, 1.000, 1.000, 0.875, 0.750, 0.500, 0.000, 0.000, 0.000, 0.000};
    dc.epsilon = {0, 0, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500};
    dc.beta = {0, 0, 1.000, 1.000, 1.000, 1.250, 1.500, 2.000, 3.000, 3.000, 3.000, 3.000};
    dc.gamma = {0, 0, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:848-857 (there stored as sortpair("propane","methane")).
DepartureCoeffs departure_methane_propane() {
    DepartureCoeffs dc;
    dc.n = {0.13746429958576e-1, -0.74425012129552e-2, -0.45516600213685e-2, -0.54546603350237e-2, 0.23682016824471e-2,
            0.18007763721438,    -0.44773942932486,    0.19327374888200e-1,  -0.30632197804624};
    dc.t = {1.850, 3.950, 0.000, 1.850, 3.850, 5.250, 3.850, 0.200, 6.500};
    dc.d = {3, 3, 4, 4, 4, 1, 1, 1, 2};
    dc.eta = {0, 0, 0, 0, 0, 0.250, 0.250, 0.000, 0.000};
    dc.epsilon = {0, 0, 0, 0, 0, 0.500, 0.500, 0.500, 0.500};
    dc.beta = {0, 0, 0, 0, 0, 0.750, 1.000, 2.000, 3.000};
    dc.gamma = {0, 0, 0, 0, 0, 0.500, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:858-867.
DepartureCoeffs departure_nitrogen_carbondioxide() {
    DepartureCoeffs dc;
    dc.n = {0.28661625028399, -0.10919833861247, -0.11374032082270e1, 0.76580544237358, 0.42638000926819e-2, 0.17673538204534};
    dc.d = {2, 3, 1, 1, 1, 2};
    dc.t = {1.850, 1.400, 3.200, 2.500, 8.000, 3.750};
    dc.eta = {0, 0, 0.250, 0.250, 0.000, 0.000};
    dc.epsilon = {0, 0, 0.500, 0.500, 0.500, 0.500};
    dc.beta = {0, 0, 0.750, 1.000, 2.000, 3.000};
    dc.gamma = {0, 0, 0.500, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:868-877.
DepartureCoeffs departure_nitrogen_ethane() {
    DepartureCoeffs dc;
    dc.n = {-0.47376518126608, 0.48961193461001, -0.57011062090535e-2, -0.19966820041320, -0.69411103101723, 0.69226192739021};
    dc.d = {2, 2, 3, 1, 2, 2};
    dc.t = {0.000, 0.050, 0.000, 3.650, 4.900, 4.450};
    dc.eta = {0, 0, 0, 1.000, 1.000, 0.875};
    dc.epsilon = {0, 0, 0, 0.500, 0.500, 0.500};
    dc.beta = {0, 0, 0, 1.000, 1.000, 1.250};
    dc.gamma = {0, 0, 0, 0.500, 0.500, 0.500};
    return dc;
}

/// teqp GERG.hpp:878-887.
DepartureCoeffs departure_methane_hydrogen() {
    DepartureCoeffs dc;
    dc.n = {-0.25157134971934, -0.62203841111983e-2, 0.88850315184396e-1, -0.35592212573239e-1};
    dc.t = {2.000, -1.000, 1.750, 1.400};
    dc.d = {1, 3, 3, 4};
    dc.eta = {0, 0, 0, 0};
    dc.epsilon = {0, 0, 0, 0};
    dc.beta = {0, 0, 0, 0};
    dc.gamma = {0, 0, 0, 0};
    return dc;
}

/// The 7 pairs above, each with its own fluid-pair-specific departure
/// coefficients.  Looked up in either order by get_departurecoeffs; no
/// reciprocation applies to any field here (unlike BetasGammas).
const std::map<BIPKey, DepartureCoeffs>& departure_specific_table() {
    static const std::map<BIPKey, DepartureCoeffs> data = {
      {{"methane", "nitrogen"}, departure_methane_nitrogen()},
      {{"methane", "carbondioxide"}, departure_methane_carbondioxide()},
      {{"methane", "ethane"}, departure_methane_ethane()},
      {{"methane", "propane"}, departure_methane_propane()},
      {{"nitrogen", "carbondioxide"}, departure_nitrogen_carbondioxide()},
      {{"nitrogen", "ethane"}, departure_nitrogen_ethane()},
      {{"methane", "hydrogen"}, departure_methane_hydrogen()},
    };
    return data;
}

/// The "generalized" departure function, teqp GERG.hpp:888-896.  Shared,
/// unmodified, by all 8 pairs in generalized_departure_pairs() below; each
/// pair scales it by its own F_ij (fij_table() above).  Every term here is
/// polynomial (eta == beta == 0 for all 10 terms): this departure function
/// has no Gaussian block at all.
DepartureCoeffs generalized_departure() {
    DepartureCoeffs dc;
    dc.n = {0.25574776844118e1, -0.79846357136353e1, 0.47859131465806e1, -0.73265392369587, 1.3805471345312,
            0.28349603476365,   -0.49087385940425,   -0.10291888921447,  0.11836314681968,  0.55527385721943e-4};
    dc.d = {1, 1, 1, 2, 2, 3, 3, 4, 4, 4};
    dc.t = {1.000, 1.550, 1.700, 0.250, 1.350, 0.000, 1.250, 0.000, 0.700, 5.400};
    dc.eta = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    dc.epsilon = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    dc.beta = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    dc.gamma = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    return dc;
}

/// teqp GERG.hpp:811-815 (the `generalized` std::set, sorted-pair form).
/// 8 pairs; 7 of them have F_ij != 1 (fij_table() above), the 8th
/// (methane/n-butane) has F_ij == 1 but still uses this generalized form
/// rather than the corresponding-states term alone.
const std::vector<BIPKey>& generalized_departure_pairs() {
    static const std::vector<BIPKey> pairs = {
      {"methane", "n-butane"}, {"methane", "isobutane"}, {"ethane", "propane"},    {"ethane", "n-butane"},
      {"ethane", "isobutane"}, {"propane", "n-butane"},  {"propane", "isobutane"}, {"n-butane", "isobutane"},
    };
    return pairs;
}

bool is_generalized_departure_pair(const std::string& f1, const std::string& f2) {
    const auto& pairs = generalized_departure_pairs();
    return std::find(pairs.begin(), pairs.end(), BIPKey{f1, f2}) != pairs.end()
           || std::find(pairs.begin(), pairs.end(), BIPKey{f2, f1}) != pairs.end();
}

}  // namespace

PureCoeffs get_pure_coeffs(GERGModel model, const std::string& gerg_name) {
    if (model == GERGModel::GERG_2008) {
        const auto& overrides = n_main_2008_overrides();
        auto it = overrides.find(gerg_name);
        if (it != overrides.end()) {
            const auto& e = main12_exponents();
            return PureCoeffs{it->second, e.t, e.d, e.c, e.l};
        }
    }
    // Fall through to GERG-2004 (also GERG-2004's own lookup path).
    if (gerg_name == "carbondioxide") {
        return carbondioxide_2004();
    }
    if (gerg_name == "hydrogen") {
        return hydrogen_2004();
    }
    if (gerg_name == "water") {
        return water_2004();
    }
    if (gerg_name == "helium") {
        return helium_2004();
    }
    {
        const auto& mne = n_mne_2004();
        auto it = mne.find(gerg_name);
        if (it != mne.end()) {
            const auto& e = mne24_exponents();
            return PureCoeffs{it->second, e.t, e.d, e.c, e.l};
        }
    }
    {
        const auto& main = n_main_2004();
        auto it = main.find(gerg_name);
        if (it != main.end()) {
            const auto& e = main12_exponents();
            return PureCoeffs{it->second, e.t, e.d, e.c, e.l};
        }
    }
    throw ValueError(format("Unable to load GERG pure residual coefficients for [%s]", gerg_name.c_str()));
}

std::pair<double, double> recalc_integration_constants(const AlphaigCoeffs& c, double T0, double Tci, double rho0, double rhoci, double Rstar_R) {
    // Faithful port of teqp GERG.hpp:49-74.  teqp builds two 3-element rows
    // {coefficient of n0[1], coefficient of n0[2], everything else} for the
    // reduced ideal-gas Helmholtz energy Aig00 and its tau-derivative Aig10,
    // then solves the resulting 2x2 system with Eigen.  A 2x2 does not need
    // Eigen, so it is written out by hand below; the row construction is kept
    // literally identical so the two implementations stay diffable.
    const double th = Tci / T0;
    auto sinh_term = [&](std::size_t i) { return (c.n0[i] != 0) ? c.n0[i] * std::log(std::abs(std::sinh(c.theta0[i] * th))) : 0.0; };
    auto cosh_term = [&](std::size_t i) { return (c.n0[i] != 0) ? c.n0[i] * std::log(std::abs(std::cosh(c.theta0[i] * th))) : 0.0; };
    auto sinh_dterm = [&](std::size_t i) { return (c.n0[i] != 0) ? c.n0[i] * c.theta0[i] * th / std::tanh(c.theta0[i] * th) : 0.0; };

    const double a00_0 = Rstar_R;
    const double a00_1 = Rstar_R * th;
    const double a00_2 = std::log(rho0 / rhoci) + Rstar_R * (c.n0[3] * std::log(th) + sinh_term(4) + sinh_term(6) - cosh_term(5) - cosh_term(7));

    // NOTE: teqp guards the two sinh terms with `!= 0` but deliberately does
    // NOT guard the two cosh terms (GERG.hpp:60-61).  Where n0[5] or n0[7] is
    // zero the unguarded term evaluates to zero anyway, so the behaviour is
    // identical; the asymmetry is preserved so the two sources stay diffable.
    const double a10_0 = 0.0;
    const double a10_1 = Rstar_R * th;
    const double a10_2 = Rstar_R
                         * (c.n0[3] + sinh_dterm(4) + sinh_dterm(6) - c.n0[5] * c.theta0[5] * th * std::tanh(c.theta0[5] * th)
                            - c.n0[7] * c.theta0[7] * th * std::tanh(c.theta0[7] * th));

    // Row 0: h0/(R*T0) = 1 + Aig10 = 0, so Aig10 = -1.
    // Row 1: s0/R = Aig10 - Aig00 = 0.
    const double A00 = a10_0, A01 = a10_1, b0 = -1.0 - a10_2;
    const double A10 = a10_0 - a00_0, A11 = a10_1 - a00_1, b1 = -a10_2 + a00_2;

    const double det = A00 * A11 - A01 * A10;
    if (std::abs(det) < 1e-300) {
        throw ValueError("GERG ideal-gas integration constants: singular 2x2 system");
    }
    const double n1 = (b0 * A11 - A01 * b1) / det;
    const double n2 = (A00 * b1 - b0 * A10) / det;
    return {n1, n2};
}

AlphaigCoeffs get_alphaig_coeffs(GERGModel model, const std::string& gerg_name) {
    // Throws if gerg_name is not a component of this model, and gives us the
    // reducing Tc/rhoc that the integration constants are solved against.
    const PureInfo info = get_pure_info(model, gerg_name);

    AlphaigCoeffs c;
    bool found = false;
    if (model == GERGModel::GERG_2008) {
        const auto& ov = alphaig_2008_overrides();
        auto it = ov.find(gerg_name);
        if (it != ov.end()) {
            c = pad_alphaig(gerg_name, it->second);
            found = true;
        }
    }
    if (!found) {
        const auto& base = alphaig_2004();
        auto it = base.find(gerg_name);
        if (it == base.end()) {
            // Unreachable with the tables as shipped: get_pure_info above has
            // already rejected anything outside component_names(model), and
            // the two ideal-gas tables cover every name in it.  Kept as a
            // drift guard -- if a component is ever added to component_names
            // without an ideal-gas row, this fires (the "padded monograph
            // shape" test sweeps every component of both models, so it fires
            // in the test suite rather than in user code).
            throw ValueError(format("Unable to load GERG ideal-gas coefficients for [%s]", gerg_name.c_str()));
        }
        c = pad_alphaig(gerg_name, it->second);
    }

    // Discard the published integration constants and re-solve them so that
    // h = s = 0 for the IDEAL GAS at 298.15 K and 101325 Pa.  teqp
    // GERG.hpp:370-382.  Note that rho0 uses R, not R*.
    const double T0 = 298.15;    // K
    const double p0 = 101325.0;  // Pa
    const double rho0 = p0 / (R_GERG * T0);
    const auto n12 = recalc_integration_constants(c, T0, info.Tc_K, rho0, info.rhoc_molm3, RSTAR_GERG / R_GERG);
    c.n0[1] = n12.first;
    c.n0[2] = n12.second;
    return c;
}

std::string resolve_component(GERGModel model, const std::string& user_name) {
    // A canonical GERG component name resolves to itself.  Several of them
    // ("n-butane", "carbondioxide", "hydrogensulfide", ...) are NOT spellings
    // CoolProp's own alias machinery recognises, so without this fast path
    // component_names(model) would not round-trip through the factory --
    // i.e. AbstractState::factory("GERG2008", {component_names(...)[5]})
    // would throw.
    {
        const auto& names = component_names(model);
        if (std::find(names.begin(), names.end(), user_name) != names.end()) {
            return user_name;
        }
    }

    // Otherwise resolve through CoolProp's normal alias/CAS machinery, so
    // users can spell components the way they do everywhere else in CoolProp.
    std::string cas;
    try {
        cas = get_fluid_param_string(user_name, "CAS");
    } catch (const std::exception&) {
        throw ValueError(format("[%s] is not a fluid CoolProp recognises, so it cannot be a GERG component", user_name.c_str()));
    }
    const auto& table = detail::cas_to_gerg();
    auto it = table.find(cas);
    if (it == table.end()) {
        throw ValueError(format("[%s] (CAS %s) is not a component of this GERG model", user_name.c_str(), cas.c_str()));
    }
    const std::string& gerg_name = it->second;
    const auto& names = component_names(model);
    if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
        throw ValueError(format("[%s] is a GERG-2008 component but not a GERG-2004 component", user_name.c_str()));
    }
    return gerg_name;
}

BetasGammas get_betasgammas(GERGModel model, const std::string& f1, const std::string& f2) {
    if (model == GERGModel::GERG_2008) {
        const auto& overrides = betasgammas_2008_overrides();
        auto it = overrides.find({f1, f2});
        if (it != overrides.end()) {
            return it->second;
        }
        auto rit = overrides.find({f2, f1});
        if (rit != overrides.end()) {
            BetasGammas bg = rit->second;
            bg.betaV = 1.0 / bg.betaV;
            bg.betaT = 1.0 / bg.betaT;
            return bg;
        }
    }
    // Fall through to GERG-2004 (also GERG-2004's own lookup path).
    const auto& base = betasgammas_2004();
    auto it = base.find({f1, f2});
    if (it != base.end()) {
        return it->second;
    }
    auto rit = base.find({f2, f1});
    if (rit != base.end()) {
        BetasGammas bg = rit->second;
        bg.betaV = 1.0 / bg.betaV;
        bg.betaT = 1.0 / bg.betaT;
        return bg;
    }
    throw ValueError(format("Unable to obtain GERG binary reducing parameters for the pair [%s, %s]", f1.c_str(), f2.c_str()));
}

bool get_Fij(GERGModel model, const std::string& f1, const std::string& f2, double& F) {
    // GERG-2008 reuses GERG-2004's F_ij unchanged, so model does not select
    // between two tables here (unlike get_pure_info/get_betasgammas above).
    (void)model;
    const auto& table = fij_table();
    auto it = table.find({f1, f2});
    if (it != table.end()) {
        F = it->second;
        return true;
    }
    auto rit = table.find({f2, f1});
    if (rit != table.end()) {
        F = rit->second;
        return true;
    }
    return false;
}

DepartureCoeffs get_departurecoeffs(GERGModel model, const std::string& f1, const std::string& f2) {
    // GERG-2008 reuses GERG-2004's departure coefficients unchanged.
    (void)model;
    if (is_generalized_departure_pair(f1, f2)) {
        return generalized_departure();
    }
    const auto& table = departure_specific_table();
    auto it = table.find({f1, f2});
    if (it != table.end()) {
        return it->second;
    }
    auto rit = table.find({f2, f1});
    if (rit != table.end()) {
        return rit->second;
    }
    throw ValueError(format("Unable to obtain GERG departure coefficients for the pair [%s, %s]", f1.c_str(), f2.c_str()));
}

std::size_t departure_Npower(const DepartureCoeffs& dc) {
    std::size_t np = 0;
    while (np < dc.n.size() && dc.eta[np] == 0.0 && dc.beta[np] == 0.0) {
        ++np;
    }
    // CoolProp's GERG2008DepartureFunction assumes the polynomial block is a
    // contiguous prefix.  If any later term is also polynomial, that
    // assumption is violated and the split would silently drop terms.
    for (std::size_t k = np; k < dc.n.size(); ++k) {
        if (dc.eta[k] == 0.0 && dc.beta[k] == 0.0) {
            throw ValueError("GERG departure coefficients: polynomial terms are not a contiguous prefix");
        }
    }
    return np;
}

namespace {

/// The generated ancillary row for one component, resolved through the same
/// base-plus-overrides shape get_pure_info uses.  Throws ValueError if the
/// component is not part of the given model, so a GERG-2004 caller cannot
/// silently receive a GERG-2008-only fluid's ancillary.
const detail::AncillarySet& lookup_ancillary_set(GERGModel model, const std::string& gerg_name) {
    if (model == GERGModel::GERG_2008) {
        const auto& ov = detail::ancillaries_2008_overrides();
        auto it = ov.find(gerg_name);
        if (it != ov.end()) {
            return it->second;
        }
    }
    // GERG-2004 does not contain the three fluids added in GERG-2008; reject
    // them here rather than falling through to the 2004 base table, which does
    // not have them either but would produce a less specific error.
    const auto& names = component_names(model);
    if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
        throw ValueError(format("[%s] is not a component of this GERG model", gerg_name.c_str()));
    }
    const auto& base = detail::ancillaries_2004();
    auto it = base.find(gerg_name);
    if (it == base.end()) {
        throw ValueError(format("Unable to load GERG saturation ancillaries for [%s]", gerg_name.c_str()));
    }
    return it->second;
}

}  // namespace

AncillaryCoeffs get_ancillary(GERGModel model, const std::string& gerg_name, const std::string& which) {
    const detail::AncillarySet& set = lookup_ancillary_set(model, gerg_name);
    AncillaryCoeffs anc;
    if (which == "rhoL") {
        anc = set.rhoL;
        // The strings are the ones cpjson::make_saturation_ancillary accepts
        // (FluidLibraryFactories.h:57-60): "rhoLnoexp" selects
        // TYPE_NOT_EXPONENTIAL, anything else TYPE_EXPONENTIAL.
        anc.type = "rhoLnoexp";
    } else if (which == "rhoV") {
        anc = set.rhoV;
        anc.type = "rhoV";
    } else if (which == "pV") {
        anc = set.pV;
        anc.type = "pV";
    } else {
        throw ValueError(format("GERG ancillary [%s] is not one of rhoL, rhoV, pV", which.c_str()));
    }
    return anc;
}

SatEndState get_sat_min_state(GERGModel model, const std::string& gerg_name) {
    return lookup_ancillary_set(model, gerg_name).sat_min;
}

double get_acentric_factor(GERGModel model, const std::string& gerg_name) {
    // SINGLE EXIT POINT, deliberately.  The table lookup below has two
    // successful paths -- the GERG-2008 override table and the GERG-2004 base
    // table -- and the finiteness check at the bottom is worthless to any path
    // that returns around it.  Returning `it->second` straight out of the
    // override branch would leave the guard inert for exactly the five rows
    // that live in that table (carbonmonoxide, isopentane, hydrogensulfide,
    // n-nonane, n-decane):
    // injecting an infinity into n-nonane's row built a fluid and answered
    // `inf` from acentric_factor(), while the same injection into water threw.
    // The result is assigned into `omega` here and validated once, on the way
    // out, so a new lookup path cannot silently acquire the same hole.
    double omega = 0.0;
    bool found = false;

    // Same base-plus-overrides resolution as lookup_ancillary_set /
    // get_pure_info: GERG-2008 moved carbon monoxide's and isopentane's
    // reducing parameters, which moves their saturation curve and therefore
    // their acentric factor, and the three fluids GERG-2008 adds must NOT
    // fall through to the GERG-2004 table.
    if (model == GERGModel::GERG_2008) {
        const auto& ov = detail::acentric_2008_overrides();
        auto it = ov.find(gerg_name);
        if (it != ov.end()) {
            omega = it->second;
            found = true;
        }
    }
    if (!found) {
        const auto& names = component_names(model);
        if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
            throw ValueError(format("[%s] is not a component of this GERG model", gerg_name.c_str()));
        }
        const auto& base = detail::acentric_2004();
        auto it = base.find(gerg_name);
        if (it == base.end()) {
            throw ValueError(format("Unable to load a GERG acentric factor for [%s]", gerg_name.c_str()));
        }
        omega = it->second;
    }

    // A table row that is not a finite number would propagate silently: every
    // consumer (Wilson K-factors, solver_rho_Tp_SRK, T_DP_PengRobinson) reads
    // EquationOfState::acentric directly and none of them checks.  That is the
    // exact failure mode the _HUGE sentinel caused, so it is rejected here
    // rather than re-created by a hand-edit of GERGAcentric.h.  (The generator
    // also refuses to emit a non-finite omega, but that is not what this guard
    // is for -- the generator cannot police an edit made after it ran.)
    if (!ValidNumber(omega)) {
        throw ValueError(format("The GERG acentric factor for [%s] is not a finite number", gerg_name.c_str()));
    }
    return omega;
}

double evaluate_ancillary(const AncillaryCoeffs& anc, double T) {
    // Deliberately mirrors SaturationAncillaryFunction::evaluate
    // (Ancillaries.cpp:43-79) term for term, INCLUDING the quiet NaN for
    // T > T_r: pow(negative, fractional) is otherwise a NaN-or-SIGFPE
    // (GitHub #1611), and internal callers rely on the NaN.  If this drifts
    // from CoolProp's evaluator, the tests that compare a fitted ancillary to
    // a converged saturation state stop testing the code that actually runs.
    const double theta = 1.0 - T / anc.T_r;
    if (theta < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double summer = 0;
    for (std::size_t i = 0; i < anc.n.size(); ++i) {
        summer += anc.n[i] * std::pow(theta, anc.t[i]);
    }
    if (anc.type == "rhoLnoexp") {
        return anc.reducing_value * (1 + summer);
    }
    // "rhoV" and "pV" both carry the T_r/T factor (using_tau_r = true).
    return anc.reducing_value * std::exp(anc.T_r / T * summer);
}

namespace {

/// Build a CoolProp SaturationAncillaryFunction from a GERG-fitted row.
///
/// This goes through SaturationAncillaryFunction::Values -- the same
/// plain-typed bundle cpjson::make_saturation_ancillary fills from JSON
/// (FluidLibraryFactories.h:41-73) -- so the ancillary that ends up on a GERG
/// fluid is byte-for-byte the same kind of object every other CoolProp fluid
/// carries, and SaturationAncillaryFunction::evaluate is the single evaluator.
/// Reimplementing the evaluation here instead would leave the flash routines
/// reading one function and the tests another.
SaturationAncillaryFunction make_anc(const AncillaryCoeffs& anc) {
    SaturationAncillaryFunction::Values v;
    v.type = (anc.type == "rhoLnoexp") ? SaturationAncillaryFunction::TYPE_NOT_EXPONENTIAL : SaturationAncillaryFunction::TYPE_EXPONENTIAL;
    v.using_tau_r = (anc.type != "rhoLnoexp");
    v.n = anc.n;
    v.t = anc.t;
    v.Tmin = anc.Tmin;
    v.Tmax = anc.Tmax;
    v.reducing_value = anc.reducing_value;
    v.T_r = anc.T_r;
    return SaturationAncillaryFunction(v);
}

}  // namespace

CoolPropFluid make_gerg_fluid(GERGModel model, const std::string& gerg_name) {
    const PureInfo info = get_pure_info(model, gerg_name);
    const PureCoeffs pc = get_pure_coeffs(model, gerg_name);
    // get_alphaig_coeffs re-solves a 2x2 system on every call.  It is called
    // exactly once per component here, at construction; nothing in a flash
    // loop may call it.
    const AlphaigCoeffs ac = get_alphaig_coeffs(model, gerg_name);

    CoolPropFluid fluid;
    fluid.name = gerg_name;
    fluid.EOSVector.resize(1);
    EquationOfState& EOS = fluid.EOS();

    EOS.R_u = R_GERG;  // 8.314472, NOT CoolProp's per-fluid R and NOT R*
    EOS.molar_mass = info.M_kgmol;
    EOS.pseudo_pure = false;
    // EquationOfState has no default constructor, so every scalar it owns has
    // to be written here or it stays indeterminate.  GERG publishes no triple
    // point, so Ttriple is a sentinel: Ttriple = 0 keeps the
    // (T > Tc && T > Ttriple) supercritical branch of
    // FlashRoutines::DHSU_T_flash reachable at every temperature the EOS is
    // evaluated at, including helium's 3.6 K.
    //
    // GERG publishes no acentric factor either, but unlike the triple point
    // that one cannot be left at a sentinel: three separate consumers read
    // EquationOfState::acentric DIRECTLY, without going through the
    // overridable calc_acentric_factor accessor, and none of them checks it --
    // SaturationSolvers' Wilson K-factor seed (so every mixture QT/PQ/DQ flash
    // and build_phase_envelope), FlashRoutines::T_DP_PengRobinson (every
    // gas-phase DmolarP flash), and solver_rho_Tp_SRK, which is what
    // solver_rho_Tp seeds from (the pure-fluid HmolarP/PSmolar/PUmolar
    // flashes).  With _HUGE in place they produced NaN guesses and failed
    // several solver frames later with unrelated-looking messages.  So it is
    // DERIVED from GERG's own equation instead of borrowed from CoolProp's
    // per-fluid library, which belongs to a different equation of state:
    // omega = -1 - log10(p_sat(0.7 T_c)/p_c), evaluated offline against this
    // backend's own coefficient tables.  See GERGAcentric.h and
    // dev/gerg/compute_acentric.py.
    EOS.acentric = get_acentric_factor(model, gerg_name);
    EOS.Ttriple = 0.0;
    EOS.ptriple = _HUGE;

    // Reducing state == the GERG reducing state of Table A3.5.  For most
    // components this is the critical point; for a few (methane, water, ...)
    // the tabulated rho_c is a fitted reducing density rather than the true
    // critical density.  Both `reduce` and `crit` are set from it because
    // calc_alpha0_deriv_nocache reads iT_reducing on the pure branch and
    // iT_critical on the mixture branch, and they must agree.
    EOS.reduce.T = info.Tc_K;
    EOS.reduce.rhomolar = info.rhoc_molm3;
    // fluid.crit is assigned from EOS.reduce AFTER the residual terms exist,
    // because reduce.p is evaluated from them (below).

    // GERG-2008 extended range of validity: 60-700 K, p <= 70 MPa
    // (Kunz & Wagner 2012, section 4.1).  Applied verbatim, with no
    // per-component adjustment: GERG publishes no per-component lower
    // temperature limit, and CoolProp's own triple-point data is deliberately
    // not consulted here (a strict GERG backend must not borrow it).
    //
    // For helium (Tc = 5.1953 K) and hydrogen (Tc = 33.19 K) the whole valid
    // range therefore lies ABOVE the critical temperature.  That is not a
    // contradiction, it is the physics: at 60 K and above these two are
    // supercritical, so as PURE fluids they have no subcritical region here
    // and no saturation state can be requested.  The range check refuses such
    // a request with the range in the message, which is a better answer than
    // inventing a lower limit the model does not publish.  Both are ordinary
    // components in a MIXTURE, where the mixture reducing temperature governs.
    EOS.limits.Tmin = 60.0;
    EOS.limits.Tmax = 700.0;
    EOS.limits.pmax = 70e6;
    EOS.limits.rhomax = 1e6;

    // --- Residual part -------------------------------------------------
    // alpha^r = sum_i n_i delta^{d_i} tau^{t_i} exp(-c_i delta^{l_i}), with
    // c_i == 1 exactly when l_i > 0.  add_Power derives c from l that same
    // way, so PureCoeffs::c is not passed (it exists only for diffing against
    // teqp).
    {
        const std::vector<CoolPropDbl> n(pc.n.begin(), pc.n.end());
        const std::vector<CoolPropDbl> d(pc.d.begin(), pc.d.end());
        const std::vector<CoolPropDbl> t(pc.t.begin(), pc.t.end());
        const std::vector<CoolPropDbl> l(pc.l.begin(), pc.l.end());
        EOS.alphar.GenExp.add_Power(n, d, t, l);
        EOS.alphar.GenExp.finish();
    }

    // Pressure at the reducing state, from the EOS that was just assembled:
    // p = rho*R*T*(1 + delta*dalphar/ddelta) evaluated at tau = delta = 1.
    // Without this, SimpleState's _HUGE default propagates and
    // AbstractState::p_critical() returns inf (calc_p_critical returns
    // components[0].crit.p for a pure fluid).  Note this is the pressure at
    // the GERG REDUCING point, which for a handful of components is a fitted
    // point rather than the measured critical point -- it is the only
    // pressure the GERG tables can produce without importing outside data.
    EOS.reduce.p = info.rhoc_molm3 * R_GERG * info.Tc_K * (1.0 + EOS.dalphar_dDelta(1.0, 1.0));

    // get_fluid_constant maps iT_critical/irhomolar_critical to fluid.crit and
    // iT_reducing/irhomolar_reducing to EOS().reduce; calc_alpha0_deriv_nocache
    // reads the reducing pair on the pure branch and the critical pair on the
    // mixture branch, so the two MUST agree.  AbstractState::T_critical() /
    // rhomolar_critical() / p_critical() also read fluid.crit, which is what
    // the "GERG pure fluid reports the GERG reducing point as its critical
    // point" test pins.
    fluid.crit = EOS.reduce;

    // --- Ideal-gas part ------------------------------------------------
    //
    // The expression the stored table is defined against (GERGData.h):
    //
    //   alpha^o_i = ln(rho/rho_c,i)
    //             + (R*/R) * [ n0[1] + n0[2]*(Tc,i/T) + n0[3]*ln(Tc,i/T)
    //                        + n0[4]*ln|sinh(theta0[4]*Tc,i/T)|
    //                        + n0[6]*ln|sinh(theta0[6]*Tc,i/T)|
    //                        - n0[5]*ln|cosh(theta0[5]*Tc,i/T)|
    //                        - n0[7]*ln|cosh(theta0[7]*Tc,i/T)| ]
    //
    // Four things have to line up, and three of them are invisible in p, c_v
    // and w -- they show up only in h and s (equivalently, in alphaig).
    //
    // 1. SIGN.  IdealHelmholtzGERG2004Sinh accumulates +n[i]*log|sinh(...)|
    //    and IdealHelmholtzGERG2004Cosh accumulates -n[i]*log|cosh(...)|
    //    (src/Helmholtz.cpp:1260-1341): the minus in front of the two cosh
    //    terms above is applied by the TERM, not carried in the coefficient.
    //    So all four of n0[4..7] go in exactly as stored, un-negated.  Some
    //    published n0 are themselves negative (methane's n0[7] is -4.46921);
    //    "as published" is the rule, not "positive".
    //
    // 2. R*/R.  CoolProp's terms know nothing about it, so it is applied here
    //    to the WHOLE bracket, i.e. to all seven of n0[1..7].  It must NOT
    //    touch the ln(rho/rho_c) term -- and it does not, because that term is
    //    generated by IdealHelmholtzLead as an unscaled log(delta) alongside
    //    a1 + a2*tau (Helmholtz.cpp:1083-1093).  n0[1] and n0[2] as returned
    //    by get_alphaig_coeffs were re-solved on exactly this basis, so they
    //    stay consistent under the uniform scaling.
    //
    // 3. Tc-vs-T_red.  NO Tc/T_red folding is applied to n0[2]/n0[3] here,
    //    because CoolProp already does it, and doing it again would
    //    double-apply.  Read HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache:
    //
    //      * pure branch (HelmholtzEOSMixtureBackend.cpp:3586-3600) evaluates
    //        the terms at taustar = Tc/Tr * tau == Tc/T, having first called
    //        E.alpha0.set_Tred(Tc) with the COMPONENT's Tc;
    //      * mixture branch (:3636-3644) evaluates component i's terms at
    //        tau_i = T_ci * tau / Tr == Tc,i/T, with T_ci read from
    //        iT_critical.
    //
    //    In BOTH branches the argument handed to alpha0 is the component's own
    //    Tc,i/T -- never Tr/T.  So IdealHelmholtzLead(a1, a2) already yields
    //    n0[1] + n0[2]*Tc,i/T and IdealHelmholtzLogTau(a1) already yields
    //    n0[3]*ln(Tc,i/T), for pure fluids AND for mixtures, with the raw
    //    coefficients.  Folding Tc/T_red into n0[2] and adding
    //    n0[3]*ln(Tc/T_red) would apply the rescale a second time.
    //
    //    The residual risk moved to the sinh/cosh terms, and it was a MIXTURE
    //    risk only.  CoolProp's IdealHelmholtzGERG2004Sinh/Cosh compute
    //    t = Tc_ctor/T_red internally and evaluate at t*tau_arg
    //    (src/Helmholtz.cpp:1260-1341).  On the pure branch the base class
    //    calls set_Tred(Tc) with the component's own Tc, so t == 1 and the
    //    argument is theta*Tc/T -- correct.  On the MIXTURE branch it calls
    //    set_Tred(Tr) with the mixture reducing temperature while still
    //    passing tau_i = Tc,i/T, leaving a spurious extra factor Tc,i/Tr
    //    inside every sinh/cosh.  No choice of Tc_ctor fixes that, because the
    //    value it would need is Tr, which is composition-dependent and not
    //    known at construction.
    //
    // 4. WHY THIS BACKEND DOES NOT USE GERG2004Sinh/GERG2004Cosh AT ALL.
    //    Rather than route around the above with a virtual override of
    //    calc_alpha0_deriv_nocache / calc_all_alpha0_derivs_nocache (neither
    //    of which is virtual in this tree -- HelmholtzEOSMixtureBackend.h:641
    //    and :644 -- and which in any case is behaviour, so it would not
    //    survive the copy into a linked saturation state), the two terms are
    //    re-expressed exactly, in a form that has no T_red dependence to get
    //    wrong.  With x = theta*tau and tau the argument CoolProp actually
    //    passes (Tc,i/T on BOTH branches, per point 3):
    //
    //      ln|sinh x| =  x + ln(1 - e^{-2x}) - ln 2
    //      ln|cosh x| =  x + ln(1 + e^{-2x}) - ln 2      (x >= 0)
    //
    //    so  +n*ln|sinh(theta*tau)| becomes
    //          a2 += n*theta ; a1 -= n*ln2 ; PlanckEinsteinGeneralized term
    //          (n, theta_term = -2*theta, c = 1, d = -1)
    //    and -n*ln|cosh(theta*tau)| becomes
    //          a2 -= n*theta ; a1 += n*ln2 ; PlanckEinsteinGeneralized term
    //          (-n, theta_term = -2*theta, c = 1, d = +1),
    //    where IdealHelmholtzPlanckEinsteinGeneralized contributes
    //    n*log(c + d*exp(theta_term*tau)) (Helmholtz.cpp:1152-1177).
    //
    //    None of IdealHelmholtzLead / LogTau / PlanckEinsteinGeneralized reads
    //    _Tr, so set_Tred becomes irrelevant and the pure and mixture branches
    //    agree by construction -- including in SatL/SatV, which carry data
    //    rather than behaviour.  The identity is exact, so the 3070 pure-fluid
    //    reference comparisons (which pass at 1e-12 on alphaig) are a direct
    //    check on the rewrite; the mixture alphaig column is what checks that
    //    the mixture branch now agrees too.  As a side benefit e^{-2x}
    //    underflows to 0 for large x where sinh(x) would overflow.
    //
    //    cosh is EVEN and two published theta0 are negative (nitrogen's
    //    theta0[5] = -5.3931, carbon dioxide's = -2.8444), so |theta| is used:
    //    ln|cosh(theta*tau)| == ln|cosh(|theta|*tau)|.  Every theta0 on a SINH
    //    index (4, 6) with a nonzero n0 is positive as published, so the same
    //    absolute value is a no-op there; it is applied to both for symmetry.
    {
        const double RR = RSTAR_GERG / R_GERG;
        const double LN2 = std::log(2.0);

        double a1 = RR * ac.n0[1];
        double a2 = RR * ac.n0[2];
        std::vector<CoolPropDbl> pe_n, pe_theta, pe_c, pe_d;

        // AlphaigCoeffs vectors are length 8 and 1-based: n0[0] and
        // theta0[0..3] are padding.  Never iterate from 0.
        //
        // Guard on n0, exactly as teqp's recalc_integration_constants does:
        // helium and argon have all four of n0[4..7] AND all four theta0
        // zero, and an unguarded n*log|sinh(0*tau)| would be 0*(-inf) == NaN
        // rather than 0.  A nonzero n0 paired with a zero theta0 would be a
        // genuine singularity in the published expression (ln|sinh 0|), not
        // something the rewrite introduces; no shipped component has one, and
        // the throw below is a drift guard for a future table edit.
        auto emit = [&](std::size_t k, bool is_sinh) {
            if (ac.n0[k] == 0.0) {
                return;
            }
            const double th = std::abs(ac.theta0[k]);
            if (th == 0.0) {
                throw ValueError(format("GERG ideal-gas coefficients for [%s]: n0[%d] is nonzero but theta0[%d] is zero", gerg_name.c_str(),
                                        static_cast<int>(k), static_cast<int>(k)));
            }
            const double n = (is_sinh ? 1.0 : -1.0) * RR * ac.n0[k];
            a2 += n * th;
            a1 -= n * LN2;
            pe_n.push_back(n);
            pe_theta.push_back(-2.0 * th);
            pe_c.push_back(1.0);
            pe_d.push_back(is_sinh ? -1.0 : 1.0);
        };
        emit(4, true);
        emit(6, true);
        emit(5, false);
        emit(7, false);

        EOS.alpha0.Lead = IdealHelmholtzLead(a1, a2);
        EOS.alpha0.LogTau = IdealHelmholtzLogTau(RR * ac.n0[3]);
        if (!pe_n.empty()) {
            EOS.alpha0.PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(pe_n, pe_theta, pe_c, pe_d);
        }
    }

    // --- Saturation ancillaries and the saturation-curve end states -----
    //
    // The ancillaries are FITTED AGAINST THIS EOS (dev/gerg/fit_ancillaries.py,
    // table in GERGAncillaries.h), not borrowed from CoolProp's library.
    // CoolProp's shipped ancillaries for these same fluids belong to each
    // fluid's REFERENCE equation of state -- Setzmann-Wagner for methane,
    // IAPWS-95 for water -- and would work perfectly well as seeds, which is
    // exactly the trap: a backend named GERG-2008 would then carry data
    // traceable to a different equation with nothing anywhere to say so.
    //
    // A SUPERANCILLARY IS A DIFFERENT THING AND MUST NEVER BE ATTACHED.  An
    // ancillary is a seed that saturation_T_pure_Maxwell iterates away from;
    // FlashRoutines::sat_superanc_path_applies (FlashRoutines.cpp:558) takes
    // the Chebyshev superancillary and RETURNS IT as the answer.  Nothing
    // below calls EquationOfState::set_superancillaries_str, so get_superanc()
    // stays null and that path is unreachable -- pinned by the test
    // "GERG fluids carry no superancillary".
    {
        fluid.ancillaries.rhoL = make_anc(get_ancillary(model, gerg_name, "rhoL"));
        fluid.ancillaries.rhoV = make_anc(get_ancillary(model, gerg_name, "rhoV"));
        // CoolProp keeps separate saturated-liquid and saturated-vapour
        // pressure ancillaries so that pseudo-pure fluids can have a glide.
        // GERG's components are all true pure fluids, so both slots get the
        // same curve -- the same thing FluidLibrary.h:1136-1137 does when a
        // fluid's JSON has one "pS" block instead of a "pL"/"pV" pair.
        const AncillaryCoeffs pV = get_ancillary(model, gerg_name, "pV");
        fluid.ancillaries.pL = make_anc(pV);
        fluid.ancillaries.pV = make_anc(pV);

        // The saturation state at the low-temperature end of the fitted range,
        // traced with teqp rather than guessed.  calc_Tmin_sat/calc_pmin_sat
        // return these, and QT_flash uses them as its lower bound
        // (FlashRoutines.cpp:920-934), so a wrong value here shows up as
        // either a spurious out-of-range throw or a saturation call attempted
        // where there is no data.
        const SatEndState end = get_sat_min_state(model, gerg_name);
        EOS.sat_min_liquid.T = end.T_K;
        EOS.sat_min_liquid.p = end.p_Pa;
        EOS.sat_min_liquid.rhomolar = end.rhoL_molm3;
        EOS.sat_min_vapor.T = end.T_K;
        EOS.sat_min_vapor.p = end.p_Pa;
        EOS.sat_min_vapor.rhomolar = end.rhoV_molm3;

        // triple_liquid / triple_vapor get the SAME state, and the name is
        // CoolProp's, not a claim about GERG.  GERG publishes no triple point
        // and EOS.Ttriple stays 0 above; saturation_T_pure_Maxwell
        // (VLERoutines.cpp:965-980) reads these two slots purely as "the
        // low-temperature end of the saturation curve", to sanity-band the
        // ancillary seed (rhoL > 1.2*tripleL.rhomolar rejects it) and to build
        // a linear fallback seed through (tripleL.T, tripleL.rhomolar).  With
        // both left at the SimpleState default of _HUGE, every seed would fail
        // that band and every pure-fluid saturation call would take the
        // fallback path.
        fluid.triple_liquid.T = end.T_K;
        fluid.triple_liquid.p = end.p_Pa;
        fluid.triple_liquid.rhomolar = end.rhoL_molm3;
        fluid.triple_vapor.T = end.T_K;
        fluid.triple_vapor.p = end.p_Pa;
        fluid.triple_vapor.rhomolar = end.rhoV_molm3;
    }

    // NO hs_anchor.  CoolProp's convention is a single-phase point at
    // (1.1*Tc, 0.9*rhoc) whose h and s serve as the reference offset for
    // h/s-input flashes.  It is dead weight for this backend: every site that
    // reads hs_anchor.hmolar/smolar (VLERoutines.cpp:200-294) pairs it with
    // components[0].ancillaries.hL or .sL, and GERG populates only rhoL, rhoV,
    // pL and pV -- so those paths fail on the empty h/s ancillary whatever
    // hs_anchor holds.  They are the saturation_PHSU_pure family, which is
    // unsupported for pure GERG fluids anyway (bd CoolProp-kr4o).
    //
    // Leaving it unset also removes the only reason this function had to
    // evaluate at 1.1*Tc, which for helium and hydrogen is below the published
    // 60 K lower limit.  The repo-wide hs_anchor machinery is untouched;
    // FluidLibrary still fills it for every JSON fluid.

    // NOTE: EquationOfState::validate() (CoolPropFluid.h:450-453) is two bare
    // assert()s on R_u and molar_mass, so it is compiled out entirely under
    // NDEBUG -- i.e. in the mandated Release build it does nothing.  It is
    // called because it is the conventional end of a fluid build, NOT because
    // it guards anything here.  The real coverage for R_u and molar_mass is
    // the reference `w` column (w = sqrt(-R T/M * ...)) and the
    // "uses the GERG gas constant and reducing state" test.
    EOS.validate();

    // reduce.hmolar/smolar, and the same two on triple_liquid/triple_vapor,
    // can only be filled by EVALUATING the assembled EOS, which needs a
    // backend.  Nothing in CoolProp READS these four; they surface only
    // through get_state("reducing") / get_state("critical") /
    // get_state("triple_liquid"), and a public accessor returning _HUGE is a
    // trap for the next reader, so they are filled here.
    //
    // The temporary is GERG-typed, with generate_SatL_and_SatV = false so the
    // recursion stops: a base HelmholtzEOSMixtureBackend would use CODATA R
    // rather than R_GERG for anything but a pure fluid, and h/s offsets wrong
    // in the 7th digit are exactly the kind of error that never shows up in p
    // or w.
    //
    // update_DmolarT_direct, NOT update_states(): update_states() evaluates at
    // the reducing state through update(), and for helium and hydrogen the
    // reducing temperature is below the published 60 K lower limit, so
    // check_gerg_range_of_validity would refuse to let the fluid be built at
    // all.  The bypass is unambiguously right here for the same reason it is
    // right just below -- the point being evaluated is one this backend
    // computed itself and stored.
    {
        GERGMixtureBackend probe(model, std::vector<CoolPropFluid>(1, fluid), false);
        probe.set_mole_fractions(std::vector<CoolPropDbl>(1, 1.0));

        // specify_phase before every update_DmolarT_direct.  The direct
        // updater skips phase determination by design, so calc_hmolar throws
        // "phase is invalid" unless the phase is imposed first.  This used to
        // be masked: update_states() ran an ordinary update() beforehand and
        // left _phase populated, which the triple-state loop below then
        // silently inherited.  Removing update_states() exposed that; imposing
        // the phase explicitly is what it should always have done.  The phase
        // is not used to select a branch of the EOS here -- alpha0 and alphar
        // are evaluated at the given (rho, T) regardless -- it only has to be
        // a valid single-phase value.
        //
        // update_DmolarT_direct, NOT update: for helium and hydrogen the
        // reducing temperature is below the published 60 K lower limit, and
        // for the light fluids the traced end of the saturation curve is too
        // (methane's is 57.17 K), so update() would throw OutOfRangeError from
        // check_gerg_range_of_validity and the fluid could not be built at
        // all.  The bypass is unambiguously right here: every point evaluated
        // is one this backend computed itself and stored.
        probe.specify_phase(iphase_supercritical);
        probe.update_DmolarT_direct(EOS.reduce.rhomolar, EOS.reduce.T);
        EOS.reduce.hmolar = probe.hmolar();
        EOS.reduce.smolar = probe.smolar();

        for (const auto& item : {std::make_pair(&fluid.triple_liquid, iphase_liquid), std::make_pair(&fluid.triple_vapor, iphase_gas)}) {
            SimpleState* st = item.first;
            probe.specify_phase(item.second);
            probe.update_DmolarT_direct(st->rhomolar, st->T);
            st->hmolar = probe.hmolar();
            st->smolar = probe.smolar();
        }
        probe.unspecify_phase();
    }
    // fluid.crit was copied from EOS.reduce before h/s existed, so refresh the
    // two fields just filled in above.  Leaving crit.hmolar at
    // _HUGE while reduce.hmolar is finite would make get_state("critical") and
    // get_state("reducing") disagree about the same point.
    fluid.crit.hmolar = EOS.reduce.hmolar;
    fluid.crit.smolar = EOS.reduce.smolar;

    return fluid;
}

}  // namespace GERG

class GERG2004Generator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) override {
        return new GERGMixtureBackend(GERGModel::GERG_2004, fluid_names);
    }
};
// This static initialization will cause the generator to register
// NOLINTNEXTLINE(cert-err58-cpp)
static GeneratorInitializer<GERG2004Generator> gerg2004_gen(GERG2004_BACKEND_FAMILY);

class GERG2008Generator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) override {
        return new GERGMixtureBackend(GERGModel::GERG_2008, fluid_names);
    }
};
// This static initialization will cause the generator to register
// NOLINTNEXTLINE(cert-err58-cpp)
static GeneratorInitializer<GERG2008Generator> gerg2008_gen(GERG2008_BACKEND_FAMILY);

} /* namespace CoolProp */
