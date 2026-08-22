#ifndef GERGBACKEND_H_
#define GERGBACKEND_H_

#include <string>
#include <vector>

#include "../Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/CoolPropFluid.h"
#include "CoolProp/DataStructures.h"
#include "GERGData.h"

namespace CoolProp {

/// Which GERG model year this backend instance represents.  Defined in
/// GERGData.h (CoolProp::GERG::GERGModel); aliased here so backend code can
/// keep referring to it as CoolProp::GERGModel.
using GERGModel = GERG::GERGModel;

namespace GERG {

/// Assemble a CoolPropFluid carrying ONLY the published GERG parameters for
/// one component: the residual Helmholtz terms of Table A3.2, the ideal-gas
/// terms of Table A3.1, the reducing state of Table A3.5, and R = 8.314472
/// J/mol/K.  Nothing is taken from CoolProp's own fluid library -- the name is
/// used purely as a table key -- so a GERG backend can never silently mix a
/// CoolProp EOS into a GERG calculation.
///
/// The ideal-gas wiring is the delicate part; see the long comment on the
/// definition in GERGBackend.cpp for the sign, R*/R and Tc-vs-T_red contract.
///
/// @param model      Which GERG model's tables to read
/// @param gerg_name  A GERG component name (as returned by resolve_component)
CoolPropFluid make_gerg_fluid(GERGModel model, const std::string& gerg_name);

}  // namespace GERG

/**
 * \brief Strict GERG-2004 / GERG-2008 backend.
 *
 * Admits only the components published with the selected model (18 for
 * GERG-2004, 21 for GERG-2008) and uses only parameters published with it.
 * See docs/superpowers/specs/2026-07-25-gerg-strict-backend-design.md.
 */
class GERGMixtureBackend : public HelmholtzEOSMixtureBackend
{
   public:
    GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names);
    ~GERGMixtureBackend() override = default;

    std::string backend_name() override {
        return get_backend_string(m_model == GERGModel::GERG_2004 ? GERG2004_BACKEND : GERG2008_BACKEND);
    }

    GERGModel model() const {
        return m_model;
    }

    /// Copy this state.  Overridden so that TPD_state / critical_state /
    /// transient_pure_state (HelmholtzEOSMixtureBackend.h:84, 92, 100, which
    /// call get_copy UNqualified) are GERG-typed and therefore reach the
    /// set_mixture_parameters override below.  A plain
    /// HelmholtzEOSMixtureBackend copy of a GERG mixture would run the default
    /// binary-pair lookup instead.
    HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true) override;

    /// Enforce make_gerg_fluid's EOS.limits.Tmin/Tmax (GERGBackend.cpp) on
    /// every update, regardless of input pair.  (Deliberately T only, not
    /// pmax too -- see check_gerg_range_of_validity's definition for why a
    /// pressure check does not belong here.)
    ///
    /// This is NOT redundant with CoolProp's ordinary flash machinery: for a
    /// pure fluid at PT_INPUTS with T above Tmax, FlashRoutines::PT_flash
    /// (FlashRoutines.cpp:298) determines the phase from p
    /// (T_phase_determination_pure_or_pseudopure), solves rho_Tp, and
    /// RETURNS -- nothing in that path compares the resulting T to
    /// EOS.limits.Tmax.  (Historical note: when this override was written, T
    /// below Tmin happened to throw as an ACCIDENT -- solver_rho_Tp's liquid
    /// branch evaluated components[0].ancillaries.rhoL on an ancillary GERG
    /// fluids did not populate.  Task 11 populates rhoL/rhoV/pL/pV, so that
    /// accidental throw is GONE; the override is now the only thing enforcing
    /// EITHER bound, not just Tmax.)  Verified empirically while writing
    /// this override, `update(PT_INPUTS, 1e5, 900.0)` on GERG2008 methane
    /// (Tmax = 700 K) returns a state instead of throwing.
    ///
    /// `update` is the one AbstractState entry point every input pair
    /// (PT/DmolarT/HmolarP/...) funnels through, so checking _T here after
    /// the base class has finished the flash catches all of them, including
    /// pairs where T is SOLVED FOR rather than given directly.  Honors
    /// DONT_CHECK_PROPERTY_LIMITS, matching every other property-limit guard
    /// in HelmholtzEOSMixtureBackend.cpp (melting-line, Tmax_sat, ...), so
    /// `update` is overridden rather than post_update: the latter is not
    /// virtual (HelmholtzEOSMixtureBackend.h:67), so
    /// HelmholtzEOSMixtureBackend::update's unqualified `post_update()` call
    /// binds to the base implementation regardless of the dynamic type and a
    /// post_update override here would never run.
    void update(CoolProp::input_pairs input_pair, double value1, double value2) override;

    /// Same guard as `update` above, on the SEPARATE virtual entry point a
    /// caller reaches by supplying an explicit `GuessesStructure`
    /// (AbstractState.h:885, HelmholtzEOSMixtureBackend.h:397).  Task 10's
    /// first review round left this open: it is virtual on `AbstractState`
    /// itself, so it is reachable with NO downcast from the plain
    /// `AbstractState*` the factory returns -- confirmed to accept T = 900 K
    /// before this override existed.  Three lines, identical shape to
    /// `update`: run the inherited flash, then apply the same range check.
    void update_with_guesses(CoolProp::input_pairs input_pair, double value1, double value2, const GuessesStructure& guesses) override {
        HelmholtzEOSMixtureBackend::update_with_guesses(input_pair, value1, value2, guesses);
        check_gerg_range_of_validity();
    }

    /// GERG-2004/GERG-2008 publish no transport correlations at all -- there is
    /// nothing "GERG" a viscosity number could be checked against, unlike the
    /// residual/ideal-gas/reducing terms above.  HelmholtzEOSMixtureBackend's
    /// inherited calc_viscosity()/calc_conductivity() would instead run
    /// CoolProp's own transport models (fit against different, possibly
    /// non-GERG, reference EOS) and return a plausible-looking number
    /// labelled GERG.  AbstractState::viscosity()/conductivity()/
    /// surface_tension() (AbstractState.cpp:793-813) cache the result but
    /// only ever populate that cache by calling calc_viscosity() etc, so
    /// overriding the calc_* hook is sufficient -- there is no separate
    /// caching path that could return a value without consulting it.
    ///
    /// Signatures copied verbatim from AbstractState.h:221,225,229 /
    /// HelmholtzEOSMixtureBackend.h:564,565,569 (return type, argument list,
    /// and the ABSENCE of a const qualifier all have to match or `override`
    /// silently declares an unrelated non-virtual function instead of
    /// overriding).
    CoolPropDbl calc_viscosity() override {
        throw NotImplementedError(
          "Transport properties are not part of the GERG-2004/GERG-2008 models. Use the HEOS backend if you want CoolProp's transport correlations.");
    }
    CoolPropDbl calc_conductivity() override {
        throw NotImplementedError(
          "Transport properties are not part of the GERG-2004/GERG-2008 models. Use the HEOS backend if you want CoolProp's transport correlations.");
    }
    CoolPropDbl calc_surface_tension() override {
        throw NotImplementedError("Surface tension is not part of the GERG-2004/GERG-2008 models.");
    }

    // NOTE: there is deliberately NO calc_acentric_factor override here any
    // more.  It used to throw NotImplementedError, on the grounds that GERG
    // publishes no acentric factor and the inherited accessor would hand back
    // make_gerg_fluid's _HUGE sentinel -- +inf presented as an answer.  That
    // reasoning ended when the sentinel did: make_gerg_fluid now writes a
    // value DERIVED from GERG's own equation (omega = -1 - log10(p_sat(0.7
    // T_c)/p_c), GERGAcentric.h), so the inherited
    // HelmholtzEOSMixtureBackend::calc_acentric_factor returns a real,
    // GERG-traceable number for a pure fluid and CoolProp's own
    // "acentric factor cannot be calculated for mixtures" ValueError for a
    // mixture.  Keeping the throw would now be the strictness violation: it
    // would hide a quantity this backend genuinely has.  The override was
    // REMOVED rather than narrowed because there is no longer any case it
    // would be right in -- the mixture case is already handled, correctly and
    // identically, by the base class.

    /// `AbstractState::change_EOS(i, name)` (AbstractState.h:1547) is a
    /// non-virtual wrapper over this virtual hook (AbstractState.h:719) --
    /// so it is reachable directly on a plain `AbstractState*`, no downcast
    /// needed, MORE exposed than the `Reducing` bypass documented below.  It
    /// violates model purity more severely than a BIP tweak too: the
    /// inherited HelmholtzEOSMixtureBackend::calc_change_EOS
    /// (HelmholtzEOSMixtureBackend.cpp:483) calls
    /// `EOS.alphar.empty_the_EOS()` and installs an entirely different
    /// residual model (SRK / Peng-Robinson / Xiang-Deiters cubic) in place of
    /// the GERG residual term -- there would be nothing "GERG" left in the
    /// EOS at all.
    ///
    /// At task-10 review this failed on its own, by accident.  The cubic
    /// path reads `EOS.reduce.T`, `EOS.reduce.p`, `EOS.reduce.rhomolar` and
    /// `EOS.acentric` straight off the stored component
    /// (HelmholtzEOSMixtureBackend.cpp:491-495), and `reduce.p` was not
    /// populated for a GERG fluid then, so the cubic got a NaN pc and the
    /// very next `update()` threw "p is not a valid number".
    ///
    /// THAT ACCIDENT IS GONE.  `make_gerg_fluid` now fills all four fields
    /// the cubic reads (GERGBackend.cpp:1427 `acentric`, :1437-1438
    /// `reduce.T`/`reduce.rhomolar`, :1481 `reduce.p`), so the swap would now
    /// succeed silently and leave an SRK or Peng-Robinson cubic sitting
    /// behind a backend still called GERG2004/GERG2008.  This throw is the
    /// only thing closing that -- which is precisely why task 10 declined to
    /// rely on the accident while it still held (compare the T-below-Tmin
    /// case, see check_gerg_range_of_validity).
    void calc_change_EOS(const std::size_t, const std::string&) override {
        throw ValueError("GERG fluids' equation of state is fixed by the published model and cannot be replaced. "
                         "Use the HEOS or cubic backends directly if you want a different EOS.");
    }

   public:
    /// GERG's betas/gammas and departure functions are the published model,
    /// not a fit that a caller is meant to adjust -- set_mixture_parameters
    /// (declared below) exists specifically so this backend does NOT answer from
    /// CoolProp's mutable global BIP library.  Overriding
    /// set_mixture_parameters alone does not close this: AbstractState
    /// declares FOUR public mutator routes to the same underlying data --
    /// index-keyed and CAS-keyed double, index-keyed and CAS-keyed string
    /// (AbstractState.h:950-965) -- and HelmholtzEOSMixtureBackend overrides
    /// only the two INDEX-keyed ones (HelmholtzEOSMixtureBackend.h:217,223).
    /// The two CAS-keyed overloads are therefore inherited straight from
    /// AbstractState's default, which already throws NotImplementedError
    /// unconditionally ("GERG refuses binary-interaction mutation by every
    /// public route" pins that this stays true without a GERG-specific
    /// override of its own -- if HelmholtzEOSMixtureBackend ever grows a
    /// CAS-keyed override, that test starts failing instead of silently
    /// reopening the hole).  Only the two index-keyed overloads need a
    /// GERG-specific override here.
    ///
    /// apply_simple_mixing_rule (HelmholtzEOSMixtureBackend.cpp:383) is a
    /// FIFTH nominal route, but it is not overridden separately: it calls
    /// set_binary_interaction_double(i, j, ...) UNQUALIFIED on `this`, so
    /// virtual dispatch already sends that call to the override below
    /// whenever `this` is a GERGMixtureBackend.
    ///
    /// ==================================================================
    /// KNOWN BYPASSES.  Full enumeration, task-10 review round 2.  Two tiers:
    /// reachable from a plain `AbstractState*` with NO downcast (the shape
    /// the factory hands back and the only one a strictness claim about
    /// "GERG cannot be mutated" can honestly ignore), and reachable only by
    /// downcasting to `HelmholtzEOSMixtureBackend*` or reaching through a
    /// public data member (a much narrower, opt-in exposure).
    ///
    /// NOT closed, NO downcast required:
    ///
    /// * `Reducing` (HelmholtzEOSMixtureBackend.h:147) is a PUBLIC
    ///   shared_ptr<ReducingFunction> data member.  A caller with a
    ///   GERGMixtureBackend*/HelmholtzEOSMixtureBackend* can do
    ///   `heos->Reducing->set_binary_interaction_double(i, j, param, value)`
    ///   directly and mutate the GERG2008ReducingFunction in place, bypassing
    ///   every guard above entirely.  This is a CHOICE not to close, not an
    ///   impossibility: `ConstantReducingFunction::set_binary_interaction_double`
    ///   (ReducingFunctions.h:558, `{ return; }`) already establishes the
    ///   pattern of a ReducingFunction subclass that refuses mutation itself,
    ///   so a GERG-specific ReducingFunction wrapper that overrides
    ///   set_binary_interaction_double to throw would close this route
    ///   without touching HelmholtzEOSMixtureBackend at all. Not built here:
    ///   it is a second class (wrapping or subclassing GERG2008ReducingFunction)
    ///   for one hard-to-reach hole, and `residual_helmholtz` below has the
    ///   identical problem and would need the identical treatment, doubling
    ///   the scope.  Left as a documented, tested (CHECK_NOTHROW) limitation.
    /// * `residual_helmholtz` (HelmholtzEOSMixtureBackend.h:148) is equally
    ///   public, and `residual_helmholtz->Excess.F[i][j]` /
    ///   `.DepartureFunctionMatrix[i][j]` are themselves public members of
    ///   `ExcessTerm` (ExcessHEFunction.h) with no setter to intercept at all
    ///   -- confirmed live: `h->residual_helmholtz->Excess.F[0][1] = 0` on
    ///   GERG2008 Methane/Nitrogen changes alphar (-0.259933 -> -0.257183)
    ///   with no error.  This is the exact class of silent substitution the
    ///   sanctioned `set_binary_interaction_string(..., "function", ...)`
    ///   guard above exists to prevent -- reachable here without going
    ///   through any setter whatsoever.  Same "known limitation, not hidden"
    ///   treatment as `Reducing`.
    ///
    /// CLOSED by this task (round 2): `AbstractState::change_EOS` (see
    /// calc_change_EOS override above) and `update_with_guesses` (see
    /// override above) were BOTH reachable with no downcast and are now
    /// guarded.
    ///
    /// NOT closed, downcast to HelmholtzEOSMixtureBackend* REQUIRED (a much
    /// narrower exposure -- documented, not contorted around): `update_TP_guessrho`,
    /// `update_DmolarT_direct` (confirmed live: accepts T = 900 K, returns
    /// p = 4.29e7 Pa with no error), and `update_TDmolarP_unchecked`
    /// (HelmholtzEOSMixtureBackend.h:408-410) are public, non-virtual, and
    /// not declared anywhere on `AbstractState` -- so calling any of them
    /// needs source code that already knows it is holding a GERG/HEOS-family
    /// object, not merely an `AbstractState`.  Overriding them individually
    /// would mean re-deriving three separate flash entry points' worth of
    /// range-check plumbing for a caller sophistication level ("I have
    /// downcast to the concrete backend type") that has already opted out of
    /// the ordinary public API's guarantees.
    /// ==================================================================
    void set_binary_interaction_double(const std::size_t, const std::size_t, const std::string&, const double) override {
        throw ValueError("GERG binary interaction parameters are fixed by the published model and cannot be modified. A mixture with altered "
                         "beta/gamma is not GERG.");
    }
    void set_binary_interaction_string(const std::size_t, const std::size_t, const std::string&, const std::string&) override {
        throw ValueError("GERG binary interaction parameters are fixed by the published model and cannot be modified.");
    }

   protected:
    /// R = 8.314472 J/mol/K, ALWAYS -- for mixtures as well as pure fluids.
    ///
    /// This is not redundant with `EquationOfState::R_u`.
    /// HelmholtzEOSMixtureBackend::calc_gas_constant returns the component's
    /// own R only when the state is pure; for a MIXTURE it returns
    /// `get_config_double(R_U_CODATA)` (8.31446261815324) whenever the global
    /// NORMALIZE_GAS_CONSTANTS config flag is set, which is CoolProp's default
    /// (HelmholtzEOSMixtureBackend.cpp:599-606).  That silently rescales p, the
    /// ideal-gas R_i/R_mix ratios in calc_alpha0_deriv_nocache, c_v and w by
    /// ~1.1e-6 relative -- large enough to be wrong, small enough to look like
    /// a plausible answer -- while leaving alpha^r untouched.  GERG defines its
    /// own R, so the backend must not be at the mercy of a global setting.
    CoolPropDbl calc_gas_constant() override;

    /// Build the reducing function and excess (departure) term from the
    /// published GERG tables ONLY.  Overriding this is the whole point of the
    /// backend: the inherited implementation resolves each binary pair through
    /// CoolProp's global BIP library.  All 210 GERG-2008 pairs are in that
    /// library, but 16 of them carry a LATER refit (15 Gernert-Thesis-2013,
    /// 1 Tkaczuk-JPCRD-2020) instead of the Kunz-JCED-2012 row GERG publishes
    /// -- and the departure function attached to a pair can differ even where
    /// the betas and gammas agree.  Reaching the inherited path would answer
    /// with those, silently, under a backend named GERG2004/GERG2008.
    void set_mixture_parameters() override;

    /// Overridden so the linked SatL/SatV states are GERGMixtureBackend
    /// objects.  The base implementation builds them with an explicitly
    /// qualified `HelmholtzEOSMixtureBackend::get_copy(false)` call
    /// (HelmholtzEOSMixtureBackend.cpp:142, :146), which constructs a
    /// base-class object whose own set_components would run the INHERITED
    /// set_mixture_parameters -- i.e. it would either throw (GERG fluids carry
    /// no CAS) or, worse, silently answer from the global BIP library.
    void set_components(const std::vector<CoolPropFluid>& components, bool generate_SatL_and_SatV = true) override;

   private:
    /// Throw if the current state's T is outside make_gerg_fluid's published
    /// range.  See the `update` override above for why this exists as a
    /// check here rather than relying on the inherited flash machinery, and
    /// for why it is T-only.
    ///
    /// The bounds are the INTERSECTION of the components' own
    /// EOS.limits.Tmin/Tmax (max of the Tmin values, min of the Tmax
    /// values), NOT AbstractState::Tmin()/Tmax() -- which for a mixture are
    /// un-normalised mole-fraction-weighted sums and therefore scale with
    /// sum(x).  The rationale for both halves of that choice is written out
    /// at the definition in GERGBackend.cpp.
    void check_gerg_range_of_validity();

    /// Build directly from already-assembled GERG CoolPropFluid objects.
    /// Used by set_components/get_copy to make the linked states; private
    /// because a caller outside this class has no way to obtain GERG fluids
    /// that were not built by make_gerg_fluid.
    GERGMixtureBackend(GERGModel model, const std::vector<CoolPropFluid>& fluids, bool generate_SatL_and_SatV);

    /// make_gerg_fluid needs this constructor for one specific job: filling in
    /// EOS.hs_anchor.hmolar/smolar and EOS.reduce.hmolar/smolar, which can only
    /// be obtained by EVALUATING the fluid it has just assembled (via
    /// update_states()).  The friendship is narrow on purpose -- it grants
    /// make_gerg_fluid, and nothing else, the ability to build a backend from
    /// raw fluids, and make_gerg_fluid is the one function whose whole job is
    /// producing those fluids in the first place.  Widening the constructor to
    /// public instead would let any caller hand in arbitrary CoolPropFluid
    /// objects -- e.g. CoolProp's own reference-EOS fluids -- and get them
    /// treated as GERG components with GERG's gas constant.
    friend CoolPropFluid GERG::make_gerg_fluid(GERGModel model, const std::string& gerg_name);

    GERGModel m_model;
};

} /* namespace CoolProp */
#endif /* GERGBACKEND_H_ */
