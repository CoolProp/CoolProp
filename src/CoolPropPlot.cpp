#include "CoolProp/plotting/CoolPropPlot.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/numerics/numerics.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>

// The numerics of CoolPlot's thermo layer (thermo/fluid.py, tracer.py,
// stepping.py, isolines.py of CoolPlot 0.7), written for readers who know
// those files: the names and the order of the steps are kept.

namespace CoolProp {
namespace Plot {

namespace Detail {

const double NaN = std::numeric_limits<double>::quiet_NaN();

inline bool finite(double v) {
    return std::isfinite(v);
}

Scale default_scale(CoolProp::parameters key) {
    switch (key) {
        case CoolProp::iDmass:
        case CoolProp::iDmolar:
        case CoolProp::iP:
            return Scale::Log;
        default:
            return Scale::Lin;
    }
}

// ---------------------------------------------------------------------------
// The fluid: a CoolProp state with the caches CoolPlot's Fluid keeps
// ---------------------------------------------------------------------------

/// The two-phase boundary of a mixture from CoolProp's phase envelope
/// (CoolPlot's PhaseEnvelope). The envelope runs along the dew line
/// (Q = 1) up through the critical region and back down the bubble line
/// (Q = 0); the bulk phase sits in the "vap" fields on both branches.
class Envelope
{
   public:
    std::vector<double> T, p, Q, rhomolar;
    bool closed = false;  ///< the two densities meet at the hottest point: it ran to the critical point
    size_t hot = 0;       ///< index of the hottest point

    Envelope(const std::shared_ptr<CoolProp::AbstractState>& work) : work_(work) {
        const double CLOSED_RTOL = 0.35;
        const CoolProp::PhaseEnvelopeData& data = work_->get_phase_envelope_data();
        T = data.T;
        p = data.p;
        Q = data.Q;
        rhomolar = data.rhomolar_vap;
        if (T.empty() || p.size() != T.size() || Q.size() != T.size() || rhomolar.size() != T.size() || data.rhomolar_liq.size() != T.size())
            throw CoolProp::ValueError("The phase envelope is empty or its columns differ in length.");
        hot = static_cast<size_t>(std::max_element(T.begin(), T.end()) - T.begin());
        double rho_liq = data.rhomolar_liq[hot], rho_vap = rhomolar[hot];
        closed = std::fabs(rho_liq - rho_vap) <= CLOSED_RTOL * std::max(rho_liq, rho_vap);
        bool has_bubble = false, has_dew = false;
        for (double q : Q) {
            has_bubble = has_bubble || q == 0.0;
            has_dew = has_dew || q == 1.0;
        }
        if (!(has_bubble && has_dew)) throw CoolProp::ValueError("The phase envelope has no bubble or no dew line.");
    }

    /// A property along the envelope in SI units, NaN where it fails: T, p and
    /// density from the envelope itself, the rest evaluated at (T, rho) with
    /// the phase imposed (no flash).
    const std::vector<double>& values(CoolProp::parameters key) {
        if (key == CoolProp::iT) return T;
        if (key == CoolProp::iP) return p;
        if (key == CoolProp::iQ) return Q;
        auto found = cache_.find(key);
        if (found != cache_.end()) return found->second;
        std::vector<double> out(T.size(), NaN);
        double M = work_->molar_mass();
        for (size_t i = 0; i < T.size(); ++i) {
            if (key == CoolProp::iDmass) {
                out[i] = rhomolar[i] * M;
                continue;
            }
            try {
                work_->specify_phase(Q[i] == 0.0 ? CoolProp::iphase_liquid : CoolProp::iphase_gas);
                work_->update(CoolProp::DmolarT_INPUTS, rhomolar[i], T[i]);
                out[i] = work_->keyed_output(key);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // no value at this point: it stays NaN
            }
        }
        work_->unspecify_phase();
        return cache_[key] = out;
    }
    double top(CoolProp::parameters key) {
        return values(key)[hot];
    }

   private:
    std::shared_ptr<CoolProp::AbstractState> work_;
    std::map<CoolProp::parameters, std::vector<double>> cache_;
};

/// The critical point: T, p, rhomolar from CoolProp, the rest evaluated there when possible.
struct CriticalPoint
{
    double T = NaN, p = NaN, rhomolar = NaN;
    std::map<CoolProp::parameters, double> values;  // h, s, rho, u where known
    [[nodiscard]] double get(CoolProp::parameters key) const {
        if (key == CoolProp::iT) return T;
        if (key == CoolProp::iP) return p;
        auto found = values.find(key);
        return found == values.end() ? NaN : found->second;
    }
};

/**
 * A fluid model with its own CoolProp states and cached constants
 * (CoolPlot's thermo.fluid.Fluid). ``update_pair`` retries a mixture on a
 * scratch state with the phase envelope built, and refuses a mixture state
 * CoolProp got wrong (reported single-phase inside the dome).
 */
class Fluid
{
   public:
    explicit Fluid(const std::string& fluid_ref) : name_(fluid_ref) {
        std::string fluids;
        CoolProp::extract_backend(fluid_ref, backend_, fluids);
        std::vector<double> fractions;
        fluids = CoolProp::extract_fractions(fluids, fractions);
        // Not a member initializer: backend_ is only known after extract_backend above
        mass_basis_ = backend_.rfind("INCOMP", 0) == 0;  // NOLINT(cppcoreguidelines-prefer-member-initializer)
        std::shared_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory(backend_, fluids));
        if (!fractions.empty()) {
            if (mass_basis_)
                state->set_mass_fractions(fractions);
            else
                state->set_mole_fractions(fractions);
        }
        // A predefined blend such as "R407C.mix" is one name for several
        // components; describe it by those, so that fresh states of the
        // same model can be built (the tracer, the envelope).
        names_ = {fluids};
        fractions_ = fractions;
        if (!mass_basis_) {
            try {
                std::vector<std::string> components = state->fluid_names();
                if (components.size() > 1) {
                    names_ = components;
                    const std::vector<CoolPropDbl>& x = state->get_mole_fractions();
                    fractions_.assign(x.begin(), x.end());
                }
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // keep the name and fractions as given
            }
        }
        state_ = state;
        current_ = state_;
    }

    [[nodiscard]] const std::string& name() const {
        return name_;
    }
    [[nodiscard]] bool is_pure() const {
        return names_.size() <= 1;
    }
    [[nodiscard]] bool mass_basis() const {
        return mass_basis_;
    }
    [[nodiscard]] const std::shared_ptr<CoolProp::AbstractState>& state() const {
        return state_;
    }

    /// A fresh state of the same model.
    [[nodiscard]] std::shared_ptr<CoolProp::AbstractState> make_state() const {
        std::shared_ptr<CoolProp::AbstractState> fresh(CoolProp::AbstractState::factory(backend_, names_));
        if (!fractions_.empty()) {
            if (mass_basis_)
                fresh->set_mass_fractions(fractions_);
            else
                fresh->set_mole_fractions(fractions_);
        }
        return fresh;
    }

    /// Set the working state from an input pair; the state that holds the result is returned.
    CoolProp::AbstractState& update_pair(CoolProp::input_pairs pair, double v1, double v2) {
        current_ = state_;
        try {
            state_->update(pair, v1, v2);
        } catch (...) {
            std::shared_ptr<CoolProp::AbstractState> guide = is_pure() ? nullptr : guide_state();
            if (!guide) throw;
            guide->update(pair, v1, v2);
            current_ = guide;
        }
        std::string problem = inconsistency(*current_);
        if (!problem.empty()) throw CoolProp::ValueError("CoolProp returned an inconsistent state for " + name_ + ": " + problem);
        // The check used the guide state for its saturation lookups; when the
        // guide holds the result, put the requested state back on it
        if (current_ != state_) current_->update(pair, v1, v2);
        return *current_;
    }
    CoolProp::AbstractState& update(CoolProp::parameters k1, double v1, CoolProp::parameters k2, double v2) {
        double o1 = 0.0, o2 = 0.0;
        CoolProp::input_pairs pair = CoolProp::generate_update_pair(k1, v1, k2, v2, o1, o2);
        return update_pair(pair, o1, o2);
    }
    [[nodiscard]] double get(CoolProp::parameters key) const {
        return current_->keyed_output(key);
    }

    /// Why a state of a mixture cannot be right (empty if it passes): CoolProp 8
    /// reports some inputs inside the two-phase region as a single phase.
    std::string inconsistency(CoolProp::AbstractState& st) {
        if (is_pure()) return "";
        double Q = NaN, T = NaN, rho = NaN;
        try {
            Q = st.Q();
            if (Q >= 0.0 && Q <= 1.0) return "";
            T = st.T();
            rho = st.rhomolar();
        } catch (...) {
            return "";
        }
        if (!has_critical_point() || T >= critical().T * (1.0 - 1e-6)) return "";
        std::shared_ptr<CoolProp::AbstractState> guide = guide_state();
        if (!guide) return "";
        double rho_vap, rho_liq;
        try {
            guide->update(CoolProp::QT_INPUTS, 1.0, T);
            rho_vap = guide->rhomolar();
            guide->update(CoolProp::QT_INPUTS, 0.0, T);
            rho_liq = guide->rhomolar();
        } catch (...) {
            return "";
        }
        if (rho_vap * (1.0 + 1e-6) < rho && rho < rho_liq * (1.0 - 1e-6)) {
            std::ostringstream text;
            text << "reported as single phase at T = " << T << " K with rho = " << rho << " mol/m3, between the saturated vapour (" << rho_vap
                 << ") and liquid (" << rho_liq << ") densities, i.e. inside the two-phase region";
            return text.str();
        }
        return "";
    }

    /// A scratch state with the phase envelope built, made once; null if that fails.
    std::shared_ptr<CoolProp::AbstractState> guide_state() {
        if (!guide_tried_) {
            guide_tried_ = true;
            try {
                std::shared_ptr<CoolProp::AbstractState> guide = make_state();
                guide->build_phase_envelope("veryfine");
                guide_ = guide;
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // no guide state: callers check for null
            }
        }
        return guide_;
    }

    /// The phase envelope of a mixture, computed once; null for pure fluids or on failure.
    Envelope* envelope() {
        if (!envelope_tried_) {
            envelope_tried_ = true;
            if (!is_pure()) {
                for (const char* level : {"veryfine", "none"}) {
                    try {
                        std::shared_ptr<CoolProp::AbstractState> work = make_state();
                        work->build_phase_envelope(level);
                        envelope_ = std::make_unique<Envelope>(work);
                        break;
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        // try the next, coarser level
                    }
                }
            }
        }
        return envelope_.get();
    }

    /// The highest stable critical temperature CoolProp finds for the fluid
    /// (all_critical_points, seconds to minutes for a blend with a wide
    /// glide), computed once; NaN if there is none.
    double dome_critical_T() {
        if (!dome_critical_tried_) {
            dome_critical_tried_ = true;
            try {
                for (const CoolProp::CriticalState& c : state_->all_critical_points())
                    if (c.stable && !(c.T <= dome_critical_T_)) dome_critical_T_ = c.T;
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // no stable critical point: stays NaN
            }
        }
        return dome_critical_T_;
    }

    bool has_critical_point() {
        if (!critical_tried_) {
            try {
                critical();
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // has_critical_ stays false
            }
        }
        return has_critical_;
    }

    const CriticalPoint& critical() {
        if (!critical_tried_) {
            critical_tried_ = true;
            CriticalPoint found;
            try {
                found.T = state_->T_critical();
                found.p = state_->p_critical();
                found.rhomolar = state_->rhomolar_critical();
            } catch (...) {
                try {
                    for (const CoolProp::CriticalState& c : state_->all_critical_points())
                        if (c.stable && !(c.T <= found.T)) {
                            found.T = c.T;
                            found.p = c.p;
                            found.rhomolar = c.rhomolar;
                        }
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                    // found.T stays NaN, reported below
                }
            }
            if (!finite(found.T)) throw CoolProp::ValueError("CoolProp could not provide the critical point of " + name_ + ".");
            std::shared_ptr<CoolProp::AbstractState> crit = make_state();
            bool evaluated = false;
            std::vector<std::pair<CoolProp::input_pairs, std::pair<double, double>>> attempts;
            if (finite(found.rhomolar)) attempts.push_back({CoolProp::DmolarT_INPUTS, {found.rhomolar, found.T}});
            if (finite(found.p)) attempts.push_back({CoolProp::PT_INPUTS, {found.p, found.T}});
            for (auto& attempt : attempts) {
                for (bool impose : {true, false}) {
                    try {
                        if (impose)
                            crit->specify_phase(CoolProp::iphase_critical_point);
                        else
                            crit->unspecify_phase();
                        crit->update(attempt.first, attempt.second.first, attempt.second.second);
                        for (CoolProp::parameters k : {CoolProp::iHmass, CoolProp::iSmass, CoolProp::iDmass, CoolProp::iUmass})
                            found.values[k] = crit->keyed_output(k);
                        evaluated = true;
                        break;
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        // try the next attempt
                    }
                }
                if (evaluated) break;
            }
            if (!evaluated && !finite(found.p)) throw CoolProp::ValueError("Could not evaluate the critical point of " + name_ + ".");
            critical_ = found;
            has_critical_ = true;
        }
        if (!has_critical_) throw CoolProp::ValueError("CoolProp could not provide the critical point of " + name_ + ".");
        return critical_;
    }

    /// A trivial parameter of the working state, NaN if the backend has none.
    [[nodiscard]] double trivial(CoolProp::parameters key) const {
        try {
            double v = state_->trivial_keyed_output(key);
            return finite(v) ? v : NaN;
        } catch (...) {
            return NaN;
        }
    }
    /// A lower temperature limit: for a mixture the highest one of its components.
    double lowest_temperature(CoolProp::parameters key) {
        double value = trivial(key);
        if (!is_pure()) {
            for (const std::string& component : names_) {
                try {
                    std::shared_ptr<CoolProp::AbstractState> one(CoolProp::AbstractState::factory(backend_, component));
                    double v = one->trivial_keyed_output(key);
                    if (finite(v) && (!finite(value) || v > value)) value = v;
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                    // skip a component without this limit
                }
            }
        }
        return value;
    }
    double T_triple_K() {
        return lowest_temperature(CoolProp::iT_triple);
    }
    double T_min_K() {
        return lowest_temperature(CoolProp::iT_min);
    }
    [[nodiscard]] double T_max_K() const {
        return trivial(CoolProp::iT_max);
    }
    [[nodiscard]] double p_min_Pa() const {
        return trivial(CoolProp::iP_min);
    }
    [[nodiscard]] double p_max_Pa() const {
        return trivial(CoolProp::iP_max);
    }

    /// Lowest and highest saturation temperature or pressure, inside the margin.
    Range saturation_range(CoolProp::parameters key, double margin) {
        if (!has_critical_point())
            throw CoolProp::ValueError(name_
                                       + " has no saturation curve (CoolProp provides no critical point for it); give absolute "
                                         "temperature and pressure limits.");
        const CriticalPoint& crit = critical();
        double T_lo = 0.0;
        for (double v : {T_triple_K(), T_min_K()})
            if (finite(v)) T_lo = std::max(T_lo, v);
        T_lo += crit.T * margin;
        if (key == CoolProp::iT) return {T_lo, crit.T * (1.0 - margin)};
        if (key != CoolProp::iP) throw CoolProp::ValueError("Saturation ranges are defined in T or p.");
        update(CoolProp::iQ, 0.0, CoolProp::iT, T_lo);
        return {get(CoolProp::iP) * (1.0 + margin), crit.p * (1.0 - margin)};
    }

   private:
    std::string name_, backend_;
    std::vector<std::string> names_;
    std::vector<double> fractions_;
    bool mass_basis_ = false;
    std::shared_ptr<CoolProp::AbstractState> state_, current_, guide_;
    bool guide_tried_ = false, envelope_tried_ = false, critical_tried_ = false, has_critical_ = false, dome_critical_tried_ = false;
    double dome_critical_T_ = NaN;
    std::unique_ptr<Envelope> envelope_;
    CriticalPoint critical_;
};

std::shared_ptr<CoolProp::AbstractState> process_fluid_state(const std::string& fluid_ref) {
    return Fluid(fluid_ref).make_state();
}

// ---------------------------------------------------------------------------
// The tracer: walking along an isoline of a mixture (CoolPlot's thermo/tracer.py)
// ---------------------------------------------------------------------------

/// A two-phase state whose two phases differ in density by less than this,
/// relatively, is the trivial solution of the mixture solver.
const double TRIVIAL_SOLUTION_RTOL = 1e-4;

bool phases_are_distinct(CoolProp::AbstractState& st) {
    double rho_liq, rho_vap;
    try {
        rho_liq = st.saturated_liquid_keyed_output(CoolProp::iDmolar);
        rho_vap = st.saturated_vapor_keyed_output(CoolProp::iDmolar);
    } catch (...) {
        // HEOS, the cubics, GERG and REFPROP give these outputs after a
        // quality flash, so a failure there is a failed state.  A backend
        // without saturated outputs (IF97, PC-SAFT) throws every time.  Both
        // count as not distinct: that costs the caller a fallback flash or a
        // gap, never a collapsed state taken as a point.
        return false;
    }
    if (!(finite(rho_liq) && finite(rho_vap)) || rho_liq <= 0.0 || rho_vap <= 0.0) return false;
    return (rho_liq - rho_vap) > TRIVIAL_SOLUTION_RTOL * std::max(rho_liq, rho_vap);
}

/**
 * Walk along an isoline, warm-starting every point from the previous one.
 * A single-phase point is a two-by-two Newton iteration in (T, rhomolar)
 * with the phase imposed; a two-phase point is bracketed by the saturation
 * states at the constant temperature or pressure and its quality found by
 * a secant iteration from the previous quality (a regula falsi over the
 * whole dome as the fallback). Every point the tracer cannot solve with
 * confidence throws, and the caller flashes it instead.
 */
class IsoLineTracer
{
   public:
    static constexpr double NEWTON_RTOL = 1e-10;
    static constexpr double RESIDUAL_RTOL = 1e-9;
    static constexpr int NEWTON_ITMAX = 30;
    static constexpr double NEWTON_TRUST = 0.25;
    static constexpr double QUALITY_XTOL = 1e-12;
    static constexpr double QUALITY_RTOL = 1e-8;
    static constexpr int QUALITY_ITMAX = 60;
    static constexpr double BRACKET_RTOL = 1e-8;
    static constexpr double BRANCH_RTOL = 1e-6;
    static constexpr double DOME_MARGIN = 1.05;
    static constexpr double CLOSED_ENVELOPE_MARGIN = 1.005;
    static constexpr double ENVELOPE_CLOSED_RTOL = 0.35;
    static constexpr double TOP_DENSITY_MARGIN = 1.5;
    static constexpr double ABOVE_TOP_RTOL = 2e-4;
    static constexpr double TOP_GAP_RTOL = 0.05;
    static constexpr double RANGE_FACTOR = 2.0;
    static const int LIQUID = -1;
    static const int VAPOUR = 1;

    static bool supports(CoolProp::parameters i1, CoolProp::parameters i2) {
        return i1 == CoolProp::iT || i2 == CoolProp::iT || i1 == CoolProp::iP || i2 == CoolProp::iP;
    }

    IsoLineTracer(Fluid& fluid, CoolProp::parameters index1, CoolProp::parameters index2, CoolProp::parameters iso_index)
      : i1_(index1), i2_(index2), critical_T_(fluid.dome_critical_T()) {
        if (index1 == index2) throw CoolProp::ValueError("The two inputs of an isoline cannot be the same property.");
        if (!supports(index1, index2))
            throw CoolProp::ValueError("An isoline whose inputs include neither temperature nor pressure cannot be traced.");
        sat_ = fluid.make_state();
        sp_ = fluid.make_state();
        try {
            sp_->specify_phase(CoolProp::iphase_gas);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // the backend cannot impose a phase: leave it free
        }
        if ((iso_index == CoolProp::iT || iso_index == CoolProp::iP) && (iso_index == index1 || iso_index == index2))
            sat_index_ = iso_index;
        else if (index1 == CoolProp::iT || index2 == CoolProp::iT)
            sat_index_ = CoolProp::iT;
        else
            sat_index_ = CoolProp::iP;
        find_eos_bounds();
        find_saturation_bounds();
        find_dome_temperature_limit();
    }

    /// The state holding the most recently traced point.
    CoolProp::AbstractState& state() {
        if (!out_) throw CoolProp::ValueError("No point has been traced yet.");
        return *out_;
    }

    /// Adopt a single-phase state the caller solved itself as the warm start.
    void set_guess(CoolProp::AbstractState& st) {
        double T, rho;
        try {
            if (st.phase() == CoolProp::iphase_twophase) return;
            double Q = st.Q();
            if (finite(Q) && Q >= 0.0 && Q <= 1.0) return;
            T = st.T();
            rho = st.rhomolar();
        } catch (...) {
            return;
        }
        if (finite(T) && finite(rho) && T > 0.0 && rho > 0.0) {
            T_ = T;
            rhomolar_ = rho;
            side_ = 0;
        }
    }

    /// Trace the state with index1 == value1 and index2 == value2; throws if it cannot.
    void update(double value1, double value2) {
        out_ = nullptr;
        std::map<CoolProp::parameters, double> targets = {{i1_, value1}, {i2_, value2}};
        Bracket bracket;
        bool have_bracket = get_bracket(targets, bracket);
        if (have_bracket && bracket.two_phase)
            solve_quality(bracket);
        else
            solve_single_phase(targets, have_bracket ? &bracket : nullptr);
    }

   private:
    struct End
    {
        double value, T, rhomolar;
    };
    struct Bracket
    {
        CoolProp::parameters target_index = CoolProp::iundefined_parameter;
        End liq{}, vap{};
        bool degenerate = false;
        double slope = NaN;  // NaN when unavailable
        double target = NaN;
        bool two_phase = false;
        int side = 0;  // LIQUID, VAPOUR or 0
        std::function<double(double)> scale;
    };

    CoolProp::parameters i1_, i2_, sat_index_;
    double critical_T_;  // the highest stable critical temperature of the fluid, NaN if none (cached by the Fluid)
    std::shared_ptr<CoolProp::AbstractState> sat_, sp_;
    CoolProp::AbstractState* out_ = nullptr;
    std::map<CoolProp::parameters, double> eos_bounds_;
    double eos_T_min_ = NaN;
    std::map<CoolProp::parameters, double> saturation_bounds_;  // the outer bound
    double dome_limit_ = NaN;                                   // inner bound in the saturation variable
    double dome_T_limit_ = NaN;
    bool have_top_ = false;
    double top_T_ = NaN, top_p_ = NaN, top_rho_liq_ = NaN, top_rho_vap_ = NaN;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> dome_densities_;  // (T, bulk rho) per branch, sorted
    double sat_value_ = NaN;
    bool have_sat_value_ = false;
    Bracket cached_bracket_;
    bool cached_is_null_ = true;
    std::string bracket_failure_;
    double T_ = NaN, rhomolar_ = NaN;
    int side_ = 0;
    double last_quality_ = NaN;
    double last_two_phase_T_ = NaN, last_two_phase_rho_ = NaN;

    void find_eos_bounds() {
        for (auto pair : {std::make_pair(CoolProp::iT, CoolProp::iT_max), std::make_pair(CoolProp::iP, CoolProp::iP_max)}) {
            try {
                double v = sp_->trivial_keyed_output(pair.second);
                if (finite(v) && v > 0.0) eos_bounds_[pair.first] = v;
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // no limit from the equation of state
            }
        }
        try {
            double v = sp_->trivial_keyed_output(CoolProp::iT_min);
            if (finite(v) && v > 0.0) eos_T_min_ = v;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // no limit from the equation of state
        }
    }

    bool envelope_is_closed(const CoolProp::PhaseEnvelopeData& data, size_t& hot) const {
        if (data.T.empty() || data.rhomolar_liq.size() != data.T.size() || data.rhomolar_vap.size() != data.T.size()
            || data.Q.size() != data.T.size())
            return false;
        hot = static_cast<size_t>(std::max_element(data.T.begin(), data.T.end()) - data.T.begin());
        double rho_liq = data.rhomolar_liq[hot], rho_vap = data.rhomolar_vap[hot];
        return std::fabs(rho_liq - rho_vap) <= ENVELOPE_CLOSED_RTOL * std::max(rho_liq, rho_vap);
    }

    void find_saturation_bounds() {
        std::map<CoolProp::parameters, double> envelope;
        bool closed = false;
        // Building the envelope on sat_ is also what makes the mixture saturation
        // flashes on it converge (they take their guesses from it), so every
        // tracer builds its own and none is shared or cached: with a shared
        // state a line would depend on the lines computed before it, and lines
        // are cached and computed one by one in any order.  It costs about a
        // quarter of the time of a traced line (CoolPlot's tracer.py, same note).
        try {
            sat_->build_phase_envelope("none");
            const CoolProp::PhaseEnvelopeData& data = sat_->get_phase_envelope_data();
            if (data.T.empty() || data.p.size() != data.T.size()) throw CoolProp::ValueError("The phase envelope is empty.");
            envelope[CoolProp::iT] = *std::max_element(data.T.begin(), data.T.end());
            envelope[CoolProp::iP] = *std::max_element(data.p.begin(), data.p.end());
            size_t hot = 0;
            closed = envelope_is_closed(data, hot);
            if (closed) {
                have_top_ = true;
                top_T_ = data.T[hot];
                top_p_ = envelope[CoolProp::iP];
                top_rho_liq_ = data.rhomolar_liq[hot];
                top_rho_vap_ = data.rhomolar_vap[hot];
                dome_density_table(data);
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // no envelope: the bounds come from the components below
        }
        // The outer bound: the largest critical point of the components, when
        // every one of them can be looked up (Props1SI reports a failure as an
        // infinite value, not an exception; then the envelope's bound stays).
        std::map<CoolProp::parameters, double> outer = envelope;
        try {
            std::map<CoolProp::parameters, double> components;
            bool complete = true;
            for (auto pair : {std::make_pair(CoolProp::iT, "Tcrit"), std::make_pair(CoolProp::iP, "pcrit")}) {
                double largest = 0.0;
                for (const std::string& fluid : sat_->fluid_names()) {
                    double v = CoolProp::Props1SI(fluid, pair.second);
                    complete = complete && finite(v) && v < 1e300;
                    largest = std::max(largest, v);
                }
                components[pair.first] = largest;
            }
            if (complete)
                for (auto& kv : components)
                    outer[kv.first] = std::max(outer.count(kv.first) != 0 ? outer[kv.first] : 0.0, kv.second);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // a lookup failed: keep the envelope's bound
        }
        auto clean = [&](const std::map<CoolProp::parameters, double>& candidate, double margin) {
            std::map<CoolProp::parameters, double> bounds;
            for (auto& kv : candidate) {
                if (!(finite(kv.second) && kv.second > 0.0)) continue;
                double value = kv.second;
                auto limit = eos_bounds_.find(kv.first);
                if (limit != eos_bounds_.end()) value = std::min(value, limit->second);
                bounds[kv.first] = value * margin;
            }
            return bounds;
        };
        std::map<CoolProp::parameters, double> inner = clean(envelope, closed ? CLOSED_ENVELOPE_MARGIN : DOME_MARGIN);
        saturation_bounds_ = clean(outer, DOME_MARGIN);
        auto found = inner.find(sat_index_);
        dome_limit_ = found == inner.end() ? NaN : found->second;
    }

    void dome_density_table(const CoolProp::PhaseEnvelopeData& data) {
        for (double quality : {1.0, 0.0}) {
            std::vector<std::pair<double, double>> rows;
            for (size_t i = 0; i < data.T.size(); ++i)
                if (data.Q[i] == quality && finite(data.T[i]) && finite(data.rhomolar_vap[i])) rows.emplace_back(data.T[i], data.rhomolar_vap[i]);
            std::sort(rows.begin(), rows.end());
            std::vector<double> T, rho;
            for (auto& row : rows) {
                T.push_back(row.first);
                rho.push_back(row.second);
            }
            dome_densities_.emplace_back(T, rho);
        }
    }

    void find_dome_temperature_limit() {
        try {
            const CoolProp::PhaseEnvelopeData& data = sat_->get_phase_envelope_data();
            size_t hot = 0;
            if (envelope_is_closed(data, hot)) {
                dome_T_limit_ = data.T[hot] * CLOSED_ENVELOPE_MARGIN;
                return;
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // no closed envelope: fall back to the critical temperature
        }
        if (finite(critical_T_)) {
            dome_T_limit_ = critical_T_ * DOME_MARGIN;
            return;
        }
        auto found = saturation_bounds_.find(CoolProp::iT);
        dome_T_limit_ = found == saturation_bounds_.end() ? NaN : found->second;
    }

    [[nodiscard]] bool is_plausible_saturation(double T, double p) const {
        for (auto pair : {std::make_pair(CoolProp::iT, T), std::make_pair(CoolProp::iP, p)}) {
            auto bound = saturation_bounds_.find(pair.first);
            if (bound != saturation_bounds_.end() && !(pair.second <= bound->second)) return false;
        }
        return true;
    }

    void saturation_endpoint(double sat_value, double quality) {
        if (sat_index_ == CoolProp::iT)
            sat_->update(CoolProp::QT_INPUTS, quality, sat_value);
        else
            sat_->update(CoolProp::PQ_INPUTS, sat_value, quality);
    }

    [[nodiscard]] bool beyond_the_envelope_top(double sat_value) const {
        if (!have_top_) return false;
        return sat_value > (sat_index_ == CoolProp::iT ? top_T_ : top_p_);
    }

    /// The bracket at the saturation value, or false where there is no dome.
    /// Throws when the dome should be there and could not be found.
    bool get_bracket(const std::map<CoolProp::parameters, double>& targets, Bracket& bracket) {
        double sat_value = targets.at(sat_index_);
        CoolProp::parameters target_index = sat_index_ == i1_ ? i2_ : i1_;
        double target = targets.at(target_index);
        if (!have_sat_value_ || sat_value != sat_value_) {
            bracket_failure_.clear();
            cached_is_null_ = true;
            try {
                cached_is_null_ = !build_bracket(sat_value, target_index, cached_bracket_);
            } catch (CoolProp::CoolPropBaseError& e) {
                bracket_failure_ = e.what();
            } catch (std::exception& e) {
                bracket_failure_ = e.what();
            }
            sat_value_ = sat_value;
            have_sat_value_ = true;
        }
        if (!bracket_failure_.empty()) throw CoolProp::ValueError(bracket_failure_);
        if (cached_is_null_) return false;
        bracket = cached_bracket_;
        double v_liq = bracket.liq.value, v_vap = bracket.vap.value;
        bracket.target = target;
        bracket.two_phase = !bracket.degenerate && std::min(v_liq, v_vap) < target && target < std::max(v_liq, v_vap);
        if (bracket.two_phase)
            bracket.side = 0;
        else if (finite(bracket.slope))
            bracket.side = (target - 0.5 * (v_liq + v_vap)) * bracket.slope > 0 ? LIQUID : VAPOUR;
        else if (!bracket.degenerate)
            bracket.side = (target - v_liq) * (v_liq - v_vap) > 0 ? LIQUID : VAPOUR;
        else
            throw CoolProp::ValueError("Cannot tell which side of the saturation curve this point is on.");
        return true;
    }

    bool build_bracket(double sat_value, CoolProp::parameters target_index, Bracket& bracket) {
        End ends[2];
        try {
            int i = 0;
            for (double quality : {0.0, 1.0}) {
                saturation_endpoint(sat_value, quality);
                ends[i] = {sat_->keyed_output(target_index), sat_->T(), sat_->rhomolar()};
                if (!is_plausible_saturation(sat_->T(), sat_->p())) {
                    std::ostringstream text;
                    text << "The saturation solver returned a state at T=" << sat_->T() << " K, p=" << sat_->p() << " Pa, beyond any possible dome.";
                    throw CoolProp::ValueError(text.str());
                }
                ++i;
            }
        } catch (...) {
            if (finite(dome_limit_) && sat_value > dome_limit_) return false;
            if (beyond_the_envelope_top(sat_value)) return false;
            throw;
        }
        double v_liq = ends[0].value, v_vap = ends[1].value;
        if (!(finite(v_liq) && finite(v_vap))) throw CoolProp::ValueError("The saturation states are not finite.");
        if (!(ends[0].rhomolar > ends[1].rhomolar * (1.0 + TRIVIAL_SOLUTION_RTOL)))
            throw CoolProp::ValueError("The saturation states are not on their own branches.");
        bracket.target_index = target_index;
        bracket.liq = ends[0];
        bracket.vap = ends[1];
        bracket.degenerate = std::fabs(v_vap - v_liq) <= BRACKET_RTOL * std::max({std::fabs(v_liq), std::fabs(v_vap), 1e-30});
        bracket.slope = branch_slope(target_index, ends[0]);
        return true;
    }

    double branch_slope(CoolProp::parameters target_index, const End& liquid_end) {
        try {
            sp_->update(CoolProp::DmolarT_INPUTS, liquid_end.rhomolar, liquid_end.T);
            double slope = sp_->first_partial_deriv(target_index, CoolProp::iDmolar, sat_index_);
            return (finite(slope) && slope != 0.0) ? slope : NaN;
        } catch (...) {
            return NaN;
        }
    }

    // Two-phase point ----------------------------------------------------
    static std::function<double(double)> quality_scale(CoolProp::parameters target_index, double target) {
        if ((target_index == CoolProp::iDmolar || target_index == CoolProp::iDmass) && target > 0.0) return [](double v) { return 1.0 / v; };
        return [](double v) { return v; };
    }

    double quality_residual(double quality, const Bracket& bracket) {
        if (sat_index_ == CoolProp::iT)
            sat_->update(CoolProp::QT_INPUTS, quality, sat_value_);
        else
            sat_->update(CoolProp::PQ_INPUTS, sat_value_, quality);
        return bracket.scale(sat_->keyed_output(bracket.target_index)) - bracket.scale(bracket.target);
    }

    void solve_quality(Bracket& bracket) {
        bracket.scale = quality_scale(bracket.target_index, bracket.target);
        double target = bracket.scale(bracket.target);
        double fa = bracket.scale(bracket.liq.value) - target, fb = bracket.scale(bracket.vap.value) - target;
        if (fa * fb > 0.0) throw CoolProp::ValueError("The saturation bracket does not contain the target.");
        double width = std::fabs(fb - fa);
        double quality = NaN;
        if (!secant_quality(bracket, fa, fb, width, quality)) quality = bracketed_quality(bracket, fa, fb, width);
        out_ = sat_.get();
        last_quality_ = quality;
        try {
            last_two_phase_T_ = sat_->T();
            last_two_phase_rho_ = sat_->rhomolar();
        } catch (...) {
            last_two_phase_T_ = NaN;
        }
        T_ = NaN;
        rhomolar_ = NaN;
        side_ = 0;
    }

    bool secant_quality(const Bracket& bracket, double fa, double fb, double width, double& quality) {
        if (!finite(last_quality_)) return false;
        double tolerance = QUALITY_RTOL * width;
        double q0 = std::min(std::max(last_quality_, 1e-6), 1.0 - 1e-6);
        double f0 = quality_residual(q0, bracket);
        if (std::fabs(f0) <= tolerance) {
            quality = q0;
            return true;
        }
        double slope = fb - fa;
        double q1 = q0 - f0 / slope;
        for (int i = 0; i < 8; ++i) {
            if (!finite(q1) || !(q1 > 0.0 && q1 < 1.0)) return false;
            double f1 = quality_residual(q1, bracket);
            if (std::fabs(f1) <= tolerance) {
                quality = q1;
                return true;
            }
            if (f1 == f0) return false;
            slope = (f1 - f0) / (q1 - q0);
            q0 = q1;
            f0 = f1;
            q1 = q0 - f0 / slope;
        }
        return false;
    }

    double bracketed_quality(const Bracket& bracket, double fa, double fb, double width) {
        double a = 0.0, b = 1.0;
        bool converged = false;
        for (int i = 0; i < QUALITY_ITMAX; ++i) {
            double quality = b - fb * (b - a) / (fb - fa);
            double lo = std::min(a, b), hi = std::max(a, b);
            if (!finite(quality) || !(quality > lo && quality < hi)) quality = 0.5 * (lo + hi);
            double fq = quality_residual(quality, bracket);
            // NaN would compare as "same side" below and pass the residual gate
            if (!finite(fq)) throw CoolProp::ValueError("The quality solver met a residual that is not a number.");
            if (fq * fb < 0.0) {
                a = b;
                fa = fb;
            } else {
                fa *= 0.5;
            }
            b = quality;
            fb = fq;
            if (fq == 0.0 || std::fabs(b - a) < QUALITY_XTOL) {
                converged = true;
                break;
            }
        }
        if (!converged) throw CoolProp::ValueError("The quality solver did not converge.");
        if (!(std::fabs(fb) <= QUALITY_RTOL * width))
            throw CoolProp::ValueError("The quality solver ran out of bracket with the residual unresolved.");
        return b;
    }

    // Single-phase point -------------------------------------------------
    void seed(const Bracket* bracket, double& T, double& rhomolar) {
        int side = bracket ? bracket->side : 0;
        if (finite(T_) && finite(rhomolar_) && (side == 0 || side_ == 0 || side == side_)) {
            T = T_;
            rhomolar = rhomolar_;
            return;
        }
        if (bracket) {
            const End& end = side == LIQUID ? bracket->liq : bracket->vap;
            T = end.T;
            rhomolar = end.rhomolar;
            return;
        }
        if (finite(last_two_phase_T_)) {
            T = last_two_phase_T_;
            rhomolar = last_two_phase_rho_;
            return;
        }
        throw CoolProp::ValueError("No initial guess is available for the Newton iteration.");
    }

    void solve_single_phase(const std::map<CoolProp::parameters, double>& targets, const Bracket* bracket) {
        double T, rhomolar;
        seed(bracket, T, rhomolar);
        double v1 = targets.at(i1_), v2 = targets.at(i2_);
        for (int it = 0; it < NEWTON_ITMAX; ++it) {
            sp_->update(CoolProp::DmolarT_INPUTS, rhomolar, T);
            double f1 = sp_->keyed_output(i1_) - v1, f2 = sp_->keyed_output(i2_) - v2;
            double J11 = sp_->first_partial_deriv(i1_, CoolProp::iT, CoolProp::iDmolar);
            double J12 = sp_->first_partial_deriv(i1_, CoolProp::iDmolar, CoolProp::iT);
            double J21 = sp_->first_partial_deriv(i2_, CoolProp::iT, CoolProp::iDmolar);
            double J22 = sp_->first_partial_deriv(i2_, CoolProp::iDmolar, CoolProp::iT);
            double det = J11 * J22 - J12 * J21;
            if (!finite(det) || det == 0.0) throw CoolProp::ValueError("The Newton iteration hit a singular Jacobian.");
            double dT = (J12 * f2 - J22 * f1) / det, dD = (J21 * f1 - J11 * f2) / det;
            if (!(finite(dT) && finite(dD))) throw CoolProp::ValueError("The Newton iteration produced a non-finite step.");
            if (std::fabs(dT) <= NEWTON_RTOL * T && std::fabs(dD) <= NEWTON_RTOL * rhomolar) {
                double scale1 = std::fabs(J11) * T + std::fabs(J12) * rhomolar, scale2 = std::fabs(J21) * T + std::fabs(J22) * rhomolar;
                if (!(finite(scale1) && finite(scale2))) throw CoolProp::ValueError("The residual cannot be scaled at this point.");
                if (std::fabs(f1) > RESIDUAL_RTOL * scale1 || std::fabs(f2) > RESIDUAL_RTOL * scale2)
                    throw CoolProp::ValueError("The Newton iteration stalled with the residual unresolved.");
                check_in_range(T);
                int side = bracket ? bracket->side : 0;
                check_branch(bracket, side, T, rhomolar);
                T_ = T;
                rhomolar_ = rhomolar;
                side_ = side;
                out_ = sp_.get();
                return;
            }
            double damp = 1.0;
            if (std::fabs(dT) > NEWTON_TRUST * T) damp = std::min(damp, NEWTON_TRUST * T / std::fabs(dT));
            if (std::fabs(dD) > NEWTON_TRUST * rhomolar) damp = std::min(damp, NEWTON_TRUST * rhomolar / std::fabs(dD));
            T += damp * dT;
            rhomolar += damp * dD;
            if (!(T > 0.0 && rhomolar > 0.0)) throw CoolProp::ValueError("The Newton iteration left the physical domain.");
        }
        throw CoolProp::ValueError("The Newton iteration did not converge.");
    }

    void check_in_range(double T) const {
        auto upper = eos_bounds_.find(CoolProp::iT);
        if ((upper != eos_bounds_.end() && T > upper->second * RANGE_FACTOR) || (finite(eos_T_min_) && T < eos_T_min_ / RANGE_FACTOR))
            throw CoolProp::ValueError("The Newton iteration converged far outside the range of the equation of state.");
    }

    void check_branch(const Bracket* bracket, int side, double T, double rhomolar) {
        if (!bracket) {
            check_outside_the_dome(T, rhomolar);
            return;
        }
        if (side == 0) return;
        bool wrong = side == LIQUID ? rhomolar < bracket->liq.rhomolar * (1.0 - BRANCH_RTOL) : rhomolar > bracket->vap.rhomolar * (1.0 + BRANCH_RTOL);
        if (wrong) throw CoolProp::ValueError("The Newton iteration converged on the wrong side of the saturation curve.");
    }

    [[nodiscard]] bool far_from_the_dome(double T, double rhomolar) const {
        if (!have_top_ || dome_densities_.size() < 2) return false;
        bool top_is_critical = std::fabs(top_rho_liq_ - top_rho_vap_) <= TOP_GAP_RTOL * std::max(top_rho_liq_, top_rho_vap_);
        if (top_is_critical && T > top_T_ * (1.0 + ABOVE_TOP_RTOL)) return true;
        double lo = NaN, hi = NaN;
        for (auto& table : dome_densities_) {
            if (table.first.size() < 2) return false;
            double rho = interpolate_clamped(table.first, table.second, T);
            lo = finite(lo) ? std::min(lo, rho) : rho;
            hi = finite(hi) ? std::max(hi, rho) : rho;
        }
        return rhomolar > hi * TOP_DENSITY_MARGIN || rhomolar < lo / TOP_DENSITY_MARGIN;
    }

    static double interpolate_clamped(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
        if (x <= xs.front()) return ys.front();
        if (x >= xs.back()) return ys.back();
        size_t i = static_cast<size_t>(std::upper_bound(xs.begin(), xs.end(), x) - xs.begin());
        double t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
        return ys[i - 1] + t * (ys[i] - ys[i - 1]);
    }

    void check_outside_the_dome(double T, double rhomolar) {
        if (finite(dome_T_limit_) && T > dome_T_limit_) return;
        std::string undecided;
        double rho_liq = NaN, rho_vap = NaN, p_liq = NaN, p_vap = NaN;
        try {
            sat_->update(CoolProp::QT_INPUTS, 0.0, T);
            rho_liq = sat_->rhomolar();
            p_liq = sat_->p();
            sat_->update(CoolProp::QT_INPUTS, 1.0, T);
            rho_vap = sat_->rhomolar();
            p_vap = sat_->p();
        } catch (...) {
            undecided = "the saturation solver did not converge";
        }
        // A NaN would slip through both checks below: std::max drops a NaN
        // p_vap, and a NaN density makes the inside-the-dome test false
        if (undecided.empty() && !(finite(p_liq) && finite(p_vap) && finite(rho_liq) && finite(rho_vap)))
            undecided = "the saturation solver returned a value that is not a number";
        if (undecided.empty() && !is_plausible_saturation(T, std::max(p_liq, p_vap))) undecided = "the saturation solver left the dome behind";
        if (undecided.empty() && !phases_are_distinct(*sat_)) undecided = "the saturation states are not on their own branches";
        if (!undecided.empty()) {
            if (!finite(dome_T_limit_) || T <= dome_T_limit_) {
                if (far_from_the_dome(T, rhomolar)) return;
                throw CoolProp::ValueError("This root cannot be shown to be outside the two-phase dome: " + undecided + ".");
            }
            return;
        }
        if (rho_vap < rhomolar && rhomolar < rho_liq) throw CoolProp::ValueError("The Newton iteration converged inside the two-phase dome.");
    }
};

// ---------------------------------------------------------------------------
// Stepping along a line with steps chosen from its shape on the page
// (CoolPlot's thermo/stepping.py)
// ---------------------------------------------------------------------------

/// Normalised plot coordinates: each axis range mapped to 0..1, log for log axes.
struct Metric
{
    Range x, y;
    bool x_log, y_log;
    static double scale(double value, double lo, double hi, bool log) {
        if (log && lo > 0.0 && hi > 0.0 && value > 0.0) return (std::log10(value) - std::log10(lo)) / (std::log10(hi) - std::log10(lo));
        return (value - lo) / (hi - lo);
    }
    [[nodiscard]] std::pair<double, double> point(double px, double py) const {
        return {scale(px, x.min, x.max, x_log), scale(py, y.min, y.max, y_log)};
    }
};

/// One computed point of a line. ``valid`` false marks a gap in a list of samples.
struct Sample
{
    bool valid = false;
    double x = NaN, y = NaN;
    bool inside = true;
    int two_phase = -1;  // -1 unknown, 0 single phase, 1 two-phase
    double w = NaN, T = NaN, p = NaN, Q = NaN, rho = NaN;
    std::string source;
};

struct Trace
{
    std::vector<double> u;
    std::vector<Sample> samples;
    int n_evaluated = 0, n_outside = 0, n_failed = 0;
    int n_ok = 0;          // evaluations that gave a state: a failed look while scanning outside counts as outside, not failed
    int n_outside_ok = 0;  // evaluations that gave a state outside the domain (n_outside also counts the failed looks)
    std::string stopped;

    void add(double at, const Sample& sample) {
        if (!sample.valid && (samples.empty() || !samples.back().valid)) return;
        u.push_back(at);
        samples.push_back(sample);
    }
    void close() {
        while (!samples.empty() && !samples.back().valid) {
            u.pop_back();
            samples.pop_back();
        }
    }
};

const double MAX_ANGLE_DEG = 5.0;
const double STEP_FACTOR = 2.5;
const double MIN_STEP = 1e-4;
const double MAX_STEP = 1.0 / 8.0;
const double SCAN_STEP = 1.0 / 64.0;
const double SCAN_MAX = 1.0 / 32.0;
const double EDGE_TOLERANCE = 1e-3;

using Evaluate = std::function<bool(double, bool, Sample&)>;

/// The state of one call to trace().
class Walk
{
   public:
    Walk(Evaluate evaluate, const Metric& metric, int points, double u_lo, double u_hi, const Sample* start, const Sample* end, bool skip_failures)
      : evaluate_(std::move(evaluate)),
        metric_(metric),
        skip_failures_(skip_failures),
        u_lo_(u_lo),
        u_hi_(u_hi),
        span_(u_hi - u_lo),
        start_(start),
        end_(end),
        points_(std::max(points, 2)),
        max_length_(STEP_FACTOR / points_),
        max_cos_(std::cos(MAX_ANGLE_DEG * M_PI / 180.0)),
        budget_(8 * points_ + 100) {}

    Trace trace;

    void run() {
        double u = u_lo_;
        Sample sample;
        bool inside = sample_at(u, true, sample) && sample.inside;
        if (inside) accept(u, sample);
        while (u < u_hi_ && !over_budget()) {
            if (inside)
                walk_inside(u, inside);
            else
                scan_outside(u, inside);
        }
        if (u >= u_hi_ && trace.stopped.empty()) trace.stopped = "the end of the range";
        trace.close();
    }

   private:
    Evaluate evaluate_;
    Metric metric_;
    bool skip_failures_;
    double u_lo_, u_hi_, span_;
    const Sample* start_;
    const Sample* end_;
    int points_;
    double max_length_, max_cos_;
    int budget_;
    bool have_previous_ = false;
    Sample previous_;
    std::pair<double, double> previous_point_;
    bool have_direction_ = false;
    std::pair<double, double> direction_;
    double last_length_ = 0.0, last_turn_ = 1.0;

    bool sample_at(double u, bool careful, Sample& out) {
        if (start_ && u == u_lo_) {
            out = *start_;
            return true;
        }
        if (end_ && u == u_hi_) {
            out = *end_;
            return true;
        }
        ++trace.n_evaluated;
        bool ok = evaluate_(u, careful, out);
        if (ok) ++trace.n_ok;
        if (!ok) {
            if (careful)
                ++trace.n_failed;
            else
                ++trace.n_outside;
        } else if (!out.inside) {
            ++trace.n_outside;
            ++trace.n_outside_ok;
        }
        return ok;
    }

    bool over_budget() {
        if (trace.n_evaluated < budget_) return false;
        std::ostringstream text;
        text << "the evaluation budget of " << budget_ << " points is used up";
        trace.stopped = text.str();
        return true;
    }

    [[nodiscard]] bool step_ok(const std::pair<double, double>& point, const Sample& sample) const {
        if (!have_previous_) return true;
        double dx = point.first - previous_point_.first, dy = point.second - previous_point_.second;
        double length = std::hypot(dx, dy);
        if (!finite(length) || length > max_length_) return false;
        if (!have_direction_ || length < 1e-12) return true;
        if (previous_.two_phase >= 0 && sample.two_phase >= 0 && sample.two_phase != previous_.two_phase) return true;
        double cos = (dx * direction_.first + dy * direction_.second) / length;
        return cos >= max_cos_;
    }

    void accept(double u, const Sample& sample) {
        std::pair<double, double> point = metric_.point(sample.x, sample.y);
        last_length_ = 0.0;
        last_turn_ = 1.0;
        if (have_previous_) {
            double dx = point.first - previous_point_.first, dy = point.second - previous_point_.second;
            double length = std::hypot(dx, dy);
            if (length > 1e-12) {
                std::pair<double, double> direction = {dx / length, dy / length};
                if (have_direction_) last_turn_ = direction.first * direction_.first + direction.second * direction_.second;
                direction_ = direction;
                have_direction_ = true;
                last_length_ = length;
            }
        }
        previous_ = sample;
        previous_point_ = point;
        have_previous_ = true;
        trace.add(u, sample);
    }

    [[nodiscard]] double next_step(double du) const {
        double factor = 1.5;
        if (last_length_ > 0.0) factor = std::min(factor, 0.8 * max_length_ / last_length_);
        double turn = std::acos(std::min(std::max(last_turn_, -1.0), 1.0)) * 180.0 / M_PI;
        if (turn > 0.0) factor = std::min(factor, 0.8 * MAX_ANGLE_DEG / turn);
        return std::min(std::max(du * factor, du * 0.25), span_ * MAX_STEP);
    }

    void walk_inside(double& u, bool& inside) {
        double du = span_ / points_;
        double min_du = span_ * MIN_STEP;
        while (u < u_hi_ && !over_budget()) {
            double u_try = std::min(u + du, u_hi_);
            Sample sample;
            if (!sample_at(u_try, true, sample)) {
                if (skip_failures_) {
                    u = u_try;
                    du = std::min(du * 2.0, span_ * MAX_STEP);
                    continue;
                }
                if (du > min_du) {
                    du *= 0.5;
                    continue;
                }
                trace.add(u_try, Sample());
                u = u_try;
                inside = false;
                return;
            }
            if (!sample.inside) {
                leave_domain(u, u_try);
                u = u_try;
                inside = false;
                return;
            }
            std::pair<double, double> point = metric_.point(sample.x, sample.y);
            if (!step_ok(point, sample) && du > min_du) {
                du *= 0.5;
                continue;
            }
            accept(u_try, sample);
            u = u_try;
            du = next_step(du);
        }
        inside = true;
    }

    void leave_domain(double u_in, double u_out) {
        Sample edge;
        double u_edge = bisect_edge(u_in, u_out, edge);
        if (edge.valid) accept(u_edge, edge);
        trace.add(u_out, Sample());
    }

    double bisect_edge(double u_in, double u_out, Sample& found) {
        found = Sample();
        double tolerance = span_ * EDGE_TOLERANCE;
        while (std::fabs(u_out - u_in) > tolerance && !over_budget()) {
            double u_mid = 0.5 * (u_in + u_out);
            Sample sample;
            if (sample_at(u_mid, true, sample) && sample.inside) {
                u_in = u_mid;
                found = sample;
            } else {
                u_out = u_mid;
            }
        }
        return u_in;
    }

    void scan_outside(double& u, bool& inside) {
        double du = span_ * SCAN_STEP;
        while (u < u_hi_ && !over_budget()) {
            double u_try = std::min(u + du, u_hi_);
            Sample sample;
            if (sample_at(u_try, false, sample) && sample.inside) {
                enter_domain(u, u_try, sample);
                u = u_try;
                inside = true;
                return;
            }
            u = u_try;
            du = std::min(du * 2.0, span_ * SCAN_MAX);
        }
        inside = false;
    }

    void enter_domain(double u_out, double u_in, const Sample& sample) {
        have_previous_ = false;
        have_direction_ = false;
        Sample edge;
        double u_edge = bisect_edge(u_in, u_out, edge);
        if (edge.valid && u_edge != u_in) accept(u_edge, edge);
        accept(u_in, sample);
    }
};

Trace trace(const Evaluate& evaluate, const Metric& metric, int points, double u_lo = 0.0, double u_hi = 1.0, const Sample* start = nullptr,
            const Sample* end = nullptr, bool skip_failures = false) {
    Walk walk(evaluate, metric, points, u_lo, u_hi, start, end, skip_failures);
    walk.run();
    return walk.trace;
}

// ---------------------------------------------------------------------------
// Lines of constant property (CoolPlot's thermo/isolines.py)
// ---------------------------------------------------------------------------

const double EDGE_TOLERANCE_AXES = 1e-7;

bool inside_axes(const Metric* metric, double x, double y) {
    if (!metric) return true;
    std::pair<double, double> uv = metric->point(x, y);
    return uv.first >= -EDGE_TOLERANCE_AXES && uv.first <= 1.0 + EDGE_TOLERANCE_AXES && uv.second >= -EDGE_TOLERANCE_AXES
           && uv.second <= 1.0 + EDGE_TOLERANCE_AXES;
}

int two_phase_of(double Q) {
    return finite(Q) ? ((Q >= 0.0 && Q <= 1.0) ? 1 : 0) : -1;
}

double quality_or_nan(CoolProp::AbstractState& st) {
    try {
        return st.keyed_output(CoolProp::iQ);
    } catch (...) {
        return NaN;
    }
}

bool use_tracer(const Fluid& fluid, const std::string& solver) {
    if (solver != "auto" && solver != "tracer" && solver != "flash")
        throw CoolProp::ValueError("Unknown solver '" + solver + "', expected auto, tracer or flash.");
    return solver == "tracer" || (solver == "auto" && !fluid.is_pure());
}

/// The stepped quantity of a line mapped to a parameter u in 0..1 (log for pressure).
struct Sweep
{
    CoolProp::parameters key;
    double lo, hi;
    bool log;
    Sweep(CoolProp::parameters k, double low, double high)
      : key(k), lo(low), hi(high), log(default_scale(k) == Scale::Log && low > 0.0 && high > 0.0) {}
    [[nodiscard]] double value(double u) const {
        return log ? lo * std::pow(hi / lo, u) : lo + (hi - lo) * u;
    }
    [[nodiscard]] double parameter(double v) const {
        return log ? std::log(v / lo) / std::log(hi / lo) : (v - lo) / (hi - lo);
    }
};

/// Computes the states of one isoline: traced with a flash as fallback, or flashed.
class LineSolver
{
   public:
    LineSolver(Fluid& fluid, CoolProp::parameters iso_key, double value, CoolProp::parameters sweep_key, CoolProp::parameters xkey,
               CoolProp::parameters ykey, const Metric* metric, const std::string& solver)
      : fluid_(fluid), value_(value), xkey_(xkey), ykey_(ykey), metric_(metric) {
        double o1 = 0.0, o2 = 0.0;
        pair_ = CoolProp::generate_update_pair(iso_key, 0.0, sweep_key, 1.0, o1, o2);
        constant_first_ = (o1 == 0.0);
        CoolProp::parameters first = constant_first_ ? iso_key : sweep_key, second = constant_first_ ? sweep_key : iso_key;
        if (use_tracer(fluid, solver) && IsoLineTracer::supports(first, second)) {
            try {
                tracer_ = std::make_unique<IsoLineTracer>(fluid, first, second, iso_key);
            } catch (...) {
                tracer_.reset();
            }
        }
        method_ = tracer_ ? "tracer" : "flash";
    }

    [[nodiscard]] const std::string& method() const {
        return method_;
    }
    [[nodiscard]] int n_fallback() const {
        return n_fallback_;
    }
    IsoLineTracer* tracer() {
        return tracer_.get();
    }

    bool sample(double w, bool careful, Sample& out) {
        double v1 = constant_first_ ? value_ : w, v2 = constant_first_ ? w : value_;
        if (tracer_) {
            bool ok = false;
            try {
                tracer_->update(v1, v2);
                ok = sample_of(tracer_->state(), w, "tracer", out);
            } catch (...) {
                ok = false;
            }
            if (ok) return true;
            if (!careful) return false;
            ++n_fallback_;
        }
        try {
            CoolProp::AbstractState& st = fluid_.update_pair(pair_, v1, v2);
            if (tracer_) tracer_->set_guess(st);
            return sample_of(st, w, "flash", out);
        } catch (...) {
            return false;
        }
    }

   private:
    // A short-lived helper of one calc_isolines call: it refers to the plot's fluid by design
    Fluid& fluid_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
    double value_;
    CoolProp::parameters xkey_, ykey_;
    const Metric* metric_;
    CoolProp::input_pairs pair_;
    bool constant_first_;
    std::unique_ptr<IsoLineTracer> tracer_;
    std::string method_;
    int n_fallback_ = 0;

    bool sample_of(CoolProp::AbstractState& st, double w, const char* source, Sample& out) {
        double x = st.keyed_output(xkey_), y = st.keyed_output(ykey_), T = st.T(), p = st.p();
        if (!(finite(x) && finite(y) && finite(T) && finite(p))) return false;
        out.valid = true;
        out.x = x;
        out.y = y;
        out.inside = inside_axes(metric_, x, y);
        out.Q = quality_or_nan(st);
        out.two_phase = two_phase_of(out.Q);
        out.w = w;
        out.T = T;
        out.p = p;
        out.rho = st.rhomass();
        out.source = source;
        return true;
    }
};

/// The saturation state at T or p = value and the given quality, false if it cannot be computed.
bool saturation_sample(Fluid& fluid, CoolProp::parameters key, double value, double quality, CoolProp::parameters xkey, CoolProp::parameters ykey,
                       const Metric* metric, CoolProp::parameters sweep_key, Sample& out) {
    try {
        CoolProp::AbstractState& st = fluid.update(key, value, CoolProp::iQ, quality);
        if (!fluid.is_pure() && !phases_are_distinct(st)) return false;
        double x = st.keyed_output(xkey), y = st.keyed_output(ykey), T = st.T(), p = st.p();
        if (!(finite(x) && finite(y) && finite(T) && finite(p))) return false;
        out.valid = true;
        out.x = x;
        out.y = y;
        out.inside = inside_axes(metric, x, y);
        out.two_phase = 1;
        out.w = sweep_key == CoolProp::iT ? T : p;
        out.T = T;
        out.p = p;
        out.Q = quality;
        out.rho = st.rhomass();
        out.source = "flash";
        return true;
    } catch (...) {
        return false;
    }
}

/// (iso value, x, y, T, p, rho) of the saturation state at T or p = sat_value, false if none.
bool saturation_state(Fluid& fluid, CoolProp::parameters sat_key, double sat_value, double quality, CoolProp::parameters xkey,
                      CoolProp::parameters ykey, CoolProp::parameters iso_key, std::vector<double>& out) {
    try {
        CoolProp::AbstractState& st = fluid.update(sat_key, sat_value, CoolProp::iQ, quality);
        if (!fluid.is_pure() && !phases_are_distinct(st)) return false;
        out = {st.keyed_output(iso_key), st.keyed_output(xkey), st.keyed_output(ykey), st.T(), st.p(), st.rhomass()};
        for (double v : out)
            if (!finite(v)) return false;
        return true;
    } catch (...) {
        return false;
    }
}

struct Crossing
{
    double w = NaN, quality = NaN;
    Sample sample;
};

Sample crossing_sample(const std::vector<double>& found, double w, double quality, const Metric* metric, const char* source) {
    Sample s;
    s.valid = true;
    s.x = found[1];
    s.y = found[2];
    s.inside = inside_axes(metric, s.x, s.y);
    s.two_phase = -1;
    s.w = w;
    s.T = found[3];
    s.p = found[4];
    s.Q = quality;
    s.rho = found[5];
    s.source = source;
    return s;
}

/// Whether two sweep values are the same up to rounding (relative 1e-9).
bool same_sweep_value(double a, double b) {
    return std::fabs(a - b) <= 1e-9 * std::max(std::fabs(a), std::fabs(b));
}

/// Whether the dome end at sweep value w_a and quality q_a comes before the one at (w_b, q_b).
/// The bubble and dew states of a pure fluid have the same sweep value, but computed apart
/// they can differ in the last bits; a bit-exact compare would then order them by that noise.
/// Equal up to rounding, the quality decides: liquid first when the sweep is in temperature
/// (liquid_first, an isobar), vapour first when it is in pressure (an isotherm).
bool dome_end_comes_first(double w_a, double q_a, double w_b, double q_b, bool liquid_first) {
    if (!same_sweep_value(w_a, w_b)) return w_a < w_b;
    return liquid_first ? q_a < q_b : q_a > q_b;
}

/// The sweep value in [a, b] where the saturation state has iso_key == value: regula falsi (Illinois).
bool solve_crossing(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, CoolProp::parameters iso_key, double value,
                    CoolProp::parameters sweep_key, double a, double b, double quality, const Metric* metric, Crossing& out) {
    std::vector<double> end_a, end_b;
    if (!saturation_state(fluid, sweep_key, a, quality, xkey, ykey, iso_key, end_a)) return false;
    if (!saturation_state(fluid, sweep_key, b, quality, xkey, ykey, iso_key, end_b)) return false;
    double fa = end_a[0] - value, fb = end_b[0] - value;
    if (fa == 0.0) {
        out = {a, quality, crossing_sample(end_a, a, quality, metric, "flash")};
        return true;
    }
    if (fb == 0.0) {
        out = {b, quality, crossing_sample(end_b, b, quality, metric, "flash")};
        return true;
    }
    if (fa * fb > 0.0) return false;
    for (int i = 0; i < 60; ++i) {
        double c = b - fb * (b - a) / (fb - fa);
        std::vector<double> state;
        if (!saturation_state(fluid, sweep_key, c, quality, xkey, ykey, iso_key, state)) return false;
        double fc = state[0] - value;
        out = {c, quality, crossing_sample(state, c, quality, metric, "flash")};
        if (fc == 0.0 || std::fabs(b - a) <= 1e-12 * std::fabs(c)) return true;
        if (fc * fb < 0.0) {
            a = b;
            fa = fb;
        } else {
            fa *= 0.5;
        }
        b = c;
        fb = fc;
        if (std::fabs(fc) <= 1e-10 * std::max(std::fabs(value), 1.0)) return true;
    }
    // Not converged in 60 steps: the last iterate is not a crossing to put on
    // the line.  A mixture line then takes the phase-envelope crossing; a
    // pure fluid's line goes on without a point at the dome.
    return false;
}

/// Where an isoline crosses the phase envelope of a mixture, interpolated linearly.
std::vector<Crossing> envelope_crossings(Envelope& envelope, CoolProp::parameters xkey, CoolProp::parameters ykey, CoolProp::parameters iso_key,
                                         double value, CoolProp::parameters sweep_key, const Metric* metric, const Range* sweep_range) {
    const std::vector<double>& k = envelope.values(iso_key);
    std::vector<const std::vector<double>*> columns = {&envelope.values(sweep_key),    &envelope.values(xkey),
                                                       &envelope.values(ykey),         &envelope.values(CoolProp::iT),
                                                       &envelope.values(CoolProp::iP), &envelope.values(CoolProp::iDmass)};
    std::vector<Crossing> found;
    for (size_t i = 0; i + 1 < k.size(); ++i) {
        double a = k[i] - value, b = k[i + 1] - value;
        if (!(finite(a) && finite(b)) || a * b > 0.0 || a == b) continue;
        if (b == 0.0 && i + 1 < k.size() - 1) continue;
        double t = a / (a - b);
        std::vector<double> point(6);
        bool ok = true;
        for (size_t c = 0; c < 6; ++c) {
            const std::vector<double>& v = *columns[c];
            point[c] = v[i] + t * (v[i + 1] - v[i]);
            ok = ok && finite(point[c]);
        }
        if (!ok) continue;
        if (sweep_range && !(sweep_range->min <= point[0] && point[0] <= sweep_range->max)) continue;
        double quality = t < 0.5 ? envelope.Q[i] : envelope.Q[i + 1];
        std::vector<double> as_state = {value, point[1], point[2], point[3], point[4], point[5]};
        found.push_back({point[0], quality, crossing_sample(as_state, point[0], quality, metric, "envelope")});
    }
    return found;
}

/// Why a line has fewer than two points, from the counts of its evaluations; empty for two or more points.
/// A state computed outside the axes is not a point of the line, so it is told apart from a failure.
std::string short_line_reason(size_t n, int n_evaluated, int n_ok, int n_outside_ok) {
    if (n >= 2) return "";
    std::ostringstream text;
    if (n_evaluated > 0 && n_ok == 0)
        text << "none of the " << n_evaluated << " states of the line could be computed";
    else if (n_outside_ok > 0)
        text << n << " point(s) of the line lie in the area the axes show";
    else
        text << "only " << n << " point(s) of the line could be computed";
    return text.str();
}

/// Assembles one isoline from stepped parts.
class Line
{
   public:
    Line(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, CoolProp::parameters iso_key, double value, const Sweep& sweep,
         const Metric& metric, int points, LineSolver& solver)
      : fluid_(fluid), xkey_(xkey), ykey_(ykey), iso_key_(iso_key), value_(value), sweep_(sweep), metric_(metric), points_(points), solver_(solver) {}

    std::vector<Sample> samples;
    int n_inserted = 0, n_evaluated = 0, n_outside = 0, n_failed = 0;
    int n_ok = 0;          // evaluations that gave a state, see Trace::n_ok
    int n_outside_ok = 0;  // of those, the states outside the domain
    std::string stopped;

    void through_the_dome() {
        std::vector<std::pair<double, Sample>> ends;  // (quality, sample) in order of w
        if (!dome_interval(ends)) {
            part(0.0, 1.0, nullptr, nullptr);
            return;
        }
        double w_a = ends[0].second.w, w_b = ends[1].second.w;
        if (w_b < sweep_.lo || w_a > sweep_.hi) {
            part(0.0, 1.0, nullptr, nullptr);
            return;
        }
        double u_a = std::min(std::max(sweep_.parameter(w_a), 0.0), 1.0), u_b = std::min(std::max(sweep_.parameter(w_b), 0.0), 1.0);
        Sample sample_a = in_sweep_range(ends[0].second), sample_b = in_sweep_range(ends[1].second);
        part(0.0, u_a, nullptr, &sample_a);
        quality_branch(ends[0].first, sample_a, ends[1].first, sample_b);
        part(u_b, 1.0, &sample_b, nullptr);
    }

    void across_the_dome() {
        part(0.0, 1.0, nullptr, nullptr);
        if (!fluid_.has_critical_point()) return;
        Envelope* envelope = fluid_.envelope();
        if (envelope) {
            Range range = {sweep_.lo, sweep_.hi};
            merge(envelope_crossings(*envelope, xkey_, ykey_, iso_key_, value_, sweep_.key, &metric_, &range));
        } else {
            insert_pure_crossings();
        }
    }

    [[nodiscard]] Isoline curve(const std::string& method, int n_fallback) const {
        Isoline out;
        out.value = value_;
        size_t n = samples.size();
        while (n > 0 && !samples[n - 1].valid)
            --n;
        for (size_t i = 0; i < n; ++i) {
            out.x.push_back(samples[i].valid ? samples[i].x : NaN);
            out.y.push_back(samples[i].valid ? samples[i].y : NaN);
        }
        out.n_failed = n_failed;
        out.n_fallback = n_fallback;
        out.n_inserted = n_inserted;
        out.method = method;
        out.n_evaluated = n_evaluated;
        out.n_outside = n_outside;
        out.reason = stopped;
        if (out.reason.empty()) out.reason = short_line_reason(n, n_evaluated, n_ok, n_outside_ok);
        return out;
    }

   private:
    // A short-lived helper of one calc_isolines call: it refers to the plot's fluid and solver by design
    Fluid& fluid_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
    CoolProp::parameters xkey_, ykey_, iso_key_;
    double value_;
    Sweep sweep_;
    Metric metric_;
    int points_;
    LineSolver& solver_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

    void count(const Trace& found) {
        n_evaluated += found.n_evaluated;
        n_outside += found.n_outside;
        n_failed += found.n_failed;
        n_ok += found.n_ok;
        n_outside_ok += found.n_outside_ok;
        if (!found.stopped.empty() && found.stopped.rfind("the end of the range", 0) != 0) stopped = found.stopped;
    }

    void extend(const std::vector<Sample>& found) {
        for (const Sample& s : found) {
            if (s.valid && !samples.empty() && samples.back().valid && samples.back().source == s.source && samples.back().w == s.w
                && samples.back().x == s.x && samples.back().y == s.y)
                continue;  // the end of one part is the start of the next
            if (!s.valid && (samples.empty() || !samples.back().valid)) continue;
            samples.push_back(s);
        }
    }

    void part(double u_lo, double u_hi, const Sample* start, const Sample* end) {
        if (u_hi <= u_lo) return;
        Trace found = trace([this](double u, bool careful, Sample& out) { return solver_.sample(sweep_.value(u), careful, out); }, metric_, points_,
                            u_lo, u_hi, start, end, false);
        count(found);
        extend(found.samples);
    }

    [[nodiscard]] Sample in_sweep_range(const Sample& s) const {
        Sample copy = s;
        if (!(sweep_.lo <= s.w && s.w <= sweep_.hi)) copy.inside = false;
        return copy;
    }

    bool dome_interval(std::vector<std::pair<double, Sample>>& ends) {
        if (!fluid_.has_critical_point()) return false;
        for (double q : {0.0, 1.0}) {
            Sample s;
            ++n_evaluated;
            if (saturation_sample(fluid_, iso_key_, value_, q, xkey_, ykey_, &metric_, sweep_.key, s)) {
                ++n_ok;
                if (!s.inside) {
                    ++n_outside;
                    ++n_outside_ok;
                }
                ends.emplace_back(q, s);
            }
        }
        if (ends.size() < 2) {
            ends.clear();
            Envelope* envelope = fluid_.envelope();
            if (!envelope) return false;
            std::map<double, Crossing> found;
            for (const Crossing& c : envelope_crossings(*envelope, xkey_, ykey_, iso_key_, value_, sweep_.key, &metric_, nullptr))
                if (found.count(c.quality) == 0) found[c.quality] = c;
            for (double q : {0.0, 1.0})
                if (found.count(q) != 0) ends.emplace_back(q, found[q].sample);
        }
        if (ends.size() != 2) return false;
        // Two ends, so a swap orders them (see dome_end_comes_first)
        if (!dome_end_comes_first(ends[0].second.w, ends[0].first, ends[1].second.w, ends[1].first, sweep_.key == CoolProp::iT))
            std::swap(ends[0], ends[1]);
        return true;
    }

    void quality_branch(double q_a, const Sample& sample_a, double q_b, const Sample& sample_b) {
        Evaluate evaluate = [&](double u, bool, Sample& out) {
            double q = q_a + (q_b - q_a) * u;
            if (!saturation_sample(fluid_, iso_key_, value_, q, xkey_, ykey_, &metric_, sweep_.key, out)) return false;
            out = in_sweep_range(out);
            return true;
        };
        Trace found = trace(evaluate, metric_, points_, 0.0, 1.0, &sample_a, &sample_b, true);
        count(found);
        for (const Sample& s : found.samples)
            if (s.valid && ((s.w == sample_a.w && s.Q == sample_a.Q) || (s.w == sample_b.w && s.Q == sample_b.Q))) ++n_inserted;
        extend(found.samples);
    }

    void insert_pure_crossings() {
        std::vector<Sample> merged;
        for (size_t i = 0; i < samples.size(); ++i) {
            const Sample& s = samples[i];
            merged.push_back(s);
            if (i + 1 >= samples.size()) continue;
            const Sample& after = samples[i + 1];
            if (!s.valid || !after.valid || s.two_phase < 0 || after.two_phase < 0 || s.two_phase == after.two_phase) continue;
            double q_inside = s.two_phase != 0 ? s.Q : after.Q;
            std::vector<double> order = q_inside <= 0.5 ? std::vector<double>{0.0, 1.0} : std::vector<double>{1.0, 0.0};
            for (double q : order) {
                Crossing crossing;
                if (solve_crossing(fluid_, xkey_, ykey_, iso_key_, value_, sweep_.key, s.w, after.w, q, &metric_, crossing)) {
                    merged.push_back(crossing.sample);
                    ++n_inserted;
                    break;
                }
            }
        }
        samples = merged;
    }

    void merge(std::vector<Crossing> crossings) {
        std::sort(crossings.begin(), crossings.end(), [](const Crossing& a, const Crossing& b) { return a.w < b.w; });
        for (Crossing& c : crossings) {
            size_t position = samples.size();
            for (size_t i = 0; i < samples.size(); ++i)
                if (samples[i].valid && samples[i].w > c.w) {
                    position = i;
                    break;
                }
            const Sample* before = position > 0 ? &samples[position - 1] : nullptr;
            const Sample* after = position < samples.size() ? &samples[position] : nullptr;
            Sample inserted = c.sample;
            if (before && after && before->valid && after->valid) {
                Crossing exact;
                if (solve_crossing(fluid_, xkey_, ykey_, iso_key_, value_, sweep_.key, before->w, after->w, c.quality, &metric_, exact))
                    inserted = exact.sample;
            }
            samples.insert(samples.begin() + static_cast<long>(position), inserted);
            ++n_inserted;
        }
    }
};

// ---------------------------------------------------------------------------
// Lines of constant quality
// ---------------------------------------------------------------------------

Isoline envelope_line(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, double quality, bool& ok) {
    ok = false;
    Isoline out;
    Envelope* envelope = fluid.envelope();
    if (!envelope) return out;
    std::vector<size_t> dew, bubble;
    for (size_t i = 0; i < envelope->Q.size(); ++i) {
        if (envelope->Q[i] == 1.0) dew.push_back(i);
        if (envelope->Q[i] == 0.0) bubble.push_back(i);
    }
    if (dew.size() < 2 || bubble.size() < 2) return out;
    std::vector<size_t> rows;
    if (quality == 1.0)
        rows = dew;
    else if (dew.back() + 1 == bubble.front()) {
        rows.push_back(dew.back());
        rows.insert(rows.end(), bubble.begin(), bubble.end());
    } else if (bubble.back() + 1 == dew.front()) {
        rows = bubble;
        rows.push_back(dew.front());
    } else
        rows = bubble;
    const std::vector<double>& xs = envelope->values(xkey);
    const std::vector<double>& ys = envelope->values(ykey);
    int n_failed = 0;
    for (size_t i : rows) {
        out.x.push_back(xs[i]);
        out.y.push_back(ys[i]);
        if (!(finite(xs[i]) && finite(ys[i]))) ++n_failed;
    }
    if (static_cast<size_t>(n_failed) > rows.size() / 2) return {};
    out.value = quality;
    out.n_failed = n_failed;
    out.method = "envelope";
    ok = true;
    return out;
}

Isoline mixture_quality_line(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, double quality, int points, double margin,
                             Envelope& envelope) {
    double T_lo = fluid.saturation_range(CoolProp::iT, margin).min;
    double T_hi = envelope.top(CoolProp::iT);
    Metric metric{};
    for (int axis = 0; axis < 2; ++axis) {
        const std::vector<double>& v = envelope.values(axis == 0 ? xkey : ykey);
        double lo = NaN, hi = NaN;
        for (double value : v) {
            if (!finite(value)) continue;
            lo = finite(lo) ? std::min(lo, value) : value;
            hi = finite(hi) ? std::max(hi, value) : value;
        }
        if (!(lo < hi)) {
            lo = (finite(lo) ? lo : 0.0) - 1.0;
            hi = lo + 2.0;
        }
        if (axis == 0)
            metric.x = {lo, hi};
        else
            metric.y = {lo, hi};
    }
    metric.x_log = default_scale(xkey) == Scale::Log;
    metric.y_log = default_scale(ykey) == Scale::Log;
    Sample top;
    top.valid = true;
    top.x = envelope.top(xkey);
    top.y = envelope.top(ykey);
    top.inside = true;
    top.two_phase = 1;
    top.w = T_hi;
    top.T = T_hi;
    top.p = envelope.top(CoolProp::iP);
    top.Q = quality;
    top.rho = envelope.top(CoolProp::iDmass);
    top.source = "envelope";
    Evaluate evaluate = [&](double u, bool, Sample& out) {
        return saturation_sample(fluid, CoolProp::iT, T_lo + (T_hi - T_lo) * u, quality, xkey, ykey, nullptr, CoolProp::iT, out);
    };
    Trace found = trace(evaluate, metric, points, 0.0, 1.0, nullptr, &top, true);
    Isoline out;
    out.value = quality;
    for (const Sample& s : found.samples) {
        out.x.push_back(s.valid ? s.x : NaN);
        out.y.push_back(s.valid ? s.y : NaN);
    }
    out.n_failed = found.n_failed;
    out.n_evaluated = found.n_evaluated;
    out.reason = found.stopped.rfind("the end of the range", 0) == 0 ? "" : found.stopped;
    // The top of the envelope is always on the line, so count the drawable points
    if (out.reason.empty()) out.reason = short_line_reason(out.n_valid(), found.n_evaluated, found.n_ok, found.n_outside_ok);
    return out;
}

Isoline compute_quality_line(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, double quality, int points, double margin,
                             const std::string& spacing) {
    if (!(quality >= 0.0 && quality <= 1.0)) throw CoolProp::ValueError("The vapour quality has to be between 0 and 1.");
    if (!fluid.has_critical_point())
        throw CoolProp::ValueError(fluid.name()
                                   + " has no two-phase region (CoolProp provides no critical point for it), so lines of constant "
                                     "quality cannot be drawn.");
    if (!fluid.is_pure()) {
        if (quality == 0.0 || quality == 1.0) {
            bool ok = false;
            Isoline line = envelope_line(fluid, xkey, ykey, quality, ok);
            if (ok) return line;
        }
        Envelope* envelope = fluid.envelope();
        if (envelope && envelope->closed) return mixture_quality_line(fluid, xkey, ykey, quality, points, margin, *envelope);
    }
    Range T_range = fluid.saturation_range(CoolProp::iT, margin);
    std::vector<double> u = linspace(0.0, 1.0, static_cast<size_t>(std::max(points, 2)));
    if (spacing == "critical")
        for (double& v : u)
            v = 1.0 - (1.0 - v) * (1.0 - v);
    else if (spacing != "linear")
        throw CoolProp::ValueError("Unknown saturation spacing '" + spacing + "'.");
    Isoline out;
    out.value = quality;
    for (double v : u) {
        double T = T_range.min + (T_range.max - T_range.min) * v;
        try {
            CoolProp::AbstractState& st = fluid.update(CoolProp::iT, T, CoolProp::iQ, quality);
            if (!fluid.is_pure() && !phases_are_distinct(st)) throw CoolProp::ValueError("trivial solution");
            out.x.push_back(st.keyed_output(xkey));
            out.y.push_back(st.keyed_output(ykey));
        } catch (...) {
            out.x.push_back(NaN);
            out.y.push_back(NaN);
            ++out.n_failed;
        }
    }
    out.n_evaluated = static_cast<int>(u.size());
    // Judged before the critical point is appended: a line of that one point says nothing
    out.reason = short_line_reason(out.n_valid(), out.n_evaluated, out.n_evaluated - out.n_failed, 0);
    if (fluid.is_pure()) {
        const CriticalPoint& crit = fluid.critical();
        double xc = crit.get(xkey), yc = crit.get(ykey);
        if (finite(xc) && finite(yc)) {
            out.x.push_back(xc);
            out.y.push_back(yc);
        }
    }
    return out;
}

/// The quantity stepped along a line of constant key: pressure for T, h, s and u, temperature for p, rho and Q.
CoolProp::parameters sweep_of(CoolProp::parameters key) {
    switch (key) {
        case CoolProp::iT:
        case CoolProp::iHmass:
        case CoolProp::iSmass:
        case CoolProp::iUmass:
            return CoolProp::iP;
        case CoolProp::iP:
        case CoolProp::iDmass:
        case CoolProp::iQ:
            return CoolProp::iT;
        default:
            return CoolProp::iundefined_parameter;
    }
}

Isoline compute_isoline(Fluid& fluid, CoolProp::parameters xkey, CoolProp::parameters ykey, CoolProp::parameters iso_key, double value,
                        const std::vector<double>& tp_box, const Metric& metric, const CalcSettings& settings) {
    CoolProp::parameters sweep_key = sweep_of(iso_key);
    if (sweep_key == CoolProp::iundefined_parameter || iso_key == xkey || iso_key == ykey)
        throw CoolProp::ValueError("Lines of constant " + CoolProp::get_parameter_information(iso_key, "short")
                                   + " cannot be drawn in this diagram.");
    Sweep sweep(sweep_key, sweep_key == CoolProp::iT ? tp_box[0] : tp_box[2], sweep_key == CoolProp::iT ? tp_box[1] : tp_box[3]);
    LineSolver solver(fluid, iso_key, value, sweep_key, xkey, ykey, &metric, settings.solver);
    Line line(fluid, xkey, ykey, iso_key, value, sweep, metric, settings.points, solver);
    if (iso_key == CoolProp::iT || iso_key == CoolProp::iP)
        line.through_the_dome();
    else
        line.across_the_dome();
    return line.curve(solver.method(), solver.n_fallback());
}

}  // namespace Detail

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

size_t Isoline::n_valid() const {
    size_t n = 0;
    for (size_t i = 0; i < x.size(); ++i)
        if (std::isfinite(x[i]) && std::isfinite(y[i])) ++n;
    return n;
}

size_t Isoline::n_gaps() const {
    size_t first = x.size(), last = 0, valid = 0;
    for (size_t i = 0; i < x.size(); ++i)
        if (std::isfinite(x[i]) && std::isfinite(y[i])) {
            first = std::min(first, i);
            last = i;
            ++valid;
        }
    if (valid < 2) return 0;
    return last - first + 1 - valid;
}

std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count) {
    if (count <= 1) return {range.min};
    auto n = static_cast<size_t>(count);
    if (scale == Scale::Log) return logspace(range.min, range.max, n);
    return linspace(range.min, range.max, n);
}

std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count) {
    return generate_values_in_range(Detail::default_scale(type), range, count);
}

PropertyPlot::PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, TPLimits tp_limits)
  : xkey_(xkey), ykey_(ykey) {
    init(fluid_name);
    const double HI = 2.25, LO = 1.01;
    switch (tp_limits) {
        case TPLimits::None:
            finish_limits(Detail::NaN, Detail::NaN, Detail::NaN, Detail::NaN, false);
            break;
        case TPLimits::Def:
            finish_limits(LO, HI, LO, HI, true);
            break;
        case TPLimits::Achp:
            finish_limits(173.15, 493.15, 0.25e5, HI, true);
            break;
        case TPLimits::Orc:
            finish_limits(273.15, 673.15, 0.25e5, HI, true);
            break;
    }
}

PropertyPlot::PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, double T_min_K, double T_max_K,
                           double p_min_Pa, double p_max_Pa)
  : xkey_(xkey), ykey_(ykey) {
    // An infinite limit means "the limit of the equation of state"; NaN is a
    // mistake in the caller and would silently be read as the same
    for (double limit : {T_min_K, T_max_K, p_min_Pa, p_max_Pa})
        if (std::isnan(limit))
            throw CoolProp::ValueError("A limit of the plotting region is not a number; give a value, or an infinite one for the limit of "
                                       "the equation of state.");
    init(fluid_name);
    finish_limits(T_min_K, T_max_K, p_min_Pa, p_max_Pa, false);
}

PropertyPlot::~PropertyPlot() = default;
PropertyPlot::PropertyPlot(PropertyPlot&&) noexcept = default;
PropertyPlot& PropertyPlot::operator=(PropertyPlot&&) noexcept = default;

void PropertyPlot::init(const std::string& fluid_name) {
    fluid_ = std::make_unique<Detail::Fluid>(fluid_name);
    xaxis.scale = Detail::default_scale(xkey_);
    yaxis.scale = Detail::default_scale(ykey_);
    double o1 = 0.0, o2 = 0.0;
    axis_pair_ = CoolProp::generate_update_pair(xkey_, 0.0, ykey_, 1.0, o1, o2);
    if (axis_pair_ == CoolProp::INPUT_PAIR_INVALID)
        throw CoolProp::ValueError("CoolProp has no flash for the axes " + CoolProp::get_parameter_information(ykey_, "short") + " over "
                                   + CoolProp::get_parameter_information(xkey_, "short") + ".");
    swap_axis_inputs_for_update_ = (o1 == 1.0);
}

/// Resolve the limits: a value below 10 with ``relative`` is a factor of the
/// saturation range (CoolPlot's Relative bounds), NaN the limit of the
/// equation of state; then clamp to the equation of state.
void PropertyPlot::finish_limits(double T_min, double T_max, double p_min, double p_max, bool relative) {
    const double ID_FACTOR = 10.0;
    auto absolute = [&](double bound, CoolProp::parameters key, bool upper, double fallback) {
        if (!std::isfinite(bound)) return fallback;
        if (relative && bound < ID_FACTOR) {
            Range sat = fluid_->saturation_range(key, 1e-7);
            return bound * (upper ? sat.max : sat.min);
        }
        return bound;
    };
    double T_lo = absolute(T_min, CoolProp::iT, false, 0.0), T_hi = absolute(T_max, CoolProp::iT, true, std::numeric_limits<double>::infinity());
    double p_lo = absolute(p_min, CoolProp::iP, false, 0.0), p_hi = absolute(p_max, CoolProp::iP, true, std::numeric_limits<double>::infinity());
    const double eos_T_min = fluid_->T_min_K();
    const double eos_T_max = fluid_->T_max_K();
    const double eos_p_min = fluid_->p_min_Pa();
    const double eos_p_max = fluid_->p_max_Pa();
    if (std::isfinite(eos_T_min)) T_lo = std::max(T_lo, eos_T_min);
    if (std::isfinite(eos_T_max)) T_hi = std::min(T_hi, eos_T_max);
    if (std::isfinite(eos_p_min)) p_lo = std::max(p_lo, eos_p_min);
    if (std::isfinite(eos_p_max)) p_hi = std::min(p_hi, eos_p_max);
    if (!(T_lo < T_hi && p_lo < p_hi)) {
        std::ostringstream text;
        text << "Empty plotting region: T = " << T_lo << ".." << T_hi << " K, p = " << p_lo << ".." << p_hi
             << " Pa after clamping to the equation of state.";
        throw CoolProp::ValueError(text.str());
    }
    T_min_K_ = T_lo;
    T_max_K_ = T_hi;
    p_min_Pa_ = p_lo;
    p_max_Pa_ = p_hi;
    // Axis is a public union (min and max alias range), kept for compatibility with existing callers
    xaxis.range = axis_range(xkey_);  // NOLINT(cppcoreguidelines-pro-type-union-access)
    yaxis.range = axis_range(ykey_);  // NOLINT(cppcoreguidelines-pro-type-union-access)
}

std::vector<double> PropertyPlot::tp_box() const {
    return {T_min_K_, T_max_K_, p_min_Pa_, p_max_Pa_};
}

bool PropertyPlot::is_mixture() const {
    return !fluid_->is_pure();
}

/// Smallest and largest value of a quantity on the corners of the T-p box.
Range PropertyPlot::axis_range(CoolProp::parameters key) const {
    if (key == CoolProp::iT) return {T_min_K_, T_max_K_};
    if (key == CoolProp::iP) return {p_min_Pa_, p_max_Pa_};
    if (key == CoolProp::iQ) return {0.0, 1.0};
    std::vector<double> values;
    for (double T : {T_min_K_, T_max_K_})
        for (double p : {p_min_Pa_, p_max_Pa_}) {
            try {
                fluid_->update(CoolProp::iP, p, CoolProp::iT, T);
                double v = fluid_->get(key);
                if (std::isfinite(v)) values.push_back(v);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // a corner without a state does not count
            }
        }
    if (values.size() < 2) throw CoolProp::ValueError("Could not determine the range of " + CoolProp::get_parameter_information(key, "short") + ".");
    return {*std::min_element(values.begin(), values.end()), *std::max_element(values.begin(), values.end())};
}

Range PropertyPlot::isoline_range(CoolProp::parameters key) const {
    std::lock_guard<std::mutex> guard(*lock_);
    return axis_range(key);
}

std::vector<CoolProp::parameters> PropertyPlot::supported_isoline_keys() const {
    std::vector<CoolProp::parameters> keys;
    for (CoolProp::parameters key :
         {CoolProp::iT, CoolProp::iP, CoolProp::iHmass, CoolProp::iSmass, CoolProp::iDmass, CoolProp::iUmass, CoolProp::iQ})
        if (key != xkey_ && key != ykey_) keys.push_back(key);
    return keys;
}

CoolProp::parameters PropertyPlot::sweep_key(CoolProp::parameters key) const {
    if (key == xkey_ || key == ykey_) return CoolProp::iundefined_parameter;
    return Detail::sweep_of(key);
}

Isolines PropertyPlot::calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const {
    CalcSettings settings;
    settings.points = points;
    settings.saturation_points = points;
    return calc_isolines(key, values, settings);
}

Isolines PropertyPlot::calc_isolines(CoolProp::parameters key, const std::vector<double>& values, const CalcSettings& given) const {
    // Invalid arguments throw before any line is computed
    if (given.solver != "auto" && given.solver != "tracer" && given.solver != "flash")
        throw CoolProp::ValueError("Unknown solver '" + given.solver + "'; use 'auto', 'tracer' or 'flash'.");
    if (given.saturation_spacing != "critical" && given.saturation_spacing != "linear")
        throw CoolProp::ValueError("Unknown saturation spacing '" + given.saturation_spacing + "'; use 'critical' or 'linear'.");
    if (!(given.critical_margin >= 0.0 && given.critical_margin < 0.5))
        throw CoolProp::ValueError("The critical margin has to be at least 0 and below 0.5.");
    std::lock_guard<std::mutex> guard(*lock_);
    if (key == CoolProp::iQ) {
        for (double value : values)
            if (!(value >= 0.0 && value <= 1.0)) throw CoolProp::ValueError("The vapour quality has to be between 0 and 1.");
        if (!fluid_->has_critical_point())
            throw CoolProp::ValueError(fluid_->name()
                                       + " has no two-phase region (CoolProp provides no critical point for it), so lines of constant "
                                         "quality cannot be drawn.");
        // For a pure fluid the margin is taken off both ends of the saturation
        // range (from the triple point up and from the critical point down): it
        // must leave some range, or the lines would run backwards. A mixture's
        // lines follow its phase envelope up to the envelope's hottest point,
        // which the margin does not bound, so they are not checked here.
        if (fluid_->is_pure()) {
            Range T_sat = fluid_->saturation_range(CoolProp::iT, given.critical_margin);
            if (!(T_sat.min < T_sat.max)) throw CoolProp::ValueError("The critical margin leaves no saturation range for " + fluid_->name() + ".");
        }
    } else if (sweep_key(key) == CoolProp::iundefined_parameter) {
        throw CoolProp::ValueError("Lines of constant " + CoolProp::get_parameter_information(key, "short") + " cannot be drawn in this diagram.");
    }
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
    Detail::Metric metric = {xaxis.range, yaxis.range, xaxis.scale == Scale::Log, yaxis.scale == Scale::Log};
    // At least two points per line, at most a million: the stepper's budget
    // and the quality grid are sized from it
    const int MAX_POINTS = 1000000;
    CalcSettings settings = given;
    settings.points = std::min(std::max(settings.points, 2), MAX_POINTS);
    settings.saturation_points = std::min(std::max(settings.saturation_points, 2), MAX_POINTS);
    Isolines lines;
    for (double value : values) {
        // A line that fails comes back with its reason, so that one bad value
        // does not cost the other lines of the call
        try {
            if (key == CoolProp::iQ)
                lines.push_back(Detail::compute_quality_line(*fluid_, xkey_, ykey_, value, settings.saturation_points, settings.critical_margin,
                                                             settings.saturation_spacing));
            else
                lines.push_back(Detail::compute_isoline(*fluid_, xkey_, ykey_, key, value, tp_box(), metric, settings));
        } catch (const std::exception& e) {
            Isoline failed;
            failed.value = value;
            failed.reason = e.what();
            lines.push_back(failed);
        }
    }
    return lines;
}

double PropertyPlot::value_at(CoolProp::parameters key, double xvalue, double yvalue, CoolProp::phases phase) const {
    if (key == xkey_) return xvalue;
    if (key == ykey_) return yvalue;
    std::lock_guard<std::mutex> guard(*lock_);
    try {
        if (swap_axis_inputs_for_update_) std::swap(xvalue, yvalue);
        CoolProp::AbstractState& st = *fluid_->state();
        // Backends without a phase to impose (incompressibles) throw on specify_phase
        if (phase != CoolProp::iphase_not_imposed) st.specify_phase(phase);
        st.update(axis_pair_, xvalue, yvalue);
        double v = st.keyed_output(key);
        if (phase != CoolProp::iphase_not_imposed) st.unspecify_phase();
        return v;
    } catch (...) {
        if (phase != CoolProp::iphase_not_imposed) {
            try {
                fluid_->state()->unspecify_phase();
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // already failed; nothing more to undo
            }
        }
        return Detail::NaN;
    }
}

}  // namespace Plot
}  // namespace CoolProp

#ifdef ENABLE_CATCH

#    include <catch2/catch_all.hpp>
#    include <type_traits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

/// Index of the point of a line closest to x.
size_t closest_x(const CoolProp::Plot::Isoline& line, double x) {
    size_t best = 0;
    double distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < line.x.size(); ++i) {
        double d = std::fabs(line.x[i] - x);
        if (d < distance) {
            distance = d;
            best = i;
        }
    }
    return best;
}

double smallest_distance_y(const CoolProp::Plot::Isoline& line, double y) {
    double distance = std::numeric_limits<double>::infinity();
    for (double v : line.y)
        distance = std::min(distance, std::fabs(v - y));
    return distance;
}

}  // namespace

TEST_CASE("An isobar runs horizontally through the dome in a T-s chart", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R290", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    const double p_Pa = 5e5;
    CoolProp::Plot::Isoline line = plot.calc_isolines(CoolProp::iP, {p_Pa}, 50).front();
    CHECK(line.n_inserted == 2);
    CHECK(line.n_failed == 0);
    CHECK(line.method == "flash");
    double T_sat_K = CoolProp::PropsSI("T", "P", p_Pa, "Q", 0, "R290");
    double s_liq = CoolProp::PropsSI("S", "P", p_Pa, "Q", 0, "R290"), s_vap = CoolProp::PropsSI("S", "P", p_Pa, "Q", 1, "R290");
    size_t i = closest_x(line, s_liq), j = closest_x(line, s_vap);
    // The bubble and dew states are on the line, exactly, with the tie line
    // between them at the saturation temperature
    CHECK_THAT(line.x[i], WithinRel(s_liq, 1e-9));
    CHECK_THAT(line.x[j], WithinRel(s_vap, 1e-9));
    REQUIRE(j > i + 1);
    for (size_t k = i; k <= j; ++k)
        CHECK_THAT(line.y[k], WithinRel(T_sat_K, 1e-9));
    // T rises monotonically along the isobar; s too
    for (size_t k = 1; k < line.x.size(); ++k) {
        CHECK(line.y[k] - line.y[k - 1] >= -1e-9);
        CHECK(line.x[k] > line.x[k - 1]);
    }
}

TEST_CASE("An isotherm in a p-h chart has its kinks on the dome", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R290", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::Isoline line = plot.calc_isolines(CoolProp::iT, {280.0}, 30).front();
    CHECK(line.n_inserted == 2);
    for (double q : {0.0, 1.0}) {
        double h = CoolProp::PropsSI("H", "T", 280.0, "Q", q, "R290");
        CHECK(std::fabs(line.x[closest_x(line, h)] - h) < 1e-6 * h);
    }
}

TEST_CASE("An isenthalp meets the bubble line exactly", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R290", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::Isoline line = plot.calc_isolines(CoolProp::iHmass, {3e5}, 40).front();
    CHECK(line.n_inserted >= 1);
    // Some point lies on the bubble line: its (T, s) is a saturated liquid state
    double closest = std::numeric_limits<double>::infinity(), largest = 0.0;
    for (size_t k = 0; k < line.x.size(); ++k) {
        double T = line.y[k];
        if (!(170.0 < T && T < 369.0)) continue;
        double s_liq = CoolProp::PropsSI("S", "T", T, "Q", 0, "R290");
        closest = std::min(closest, std::fabs(line.x[k] - s_liq));
        largest = std::max(largest, s_liq);
    }
    CHECK(closest < 1e-6 * largest);
}

TEST_CASE("The tracer agrees with the flash for a blend", "[Plot]") {
    const std::string blend = "HEOS::R32[0.7]&R125[0.3]";
    CoolProp::Plot::PropertyPlot plot(blend, CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::CalcSettings settings;
    settings.points = 25;
    settings.solver = "tracer";
    CoolProp::Plot::Isoline traced = plot.calc_isolines(CoolProp::iHmass, {4e5}, settings).front();
    settings.solver = "flash";
    CoolProp::Plot::Isoline flashed = plot.calc_isolines(CoolProp::iHmass, {4e5}, settings).front();
    CHECK(traced.method == "tracer");
    CHECK(flashed.method == "flash");
    CHECK(traced.n_valid() >= 20);
    CHECK(flashed.n_valid() >= 20);
    // Which of the two is faster is not checked here: a timing check is
    // flaky in a sharded parallel test run
    // The two lines step the pressure differently, so compare the two ways
    // of computing a point at the same pressures: inside and outside the dome
    CoolProp::Plot::Detail::Fluid fluid(blend);
    CoolProp::Plot::Detail::LineSolver tracer(fluid, CoolProp::iHmass, 4e5, CoolProp::iP, CoolProp::iSmass, CoolProp::iT, nullptr, "tracer");
    CoolProp::Plot::Detail::LineSolver flash(fluid, CoolProp::iHmass, 4e5, CoolProp::iP, CoolProp::iSmass, CoolProp::iT, nullptr, "flash");
    CHECK(tracer.method() == "tracer");
    CHECK(flash.method() == "flash");
    // CoolProp's H-p flash of a mixture refuses some two-phase states (its
    // residual check after the temperature search); compare where it converges
    int compared = 0;
    for (double p_Pa : logspace(1e5, 3e6, 7)) {
        CoolProp::Plot::Detail::Sample a, b;
        REQUIRE(tracer.sample(p_Pa, true, a));
        CHECK(a.source == "tracer");
        if (!flash.sample(p_Pa, true, b)) continue;
        ++compared;
        CHECK_THAT(a.x, WithinRel(b.x, 1e-8));
        CHECK_THAT(a.y, WithinRel(b.y, 1e-8));
    }
    CHECK(compared >= 5);
    CHECK(tracer.n_fallback() == 0);
}

TEST_CASE("The automatic solver traces mixtures only", "[Plot]") {
    CoolProp::Plot::PropertyPlot pure("HEOS::R290", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    CHECK(pure.calc_isolines(CoolProp::iP, {5e5}, 10).front().method == "flash");
    CoolProp::Plot::PropertyPlot blend("HEOS::R32[0.7]&R125[0.3]", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    CHECK(blend.calc_isolines(CoolProp::iP, {5e5}, 10).front().method == "tracer");
    CoolProp::Plot::CalcSettings settings;
    settings.solver = "newton";
    CHECK_THROWS_AS(blend.calc_isolines(CoolProp::iP, {5e5}, settings), CoolProp::ValueError);
}

TEST_CASE("The dome of a mixture comes from the phase envelope", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R32[0.7]&R125[0.3]", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CHECK(plot.is_mixture());
    CoolProp::Plot::Isolines lines = plot.calc_isolines(CoolProp::iQ, {0.0, 1.0}, 50);
    const CoolProp::Plot::Isoline& bubble = lines[0];
    const CoolProp::Plot::Isoline& dew = lines[1];
    CHECK(bubble.method == "envelope");
    CHECK(dew.method == "envelope");
    CHECK(bubble.n_failed == 0);
    CHECK(dew.n_failed == 0);
    // The dome closes: the bubble line starts where the dew line ends
    CHECK(bubble.x.front() == dew.x.back());
    CHECK(bubble.y.front() == dew.y.back());
    // Every point is an envelope state: in a T-rho chart the bubble line
    // runs through the envelope's temperatures and bulk densities (CoolProp
    // keeps the bulk phase in the "vap" fields on both branches)
    std::shared_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("HEOS", std::vector<std::string>{"R32", "R125"}));
    state->set_mole_fractions({0.7, 0.3});
    state->build_phase_envelope("veryfine");
    const CoolProp::PhaseEnvelopeData& data = state->get_phase_envelope_data();
    CoolProp::Plot::PropertyPlot trho_plot("HEOS::R32[0.7]&R125[0.3]", CoolProp::iT, CoolProp::iDmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::Isoline trho = trho_plot.calc_isolines(CoolProp::iQ, {0.0}, 50).front();
    size_t k = 1;
    for (size_t i = 0; i < data.T.size(); ++i) {
        if (data.Q[i] != 0.0) continue;
        REQUIRE(k < trho.x.size());
        CHECK_THAT(trho.y[k], WithinRel(data.T[i], 1e-12));
        CHECK_THAT(trho.x[k], WithinRel(data.rhomolar_vap[i] * state->molar_mass(), 1e-12));
        CHECK_THAT(bubble.y[k], WithinRel(data.p[i], 1e-12));
        ++k;
    }
    CHECK(k == trho.x.size());
    // The bulk on the bubble line is liquid: at the same temperature the
    // bubble state agrees with a QT flash (away from the critical point)
    size_t i = bubble.x.size() * 3 / 4;
    state->update(CoolProp::QT_INPUTS, 0.0, trho.y[i]);
    CHECK_THAT(state->hmass(), WithinRel(bubble.x[i], 1e-4));
    CHECK_THAT(state->p(), WithinRel(bubble.y[i], 1e-4));
}

TEST_CASE("The lines of a blend meet the dome near the critical point", "[Plot]") {
    // The saturation flash of the blend fails above about 340 K, so these
    // crossings can only come from the phase envelope
    CoolProp::Plot::PropertyPlot plot("HEOS::R32[0.7]&R125[0.3]", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::Isoline isochore = plot.calc_isolines(CoolProp::iDmass, {300.0}, 40).front();
    CHECK(isochore.n_inserted == 1);
    CoolProp::Plot::Isoline isotherm = plot.calc_isolines(CoolProp::iT, {340.0}, 40).front();
    CHECK(isotherm.n_inserted == 2);
    // The two inserted states of the isotherm sit at 340 K on the envelope
    std::shared_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("HEOS", std::vector<std::string>{"R32", "R125"}));
    state->set_mole_fractions({0.7, 0.3});
    state->build_phase_envelope("veryfine");
    const CoolProp::PhaseEnvelopeData& data = state->get_phase_envelope_data();
    for (double quality : {0.0, 1.0}) {
        std::vector<std::pair<double, double>> branch;
        for (size_t i = 0; i < data.T.size(); ++i)
            if (data.Q[i] == quality) branch.emplace_back(data.T[i], data.p[i]);
        std::sort(branch.begin(), branch.end());
        REQUIRE(branch.size() >= 2);
        size_t j = 1;
        while (j + 1 < branch.size() && branch[j].first < 340.0)
            ++j;
        // cppcheck-suppress containerOutOfBounds // REQUIRE above stops the test when a branch has fewer than two points
        double t = (340.0 - branch[j - 1].first) / (branch[j].first - branch[j - 1].first);
        // cppcheck-suppress containerOutOfBounds // as above
        double p_sat = branch[j - 1].second + t * (branch[j].second - branch[j - 1].second);
        CHECK(smallest_distance_y(isotherm, p_sat) < 1e-4 * p_sat);
    }
}

TEST_CASE("The quality lines of a pure fluid end at the critical point", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R290", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Def);
    CHECK_FALSE(plot.is_mixture());
    CoolProp::Plot::Isolines lines = plot.calc_isolines(CoolProp::iQ, {0.0, 0.5, 1.0}, 40);
    double p_crit = CoolProp::PropsSI("pcrit", "T", 0, "P", 0, "R290");
    double h_crit = CoolProp::PropsSI("H", "T", CoolProp::PropsSI("Tcrit", "T", 0, "P", 0, "R290"), "Dmass",
                                      CoolProp::PropsSI("rhocrit", "T", 0, "P", 0, "R290"), "R290");
    for (const CoolProp::Plot::Isoline& line : lines) {
        CHECK(line.n_failed == 0);
        CHECK(line.size() == 41);
        CHECK_THAT(line.y.back(), WithinRel(p_crit, 1e-9));
        CHECK_THAT(line.x.back(), WithinRel(h_crit, 1e-9));
        CHECK(line.n_gaps() == 0);
    }
    // The bubble line's enthalpy rises with the pressure
    for (size_t k = 1; k < lines[0].x.size(); ++k)
        CHECK(lines[0].x[k] > lines[0].x[k - 1]);
}

TEST_CASE("States computed outside the axes are told apart from failures", "[Plot]") {
    using namespace CoolProp::Plot::Detail;
    // The first, careful look fails; every other state lies above the chart
    Evaluate evaluate = [](double u, bool /*careful*/, Sample& out) {
        if (u == 0.0) return false;
        out.valid = true;
        out.x = u;
        out.y = 2.0;
        out.inside = false;
        return true;
    };
    Metric metric = {{0.0, 1.0}, {0.0, 1.0}, false, false};
    Trace found = trace(evaluate, metric, 20);
    CHECK(found.n_failed >= 1);
    CHECK(found.n_outside_ok >= 1);
    CHECK(std::none_of(found.samples.begin(), found.samples.end(), [](const Sample& s) { return s.valid; }));
    CHECK(short_line_reason(0, found.n_evaluated, found.n_ok, found.n_outside_ok).find("lie in the area the axes show") != std::string::npos);
    CHECK(short_line_reason(0, 5, 0, 0).rfind("none of the 5 states", 0) == 0);
    CHECK(short_line_reason(1, 5, 1, 0).rfind("only 1 point(s)", 0) == 0);
    CHECK(short_line_reason(2, 5, 2, 0).empty());
}

TEST_CASE("Review fixes: limits, sweep values, trivial solutions", "[Plot]") {
    using namespace CoolProp::Plot::Detail;
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    // A NaN limit is a mistake and throws; an infinite one is the limit of the equation of state
    CHECK_THROWS_AS(CoolProp::Plot::PropertyPlot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, nan, 400.0, 1e5, 1e7), CoolProp::ValueError);
    CHECK_NOTHROW(CoolProp::Plot::PropertyPlot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, 200.0, inf, 1e5, 1e7));
    // Sweep values equal up to rounding order the dome ends by quality, not by noise
    CHECK(same_sweep_value(300.0, 300.0 * (1.0 + 1e-15)));
    CHECK_FALSE(same_sweep_value(300.0, 300.001));
    const double w_noisy = 300.0 * (1.0 + 1e-15);                        // the liquid end, a few bits above the vapour end
    CHECK(dome_end_comes_first(w_noisy, 0.0, 300.0, 1.0, true));         // sweep in T: liquid first
    CHECK_FALSE(dome_end_comes_first(w_noisy, 0.0, 300.0, 1.0, false));  // sweep in p: vapour first
    CHECK(dome_end_comes_first(300.0, 1.0, 300.001, 0.0, true));         // distinct sweep values decide
    // A state whose saturated outputs throw is not a pair of distinct phases (the check fails closed;
    // the incompressible backend has no saturated outputs at all, so its call throws)
    Fluid liquid("INCOMP::MEG[0.3]");
    CoolProp::AbstractState& brine = liquid.update(CoolProp::iT, 300.0, CoolProp::iP, 1e5);
    CHECK_THROWS(brine.saturated_liquid_keyed_output(CoolProp::iDmolar));
    CHECK_FALSE(phases_are_distinct(brine));
}

TEST_CASE("Gaps count only inside a line", "[Plot]") {
    const double nan = std::nan("");
    CoolProp::Plot::Isoline line;
    line.x = {nan, 1, 2, nan, 4, nan, nan};
    line.y = line.x;
    CHECK(line.n_valid() == 3);
    CHECK(line.n_gaps() == 1);
    line.x = {1, 2, nan};
    line.y = line.x;
    CHECK(line.n_gaps() == 0);
}

TEST_CASE("Check value_at for p-h plots", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CHECK_THAT(plot.value_at(CoolProp::iP, 300000 /*J/kg*/, 200000 /*Pa*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iT, 300000, 200000), WithinRel(CoolProp::PropsSI("T", "P", 200000, "H", 300000, "R134a"), 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iQ, 300000, 200000), WithinRel(CoolProp::PropsSI("Q", "P", 200000, "H", 300000, "R134a"), 1e-10));
    CHECK(std::isnan(plot.value_at(CoolProp::iT, -1.0, -1.0)));
}

TEST_CASE("The axes and the supported isolines of a plot", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Log);
    std::vector<double> box = plot.tp_box();
    CHECK_THAT(box[0], WithinAbs(173.15, 1e-9));
    // The equation of state of R134a ends below the 493.15 K of the ACHP limits
    CHECK_THAT(box[1], WithinAbs(CoolProp::PropsSI("Tmax", "T", 0, "P", 0, "R134a"), 1e-9));
    CHECK(box[1] < 493.15);
    CHECK_THAT(box[2], WithinAbs(0.25e5, 1e-9));
    CHECK_THAT(box[3], WithinRel(2.25 * CoolProp::PropsSI("pcrit", "T", 0, "P", 0, "R134a"), 1e-6));
    CHECK(plot.yaxis.min == box[2]);         // NOLINT(cppcoreguidelines-pro-type-union-access)
    CHECK(plot.yaxis.max == box[3]);         // NOLINT(cppcoreguidelines-pro-type-union-access)
    CHECK(plot.xaxis.min < plot.xaxis.max);  // NOLINT(cppcoreguidelines-pro-type-union-access)
    std::vector<CoolProp::parameters> keys = plot.supported_isoline_keys();
    CHECK(keys == std::vector<CoolProp::parameters>{CoolProp::iT, CoolProp::iSmass, CoolProp::iDmass, CoolProp::iUmass, CoolProp::iQ});
    CHECK(plot.sweep_key(CoolProp::iT) == CoolProp::iP);
    CHECK(plot.sweep_key(CoolProp::iDmass) == CoolProp::iT);
    CHECK(plot.sweep_key(CoolProp::iP) == CoolProp::iundefined_parameter);
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iP, {1e5}, 10), CoolProp::ValueError);
    CoolProp::Plot::Range range = plot.isoline_range(CoolProp::iT);
    CHECK(range.min == box[0]);
    CHECK(range.max == box[1]);
    std::vector<double> values = CoolProp::Plot::generate_values_in_range(CoolProp::iP, {1e5, 1e7}, 3);
    REQUIRE(values.size() == 3);
    CHECK_THAT(values[1], WithinRel(1e6, 1e-9));
}

TEST_CASE("Settings are checked up front and a failing line keeps its reason", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    CoolProp::Plot::CalcSettings settings;
    settings.solver = "newton";  // checked for lines of constant quality as well
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iQ, {0.5}, settings), CoolProp::ValueError);
    settings = CoolProp::Plot::CalcSettings();
    settings.saturation_spacing = "cubic";
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iT, {300.0}, settings), CoolProp::ValueError);
    settings = CoolProp::Plot::CalcSettings();
    settings.critical_margin = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iQ, {0.5}, settings), CoolProp::ValueError);
    settings.critical_margin = 0.5;
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iQ, {0.5}, settings), CoolProp::ValueError);
    settings.critical_margin = 0.4;  // below 0.5, but more than half the saturation range of R134a
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iQ, {0.5}, settings), CoolProp::ValueError);
    CHECK_THROWS_AS(plot.calc_isolines(CoolProp::iQ, {0.5, 1.5}, 10), CoolProp::ValueError);
    // A value no state can have: that line comes back empty with its reason, the other one is drawn
    CoolProp::Plot::Isolines lines = plot.calc_isolines(CoolProp::iT, {300.0, -50.0}, 20);
    REQUIRE(lines.size() == 2);
    CHECK(lines[0].n_valid() >= 10);
    CHECK(lines[1].n_valid() == 0);
    CHECK_FALSE(lines[1].reason.empty());
    // A plot can be moved (into a container, out of a function)
    static_assert(std::is_nothrow_move_constructible_v<CoolProp::Plot::PropertyPlot>, "PropertyPlot is movable");
    std::vector<CoolProp::Plot::PropertyPlot> plots;
    plots.push_back(std::move(plot));
    CHECK(plots.back().calc_isolines(CoolProp::iT, {300.0}, 20).front().n_valid() >= 10);
}

TEST_CASE("Absolute limits and a fluid without a critical point", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::Water", CoolProp::iP, CoolProp::iHmass, 300.0, 600.0, 1e5, 1e7);
    std::vector<double> box = plot.tp_box();
    CHECK(box == std::vector<double>{300.0, 600.0, 1e5, 1e7});
    CHECK(plot.calc_isolines(CoolProp::iT, {400.0}, 20).front().n_valid() >= 10);
    // An incompressible fluid has no saturation curve: relative limits cannot be resolved
    CHECK_THROWS_AS(CoolProp::Plot::PropertyPlot("INCOMP::MEG[0.5]", CoolProp::iP, CoolProp::iT, CoolProp::Plot::TPLimits::Def),
                    CoolProp::ValueError);
    CoolProp::Plot::PropertyPlot incomp("INCOMP::MEG[0.5]", CoolProp::iHmass, CoolProp::iT, 270.0, 350.0, 1e5, 1e6);
    CHECK(incomp.calc_isolines(CoolProp::iP, {5e5}, 20).front().n_valid() >= 10);
    CHECK_THROWS_AS(incomp.calc_isolines(CoolProp::iQ, {0.5}, 20), CoolProp::ValueError);
    // The cursor readout works without a phase to impose (a p-T chart: x = T, y = p)
    CoolProp::Plot::PropertyPlot incomp_pT("INCOMP::MEG[0.5]", CoolProp::iP, CoolProp::iT, 270.0, 350.0, 1e5, 1e6);
    CHECK_THAT(incomp_pT.value_at(CoolProp::iDmass, 300.0, 5e5), WithinRel(CoolProp::PropsSI("D", "T", 300.0, "P", 5e5, "INCOMP::MEG[0.5]"), 1e-9));
    CHECK_THAT(incomp_pT.value_at(CoolProp::iHmass, 300.0, 5e5), WithinRel(CoolProp::PropsSI("H", "T", 300.0, "P", 5e5, "INCOMP::MEG[0.5]"), 1e-9));
    // Axes CoolProp has no flash for are refused at construction
    CHECK_THROWS_AS(CoolProp::Plot::PropertyPlot("HEOS::R134a", CoolProp::iHmass, CoolProp::iUmass, CoolProp::Plot::TPLimits::Achp),
                    CoolProp::ValueError);
    CHECK(CoolProp::Plot::generate_values_in_range(CoolProp::iP, {1e5, 1e7}, 0) == std::vector<double>{1e5});
    CHECK(CoolProp::Plot::generate_values_in_range(CoolProp::iP, {1e5, 1e7}, 1) == std::vector<double>{1e5});
}

#endif  // ENABLE_CATCH
