/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INCOMPRESSIBLEFLUID_H_
#define INCOMPRESSIBLEFLUID_H_

#include "CoolProp/DataStructures.h"
#include "CoolProp/fluids/Helmholtz.h"
#include "CoolProp/numerics/Solvers.h"

#include <cmath>
#include <exception>
#include <cstdint>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cassert>
#include <iterator>

#include <Eigen/Core>
#include "CoolProp/numerics/PolyMath.h"
#include "CoolProp/numerics/MatrixMath.h"

namespace CoolProp {

struct IncompressibleData
{
    enum IncompressibleTypeEnum
    {
        INCOMPRESSIBLE_NOT_SET,
        INCOMPRESSIBLE_POLYNOMIAL,
        INCOMPRESSIBLE_EXPPOLYNOMIAL,
        INCOMPRESSIBLE_EXPONENTIAL,
        INCOMPRESSIBLE_LOGEXPONENTIAL,
        INCOMPRESSIBLE_POLYOFFSET,
        INCOMPRESSIBLE_CHEBYSHEV
    };
    IncompressibleTypeEnum type;
    Eigen::MatrixXd coeffs;  //TODO: Can we store the Eigen::Matrix objects more efficiently?

    /// Chebyshev fits only (type == INCOMPRESSIBLE_CHEBYSHEV). The rows of
    /// coeffs are Chebyshev-in-T coefficients on [cheb_Tmin, cheb_Tmax]; the
    /// columns multiply powers of (x - cheb_xbase).
    double cheb_Tmin = 0.0, cheb_Tmax = 0.0, cheb_xbase = 0.0;

    /// Derived coefficient matrices for the caloric fits (density and
    /// specific heat, type INCOMPRESSIBLE_POLYNOMIAL or INCOMPRESSIBLE_CHEBYSHEV).
    /// They are built once at load time by IncompressibleFluid::validate(), so
    /// that evaluating h and s is a plain evaluation of a stored matrix, in the
    /// same basis as coeffs:
    ///  - ddT:        coefficients of d/dT (the density fit needs this for
    ///                the pressure terms of h and s; both fits use it as the
    ///                slope when T_rho and T_c invert them with Newton)
    ///  - ddT2:       density fits only, coefficients of d2/dT2 (the exact
    ///                slope of h and s at constant p in the HmassP and PSmass
    ///                flashes)
    ///  - int_dT:     coefficients of the indefinite Int f dT (enthalpy)
    ///  - int_dTdivT: coefficients of the indefinite Int f/T dT (entropy)
    ///  - int_dTdivT_log: polynomial fits only. Int f/T dT of a polynomial
    ///                in (T - Tbase) is a polynomial plus a ln(T) term; this
    ///                single row holds the factor of ln(T) as a polynomial in
    ///                (x - xbase). The Chebyshev fits need no log term: they
    ///                re-expand f/T itself (see chebDivideByTMatrix), which
    ///                converges geometrically since T > 0.
    /// Integration constants are irrelevant: the backend always subtracts the
    /// same expression evaluated at the reference state.
    Eigen::MatrixXd ddT, ddT2, int_dT, int_dTdivT, int_dTdivT_log;

    IncompressibleData()
      : type(INCOMPRESSIBLE_NOT_SET) {

        };
};

#if !defined(NO_FMTLIB) && FMT_VERSION >= 90000
inline int format_as(IncompressibleData::IncompressibleTypeEnum type) {
    return fmt::underlying(type);
}
#endif

/// One fit with the composition already put in: a series in T alone.
/**
The polynomial and Chebyshev fits are matrices, rows in T and columns in
(x - xbase). For a state the composition is fixed, so each row can be
reduced to a single number once, and every later evaluation is a short sum
over T. The reduced coefficients are computed exactly as the full evaluation
computes them (Horner over the columns, row by row), and the sum over T runs
in the same order, so the result is bit for bit the same (on x86-64; where a
compiler fuses multiply-adds, as on aarch64, the two inlined copies can differ
by a few units in the last place).
*/
struct FitAtComposition
{
    enum Kind : std::uint8_t
    {
        NONE,       ///< not reduced (another fit type): evaluate the full fit instead
        POWER,      ///< sum a[k] * (T - T_offset)^k, Horner's scheme
        CHEBYSHEV,  ///< sum a[k] * T_k(u), u = (2 T - T_offset) / T_span, Clenshaw
    };
    Kind kind = NONE;
    bool exponential = false;  ///< exppolynomial fits: the value is exp(series)
    std::vector<double> a;
    double T_offset = 0.0, T_span = 1.0;
    [[nodiscard]] double eval(double T) const;
};

/// All fits one state evaluates, at one composition and for one choice of
/// caloric fit (classic or Chebyshev). Built by IncompressibleFluid::fixComposition.
struct CompositionFits
{
    double x = NAN;
    bool chebyshev = true;
    unsigned long data_version = 0;  ///< IncompressibleFluid::data_version these were built from
    FitAtComposition rho, drhodT, d2rhodT2, c, int_c, int_c_divT, int_c_divT_log, visc, cond;
};

/// A property provider for incompressible solutions and pure fluids
/**
This fluid instance is populated using an entry from a JSON file
*/
class IncompressibleFluid
{

   protected:
    bool strict;
    /// Changes every time validate() runs or a fit, Tbase or xbase is set, so
    /// that a backend holding CompositionFits of this fluid sees the change
    /// (also when add_fluids_as_JSON replaces the fluid in place) and builds
    /// them again.
    unsigned long data_version = 0;
    /// A new data_version, unique across all fluids (see validate()). Every
    /// setter of data that a backend caches calls it, as validate() does.
    static unsigned long nextDataVersion();

    std::string name;
    std::string description;
    std::string reference;

    double Tmin, Tmax;
    double xmin, xmax;
    composition_types xid;

    double TminPsat;
    double xbase, Tbase;

    /// These are the objects that hold the coefficients
    /** Note that all polynomials require a 2-dimensional array
     *  of coefficients. This array may have only one row or
     *  column, but the structure should be 2D. This behaviour is
     *  hard-coded in the JSON file reader that resides inside
     *  the IncompressibleLibrary.cpp
     *  All other functions, also polyoffset, can only handle 1D
     *  input and throw an error if you feed them other coefficients.
     */

    /// Density coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData density;
    /// Specific heat coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     *  Fails for all other forms than polynomial due to the automatic
     *  integration for internal energy and entropy.
     */
    IncompressibleData specific_heat;
    /// Optional Chebyshev alternatives to density and specific_heat
    /** The fluid data may ship both a classic fit ("density") and a
     *  Chebyshev fit ("density_cheb") of the same property; both are kept.
     *  Which one is evaluated is chosen per call, see useChebyshevFit() and the
     *  INCOMP_PREFER_CHEBYSHEV_CALORIC configuration key. Either may be
     *  missing (type INCOMPRESSIBLE_NOT_SET), in which case the other is used.
     */
    IncompressibleData density_cheb, specific_heat_cheb;
    /// Viscosity coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData viscosity;
    /// Conductivity coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData conductivity;
    /// Saturation pressure coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData p_sat;
    /// Freezing temperature coefficients
    /** If 2D, the rows are concentration and the columns are pressure.
     *  If 1D, should be a column vector of concentration coefficients
     */
    IncompressibleData T_freeze;

    /// Mass fraction conversion coefficients
    /** If the fluid type is mass-based, it does not do anything. Otherwise,
     *  it converts the mass fraction to the required input.
     */
    IncompressibleData mass2input;
    /// Volume fraction conversion coefficients
    /** If the fluid type is volume-based, it does not do anything. Otherwise,
     *  it converts the volume fraction to the required input.
     */
    IncompressibleData volume2input;
    /// Mole fraction conversion coefficients
    /** If the fluid type is mole-based, it does not do anything. Otherwise,
     *  it converts the mole fraction to the required input.
     */
    IncompressibleData mole2input;

    Polynomial2DFrac poly;

    /// Decide which fit to evaluate for a caloric property.
    /** True (use the Chebyshev fit) if chebyshev is true and that fit exists,
     *  or if the classic fit does not exist; false (use the classic fit)
     *  otherwise. With neither fit present the classic one is chosen; it is
     *  unset, and evaluating it reports the missing property.
     */
    static bool useChebyshevFit(const IncompressibleData& classic_fit, const IncompressibleData& chebyshev_fit, bool chebyshev) {
        const bool has_classic = classic_fit.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET;
        const bool has_chebyshev = chebyshev_fit.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET;
        return has_chebyshev && (chebyshev || !has_classic);
    }
    [[nodiscard]] const IncompressibleData& densityFit(bool chebyshev) const {
        return useChebyshevFit(density, density_cheb, chebyshev) ? density_cheb : density;
    }
    [[nodiscard]] const IncompressibleData& specificHeatFit(bool chebyshev) const {
        return useChebyshevFit(specific_heat, specific_heat_cheb, chebyshev) ? specific_heat_cheb : specific_heat;
    }

    /// Throw unless data is a Chebyshev fit or unset; see setDensityChebyshev.
    static void requireChebyshevOrUnset(const IncompressibleData& data, const char* what) {
        if (data.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET && data.type != IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV) {
            throw ValueError(
              format("The Chebyshev %s fit slot only takes a Chebyshev fit, not one of type [%d].", what, static_cast<int>(data.type)));
        }
    }

    /// Build the error for a property this fluid carries no fit for.
    /** Returned, not thrown, so each call site keeps its own throw statement.
     *  See the definition in IncompressibleFluid.cpp for why.
     */
    ValueError notDefinedError(const std::string& property) const;
    /// Evaluate the matrix M of a caloric fit (fit.coeffs or one derived from
    /// it in validate()) for either kind, polynomial or Chebyshev. Raises the
    /// "not defined" error naming property for a missing fit, and an error
    /// naming quantity for a fit type that has no such matrix.
    double evalCaloricFit(const IncompressibleData& fit, const Eigen::MatrixXd& M, double T, double x, const char* property,
                          const char* quantity) const;

    // Forward declaration of the some internal functions
    //double h_u(double T, double p, double x);
    //double u_h(double T, double p, double x);

   public:
    IncompressibleFluid()
      : strict(true), Tmin(_HUGE), Tmax(_HUGE), xmin(_HUGE), xmax(_HUGE), xid(IFRAC_UNDEFINED), TminPsat(_HUGE), xbase(_HUGE), Tbase(_HUGE) {

        };
    // The virtual destructor below suppresses the implicitly-declared move
    // constructor/assignment operator (even though it is only '= default'),
    // so state them explicitly - otherwise callers that std::move() a fluid
    // (e.g. JSONIncompressibleLibrary::add_one) silently fall back to a copy.
    IncompressibleFluid(const IncompressibleFluid&) = default;
    IncompressibleFluid(IncompressibleFluid&&) = default;
    IncompressibleFluid& operator=(const IncompressibleFluid&) = default;
    IncompressibleFluid& operator=(IncompressibleFluid&&) = default;
    virtual ~IncompressibleFluid() = default;

    std::string getName() const {
        return name;
    }
    std::string get_name() const {
        return getName();
    }  // For backwards-compatibility.
    std::string getDescription() const {
        return description;
    }
    std::string getReference() const {
        return reference;
    }

    double getTmax() const {
        return Tmax;
    }
    double getTmin() const {
        return Tmin;
    }
    /// Lowest valid temperature at this pressure and composition: Tmin, or
    /// the freezing point where that is higher (the limit checkT() applies).
    /// The inverse functions search from here, so that a fit extrapolated
    /// below the freezing point cannot offer a second, frozen solution.
    double TminAt(double p, double x);
    /// Highest valid temperature at this pressure and composition: Tmax, or
    /// the saturation temperature where the liquid boils below Tmax (the limit
    /// checkP() applies, where p must not be below psat). The inverse functions
    /// search up to here, so that a fit extrapolated above the boiling point,
    /// where h or s can fall with T again, cannot hide the liquid solution.
    double TmaxAt(double p, double x);
    /// Temperature margin at the boiling point, in K. TmaxAt() returns the
    /// temperature where the liquid ends (the saturation temperature, or
    /// TminPsat if the liquid boils just above it) plus 0.3 to 0.7 times this
    /// margin, and checkP() accepts a state within the full margin, so that a
    /// state solved for at saturation, which can land a few ulp above it, is
    /// still accepted as liquid.
    [[nodiscard]] double saturationMargin() const {
        return 2e-9 * Tmax;
    }
    /// True if (T, p) is a liquid state: below TminPsat (no saturation curve
    /// there), at or above psat(T), or liquid within saturationMargin(), that
    /// is at T minus the margin. checkP() and solveInLiquidRange() both use it,
    /// so a state one accepts the other does too.
    bool isLiquidAt(double T, double p, double x) {
        if (p_sat.type == IncompressibleData::INCOMPRESSIBLE_NOT_SET || T <= TminPsat) return true;
        if (psat(T, x) <= p) return true;
        const double T_lower = T - saturationMargin();
        return T_lower <= TminPsat || psat(T_lower, x) <= p;
    }
    /// Solve for a temperature inside the valid liquid range at (p, x).
    /** solve(T_lo, T_hi) must return the root on [T_lo, T_hi] or throw. The
     *  search runs on [TminAt(p, x), Tmax] first. When that gives a liquid
     *  state (p not below psat), it is the answer, and the boiling point is
     *  never computed, which keeps this cheap. When it gives a boiling state,
     *  or finds no root (the fit can turn over above the boiling point), the
     *  search is repeated on [TminAt(p, x), TmaxAt(p, x)].
     *  A template, so that the solver call is inlined: Solve is any callable,
     *  typically a lambda around BoundedNewton.
     */
    template <typename Solve>
    double solveInLiquidRange(double p, double x, Solve solve) {
        const double T_lo = TminAt(p, x);
        if (p_sat.type == IncompressibleData::INCOMPRESSIBLE_NOT_SET || !(TminPsat < Tmax)) return solve(T_lo, Tmax);
        // First pass on the full range. Keep its error, if any: when the second
        // pass fails as well, the first error is the one that explains why (a
        // target below freezing, a NaN input, a density reached twice).
        std::exception_ptr first_error;
        bool boiling = false;
        try {
            const double T = solve(T_lo, Tmax);
            if (isLiquidAt(T, p, x)) return T;  // a liquid state: done
            boiling = true;                     // the solution boils at p
        } catch (const CoolPropBaseError&) {
            first_error = std::current_exception();  // maybe only the extrapolation above boiling
        }
        double T_boil = Tmax;
        try {
            T_boil = TmaxAt(p, x);
        } catch (const CoolPropBaseError&) {
            if (first_error) std::rethrow_exception(first_error);
            throw;
        }
        if (!(T_boil < Tmax)) {
            // No boiling point inside the range; with a monotonic psat a
            // boiling first result cannot happen here
            if (first_error) std::rethrow_exception(first_error);
            throw ValueError(format("%s: no liquid state at p = %g Pa", name.c_str(), p));
        }
        // T_boil lies 0.3 to 0.7 saturationMargin() above the temperature where
        // the liquid ends, while isLiquidAt() accepts up to one margin above it.
        // The second pass searches 0.75 margin further, to 1.05 to 1.45
        // margins above saturation, so it reaches every state isLiquidAt()
        // accepts (a fit that turns over can have its only liquid root there),
        // and a root it finds beyond that band is rejected as boiling by the
        // same test as in the first pass.
        const double T_hi = std::min(Tmax, T_boil + 0.75 * saturationMargin());
        double T = 0.0;
        try {
            T = solve(T_lo, T_hi);
        } catch (const CoolPropBaseError& e) {
            if (boiling) {
                throw ValueError(format("%s: no liquid state at p = %g Pa, where the liquid boils at %g K (%s)", name.c_str(), p, T_boil, e.what()));
            }
            std::rethrow_exception(first_error);
        }
        if (!isLiquidAt(T, p, x)) {
            throw ValueError(
              format("%s: no liquid state at p = %g Pa, where the liquid boils at %g K (solution %.17g K)", name.c_str(), p, T_boil, T));
        }
        return T;
    }
    /// Reduce every fit a state evaluates to a series in T at composition x
    /// (see FitAtComposition). Fits of other types are left as NONE.
    void fixComposition(double x, bool chebyshev, CompositionFits& out) const;
    [[nodiscard]] unsigned long getDataVersion() const {
        return data_version;
    }
    double getxmax() const {
        return xmax;
    }
    double getxmin() const {
        return xmin;
    }
    composition_types getxid() const {
        return xid;
    }
    double getTminPsat() const {
        return TminPsat;
    }
    double getTbase() const {
        return Tbase;
    }
    double getxbase() const {
        return xbase;
    }

    void setName(const std::string& name) {
        this->name = name;
    }
    void setDescription(const std::string& description) {
        this->description = description;
    }
    void setReference(const std::string& reference) {
        this->reference = reference;
    }
    void setTmax(double Tmax) {
        this->Tmax = Tmax;
    }
    void setTmin(double Tmin) {
        this->Tmin = Tmin;
    }
    void setxmax(double xmax) {
        this->xmax = xmax;
    }
    void setxmin(double xmin) {
        this->xmin = xmin;
    }
    void setxid(composition_types xid) {
        this->xid = xid;
    }
    void setTminPsat(double TminPsat) {
        this->TminPsat = TminPsat;
    }
    void setTbase(double Tbase) {
        this->Tbase = Tbase;
        data_version = nextDataVersion();
        // The entropy integral of a polynomial fit has the value of Tbase
        // built in (see int_dTdivT_log), so drop it: evaluating s then throws
        // until validate() rebuilds it, instead of silently mixing the old
        // and the new Tbase. ddT and int_dT do not depend on the value of
        // Tbase (they are evaluated at the current T - Tbase) and stay valid.
        for (IncompressibleData* data : {&density, &specific_heat}) {
            if (data->type == IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL) {
                data->int_dTdivT.resize(0, 0);
                data->int_dTdivT_log.resize(0, 0);
            }
        }
    }
    void setxbase(double xbase) {
        this->xbase = xbase;
        data_version = nextDataVersion();
    }

    /// Setters for the coefficients
    void setDensity(IncompressibleData density) {
        this->density = std::move(density);
        data_version = nextDataVersion();
    }
    void setSpecificHeat(IncompressibleData specific_heat) {
        this->specific_heat = std::move(specific_heat);
        data_version = nextDataVersion();
    }
    /// Setters for the optional Chebyshev alternatives, see density_cheb
    /** These slots only take Chebyshev fits (or an unset fit): anything else
     *  there would be evaluated as the "Chebyshev" choice without being one.
     *  Checked here, not only in validate(), because rho() and c() also work
     *  on a fluid that never went through validate().
     */
    void setDensityChebyshev(IncompressibleData density_cheb) {
        requireChebyshevOrUnset(density_cheb, "density");
        this->density_cheb = std::move(density_cheb);
        data_version = nextDataVersion();
    }
    void setSpecificHeatChebyshev(IncompressibleData specific_heat_cheb) {
        requireChebyshevOrUnset(specific_heat_cheb, "specific heat");
        this->specific_heat_cheb = std::move(specific_heat_cheb);
        data_version = nextDataVersion();
    }
    void setViscosity(IncompressibleData viscosity) {
        this->viscosity = std::move(viscosity);
        data_version = nextDataVersion();
    }
    void setConductivity(IncompressibleData conductivity) {
        this->conductivity = std::move(conductivity);
        data_version = nextDataVersion();
    }
    void setPsat(IncompressibleData p_sat) {
        this->p_sat = std::move(p_sat);
    }
    void setTfreeze(IncompressibleData T_freeze) {
        this->T_freeze = std::move(T_freeze);
    }

    /// Setters for the concentration conversion coefficients
    void setMass2input(IncompressibleData mass2input) {
        this->mass2input = std::move(mass2input);
    }
    void setVolume2input(IncompressibleData volume2input) {
        this->volume2input = std::move(volume2input);
    }
    void setMole2input(IncompressibleData mole2input) {
        this->mole2input = std::move(mole2input);
    }

    /// A function to check coefficients and equation types.
    void validate();
    /// A function to test the density coefficients for 1D or 2D
    bool is_pure();

    /// The current value of the INCOMP_PREFER_CHEBYSHEV_CALORIC configuration
    /// key. It is the default for the chebyshev argument of the caloric
    /// functions below; the backend reads it once per state instead.
    static bool chebyshevCaloricPreferred();
    /// True if this fluid has Chebyshev fits for both density and specific heat
    [[nodiscard]] bool hasChebyshevCaloric() const {
        return density_cheb.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET
               && specific_heat_cheb.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET;
    }
    /// True if this fluid has classic fits for both density and specific heat
    /** A Chebyshev fit stored under the classic name does not count: with
     *  the setting at false it would be evaluated as the "classic" choice.
     */
    [[nodiscard]] bool hasClassicCaloric() const {
        auto is_classic = [](const IncompressibleData& data) {
            return data.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET && data.type != IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV;
        };
        return is_classic(density) && is_classic(specific_heat);
    }

   protected:
    /// Base functions that handle the custom function types
    double baseExponential(const IncompressibleData& data, double y, double ybase);
    double baseLogexponential(const IncompressibleData& data, double y, double ybase);
    double baseExponentialOffset(IncompressibleData data, double y);
    double basePolyOffset(const IncompressibleData& data, double y, double z = 0.0);

   public:
    /* All functions need T and p as input. Might not
     * be necessary, but gives a clearer structure.
     */
    /// Density as a function of temperature, pressure and composition.
    /** The caloric functions (rho, c, their derivatives and integrals, and
     *  their inverses) take a trailing chebyshev flag: true evaluates the
     *  Chebyshev fit where the fluid ships one, false the classic fit where
     *  the fluid ships one (see useChebyshevFit). It defaults to the current value
     *  of the INCOMP_PREFER_CHEBYSHEV_CALORIC configuration key.
     */
    double rho(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    /// Heat capacities as a function of temperature, pressure and composition.
    double c(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    double cp(double T, double p, double x) {
        throw ValueError(format("%s (%d): Please use the c-function instead.", __FILE__, __LINE__));
    }
    double cv(double T, double p, double x) {
        throw ValueError(format("%s (%d): Please use the c-function instead.", __FILE__, __LINE__));
    }
    /// Entropy as a function of temperature, pressure and composition.
    double s(double T, double p, double x) {
        throw ValueError(format("%s (%d): The internal calculations have changed, use the backend to calculate entropy from the partial derivatives.",
                                __FILE__, __LINE__));
    }
    /// Internal energy as a function of temperature, pressure and composition.
    double u(double T, double p, double x) {
        throw ValueError(
          format("%s (%d): The internal calculations have changed, use the backend to calculate internal energy from enthalpy.", __FILE__, __LINE__));
    }
    /// Enthalpy as a function of temperature, pressure and composition.
    double h(double T, double p, double x) {
        throw ValueError(
          format("%s (%d): The internal calculations have changed, use the backend to calculate enthalpy from the partial derivatives.", __FILE__,
                 __LINE__));
    }
    /// Viscosity as a function of temperature, pressure and composition.
    double visc(double T, double p, double x);
    /// Thermal conductivity as a function of temperature, pressure and composition.
    double cond(double T, double p, double x);
    /// Saturation pressure as a function of temperature and composition.
    double psat(double T, double x);
    /// Freezing temperature as a function of pressure and composition.
    double Tfreeze(double p, double x);

    /* Below are direct calculations of the derivatives. Nothing
     * special is going on, we simply use the polynomial class to
     * derive the different functions with respect to temperature.
     */
    /// Partial derivative of density
    //  with respect to temperature at constant pressure and composition
    double drhodTatPx(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    /// Second partial derivative of density with respect to temperature at constant pressure and composition
    double d2rhodT2atPx(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    ///// Partial derivative of entropy
    ////  with respect to temperature at constant pressure and composition
    //double dsdTatPx  (double T, double p, double x){return c(T,p,x)/T;};
    ///// Partial derivative of enthalpy
    ////  with respect to temperature at constant pressure and composition
    //double dhdTatPx  (double T, double p, double x){return c(T,p,x);};
    /// Partial derivative of entropy
    //  with respect to temperature at constant pressure and composition
    //  integrated in temperature
    double dsdTatPxdT(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    /// Partial derivative of enthalpy
    //  with respect to temperature at constant pressure and composition
    //  integrated in temperature
    double dhdTatPxdT(double T, double p, double x, bool chebyshev = chebyshevCaloricPreferred());

    /// Mass fraction conversion function
    /** If the fluid type is mass-based, it does not do anything. Otherwise,
     *  it converts the mass fraction to the required input. */
    double inputFromMass(double T, double x);
    /// Volume fraction conversion function
    /** If the fluid type is volume-based, it does not do anything. Otherwise,
     *  it converts the volume fraction to the required input. */
    double inputFromVolume(double T, double x);
    /// Mole fraction conversion function
    /** If the fluid type is mole-based, it does not do anything. Otherwise,
     *  it converts the mole fraction to the required input. */
    double inputFromMole(double T, double x);

    /* Some functions can be inverted directly, those are listed
     * here. It is also possible to solve for other quantities, but
     * that involves some more sophisticated processing and is not
     * done here, but in the backend, T(h,p) for example.
     */
    /// Temperature as a function of density, pressure and composition.
    double T_rho(double Dmass, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    /// Temperature as a function of heat capacities as a function of temperature, pressure and composition.
    double T_c(double Cmass, double p, double x, bool chebyshev = chebyshevCaloricPreferred());
    /// Temperature as a function of entropy as a function of temperature, pressure and composition.
    double T_s(double Smass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from entropy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Temperature as a function of internal energy as a function of temperature, pressure and composition.
    double T_u(double Umass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from internal energy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Temperature as a function of enthalpy, pressure and composition.
    double T_h(double Hmass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from enthalpy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Viscosity as a function of temperature, pressure and composition.
    double T_visc(double visc, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from viscosity is not implemented.", __FILE__, __LINE__));
    }
    /// Thermal conductivity as a function of temperature, pressure and composition.
    double T_cond(double cond, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from conductivity is not implemented.", __FILE__, __LINE__));
    }
    /// Saturation pressure as a function of temperature and composition.
    double T_psat(double psat, double x) {
        throw NotImplementedError(format("%s (%d): T from psat is not implemented.", __FILE__, __LINE__));
    }
    /// Composition as a function of freezing temperature and pressure.
    double x_Tfreeze(double Tfreeze, double p) {
        throw NotImplementedError(format("%s (%d): x from T_freeze is not implemented.", __FILE__, __LINE__));
    }

   protected:
    /* Define internal energy and enthalpy as functions of the
     * other properties to provide data in case there are no
     * coefficients.
     */
    /// Enthalpy from u, p and rho.
    /** Calculate enthalpy as a function of temperature and
     *  pressure employing functions for internal energy and
     *  density. Provides consistent formulations. */
    double h_u(double T, double p, double x) {
        return u(T, p, x) + p / rho(T, p, x);
    };

    /// Internal energy from h, p and rho.
    /** Calculate internal energy as a function of temperature
     *  and pressure employing functions for enthalpy and
     *  density. Provides consistent formulations. */
    double u_h(double T, double p, double x) {
        return h(T, p, x) - p / rho(T, p, x);
    };

    /*
     * Some more functions to provide a single implementation
     * of important routines.
     * We start with the check functions that can validate input
     * in terms of pressure p, temperature T and composition x.
     */
    /// Check validity of temperature input.
    /** Compares the given temperature T to the result of a
     *  freezing point calculation. This is not necessarily
     *  defined for all fluids, default values do not cause errors. */
    bool checkT(double T, double p, double x);

    /// Check validity of pressure input.
    /** Compares the given pressure p to the saturation pressure at
     *  temperature T and throws and exception if p is lower than
     *  the saturation conditions.
     *  The default value for psat is -1 yielding true if psat
     *  is not redefined in the subclass.
     *  */
    bool checkP(double T, double p, double x);

   public:
    /// Check validity of composition input.
    /** Compares the given composition x to a stored minimum and
     *  maximum value. Enforces the redefinition of xmin and
     *  xmax since the default values cause an error. */
    bool checkX(double x);

    /// Check validity of temperature, pressure and composition input.
    bool checkTPX(double T, double p, double x) {
        return (checkT(T, p, x) && checkP(T, p, x) && checkX(x));
    };
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEFLUID_H_ */
