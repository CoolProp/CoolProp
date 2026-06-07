#include "CoolProp/fluids/Ancillaries.h"

#include <cmath>
#include <limits>
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/AbstractState.h"
#include "CoolProp/Configuration.h"

#if defined(ENABLE_CATCH)

#    include <memory>
using std::shared_ptr;
#    include <catch2/catch_all.hpp>

#endif

namespace CoolProp {

// The anonymous union (max_abs_error vs. {using_tau_r, reducing_value, T_r, N})
// is intentionally written one active member per branch; the `type` discriminant
// gates every read, so the inactive members are never accessed. clang-tidy's
// union-access / member-init guidance does not model discriminated unions, so it
// is suppressed across this constructor.
// NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,cppcoreguidelines-pro-type-union-access)
SaturationAncillaryFunction::SaturationAncillaryFunction(const Values& v) : type(v.type), Tmin(v.Tmin), Tmax(v.Tmax) {
    if (type == TYPE_RATIONAL_POLYNOMIAL) {
        num_coeffs = v.num_coeffs;
        den_coeffs = v.den_coeffs;
        max_abs_error = v.max_abs_error;
    } else {
        n = v.n;
        t = v.t;
        N = n.size();
        s = n;
        reducing_value = v.reducing_value;
        using_tau_r = v.using_tau_r;
        T_r = v.T_r;
    }
}
// NOLINTEND(cppcoreguidelines-pro-type-member-init,cppcoreguidelines-pro-type-union-access)

double SaturationAncillaryFunction::evaluate(double T) {
    if (type == TYPE_NOT_SET) {
        throw ValueError(format("type not set"));
    } else if (type == TYPE_RATIONAL_POLYNOMIAL) {
        Polynomial2D poly;
        return poly.evaluate(num_coeffs, T) / poly.evaluate(den_coeffs, T);
    } else {
        double THETA = 1 - T / T_r;
        // Saturation ancillaries are only defined at or below the reducing
        // temperature. For T > T_r, pow(THETA, t) with non-integer t was
        // pow(negative, fractional), producing NaN and (depending on FP
        // trap settings) sometimes a SIGFPE that escaped C++ try/catch
        // (#1611). Many internal callers (mixture critical-point search,
        // VLE initialisers, etc.) legitimately probe ancillaries at T
        // ranges where the pure component is supercritical and rely on
        // the result being well-behaved-but-invalid (NaN) rather than an
        // exception. Return NaN explicitly so the caller's existing
        // sanity check / fallback runs without a SIGFPE.
        if (THETA < 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        for (std::size_t i = 0; i < N; ++i) {
            s[i] = n[i] * pow(THETA, t[i]);
        }
        double summer = std::accumulate(s.begin(), s.end(), 0.0);

        if (type == TYPE_NOT_EXPONENTIAL) {
            return reducing_value * (1 + summer);
        } else {
            double tau_r_value = NAN;
            if (using_tau_r)
                tau_r_value = T_r / T;
            else
                tau_r_value = 1.0;
            return reducing_value * exp(tau_r_value * summer);
        }
    }
}
double SaturationAncillaryFunction::invert(double value, double min_bound, double max_bound) {
    // Invert the ancillary curve to get the temperature as a function of the output variable
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
       public:
        SaturationAncillaryFunction* anc;
        CoolPropDbl value;

        solver_resid(SaturationAncillaryFunction* anc, CoolPropDbl value) : anc(anc), value(value) {}

        double call(double T) override {
            CoolPropDbl current_value = anc->evaluate(T);
            return current_value - value;
        }
    };
    solver_resid resid(this, value);
    if (min_bound < 0) {
        min_bound = Tmin - 0.01;
    }
    if (max_bound < 0) {
        max_bound = Tmax;
    }

    try {
        // Safe to expand the domain a little bit to lower temperature, absolutely cannot exceed Tmax
        // because then you get (negative number)^(double) which is undefined.
        return Brent(resid, min_bound, max_bound, DBL_EPSILON, 1e-10, 100);
    } catch (...) {
        return ExtrapolatingSecant(resid, max_bound, -0.01, 1e-12, 100);
    }
}

void MeltingLineVariables::set_limits() {
    if (type == MELTING_LINE_SIMON_TYPE) {

        // Fill in the min and max pressures for each part
        for (auto& part : simon.parts) {
            part.p_min = part.p_0 + part.a * (pow(part.T_min / part.T_0, part.c) - 1);
            part.p_max = part.p_0 + part.a * (pow(part.T_max / part.T_0, part.c) - 1);
        }
        pmin = simon.parts.front().p_min;
        pmax = simon.parts.back().p_max;
        Tmin = simon.parts.front().T_min;
        Tmax = simon.parts.back().T_max;
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
        // Fill in the min and max pressures for each part
        for (auto& part : polynomial_in_Tr.parts) {
            part.p_min = part.evaluate(part.T_min);
            part.p_max = part.evaluate(part.T_max);
        }
        Tmin = polynomial_in_Tr.parts.front().T_min;
        pmin = polynomial_in_Tr.parts.front().p_min;
        Tmax = polynomial_in_Tr.parts.back().T_max;
        pmax = polynomial_in_Tr.parts.back().p_max;
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {
        // Fill in the min and max pressures for each part
        for (auto& part : polynomial_in_Theta.parts) {
            part.p_min = part.evaluate(part.T_min);
            part.p_max = part.evaluate(part.T_max);
        }
        Tmin = polynomial_in_Theta.parts.front().T_min;
        pmin = polynomial_in_Theta.parts.front().p_min;
        Tmax = polynomial_in_Theta.parts.back().T_max;
        pmax = polynomial_in_Theta.parts.back().p_max;
    } else {
        throw ValueError("only Simon supported now");
    }
}

CoolPropDbl MeltingLineVariables::evaluate(int OF, int GIVEN, CoolPropDbl value) {
    if (type == MELTING_LINE_NOT_SET) {
        throw ValueError("Melting line curve not set");
    }
    if (OF == iP_max) {
        return pmax;
    } else if (OF == iP_min) {
        return pmin;
    } else if (OF == iT_max) {
        return Tmax;
    } else if (OF == iT_min) {
        return Tmin;
    } else if (OF == iP && GIVEN == iT) {
        CoolPropDbl T = value;
        if (type == MELTING_LINE_SIMON_TYPE) {
            // Need to find the right segment
            for (auto& part : simon.parts) {
                if (is_in_closed_range(part.T_min, part.T_max, T)) {
                    return part.p_0 + part.a * (pow(T / part.T_0, part.c) - 1);
                }
            }
            throw ValueError("unable to calculate melting line (p,T) for Simon curve");
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
            // Need to find the right segment
            for (auto& part : polynomial_in_Tr.parts) {
                if (is_in_closed_range(part.T_min, part.T_max, T)) {
                    return part.evaluate(T);
                }
            }
            throw ValueError("unable to calculate melting line (p,T) for polynomial_in_Tr curve");
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {
            // Need to find the right segment
            for (auto& part : polynomial_in_Theta.parts) {
                if (is_in_closed_range(part.T_min, part.T_max, T)) {
                    return part.evaluate(T);
                }
            }
            throw ValueError("unable to calculate melting line (p,T) for polynomial_in_Theta curve");
        } else {
            throw ValueError(format("Invalid melting line type [%d]", type));
        }
    } else {
        if (type == MELTING_LINE_SIMON_TYPE) {
            // Need to find the right segment
            for (auto& part : simon.parts) {
                //  p = part.p_0 + part.a*(pow(T/part.T_0,part.c)-1);
                CoolPropDbl T = pow((value - part.p_0) / part.a + 1, 1 / part.c) * part.T_0;
                if (get_config_bool(DONT_CHECK_PROPERTY_LIMITS) || (T >= part.T_0 && T <= part.T_max)) {
                    return T;
                }
            }
            throw ValueError(format("unable to calculate melting line T(p) for Simon curve for p=%Lg; bounds are %Lg,%Lg Pa", value, pmin, pmax));
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
            class solver_resid : public FuncWrapper1D
            {
               public:
                MeltingLinePiecewisePolynomialInTrSegment* part;
                CoolPropDbl given_p;
                solver_resid(MeltingLinePiecewisePolynomialInTrSegment* part, CoolPropDbl p) : part(part), given_p(p) {};
                double call(double T) override {

                    CoolPropDbl calc_p = part->evaluate(T);

                    // Difference between the two is to be driven to zero
                    return given_p - calc_p;
                };
            };

            // Need to find the right segment
            for (auto& part : polynomial_in_Tr.parts) {
                if (is_in_closed_range(part.p_min, part.p_max, value)) {
                    solver_resid resid(&part, value);
                    double T = Brent(resid, part.T_min, part.T_max, DBL_EPSILON, 1e-12, 100);
                    return T;
                }
            }
            throw ValueError(
              format("unable to calculate melting line T(p) for polynomial_in_Theta curve for p=%Lg; bounds are %Lg,%Lg Pa", value, pmin, pmax));
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {

            class solver_resid : public FuncWrapper1D
            {
               public:
                MeltingLinePiecewisePolynomialInThetaSegment* part;
                CoolPropDbl given_p;
                solver_resid(MeltingLinePiecewisePolynomialInThetaSegment* part, CoolPropDbl p) : part(part), given_p(p) {};
                double call(double T) override {

                    CoolPropDbl calc_p = part->evaluate(T);

                    // Difference between the two is to be driven to zero
                    return given_p - calc_p;
                };
            };

            // Need to find the right segment
            for (auto& part : polynomial_in_Theta.parts) {
                if (is_in_closed_range(part.p_min, part.p_max, value)) {
                    solver_resid resid(&part, value);
                    double T = Brent(resid, part.T_min, part.T_max, DBL_EPSILON, 1e-12, 100);
                    return T;
                }
            }

            throw ValueError(
              format("unable to calculate melting line T(p) for polynomial_in_Theta curve for p=%Lg; bounds are %Lg,%Lg Pa", value, pmin, pmax));
        } else {
            throw ValueError(format("Invalid melting line type T(p) [%d]", type));
        }
    }
}

std::vector<std::pair<CoolPropDbl, CoolPropDbl>> MeltingLineVariables::get_parts_pranges() const {
    std::vector<std::pair<CoolPropDbl, CoolPropDbl>> out;
    if (type == MELTING_LINE_SIMON_TYPE) {
        for (const auto& p : simon.parts)
            out.emplace_back(p.p_min, p.p_max);
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
        for (const auto& p : polynomial_in_Tr.parts)
            out.emplace_back(p.p_min, p.p_max);
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {
        for (const auto& p : polynomial_in_Theta.parts)
            out.emplace_back(p.p_min, p.p_max);
    }
    return out;
}

}; /* namespace CoolProp */

#if defined(ENABLE_CATCH)
TEST_CASE("Water melting line", "[melting]") {
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "water"));
    int iT = CoolProp::iT, iP = CoolProp::iP;
    SECTION("Ice Ih-liquid") {
        double actual = AS->melting_line(iT, iP, 138.268e6);
        double expected = 260.0;
        CAPTURE(actual);
        CAPTURE(expected);
        CHECK(std::abs(actual - expected) < 0.01);
    }
    SECTION("Ice III-liquid") {
        double actual = AS->melting_line(iT, iP, 268.685e6);
        double expected = 254;
        CAPTURE(actual);
        CAPTURE(expected);
        CHECK(std::abs(actual - expected) < 0.01);
    }
    SECTION("Ice V-liquid") {
        double actual = AS->melting_line(iT, iP, 479.640e6);
        double expected = 265;
        CAPTURE(actual);
        CAPTURE(expected);
        CHECK(std::abs(actual - expected) < 0.01);
    }
    SECTION("Ice VI-liquid") {
        double actual = AS->melting_line(iT, iP, 1356.76e6);
        double expected = 320;
        CAPTURE(actual);
        CAPTURE(expected);
        CHECK(std::abs(actual - expected) < 1);
    }
}

TEST_CASE("Tests for values from melting lines", "[melting]") {
    int iT = CoolProp::iT, iP = CoolProp::iP;
    std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    for (const auto& fluid : fluids) {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));

        // Skip fluids without a melting line. Also skip Water and HeavyWater:
        // their melting curves are anomalous (the ice-Ih branch folds back BELOW
        // the triple temperature), which violates this generic test's
        // actual_T > Tmin / actual_T < Tmax bounds. Both have dedicated
        // check-value tests instead (see the "[melting]" cases in CoolProp-Tests).
        if (!AS->has_melting_line() || !fluid.compare("Water") || !fluid.compare("HeavyWater")) {
            continue;
        }

        double pmax = AS->melting_line(CoolProp::iP_max, iT, 0);
        double pmin = AS->melting_line(CoolProp::iP_min, iT, 0);
        double Tmax = AS->melting_line(CoolProp::iT_max, iT, 0);
        double Tmin = AS->melting_line(CoolProp::iT_min, iT, 0);

        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss0;
        ss0 << "Check melting line limits for fluid " << fluid;
        SECTION(ss0.str(), "") {
            CAPTURE(Tmin);
            CAPTURE(Tmax);
            CAPTURE(pmin);
            CAPTURE(pmax);
            CHECK(Tmax > Tmin);
            CHECK(pmax > pmin);
            CHECK(pmin > 0);
        }
        // Integer-indexed grid (cert-flp30-c): the original
        //   for (p = 0.1*range + pmin; p < pmax; p += 0.2*range)
        // intends 5 samples at p_i = pmin + (0.1 + 0.2*i)*range for
        // i in {0..4}; preserve that count exactly.
        const double p_range = pmax - pmin;
        for (int p_i = 0; p_i < 5; ++p_i) {
            const double p = pmin + (0.1 + 0.2 * p_i) * p_range;
            // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
            std::ostringstream ss1;
            ss1 << "Melting line for " << fluid << " at p=" << p;
            SECTION(ss1.str(), "") {
                double actual_T = AS->melting_line(iT, iP, p);
                CAPTURE(Tmin);
                CAPTURE(Tmax);
                CAPTURE(actual_T);
                CHECK(actual_T > Tmin);
                CHECK(actual_T < Tmax);
            }
        }
        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss2;
        ss2 << "Ensure melting line valid for " << fluid << " @ EOS pmax";
        SECTION(ss2.str(), "") {
            double actual_T = -_HUGE;
            double EOS_pmax = AS->pmax();
            double T_pmax_required = -1;
            try {
                CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, true);
                T_pmax_required = AS->melting_line(iT, iP, EOS_pmax);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // Best-effort probe: T_pmax_required stays at its -1
                // sentinel and shows up in the CAPTURE below.  The real
                // assertion is CHECK_NOTHROW on the next melting_line
                // call.
            }
            CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, false);
            CAPTURE(T_pmax_required);
            CAPTURE(EOS_pmax);
            CHECK_NOTHROW(actual_T = AS->melting_line(iT, iP, EOS_pmax));
            CAPTURE(actual_T);
        }
    }
}

TEST_CASE("Test that hs_anchor enthalpy/entropy agrees with EOS", "[ancillaries]") {
    std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    for (const auto& fluid : fluids) {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));

        CoolProp::SimpleState hs_anchor = AS->get_state("hs_anchor");

        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "Check hs_anchor for " << fluid;
        SECTION(ss1.str(), "") {
            INFO("The enthalpy and entropy are hardcoded in the fluid JSON files.  They MUST agree with the values calculated by the EOS");
            AS->update(CoolProp::DmolarT_INPUTS, hs_anchor.rhomolar, hs_anchor.T);
            double EOS_hmolar = AS->hmolar();
            double EOS_smolar = AS->smolar();
            CAPTURE(hs_anchor.hmolar);
            CAPTURE(hs_anchor.smolar);
            CAPTURE(EOS_hmolar);
            CAPTURE(EOS_smolar);
            CHECK(std::abs(EOS_hmolar - hs_anchor.hmolar) < 1e-3);
            CHECK(std::abs(EOS_smolar - hs_anchor.smolar) < 1e-3);
        }
    }
}

TEST_CASE("Surface tension", "[surface_tension]") {
    SECTION("from PropsSI") {
        CHECK(ValidNumber(CoolProp::PropsSI("surface_tension", "T", 300, "Q", 0, "Water")));
    }
    SECTION("from saturation_ancillary") {
        CHECK(ValidNumber(CoolProp::saturation_ancillary("Water", "surface_tension", 0, "T", 300)));
    }
    SECTION("from AbstractState") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(CoolProp::QT_INPUTS, 0, 300);
        CHECK_NOTHROW(AS->surface_tension());
    }
}
#endif
