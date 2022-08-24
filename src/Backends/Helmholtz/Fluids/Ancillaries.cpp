#include "Ancillaries.h"
#include "DataStructures.h"
#include "AbstractState.h"

#if defined(ENABLE_CATCH)

#    include "crossplatform_shared_ptr.h"
#    include <catch2/catch_all.hpp>

#endif

namespace CoolProp {

SaturationAncillaryFunction::SaturationAncillaryFunction(rapidjson::Value& json_code) {
    std::string type = cpjson::get_string(json_code, "type");
    if (!type.compare("rational_polynomial")) {
        this->type = TYPE_RATIONAL_POLYNOMIAL;
        num_coeffs = vec_to_eigen(cpjson::get_double_array(json_code["A"]));
        den_coeffs = vec_to_eigen(cpjson::get_double_array(json_code["B"]));
        max_abs_error = cpjson::get_double(json_code, "max_abs_error");
        try {
            Tmin = cpjson::get_double(json_code, "Tmin");
            Tmax = cpjson::get_double(json_code, "Tmax");
        } catch (...) {
            Tmin = _HUGE;
            Tmax = _HUGE;
        }
    } else {
        if (!type.compare("rhoLnoexp"))
            this->type = TYPE_NOT_EXPONENTIAL;
        else
            this->type = TYPE_EXPONENTIAL;
        n = cpjson::get_double_array(json_code["n"]);
        N = n.size();
        s = n;
        t = cpjson::get_double_array(json_code["t"]);
        Tmin = cpjson::get_double(json_code, "Tmin");
        Tmax = cpjson::get_double(json_code, "Tmax");
        reducing_value = cpjson::get_double(json_code, "reducing_value");
        using_tau_r = cpjson::get_bool(json_code, "using_tau_r");
        T_r = cpjson::get_double(json_code, "T_r");
    }
};

double SaturationAncillaryFunction::evaluate(double T) {
    if (type == TYPE_NOT_SET) {
        throw ValueError(format("type not set"));
    } else if (type == TYPE_RATIONAL_POLYNOMIAL) {
        Polynomial2D poly;
        return poly.evaluate(num_coeffs, T) / poly.evaluate(den_coeffs, T);
    } else {
        double THETA = 1 - T / T_r;

        for (std::size_t i = 0; i < N; ++i) {
            s[i] = n[i] * pow(THETA, t[i]);
        }
        double summer = std::accumulate(s.begin(), s.end(), 0.0);

        if (type == TYPE_NOT_EXPONENTIAL) {
            return reducing_value * (1 + summer);
        } else {
            double tau_r_value;
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

        double call(double T) {
            CoolPropDbl current_value = anc->evaluate(T);
            return current_value - value;
        }
    };
    solver_resid resid(this, value);
    std::string errstring;
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
        return Secant(resid, max_bound, -0.01, 1e-12, 100);
    }
}

void MeltingLineVariables::set_limits(void) {
    if (type == MELTING_LINE_SIMON_TYPE) {

        // Fill in the min and max pressures for each part
        for (std::size_t i = 0; i < simon.parts.size(); ++i) {
            MeltingLinePiecewiseSimonSegment& part = simon.parts[i];
            part.p_min = part.p_0 + part.a * (pow(part.T_min / part.T_0, part.c) - 1);
            part.p_max = part.p_0 + part.a * (pow(part.T_max / part.T_0, part.c) - 1);
        }
        pmin = simon.parts.front().p_min;
        pmax = simon.parts.back().p_max;
        Tmin = simon.parts.front().T_min;
        Tmax = simon.parts.back().T_max;
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
        // Fill in the min and max pressures for each part
        for (std::size_t i = 0; i < polynomial_in_Tr.parts.size(); ++i) {
            MeltingLinePiecewisePolynomialInTrSegment& part = polynomial_in_Tr.parts[i];
            part.p_min = part.evaluate(part.T_min);
            part.p_max = part.evaluate(part.T_max);
        }
        Tmin = polynomial_in_Tr.parts.front().T_min;
        pmin = polynomial_in_Tr.parts.front().p_min;
        Tmax = polynomial_in_Tr.parts.back().T_max;
        pmax = polynomial_in_Tr.parts.back().p_max;
    } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {
        // Fill in the min and max pressures for each part
        for (std::size_t i = 0; i < polynomial_in_Theta.parts.size(); ++i) {
            MeltingLinePiecewisePolynomialInThetaSegment& part = polynomial_in_Theta.parts[i];
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
            for (std::size_t i = 0; i < simon.parts.size(); ++i) {
                MeltingLinePiecewiseSimonSegment& part = simon.parts[i];
                if (is_in_closed_range(part.T_min, part.T_max, T)) {
                    return part.p_0 + part.a * (pow(T / part.T_0, part.c) - 1);
                }
            }
            throw ValueError("unable to calculate melting line (p,T) for Simon curve");
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE) {
            // Need to find the right segment
            for (std::size_t i = 0; i < polynomial_in_Tr.parts.size(); ++i) {
                MeltingLinePiecewisePolynomialInTrSegment& part = polynomial_in_Tr.parts[i];
                if (is_in_closed_range(part.T_min, part.T_max, T)) {
                    return part.evaluate(T);
                }
            }
            throw ValueError("unable to calculate melting line (p,T) for polynomial_in_Tr curve");
        } else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE) {
            // Need to find the right segment
            for (std::size_t i = 0; i < polynomial_in_Theta.parts.size(); ++i) {
                MeltingLinePiecewisePolynomialInThetaSegment& part = polynomial_in_Theta.parts[i];
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
            for (std::size_t i = 0; i < simon.parts.size(); ++i) {
                MeltingLinePiecewiseSimonSegment& part = simon.parts[i];
                //  p = part.p_0 + part.a*(pow(T/part.T_0,part.c)-1);
                CoolPropDbl T = pow((value - part.p_0) / part.a + 1, 1 / part.c) * part.T_0;
                if (T >= part.T_0 && T <= part.T_max) {
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
                solver_resid(MeltingLinePiecewisePolynomialInTrSegment* part, CoolPropDbl p) : part(part), given_p(p){};
                double call(double T) {

                    CoolPropDbl calc_p = part->evaluate(T);

                    // Difference between the two is to be driven to zero
                    return given_p - calc_p;
                };
            };

            // Need to find the right segment
            for (std::size_t i = 0; i < polynomial_in_Tr.parts.size(); ++i) {
                MeltingLinePiecewisePolynomialInTrSegment& part = polynomial_in_Tr.parts[i];
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
                solver_resid(MeltingLinePiecewisePolynomialInThetaSegment* part, CoolPropDbl p) : part(part), given_p(p){};
                double call(double T) {

                    CoolPropDbl calc_p = part->evaluate(T);

                    // Difference between the two is to be driven to zero
                    return given_p - calc_p;
                };
            };

            // Need to find the right segment
            for (std::size_t i = 0; i < polynomial_in_Theta.parts.size(); ++i) {
                MeltingLinePiecewisePolynomialInThetaSegment& part = polynomial_in_Theta.parts[i];
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
    for (std::size_t i = 0; i < fluids.size(); ++i) {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluids[i]));

        // Water has its own better tests; skip fluids without melting line
        if (!AS->has_melting_line() || !fluids[i].compare("Water")) {
            continue;
        }

        double pmax = AS->melting_line(CoolProp::iP_max, iT, 0);
        double pmin = AS->melting_line(CoolProp::iP_min, iT, 0);
        double Tmax = AS->melting_line(CoolProp::iT_max, iT, 0);
        double Tmin = AS->melting_line(CoolProp::iT_min, iT, 0);

        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss0;
        ss0 << "Check melting line limits for fluid " << fluids[i];
        SECTION(ss0.str(), "") {
            CAPTURE(Tmin);
            CAPTURE(Tmax);
            CAPTURE(pmin);
            CAPTURE(pmax);
            CHECK(Tmax > Tmin);
            CHECK(pmax > pmin);
            CHECK(pmin > 0);
        }
        for (double p = 0.1 * (pmax - pmin) + pmin; p < pmax; p += 0.2 * (pmax - pmin)) {
            // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
            std::ostringstream ss1;
            ss1 << "Melting line for " << fluids[i] << " at p=" << p;
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
        ss2 << "Ensure melting line valid for " << fluids[i] << " @ EOS pmax";
        SECTION(ss2.str(), "") {
            double actual_T;
            double EOS_pmax = AS->pmax();
            CAPTURE(EOS_pmax);
            CHECK_NOTHROW(actual_T = AS->melting_line(iT, iP, EOS_pmax));
            CAPTURE(actual_T);
        }
    }
}

TEST_CASE("Test that hs_anchor enthalpy/entropy agrees with EOS", "[ancillaries]") {
    std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    for (std::size_t i = 0; i < fluids.size(); ++i) {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluids[i]));

        CoolProp::SimpleState hs_anchor = AS->get_state("hs_anchor");

        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "Check hs_anchor for " << fluids[i];
        SECTION(ss1.str(), "") {
            std::string note =
              "The enthalpy and entropy are hardcoded in the fluid JSON files.  They MUST agree with the values calculated by the EOS";
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
