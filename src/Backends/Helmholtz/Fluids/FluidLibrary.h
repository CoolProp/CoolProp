
#ifndef FLUIDLIBRARY_H
#define FLUIDLIBRARY_H

#include "CoolPropFluid.h"

#include "rapidjson_include.h"
#include "crossplatform_shared_ptr.h"

#include <map>
#include <algorithm>
#include "Configuration.h"
#include "Backends/Cubics/CubicsLibrary.h"
#include "Helmholtz.h"

namespace CoolProp {

// Forward declaration of the necessary debug function to avoid including the whole header
extern int get_debug_level();

/// A container for the fluid parameters for the CoolProp fluids
/**
This container holds copies of all of the fluid instances for the fluids that are loaded in CoolProp.
New fluids can be added by passing in a rapidjson::Value instance to the add_one function, or
a rapidjson array of fluids to the add_many function.
*/
class JSONFluidLibrary
{
    /// Map from CAS code to JSON instance.  For pseudo-pure fluids, use name in place of CAS code since no CAS number is defined for mixtures
    std::map<std::size_t, CoolPropFluid> fluid_map;
    /// Map from index of fluid to a string
    std::map<std::size_t, std::string> JSONstring_map;
    std::vector<std::string> name_vector;
    std::map<std::string, std::size_t> string_to_index_map;
    bool _is_empty;

   public:
    /// Parse the contributions to the residual Helmholtz energy
    static ResidualHelmholtzContainer parse_alphar(rapidjson::Value& jsonalphar) {
        ResidualHelmholtzContainer alphar;

        for (rapidjson::Value::ValueIterator itr = jsonalphar.Begin(); itr != jsonalphar.End(); ++itr) {
            // A reference for code cleanness
            rapidjson::Value& contribution = *itr;

            // Get the type (required!)
            std::string type = contribution["type"].GetString();

            if (!type.compare("ResidualHelmholtzPower")) {
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<CoolPropDbl> l = cpjson::get_long_double_array(contribution["l"]);
                assert(n.size() == d.size());
                assert(n.size() == t.size());
                assert(n.size() == l.size());

                alphar.GenExp.add_Power(n, d, t, l);
            } else if (!type.compare("ResidualHelmholtzGaussian")) {
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<CoolPropDbl> eta = cpjson::get_long_double_array(contribution["eta"]);
                std::vector<CoolPropDbl> epsilon = cpjson::get_long_double_array(contribution["epsilon"]);
                std::vector<CoolPropDbl> beta = cpjson::get_long_double_array(contribution["beta"]);
                std::vector<CoolPropDbl> gamma = cpjson::get_long_double_array(contribution["gamma"]);
                assert(n.size() == d.size());
                assert(n.size() == t.size());
                assert(n.size() == eta.size());
                assert(n.size() == epsilon.size());
                assert(n.size() == beta.size());
                assert(n.size() == gamma.size());
                alphar.GenExp.add_Gaussian(n, d, t, eta, epsilon, beta, gamma);
            } else if (!type.compare("ResidualHelmholtzGaoB")) {
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<CoolPropDbl> eta = cpjson::get_long_double_array(contribution["eta"]);
                std::vector<CoolPropDbl> beta = cpjson::get_long_double_array(contribution["beta"]);
                std::vector<CoolPropDbl> gamma = cpjson::get_long_double_array(contribution["gamma"]);
                std::vector<CoolPropDbl> epsilon = cpjson::get_long_double_array(contribution["epsilon"]);
                std::vector<CoolPropDbl> b = cpjson::get_long_double_array(contribution["b"]);
                assert(n.size() == t.size());
                assert(n.size() == d.size());
                assert(n.size() == eta.size());
                assert(n.size() == epsilon.size());
                assert(n.size() == beta.size());
                assert(n.size() == gamma.size());
                assert(n.size() == b.size());
                alphar.GaoB = ResidualHelmholtzGaoB(n, t, d, eta, beta, gamma, epsilon, b);
            } else if (!type.compare("ResidualHelmholtzNonAnalytic")) {
                if (alphar.NonAnalytic.N > 0) {
                    throw ValueError("Cannot add ");
                }
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> a = cpjson::get_long_double_array(contribution["a"]);
                std::vector<CoolPropDbl> b = cpjson::get_long_double_array(contribution["b"]);
                std::vector<CoolPropDbl> beta = cpjson::get_long_double_array(contribution["beta"]);
                std::vector<CoolPropDbl> A = cpjson::get_long_double_array(contribution["A"]);
                std::vector<CoolPropDbl> B = cpjson::get_long_double_array(contribution["B"]);
                std::vector<CoolPropDbl> C = cpjson::get_long_double_array(contribution["C"]);
                std::vector<CoolPropDbl> D = cpjson::get_long_double_array(contribution["D"]);
                assert(n.size() == a.size());
                assert(n.size() == b.size());
                assert(n.size() == beta.size());
                assert(n.size() == A.size());
                assert(n.size() == B.size());
                assert(n.size() == C.size());
                assert(n.size() == D.size());
                alphar.NonAnalytic = ResidualHelmholtzNonAnalytic(n, a, b, beta, A, B, C, D);
            } else if (!type.compare("ResidualHelmholtzLemmon2005")) {
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<CoolPropDbl> l = cpjson::get_long_double_array(contribution["l"]);
                std::vector<CoolPropDbl> m = cpjson::get_long_double_array(contribution["m"]);
                assert(n.size() == d.size());
                assert(n.size() == t.size());
                assert(n.size() == l.size());
                assert(n.size() == m.size());
                alphar.GenExp.add_Lemmon2005(n, d, t, l, m);
            } else if (!type.compare("ResidualHelmholtzExponential")) {
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<CoolPropDbl> g = cpjson::get_long_double_array(contribution["g"]);
                std::vector<CoolPropDbl> l = cpjson::get_long_double_array(contribution["l"]);
                assert(n.size() == d.size());
                assert(n.size() == t.size());
                assert(n.size() == g.size());
                assert(n.size() == l.size());
                alphar.GenExp.add_Exponential(n, d, t, g, l);
            } else if (!type.compare("ResidualHelmholtzAssociating")) {
                if (alphar.SAFT.disabled == false) {
                    throw ValueError("Cannot add ");
                }
                CoolPropDbl a = cpjson::get_double(contribution, "a");
                CoolPropDbl m = cpjson::get_double(contribution, "m");
                CoolPropDbl epsilonbar = cpjson::get_double(contribution, "epsilonbar");
                CoolPropDbl vbarn = cpjson::get_double(contribution, "vbarn");
                CoolPropDbl kappabar = cpjson::get_double(contribution, "kappabar");
                alphar.SAFT = ResidualHelmholtzSAFTAssociating(a, m, epsilonbar, vbarn, kappabar);
            } else {
                throw ValueError(format("Unsupported Residual helmholtz type: %s", type.c_str()));
            }
        }

        // Finish adding parts to the Generalized Exponential term, build other vectors
        alphar.GenExp.finish();

        return alphar;
    };

    /// Parse the contributions to the ideal-gas Helmholtz energy
    static IdealHelmholtzContainer parse_alpha0(rapidjson::Value& jsonalpha0) {
        if (!jsonalpha0.IsArray()) {
            throw ValueError();
        }

        IdealHelmholtzContainer alpha0;

        for (rapidjson::Value::ConstValueIterator itr = jsonalpha0.Begin(); itr != jsonalpha0.End(); ++itr) {
            // A reference for code cleanness
            const rapidjson::Value& contribution = *itr;

            // Get the type (required!)
            std::string type = contribution["type"].GetString();

            if (!type.compare("IdealGasHelmholtzLead")) {
                if (alpha0.Lead.is_enabled() == true) {
                    throw ValueError("Cannot add ");
                }
                CoolPropDbl a1 = cpjson::get_double(contribution, "a1");
                CoolPropDbl a2 = cpjson::get_double(contribution, "a2");

                alpha0.Lead = IdealHelmholtzLead(a1, a2);
            } else if (!type.compare("IdealGasHelmholtzPower")) {
                if (alpha0.Power.is_enabled() == true) {
                    throw ValueError("Cannot add ");
                }
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);

                alpha0.Power = IdealHelmholtzPower(n, t);
            } else if (!type.compare("IdealGasHelmholtzLogTau")) {
                if (alpha0.LogTau.is_enabled() == true) {
                    throw ValueError("Cannot add ");
                }
                CoolPropDbl a = cpjson::get_double(contribution, "a");

                alpha0.LogTau = IdealHelmholtzLogTau(a);
            } else if (!type.compare("IdealGasHelmholtzPlanckEinsteinGeneralized")) {
                // Retrieve the values
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);

                std::vector<CoolPropDbl> c = cpjson::get_long_double_array(contribution["c"]);
                std::vector<CoolPropDbl> d = cpjson::get_long_double_array(contribution["d"]);

                if (alpha0.PlanckEinstein.is_enabled() == true) {
                    alpha0.PlanckEinstein.extend(n, t, c, d);
                } else {
                    alpha0.PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
                }
            } else if (!type.compare("IdealGasHelmholtzPlanckEinstein")) {
                // Retrieve the values
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                // Flip the sign of theta
                for (std::size_t i = 0; i < t.size(); ++i) {
                    t[i] *= -1;
                }
                std::vector<CoolPropDbl> c(n.size(), 1);
                std::vector<CoolPropDbl> d(c.size(), -1);

                if (alpha0.PlanckEinstein.is_enabled() == true) {
                    alpha0.PlanckEinstein.extend(n, t, c, d);
                } else {
                    alpha0.PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
                }
            } else if (!type.compare("IdealGasHelmholtzPlanckEinsteinFunctionT")) {
                // Retrieve the values
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> v = cpjson::get_long_double_array(contribution["v"]), theta(n.size(), 0.0);
                // Calculate theta
                double Tc = cpjson::get_double(contribution, "Tcrit");
                for (std::size_t i = 0; i < v.size(); ++i) {
                    theta[i] = -v[i] / Tc;
                }
                std::vector<CoolPropDbl> c(n.size(), 1);
                std::vector<CoolPropDbl> d(c.size(), -1);

                if (alpha0.PlanckEinstein.is_enabled() == true) {
                    alpha0.PlanckEinstein.extend(n, theta, c, d);
                } else {
                    alpha0.PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, theta, c, d);
                }
            } else if (!type.compare("IdealGasHelmholtzGERG2004Cosh")) {
                // Retrieve the values
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> theta = cpjson::get_long_double_array(contribution["theta"]);
                double Tc = cpjson::get_double(contribution, "Tcrit");
                if (alpha0.GERG2004Cosh.is_enabled() == true) {
                    alpha0.GERG2004Cosh.extend(n, theta);
                } else {
                    alpha0.GERG2004Cosh = IdealHelmholtzGERG2004Cosh(n, theta, Tc);
                }
            } else if (!type.compare("IdealGasHelmholtzGERG2004Sinh")) {
                // Retrieve the values
                std::vector<CoolPropDbl> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<CoolPropDbl> theta = cpjson::get_long_double_array(contribution["theta"]);
                double Tc = cpjson::get_double(contribution, "Tcrit");
                if (alpha0.GERG2004Sinh.is_enabled() == true) {
                    alpha0.GERG2004Sinh.extend(n, theta);
                } else {
                    alpha0.GERG2004Sinh = IdealHelmholtzGERG2004Sinh(n, theta, Tc);
                }
            } else if (!type.compare("IdealGasHelmholtzCP0Constant")) {
                if (alpha0.CP0Constant.is_enabled() == true) {
                    throw ValueError("Cannot add another IdealGasHelmholtzCP0Constant term; join them together");
                }
                CoolPropDbl cp_over_R = cpjson::get_double(contribution, "cp_over_R");
                CoolPropDbl Tc = cpjson::get_double(contribution, "Tc");
                CoolPropDbl T0 = cpjson::get_double(contribution, "T0");
                alpha0.CP0Constant = IdealHelmholtzCP0Constant(cp_over_R, Tc, T0);
            } else if (!type.compare("IdealGasHelmholtzCP0PolyT")) {
                if (alpha0.CP0PolyT.is_enabled() == true) {
                    throw ValueError("Cannot add another CP0PolyT term; join them together");
                }
                std::vector<CoolPropDbl> c = cpjson::get_long_double_array(contribution["c"]);
                std::vector<CoolPropDbl> t = cpjson::get_long_double_array(contribution["t"]);
                CoolPropDbl Tc = cpjson::get_double(contribution, "Tc");
                CoolPropDbl T0 = cpjson::get_double(contribution, "T0");
                alpha0.CP0PolyT = IdealHelmholtzCP0PolyT(c, t, Tc, T0);
            } else if (!type.compare("IdealGasHelmholtzCP0AlyLee")) {

                std::vector<CoolPropDbl> constants = cpjson::get_long_double_array(contribution["c"]);
                CoolPropDbl Tc = cpjson::get_double(contribution, "Tc");
                CoolPropDbl T0 = cpjson::get_double(contribution, "T0");

                // Take the constant term if nonzero and set it as a polyT term
                if (std::abs(constants[0]) > 1e-14) {
                    std::vector<CoolPropDbl> c(1, constants[0]), t(1, 0);
                    if (alpha0.CP0PolyT.is_enabled() == true) {
                        alpha0.CP0PolyT.extend(c, t);
                    } else {
                        alpha0.CP0PolyT = IdealHelmholtzCP0PolyT(c, t, Tc, T0);
                    }
                }
                std::vector<CoolPropDbl> n, c, d, t;
                if (std::abs(constants[1]) > 1e-14) {
                    // sinh term can be converted by setting  a_k = C, b_k = 2*D, c_k = -1, d_k = 1
                    n.push_back(constants[1]);
                    t.push_back(-2 * constants[2] / Tc);
                    c.push_back(1);
                    d.push_back(-1);
                }
                if (std::abs(constants[3]) > 1e-14) {
                    // cosh term can be converted by setting  a_k = C, b_k = 2*D, c_k = 1, d_k = 1
                    n.push_back(-constants[3]);
                    t.push_back(-2 * constants[4] / Tc);
                    c.push_back(1);
                    d.push_back(1);
                }
                if (alpha0.PlanckEinstein.is_enabled() == true) {
                    alpha0.PlanckEinstein.extend(n, t, c, d);
                } else {
                    alpha0.PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
                }
            } else if (!type.compare("IdealGasHelmholtzEnthalpyEntropyOffset")) {
                CoolPropDbl a1 = cpjson::get_double(contribution, "a1");
                CoolPropDbl a2 = cpjson::get_double(contribution, "a2");
                std::string reference = cpjson::get_string(contribution, "reference");
                alpha0.EnthalpyEntropyOffsetCore = IdealHelmholtzEnthalpyEntropyOffset(a1, a2, reference);
            } else {
                std::cout << format("Unsupported ideal-gas Helmholtz type: %s\n", type.c_str());
                //throw ValueError(format("Unsupported ideal-gas Helmholtz type: %s",type.c_str()));
            }
        }
        return alpha0;
    };

   protected:
    /// Parse the environmental parameters (ODP, GWP, etc.)
    void parse_environmental(rapidjson::Value& json, CoolPropFluid& fluid) {
        fluid.environment.ASHRAE34 = cpjson::get_string(json, "ASHRAE34");
        fluid.environment.GWP20 = cpjson::get_double(json, "GWP20");
        fluid.environment.GWP100 = cpjson::get_double(json, "GWP100");
        fluid.environment.GWP500 = cpjson::get_double(json, "GWP500");
        fluid.environment.HH = cpjson::get_double(json, "HH");
        fluid.environment.FH = cpjson::get_double(json, "FH");
        fluid.environment.PH = cpjson::get_double(json, "PH");
        fluid.environment.ODP = cpjson::get_double(json, "ODP");
    }

    /// Parse the Equation of state JSON entry
    void parse_EOS(rapidjson::Value& EOS_json, CoolPropFluid& fluid) {
        EquationOfState E;
        fluid.EOSVector.push_back(E);

        EquationOfState& EOS = fluid.EOSVector.at(fluid.EOSVector.size() - 1);

        // Universal gas constant [J/mol/K]
        EOS.R_u = cpjson::get_double(EOS_json, "gas_constant");
        EOS.molar_mass = cpjson::get_double(EOS_json, "molar_mass");
        EOS.acentric = cpjson::get_double(EOS_json, "acentric");

        EOS.pseudo_pure = cpjson::get_bool(EOS_json, "pseudo_pure");
        EOS.limits.Tmax = cpjson::get_double(EOS_json, "T_max");
        EOS.limits.pmax = cpjson::get_double(EOS_json, "p_max");

        rapidjson::Value& reducing_state = EOS_json["STATES"]["reducing"];
        rapidjson::Value& satminL_state = EOS_json["STATES"]["sat_min_liquid"];
        rapidjson::Value& satminV_state = EOS_json["STATES"]["sat_min_vapor"];

        // Reducing state
        EOS.reduce.T = cpjson::get_double(reducing_state, "T");
        EOS.reduce.rhomolar = cpjson::get_double(reducing_state, "rhomolar");
        EOS.reduce.p = cpjson::get_double(reducing_state, "p");
        EOS.reduce.hmolar = cpjson::get_double(reducing_state, "hmolar");
        EOS.reduce.smolar = cpjson::get_double(reducing_state, "smolar");

        EOS.sat_min_liquid.T = cpjson::get_double(satminL_state, "T");
        EOS.sat_min_liquid.p = cpjson::get_double(satminL_state, "p");
        EOS.sat_min_liquid.rhomolar = cpjson::get_double(satminL_state, "rhomolar");
        EOS.sat_min_vapor.T = cpjson::get_double(satminV_state, "T");
        EOS.sat_min_vapor.p = cpjson::get_double(satminV_state, "p");
        EOS.sat_min_vapor.rhomolar = cpjson::get_double(satminV_state, "rhomolar");

        /// \todo: define limits of EOS better
        EOS.limits.Tmin = cpjson::get_double(satminL_state, "T");
        EOS.ptriple = cpjson::get_double(satminL_state, "p");
        EOS.Ttriple = EOS.limits.Tmin;

        // BibTex keys
        EOS.BibTeX_EOS = cpjson::get_string(EOS_json, "BibTeX_EOS");
        EOS.BibTeX_CP0 = cpjson::get_string(EOS_json, "BibTeX_CP0");

        EOS.alphar = parse_alphar(EOS_json["alphar"]);
        EOS.alpha0 = parse_alpha0(EOS_json["alpha0"]);

        // Store the prefactor multipliying alpha0 if present
        if (EOS_json.HasMember("alpha0_prefactor")) {
            EOS.alpha0.set_prefactor(cpjson::get_double(EOS_json, "alpha0_prefactor"));
        }
        if (EOS_json["STATES"].HasMember("hs_anchor")) {
            rapidjson::Value& hs_anchor = EOS_json["STATES"]["hs_anchor"];
            EOS.hs_anchor.T = cpjson::get_double(hs_anchor, "T");
            EOS.hs_anchor.p = cpjson::get_double(hs_anchor, "p");
            EOS.hs_anchor.rhomolar = cpjson::get_double(hs_anchor, "rhomolar");
            EOS.hs_anchor.hmolar = cpjson::get_double(hs_anchor, "hmolar");
            EOS.hs_anchor.smolar = cpjson::get_double(hs_anchor, "smolar");
        }

        if (EOS_json["STATES"].HasMember("pressure_max_sat")) {
            rapidjson::Value& s = EOS_json["STATES"]["pressure_max_sat"];
            EOS.max_sat_p.T = cpjson::get_double(s, "T");
            EOS.max_sat_p.p = cpjson::get_double(s, "p");
            EOS.max_sat_p.rhomolar = cpjson::get_double(s, "rhomolar");
            if (s.HasMember("hmolar")) {
                EOS.max_sat_p.hmolar = cpjson::get_double(s, "hmolar");
                EOS.max_sat_p.smolar = cpjson::get_double(s, "smolar");
            }
        }

        if (EOS_json["STATES"].HasMember("temperature_max_sat")) {
            rapidjson::Value& s = EOS_json["STATES"]["temperature_max_sat"];
            EOS.max_sat_T.T = cpjson::get_double(s, "T");
            EOS.max_sat_T.p = cpjson::get_double(s, "p");
            EOS.max_sat_T.rhomolar = cpjson::get_double(s, "rhomolar");
            if (s.HasMember("hmolar")) {
                EOS.max_sat_T.hmolar = cpjson::get_double(s, "hmolar");
                EOS.max_sat_T.smolar = cpjson::get_double(s, "smolar");
            }
        }

        if (EOS_json.HasMember("critical_region_splines")) {
            rapidjson::Value& spline = EOS_json["critical_region_splines"];
            EOS.critical_region_splines.T_min = cpjson::get_double(spline, "T_min");
            EOS.critical_region_splines.T_max = cpjson::get_double(spline, "T_max");
            EOS.critical_region_splines.rhomolar_min = cpjson::get_double(spline, "rhomolar_min");
            EOS.critical_region_splines.rhomolar_max = cpjson::get_double(spline, "rhomolar_max");
            EOS.critical_region_splines.cL = cpjson::get_double_array(spline["cL"]);
            EOS.critical_region_splines.cV = cpjson::get_double_array(spline["cV"]);
            EOS.critical_region_splines.enabled = true;
        }

        // Validate the equation of state that was just created
        EOS.validate();
    }

    /// Parse the list of possible equations of state
    void parse_EOS_listing(rapidjson::Value& EOS_array, CoolPropFluid& fluid) {
        for (rapidjson::Value::ValueIterator itr = EOS_array.Begin(); itr != EOS_array.End(); ++itr) {
            parse_EOS(*itr, fluid);
        }
    };

    /// Parse the transport properties
    void parse_dilute_viscosity(rapidjson::Value& dilute, CoolPropFluid& fluid) {
        if (dilute.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("Ethane")) {
                fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_ETHANE;
                return;
            } else if (!target.compare("Cyclohexane")) {
                fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_CYCLOHEXANE;
                return;
            } else {
                throw ValueError(format("hardcoded dilute viscosity [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("collision_integral")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteGasCollisionIntegralData& CI = fluid.transport.viscosity_dilute.collision_integral;

            // Set the type flag
            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);
            CI.molar_mass = cpjson::get_double(dilute, "molar_mass");
            CI.C = cpjson::get_double(dilute, "C");
        } else if (!type.compare("kinetic_theory")) {
            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_KINETIC_THEORY;
        } else if (!type.compare("powers_of_T")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteGasPowersOfT& CI = fluid.transport.viscosity_dilute.powers_of_T;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);

            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_POWERS_OF_T;
        } else if (!type.compare("powers_of_Tr")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteGasPowersOfTr& CI = fluid.transport.viscosity_dilute.powers_of_Tr;
            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);
            CI.T_reducing = cpjson::get_double(dilute, "T_reducing");
            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_POWERS_OF_TR;
        } else if (!type.compare("collision_integral_powers_of_Tstar")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteCollisionIntegralPowersOfTstarData& CI = fluid.transport.viscosity_dilute.collision_integral_powers_of_Tstar;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);
            CI.T_reducing = cpjson::get_double(dilute, "T_reducing");
            CI.C = cpjson::get_double(dilute, "C");

            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL_POWERS_OF_TSTAR;
        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_initial_density_viscosity(rapidjson::Value& initial_density, CoolPropFluid& fluid) {
        std::string type = cpjson::get_string(initial_density, "type");
        if (!type.compare("Rainwater-Friend")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityRainWaterFriendData& RF = fluid.transport.viscosity_initial.rainwater_friend;

            // Load up the values
            RF.b = cpjson::get_long_double_array(initial_density["b"]);
            RF.t = cpjson::get_long_double_array(initial_density["t"]);

            // Set the type flag
            fluid.transport.viscosity_initial.type = CoolProp::ViscosityInitialDensityVariables::VISCOSITY_INITIAL_DENSITY_RAINWATER_FRIEND;
        } else if (!type.compare("empirical")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityInitialDensityEmpiricalData& EM = fluid.transport.viscosity_initial.empirical;

            // Load up the values
            EM.n = cpjson::get_long_double_array(initial_density["n"]);
            EM.d = cpjson::get_long_double_array(initial_density["d"]);
            EM.t = cpjson::get_long_double_array(initial_density["t"]);
            EM.T_reducing = cpjson::get_double(initial_density, "T_reducing");
            EM.rhomolar_reducing = cpjson::get_double(initial_density, "rhomolar_reducing");

            // Set the type flag
            fluid.transport.viscosity_initial.type = CoolProp::ViscosityInitialDensityVariables::VISCOSITY_INITIAL_DENSITY_EMPIRICAL;
        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_higher_order_viscosity(rapidjson::Value& higher, CoolPropFluid& fluid) {
        // First check for hardcoded higher-order term
        if (higher.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(higher, "hardcoded");
            if (!target.compare("Hydrogen")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HYDROGEN;
                return;
            } else if (!target.compare("n-Hexane")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HEXANE;
                return;
            } else if (!target.compare("n-Heptane")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HEPTANE;
                return;
            } else if (!target.compare("Toluene")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_TOLUENE;
                return;
            } else if (!target.compare("Ethane")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_ETHANE;
                return;
            } else if (!target.compare("Benzene")) {
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_BENZENE;
                return;
            } else {
                throw ValueError(
                  format("hardcoded higher order viscosity term [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }

        std::string type = cpjson::get_string(higher, "type");
        if (!type.compare("modified_Batschinski_Hildebrand")) {
            // Get a reference to the entry in the fluid instance to simplify the code that follows
            CoolProp::ViscosityModifiedBatschinskiHildebrandData& BH = fluid.transport.viscosity_higher_order.modified_Batschinski_Hildebrand;

            // Set the flag for the type of this model
            fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_BATSCHINKI_HILDEBRAND;

            BH.T_reduce = cpjson::get_double(higher, "T_reduce");
            BH.rhomolar_reduce = cpjson::get_double(higher, "rhomolar_reduce");
            // Load up the values
            BH.a = cpjson::get_long_double_array(higher["a"]);
            BH.t1 = cpjson::get_long_double_array(higher["t1"]);
            BH.d1 = cpjson::get_long_double_array(higher["d1"]);
            BH.gamma = cpjson::get_long_double_array(higher["gamma"]);
            BH.l = cpjson::get_long_double_array(higher["l"]);
            assert(BH.a.size() == BH.t1.size());
            assert(BH.a.size() == BH.d1.size());
            assert(BH.a.size() == BH.gamma.size());
            assert(BH.a.size() == BH.l.size());
            BH.f = cpjson::get_long_double_array(higher["f"]);
            BH.t2 = cpjson::get_long_double_array(higher["t2"]);
            BH.d2 = cpjson::get_long_double_array(higher["d2"]);
            assert(BH.f.size() == BH.t2.size());
            assert(BH.f.size() == BH.d2.size());
            BH.g = cpjson::get_long_double_array(higher["g"]);
            BH.h = cpjson::get_long_double_array(higher["h"]);
            assert(BH.g.size() == BH.h.size());
            BH.p = cpjson::get_long_double_array(higher["p"]);
            BH.q = cpjson::get_long_double_array(higher["q"]);
            assert(BH.p.size() == BH.q.size());
        } else if (!type.compare("friction_theory")) {
            // Get a reference to the entry in the fluid instance to simplify the code that follows
            CoolProp::ViscosityFrictionTheoryData& F = fluid.transport.viscosity_higher_order.friction_theory;

            // Set the flag for the type of this model
            fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_FRICTION_THEORY;

            // Always need these terms
            F.Ai = cpjson::get_long_double_array(higher["Ai"]);
            F.Aa = cpjson::get_long_double_array(higher["Aa"]);
            F.Aaa = cpjson::get_long_double_array(higher["Aaa"]);
            F.Ar = cpjson::get_long_double_array(higher["Ar"]);

            F.Na = cpjson::get_integer(higher, "Na");
            F.Naa = cpjson::get_integer(higher, "Naa");
            F.Nr = cpjson::get_integer(higher, "Nr");
            F.Nrr = cpjson::get_integer(higher, "Nrr");
            F.c1 = cpjson::get_double(higher, "c1");
            F.c2 = cpjson::get_double(higher, "c2");
            assert(F.Aa.size() == 3);
            assert(F.Aaa.size() == 3);
            assert(F.Ar.size() == 3);

            F.T_reduce = cpjson::get_double(higher, "T_reduce");

            if (higher.HasMember("Arr") && !higher.HasMember("Adrdr")) {
                F.Arr = cpjson::get_long_double_array(higher["Arr"]);
                assert(F.Arr.size() == 3);
            } else if (higher.HasMember("Adrdr") && !higher.HasMember("Arr")) {
                F.Adrdr = cpjson::get_long_double_array(higher["Adrdr"]);
                assert(F.Adrdr.size() == 3);
            } else {
                throw ValueError(format("can only provide one of Arr or Adrdr for fluid %s", fluid.name.c_str()));
            }
            if (higher.HasMember("Aii")) {
                F.Aii = cpjson::get_long_double_array(higher["Aii"]);
                F.Nii = cpjson::get_integer(higher, "Nii");
            }
            if (higher.HasMember("Aaaa") && higher.HasMember("Arrr")) {
                F.Aaaa = cpjson::get_long_double_array(higher["Aaaa"]);
                F.Arrr = cpjson::get_long_double_array(higher["Arrr"]);
                F.Naaa = cpjson::get_integer(higher, "Naaa");
                F.Nrrr = cpjson::get_integer(higher, "Nrrr");
            }

        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    void parse_ECS_conductivity(rapidjson::Value& conductivity, CoolPropFluid& fluid) {
        fluid.transport.conductivity_ecs.reference_fluid = cpjson::get_string(conductivity, "reference_fluid");

        // Parameters for correction polynomials
        fluid.transport.conductivity_ecs.psi_a = cpjson::get_long_double_array(conductivity["psi"]["a"]);
        fluid.transport.conductivity_ecs.psi_t = cpjson::get_long_double_array(conductivity["psi"]["t"]);
        fluid.transport.conductivity_ecs.psi_rhomolar_reducing = cpjson::get_double(conductivity["psi"], "rhomolar_reducing");
        fluid.transport.conductivity_ecs.f_int_a = cpjson::get_long_double_array(conductivity["f_int"]["a"]);
        fluid.transport.conductivity_ecs.f_int_t = cpjson::get_long_double_array(conductivity["f_int"]["t"]);
        fluid.transport.conductivity_ecs.f_int_T_reducing = cpjson::get_double(conductivity["f_int"], "T_reducing");

        fluid.transport.conductivity_using_ECS = true;
    }

    void parse_ECS_viscosity(rapidjson::Value& viscosity, CoolPropFluid& fluid) {
        fluid.transport.viscosity_ecs.reference_fluid = cpjson::get_string(viscosity, "reference_fluid");

        // Parameters for correction polynomial
        fluid.transport.viscosity_ecs.psi_a = cpjson::get_long_double_array(viscosity["psi"]["a"]);
        fluid.transport.viscosity_ecs.psi_t = cpjson::get_long_double_array(viscosity["psi"]["t"]);
        fluid.transport.viscosity_ecs.psi_rhomolar_reducing = cpjson::get_double(viscosity["psi"], "rhomolar_reducing");

        fluid.transport.viscosity_using_ECS = true;
    }

    void parse_Chung_viscosity(rapidjson::Value& viscosity, CoolPropFluid& fluid) {
        // These in base SI units
        fluid.transport.viscosity_Chung.rhomolar_critical = cpjson::get_double(viscosity, "rhomolar_critical");
        fluid.transport.viscosity_Chung.T_critical = cpjson::get_double(viscosity, "T_critical");
        fluid.transport.viscosity_Chung.molar_mass = cpjson::get_double(viscosity, "molar_mass");
        fluid.transport.viscosity_Chung.dipole_moment_D = cpjson::get_double(viscosity, "dipole_moment_D");
        fluid.transport.viscosity_Chung.acentric = cpjson::get_double(viscosity, "acentric");
        fluid.transport.viscosity_using_Chung = true;
    }

    void parse_rhosr_viscosity(rapidjson::Value& viscosity, CoolPropFluid& fluid) {
        fluid.transport.viscosity_rhosr.C = cpjson::get_double(viscosity, "C");
        fluid.transport.viscosity_rhosr.c_liq = cpjson::get_double_array(viscosity, "c_liq");
        fluid.transport.viscosity_rhosr.c_vap = cpjson::get_double_array(viscosity, "c_vap");
        fluid.transport.viscosity_rhosr.rhosr_critical = cpjson::get_double(viscosity, "rhosr_critical");
        fluid.transport.viscosity_rhosr.x_crossover = cpjson::get_double(viscosity, "x_crossover");
        fluid.transport.viscosity_using_rhosr = true;
    }

    /// Parse the transport properties
    void parse_viscosity(rapidjson::Value& viscosity, CoolPropFluid& fluid) {
        // If an array, use the first one, and then stop;
        if (viscosity.IsArray()) {
            rapidjson::Value::ValueIterator itr = viscosity.Begin();
            parse_viscosity(*itr, fluid);
            return;
        }

        // Load the BibTeX key
        fluid.transport.BibTeX_viscosity = cpjson::get_string(viscosity, "BibTeX");

        // Set the Lennard-Jones 12-6 potential variables, or approximate them from method of Chung
        if (!viscosity.HasMember("sigma_eta") || !viscosity.HasMember("epsilon_over_k")) {
            default_transport(fluid);
        } else {
            fluid.transport.sigma_eta = cpjson::get_double(viscosity, "sigma_eta");
            fluid.transport.epsilon_over_k = cpjson::get_double(viscosity, "epsilon_over_k");
        }

        // If it is using ECS, set ECS parameters and quit
        if (viscosity.HasMember("type") && !cpjson::get_string(viscosity, "type").compare("ECS")) {
            parse_ECS_viscosity(viscosity, fluid);
            return;
        }

        // If it is using rho*sr CS, set parameters and quit
        if (viscosity.HasMember("type") && !cpjson::get_string(viscosity, "type").compare("rhosr-CS")) {
            parse_rhosr_viscosity(viscosity, fluid);
            return;
        }

        // Use the method of Chung
        // If it is using ECS, set ECS parameters and quit
        if (viscosity.HasMember("type") && !cpjson::get_string(viscosity, "type").compare("Chung")) {
            parse_Chung_viscosity(viscosity, fluid);
            return;
        }

        if (viscosity.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(viscosity, "hardcoded");
            if (!target.compare("Water")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_WATER;
                return;
            } else if (!target.compare("HeavyWater")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_HEAVYWATER;
                return;
            } else if (!target.compare("Helium")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_HELIUM;
                return;
            } else if (!target.compare("R23")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_R23;
                return;
            } else if (!target.compare("Methanol")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_METHANOL;
                return;
            } else if (!target.compare("m-Xylene")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_M_XYLENE;
                return;
            } else if (!target.compare("o-Xylene")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_O_XYLENE;
                return;
            } else if (!target.compare("p-Xylene")) {
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_P_XYLENE;
                return;
            } else {
                throw ValueError(format("hardcoded viscosity [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }

        // Load dilute viscosity term
        if (viscosity.HasMember("dilute")) {
            parse_dilute_viscosity(viscosity["dilute"], fluid);
        }
        // Load initial density term
        if (viscosity.HasMember("initial_density")) {
            parse_initial_density_viscosity(viscosity["initial_density"], fluid);
        }
        // Load higher_order term
        if (viscosity.HasMember("higher_order")) {
            parse_higher_order_viscosity(viscosity["higher_order"], fluid);
        }
    };

    /// Parse the transport properties
    void parse_dilute_conductivity(rapidjson::Value& dilute, CoolPropFluid& fluid) {
        if (dilute.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("CO2")) {
                fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_CO2;
                return;
            } else if (!target.compare("Ethane")) {
                fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETHANE;
                return;
            } else if (!target.compare("none")) {
                fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_NONE;
                return;
            } else {
                throw ValueError(
                  format("hardcoded dilute conductivity term [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("ratio_of_polynomials")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityDiluteRatioPolynomialsData& data = fluid.transport.conductivity_dilute.ratio_polynomials;

            // Set the type flag
            fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_RATIO_POLYNOMIALS;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.B = cpjson::get_long_double_array(dilute["B"]);
            data.n = cpjson::get_long_double_array(dilute["n"]);
            data.m = cpjson::get_long_double_array(dilute["m"]);
            data.T_reducing = cpjson::get_double(dilute, "T_reducing");
        } else if (!type.compare("eta0_and_poly")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityDiluteEta0AndPolyData& data = fluid.transport.conductivity_dilute.eta0_and_poly;

            // Set the type flag
            fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETA0_AND_POLY;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);
        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_residual_conductivity(rapidjson::Value& dilute, CoolPropFluid& fluid) {
        if (dilute.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("CO2")) {
                fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_CO2;
                return;
            } else {
                throw ValueError(
                  format("hardcoded residual conductivity term [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("polynomial")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityResidualPolynomialData& data = fluid.transport.conductivity_residual.polynomials;

            // Set the type flag
            fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL;

            // Load up the values
            data.B = cpjson::get_long_double_array(dilute["B"]);
            data.d = cpjson::get_long_double_array(dilute["d"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);
            data.T_reducing = cpjson::get_double(dilute, "T_reducing");
            data.rhomass_reducing = cpjson::get_double(dilute, "rhomass_reducing");
        } else if (!type.compare("polynomial_and_exponential")) {
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityResidualPolynomialAndExponentialData& data = fluid.transport.conductivity_residual.polynomial_and_exponential;

            // Set the type flag
            fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL_AND_EXPONENTIAL;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.d = cpjson::get_long_double_array(dilute["d"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);
            data.gamma = cpjson::get_long_double_array(dilute["gamma"]);
            data.l = cpjson::get_long_double_array(dilute["l"]);
        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    void parse_critical_conductivity(rapidjson::Value& critical, CoolPropFluid& fluid) {
        if (critical.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(critical, "hardcoded");
            if (!target.compare("R123")) {
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_R123;
                return;
            } else if (!target.compare("Ammonia")) {
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_AMMONIA;
                return;
            } else if (!target.compare("CarbonDioxideScalabrinJPCRD2006")) {
                fluid.transport.conductivity_critical.type =
                  CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_CARBONDIOXIDE_SCALABRIN_JPCRD_2006;
                return;
            } else if (!target.compare("None")) {
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_NONE;
                return;
            } else {
                throw ValueError(format("critical conductivity term [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(critical, "type");
        if (!type.compare("simplified_Olchowy_Sengers")) {
            //// Get a reference to the entry in the fluid instance
            CoolProp::ConductivityCriticalSimplifiedOlchowySengersData& data = fluid.transport.conductivity_critical.Olchowy_Sengers;

            // Set the type flag
            fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_SIMPLIFIED_OLCHOWY_SENGERS;

            // Set values if they are found - otherwise fall back to default values
            if (critical.HasMember("qD")) {
                data.qD = cpjson::get_double(critical, "qD");
            }
            if (critical.HasMember("zeta0")) {
                data.zeta0 = cpjson::get_double(critical, "zeta0");
            }
            if (critical.HasMember("GAMMA")) {
                data.GAMMA = cpjson::get_double(critical, "GAMMA");
            }
            if (critical.HasMember("gamma")) {
                data.gamma = cpjson::get_double(critical, "gamma");
            }
            if (critical.HasMember("R0")) {
                data.R0 = cpjson::get_double(critical, "R0");
            }
            if (critical.HasMember("T_ref")) {
                data.T_ref = cpjson::get_double(critical, "T_ref");
            }
        } else {
            throw ValueError(format("type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
        }
    };

    /// Parse the thermal conductivity data
    void parse_thermal_conductivity(rapidjson::Value& conductivity, CoolPropFluid& fluid) {
        // Load the BibTeX key
        fluid.transport.BibTeX_conductivity = cpjson::get_string(conductivity, "BibTeX");

        // If it is using ECS, set ECS parameters and quit
        if (conductivity.HasMember("type") && !cpjson::get_string(conductivity, "type").compare("ECS")) {
            parse_ECS_conductivity(conductivity, fluid);
            return;
        }

        if (conductivity.HasMember("hardcoded")) {
            std::string target = cpjson::get_string(conductivity, "hardcoded");
            if (!target.compare("Water")) {
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_WATER;
                return;
            } else if (!target.compare("HeavyWater")) {
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_HEAVYWATER;
                return;
            } else if (!target.compare("Methane")) {
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_METHANE;
                return;
            } else if (!target.compare("R23")) {
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_R23;
                return;
            } else if (!target.compare("Helium")) {
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_HELIUM;
                return;
            } else {
                throw ValueError(
                  format("hardcoded residual conductivity term [%s] is not understood for fluid %s", target.c_str(), fluid.name.c_str()));
            }
        }

        // Load dilute conductivity term
        if (conductivity.HasMember("dilute")) {
            parse_dilute_conductivity(conductivity["dilute"], fluid);
        }
        // Load residual conductivity term
        if (conductivity.HasMember("residual")) {
            parse_residual_conductivity(conductivity["residual"], fluid);
        }
        // Load critical conductivity term
        if (conductivity.HasMember("critical")) {
            parse_critical_conductivity(conductivity["critical"], fluid);
        }
    };

    /// Parse the transport properties
    void parse_transport(rapidjson::Value& transport, CoolPropFluid& fluid) {

        // Parse viscosity
        if (transport.HasMember("viscosity")) {
            parse_viscosity(transport["viscosity"], fluid);
            fluid.transport.viscosity_model_provided = true;
        }

        // Parse thermal conductivity
        if (transport.HasMember("conductivity")) {
            parse_thermal_conductivity(transport["conductivity"], fluid);
            fluid.transport.conductivity_model_provided = true;
        }
    };

    void default_transport(CoolPropFluid& fluid) {
        // Use the method of Chung to approximate the values for epsilon_over_k and sigma_eta
        // Chung, T.-H.; Ajlan, M.; Lee, L. L.; Starling, K. E. Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties. Ind. Eng. Chem. Res. 1988, 27, 671-679.
        // rhoc needs to be in mol/L to yield a sigma in nm,
        CoolPropDbl rho_crit_molar = fluid.EOS().reduce.rhomolar / 1000.0;  // [mol/m3 to mol/L]
        CoolPropDbl Tc = fluid.EOS().reduce.T;
        fluid.transport.sigma_eta = 0.809 / pow(rho_crit_molar, static_cast<CoolPropDbl>(1.0 / 3.0)) / 1e9;  // 1e9 is to convert from nm to m
        fluid.transport.epsilon_over_k = Tc / 1.2593;                                                        // [K]
    }

    void parse_melting_line(rapidjson::Value& melting_line, CoolPropFluid& fluid) {
        fluid.ancillaries.melting_line.T_m = cpjson::get_double(melting_line, "T_m");
        fluid.ancillaries.melting_line.BibTeX = cpjson::get_string(melting_line, "BibTeX");

        if (melting_line.HasMember("type")) {
            std::string type = cpjson::get_string(melting_line, "type");
            if (!type.compare("Simon")) {
                rapidjson::Value& parts = melting_line["parts"];
                fluid.ancillaries.melting_line.type = MeltingLineVariables::MELTING_LINE_SIMON_TYPE;
                for (rapidjson::Value::ValueIterator itr = parts.Begin(); itr != parts.End(); ++itr) {
                    MeltingLinePiecewiseSimonSegment data;
                    data.a = cpjson::get_double((*itr), "a");
                    data.c = cpjson::get_double((*itr), "c");
                    data.T_min = cpjson::get_double((*itr), "T_min");
                    data.T_max = cpjson::get_double((*itr), "T_max");
                    data.T_0 = cpjson::get_double((*itr), "T_0");
                    data.p_0 = cpjson::get_double((*itr), "p_0");
                    fluid.ancillaries.melting_line.simon.parts.push_back(data);
                }
            } else if (!type.compare("polynomial_in_Tr")) {
                rapidjson::Value& parts = melting_line["parts"];
                fluid.ancillaries.melting_line.type = MeltingLineVariables::MELTING_LINE_POLYNOMIAL_IN_TR_TYPE;
                for (rapidjson::Value::ValueIterator itr = parts.Begin(); itr != parts.End(); ++itr) {
                    MeltingLinePiecewisePolynomialInTrSegment data;
                    data.a = cpjson::get_long_double_array((*itr), "a");
                    data.t = cpjson::get_long_double_array((*itr), "t");
                    data.T_min = cpjson::get_double((*itr), "T_min");
                    data.T_max = cpjson::get_double((*itr), "T_max");
                    data.T_0 = cpjson::get_double((*itr), "T_0");
                    data.p_0 = cpjson::get_double((*itr), "p_0");
                    fluid.ancillaries.melting_line.polynomial_in_Tr.parts.push_back(data);
                }
            } else if (!type.compare("polynomial_in_Theta")) {
                rapidjson::Value& parts = melting_line["parts"];
                fluid.ancillaries.melting_line.type = MeltingLineVariables::MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE;
                for (rapidjson::Value::ValueIterator itr = parts.Begin(); itr != parts.End(); ++itr) {
                    MeltingLinePiecewisePolynomialInThetaSegment data;
                    data.a = cpjson::get_long_double_array((*itr), "a");
                    data.t = cpjson::get_long_double_array((*itr), "t");
                    data.T_min = cpjson::get_double((*itr), "T_min");
                    data.T_max = cpjson::get_double((*itr), "T_max");
                    data.T_0 = cpjson::get_double((*itr), "T_0");
                    data.p_0 = cpjson::get_double((*itr), "p_0");
                    fluid.ancillaries.melting_line.polynomial_in_Theta.parts.push_back(data);
                }
            } else {
                throw ValueError(format("melting line type [%s] is not understood for fluid %s", type.c_str(), fluid.name.c_str()));
            }
            // Set the limits for the melting line curve
            fluid.ancillaries.melting_line.set_limits();
        } else {
            throw ValueError(format("melting line does not have \"type\" for fluid %s", fluid.name.c_str()));
        }
    };

    /// Parse the critical state for the given EOS
    void parse_states(rapidjson::Value& states, CoolPropFluid& fluid) {
        if (!states.HasMember("critical")) {
            throw ValueError(format("fluid[\"STATES\"] [%s] does not have \"critical\" member", fluid.name.c_str()));
        }
        rapidjson::Value& crit = states["critical"];
        fluid.crit.T = cpjson::get_double(crit, "T");
        fluid.crit.p = cpjson::get_double(crit, "p");
        fluid.crit.rhomolar = cpjson::get_double(crit, "rhomolar");
        fluid.crit.hmolar = cpjson::get_double(crit, "hmolar");
        fluid.crit.smolar = cpjson::get_double(crit, "smolar");

        if (!states.HasMember("triple_liquid")) {
            throw ValueError(format("fluid[\"STATES\"] [%s] does not have \"triple_liquid\" member", fluid.name.c_str()));
        }
        rapidjson::Value& triple_liquid = states["triple_liquid"];
        if (triple_liquid.ObjectEmpty()) {
            // State is empty - probably because the triple point temperature is below the minimum saturation temperature
            fluid.triple_liquid.T = -1;
            fluid.triple_liquid.p = -1;
            fluid.triple_liquid.rhomolar = -1;
            fluid.triple_liquid.hmolar = _HUGE;
            fluid.triple_liquid.smolar = _HUGE;
        } else {
            fluid.triple_liquid.T = cpjson::get_double(triple_liquid, "T");
            fluid.triple_liquid.p = cpjson::get_double(triple_liquid, "p");
            fluid.triple_liquid.rhomolar = cpjson::get_double(triple_liquid, "rhomolar");
            fluid.triple_liquid.hmolar = cpjson::get_double(triple_liquid, "hmolar");
            fluid.triple_liquid.smolar = cpjson::get_double(triple_liquid, "smolar");
        }

        if (!states.HasMember("triple_vapor")) {
            throw ValueError(format("fluid[\"STATES\"] [%s] does not have \"triple_vapor\" member", fluid.name.c_str()));
        }
        rapidjson::Value& triple_vapor = states["triple_vapor"];
        if (triple_vapor.ObjectEmpty()) {
            // State is empty - probably because the triple point temperature is below the minimum saturation temperature
            fluid.triple_vapor.T = -1;
            fluid.triple_vapor.p = -1;
            fluid.triple_vapor.rhomolar = -1;
            fluid.triple_vapor.hmolar = _HUGE;
            fluid.triple_vapor.smolar = _HUGE;
        } else {
            fluid.triple_vapor.T = cpjson::get_double(triple_vapor, "T");
            fluid.triple_vapor.p = cpjson::get_double(triple_vapor, "p");
            fluid.triple_vapor.rhomolar = cpjson::get_double(triple_vapor, "rhomolar");
            fluid.triple_vapor.hmolar = cpjson::get_double(triple_vapor, "hmolar");
            fluid.triple_vapor.smolar = cpjson::get_double(triple_vapor, "smolar");
        }
    };

    /// Parse the critical state for the given EOS
    void parse_ancillaries(rapidjson::Value& ancillaries, CoolPropFluid& fluid) {
        if (!ancillaries.HasMember("rhoL") || !ancillaries.HasMember("rhoV")) {
            throw ValueError("Ancillary curves for either rhoL or rhoV are missing");
        }
        fluid.ancillaries.rhoL = SaturationAncillaryFunction(ancillaries["rhoL"]);
        fluid.ancillaries.rhoV = SaturationAncillaryFunction(ancillaries["rhoV"]);

        // If a pseudo-pure fluid, has pL and pV curves
        if (ancillaries.HasMember("pL") && ancillaries.HasMember("pV")) {
            fluid.ancillaries.pL = SaturationAncillaryFunction(ancillaries["pL"]);
            fluid.ancillaries.pV = SaturationAncillaryFunction(ancillaries["pV"]);
        }
        // Otherwise has a single pS curve and not pL and not pV
        else if (!ancillaries.HasMember("pL") && !ancillaries.HasMember("pV") && ancillaries.HasMember("pS")) {
            fluid.ancillaries.pL = SaturationAncillaryFunction(ancillaries["pS"]);
            fluid.ancillaries.pV = SaturationAncillaryFunction(ancillaries["pS"]);
        } else {
            throw ValueError("Pressure ancillary curves are missing or invalid");
        }

        if (ancillaries.HasMember("hL")) {
            fluid.ancillaries.hL = SaturationAncillaryFunction(ancillaries["hL"]);
        } else {
            if (get_debug_level() > 0) {
                std::cout << "Missing hL ancillary for fluid " << fluid.name;
            }
        }
        if (ancillaries.HasMember("hLV")) {
            fluid.ancillaries.hLV = SaturationAncillaryFunction(ancillaries["hLV"]);
        } else {
            if (get_debug_level() > 0) {
                std::cout << "Missing hLV ancillary for fluid " << fluid.name;
            }
        }

        if (ancillaries.HasMember("sL")) {
            fluid.ancillaries.sL = SaturationAncillaryFunction(ancillaries["sL"]);
        } else {
            if (get_debug_level() > 0) {
                std::cout << "Missing sL ancillary for fluid " << fluid.name;
            }
        }
        if (ancillaries.HasMember("sLV")) {
            fluid.ancillaries.sLV = SaturationAncillaryFunction(ancillaries["sLV"]);
        } else {
            if (get_debug_level() > 0) {
                std::cout << "Missing sLV ancillary for fluid " << fluid.name;
            }
        }
        if (!ValidNumber(fluid.ancillaries.sL.get_Tmin()) && get_debug_level() > 0) {
            std::cout << "Tmin invalid for sL for " << fluid.name << std::endl;
        }
    };

    /// Parse the surface_tension
    void parse_surface_tension(rapidjson::Value& surface_tension, CoolPropFluid& fluid) {
        fluid.ancillaries.surface_tension = SurfaceTensionCorrelation(surface_tension);
    };

    /// Validate the fluid file that was just constructed
    void validate(CoolPropFluid& fluid) {
        assert(fluid.EOSVector.size() > 0);
        assert(fluid.CAS.length() > 0);
        assert(fluid.name.length() > 0);
    }

   public:
    // Default constructor;
    JSONFluidLibrary() {
        _is_empty = true;
    };
    bool is_empty(void) {
        return _is_empty;
    };

    /// Add all the fluid entries in the JSON-encoded string passed in
    static void add_many(const std::string& JSON_string);

    /// Add all the fluid entries in the rapidjson::Value instance passed in
    void add_many(rapidjson::Value& listing);

    void add_one(rapidjson::Value& fluid_json);

    std::string get_JSONstring(const std::string& key) {
        // Try to find it
        std::map<std::string, std::size_t>::const_iterator it = string_to_index_map.find(key);
        if (it != string_to_index_map.end()) {

            std::map<std::size_t, std::string>::const_iterator it2 = JSONstring_map.find(it->second);
            if (it2 != JSONstring_map.end()) {
                // Then, load the fluids we would like to add
                rapidjson::Document doc;
                cpjson::JSON_string_to_rapidjson(it2->second, doc);
                rapidjson::Document doc2;
                doc2.SetArray();
                doc2.PushBack(doc, doc.GetAllocator());
                return cpjson::json2string(doc2);
            } else {
                throw ValueError(format("Unable to obtain JSON string for this identifier [%d]", it->second));
            }
        } else {
            throw ValueError(format("Unable to obtain index for this identifier [%s]", key.c_str()));
        }
    }

    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key Either a CAS number or the name (CAS number should be preferred)
    */
    CoolPropFluid get(const std::string& key) {
        // Try to find it
        std::map<std::string, std::size_t>::const_iterator it = string_to_index_map.find(key);
        // If it is found
        if (it != string_to_index_map.end()) {
            return get(it->second);
        } else {
            // Here we check for the use of a cubic Helmholtz energy transformation for a multi-fluid model
            std::vector<std::string> endings;
            endings.push_back("-SRK");
            endings.push_back("-PengRobinson");
            for (std::vector<std::string>::const_iterator end = endings.begin(); end != endings.end(); ++end) {
                if (endswith(key, *end)) {
                    std::string used_name = key.substr(0, key.size() - (*end).size());
                    it = string_to_index_map.find(used_name);
                    if (it != string_to_index_map.end()) {
                        // We found the name of the fluid within the library of multiparameter
                        // Helmholtz-explicit models.  We will load its parameters from the
                        // multiparameter EOS
                        //
                        CoolPropFluid fluid = get(it->second);
                        // Remove all the residual contributions to the Helmholtz energy
                        fluid.EOSVector[0].alphar.empty_the_EOS();
                        // Get the parameters for the cubic EOS
                        CoolPropDbl Tc = fluid.EOSVector[0].reduce.T;
                        CoolPropDbl pc = fluid.EOSVector[0].reduce.p;
                        CoolPropDbl rhomolarc = fluid.EOSVector[0].reduce.rhomolar;
                        CoolPropDbl acentric = fluid.EOSVector[0].acentric;
                        CoolPropDbl R = 8.3144598;  // fluid.EOSVector[0].R_u;
                        // Set the cubic contribution to the residual Helmholtz energy
                        shared_ptr<AbstractCubic> ac;
                        if (*end == "-SRK") {
                            ac.reset(new SRK(Tc, pc, acentric, R));
                        } else if (*end == "-PengRobinson") {
                            ac.reset(new PengRobinson(Tc, pc, acentric, R));
                        } else {
                            throw CoolProp::ValueError(format("Unable to match this ending [%s]", (*end).c_str()));
                        }
                        ac->set_Tr(Tc);
                        ac->set_rhor(rhomolarc);
                        fluid.EOSVector[0].alphar.cubic = ResidualHelmholtzGeneralizedCubic(ac);
                        return fluid;
                    } else {
                        // Let's look in the library of cubic EOS
                        CubicLibrary::CubicsValues vals = CubicLibrary::get_cubic_values(used_name);
                        // Set the cubic contribution to the residual Helmholtz energy
                        shared_ptr<AbstractCubic> ac;
                        if (*end == "-SRK") {
                            ac.reset(new SRK(vals.Tc, vals.pc, vals.acentric, get_config_double(R_U_CODATA)));
                        } else if (*end == "-PengRobinson") {
                            ac.reset(new PengRobinson(vals.Tc, vals.pc, vals.acentric, get_config_double(R_U_CODATA)));
                        } else {
                            throw CoolProp::ValueError(format("Unable to match this ending [%s]", (*end).c_str()));
                        }
                        ac->set_Tr(vals.Tc);
                        if (vals.rhomolarc > 0) {
                            ac->set_rhor(vals.rhomolarc);
                        } else {
                            // Curve fit from all the pure fluids in CoolProp (thanks to recommendation of A. Kazakov)
                            double v_c_Lmol = 2.14107171795 * (vals.Tc / vals.pc * 1000) + 0.00773144012514;  // [L/mol]
                            ac->set_rhor(1 / (v_c_Lmol / 1000.0));
                        }
                        if (vals.alpha_type == "Twu") {
                            std::vector<double>& c = vals.alpha_coeffs;
                            ac->set_C_Twu(0, c[0], c[1], c[2]);
                        }
                        CoolPropFluid fluid;
                        fluid.CAS = vals.CAS;
                        EquationOfState E;
                        E.acentric = vals.acentric;
                        E.sat_min_liquid.T = _HUGE;
                        E.sat_min_liquid.p = _HUGE;
                        E.reduce.T = vals.Tc;
                        E.reduce.p = vals.pc;
                        E.reduce.rhomolar = ac->get_rhor();
                        fluid.EOSVector.push_back(E);
                        fluid.EOS().alphar.cubic = ResidualHelmholtzGeneralizedCubic(ac);
                        fluid.EOS().alpha0 = vals.alpha0;
                        fluid.crit.T = vals.Tc;
                        fluid.crit.p = vals.pc;
                        fluid.crit.rhomolar = ac->get_rhor();

                        return fluid;
                    }
                }
            }
            throw ValueError(format("key [%s] was not found in string_to_index_map in JSONFluidLibrary", key.c_str()));
        }
    };

    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key The index of the fluid in the map
    */
    CoolPropFluid get(std::size_t key) {
        // Try to find it
        std::map<std::size_t, CoolPropFluid>::iterator it = fluid_map.find(key);
        // If it is found
        if (it != fluid_map.end()) {
            return it->second;
        } else {
            throw ValueError(format("key [%d] was not found in JSONFluidLibrary", key));
        }
    };
    void set_fluid_enthalpy_entropy_offset(const std::string& fluid, double delta_a1, double delta_a2, const std::string& ref);
    /// Return a comma-separated list of fluid names
    std::string get_fluid_list(void) {
        return strjoin(name_vector, get_config_string(LIST_STRING_DELIMITER));
    };
};

/// Get a reference to the library instance
JSONFluidLibrary& get_library(void);

/// Get a comma-separated-list of fluids that are included
std::string get_fluid_list(void);

/// Get the fluid structure
CoolPropFluid get_fluid(const std::string& fluid_string);

/// Get the fluid as a JSON string, suitable for modification and reloading
std::string get_fluid_as_JSONstring(const std::string& indentifier);

/// Set the internal enthalpy and entropy offset variables
void set_fluid_enthalpy_entropy_offset(const std::string& fluid, double delta_a1, double delta_a2, const std::string& ref);

} /* namespace CoolProp */
#endif
