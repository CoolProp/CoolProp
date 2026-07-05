#include "IncompressibleLibrary.h"
#include "CoolProp/numerics/MatrixMath.h"
#include "CoolProp/DataStructures.h"
//#include "crossplatform_shared_ptr.h"
#include "CoolProp/detail/json.h"
#include "all_incompressibles_JSON.h"  // Makes a std::string variable called all_incompressibles_JSON

#include <mutex>

namespace CoolProp {

///// Class to access Lithium-Bromide solutions
///** Employs some basic wrapper-like functionality
// *  to bridge the gap between the solution functions
// *  used in the paper by Pátek and Klomfar:
// *  http://dx.doi.org/10.1016/j.ijrefrig.2005.10.007
// *
// *  We owe gratitude to the authors for providing
// *  both access to the paper as well as the equations
// *  in the form of C source code. */
//
//double const LiBrSolution::M_H2O  = 0.018015268; /* kg/mol, molar mass of H2O */
//double const LiBrSolution::M_LiBr = 0.08685; /* kg/mol, molar mass of LiBr */
//double const LiBrSolution::T0     = 221; /* K, constant */
//
///* Critical point of H2O */
//double const LiBrSolution::Tc_H2O   = 647.096; /* K, temperature  */
//double const LiBrSolution::pc_H2O   = 22.064; /* MPa, pressure */
//double const LiBrSolution::rhoc_H2O = 17873; /* mol/m^3 (322 kg/m^3), molar density */
//double const LiBrSolution::hc_H2O   = 37548.5; /* J/mol, molar enthalpy */
//double const LiBrSolution::sc_H2O   = 79.3933; /* J/(mol.K) molar entropy*/
//
///*Triple point of H2O */
//double const LiBrSolution::Tt_H2O   = 273.16; /* K, temperature */
//double const LiBrSolution::cpt_H2O  = 76.0226; /* J/(mol.K), molar isobaric heat capacity*/
//
//double LiBrSolution::ps_mix(double T, double x)
///* Equation (1) */
//{
//    static double m[8] = { 3.0, 4.0, 4.0, 8.0, 1.0, 1.0, 4.0, 6.0 };
//    static double n[8] = { 0.0, 5.0, 6.0, 3.0, 0.0, 2.0, 6.0, 0.0 };
//    static double t[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
//    static double a[8] = { -2.41303e2, 1.91750e7, -1.75521e8, 3.25430e7,
//            3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3 };
//    double tau, suma = 0.0;
//    int i;
//
//    tau = T / Tc_H2O;
//    for (i = 0; i <= 7; i++)
//        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
//    return (ps_H2O(T - suma));
//
//} /* end function ps_mix */
//
//double LiBrSolution::rho_mix(double T, double x)
///* Equation (2) */
//{
//    static double m[2] = { 1.0, 1.0 };
//    static double n[2] = { 0.0, 6.0 };
//    static double a[2] = { 1.746, 4.709 };
//
//    double tau, suma = 0.0;
//    int i;
//
//    tau = T / Tc_H2O;
//    for (i = 0; i <= 1; i++)
//        suma += a[i] * pow(x, m[i]) * pow(tau, n[i]);
//
//    return ((1.0 - x) * rho_H2O(T) + rhoc_H2O * suma);
//
//} /* end function rho_mix */
//
//double LiBrSolution::cp_mix(double T, double x)
///* Equation (3) */
//{
//    static double m[8] = { 2.0, 3.0, 3.0, 3.0, 3.0, 2.0, 1.0, 1.0 };
//    static double n[8] = { 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 3.0, 2.0 };
//    static double t[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0 };
//    static double a[8] = { -1.42094e1, 4.04943e1, 1.11135e2, 2.29980e2,
//            1.34526e3, -1.41010e-2, 1.24977e-2, -6.83209e-4 };
//
//    double tau, suma = 0.0;
//    int i;
//
//    tau = Tc_H2O / (T - T0);
//    for (i = 0; i <= 7; i++)
//        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
//
//    return ((1.0 - x) * cp_H2O(T) + cpt_H2O * suma);
//
//} /* end function cp_mix */
//
//double LiBrSolution::h_mix(double T, double x)
///* Equation (4) */
//{
//    static double m[30] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 4.0,
//            5.0, 5.0, 6.0, 6.0, 1.0, 2.0, 2.0, 2.0, 5.0, 6.0, 7.0, 1.0, 1.0,
//            2.0, 2.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };
//
//    static double n[30] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
//            4.0, 5.0, 5.0, 6.0, 0.0, 3.0, 5.0, 7.0, 0.0, 3.0, 1.0, 0.0, 4.0,
//            2.0, 6.0, 7.0, 0.0, 0.0, 1.0, 2.0, 3.0 };
//
//    static double t[30] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
//            2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
//            4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };
//
//    static double a[30] = { 2.27431, -7.99511, 3.85239e2, -1.63940e4,
//            -4.22562e2, 1.13314e-1, -8.33474, -1.73833e4, 6.49763,
//            3.24552e3, -1.34643e4, 3.99322e4, -2.58877e5, -1.93046e-3,
//            2.80616, -4.04479e1, 1.45342e2, -2.74873, -4.49743e2,
//            -1.21794e1, -5.83739e-3, 2.33910e-1, 3.41888e-1, 8.85259,
//            -1.78731e1, 7.35179e-2, -1.79430e-4, 1.84261e-3, -6.24282e-3,
//            6.84765e-3 };
//
//    double tau, suma = 0.0;
//    int i;
//
//    tau = Tc_H2O / (T - T0);
//    for (i = 0; i <= 29; i++)
//        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
//
//    return ((1.0 - x) * h_H2O(T) + hc_H2O * suma);
//
//} /* end function h_mix */
//
//double LiBrSolution::s_mix(double T, double x)
///* Equation (5) */
//{
//    static double m[29] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 1.0,
//            2.0, 2.0, 4.0, 5.0, 5.0, 6.0, 6.0, 1.0, 3.0, 5.0, 7.0, 1.0, 1.0,
//            1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };
//
//    static double n[29] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
//            0.0, 4.0, 0.0, 4.0, 5.0, 2.0, 5.0, 0.0, 4.0, 0.0, 1.0, 0.0, 2.0,
//            4.0, 7.0, 1.0, 0.0, 1.0, 2.0, 3.0 };
//
//    static double t[29] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
//            2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
//            4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };
//
//    static double a[29] = { 1.53091, -4.52564, 6.98302e+2, -2.16664e+4,
//            -1.47533e+3, 8.47012e-2, -6.59523, -2.95331e+4, 9.56314e-3,
//            -1.88679e-1, 9.31752, 5.78104, 1.38931e+4, -1.71762e+4,
//            4.15108e+2, -5.55647e+4, -4.23409e-3, 3.05242e+1, -1.67620,
//            1.48283e+1, 3.03055e-3, -4.01810e-2, 1.49252e-1, 2.59240,
//            -1.77421e-1, -6.99650e-5, 6.05007e-4, -1.65228e-3, 1.22966e-3 };
//
//    double tau, suma = 0.0;
//    int i;
//
//    tau = Tc_H2O / (T - T0);
//    for (i = 0; i <= 28; i++)
//        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
//
//    return ((1.0 - x) * s_H2O(T) + sc_H2O * suma);
//
//} /* end function s_mix */
//
//double LiBrSolution::ps_H2O(double T)
///* Equation (28) */
//{
//    static double a[7] = { 0.0, -7.85951783, 1.84408259, -11.7866497,
//            22.6807411, -15.9618719, 1.80122502 };
//
//    double tau, ps;
//
//    tau = 1 - T / Tc_H2O;
//
//    ps = pc_H2O
//            * exp(
//                    Tc_H2O / T
//                            * (a[1] * tau + a[2] * pow(tau, 1.5)
//                                    + a[3] * pow(tau, 3.0)
//                                    + a[4] * pow(tau, 3.5)
//                                    + a[5] * pow(tau, 4.0)
//                                    + a[6] * pow(tau, 7.5)));
//
//    return (ps * 1.0e6);
//
//} /* end function ps_H2O */
//
//double LiBrSolution::rho_H2O(double T)
///* Equation (29) */
//{
//    static double b[7] = { 0.0, 1.99274064, 1.09965342, -0.510839303,
//            -1.75493479, -45.5170352, -6.7469445e5 };
//    double theta, rho;
//
//    theta = 1.0 - T / Tc_H2O;
//
//    rho = rhoc_H2O
//            * (1.0 + b[1] * pow(theta, 1.0 / 3.0)
//                    + b[2] * pow(theta, 2.0 / 3.0)
//                    + b[3] * pow(theta, 5.0 / 3.0)
//                    + b[4] * pow(theta, 16.0 / 3.0)
//                    + b[5] * pow(theta, 43.0 / 3.0)
//                    + b[6] * pow(theta, 110.0 / 3.0));
//
//    return (rho);
//} /* end function rho_H2O */
//
//double LiBrSolution::cp_H2O(double T)
///* Equation (30) */
//{
//    static double a[5] =
//            { 1.38801, -2.95318, 3.18721, -0.645473, 9.18946e5 };
//    static double b[5] = { 0.0, 2.0, 3.0, 6.0, 34.0 };
//    static double c[5] = { 0.0, 2.0, 3.0, 5.0, 0.0 };
//
//    double suma = 0.0;
//    int i;
//
//    for (i = 0; i <= 4; i++)
//        suma += a[i] * exp(b[i] * log(1.0 - T / Tc_H2O))
//                * exp(c[i] * log(T / Tt_H2O));
//
//    return (cpt_H2O * suma);
//
//} /* end function cp_H2O */
//
//double LiBrSolution::h_H2O(double T)
///* Equation (31) */
//{
//    static double a[4] = { -4.37196e-1, 3.03440e-1, -1.29582, -1.76410e-1 };
//    static double alpha[4] = { 1.0 / 3.0, 2.0 / 3.0, 5.0 / 6.0, 21.0 / 6.0 };
//
//    double suma = 0.0;
//    int i;
//
//    for (i = 0; i <= 3; i++)
//        suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));
//
//    return (hc_H2O * (1.0 + suma));
//
//} /* end function h_H2O */
//
//double LiBrSolution::s_H2O(double T)
///* Equation (32)  */
//{
//    static double a[4] = { -3.34112e-1, -8.47987e-1, -9.11980e-1, -1.64046 };
//    static double alpha[4] = { 1.0 / 3.0, 3.0 / 3.0, 8.0 / 3.0, 24.0 / 3.0 };
//
//    double suma = 0.0;
//    int i;
//
//    for (i = 0; i <= 3; i++)
//        suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));
//
//    return (sc_H2O * (1.0 + suma));
//
//} /* end function s_H2O */
//
//
///** Finished with the code from the paper. Now we need to
// *  convert the molar values to mass-based units. */
//double LiBrSolution::massToMole(double w)
///* Equation (7)  */
//{
//    return (w/M_LiBr)/(w/M_LiBr+(1.-w)/M_H2O);
//    //return (w*M_LiBr)/(w*M_LiBr+(1.-w)*M_H2O);
//}
//
//double LiBrSolution::molarToSpecific(double w, double value)
///* Equation (7,8)  */
//{
//    double x = massToMole(w);
//    //return w/(x*M_LiBr) * value;
//    return 1. / ( x*M_LiBr + (1.-x)*M_H2O ) * value;
//}
//
//bool const LiBrSolution::debug = false;
//
//
//
//LiBrSolution::LiBrSolution():IncompressibleFluid(){
//    name = std::string("LiBr");
//    description = std::string("Lithium-Bromide solution from Patek2006");
//    reference = std::string("Patek2006");
//
//    Tmin     = 273.00;
//    Tmax     = 500.00;
//    TminPsat = Tmin;
//
//    xmin     = 0.0;
//    xmax     = 1.0;
//
//    xbase  = 0.0;
//    Tbase  = 0.0;
//
//};
//
//double LiBrSolution::rho(double T, double p, double x){
//    checkTPX(T, p, x);
//    return 1./molarToSpecific(x, 1./rho_mix(T,massToMole(x)));
//}
//double LiBrSolution::c(double T, double p, double x){
//    checkTPX(T, p, x);
//    return molarToSpecific(x, cp_mix(T,massToMole(x)));
//}
////double h(double T, double p, double x){
////    return h_u(T,p,x);
////}
//double LiBrSolution::s(double T, double p, double x){
//    checkTPX(T, p, x);
//    return molarToSpecific(x, s_mix(T,massToMole(x)));
//}
//double LiBrSolution::visc(double T, double p, double x){
//    throw ValueError("Viscosity is not defined for LiBr-solutions.");
//}
//double LiBrSolution::cond(double T, double p, double x){
//    throw ValueError("Thermal conductivity is not defined for LiBr-solutions.");
//}
//double LiBrSolution::u(double T, double p, double x){
//    checkTPX(T, p, x);
//    return molarToSpecific(x, h_mix(T,massToMole(x)));
//}
//double LiBrSolution::psat(double T, double x){
//    //checkT(T,p,x);
//    if (debug) throw ValueError(format("Your concentration is %f in kg/kg and %f in mol/mol.",x,massToMole(x)));
//    return ps_mix(T,massToMole(x));
//};
//double LiBrSolution::Tfreeze(double p, double x){
//    if (debug) throw ValueError(format("No freezing point data available for Lithium-Bromide: p=%f, x=%f",p,x));
//    return Tmin;
//}

/// Default constructor
JSONIncompressibleLibrary::JSONIncompressibleLibrary()
  : _is_empty(true) {

        //    fluid_map.clear();
        //    name_vector.clear();
        //    string_to_index_map.clear();
        //
        //    //shared_ptr<double> array (new double [256], ArrayDeleter<double> ());
    };

/// Default destructor
JSONIncompressibleLibrary::~JSONIncompressibleLibrary() {
    //    freeClear(fluid_map);
    //      fluid_map.clear();
    //    name_vector.clear();
    //    string_to_index_map.clear();
};

/// A general function to parse the json files that hold the coefficient matrices
IncompressibleData JSONIncompressibleLibrary::parse_coefficients(const nlohmann::json& obj, const std::string& id, bool vital) {
    IncompressibleData fluidData;
    if (obj.contains(id)) {
        const nlohmann::json& entry = obj.at(id);
        if (entry.contains("type")) {
            if (entry.contains("coeffs")) {
                std::string type = cpjson::get_string(entry, "type");
                if (!type.compare("polynomial")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("exponential")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("logexponential")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("exppolynomial")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("polyoffset")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYOFFSET;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("chebyshev")) {
                    // Chebyshev-in-T x monomial-in-(x - xbase) caloric fit; the
                    // fit domain is explicit per entry (it is the data range,
                    // which need not equal the fluid-level Tmin/Tmax exactly).
                    // The derived integral/derivative matrices are built in
                    // IncompressibleFluid::validate() once the fluid is
                    // assembled.
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(entry.at("coeffs")));
                    std::vector<double> Trange = cpjson::get_double_array(entry, "Trange");
                    if (Trange.size() != 2) {
                        throw ValueError(
                          format("The \"Trange\" of [%s] must have exactly 2 entries, got %d.", id.c_str(), static_cast<int>(Trange.size())));
                    }
                    fluidData.cheb_Tmin = Trange[0];
                    fluidData.cheb_Tmax = Trange[1];
                    fluidData.cheb_xbase = entry.contains("xbase") ? cpjson::get_double(entry, "xbase") : 0.0;
                    return fluidData;
                } else if (vital) {
                    throw ValueError(format("The type [%s] is not understood for [%s] of incompressible fluids. Please check your JSON file.",
                                            type.c_str(), id.c_str()));
                }
            } else {
                throw ValueError(format("Your file does not have an entry for \"coeffs\" in [%s], which is vital for this function.", id.c_str()));
            }
        } else {
            throw ValueError(format("Your file does not have an entry for \"type\" in [%s], which is vital for this function.", id.c_str()));
        }
    } else {
        if (vital) {
            throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
        }
    }
    return fluidData;
}

/// Get a double from the JSON storage if it is defined, otherwise return def
double JSONIncompressibleLibrary::parse_value(const nlohmann::json& obj, const std::string& id, bool vital, double def = 0.0) {
    if (obj.contains(id)) {
        return cpjson::get_double(obj, id);
    } else {
        if (vital) {
            throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
        } else {
            return def;
        }
    }
}

/// Get an integer from the JSON storage to identify the composition
composition_types JSONIncompressibleLibrary::parse_ifrac(const nlohmann::json& obj, const std::string& id) {
    std::string res = cpjson::get_string(obj, id);
    if (!res.compare("mass")) return IFRAC_MASS;
    if (!res.compare("mole")) return IFRAC_MOLE;
    if (!res.compare("volume")) return IFRAC_VOLUME;
    if (!res.compare("not defined")) return IFRAC_UNDEFINED;
    if (!res.compare("pure")) return IFRAC_PURE;

    throw ValueError(format("Cannot recognise the entry for [%s], [%s] is unknown for incompressible fluids.", id.c_str(), res.c_str()));
    return IFRAC_UNDEFINED;
}

/// Add all the fluid entries in the nlohmann::json array passed in
void JSONIncompressibleLibrary::add_many(const nlohmann::json& listing) {
    for (const auto& fluid_json : listing) {
        add_one(fluid_json);
    }
};

void JSONIncompressibleLibrary::add_one(const nlohmann::json& fluid_json) {
    _is_empty = false;

    // Build the fluid locally first: nothing is registered until parsing and
    // validation succeed, so a malformed definition (now reachable at runtime
    // through add_fluids_as_JSON) cannot leave a half-initialized entry in
    // the maps.
    IncompressibleFluid fluid;
    fluid.setName("unloaded");
    try {
        fluid.setName(cpjson::get_string(fluid_json, "name"));
        if (get_debug_level() >= 20) std::cout << format("Incompressible library: Loading base values for %s ", fluid.getName().c_str()) << '\n';
        fluid.setDescription(cpjson::get_string(fluid_json, "description"));
        fluid.setReference(cpjson::get_string(fluid_json, "reference"));
        fluid.setTmax(parse_value(fluid_json, "Tmax", true, 0.0));
        fluid.setTmin(parse_value(fluid_json, "Tmin", true, 0.0));
        fluid.setxmax(parse_value(fluid_json, "xmax", false, 1.0));
        fluid.setxmin(parse_value(fluid_json, "xmin", false, 0.0));
        fluid.setxid(parse_ifrac(fluid_json, "xid"));
        fluid.setTminPsat(parse_value(fluid_json, "TminPsat", false, 0.0));

        fluid.setTbase(parse_value(fluid_json, "Tbase", false, 0.0));
        fluid.setxbase(parse_value(fluid_json, "xbase", false, 0.0));

        /// Setters for the coefficients
        if (get_debug_level() >= 20) std::cout << format("Incompressible library: Loading coefficients for %s ", fluid.getName().c_str()) << '\n';
        // The caloric properties prefer the Chebyshev entries when the JSON
        // carries them (exact singularity-free enthalpy/entropy integrals,
        // see dev/incompressible_liquids/NOTES_thermodynamic_consistency.md);
        // the classic polynomial entries remain the fallback and the format
        // for every other property. Flip this to false to A/B against the
        // polynomial caloric path with the same library.
        static constexpr bool prefer_chebyshev_caloric = true;
        auto parse_caloric = [this, &fluid_json](const std::string& id) {
            if (prefer_chebyshev_caloric && fluid_json.contains(id + "_cheb")) {
                IncompressibleData data = parse_coefficients(fluid_json, id + "_cheb", false);
                if (data.type != CoolProp::IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV) {
                    // The entry exists but did not parse as chebyshev (e.g. a
                    // typo in "type"): failing loudly beats silently running
                    // on the polynomial fallback while the author believes
                    // the Chebyshev fit is in use.
                    throw ValueError(format("The entry [%s_cheb] exists but its type is not \"chebyshev\"; fix or remove it.", id.c_str()));
                }
                return data;
            }
            return parse_coefficients(fluid_json, id, true);
        };
        fluid.setDensity(parse_caloric("density"));
        fluid.setSpecificHeat(parse_caloric("specific_heat"));
        fluid.setViscosity(parse_coefficients(fluid_json, "viscosity", false));
        fluid.setConductivity(parse_coefficients(fluid_json, "conductivity", false));
        fluid.setPsat(parse_coefficients(fluid_json, "saturation_pressure", false));
        fluid.setTfreeze(parse_coefficients(fluid_json, "T_freeze", false));
        fluid.setMass2input(parse_coefficients(fluid_json, "mass2input", false));
        fluid.setVolume2input(parse_coefficients(fluid_json, "volume2input", false));
        fluid.setMole2input(parse_coefficients(fluid_json, "mole2input", false));

        //if (get_debug_level()>=20) std::cout << format("Incompressible library: Loading reference state for %s ",fluid.getName().c_str()) << std::endl;
        //fluid.set_reference_state(
        //        parse_value(fluid_json, "Tref", false, 20+273.15) ,
        //        parse_value(fluid_json, "pref", false, 1.01325e5) ,
        //        parse_value(fluid_json, "xref", false, 0.0) ,
        //        parse_value(fluid_json, "href", false, 0.0) ,
        //        parse_value(fluid_json, "sref", false, 0.0)
        //        );

        /// A function to check coefficients and equation types.
        fluid.validate();
    } catch (std::exception& e) {
        std::cout << format("Unable to load fluid: %s; error was %s\n", fluid.getName().c_str(), e.what());
        throw;
    }

    const std::string name = fluid.getName();
    // These characters would corrupt the comma-joined fluid lists or the
    // "INCOMP::Name" / "Name[x]" fluid-string parsing.
    if (name.find_first_of(",|[]:&") != std::string::npos) {
        throw ValueError(format("Invalid incompressible fluid name [%s]: must not contain any of ',|[]:&'.", name.c_str()));
    }

    std::map<std::string, std::size_t>::const_iterator it = string_to_index_map.find(name);
    if (it != string_to_index_map.end()) {
        // Re-adding an existing name replaces the fluid in place (idempotent
        // re-registration and edit flows); the name vectors must not grow
        // duplicates. A pure/solution flip would leave the name in the wrong
        // list, so reject that instead of silently misfiling it.
        if (fluid_map[it->second].is_pure() != fluid.is_pure()) {
            throw ValueError(
              format("Cannot replace incompressible fluid [%s]: pure/solution classification differs from the existing entry.", name.c_str()));
        }
        fluid_map[it->second] = std::move(fluid);
        return;
    }

    const std::size_t index = fluid_map.size();
    string_to_index_map[name] = index;
    if (fluid.is_pure()) {
        this->name_vector_pure.push_back(name);
    } else {
        this->name_vector_solution.push_back(name);
    }
    fluid_map[index] = std::move(fluid);
};

void JSONIncompressibleLibrary::add_obj(const IncompressibleFluid& fluid_obj) {
    _is_empty = false;

    // Get the next index for this fluid
    std::size_t index = fluid_map.size();

    // Add index->fluid mapping
    fluid_map[index] = fluid_obj;

    // Create an instance of the fluid
    IncompressibleFluid& fluid = fluid_map[index];

    /// A function to check coefficients and equation types.
    fluid.validate();

    // Add name->index mapping
    string_to_index_map[fluid.getName()] = index;
}

// Get an IncompressibleFluid instance stored in this library
IncompressibleFluid& JSONIncompressibleLibrary::get(const std::string& key) {
    // Try to find it
    auto it = string_to_index_map.find(key);
    // If it is found
    if (it != string_to_index_map.end()) {
        return get(it->second);
    } else {
        throw ValueError(format("key [%s] was not found in string_to_index_map in JSONIncompressibleLibrary", key.c_str()));
    }
};

/// Get a IncompressibleFluid instance stored in this library
/**
 @param key The index of the fluid in the map
 */
IncompressibleFluid& JSONIncompressibleLibrary::get(std::size_t key) {
    // Try to find it
    auto it = fluid_map.find(key);
    // If it is found
    if (it != fluid_map.end()) {
        return it->second;
    } else {
        throw ValueError(format("key [%d] was not found in JSONIncompressibleLibrary", key));
    }
};

static JSONIncompressibleLibrary library;

void load_incompressible_library();

// Thread-safe lazy initialization — see FluidLibrary.cpp for the same pattern
// and rationale (gh-2787).
static std::once_flag library_load_flag;

static void ensure_library_loaded() {
    std::call_once(library_load_flag, &load_incompressible_library);
}

void load_incompressible_library() {
    // This json formatted string comes from the all_incompressibles_JSON.h header
    nlohmann::json dd = cpjson::parse(all_incompressibles_JSON);
    try {
        library.add_many(dd);
    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }
    // TODO: Implement LiBr in the source code!
    //library.add_obj(LiBrSolution());
}

JSONIncompressibleLibrary& get_incompressible_library() {
    ensure_library_loaded();
    return library;
}

IncompressibleFluid& get_incompressible_fluid(const std::string& fluid_string) {
    ensure_library_loaded();
    return library.get(fluid_string);
}

std::string get_incompressible_list_pure() {
    ensure_library_loaded();
    return library.get_incompressible_list_pure();
};
std::string get_incompressible_list_solution() {
    ensure_library_loaded();
    return library.get_incompressible_list_solution();
};

} /* namespace CoolProp */
