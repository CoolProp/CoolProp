
#ifndef INCOMPRESSIBLELIBRARY_H
#define INCOMPRESSIBLELIBRARY_H

#include "DataStructures.h"
#include "IncompressibleFluid.h"
//#include "crossplatform_shared_ptr.h"

#include "rapidjson_include.h"

#include <map>
#include <algorithm>

namespace CoolProp {

// Forward declaration of the necessary debug function to avoid including the whole header
extern int get_debug_level();

///// Class to access Lithium-Bromide solutions
///** Employs some basic wrapper-like functionality
// *  to bridge the gap between the solution functions
// *  used in the paper by Pátek and Klomfar:
// *  http://dx.doi.org/10.1016/j.ijrefrig.2005.10.007
// *
// *  We owe gratitude to the authors for providing
// *  both access to the paper as well as the equations
// *  in the form of C source code. */
//class LiBrSolution : public IncompressibleFluid{
//
//protected:
//    static double const M_H2O; /* kg/mol, molar mass of H2O */
//    static double const M_LiBr; /* kg/mol, molar mass of LiBr */
//    static double const T0; /* K, constant */
//
//    /* Critical point of H2O */
//    static double const Tc_H2O; /* K, temperature  */
//    static double const pc_H2O; /* MPa, pressure */
//    static double const rhoc_H2O; /* mol/m^3 (322 kg/m^3), molar density */
//    static double const hc_H2O; /* J/mol, molar enthalpy */
//    static double const sc_H2O; /* J/(mol.K) molar entropy*/
//
//    /*Triple point of H2O */
//    static double const Tt_H2O; /* K, temperature */
//    static double const cpt_H2O; /* J/(mol.K), molar isobaric heat capacity*/
//
//    double ps_mix(double T, double x);
//    double rho_mix(double T, double x);
//    double cp_mix(double T, double x);
//    double h_mix(double T, double x);
//    double s_mix(double T, double x);
//    double ps_H2O(double T);
//    double rho_H2O(double T);
//    double cp_H2O(double T);
//    double h_H2O(double T);
//    double s_H2O(double T);
//
//    /** Finished with the code from the paper. Now we need to
//     *  convert the molar values to mass-based units. */
//    double massToMole(double w);
//    double molarToSpecific(double w, double value);
//
//    static const bool debug;
//
//public:
//
//    LiBrSolution();
//
//    double rho(double T, double p, double x);
//    double c(double T, double p, double x);
//    //double h(double T_K, double p, double x);
//    double s(double T, double p, double x);
//    double visc(double T, double p, double x);
//    double cond(double T, double p, double x);
//    double u(double T, double p, double x);
//    double psat(double T, double x);
//    double Tfreeze(double p, double x);
//
//    /* Some functions can be inverted directly, those are listed
//     * here. It is also possible to solve for other quantities, but
//     * that involves some more sophisticated processing and is not
//     * done here, but in the backend, T(h,p) for example.
//     */
//    /// Temperature as a function of density, pressure and composition.
//    double T_rho (double Dmass, double p, double x){throw NotImplementedError(format("%s (%d): T from density is not implemented for LiBr.",__FILE__,__LINE__));}
//    /// Temperature as a function of heat capacities as a function of temperature, pressure and composition.
//    double T_c   (double Cmass, double p, double x){throw NotImplementedError(format("%s (%d): T from heat capacity is not implemented for LiBr.",__FILE__,__LINE__));}
//    /// Temperature as a function of entropy as a function of temperature, pressure and composition.
//    double T_s   (double Smass, double p, double x){throw NotImplementedError(format("%s (%d): T from entropy is not implemented for LiBr.",__FILE__,__LINE__));}
//    /// Temperature as a function of internal energy as a function of temperature, pressure and composition.
//    double T_u   (double Umass, double p, double x){throw NotImplementedError(format("%s (%d): T from internal energy is not implemented for LiBr.",__FILE__,__LINE__));}
//    /// Temperature as a function of enthalpy, pressure and composition.
//    //double T_h   (double Hmass, double p, double x){throw NotImplementedError(format("%s (%d): T from enthalpy is not implemented in the fluid, use the backend.",__FILE__,__LINE__));}
//    /// Viscosity as a function of temperature, pressure and composition.
//    double T_visc(double  visc, double p, double x){throw NotImplementedError(format("%s (%d): T from viscosity is not implemented.",__FILE__,__LINE__));}
//    /// Thermal conductivity as a function of temperature, pressure and composition.
//    double T_cond(double  cond, double p, double x){throw NotImplementedError(format("%s (%d): T from conductivity is not implemented.",__FILE__,__LINE__));}
//    /// Saturation pressure as a function of temperature and composition.
//    double T_psat(double  psat,           double x){throw NotImplementedError(format("%s (%d): T from psat is not implemented.",__FILE__,__LINE__));}
//    /// Composition as a function of freezing temperature and pressure.
//    double x_Tfreeze(       double Tfreeze, double p){throw NotImplementedError(format("%s (%d): x from T_freeze is not implemented.",__FILE__,__LINE__));}
//
//
//    /// Overwrite some standard functions that cannot be used with LiBr
//    void setName(std::string name){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setDescription(std::string description){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setReference(std::string reference){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setTmax(double Tmax){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setTmin(double Tmin){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setxmax(double xmax){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setxmin(double xmin){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setTminPsat(double TminPsat){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//
//    void setTbase(double Tbase){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setxbase(double xbase){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//
//    void setDensity(IncompressibleData density){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setSpecificHeat(IncompressibleData specific_heat){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setViscosity(IncompressibleData viscosity){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setConductivity(IncompressibleData conductivity){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setPsat(IncompressibleData p_sat){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setTfreeze(IncompressibleData T_freeze){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setVolToMass(IncompressibleData volToMass){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//    void setMassToMole(IncompressibleData massToMole){throw ValueError(format("%s (%d): Cannot change property of LiBr class",__FILE__,__LINE__));}
//
//    bool is_pure() {return false;};
//
//};
//

/// A container for the fluid parameters for the incompressible fluids
/**
This container holds copies of all of the fluid instances for the fluids that are loaded in incompressible.
New fluids can be added by passing in a rapidjson::Value instance to the add_one function, or
a rapidjson array of fluids to the add_many function.
*/

//typedef shared_ptr<IncompressibleFluid> IncompressibleFluidPointer;

class JSONIncompressibleLibrary
{
    /// Map from CAS code to JSON instance.
    /** This is not practical for the incompressibles, the CAS may not be
     *  defined for blends of heat transfer fluids and solutions.
     */
    std::map<std::size_t, IncompressibleFluid> fluid_map;
    std::vector<std::string> name_vector_pure, name_vector_solution;
    std::map<std::string, std::size_t> string_to_index_map;
    bool _is_empty;

   protected:
    /// A general function to parse the json files that hold the coefficient matrices
    IncompressibleData parse_coefficients(rapidjson::Value& obj, const std::string& id, bool vital);
    double parse_value(rapidjson::Value& obj, const std::string& id, bool vital, double def);
    composition_types parse_ifrac(rapidjson::Value& obj, const std::string& id);

   public:
    // Default constructor;
    JSONIncompressibleLibrary();
    ~JSONIncompressibleLibrary();

    bool is_empty(void) {
        return _is_empty;
    };

    /// Add all the fluid entries in the rapidjson::Value instance passed in
    void add_many(rapidjson::Value& listing);
    void add_one(rapidjson::Value& fluid_json);
    void add_obj(const IncompressibleFluid& fluid_obj);

    /** \brief Get an IncompressibleFluid instance stored in this library
     *
     * @param name Name of the fluid
     */
    IncompressibleFluid& get(const std::string& name);

    /** \brief Get a CoolPropFluid instance stored in this library
     *
     * @param key The index of the fluid in the map
     */
    IncompressibleFluid& get(std::size_t key);

    /// Return a comma-separated list of incompressible pure fluid names
    std::string get_incompressible_list_pure(void) {
        return strjoin(name_vector_pure, ",");
    };
    /// Return a comma-separated list of solution names
    std::string get_incompressible_list_solution(void) {
        return strjoin(name_vector_solution, ",");
    };
};

/// Get a reference to the library instance
JSONIncompressibleLibrary& get_incompressible_library(void);

/// Return a comma-separated list of incompressible pure fluid names
std::string get_incompressible_list_pure(void);

/// Return a comma-separated list of solution names
std::string get_incompressible_list_solution(void);

/// Get the fluid structure returned as a reference
IncompressibleFluid& get_incompressible_fluid(const std::string& fluid_string);

} /* namespace CoolProp */
#endif
