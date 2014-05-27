
#ifndef FLUIDLIBRARY_H
#define FLUIDLIBRARY_H

#include "CoolPropFluid.h"

#include "rapidjson/rapidjson_include.h"

#include <map>
#include <algorithm>

namespace CoolProp{

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
    /// Map from CAS code to JSON instance.  For pseudo-pure fluids, use name in place of CAS code since no CASE number is defined for mixtures
    std::map<std::size_t, CoolPropFluid> fluid_map;
    std::vector<std::string> name_vector;
    std::map<std::string, std::size_t> string_to_index_map;
    bool _is_empty;
protected:

    /// Parse the contributions to the residual Helmholtz energy
    void parse_alphar(rapidjson::Value &alphar, EquationOfState &EOS)
    {	
        for (rapidjson::Value::ValueIterator itr = alphar.Begin(); itr != alphar.End(); ++itr)
        {	
            // A reference for code cleanness
            rapidjson::Value &contribution = *itr;

            // Get the type (required!)
            std::string type = contribution["type"].GetString();

            if (!type.compare("ResidualHelmholtzPower"))
            {
                if (EOS.alphar.Power.N > 0){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<long double> l = cpjson::get_long_double_array(contribution["l"]);
                EOS.alphar.Power = ResidualHelmholtzPower(n,d,t,l);
            }
            else if (!type.compare("ResidualHelmholtzGaussian"))
            {
                if (EOS.alphar.Gaussian.N > 0){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<long double> eta = cpjson::get_long_double_array(contribution["eta"]);
                std::vector<long double> epsilon = cpjson::get_long_double_array(contribution["epsilon"]);
                std::vector<long double> beta = cpjson::get_long_double_array(contribution["beta"]);
                std::vector<long double> gamma = cpjson::get_long_double_array(contribution["gamma"]);
                EOS.alphar.Gaussian = ResidualHelmholtzGaussian(n,d,t,eta,epsilon,beta,gamma);
            }
            else if (!type.compare("ResidualHelmholtzNonAnalytic"))
            {
                if (EOS.alphar.NonAnalytic.N > 0){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> a = cpjson::get_long_double_array(contribution["a"]);
                std::vector<long double> b = cpjson::get_long_double_array(contribution["b"]);
                std::vector<long double> beta = cpjson::get_long_double_array(contribution["beta"]);
                std::vector<long double> A = cpjson::get_long_double_array(contribution["A"]);
                std::vector<long double> B = cpjson::get_long_double_array(contribution["B"]);
                std::vector<long double> C = cpjson::get_long_double_array(contribution["C"]);
                std::vector<long double> D = cpjson::get_long_double_array(contribution["D"]);
                EOS.alphar.NonAnalytic = ResidualHelmholtzNonAnalytic(n,a,b,beta,A,B,C,D);
            }
            else if (!type.compare("ResidualHelmholtzLemmon2005"))
            {
                if (EOS.alphar.Lemmon2005.N > 0){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<long double> l = cpjson::get_long_double_array(contribution["l"]);
                std::vector<long double> m = cpjson::get_long_double_array(contribution["m"]);
                EOS.alphar.Lemmon2005 = ResidualHelmholtzLemmon2005(n,d,t,l,m);
            }
            else if (!type.compare("ResidualHelmholtzExponential"))
            {
                if (EOS.alphar.Exponential.N > 0){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> d = cpjson::get_long_double_array(contribution["d"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<long double> g = cpjson::get_long_double_array(contribution["g"]);
                std::vector<long double> l = cpjson::get_long_double_array(contribution["l"]);
                EOS.alphar.Exponential = ResidualHelmholtzExponential(n,d,t,g,l);
            }
            else if (!type.compare("ResidualHelmholtzAssociating"))
            {
                if (EOS.alphar.SAFT.disabled == false){throw ValueError("Cannot add ");}
                long double a = cpjson::get_double(contribution,"a");
                long double m = cpjson::get_double(contribution,"m");
                long double epsilonbar = cpjson::get_double(contribution,"epsilonbar");
                long double vbarn = cpjson::get_double(contribution,"vbarn");
                long double kappabar = cpjson::get_double(contribution,"kappabar");
                EOS.alphar.SAFT = ResidualHelmholtzSAFTAssociating(a,m,epsilonbar,vbarn,kappabar);
            }
            else
            {
                throw ValueError(format("Unsupported Residual helmholtz type: %s",type.c_str()));
            }
        }
    };

    /// Parse the contributions to the ideal-gas Helmholtz energy
    void parse_alpha0(rapidjson::Value &alpha0, EquationOfState &EOS)
    {
        if (!alpha0.IsArray()){throw ValueError();}
        for (rapidjson::Value::ValueIterator itr = alpha0.Begin(); itr != alpha0.End(); ++itr)
        {	
            // A reference for code cleanness
            rapidjson::Value &contribution = *itr;

            // Get the type (required!)
            std::string type = contribution["type"].GetString();

            if (!type.compare("IdealGasHelmholtzLead"))
            {
                if (EOS.alpha0.Lead.is_enabled() == true){throw ValueError("Cannot add ");}
                long double a1 = cpjson::get_double(contribution,"a1");
                long double a2 = cpjson::get_double(contribution,"a2");
                EOS.alpha0.Lead = IdealHelmholtzLead(a1, a2);
            }
            else if (!type.compare("IdealGasHelmholtzPower"))
            {
                if (EOS.alpha0.Power.is_enabled() == true){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                EOS.alpha0.Power = IdealHelmholtzPower(n, t);
            }
            else if (!type.compare("IdealGasHelmholtzLogTau"))
            {
                if (EOS.alpha0.LogTau.is_enabled() == true){throw ValueError("Cannot add ");}
                long double a = cpjson::get_double(contribution,"a");
                EOS.alpha0.LogTau = IdealHelmholtzLogTau(a);
            }
            else if (!type.compare("IdealGasHelmholtzPlanckEinstein"))
            {
                if (EOS.alpha0.PlanckEinstein.is_enabled() == true){throw ValueError("Cannot add ");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                EOS.alpha0.PlanckEinstein = IdealHelmholtzPlanckEinstein(n, t);
            }
            else if (!type.compare("IdealGasHelmholtzPlanckEinstein2"))
            {
                if (EOS.alpha0.PlanckEinstein2.is_enabled() == true){throw ValueError("Cannot add");}
                std::vector<long double> n = cpjson::get_long_double_array(contribution["n"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                std::vector<long double> c = cpjson::get_long_double_array(contribution["c"]);
                EOS.alpha0.PlanckEinstein2 = IdealHelmholtzPlanckEinstein2(n, t, c);
            }
            else if (!type.compare("IdealGasHelmholtzEnthalpyEntropyOffset"))
            {
                if (EOS.alpha0.EnthalpyEntropyOffset.is_enabled() == true){throw ValueError("Cannot add ");}
                long double a1 = cpjson::get_double(contribution, "a1");
                long double a2 = cpjson::get_double(contribution, "a2");
                EOS.alpha0.EnthalpyEntropyOffset = IdealHelmholtzEnthalpyEntropyOffset(a1, a2);
            }
            else if (!type.compare("IdealGasHelmholtzCP0Constant"))
            {
                if (EOS.alpha0.CP0Constant.is_enabled() == true){throw ValueError("Cannot add ");}
                long double cp_over_R = cpjson::get_double(contribution, "cp_over_R");
                long double Tc = cpjson::get_double(contribution, "Tc");
                long double T0 = cpjson::get_double(contribution, "T0");
                EOS.alpha0.CP0Constant = IdealHelmholtzCP0Constant(cp_over_R, Tc, T0);
            }
            else if (!type.compare("IdealGasHelmholtzCP0PolyT"))
            {
                if (EOS.alpha0.CP0PolyT.is_enabled() == true){throw ValueError("Cannot add ");}
                std::vector<long double> c = cpjson::get_long_double_array(contribution["c"]);
                std::vector<long double> t = cpjson::get_long_double_array(contribution["t"]);
                long double Tc = cpjson::get_double(contribution, "Tc");
                long double T0 = cpjson::get_double(contribution, "T0");
                EOS.alpha0.CP0PolyT = IdealHelmholtzCP0PolyT(c, t, Tc, T0);
            }
            else if (!type.compare("IdealGasHelmholtzCP0AlyLee"))
            {
                if (EOS.alpha0.CP0AlyLee.is_enabled() == true){std::cout << "Cannot add IdealGasHelmholtzCP0AlyLee\n";}
                std::vector<long double> c = cpjson::get_long_double_array(contribution["c"]);
                long double Tc = cpjson::get_double(contribution, "Tc");
                long double T0 = cpjson::get_double(contribution, "T0");
                EOS.alpha0.CP0AlyLee = IdealHelmholtzCP0AlyLee(c, Tc, T0);
            }
            else
            {
                std::cout << format("Unsupported ideal-gas Helmholtz type: %s\n",type.c_str());
                //throw ValueError(format("Unsupported ideal-gas Helmholtz type: %s",type.c_str()));
            }
        }
    };

    /// Parse the environmental parameters (ODP, GWP, etc.)
    void parse_environmental(rapidjson::Value &json, CoolPropFluid &fluid)
    {
        fluid.environment.ASHRAE34 = cpjson::get_string(json,"ASHRAE34");
        fluid.environment.GWP20 = cpjson::get_double(json,"GWP20");
        fluid.environment.GWP100 = cpjson::get_double(json,"GWP100");
        fluid.environment.GWP500 = cpjson::get_double(json,"GWP500");
        fluid.environment.HH = cpjson::get_double(json,"HH");
        fluid.environment.FH = cpjson::get_double(json,"FH");
        fluid.environment.PH = cpjson::get_double(json,"PH");
        fluid.environment.ODP = cpjson::get_double(json,"ODP");
    }

    /// Parse the Equation of state JSON entry
    void parse_EOS(rapidjson::Value &EOS_json, CoolPropFluid &fluid)
    {
        EquationOfState E;
        fluid.EOSVector.push_back(E);

        EquationOfState &EOS = fluid.EOSVector.at(fluid.EOSVector.size()-1);

        // Universal gas constant [J/mol/K]
        EOS.R_u = cpjson::get_double(EOS_json,"gas_constant");
        EOS.molar_mass = cpjson::get_double(EOS_json,"molar_mass");
        EOS.accentric = cpjson::get_double(EOS_json,"accentric");
        EOS.Ttriple = cpjson::get_double(EOS_json, "Ttriple");
        EOS.ptriple = cpjson::get_double(EOS_json, "ptriple");
        EOS.rhoLtriple = cpjson::get_double(EOS_json, "rhoLtriple");
        EOS.rhoVtriple = cpjson::get_double(EOS_json, "rhoVtriple");
        EOS.pseudo_pure = cpjson::get_bool(EOS_json, "pseudo_pure");
        EOS.limits.Tmax = cpjson::get_double(EOS_json, "T_max");
        EOS.limits.pmax = cpjson::get_double(EOS_json, "p_max");

        rapidjson::Value &reducing_state = EOS_json["reducing_state"];
        
        // Reducing state
        EOS.reduce.T = cpjson::get_double(reducing_state,"T");
        EOS.reduce.rhomolar = cpjson::get_double(reducing_state,"rhomolar");
        EOS.reduce.p = cpjson::get_double(reducing_state,"p");

        // BibTex keys
        EOS.BibTeX_EOS = cpjson::get_string(EOS_json,"BibTeX_EOS");
        EOS.BibTeX_CP0 = cpjson::get_string(EOS_json,"BibTeX_CP0");
        
        parse_alphar(EOS_json["alphar"], EOS);
        parse_alpha0(EOS_json["alpha0"], EOS);
        
        // Validate the equation of state that was just created
        EOS.validate();
        
    }

    /// Parse the list of possible equations of state
    void parse_EOS_listing(rapidjson::Value &EOS_array, CoolPropFluid & fluid)
    {
        for (rapidjson::Value::ValueIterator itr = EOS_array.Begin(); itr != EOS_array.End(); ++itr)
        {	
            parse_EOS(*itr,fluid);
        }

        // Set the EOS pointer to the first EOS
        fluid.pEOS = &(fluid.EOSVector[0]);
    };

    /// Parse the transport properties
    void parse_dilute_viscosity(rapidjson::Value &dilute, CoolPropFluid & fluid)
    {
        if (dilute.HasMember("hardcoded")){
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("Ethane")){
                fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_ETHANE; return;
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("collision_integral")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteGasCollisionIntegralData &CI = fluid.transport.viscosity_dilute.collision_integral;

            // Set the type flag
            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);
            CI.molar_mass = cpjson::get_double(dilute, "molar_mass");
            CI.C = cpjson::get_double(dilute, "C");
        }
        else if (!type.compare("kinetic_theory")){
            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_KINETIC_THEORY;
        }
        else if (!type.compare("powers_of_T")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteGasPowersOfT &CI = fluid.transport.viscosity_dilute.powers_of_T;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);

            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_POWERS_OF_T;
        }
        else if (!type.compare("collision_integral_powers_of_Tstar")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityDiluteCollisionIntegralPowersOfTstarData &CI = fluid.transport.viscosity_dilute.collision_integral_powers_of_Tstar;

            // Load up the values
            CI.a = cpjson::get_long_double_array(dilute["a"]);
            CI.t = cpjson::get_long_double_array(dilute["t"]);
            CI.T_reducing = cpjson::get_double(dilute,"T_reducing");
            CI.C = cpjson::get_double(dilute,"C");

            fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL_POWERS_OF_TSTAR ;
        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(),fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_initial_density_viscosity(rapidjson::Value &dilute, CoolPropFluid & fluid)
    {
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("Rainwater-Friend")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ViscosityRainWaterFriendData &RF = fluid.transport.viscosity_initial.rainwater_friend;

            // Load up the values
            RF.b = cpjson::get_long_double_array(dilute["b"]);
            RF.t = cpjson::get_long_double_array(dilute["t"]);
        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(),fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_higher_order_viscosity(rapidjson::Value &higher, CoolPropFluid & fluid)
    {
        // First check for hardcoded higher-order term
        if (higher.HasMember("hardcoded")){
            std::string target = cpjson::get_string(higher,"hardcoded");
            if (!target.compare("Hydrogen")){
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HYDROGEN; return;
            }
            else if (!target.compare("n-Hexane")){
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HEXANE; return;
            }
            else if (!target.compare("Ethane")){
                fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_ETHANE; return;
            }
            else{
                throw ValueError(format("hardcoded higher order viscosity term [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }

        std::string type = cpjson::get_string(higher, "type");
        if (!type.compare("modified_Batschinski_Hildebrand")){
            // Get a reference to the entry in the fluid instance to simplify the code that follows
            CoolProp::ViscosityModifiedBatschinskiHildebrandData &BH = fluid.transport.viscosity_higher_order.modified_Batschinski_Hildebrand;
            
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
        }
        else if (!type.compare("friction_theory")){
            // Get a reference to the entry in the fluid instance to simplify the code that follows
            CoolProp::ViscosityFrictionTheoryData &F = fluid.transport.viscosity_higher_order.friction_theory;
            
            // Set the flag for the type of this model
            fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_FRICTION_THEORY;

            // Always need these terms
            F.Ai = cpjson::get_long_double_array(higher["Ai"]);
            F.Aa = cpjson::get_long_double_array(higher["Aa"]);
            F.Aaa = cpjson::get_long_double_array(higher["Aaa"]);
            F.Ar = cpjson::get_long_double_array(higher["Ar"]);
            

            F.Na = cpjson::get_integer(higher,"Na");
            F.Naa = cpjson::get_integer(higher,"Naa");
            F.Nr = cpjson::get_integer(higher,"Nr");
            F.Nrr = cpjson::get_integer(higher,"Nrr");
            F.c1 = cpjson::get_double(higher,"c1");
            F.c2 = cpjson::get_double(higher,"c2");
            assert(F.Aa.size() == 3);
            assert(F.Aaa.size() == 3);
            assert(F.Ar.size() == 3);
            
            F.T_reduce = cpjson::get_double(higher,"T_reduce");

            if (higher.HasMember("Arr") && !higher.HasMember("Adrdr")){
                F.Arr = cpjson::get_long_double_array(higher["Arr"]);
                assert(F.Arr.size() == 3);
            }
            else if (higher.HasMember("Adrdr") && !higher.HasMember("Arr")){
                F.Adrdr = cpjson::get_long_double_array(higher["Adrdr"]);
                assert(F.Adrdr.size() == 3);
            }
            else{
                throw ValueError(format("can only provide one of Arr or Adrdr for fluid %s",fluid.name.c_str()));
            }

            if (higher.HasMember("Aaaa") && higher.HasMember("Arrr") && higher.HasMember("Aii")){
                F.Aaaa = cpjson::get_long_double_array(higher["Aaaa"]);
                F.Arrr = cpjson::get_long_double_array(higher["Arrr"]);
                F.Aii = cpjson::get_long_double_array(higher["Aii"]);
                F.Naaa = cpjson::get_integer(higher,"Naaa");
                F.Nrrr = cpjson::get_integer(higher,"Nrrr");
                F.Nii = cpjson::get_integer(higher,"Nii");
            }

        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(), fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_viscosity(rapidjson::Value &viscosity, CoolPropFluid & fluid)
    {
        // Load the BibTeX key
        fluid.transport.BibTeX_viscosity = cpjson::get_string(viscosity,"BibTeX");
        if (viscosity.HasMember("hardcoded")){
            std::string target = cpjson::get_string(viscosity,"hardcoded");
            if (!target.compare("Water")){
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_WATER; return;
            }
            else if (!target.compare("Helium")){
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_HELIUM; return;
            }
            else if (!target.compare("R23")){
                fluid.transport.hardcoded_viscosity = CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_R23; return;
            }
            else{
                throw ValueError(format("hardcoded viscosity [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }

        // Set the Lennard-Jones 12-6 potential variables, or approximate them from method of Chung
        if (!viscosity.HasMember("sigma_eta")|| !viscosity.HasMember("epsilon_over_k")){
            default_transport(fluid);
        }
        else{
            fluid.transport.sigma_eta = cpjson::get_double(viscosity, "sigma_eta");
            fluid.transport.epsilon_over_k = cpjson::get_double(viscosity, "epsilon_over_k");
        }

        // Load dilute viscosity term
        if (viscosity.HasMember("dilute")){
            parse_dilute_viscosity(viscosity["dilute"], fluid);
        }
        // Load initial density term
        if (viscosity.HasMember("initial_density")){
            parse_initial_density_viscosity(viscosity["initial_density"], fluid);
        }
        // Load higher_order term
        if (viscosity.HasMember("higher_order")){
            parse_higher_order_viscosity(viscosity["higher_order"], fluid);
        }
    };
    
    /// Parse the transport properties
    void parse_dilute_conductivity(rapidjson::Value &dilute, CoolPropFluid & fluid)
    {
        if (dilute.HasMember("hardcoded")){
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("CO2")){
                fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_CO2; return;
            }
            else if (!target.compare("Ethane")){
                fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETHANE; return;
            }
            else{
                throw ValueError(format("hardcoded dilute conductivity term [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("ratio_of_polynomials")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityDiluteRatioPolynomialsData &data = fluid.transport.conductivity_dilute.ratio_polynomials;

            // Set the type flag
            fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_RATIO_POLYNOMIALS;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.B = cpjson::get_long_double_array(dilute["B"]);
            data.n = cpjson::get_long_double_array(dilute["n"]);
            data.m = cpjson::get_long_double_array(dilute["m"]);
            data.T_reducing = cpjson::get_double(dilute, "T_reducing");
        }
        else if (!type.compare("eta0_and_poly")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityDiluteEta0AndPolyData &data = fluid.transport.conductivity_dilute.eta0_and_poly;

            // Set the type flag
            fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETA0_AND_POLY;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);           
        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(),fluid.name.c_str()));
        }
    };

    /// Parse the transport properties
    void parse_residual_conductivity(rapidjson::Value &dilute, CoolPropFluid & fluid)
    {
        if (dilute.HasMember("hardcoded")){
            std::string target = cpjson::get_string(dilute, "hardcoded");
            if (!target.compare("CO2")){
                fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_CO2; return;
            }
            else{
                throw ValueError(format("hardcoded residual conductivity term [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(dilute, "type");
        if (!type.compare("polynomial")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityResidualPolynomialData &data = fluid.transport.conductivity_residual.polynomials;

            // Set the type flag
            fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL;

            // Load up the values
            data.B = cpjson::get_long_double_array(dilute["B"]);
            data.d = cpjson::get_long_double_array(dilute["d"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);
            data.T_reducing = cpjson::get_double(dilute, "T_reducing");
            data.rhomass_reducing = cpjson::get_double(dilute, "rhomass_reducing");
        }
        else if (!type.compare("polynomial_and_exponential")){
            // Get a reference to the entry in the fluid instance
            CoolProp::ConductivityResidualPolynomialAndExponentialData &data = fluid.transport.conductivity_residual.polynomial_and_exponential;

            // Set the type flag
            fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL_AND_EXPONENTIAL;

            // Load up the values
            data.A = cpjson::get_long_double_array(dilute["A"]);
            data.d = cpjson::get_long_double_array(dilute["d"]);
            data.t = cpjson::get_long_double_array(dilute["t"]);
            data.gamma = cpjson::get_long_double_array(dilute["gamma"]);
            data.l = cpjson::get_long_double_array(dilute["l"]);
        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(),fluid.name.c_str()));
        }
    };

    void parse_critical_conductivity(rapidjson::Value &critical, CoolPropFluid & fluid)
    {
        if (critical.HasMember("hardcoded")){
            std::string target = cpjson::get_string(critical, "hardcoded");
            if (!target.compare("R123")){
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_R123; return;
            }
            else if (!target.compare("Ammonia")){
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_AMMONIA; return;
            }
            else if (!target.compare("None")){
                fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_NONE; return;
            }
            else{
                throw ValueError(format("critical conductivity term [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }
        std::string type = cpjson::get_string(critical, "type");
        if (!type.compare("simplified_Olchowy_Sengers")){
            //// Get a reference to the entry in the fluid instance
            CoolProp::ConductivityCriticalSimplifiedOlchowySengersData &data = fluid.transport.conductivity_critical.Olchowy_Sengers;

            // Set the type flag
            fluid.transport.conductivity_critical.type = CoolProp::ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_SIMPLIFIED_OLCHOWY_SENGERS;

            // Set values if they are found - otherwise fall back to default values
            if (critical.HasMember("qD")){ data.qD = cpjson::get_double(critical,"qD"); }
            if (critical.HasMember("zeta0")){ data.zeta0 = cpjson::get_double(critical,"zeta0"); }
            if (critical.HasMember("GAMMA")){ data.GAMMA = cpjson::get_double(critical,"GAMMA"); }
            if (critical.HasMember("gamma")){ data.gamma = cpjson::get_double(critical,"gamma"); }
            if (critical.HasMember("R0")){ data.R0 = cpjson::get_double(critical,"R0"); }
            if (critical.HasMember("T_ref")){ data.T_ref = cpjson::get_double(critical,"T_ref"); }
        }
        else{
            throw ValueError(format("type [%s] is not understood for fluid %s",type.c_str(),fluid.name.c_str()));
        }
    };
    
    /// Parse the thermal conductivity data
    void parse_thermal_conductivity(rapidjson::Value &conductivity, CoolPropFluid & fluid)
    {
        if (conductivity.HasMember("hardcoded")){
            std::string target = cpjson::get_string(conductivity, "hardcoded");
            if (!target.compare("Water")){
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_WATER; return;
            }
            else if (!target.compare("R23")){
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_R23; return;
            }
            else if (!target.compare("Helium")){
                fluid.transport.hardcoded_conductivity = CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_HELIUM; return;
            }
            else{
                throw ValueError(format("hardcoded residual conductivity term [%s] is not understood for fluid %s",target.c_str(), fluid.name.c_str()));
            }
        }

        // Load the BibTeX key
        fluid.transport.BibTeX_conductivity = cpjson::get_string(conductivity,"BibTeX");

        // Load dilute conductivity term
        if (conductivity.HasMember("dilute")){
            parse_dilute_conductivity(conductivity["dilute"], fluid);
        }
         // Load residual conductivity term
        if (conductivity.HasMember("residual")){
            parse_residual_conductivity(conductivity["residual"], fluid);
        }
        // Load critical conductivity term
        if (conductivity.HasMember("critical")){
            parse_critical_conductivity(conductivity["critical"], fluid);
        }
    };

    /// Parse the transport properties
    void parse_transport(rapidjson::Value &transport, CoolPropFluid & fluid)
    {
        //if (!fluid.name.compare("n-Hexane")){
        //    rapidjson::Document dd;
        //    dd.SetObject();

        //    dd.AddMember("core",transport,dd.GetAllocator());
        //    rapidjson::StringBuffer buffer;
        //    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);

        //    dd.Accept(writer);
        //    std::string json0 = buffer.GetString();
        //    //std::cout << json0 << std::endl;

        //    FILE *fp;
        //    fp = fopen("nHexane_transport.json","w");
        //    fprintf(fp,"%s",json0.c_str());
        //    fclose(fp);
        //}

        // Parse viscosity
        if (transport.HasMember("viscosity")){
            parse_viscosity(transport["viscosity"],fluid);
        }

        // Parse thermal conductivity 
        if (transport.HasMember("conductivity")){
            parse_thermal_conductivity(transport["conductivity"],fluid);
        }
    };

    void default_transport(CoolPropFluid & fluid)
    {
        // Use the method of Chung to approximate the values for epsilon_over_k and sigma_eta
        // Chung, T.-H.; Ajlan, M.; Lee, L. L.; Starling, K. E. Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties. Ind. Eng. Chem. Res. 1988, 27, 671-679.
        // rhoc needs to be in mol/L to yield a sigma in nm, 
        long double rho_crit_molar = fluid.pEOS->reduce.rhomolar/1000.0;// [mol/m3 to mol/L]
        long double Tc = fluid.pEOS->reduce.T;
        fluid.transport.sigma_eta = 0.809/pow(rho_crit_molar, static_cast<long double>(1.0/3.0))/1e9; // 1e9 is to convert from nm to m
        fluid.transport.epsilon_over_k = Tc/1.3593; // [K]
    }

    /// Parse the critical state for the given EOS
    void parse_crit_state(rapidjson::Value &alphar)
    {
    };

    /// Parse the critical state for the given EOS
    void parse_ancillaries(rapidjson::Value &ancillaries, CoolPropFluid & fluid)
    {
        if (!ancillaries.HasMember("pL") || !ancillaries.HasMember("pV") || !ancillaries.HasMember("rhoL") || !ancillaries.HasMember("rhoV")){throw ValueError("Ancillary curves are missing");};
        fluid.ancillaries.pL = SaturationAncillaryFunction(ancillaries["pL"]);
        fluid.ancillaries.pV = SaturationAncillaryFunction(ancillaries["pV"]);
        fluid.ancillaries.rhoL = SaturationAncillaryFunction(ancillaries["rhoL"]);
        fluid.ancillaries.rhoV = SaturationAncillaryFunction(ancillaries["rhoV"]);
    };

    /// Parse the surface_tension
    void parse_surface_tension(rapidjson::Value &surface_tension, CoolPropFluid & fluid)
    {
        fluid.ancillaries.surface_tension = SurfaceTensionCorrelation(surface_tension);
    };

    /// Validate the fluid file that was just constructed
    void validate(CoolPropFluid & fluid)
    {
        assert(fluid.EOSVector.size() > 0);
        assert(fluid.CAS.length() > 0);
        assert(fluid.name.length() > 0);
    }
public:
    
    // Default constructor;
    JSONFluidLibrary(){
        _is_empty = true;
    };
    bool is_empty(void){ return _is_empty;};

    /// Add all the fluid entries in the rapidjson::Value instance passed in
    void add_many(rapidjson::Value &listing)
    {
        for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr)
        {	
            add_one(*itr);
        }
    };
    void add_one(rapidjson::Value &fluid_json)
    {
        _is_empty = false;
        
        // Get the next index for this fluid
        std::size_t index = fluid_map.size();

        // Add index->fluid mapping
        fluid_map[index] = CoolPropFluid();

        // Create an instance of the fluid
        CoolPropFluid &fluid = fluid_map[index];

        // Fluid name
        fluid.name = fluid_json["NAME"].GetString(); name_vector.push_back(fluid.name);
        // CAS number
        fluid.CAS = fluid_json["CAS"].GetString();

        // Aliases
        fluid.aliases = cpjson::get_string_array(fluid_json["ALIASES"]);
        
        // EOS
        parse_EOS_listing(fluid_json["EOS"],fluid);

        // Validate the fluid
        validate(fluid);

        // Ancillaries for saturation 
        if (!fluid_json.HasMember("ANCILLARIES")){throw ValueError(format("Ancillary curves are missing for fluid [%s]",fluid.name.c_str()));};
        parse_ancillaries(fluid_json["ANCILLARIES"],fluid);

        // Surface tension
        if (!(fluid_json["ANCILLARIES"].HasMember("surface_tension"))){
            if (get_debug_level() > 0){
                std::cout << format("Surface tension curves are missing for fluid [%s]\n", fluid.name.c_str()) ;
            }
        }
        else{
            parse_surface_tension(fluid_json["ANCILLARIES"]["surface_tension"], fluid);
        }

        // Parse the environmental parameters
        if (!(fluid_json.HasMember("ENVIRONMENTAL"))){
            if (get_debug_level() > 0){
                std::cout << format("Environmental data are missing for fluid [%s]\n", fluid.name.c_str()) ;
            }
        }
        else{
            parse_environmental(fluid_json["ENVIRONMENTAL"], fluid);
        }

        // Parse the environmental parameters
        if (!(fluid_json.HasMember("TRANSPORT"))){
            default_transport(fluid);
        }
        else{
            parse_transport(fluid_json["TRANSPORT"], fluid);
        }
        
        // If the fluid is ok...

        // Add CAS->index mapping
        string_to_index_map[fluid.CAS] = index;

        // Add name->index mapping
        string_to_index_map[fluid.name] = index;
        
        // Add the aliases
        for (std::size_t i = 0; i < fluid.aliases.size(); ++i)
        {
            string_to_index_map[fluid.aliases[i]] = index;
        }

    };
    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key Either a CAS number or the name (CAS number should be preferred)
    */
    CoolPropFluid& get(std::string key)
    {
        std::map<std::string, std::size_t>::iterator it;
        // Try to find it
        it = string_to_index_map.find(key);
        // If it is found
        if (it != string_to_index_map.end()){
            return get(it->second);
        }
        else{
            throw ValueError(format("key [%s] was not found in string_to_index_map in JSONFluidLibrary",key.c_str()));
        }
    };
    /// Get a CoolPropFluid instance stored in this library
    /**
    @param key The index of the fluid in the map
    */
    CoolPropFluid& get(std::size_t key)
    {
        std::map<std::size_t,CoolPropFluid>::iterator it;
        // Try to find it
        it = fluid_map.find(key);
        // If it is found
        if (it != fluid_map.end()){
            return it->second;
        }
        else{
            throw ValueError(format("key [%d] was not found in JSONFluidLibrary",key));
        }
    };
    /// Return a comma-separated list of fluid names
    std::string get_fluid_list(void)
    {
        return strjoin(name_vector, ",");
    };
};

/// Get a reference to the library instance
JSONFluidLibrary & get_library(void);

/// Get a comma-separated-list of fluids that are included
std::string get_fluid_list(void);

/// Get the fluid structure returned as a reference
CoolPropFluid& get_fluid(std::string fluid_string);

} /* namespace CoolProp */
#endif