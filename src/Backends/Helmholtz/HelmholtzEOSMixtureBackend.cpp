/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#include <memory>

#if defined(_MSC_VER)
#    define _CRTDBG_MAP_ALLOC
#    ifndef _CRT_SECURE_NO_WARNINGS
#        define _CRT_SECURE_NO_WARNINGS
#    endif
#    include <crtdbg.h>
#    include <sys/stat.h>
#else
#    include <sys/stat.h>
#endif

#include <string>
//#include "CoolProp.h"

#include "HelmholtzEOSMixtureBackend.h"
#include "HelmholtzEOSBackend.h"
#include "Fluids/FluidLibrary.h"
#include "Solvers.h"
#include "MatrixMath.h"
#include "VLERoutines.h"
#include "FlashRoutines.h"
#include "TransportRoutines.h"
#include "MixtureDerivatives.h"
#include "PhaseEnvelopeRoutines.h"
#include "ReducingFunctions.h"
#include "MixtureParameters.h"
#include "IdealCurves.h"
#include "MixtureParameters.h"
#include <stdlib.h>

static int deriv_counter = 0;

namespace CoolProp {

class HEOSGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        if (fluid_names.size() == 1) {
            return new HelmholtzEOSBackend(fluid_names[0]);
        } else {
            return new HelmholtzEOSMixtureBackend(fluid_names);
        }
    };
};
// This static initialization will cause the generator to register
static CoolProp::GeneratorInitializer<HEOSGenerator> heos_gen(CoolProp::HEOS_BACKEND_FAMILY);

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend() {
    imposed_phase_index = iphase_not_imposed;
    is_pure_or_pseudopure = false;
    N = 0;
    _phase = iphase_unknown;
    // Reset the residual Helmholtz energy class
    residual_helmholtz.reset(new ResidualHelmholtz());
}
HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(const std::vector<std::string>& component_names, bool generate_SatL_and_SatV) {
    std::vector<CoolPropFluid> components(component_names.size());
    for (unsigned int i = 0; i < components.size(); ++i) {
        components[i] = get_library().get(component_names[i]);
    }

    // Reset the residual Helmholtz energy class
    residual_helmholtz.reset(new ResidualHelmholtz());

    // Set the components and associated flags
    set_components(components, generate_SatL_and_SatV);

    // Set the phase to default unknown value
    _phase = iphase_unknown;
}
HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(const std::vector<CoolPropFluid>& components, bool generate_SatL_and_SatV) {

    // Reset the residual Helmholtz energy class
    residual_helmholtz.reset(new ResidualHelmholtz());

    // Set the components and associated flags
    set_components(components, generate_SatL_and_SatV);

    // Set the phase to default unknown value
    _phase = iphase_unknown;
}
void HelmholtzEOSMixtureBackend::set_components(const std::vector<CoolPropFluid>& components, bool generate_SatL_and_SatV) {

    // Copy the components
    this->components = components;
    this->N = components.size();

    is_pure_or_pseudopure = (components.size() == 1);
    if (is_pure_or_pseudopure) {
        mole_fractions = std::vector<CoolPropDbl>(1, 1);
        std::vector<std::vector<double>> ones(1, std::vector<double>(1, 1));
        Reducing = shared_ptr<ReducingFunction>(new GERG2008ReducingFunction(components, ones, ones, ones, ones));
    } else {
        // Set the mixture parameters - binary pair reducing functions, departure functions, F_ij, etc.
        set_mixture_parameters();
    }

    imposed_phase_index = iphase_not_imposed;

    // Top-level class can hold copies of the base saturation classes,
    // saturation classes cannot hold copies of the saturation classes
    if (generate_SatL_and_SatV) {
        SatL.reset(get_copy(false));
        SatL->specify_phase(iphase_liquid);
        linked_states.push_back(SatL);
        SatV.reset(get_copy(false));
        SatV->specify_phase(iphase_gas);
        linked_states.push_back(SatV);
    }
}
void HelmholtzEOSMixtureBackend::set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
    if (mole_fractions.size() != N) {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]", mole_fractions.size(), N));
    }
    // Copy values without reallocating memory
    this->mole_fractions = mole_fractions;  // Most effective copy
    this->resize(N);                        // No reallocation of this->mole_fractions happens
    // Also store the mole fractions as doubles
    this->mole_fractions_double = std::vector<double>(mole_fractions.begin(), mole_fractions.end());
    clear_comp_change();
};
void HelmholtzEOSMixtureBackend::sync_linked_states(const HelmholtzEOSMixtureBackend* const source) {
    residual_helmholtz.reset(source->residual_helmholtz->copy_ptr());
    if (source->Reducing) {
        Reducing.reset(source->Reducing->copy());
    }
    // Recurse into linked states of the class
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        it->get()->sync_linked_states(source);
    }
}
HelmholtzEOSMixtureBackend* HelmholtzEOSMixtureBackend::get_copy(bool generate_SatL_and_SatV) {
    // Set up the class with these components
    HelmholtzEOSMixtureBackend* ptr = new HelmholtzEOSMixtureBackend(components, generate_SatL_and_SatV);
    // Recursively walk into linked states, setting the departure and reducing terms
    // to be equal to the parent (this instance)
    ptr->sync_linked_states(this);
    return ptr;
};
void HelmholtzEOSMixtureBackend::set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
    if (mass_fractions.size() != N) {
        throw ValueError(format("size of mass fraction vector [%d] does not equal that of component vector [%d]", mass_fractions.size(), N));
    }
    std::vector<CoolPropDbl> moles;
    CoolPropDbl sum_moles = 0.0;
    CoolPropDbl tmp = 0.0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        tmp = mass_fractions[i] / components[i].molar_mass();
        moles.push_back(tmp);
        sum_moles += tmp;
    }
    std::vector<CoolPropDbl> mole_fractions;
    for (std::vector<CoolPropDbl>::iterator it = moles.begin(); it != moles.end(); ++it) {
        mole_fractions.push_back(*it / sum_moles);
    }
    this->set_mole_fractions(mole_fractions);
};
void HelmholtzEOSMixtureBackend::resize(std::size_t N) {
    this->mole_fractions.resize(N);
    this->mole_fractions_double.resize(N);
    this->K.resize(N);
    this->lnK.resize(N);
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        it->get()->N = N;
        it->get()->resize(N);
    }
}
void HelmholtzEOSMixtureBackend::recalculate_singlephase_phase() {
    if (p() > p_critical()) {
        if (T() > T_critical()) {
            _phase = iphase_supercritical;
        } else {
            _phase = iphase_supercritical_liquid;
        }
    } else {
        if (T() > T_critical()) {
            _phase = iphase_supercritical_gas;
        } else {
            // Liquid or vapor
            if (rhomolar() > rhomolar_critical()) {
                _phase = iphase_liquid;
            } else {
                _phase = iphase_gas;
            }
        }
    }
}
std::string HelmholtzEOSMixtureBackend::fluid_param_string(const std::string& ParamName) {
    CoolProp::CoolPropFluid cpfluid = get_components()[0];
    if (!ParamName.compare("name")) {
        return cpfluid.name;
    } else if (!ParamName.compare("aliases")) {
        return strjoin(cpfluid.aliases, get_config_string(LIST_STRING_DELIMITER));
    } else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number")) {
        return cpfluid.CAS;
    } else if (!ParamName.compare("formula")) {
        return cpfluid.formula;
    } else if (!ParamName.compare("ASHRAE34")) {
        return cpfluid.environment.ASHRAE34;
    } else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname")) {
        return cpfluid.REFPROPname;
    } else if (ParamName.find("BibTeX") == 0)  // Starts with "BibTeX"
    {
        std::vector<std::string> parts = strsplit(ParamName, '-');
        if (parts.size() != 2) {
            throw ValueError(format("Unable to parse BibTeX string %s", ParamName.c_str()));
        }
        std::string key = parts[1];
        if (!key.compare("EOS")) {
            return cpfluid.EOS().BibTeX_EOS;
        } else if (!key.compare("CP0")) {
            return cpfluid.EOS().BibTeX_CP0;
        } else if (!key.compare("VISCOSITY")) {
            return cpfluid.transport.BibTeX_viscosity;
        } else if (!key.compare("CONDUCTIVITY")) {
            return cpfluid.transport.BibTeX_conductivity;
        } else if (!key.compare("ECS_LENNARD_JONES")) {
            throw NotImplementedError();
        } else if (!key.compare("ECS_VISCOSITY_FITS")) {
            throw NotImplementedError();
        } else if (!key.compare("ECS_CONDUCTIVITY_FITS")) {
            throw NotImplementedError();
        } else if (!key.compare("SURFACE_TENSION")) {
            return cpfluid.ancillaries.surface_tension.BibTeX;
        } else if (!key.compare("MELTING_LINE")) {
            return cpfluid.ancillaries.melting_line.BibTeX;
        } else {
            throw CoolProp::KeyError(format("Bad key to get_BibTeXKey [%s]", key.c_str()));
        }
    } else if (ParamName.find("pure") == 0) {
        if (is_pure()) {
            return "true";
        } else {
            return "false";
        }
    } else if (ParamName == "INCHI" || ParamName == "InChI" || ParamName == "INCHI_STRING") {
        return cpfluid.InChI;
    } else if (ParamName == "INCHI_Key" || ParamName == "InChIKey" || ParamName == "INCHIKEY") {
        return cpfluid.InChIKey;
    } else if (ParamName == "2DPNG_URL") {
        return cpfluid.TwoDPNG_URL;
    } else if (ParamName == "SMILES" || ParamName == "smiles") {
        return cpfluid.smiles;
    } else if (ParamName == "CHEMSPIDER_ID") {
        return format("%d", cpfluid.ChemSpider_id);
    } else if (ParamName == "JSON") {
        return get_fluid_as_JSONstring(cpfluid.CAS);
    } else {
        throw ValueError(format("fluid parameter [%s] is invalid", ParamName.c_str()));
    }
}

void HelmholtzEOSMixtureBackend::apply_simple_mixing_rule(std::size_t i, std::size_t j, const std::string& model) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }
    if (model == "linear") {
        double Tc1 = get_fluid_constant(i, iT_critical), Tc2 = get_fluid_constant(j, iT_critical);
        double gammaT = 0.5 * (Tc1 + Tc2) / sqrt(Tc1 * Tc2);
        double rhoc1 = get_fluid_constant(i, irhomolar_critical), rhoc2 = get_fluid_constant(j, irhomolar_critical);
        double gammaV = 4.0 * (1 / rhoc1 + 1 / rhoc2) / pow(pow(rhoc1, -1.0 / 3.0) + pow(rhoc2, -1.0 / 3.0), 3);
        set_binary_interaction_double(i, j, "betaT", 1.0);
        set_binary_interaction_double(i, j, "gammaT", gammaT);
        set_binary_interaction_double(i, j, "betaV", 1.0);
        set_binary_interaction_double(i, j, "gammaV", gammaV);
    } else if (model == "Lorentz-Berthelot") {
        set_binary_interaction_double(i, j, "betaT", 1.0);
        set_binary_interaction_double(i, j, "gammaT", 1.0);
        set_binary_interaction_double(i, j, "betaV", 1.0);
        set_binary_interaction_double(i, j, "gammaV", 1.0);
    } else {
        throw ValueError(format("mixing rule [%s] is not understood", model.c_str()));
    }
}
/// Set binary mixture floating point parameter for this instance
void HelmholtzEOSMixtureBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                               const double value) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }
    if (parameter == "Fij") {
        residual_helmholtz->Excess.F[i][j] = value;
        residual_helmholtz->Excess.F[j][i] = value;
    } else {
        Reducing->set_binary_interaction_double(i, j, parameter, value);
    }
    /// Also set the parameters in the managed pointers for other states
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        it->get()->set_binary_interaction_double(i, j, parameter, value);
    }
};
/// Get binary mixture floating point parameter for this instance
double HelmholtzEOSMixtureBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }
    if (parameter == "Fij") {
        return residual_helmholtz->Excess.F[i][j];
    } else {
        return Reducing->get_binary_interaction_double(i, j, parameter);
    }
};
///// Get binary mixture string value
//std::string HelmholtzEOSMixtureBackend::get_binary_interaction_string(const std::string &CAS1, const std::string &CAS2, const std::string &parameter){
//    return CoolProp::get_mixture_binary_pair_data(CAS1, CAS2, parameter);
//}
/// Set binary mixture floating point parameter for this instance
void HelmholtzEOSMixtureBackend::set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                               const std::string& value) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }
    if (parameter == "function") {
        residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(get_departure_function(value));
        residual_helmholtz->Excess.DepartureFunctionMatrix[j][i].reset(get_departure_function(value));
    } else {
        throw ValueError(format("Cannot process this string parameter [%s] in set_binary_interaction_string", parameter.c_str()));
    }
    /// Also set the parameters in the managed pointers for other states
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        it->get()->set_binary_interaction_string(i, j, parameter, value);
    }
};

void HelmholtzEOSMixtureBackend::calc_change_EOS(const std::size_t i, const std::string& EOS_name) {

    if (i < components.size()) {
        CoolPropFluid& fluid = components[i];
        EquationOfState& EOS = fluid.EOSVector[0];

        if (EOS_name == "SRK" || EOS_name == "Peng-Robinson") {

            // Get the parameters for the cubic EOS
            CoolPropDbl Tc = EOS.reduce.T;
            CoolPropDbl pc = EOS.reduce.p;
            CoolPropDbl rhomolarc = EOS.reduce.rhomolar;
            CoolPropDbl acentric = EOS.acentric;
            CoolPropDbl R = 8.3144598;

            // Remove the residual part
            EOS.alphar.empty_the_EOS();
            // Set the contribution
            shared_ptr<AbstractCubic> ac;
            if (EOS_name == "SRK") {
                ac.reset(new SRK(Tc, pc, acentric, R));
            } else {
                ac.reset(new PengRobinson(Tc, pc, acentric, R));
            }
            ac->set_Tr(Tc);
            ac->set_rhor(rhomolarc);
            EOS.alphar.cubic = ResidualHelmholtzGeneralizedCubic(ac);
        } else if (EOS_name == "XiangDeiters") {

            // Get the parameters for the EOS
            CoolPropDbl Tc = EOS.reduce.T;
            CoolPropDbl pc = EOS.reduce.p;
            CoolPropDbl rhomolarc = EOS.reduce.rhomolar;
            CoolPropDbl acentric = EOS.acentric;
            CoolPropDbl R = 8.3144598;

            // Remove the residual part
            EOS.alphar.empty_the_EOS();
            // Set the Xiang & Deiters contribution
            EOS.alphar.XiangDeiters = ResidualHelmholtzXiangDeiters(Tc, pc, rhomolarc, acentric, R);
        }
    } else {
        throw ValueError(format("Index [%d] is invalid", i));
    }
    // Now do the same thing to the saturated liquid and vapor instances if possible
    if (this->SatL) SatL->change_EOS(i, EOS_name);
    if (this->SatV) SatV->change_EOS(i, EOS_name);
}
void HelmholtzEOSMixtureBackend::calc_phase_envelope(const std::string& type) {
    // Clear the phase envelope data
    PhaseEnvelope = PhaseEnvelopeData();
    // Build the phase envelope
    PhaseEnvelopeRoutines::build(*this, type);
    // Finalize the phase envelope
    PhaseEnvelopeRoutines::finalize(*this);
};
void HelmholtzEOSMixtureBackend::set_mixture_parameters() {
    // Build the matrix of binary-pair reducing functions
    MixtureParameters::set_mixture_parameters(*this);
}
void HelmholtzEOSMixtureBackend::update_states(void) {
    CoolPropFluid& component = components[0];
    EquationOfState& EOS = component.EOSVector[0];

    // Clear the state class
    clear();

    // Calculate the new enthalpy and entropy values
    update(DmolarT_INPUTS, EOS.hs_anchor.rhomolar, EOS.hs_anchor.T);
    EOS.hs_anchor.hmolar = hmolar();
    EOS.hs_anchor.smolar = smolar();

    // Calculate the new enthalpy and entropy values at the reducing state
    update(DmolarT_INPUTS, EOS.reduce.rhomolar, EOS.reduce.T);
    EOS.reduce.hmolar = hmolar();
    EOS.reduce.smolar = smolar();

    // Clear again just to be sure
    clear();
}
const CoolProp::SimpleState& HelmholtzEOSMixtureBackend::calc_state(const std::string& state) {
    if (is_pure_or_pseudopure) {
        if (!state.compare("hs_anchor")) {
            return components[0].EOS().hs_anchor;
        } else if (!state.compare("max_sat_T")) {
            return components[0].EOS().max_sat_T;
        } else if (!state.compare("max_sat_p")) {
            return components[0].EOS().max_sat_p;
        } else if (!state.compare("reducing")) {
            return components[0].EOS().reduce;
        } else if (!state.compare("critical")) {
            return components[0].crit;
        } else if (!state.compare("triple_liquid")) {
            return components[0].triple_liquid;
        } else if (!state.compare("triple_vapor")) {
            return components[0].triple_vapor;
        } else {
            throw ValueError(format("This state [%s] is invalid to calc_state", state.c_str()));
        }
    } else {
        if (!state.compare("critical")) {
            return _critical;
        } else {
            throw ValueError(format("calc_state not supported for mixtures"));
        }
    }
};
CoolPropDbl HelmholtzEOSMixtureBackend::calc_acentric_factor(void) {
    if (is_pure_or_pseudopure) {
        return components[0].EOS().acentric;
    } else {
        throw ValueError("acentric factor cannot be calculated for mixtures");
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_gas_constant(void) {
    if (is_pure_or_pseudopure) {
        return components[0].gas_constant();
    } else {
        if (get_config_bool(NORMALIZE_GAS_CONSTANTS)) {
            return get_config_double(R_U_CODATA);
        } else {
            // mass fraction weighted average of the components
            double summer = 0;
            for (unsigned int i = 0; i < components.size(); ++i) {
                summer += mole_fractions[i] * components[i].gas_constant();
            }
            return summer;
        }
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_molar_mass(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].molar_mass();
    }
    return summer;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_saturation_ancillary(parameters param, int Q, parameters given, double value) {
    if (is_pure_or_pseudopure) {
        if (param == iP && given == iT) {
            // p = f(T), direct evaluation
            switch (Q) {
                case 0:
                    return components[0].ancillaries.pL.evaluate(value);
                case 1:
                    return components[0].ancillaries.pV.evaluate(value);
            }
        } else if (param == iT && given == iP) {
            // T = f(p), inverse evaluation
            switch (Q) {
                case 0:
                    return components[0].ancillaries.pL.invert(value);
                case 1:
                    return components[0].ancillaries.pV.invert(value);
            }
        } else if (param == iDmolar && given == iT) {
            // rho = f(T), inverse evaluation
            switch (Q) {
                case 0:
                    return components[0].ancillaries.rhoL.evaluate(value);
                case 1:
                    return components[0].ancillaries.rhoV.evaluate(value);
            }
        } else if (param == iT && given == iDmolar) {
            // T = f(rho), inverse evaluation
            switch (Q) {
                case 0:
                    return components[0].ancillaries.rhoL.invert(value);
                case 1:
                    return components[0].ancillaries.rhoV.invert(value);
            }
        } else if (param == isurface_tension && given == iT) {
            return components[0].ancillaries.surface_tension.evaluate(value);
        } else {
            throw ValueError(format("calc of %s given %s is invalid in calc_saturation_ancillary", get_parameter_information(param, "short").c_str(),
                                    get_parameter_information(given, "short").c_str()));
        }

        throw ValueError(format("Q [%d] is invalid in calc_saturation_ancillary", Q));
    } else {
        throw NotImplementedError(format("calc_saturation_ancillary not implemented for mixtures"));
    }
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_melting_line(int param, int given, CoolPropDbl value) {
    if (is_pure_or_pseudopure) {
        return components[0].ancillaries.melting_line.evaluate(param, given, value);
    } else {
        throw NotImplementedError(format("calc_melting_line not implemented for mixtures"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_surface_tension(void) {
    if (is_pure_or_pseudopure) {
        if ((_phase == iphase_twophase) || (_phase == iphase_critical_point)) {  // if within the two phase region or at critical point
            return components[0].ancillaries.surface_tension.evaluate(T());      //    calculate surface tension and return
        } else {                                                                 // else state point not in the two phase region
            throw ValueError(format("surface tension is only defined within the two-phase region; Try PQ or QT inputs"));  // throw error
        }
    } else {
        throw NotImplementedError(format("surface tension not implemented for mixtures"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_viscosity_dilute(void) {
    if (is_pure_or_pseudopure) {
        CoolPropDbl eta_dilute;
        switch (components[0].transport.viscosity_dilute.type) {
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_KINETIC_THEORY:
                eta_dilute = TransportRoutines::viscosity_dilute_kinetic_theory(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL:
                eta_dilute = TransportRoutines::viscosity_dilute_collision_integral(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_POWERS_OF_T:
                eta_dilute = TransportRoutines::viscosity_dilute_powers_of_T(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_POWERS_OF_TR:
                eta_dilute = TransportRoutines::viscosity_dilute_powers_of_Tr(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_COLLISION_INTEGRAL_POWERS_OF_TSTAR:
                eta_dilute = TransportRoutines::viscosity_dilute_collision_integral_powers_of_T(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_ETHANE:
                eta_dilute = TransportRoutines::viscosity_dilute_ethane(*this);
                break;
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_CYCLOHEXANE:
                eta_dilute = TransportRoutines::viscosity_dilute_cyclohexane(*this);
                break;
            default:
                throw ValueError(
                  format("dilute viscosity type [%d] is invalid for fluid %s", components[0].transport.viscosity_dilute.type, name().c_str()));
        }
        return eta_dilute;
    } else {
        throw NotImplementedError(format("dilute viscosity not implemented for mixtures"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_viscosity_background() {
    CoolPropDbl eta_dilute = calc_viscosity_dilute(), initial_density = 0, residual = 0;
    return calc_viscosity_background(eta_dilute, initial_density, residual);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_viscosity_background(CoolPropDbl eta_dilute, CoolPropDbl& initial_density, CoolPropDbl& residual) {

    switch (components[0].transport.viscosity_initial.type) {
        case ViscosityInitialDensityVariables::VISCOSITY_INITIAL_DENSITY_RAINWATER_FRIEND: {
            CoolPropDbl B_eta_initial = TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend(*this);
            CoolPropDbl rho = rhomolar();
            initial_density = eta_dilute * B_eta_initial * rho;
            break;
        }
        case ViscosityInitialDensityVariables::VISCOSITY_INITIAL_DENSITY_EMPIRICAL: {
            initial_density = TransportRoutines::viscosity_initial_density_dependence_empirical(*this);
            break;
        }
        case ViscosityInitialDensityVariables::VISCOSITY_INITIAL_DENSITY_NOT_SET: {
            break;
        }
    }

    // Higher order terms
    switch (components[0].transport.viscosity_higher_order.type) {
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_BATSCHINKI_HILDEBRAND:
            residual = TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_FRICTION_THEORY:
            residual = TransportRoutines::viscosity_higher_order_friction_theory(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HYDROGEN:
            residual = TransportRoutines::viscosity_hydrogen_higher_order_hardcoded(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_TOLUENE:
            residual = TransportRoutines::viscosity_toluene_higher_order_hardcoded(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HEXANE:
            residual = TransportRoutines::viscosity_hexane_higher_order_hardcoded(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_HEPTANE:
            residual = TransportRoutines::viscosity_heptane_higher_order_hardcoded(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_ETHANE:
            residual = TransportRoutines::viscosity_ethane_higher_order_hardcoded(*this);
            break;
        case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_BENZENE:
            residual = TransportRoutines::viscosity_benzene_higher_order_hardcoded(*this);
            break;
        default:
            throw ValueError(
              format("higher order viscosity type [%d] is invalid for fluid %s", components[0].transport.viscosity_dilute.type, name().c_str()));
    }

    return initial_density + residual;
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_viscosity(void) {
    if (is_pure_or_pseudopure) {
        CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0;
        calc_viscosity_contributions(dilute, initial_density, residual, critical);
        return dilute + initial_density + residual + critical;
    } else {
        set_warning_string("Mixture model for viscosity is highly approximate");
        CoolPropDbl summer = 0;
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
            shared_ptr<HelmholtzEOSBackend> HEOS(new HelmholtzEOSBackend(components[i]));
            HEOS->update(DmolarT_INPUTS, _rhomolar, _T);
            summer += mole_fractions[i] * log(HEOS->viscosity());
        }
        return exp(summer);
    }
}
void HelmholtzEOSMixtureBackend::calc_viscosity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual,
                                                              CoolPropDbl& critical) {
    if (is_pure_or_pseudopure) {
        // Reset the variables
        dilute = 0;
        initial_density = 0;
        residual = 0;
        critical = 0;

        // Get a reference for code cleanness
        CoolPropFluid& component = components[0];

        if (!component.transport.viscosity_model_provided) {
            throw ValueError(format("Viscosity model is not available for this fluid"));
        }

        // Check if using ECS
        if (component.transport.viscosity_using_ECS) {
            // Get reference fluid name
            std::string fluid_name = component.transport.viscosity_ecs.reference_fluid;
            std::vector<std::string> names(1, fluid_name);
            // Get a managed pointer to the reference fluid for ECS
            shared_ptr<HelmholtzEOSMixtureBackend> ref_fluid(new HelmholtzEOSMixtureBackend(names));
            // Get the viscosity using ECS and stick in the critical value
            critical = TransportRoutines::viscosity_ECS(*this, *ref_fluid);
            return;
        }

        // Check if using Chung model
        if (component.transport.viscosity_using_Chung) {
            // Get the viscosity using ECS and stick in the critical value
            critical = TransportRoutines::viscosity_Chung(*this);
            return;
        }

        // Check if using rho*sr model
        if (component.transport.viscosity_using_rhosr) {
            // Get the viscosity using rho*sr model and stick in the critical value
            critical = TransportRoutines::viscosity_rhosr(*this);
            return;
        }

        if (component.transport.hardcoded_viscosity != CoolProp::TransportPropertyData::VISCOSITY_NOT_HARDCODED) {
            // Evaluate hardcoded model and stick in the critical value
            switch (component.transport.hardcoded_viscosity) {
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_WATER:
                    critical = TransportRoutines::viscosity_water_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_HEAVYWATER:
                    critical = TransportRoutines::viscosity_heavywater_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_HELIUM:
                    critical = TransportRoutines::viscosity_helium_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_R23:
                    critical = TransportRoutines::viscosity_R23_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_METHANOL:
                    critical = TransportRoutines::viscosity_methanol_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_M_XYLENE:
                    critical = TransportRoutines::viscosity_m_xylene_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_O_XYLENE:
                    critical = TransportRoutines::viscosity_o_xylene_hardcoded(*this);
                    break;
                case CoolProp::TransportPropertyData::VISCOSITY_HARDCODED_P_XYLENE:
                    critical = TransportRoutines::viscosity_p_xylene_hardcoded(*this);
                    break;
                default:
                    throw ValueError(
                      format("hardcoded viscosity type [%d] is invalid for fluid %s", component.transport.hardcoded_viscosity, name().c_str()));
            }
            return;
        }

        // -------------------------
        //     Normal evaluation
        // -------------------------

        // Dilute part
        dilute = calc_viscosity_dilute();

        // Background viscosity given by the sum of the initial density dependence and higher order terms
        calc_viscosity_background(dilute, initial_density, residual);

        // Critical part (no fluids have critical enhancement for viscosity currently)
        critical = 0;
    } else {
        throw ValueError("calc_viscosity_contributions invalid for mixtures");
    }
}
void HelmholtzEOSMixtureBackend::calc_conductivity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual,
                                                                 CoolPropDbl& critical) {
    if (is_pure_or_pseudopure) {
        // Reset the variables
        dilute = 0;
        initial_density = 0;
        residual = 0;
        critical = 0;

        // Get a reference for code cleanness
        CoolPropFluid& component = components[0];

        if (!component.transport.conductivity_model_provided) {
            throw ValueError(format("Thermal conductivity model is not available for this fluid"));
        }

        // Check if using ECS
        if (component.transport.conductivity_using_ECS) {
            // Get reference fluid name
            std::string fluid_name = component.transport.conductivity_ecs.reference_fluid;
            std::vector<std::string> name(1, fluid_name);
            // Get a managed pointer to the reference fluid for ECS
            shared_ptr<HelmholtzEOSMixtureBackend> ref_fluid(new HelmholtzEOSMixtureBackend(name));
            // Get the viscosity using ECS and store in initial_density (not normally used);
            initial_density = TransportRoutines::conductivity_ECS(*this, *ref_fluid);  // Warning: not actually initial_density
            return;
        }

        if (component.transport.hardcoded_conductivity != CoolProp::TransportPropertyData::CONDUCTIVITY_NOT_HARDCODED) {
            // Evaluate hardcoded model and deposit in initial_density variable
            // Warning: not actually initial_density
            switch (component.transport.hardcoded_conductivity) {
                case CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_WATER:
                    initial_density = TransportRoutines::conductivity_hardcoded_water(*this);
                    break;
                case CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_HEAVYWATER:
                    initial_density = TransportRoutines::conductivity_hardcoded_heavywater(*this);
                    break;
                case CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_R23:
                    initial_density = TransportRoutines::conductivity_hardcoded_R23(*this);
                    break;
                case CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_HELIUM:
                    initial_density = TransportRoutines::conductivity_hardcoded_helium(*this);
                    break;
                case CoolProp::TransportPropertyData::CONDUCTIVITY_HARDCODED_METHANE:
                    initial_density = TransportRoutines::conductivity_hardcoded_methane(*this);
                    break;
                default:
                    throw ValueError(format("hardcoded conductivity type [%d] is invalid for fluid %s",
                                            components[0].transport.hardcoded_conductivity, name().c_str()));
            }
            return;
        }

        // -------------------------
        //     Normal evaluation
        // -------------------------

        // Dilute part
        switch (component.transport.conductivity_dilute.type) {
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_RATIO_POLYNOMIALS:
                dilute = TransportRoutines::conductivity_dilute_ratio_polynomials(*this);
                break;
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETA0_AND_POLY:
                dilute = TransportRoutines::conductivity_dilute_eta0_and_poly(*this);
                break;
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_CO2:
                dilute = TransportRoutines::conductivity_dilute_hardcoded_CO2(*this);
                break;
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_ETHANE:
                dilute = TransportRoutines::conductivity_dilute_hardcoded_ethane(*this);
                break;
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_NONE:
                dilute = 0.0;
                break;
            default:
                throw ValueError(
                  format("dilute conductivity type [%d] is invalid for fluid %s", components[0].transport.conductivity_dilute.type, name().c_str()));
        }

        // Residual part
        residual = calc_conductivity_background();

        // Critical part
        switch (component.transport.conductivity_critical.type) {
            case ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_SIMPLIFIED_OLCHOWY_SENGERS:
                critical = TransportRoutines::conductivity_critical_simplified_Olchowy_Sengers(*this);
                break;
            case ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_R123:
                critical = TransportRoutines::conductivity_critical_hardcoded_R123(*this);
                break;
            case ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_AMMONIA:
                critical = TransportRoutines::conductivity_critical_hardcoded_ammonia(*this);
                break;
            case ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_NONE:
                critical = 0.0;
                break;
            case ConductivityCriticalVariables::CONDUCTIVITY_CRITICAL_CARBONDIOXIDE_SCALABRIN_JPCRD_2006:
                critical = TransportRoutines::conductivity_critical_hardcoded_CO2_ScalabrinJPCRD2006(*this);
                break;
            default:
                throw ValueError(
                  format("critical conductivity type [%d] is invalid for fluid %s", components[0].transport.viscosity_dilute.type, name().c_str()));
        }
    } else {
        throw ValueError("calc_conductivity_contributions invalid for mixtures");
    }
};

CoolPropDbl HelmholtzEOSMixtureBackend::calc_conductivity_background(void) {
    // Residual part
    CoolPropDbl lambda_residual = _HUGE;
    switch (components[0].transport.conductivity_residual.type) {
        case ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL:
            lambda_residual = TransportRoutines::conductivity_residual_polynomial(*this);
            break;
        case ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_POLYNOMIAL_AND_EXPONENTIAL:
            lambda_residual = TransportRoutines::conductivity_residual_polynomial_and_exponential(*this);
            break;
        default:
            throw ValueError(
              format("residual conductivity type [%d] is invalid for fluid %s", components[0].transport.conductivity_residual.type, name().c_str()));
    }
    return lambda_residual;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_conductivity(void) {
    if (is_pure_or_pseudopure) {
        CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0;
        calc_conductivity_contributions(dilute, initial_density, residual, critical);
        return dilute + initial_density + residual + critical;
    } else {
        set_warning_string("Mixture model for conductivity is highly approximate");
        CoolPropDbl summer = 0;
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
            shared_ptr<HelmholtzEOSBackend> HEOS(new HelmholtzEOSBackend(components[i]));
            HEOS->update(DmolarT_INPUTS, _rhomolar, _T);
            summer += mole_fractions[i] * HEOS->conductivity();
        }
        return summer;
    }
}
void HelmholtzEOSMixtureBackend::calc_conformal_state(const std::string& reference_fluid, CoolPropDbl& T, CoolPropDbl& rhomolar) {
    shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> REF(new CoolProp::HelmholtzEOSBackend(reference_fluid));

    if (T < 0 && rhomolar < 0) {
        // Collect some parameters
        CoolPropDbl Tc = T_critical(), Tc0 = REF->T_critical(), rhocmolar = rhomolar_critical(), rhocmolar0 = REF->rhomolar_critical();

        // Starting guess values for shape factors
        CoolPropDbl theta = 1;
        CoolPropDbl phi = 1;

        // The equivalent substance reducing ratios
        CoolPropDbl f = Tc / Tc0 * theta;
        CoolPropDbl h = rhocmolar0 / rhocmolar * phi;  // Must be the ratio of MOLAR densities!!

        // Starting guesses for conformal state
        T = this->T() / f;
        rhomolar = this->rhomolar() * h;
    }

    TransportRoutines::conformal_state_solver(*this, *REF, T, rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Ttriple(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].EOS().Ttriple;
    }
    return summer;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_p_triple(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].EOS().ptriple;
    }
    return summer;
}
std::string HelmholtzEOSMixtureBackend::calc_name(void) {
    if (components.size() != 1) {
        throw ValueError(format("calc_name is only valid for pure and pseudo-pure fluids, %d components", components.size()));
    } else {
        return components[0].name;
    }
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_saturated_liquid_keyed_output(parameters key) {
    if ((key == iDmolar) && _rhoLmolar) return _rhoLmolar;
    if (!SatL) throw ValueError("The saturated liquid state has not been set.");
    return SatL->keyed_output(key);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_saturated_vapor_keyed_output(parameters key) {
    if ((key == iDmolar) && _rhoVmolar) return _rhoVmolar;
    if (!SatV) throw ValueError("The saturated vapor state has not been set.");
    return SatV->keyed_output(key);
}

void HelmholtzEOSMixtureBackend::calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) {
    if (type == "Joule-Thomson") {
        JouleThomsonCurveTracer JTCT(this, 1e5, 800);
        JTCT.trace(T, p);
    } else if (type == "Joule-Inversion") {
        JouleInversionCurveTracer JICT(this, 1e5, 800);
        JICT.trace(T, p);
    } else if (type == "Ideal") {
        IdealCurveTracer ICT(this, 1e5, 800);
        ICT.trace(T, p);
    } else if (type == "Boyle") {
        BoyleCurveTracer BCT(this, 1e5, 800);
        BCT.trace(T, p);
    } else {
        throw ValueError(format("Invalid ideal curve type: %s", type.c_str()));
    }
};
std::vector<std::string> HelmholtzEOSMixtureBackend::calc_fluid_names(void) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < components.size(); ++i) {
        out.push_back(components[i].name);
    }
    return out;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_ODP(void) {
    if (components.size() != 1) {
        throw ValueError(format("For now, calc_ODP is only valid for pure and pseudo-pure fluids, %d components", components.size()));
    } else {
        CoolPropDbl v = components[0].environment.ODP;
        if (!ValidNumber(v) || v < 0) {
            throw ValueError(format("ODP value is not specified or invalid"));
        }
        return v;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_GWP20(void) {
    if (components.size() != 1) {
        throw ValueError(format("For now, calc_GWP20 is only valid for pure and pseudo-pure fluids, %d components", components.size()));
    } else {
        CoolPropDbl v = components[0].environment.GWP20;
        if (!ValidNumber(v) || v < 0) {
            throw ValueError(format("GWP20 value is not specified or invalid"));
        }
        return v;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_GWP100(void) {
    if (components.size() != 1) {
        throw ValueError(format("For now, calc_GWP100 is only valid for pure and pseudo-pure fluids, %d components", components.size()));
    } else {
        CoolPropDbl v = components[0].environment.GWP100;
        if (!ValidNumber(v) || v < 0) {
            throw ValueError(format("GWP100 value is not specified or invalid"));
        }
        return v;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_GWP500(void) {
    if (components.size() != 1) {
        throw ValueError(format("For now, calc_GWP500 is only valid for pure and pseudo-pure fluids, %d components", components.size()));
    } else {
        CoolPropDbl v = components[0].environment.GWP500;
        if (!ValidNumber(v) || v < 0) {
            throw ValueError(format("GWP500 value is not specified or invalid"));
        }
        return v;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_T_critical(void) {
    if (components.size() != 1) {
        std::vector<CriticalState> critpts = calc_all_critical_points();
        if (critpts.size() == 1) {
            //if (!critpts[0].stable){ throw ValueError(format("found one critical point but critical point is not stable")); }
            return critpts[0].T;
        } else {
            throw ValueError(format("critical point finding routine found %d critical points", critpts.size()));
        }
    } else {
        return components[0].crit.T;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_p_critical(void) {
    if (components.size() != 1) {
        std::vector<CriticalState> critpts = calc_all_critical_points();
        if (critpts.size() == 1) {
            //if (!critpts[0].stable){ throw ValueError(format("found one critical point but critical point is not stable")); }
            return critpts[0].p;
        } else {
            throw ValueError(format("critical point finding routine found %d critical points", critpts.size()));
        }
    } else {
        return components[0].crit.p;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_rhomolar_critical(void) {
    if (components.size() != 1) {
        std::vector<CriticalState> critpts = calc_all_critical_points();
        if (critpts.size() == 1) {
            //if (!critpts[0].stable){ throw ValueError(format("found one critical point but critical point is not stable")); }
            return critpts[0].rhomolar;
        } else {
            throw ValueError(format("critical point finding routine found %d critical points", critpts.size()));
        }
    } else {
        return components[0].crit.rhomolar;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_pmax_sat(void) {
    if (is_pure_or_pseudopure) {
        if (components[0].EOS().pseudo_pure) {
            return components[0].EOS().max_sat_p.p;
        } else {
            return p_critical();
        }
    } else {
        throw ValueError("calc_pmax_sat not yet defined for mixtures");
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Tmax_sat(void) {
    if (is_pure_or_pseudopure) {
        if (components[0].EOS().pseudo_pure) {
            double Tmax_sat = components[0].EOS().max_sat_T.T;
            if (!ValidNumber(Tmax_sat)) {
                return T_critical();
            } else {
                return Tmax_sat;
            }
        } else {
            return T_critical();
        }
    } else {
        throw ValueError("calc_Tmax_sat not yet defined for mixtures");
    }
}

void HelmholtzEOSMixtureBackend::calc_Tmin_sat(CoolPropDbl& Tmin_satL, CoolPropDbl& Tmin_satV) {
    if (is_pure_or_pseudopure) {
        Tmin_satL = components[0].EOS().sat_min_liquid.T;
        Tmin_satV = components[0].EOS().sat_min_vapor.T;
        return;
    } else {
        throw ValueError("calc_Tmin_sat not yet defined for mixtures");
    }
}

void HelmholtzEOSMixtureBackend::calc_pmin_sat(CoolPropDbl& pmin_satL, CoolPropDbl& pmin_satV) {
    if (is_pure_or_pseudopure) {
        pmin_satL = components[0].EOS().sat_min_liquid.p;
        pmin_satV = components[0].EOS().sat_min_vapor.p;
        return;
    } else {
        throw ValueError("calc_pmin_sat not yet defined for mixtures");
    }
}

// Minimum allowed saturation temperature the maximum of the saturation temperatures of liquid and vapor
// For pure fluids, both values are the same, for pseudo-pure they are probably the same, for mixtures they are definitely not the same

CoolPropDbl HelmholtzEOSMixtureBackend::calc_Tmax(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].EOS().limits.Tmax;
    }
    return summer;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Tmin(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].EOS().limits.Tmin;
    }
    return summer;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_pmax(void) {
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        summer += mole_fractions[i] * components[i].EOS().limits.pmax;
    }
    return summer;
}

void HelmholtzEOSMixtureBackend::update_DmolarT_direct(CoolPropDbl rhomolar, CoolPropDbl T) {
    // TODO: This is just a quick fix for #878 - should be done more systematically
    const CoolPropDbl rhomolar_min = 0;
    const CoolPropDbl T_min = 0;

    if (rhomolar < rhomolar_min) {
        throw ValueError(format("The molar density of %f mol/m3 is below the minimum of %f mol/m3", rhomolar, rhomolar_min));
    }

    if (T < T_min) {
        throw ValueError(format("The temperature of %f K is below the minimum of %f K", T, T_min));
    }

    CoolProp::input_pairs pair = DmolarT_INPUTS;
    // Set up the state
    pre_update(pair, rhomolar, T);

    _rhomolar = rhomolar;
    _T = T;
    _p = calc_pressure();

    // Cleanup
    bool optional_checks = false;
    post_update(optional_checks);
}

void HelmholtzEOSMixtureBackend::update_HmolarQ_with_guessT(CoolPropDbl hmolar, CoolPropDbl Q, CoolPropDbl Tguess) {
    CoolProp::input_pairs pair = CoolProp::HmolarQ_INPUTS;
    // Set up the state
    pre_update(pair, hmolar, Q);

    _hmolar = hmolar;
    _Q = Q;
    FlashRoutines::HQ_flash(*this, Tguess);

    // Cleanup
    post_update();
}
void HelmholtzEOSMixtureBackend::update_internal(HelmholtzEOSMixtureBackend& HEOS) {
    this->_hmolar = HEOS.hmolar();
    this->_smolar = HEOS.smolar();
    this->_T = HEOS.T();
    this->_umolar = HEOS.umolar();
    this->_p = HEOS.p();
    this->_rhomolar = HEOS.rhomolar();
    this->_Q = HEOS.Q();
    this->_phase = HEOS.phase();

    // Copy the derivatives as well
}
void HelmholtzEOSMixtureBackend::update_TP_guessrho(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_guess) {
    CoolProp::input_pairs pair = PT_INPUTS;
    // Set up the state
    pre_update(pair, p, T);

    // Do the flash call to find rho = f(T,p)
    CoolPropDbl rhomolar = solver_rho_Tp(T, p, rhomolar_guess);

    // Update the class with the new calculated density
    update_DmolarT_direct(rhomolar, T);

    // Skip the cleanup, already done in update_DmolarT_direct
}

void HelmholtzEOSMixtureBackend::pre_update(CoolProp::input_pairs& input_pair, CoolPropDbl& value1, CoolPropDbl& value2) {
    // Clear the state
    clear();

    if (is_pure_or_pseudopure == false && mole_fractions.size() == 0) {
        throw ValueError("Mole fractions must be set");
    }

    // If the inputs are in mass units, convert them to molar units
    mass_to_molar_inputs(input_pair, value1, value2);

    // Set the mole-fraction weighted gas constant for the mixture
    // (or the pure/pseudo-pure fluid) if it hasn't been set yet
    gas_constant();

    // Calculate and cache the reducing state
    calc_reducing_state();
}

void HelmholtzEOSMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    if (get_debug_level() > 10) {
        std::cout << format("%s (%d): update called with (%d: (%s), %g, %g)", __FILE__, __LINE__, input_pair,
                            get_input_pair_short_desc(input_pair).c_str(), value1, value2)
                  << std::endl;
    }

    CoolPropDbl ld_value1 = value1, ld_value2 = value2;
    pre_update(input_pair, ld_value1, ld_value2);
    value1 = ld_value1;
    value2 = ld_value2;

    switch (input_pair) {
        case PT_INPUTS:
            _p = value1;
            _T = value2;
            FlashRoutines::PT_flash(*this);
            break;
        case DmolarT_INPUTS:
            _rhomolar = value1;
            _T = value2;
            FlashRoutines::DHSU_T_flash(*this, iDmolar);
            break;
        case SmolarT_INPUTS:
            _smolar = value1;
            _T = value2;
            FlashRoutines::DHSU_T_flash(*this, iSmolar);
            break;
        //case HmolarT_INPUTS:
        //    _hmolar = value1; _T = value2; FlashRoutines::DHSU_T_flash(*this, iHmolar); break;
        //case TUmolar_INPUTS:
        //    _T = value1; _umolar = value2; FlashRoutines::DHSU_T_flash(*this, iUmolar); break;
        case DmolarP_INPUTS:
            _rhomolar = value1;
            _p = value2;
            FlashRoutines::DP_flash(*this);
            break;
        case DmolarHmolar_INPUTS:
            _rhomolar = value1;
            _hmolar = value2;
            FlashRoutines::HSU_D_flash(*this, iHmolar);
            break;
        case DmolarSmolar_INPUTS:
            _rhomolar = value1;
            _smolar = value2;
            FlashRoutines::HSU_D_flash(*this, iSmolar);
            break;
        case DmolarUmolar_INPUTS:
            _rhomolar = value1;
            _umolar = value2;
            FlashRoutines::HSU_D_flash(*this, iUmolar);
            break;
        case HmolarP_INPUTS:
            _hmolar = value1;
            _p = value2;
            FlashRoutines::HSU_P_flash(*this, iHmolar);
            break;
        case PSmolar_INPUTS:
            _p = value1;
            _smolar = value2;
            FlashRoutines::HSU_P_flash(*this, iSmolar);
            break;
        case PUmolar_INPUTS:
            _p = value1;
            _umolar = value2;
            FlashRoutines::HSU_P_flash(*this, iUmolar);
            break;
        case HmolarSmolar_INPUTS:
            _hmolar = value1;
            _smolar = value2;
            FlashRoutines::HS_flash(*this);
            break;
        case QT_INPUTS:
            _Q = value1;
            _T = value2;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            FlashRoutines::QT_flash(*this);
            break;
        case PQ_INPUTS:
            _p = value1;
            _Q = value2;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            FlashRoutines::PQ_flash(*this);
            break;
        case QSmolar_INPUTS:
            _Q = value1;
            _smolar = value2;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            FlashRoutines::QS_flash(*this);
            break;
        case HmolarQ_INPUTS:
            _hmolar = value1;
            _Q = value2;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            FlashRoutines::HQ_flash(*this);
            break;
        case DmolarQ_INPUTS:
            _rhomolar = value1;
            _Q = value2;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            FlashRoutines::DQ_flash(*this);
            break;
        default:
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
    }

    post_update();
}
const std::vector<CoolPropDbl> HelmholtzEOSMixtureBackend::calc_mass_fractions() {
    // mass fraction is mass_i/total_mass;
    CoolPropDbl mm = molar_mass();
    std::vector<CoolPropDbl>& mole_fractions = get_mole_fractions_ref();
    std::vector<CoolPropDbl> mass_fractions(mole_fractions.size());
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
        double mmi = get_fluid_constant(i, imolar_mass);
        mass_fractions[i] = mmi * (mole_fractions[i]) / mm;
    }
    return mass_fractions;
}

void HelmholtzEOSMixtureBackend::update_with_guesses(CoolProp::input_pairs input_pair, double value1, double value2,
                                                     const GuessesStructure& guesses) {
    if (get_debug_level() > 10) {
        std::cout << format("%s (%d): update called with (%d: (%s), %g, %g)", __FILE__, __LINE__, input_pair,
                            get_input_pair_short_desc(input_pair).c_str(), value1, value2)
                  << std::endl;
    }

    CoolPropDbl ld_value1 = value1, ld_value2 = value2;
    pre_update(input_pair, ld_value1, ld_value2);
    value1 = ld_value1;
    value2 = ld_value2;

    switch (input_pair) {
        case PQ_INPUTS:
            _p = value1;
            _Q = value2;
            FlashRoutines::PQ_flash_with_guesses(*this, guesses);
            break;
        case QT_INPUTS:
            _Q = value1;
            _T = value2;
            FlashRoutines::QT_flash_with_guesses(*this, guesses);
            break;
        case PT_INPUTS:
            _p = value1;
            _T = value2;
            FlashRoutines::PT_flash_with_guesses(*this, guesses);
            break;
        default:
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
    }
    post_update();
}

void HelmholtzEOSMixtureBackend::post_update(bool optional_checks) {
    // Check the values that must always be set
    //if (_p < 0){ throw ValueError("p is less than zero");}
    if (!ValidNumber(_p)) {
        throw ValueError("p is not a valid number");
    }
    //if (_T < 0){ throw ValueError("T is less than zero");}
    if (!ValidNumber(_T)) {
        throw ValueError("T is not a valid number");
    }
    if (_rhomolar < 0) {
        throw ValueError("rhomolar is less than zero");
    }
    if (!ValidNumber(_rhomolar)) {
        throw ValueError("rhomolar is not a valid number");
    }

    if (optional_checks) {
        if (!ValidNumber(_Q)) {
            throw ValueError("Q is not a valid number");
        }
        if (_phase == iphase_unknown) {
            throw ValueError("_phase is unknown");
        }
    }

    // Set the reduced variables
    _tau = _reducing.T / _T;
    _delta = _rhomolar / _reducing.rhomolar;

    // Update the terms in the excess contribution
    residual_helmholtz->Excess.update(_tau, _delta);
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_Bvirial() {
    return 1 / rhomolar_reducing() * calc_alphar_deriv_nocache(0, 1, mole_fractions, _tau, 1e-12);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dBvirial_dT() {
    SimpleState red = get_reducing_state();
    CoolPropDbl dtau_dT = -red.T / pow(_T, 2);
    return 1 / red.rhomolar * calc_alphar_deriv_nocache(1, 1, mole_fractions, _tau, 1e-12) * dtau_dT;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Cvirial() {
    return 1 / pow(rhomolar_reducing(), 2) * calc_alphar_deriv_nocache(0, 2, mole_fractions, _tau, 1e-12);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dCvirial_dT() {
    SimpleState red = get_reducing_state();
    CoolPropDbl dtau_dT = -red.T / pow(_T, 2);
    return 1 / pow(red.rhomolar, 2) * calc_alphar_deriv_nocache(1, 2, mole_fractions, _tau, 1e-12) * dtau_dT;
}
void HelmholtzEOSMixtureBackend::p_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value, bool& saturation_called) {
    /*
    Determine the phase given p and one other state variable
    */
    saturation_called = false;

    // Reference declaration to save indexing
    CoolPropFluid& component = components[0];

    // Maximum saturation temperature - Equal to critical pressure for pure fluids
    CoolPropDbl psat_max = calc_pmax_sat();

    // Check supercritical pressure
    if (_p > psat_max) {
        _Q = 1e9;
        switch (other) {
            case iT: {
                if (_T > _crit.T) {
                    this->_phase = iphase_supercritical;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iDmolar: {
                if (_rhomolar < _crit.rhomolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iSmolar: {
                if (_smolar.pt() > _crit.smolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iHmolar: {
                if (_hmolar.pt() > _crit.hmolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iUmolar: {
                if (_umolar.pt() > _crit.umolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            default: {
                throw ValueError("supercritical pressure but other invalid for now");
            }
        }
    }
    // Check between triple point pressure and psat_max
    else if (_p >= components[0].EOS().ptriple * 0.9999 && _p <= psat_max) {
        // First try the ancillaries, use them to determine the state if you can

        // Calculate dew and bubble temps from the ancillaries (everything needs them)
        _TLanc = components[0].ancillaries.pL.invert(_p);
        _TVanc = components[0].ancillaries.pV.invert(_p);

        bool definitely_two_phase = false;

        // Try using the ancillaries for P,H,S if they are there
        switch (other) {
            case iT: {

                if (has_melting_line()) {
                    double Tm = melting_line(iT, iP, _p);
                    if (get_config_bool(DONT_CHECK_PROPERTY_LIMITS)) {
                        _phase = iphase_liquid;
                    } else {
                        if (_T < Tm - 0.001) {
                            throw ValueError(format("For now, we don't support T [%g K] below Tmelt(p) [%g K]", _T, Tm));
                        }
                    }
                } else {
                    if (get_config_bool(DONT_CHECK_PROPERTY_LIMITS)) {
                        _phase = iphase_liquid;
                    } else {
                        if (_T < Tmin() - 0.001) {
                            throw ValueError(format("For now, we don't support T [%g K] below Tmin(saturation) [%g K]", _T, Tmin()));
                        }
                    }
                }

                CoolPropDbl T_vap = 0.1 + static_cast<double>(_TVanc);
                CoolPropDbl T_liq = -0.1 + static_cast<double>(_TLanc);

                if (value > T_vap) {
                    this->_phase = iphase_gas;
                    _Q = -1000;
                    return;
                } else if (value < T_liq) {
                    this->_phase = iphase_liquid;
                    _Q = 1000;
                    return;
                }
                break;
            }
            case iHmolar: {
                if (!component.ancillaries.hL.enabled()) {
                    break;
                }
                // Ancillaries are h-h_anchor, so add back h_anchor
                CoolPropDbl h_liq = component.ancillaries.hL.evaluate(_TLanc) + component.EOS().hs_anchor.hmolar;
                CoolPropDbl h_liq_error_band = component.ancillaries.hL.get_max_abs_error();
                CoolPropDbl h_vap = h_liq + component.ancillaries.hLV.evaluate(_TLanc);
                CoolPropDbl h_vap_error_band = h_liq_error_band + component.ancillaries.hLV.get_max_abs_error();

                //                HelmholtzEOSMixtureBackend HEOS(components);
                //                HEOS.update(QT_INPUTS, 1, _TLanc);
                //                double h1 = HEOS.hmolar();
                //                HEOS.update(QT_INPUTS, 0, _TLanc);
                //                double h0 = HEOS.hmolar();

                // Check if in range given the accuracy of the fit
                if (value > h_vap + h_vap_error_band) {
                    this->_phase = iphase_gas;
                    _Q = -1000;
                    return;
                } else if (value < h_liq - h_liq_error_band) {
                    this->_phase = iphase_liquid;
                    _Q = 1000;
                    return;
                } else if (value > h_liq + h_liq_error_band && value < h_vap - h_vap_error_band) {
                    definitely_two_phase = true;
                }
                break;
            }
            case iSmolar: {
                if (!component.ancillaries.sL.enabled()) {
                    break;
                }
                // Ancillaries are s-s_anchor, so add back s_anchor
                CoolPropDbl s_anchor = component.EOS().hs_anchor.smolar;
                CoolPropDbl s_liq = component.ancillaries.sL.evaluate(_TLanc) + s_anchor;
                CoolPropDbl s_liq_error_band = component.ancillaries.sL.get_max_abs_error();
                CoolPropDbl s_vap = s_liq + component.ancillaries.sLV.evaluate(_TVanc);
                CoolPropDbl s_vap_error_band = s_liq_error_band + component.ancillaries.sLV.get_max_abs_error();

                // Check if in range given the accuracy of the fit
                if (value > s_vap + s_vap_error_band) {
                    this->_phase = iphase_gas;
                    _Q = -1000;
                    return;
                } else if (value < s_liq - s_liq_error_band) {
                    this->_phase = iphase_liquid;
                    _Q = 1000;
                    return;
                } else if (value > s_liq + s_liq_error_band && value < s_vap - s_vap_error_band) {
                    definitely_two_phase = true;
                }
                break;
            }
            case iUmolar: {
                if (!component.ancillaries.hL.enabled()) {
                    break;
                }
                // u = h-p/rho

                // Ancillaries are h-h_anchor, so add back h_anchor
                CoolPropDbl h_liq = component.ancillaries.hL.evaluate(_TLanc) + component.EOS().hs_anchor.hmolar;
                CoolPropDbl h_liq_error_band = component.ancillaries.hL.get_max_abs_error();
                CoolPropDbl h_vap = h_liq + component.ancillaries.hLV.evaluate(_TLanc);
                CoolPropDbl h_vap_error_band = h_liq_error_band + component.ancillaries.hLV.get_max_abs_error();
                CoolPropDbl rho_vap = component.ancillaries.rhoV.evaluate(_TVanc);
                CoolPropDbl rho_liq = component.ancillaries.rhoL.evaluate(_TLanc);
                CoolPropDbl u_liq = h_liq - _p / rho_liq;
                CoolPropDbl u_vap = h_vap - _p / rho_vap;
                CoolPropDbl u_liq_error_band = 1.5 * h_liq_error_band;  // Most of error is in enthalpy
                CoolPropDbl u_vap_error_band = 1.5 * h_vap_error_band;  // Most of error is in enthalpy

                // Check if in range given the accuracy of the fit
                if (value > u_vap + u_vap_error_band) {
                    this->_phase = iphase_gas;
                    _Q = -1000;
                    return;
                } else if (value < u_liq - u_liq_error_band) {
                    this->_phase = iphase_liquid;
                    _Q = 1000;
                    return;
                } else if (value > u_liq + u_liq_error_band && value < u_vap - u_vap_error_band) {
                    definitely_two_phase = true;
                }
                break;
            }
            default: {
            }
        }

        // Now either density is an input, or an ancillary for h,s,u is missing
        // Always calculate the densities using the ancillaries
        if (!definitely_two_phase) {
            _rhoVanc = component.ancillaries.rhoV.evaluate(_TVanc);
            _rhoLanc = component.ancillaries.rhoL.evaluate(_TLanc);
            CoolPropDbl rho_vap = 0.95 * static_cast<double>(_rhoVanc);
            CoolPropDbl rho_liq = 1.05 * static_cast<double>(_rhoLanc);
            switch (other) {
                case iDmolar: {
                    if (value < rho_vap) {
                        this->_phase = iphase_gas;
                        return;
                    } else if (value > rho_liq) {
                        this->_phase = iphase_liquid;
                        return;
                    }
                    break;
                }
            }
        }

        if (!is_pure_or_pseudopure) {
            throw ValueError("possibly two-phase inputs not supported for mixtures for now");
        }

        // Actually have to use saturation information sadly
        // For the given pressure, find the saturation state
        // Run the saturation routines to determine the saturation densities and pressures
        HelmholtzEOSMixtureBackend HEOS(components);
        HEOS._p = this->_p;
        HEOS._Q = 0;  // ?? What is the best to do here? Doesn't matter for our purposes since pure fluid
        FlashRoutines::PQ_flash(HEOS);

        // We called the saturation routines, so HEOS.SatL and HEOS.SatV are now updated
        // with the saturated liquid and vapor values, which can therefore be used in
        // the other solvers
        saturation_called = true;

        CoolPropDbl Q;

        if (other == iT) {
            if (value < HEOS.SatL->T() - 100 * DBL_EPSILON) {
                this->_phase = iphase_liquid;
                _Q = -1000;
                return;
            } else if (value > HEOS.SatV->T() + 100 * DBL_EPSILON) {
                this->_phase = iphase_gas;
                _Q = 1000;
                return;
            } else {
                this->_phase = iphase_twophase;
            }
        }
        switch (other) {
            case iDmolar:
                Q = (1 / value - 1 / HEOS.SatL->rhomolar()) / (1 / HEOS.SatV->rhomolar() - 1 / HEOS.SatL->rhomolar());
                break;
            case iSmolar:
                Q = (value - HEOS.SatL->smolar()) / (HEOS.SatV->smolar() - HEOS.SatL->smolar());
                break;
            case iHmolar:
                Q = (value - HEOS.SatL->hmolar()) / (HEOS.SatV->hmolar() - HEOS.SatL->hmolar());
                break;
            case iUmolar:
                Q = (value - HEOS.SatL->umolar()) / (HEOS.SatV->umolar() - HEOS.SatL->umolar());
                break;
            default:
                throw ValueError(format("bad input for other"));
        }
        // TODO: Check the speed penalty of these calls
        // Update the states
        if (this->SatL) this->SatL->update(DmolarT_INPUTS, HEOS.SatL->rhomolar(), HEOS.SatL->T());
        if (this->SatV) this->SatV->update(DmolarT_INPUTS, HEOS.SatV->rhomolar(), HEOS.SatV->T());
        // Update the two-Phase variables
        _rhoLmolar = HEOS.SatL->rhomolar();
        _rhoVmolar = HEOS.SatV->rhomolar();

        //
        if (Q < -1e-9) {
            this->_phase = iphase_liquid;
            _Q = -1000;
            return;
        } else if (Q > 1 + 1e-9) {
            this->_phase = iphase_gas;
            _Q = 1000;
            return;
        } else {
            this->_phase = iphase_twophase;
        }

        _Q = Q;
        // Load the outputs
        _T = _Q * HEOS.SatV->T() + (1 - _Q) * HEOS.SatL->T();
        _rhomolar = 1 / (_Q / HEOS.SatV->rhomolar() + (1 - _Q) / HEOS.SatL->rhomolar());
        return;
    } else if (_p < components[0].EOS().ptriple * 0.9999) {
        if (other == iT) {
            if (_T > std::max(Tmin(), Ttriple())) {
                _phase = iphase_gas;
            } else {
                if (get_config_bool(DONT_CHECK_PROPERTY_LIMITS)) {
                    _phase = iphase_gas;
                } else {
                    throw NotImplementedError(format("For now, we don't support p [%g Pa] below ptriple [%g Pa] when T [%g] is less than Tmin [%g]",
                                                     _p, components[0].EOS().ptriple, _T, std::max(Tmin(), Ttriple())));
                }
            }
        } else {
            _phase = iphase_gas;
        }
    } else {
        throw ValueError(format("The pressure [%g Pa] cannot be used in p_phase_determination", _p));
    }
}
void HelmholtzEOSMixtureBackend::calc_ssat_max(void) {
    class Residual : public FuncWrapper1D
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        Residual(HelmholtzEOSMixtureBackend& HEOS) : HEOS(&HEOS){};
        double call(double T) {
            HEOS->update(QT_INPUTS, 1, T);
            // dTdp_along_sat
            double dTdp_along_sat =
              HEOS->T() * (1 / HEOS->SatV->rhomolar() - 1 / HEOS->SatL->rhomolar()) / (HEOS->SatV->hmolar() - HEOS->SatL->hmolar());
            // dsdT_along_sat;
            return HEOS->SatV->first_partial_deriv(iSmolar, iT, iP) + HEOS->SatV->first_partial_deriv(iSmolar, iP, iT) / dTdp_along_sat;
        }
    };
    if (!ssat_max.is_valid() && ssat_max.exists != SsatSimpleState::SSAT_MAX_DOESNT_EXIST) {
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS_copy(new CoolProp::HelmholtzEOSMixtureBackend(get_components()));
        Residual resid(*HEOS_copy);
        const CoolProp::SimpleState& tripleV = HEOS_copy->get_components()[0].triple_vapor;
        double v1 = resid.call(hsat_max.T);
        double v2 = resid.call(tripleV.T);
        // If there is a sign change, there is a maxima, otherwise there is no local maxima/minima
        if (v1 * v2 < 0) {
            Brent(resid, hsat_max.T, tripleV.T, DBL_EPSILON, 1e-8, 30);
            ssat_max.T = resid.HEOS->T();
            ssat_max.p = resid.HEOS->p();
            ssat_max.rhomolar = resid.HEOS->rhomolar();
            ssat_max.hmolar = resid.HEOS->hmolar();
            ssat_max.smolar = resid.HEOS->smolar();
            ssat_max.exists = SsatSimpleState::SSAT_MAX_DOES_EXIST;
        } else {
            ssat_max.exists = SsatSimpleState::SSAT_MAX_DOESNT_EXIST;
        }
    }
}
void HelmholtzEOSMixtureBackend::calc_hsat_max(void) {
    class Residualhmax : public FuncWrapper1D
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        Residualhmax(HelmholtzEOSMixtureBackend& HEOS) : HEOS(&HEOS){};
        double call(double T) {
            HEOS->update(QT_INPUTS, 1, T);
            // dTdp_along_sat
            double dTdp_along_sat =
              HEOS->T() * (1 / HEOS->SatV->rhomolar() - 1 / HEOS->SatL->rhomolar()) / (HEOS->SatV->hmolar() - HEOS->SatL->hmolar());
            // dhdT_along_sat;
            return HEOS->SatV->first_partial_deriv(iHmolar, iT, iP) + HEOS->SatV->first_partial_deriv(iHmolar, iP, iT) / dTdp_along_sat;
        }
    };
    if (!hsat_max.is_valid()) {
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS_copy(new CoolProp::HelmholtzEOSMixtureBackend(get_components()));
        Residualhmax residhmax(*HEOS_copy);
        Brent(residhmax, T_critical() - 0.1, HEOS_copy->Ttriple() + 1, DBL_EPSILON, 1e-8, 30);
        hsat_max.T = residhmax.HEOS->T();
        hsat_max.p = residhmax.HEOS->p();
        hsat_max.rhomolar = residhmax.HEOS->rhomolar();
        hsat_max.hmolar = residhmax.HEOS->hmolar();
        hsat_max.smolar = residhmax.HEOS->smolar();
    }
}
void HelmholtzEOSMixtureBackend::T_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value) {
    if (!ValidNumber(value)) {
        throw ValueError(format("value to T_phase_determination_pure_or_pseudopure is invalid"));
    };

    // T is known, another input P, T, H, S, U is given (all molar)
    if (_T < _crit.T && _p > _crit.p) {
        // Only ever true if (other = iP); otherwise _p = -HUGE
        _phase = iphase_supercritical_liquid;
    } else if (std::abs(_T - _crit.T) < 10 * DBL_EPSILON)  // Exactly at Tcrit
    {
        switch (other) {
            case iDmolar:
                if (std::abs(_rhomolar - _crit.rhomolar) < 10 * DBL_EPSILON) {
                    _phase = iphase_critical_point;
                    break;
                } else if (_rhomolar > _crit.rhomolar) {
                    _phase = iphase_supercritical_liquid;
                    break;
                } else {
                    _phase = iphase_supercritical_gas;
                    break;
                }
            case iP: {
                if (std::abs(_p - _crit.p) < 10 * DBL_EPSILON) {
                    _phase = iphase_critical_point;
                    break;
                } else if (_p > _crit.p) {
                    _phase = iphase_supercritical_liquid;
                    break;
                } else {
                    _phase = iphase_supercritical_gas;
                    break;
                }
            }
            default:
                throw ValueError(format("T=Tcrit; invalid input for other to T_phase_determination_pure_or_pseudopure"));
        }
    } else if (_T < _crit.T)  // Gas, 2-Phase, Liquid, or Supercritical Liquid Region
    {
        // Start to think about the saturation stuff
        // First try to use the ancillary equations if you are far enough away
        // You know how accurate the ancillary equations are thanks to using CoolProp code to refit them
        switch (other) {
            case iP: {
                _pLanc = components[0].ancillaries.pL.evaluate(_T);
                _pVanc = components[0].ancillaries.pV.evaluate(_T);
                CoolPropDbl p_vap = 0.98 * static_cast<double>(_pVanc);
                CoolPropDbl p_liq = 1.02 * static_cast<double>(_pLanc);

                if (value < p_vap) {
                    this->_phase = iphase_gas;
                    _Q = -1000;
                    return;
                } else if (value > p_liq) {
                    this->_phase = iphase_liquid;
                    _Q = 1000;
                    return;
                } else if (!is_pure())  // pseudo-pure
                {
                    // For pseudo-pure fluids, the ancillary pressure curves are the official
                    // arbiter of the phase
                    if (value > static_cast<CoolPropDbl>(_pLanc)) {
                        this->_phase = iphase_liquid;
                        _Q = 1000;
                        return;
                    } else if (value < static_cast<CoolPropDbl>(_pVanc)) {
                        this->_phase = iphase_gas;
                        _Q = -1000;
                        return;
                    } else {
                        throw ValueError("Two-phase inputs not supported for pseudo-pure for now");
                    }
                }
                break;
            }
            default: {
                // Always calculate the densities using the ancillaries
                _rhoVanc = components[0].ancillaries.rhoV.evaluate(_T);
                _rhoLanc = components[0].ancillaries.rhoL.evaluate(_T);
                CoolPropDbl rho_vap = 0.95 * static_cast<double>(_rhoVanc);
                CoolPropDbl rho_liq = 1.05 * static_cast<double>(_rhoLanc);
                switch (other) {
                    case iDmolar: {
                        if (value < rho_vap) {
                            this->_phase = iphase_gas;
                            return;
                        } else if (value > rho_liq) {
                            this->_phase = iphase_liquid;
                            return;
                        } else {
                            // Next we check the vapor quality based on the ancillary values
                            double Qanc = (1 / value - 1 / static_cast<double>(_rhoLanc))
                                          / (1 / static_cast<double>(_rhoVanc) - 1 / static_cast<double>(_rhoLanc));
                            // If the vapor quality is significantly inside the two-phase zone, stop, we are definitely two-phase
                            if (value > 0.95 * rho_liq || value < 1.05 * rho_vap) {
                                // Definitely single-phase
                                _phase = iphase_liquid;  // Needed for direct update call
                                _Q = -1000;              // Needed for direct update call
                                update_DmolarT_direct(value, _T);
                                CoolPropDbl pL = components[0].ancillaries.pL.evaluate(_T);
                                if (Qanc < 0.01 && _p > pL * 1.05 && first_partial_deriv(iP, iDmolar, iT) > 0
                                    && second_partial_deriv(iP, iDmolar, iT, iDmolar, iT) > 0) {
                                    _phase = iphase_liquid;
                                    _Q = -1000;
                                    return;
                                } else if (Qanc > 1.01) {
                                    break;
                                } else {
                                    _phase = iphase_unknown;
                                    _p = _HUGE;
                                }
                            }
                        }
                        break;
                    }
                    default: {
                        if (!this->SatL || !this->SatV) {
                            throw ValueError(format("The saturation properties are needed in T_phase_determination_pure_or_pseudopure"));
                        }
                        // If it is not density, update the states
                        SatV->update(DmolarT_INPUTS, rho_vap, _T);
                        SatL->update(DmolarT_INPUTS, rho_liq, _T);

                        // First we check ancillaries
                        switch (other) {
                            case iSmolar: {
                                if (value > SatV->calc_smolar()) {
                                    this->_phase = iphase_gas;
                                    return;
                                }
                                if (value < SatL->calc_smolar()) {
                                    this->_phase = iphase_liquid;
                                    return;
                                }
                                break;
                            }
                            case iHmolar: {
                                if (value > SatV->calc_hmolar()) {
                                    this->_phase = iphase_gas;
                                    return;
                                } else if (value < SatL->calc_hmolar()) {
                                    this->_phase = iphase_liquid;
                                    return;
                                }
                                break;
                            }
                            case iUmolar: {
                                if (value > SatV->calc_umolar()) {
                                    this->_phase = iphase_gas;
                                    return;
                                } else if (value < SatL->calc_umolar()) {
                                    this->_phase = iphase_liquid;
                                    return;
                                }
                                break;
                            }
                            default:
                                throw ValueError(format("invalid input for other to T_phase_determination_pure_or_pseudopure"));
                        }
                    }
                }
            }
        }

        // Actually have to use saturation information sadly
        // For the given temperature, find the saturation state
        // Run the saturation routines to determine the saturation densities and pressures
        HelmholtzEOSMixtureBackend HEOS(components);
        SaturationSolvers::saturation_T_pure_options options;
        SaturationSolvers::saturation_T_pure(HEOS, _T, options);

        CoolPropDbl Q;

        if (other == iP) {
            if (value > HEOS.SatL->p() * (1e-6 + 1)) {
                this->_phase = iphase_liquid;
                _Q = -1000;
                return;
            } else if (value < HEOS.SatV->p() * (1 - 1e-6)) {
                this->_phase = iphase_gas;
                _Q = 1000;
                return;
            } else {
                throw ValueError(
                  format("Saturation pressure [%g Pa] corresponding to T [%g K] is within 1e-4 %% of given p [%Lg Pa]", HEOS.SatL->p(), _T, value));
            }
        }

        switch (other) {
            case iDmolar:
                Q = (1 / value - 1 / HEOS.SatL->rhomolar()) / (1 / HEOS.SatV->rhomolar() - 1 / HEOS.SatL->rhomolar());
                break;
            case iSmolar:
                Q = (value - HEOS.SatL->smolar()) / (HEOS.SatV->smolar() - HEOS.SatL->smolar());
                break;
            case iHmolar:
                Q = (value - HEOS.SatL->hmolar()) / (HEOS.SatV->hmolar() - HEOS.SatL->hmolar());
                break;
            case iUmolar:
                Q = (value - HEOS.SatL->umolar()) / (HEOS.SatV->umolar() - HEOS.SatL->umolar());
                break;
            default:
                throw ValueError(format("bad input for other"));
        }

        // Update the states
        if (this->SatL) this->SatL->update(DmolarT_INPUTS, HEOS.SatL->rhomolar(), HEOS.SatL->T());
        if (this->SatV) this->SatV->update(DmolarT_INPUTS, HEOS.SatV->rhomolar(), HEOS.SatV->T());
        // Update the two-Phase variables
        _rhoLmolar = HEOS.SatL->rhomolar();
        _rhoVmolar = HEOS.SatV->rhomolar();

        if (Q < 0) {
            this->_phase = iphase_liquid;
            _Q = -1;
            return;
        } else if (Q > 1) {
            this->_phase = iphase_gas;
            _Q = 1;
            return;
        } else {
            this->_phase = iphase_twophase;
        }
        _Q = Q;
        // Load the outputs
        _p = _Q * HEOS.SatV->p() + (1 - _Q) * HEOS.SatL->p();
        _rhomolar = 1 / (_Q / HEOS.SatV->rhomolar() + (1 - _Q) / HEOS.SatL->rhomolar());
        return;
    } else if (_T > _crit.T && _T > components[0].EOS().Ttriple)  // Supercritical or Supercritical Gas Region
    {
        _Q = 1e9;
        switch (other) {
            case iP: {
                if (_p > _crit.p) {
                    this->_phase = iphase_supercritical;
                    return;
                } else {
                    this->_phase = iphase_supercritical_gas;
                    return;
                }
            }
            case iDmolar: {
                if (_rhomolar > _crit.rhomolar) {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                } else {
                    this->_phase = iphase_supercritical_gas;
                    return;
                }
            }
            case iSmolar: {
                if (_smolar.pt() > _crit.smolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iHmolar: {
                if (_hmolar.pt() > _crit.hmolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            case iUmolar: {
                if (_umolar.pt() > _crit.umolar) {
                    this->_phase = iphase_supercritical_gas;
                    return;
                } else {
                    this->_phase = iphase_supercritical_liquid;
                    return;
                }
            }
            default: {
                throw ValueError("supercritical temp but other invalid for now");
            }
        }
    } else {
        throw ValueError(format("For now, we don't support T [%g K] below Ttriple [%g K]", _T, components[0].EOS().Ttriple));
    }
}
void get_dT_drho(HelmholtzEOSMixtureBackend* HEOS, parameters index, CoolPropDbl& dT, CoolPropDbl& drho) {
    CoolPropDbl T = HEOS->T(), rho = HEOS->rhomolar(), rhor = HEOS->get_reducing_state().rhomolar, Tr = HEOS->get_reducing_state().T,
                dT_dtau = -pow(T, 2) / Tr, R = HEOS->gas_constant(), delta = rho / rhor, tau = Tr / T;

    switch (index) {
        case iT:
            dT = 1;
            drho = 0;
            break;
        case iDmolar:
            dT = 0;
            drho = 1;
            break;
        case iDmass:
            dT = 0;
            drho = HEOS->molar_mass();
            break;
        case iP: {
            // dp/drho|T
            drho = R * T * (1 + 2 * delta * HEOS->dalphar_dDelta() + pow(delta, 2) * HEOS->d2alphar_dDelta2());
            // dp/dT|rho
            dT = rho * R * (1 + delta * HEOS->dalphar_dDelta() - tau * delta * HEOS->d2alphar_dDelta_dTau());
            break;
        }
        case iHmass:
        case iHmolar: {
            // dh/dT|rho
            dT = R
                 * (-pow(tau, 2) * (HEOS->d2alpha0_dTau2() + HEOS->d2alphar_dTau2())
                    + (1 + delta * HEOS->dalphar_dDelta() - tau * delta * HEOS->d2alphar_dDelta_dTau()));
            // dh/drhomolar|T
            drho =
              T * R / rho * (tau * delta * HEOS->d2alphar_dDelta_dTau() + delta * HEOS->dalphar_dDelta() + pow(delta, 2) * HEOS->d2alphar_dDelta2());
            if (index == iHmass) {
                // dhmolar/drhomolar|T * dhmass/dhmolar where dhmass/dhmolar = 1/mole mass
                drho /= HEOS->molar_mass();
                dT /= HEOS->molar_mass();
            }
            break;
        }
        case iSmass:
        case iSmolar: {
            // ds/dT|rho
            dT = R / T * (-pow(tau, 2) * (HEOS->d2alpha0_dTau2() + HEOS->d2alphar_dTau2()));
            // ds/drho|T
            drho = R / rho * (-(1 + delta * HEOS->dalphar_dDelta() - tau * delta * HEOS->d2alphar_dDelta_dTau()));
            if (index == iSmass) {
                // ds/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
                drho /= HEOS->molar_mass();
                dT /= HEOS->molar_mass();
            }
            break;
        }
        case iUmass:
        case iUmolar: {
            // du/dT|rho
            dT = R * (-pow(tau, 2) * (HEOS->d2alpha0_dTau2() + HEOS->d2alphar_dTau2()));
            // du/drho|T
            drho = HEOS->T() * R / rho * (tau * delta * HEOS->d2alphar_dDelta_dTau());
            if (index == iUmass) {
                // du/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
                drho /= HEOS->molar_mass();
                dT /= HEOS->molar_mass();
            }
            break;
        }
        case iTau:
            dT = 1 / dT_dtau;
            drho = 0;
            break;
        case iDelta:
            dT = 0;
            drho = 1 / rhor;
            break;
        default:
            throw ValueError(format("input to get_dT_drho[%s] is invalid", get_parameter_information(index, "short").c_str()));
    }
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_pressure_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    SimpleState reducing = calc_reducing_state_nocache(mole_fractions);
    CoolPropDbl delta = rhomolar / reducing.rhomolar;
    CoolPropDbl tau = reducing.T / T;

    // Calculate derivative
    int nTau = 0, nDelta = 1;
    CoolPropDbl dalphar_dDelta = calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, tau, delta);

    // Get pressure
    return rhomolar * gas_constant() * T * (1 + delta * dalphar_dDelta);
}
HelmholtzEOSBackend::StationaryPointReturnFlag HelmholtzEOSMixtureBackend::solver_dpdrho0_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomax,
                                                                                             CoolPropDbl& light, CoolPropDbl& heavy) {

    /// The residual to be used to find the location where dpdrho=0 for given T
    class dpdrho_resid : public FuncWrapper1DWithTwoDerivs
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl T, p, delta, rhor, tau, R_u;

        dpdrho_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl T, CoolPropDbl p)
          : HEOS(HEOS),
            T(T),
            p(p),
            delta(_HUGE),
            rhor(HEOS->get_reducing_state().rhomolar),
            tau(HEOS->get_reducing_state().T / T),
            R_u(HEOS->gas_constant()) {}
        double call(double rhomolar) {
            delta = rhomolar / rhor;  // needed for derivative
            HEOS->update_DmolarT_direct(rhomolar, T);
            // dp/drho|T
            return R_u * T * (1 + 2 * delta * HEOS->dalphar_dDelta() + POW2(delta) * HEOS->d2alphar_dDelta2());
        };
        double deriv(double rhomolar) {
            // d2p/drho2|T
            return R_u * T / rhor * (2 * HEOS->dalphar_dDelta() + 4 * delta * HEOS->d2alphar_dDelta2() + POW2(delta) * HEOS->calc_d3alphar_dDelta3());
        };
        double second_deriv(double rhomolar) {
            // d3p/drho3|T
            return R_u * T / POW2(rhor)
                   * (6 * HEOS->d2alphar_dDelta2() + 6 * delta * HEOS->d3alphar_dDelta3() + POW2(delta) * HEOS->calc_d4alphar_dDelta4());
        };
    };
    dpdrho_resid resid(this, T, p);
    light = -1;
    heavy = -1;
    try {
        light = Halley(resid, 1e-6, 1e-8, 100);
        double d2pdrho2__constT = resid.deriv(light);
        if (d2pdrho2__constT > 0) {
            // Not possible since curvature should be negative
            throw CoolProp::ValueError("curvature cannot be positive");
        }
    } catch (std::exception& e) {
        if (get_debug_level() > 5) {
            std::cout << e.what() << std::endl;
        };
        light = -1;
    }

    if (light < 0) {
        try {
            // Now we are going to do something VERY slow - increase density until curvature is positive
            double rho = 1e-6;
            for (std::size_t counter = 0; counter <= 100; counter++) {
                resid.call(rho);  // Updates the state
                double curvature = resid.deriv(rho);
                if (curvature > 0) {
                    light = rho;
                    break;
                }
                rho *= 2;
            }
        } catch (...) {
        }
    }

    // First try a "normal" calculation of the stationary point on the liquid side
    for (double omega = 0.7; omega > 0; omega -= 0.2) {
        try {
            resid.options.add_number("omega", omega);
            heavy = Halley(resid, rhomax, 1e-8, 100);
            double d2pdrho2__constT = resid.deriv(heavy);
            if (d2pdrho2__constT < 0) {
                // Not possible since curvature should be positive
                throw CoolProp::ValueError("curvature cannot be negative");
            }
            break;  // Jump out, we got a good solution
        } catch (std::exception& e) {
            if (get_debug_level() > 5) {
                std::cout << e.what() << std::endl;
            };
            heavy = -1;
        }
    }

    if (heavy < 0) {
        try {
            // Now we are going to do something VERY slow - decrease density until curvature is negative or pressure is negative
            double rho = rhomax;
            for (std::size_t counter = 0; counter <= 100; counter++) {
                resid.call(rho);  // Updates the state
                double curvature = resid.deriv(rho);
                if (curvature < 0 || this->p() < 0) {
                    heavy = rho;
                    break;
                }
                rho /= 1.1;
            }
        } catch (...) {
        }
    }

    if (light > 0 && heavy > 0) {
        // Found two stationary points, done!
        return TWO_STATIONARY_POINTS_FOUND;
    }
    // If no solution is found for dpdrho|T=0 starting at high and low densities,
    // then try to do a bounded solver to see if you can find any solutions.  If you
    // can't, p = f(rho) is probably monotonic (supercritical?), and the bounds are
    else if (light < 0 && heavy < 0) {
        double dpdrho_min = resid.call(1e-10);
        double dpdrho_max = resid.call(rhomax);
        if (dpdrho_max * dpdrho_min > 0) {
            return ZERO_STATIONARY_POINTS;
        } else {
            throw CoolProp::ValueError("zero stationary points -- does this make sense?");
        }
    } else {
        return ONE_STATIONARY_POINT_FOUND;
    }
}
// Define the residual to be driven to zero
class SolverTPResid : public FuncWrapper1DWithThreeDerivs
{
   public:
    HelmholtzEOSMixtureBackend* HEOS;
    CoolPropDbl T, p, delta, rhor, tau, R_u;

    SolverTPResid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl T, CoolPropDbl p)
      : HEOS(HEOS),
        T(T),
        p(p),
        delta(_HUGE),
        rhor(HEOS->get_reducing_state().rhomolar),
        tau(HEOS->get_reducing_state().T / T),
        R_u(HEOS->gas_constant()) {}
    double call(double rhomolar) {
        delta = rhomolar / rhor;  // needed for derivative
        HEOS->update_DmolarT_direct(rhomolar, T);
        CoolPropDbl peos = HEOS->p();
        return (peos - p) / p;
    };
    double deriv(double rhomolar) {
        // dp/drho|T / pspecified
        return R_u * T * (1 + 2 * delta * HEOS->dalphar_dDelta() + POW2(delta) * HEOS->d2alphar_dDelta2()) / p;
    };
    double second_deriv(double rhomolar) {
        // d2p/drho2|T / pspecified
        return R_u * T / rhor * (2 * HEOS->dalphar_dDelta() + 4 * delta * HEOS->d2alphar_dDelta2() + POW2(delta) * HEOS->calc_d3alphar_dDelta3()) / p;
    };
    double third_deriv(double rhomolar) {
        // d3p/drho3|T / pspecified
        return R_u * T / POW2(rhor)
               * (6 * HEOS->d2alphar_dDelta2() + 6 * delta * HEOS->d3alphar_dDelta3() + POW2(delta) * HEOS->calc_d4alphar_dDelta4()) / p;
    };
};
CoolPropDbl HelmholtzEOSMixtureBackend::SRK_covolume() {
    double b = 0;
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
        // Get the parameters for the cubic EOS
        CoolPropDbl Tc = get_fluid_constant(i, iT_critical), pc = get_fluid_constant(i, iP_critical);
        CoolPropDbl R = 8.3144598;
        b += mole_fractions[i] * 0.08664 * R * Tc / pc;
    }
    return b;
}
CoolPropDbl HelmholtzEOSMixtureBackend::solver_rho_Tp_global(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_max) {
    // Find the densities along the isotherm where dpdrho|T = 0 (if you can)
    CoolPropDbl light = -1, heavy = -1;
    StationaryPointReturnFlag retval = solver_dpdrho0_Tp(T, p, rhomolar_max, light, heavy);

    // Define the solver class
    SolverTPResid resid(this, T, p);

    if (retval == ZERO_STATIONARY_POINTS) {
        // It's monotonic (no stationary points found), so do the full bounded solver
        // for the density
        double rho = Brent(resid, 1e-10, rhomolar_max, DBL_EPSILON, 1e-8, 100);
        return rho;
    } else if (retval == TWO_STATIONARY_POINTS_FOUND) {

        // Calculate the pressures at the min and max densities where dpdrho|T = 0
        double p_at_rhomin_stationary = calc_pressure_nocache(T, light);
        double p_at_rhomax_stationary = calc_pressure_nocache(T, heavy);

        double rho_liq = -1, rho_vap = -1;
        if (p > p_at_rhomax_stationary) {
            int counter = 0;
            for (/* init above, for debugging */; counter <= 10; counter++) {
                // Bump up rhomax if needed to bound the given pressure
                double p_at_rhomax = calc_pressure_nocache(T, rhomolar_max);
                if (p_at_rhomax < p) {
                    rhomolar_max *= 1.05;
                } else {
                    break;
                }
            }
            // Look for liquid root starting at stationary point density
            rho_liq = Brent(resid, heavy, rhomolar_max, DBL_EPSILON, 1e-8, 100);
        }

        if (p < p_at_rhomin_stationary) {
            // Look for vapor root starting at stationary point density
            rho_vap = Brent(resid, light, 1e-10, DBL_EPSILON, 1e-8, 100);
        }

        if (rho_vap > 0 && rho_liq > 0) {
            // Both densities are the same
            if (std::abs(rho_vap - rho_liq) < 1e-10) {
                // return one of them
                return rho_vap;
            } else {
                // Two solutions found, keep the one with lower Gibbs energy
                double gibbsmolar_vap = calc_gibbsmolar_nocache(T, rho_vap);
                double gibbsmolar_liq = calc_gibbsmolar_nocache(T, rho_liq);
                if (gibbsmolar_liq < gibbsmolar_vap) {
                    return rho_liq;
                } else {
                    return rho_vap;
                }
            }
        } else if (rho_vap < 0 && rho_liq > 0) {
            // Liquid root found, return it
            return rho_liq;
        } else if (rho_vap > 0 && rho_liq < 0) {
            // Vapor root found, return it
            return rho_vap;
        } else {
            throw CoolProp::ValueError(format("No density solutions for T=%g,p=%g,z=%s", T, p, vec_to_string(mole_fractions, "%0.12g").c_str()));
        }
    } else {
        throw CoolProp::ValueError(
          format("One stationary point (not good) for T=%g,p=%g,z=%s", T, p, vec_to_string(mole_fractions, "%0.12g").c_str()));
    }
};

CoolPropDbl HelmholtzEOSMixtureBackend::solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_guess) {
    phases phase;

    SolverTPResid resid(this, T, p);

    // Check if the phase is imposed
    if (imposed_phase_index != iphase_not_imposed)
        // Use the imposed phase index
        phase = imposed_phase_index;
    else
        // Use the phase index in the class
        phase = _phase;

    if (rhomolar_guess < 0)  // Not provided
    {
        // Calculate a guess value using SRK equation of state
        rhomolar_guess = solver_rho_Tp_SRK(T, p, phase);

        // A gas-like phase, ideal gas might not be the perfect model, but probably good enough
        if (phase == iphase_gas || phase == iphase_supercritical_gas || phase == iphase_supercritical) {
            if (rhomolar_guess < 0 || !ValidNumber(rhomolar_guess))  // If the guess is bad, probably high temperature, use ideal gas
            {
                rhomolar_guess = p / (gas_constant() * T);
            }
        } else if (phase == iphase_liquid) {
            double rhomolar;
            if (is_pure_or_pseudopure) {
                // It's liquid at subcritical pressure, we can use ancillaries as guess value
                CoolPropDbl _rhoLancval = static_cast<CoolPropDbl>(components[0].ancillaries.rhoL.evaluate(T));
                try {
                    // First we try with Halley's method starting at saturated liquid
                    rhomolar = Halley(resid, _rhoLancval, 1e-8, 100);
                    if (!ValidNumber(rhomolar) || first_partial_deriv(iP, iDmolar, iT) < 0
                        || second_partial_deriv(iP, iDmolar, iT, iDmolar, iT) < 0) {
                        throw ValueError("Liquid density is invalid");
                    }
                } catch (std::exception&) {
                    // Next we try with a Brent method bounded solver since the function is 1-1
                    rhomolar = Brent(resid, _rhoLancval * 0.9, _rhoLancval * 1.3, DBL_EPSILON, 1e-8, 100);
                    if (!ValidNumber(rhomolar)) {
                        throw ValueError();
                    }
                }
            } else {
                // Try with 4th order Householder method starting at a very high density
                rhomolar = Householder4(&resid, 3 * rhomolar_reducing(), 1e-8, 100);
            }
            return rhomolar;
        } else if (phase == iphase_supercritical_liquid) {
            CoolPropDbl rhoLancval = static_cast<CoolPropDbl>(components[0].ancillaries.rhoL.evaluate(T));
            // Next we try with a Brent method bounded solver since the function is 1-1
            double rhomolar = Brent(resid, rhoLancval * 0.99, rhomolar_critical() * 4, DBL_EPSILON, 1e-8, 100);
            if (!ValidNumber(rhomolar)) {
                throw ValueError();
            }
            return rhomolar;
        }
    }

    try {
        double rhomolar = Householder4(resid, rhomolar_guess, 1e-8, 20);
        if (!ValidNumber(rhomolar) || rhomolar < 0) {
            throw ValueError();
        }
        if (phase == iphase_liquid) {
            double dpdrho = first_partial_deriv(iP, iDmolar, iT);
            double d2pdrho2 = second_partial_deriv(iP, iDmolar, iT, iDmolar, iT);
            if (dpdrho < 0 || d2pdrho2 < 0) {
                // Try again with a larger density in order to end up at the right solution
                rhomolar = Householder4(resid, 3 * rhomolar_reducing(), 1e-8, 100);
                return rhomolar;
            }
        } else if (phase == iphase_gas) {
            double dpdrho = first_partial_deriv(iP, iDmolar, iT);
            double d2pdrho2 = second_partial_deriv(iP, iDmolar, iT, iDmolar, iT);
            if (dpdrho < 0 || d2pdrho2 > 0) {
                // Try again with a tiny density in order to end up at the right solution
                rhomolar = Householder4(resid, 1e-6, 1e-8, 100);
                return rhomolar;
            }
        }
        return rhomolar;
    } catch (std::exception& e) {
        if (phase == iphase_supercritical || phase == iphase_supercritical_gas) {
            double rhomolar = Brent(resid, 1e-10, 3 * rhomolar_reducing(), DBL_EPSILON, 1e-8, 100);
            return rhomolar;
        } else if (is_pure_or_pseudopure && T > T_critical()) {
            try {
                double rhomolar = Brent(resid, 1e-10, 5 * rhomolar_reducing(), DBL_EPSILON, 1e-8, 100);
                return rhomolar;

            } catch (...) {
                double rhomolar = Householder4(resid, 3 * rhomolar_reducing(), 1e-8, 100);
                return rhomolar;
            }
        }
        throw ValueError(format("solver_rho_Tp was unable to find a solution for T=%10Lg, p=%10Lg, with guess value %10Lg with error: %s", T, p,
                                rhomolar_guess, e.what()));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::solver_rho_Tp_SRK(CoolPropDbl T, CoolPropDbl p, phases phase) {
    CoolPropDbl rhomolar, R_u = gas_constant(), a = 0, b = 0, k_ij = 0;

    for (std::size_t i = 0; i < components.size(); ++i) {
        CoolPropDbl Tci = components[i].EOS().reduce.T, pci = components[i].EOS().reduce.p, acentric_i = components[i].EOS().acentric;
        CoolPropDbl m_i = 0.480 + 1.574 * acentric_i - 0.176 * pow(acentric_i, 2);
        CoolPropDbl b_i = 0.08664 * R_u * Tci / pci;
        b += mole_fractions[i] * b_i;

        CoolPropDbl a_i = 0.42747 * pow(R_u * Tci, 2) / pci * pow(1 + m_i * (1 - sqrt(T / Tci)), 2);

        for (std::size_t j = 0; j < components.size(); ++j) {
            CoolPropDbl Tcj = components[j].EOS().reduce.T, pcj = components[j].EOS().reduce.p, acentric_j = components[j].EOS().acentric;
            CoolPropDbl m_j = 0.480 + 1.574 * acentric_j - 0.176 * pow(acentric_j, 2);

            CoolPropDbl a_j = 0.42747 * pow(R_u * Tcj, 2) / pcj * pow(1 + m_j * (1 - sqrt(T / Tcj)), 2);

            k_ij = 0;
            //if (i == j){
            //    k_ij = 0;
            //}
            //else{
            //    k_ij = 0;
            //}

            a += mole_fractions[i] * mole_fractions[j] * sqrt(a_i * a_j) * (1 - k_ij);
        }
    }

    CoolPropDbl A = a * p / pow(R_u * T, 2);
    CoolPropDbl B = b * p / (R_u * T);

    //Solve the cubic for solutions for Z = p/(rho*R*T)
    double Z0, Z1, Z2;
    int Nsolns;
    solve_cubic(1, -1, A - B - B * B, -A * B, Nsolns, Z0, Z1, Z2);

    // Determine the guess value
    if (Nsolns == 1) {
        rhomolar = p / (Z0 * R_u * T);
    } else {
        CoolPropDbl rhomolar0 = p / (Z0 * R_u * T);
        CoolPropDbl rhomolar1 = p / (Z1 * R_u * T);
        CoolPropDbl rhomolar2 = p / (Z2 * R_u * T);

        // Check if only one solution is positive, return the solution if that is the case
        if (rhomolar0 > 0 && rhomolar1 <= 0 && rhomolar2 <= 0) {
            return rhomolar0;
        }
        if (rhomolar0 <= 0 && rhomolar1 > 0 && rhomolar2 <= 0) {
            return rhomolar1;
        }
        if (rhomolar0 <= 0 && rhomolar1 <= 0 && rhomolar2 > 0) {
            return rhomolar2;
        }

        switch (phase) {
            case iphase_liquid:
            case iphase_supercritical_liquid:
                rhomolar = max3(rhomolar0, rhomolar1, rhomolar2);
                break;
            case iphase_gas:
            case iphase_supercritical_gas:
            case iphase_supercritical:
                rhomolar = min3(rhomolar0, rhomolar1, rhomolar2);
                break;
            default:
                throw ValueError("Bad phase to solver_rho_Tp_SRK");
        };
    }
    return rhomolar;
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_pressure(void) {
    // Calculate the reducing parameters
    _delta = _rhomolar / _reducing.rhomolar;
    _tau = _reducing.T / _T;

    // Calculate derivative if needed
    CoolPropDbl dar_dDelta = dalphar_dDelta();
    CoolPropDbl R_u = gas_constant();

    // Get pressure
    _p = _rhomolar * R_u * _T * (1 + _delta.pt() * dar_dDelta);

    //std::cout << format("p: %13.12f %13.12f %10.9f %10.9f %10.9f %10.9f %g\n",_T,_rhomolar,_tau,_delta,mole_fractions[0],dar_dDelta,_p);
    //if (_p < 0){
    //    throw ValueError("Pressure is less than zero");
    //}

    return static_cast<CoolPropDbl>(_p);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_hmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    // Calculate the reducing parameters
    CoolPropDbl delta = rhomolar / _reducing.rhomolar;
    CoolPropDbl tau = _reducing.T / T;

    // Calculate derivatives if needed, or just use cached values
    // Calculate derivative if needed
    CoolPropDbl dar_dDelta = calc_alphar_deriv_nocache(0, 1, mole_fractions, tau, delta);
    CoolPropDbl dar_dTau = calc_alphar_deriv_nocache(1, 0, mole_fractions, tau, delta);
    CoolPropDbl da0_dTau = calc_alpha0_deriv_nocache(1, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    CoolPropDbl R_u = gas_constant();

    // Get molar enthalpy
    return R_u * T * (1 + tau * (da0_dTau + dar_dTau) + delta * dar_dDelta);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_hmolar(void) {
    if (get_debug_level() >= 50)
        std::cout << format("HelmholtzEOSMixtureBackend::calc_hmolar: 2phase: %d T: %g rhomomolar: %g", isTwoPhase(), _T, _rhomolar) << std::endl;
    if (isTwoPhase()) {
        if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
        if (std::abs(_Q) < DBL_EPSILON) {
            _hmolar = SatL->hmolar();
        } else if (std::abs(_Q - 1) < DBL_EPSILON) {
            _hmolar = SatV->hmolar();
        } else {
            _hmolar = _Q * SatV->hmolar() + (1 - _Q) * SatL->hmolar();
        }
        return static_cast<CoolPropDbl>(_hmolar);
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl da0_dTau = dalpha0_dTau();
        CoolPropDbl dar_dTau = dalphar_dTau();
        CoolPropDbl dar_dDelta = dalphar_dDelta();
        CoolPropDbl R_u = gas_constant();

        // Get molar enthalpy
        _hmolar = R_u * _T * (1 + _tau.pt() * (da0_dTau + dar_dTau) + _delta.pt() * dar_dDelta);

        return static_cast<CoolPropDbl>(_hmolar);
    } else {
        throw ValueError(format("phase is invalid in calc_hmolar"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_smolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    // Calculate the reducing parameters
    CoolPropDbl delta = rhomolar / _reducing.rhomolar;
    CoolPropDbl tau = _reducing.T / T;

    // Calculate derivatives if needed, or just use cached values
    // Calculate derivative if needed
    CoolPropDbl dar_dTau = calc_alphar_deriv_nocache(1, 0, mole_fractions, tau, delta);
    CoolPropDbl ar = calc_alphar_deriv_nocache(0, 0, mole_fractions, tau, delta);
    CoolPropDbl da0_dTau = calc_alpha0_deriv_nocache(1, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    CoolPropDbl a0 = calc_alpha0_deriv_nocache(0, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    CoolPropDbl R_u = gas_constant();

    // Get molar entropy
    return R_u * (tau * (da0_dTau + dar_dTau) - a0 - ar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_smolar(void) {
    if (isTwoPhase()) {
        if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
        if (std::abs(_Q) < DBL_EPSILON) {
            _smolar = SatL->smolar();
        } else if (std::abs(_Q - 1) < DBL_EPSILON) {
            _smolar = SatV->smolar();
        } else {
            _smolar = _Q * SatV->smolar() + (1 - _Q) * SatL->smolar();
        }
        return static_cast<CoolPropDbl>(_smolar);
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl da0_dTau = dalpha0_dTau();
        CoolPropDbl ar = alphar();
        CoolPropDbl a0 = alpha0();
        CoolPropDbl dar_dTau = dalphar_dTau();
        CoolPropDbl R_u = gas_constant();

        // Get molar entropy
        _smolar = R_u * (_tau.pt() * (da0_dTau + dar_dTau) - a0 - ar);

        return static_cast<CoolPropDbl>(_smolar);
    } else {
        throw ValueError(format("phase is invalid in calc_smolar"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_umolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    // Calculate the reducing parameters
    CoolPropDbl delta = rhomolar / _reducing.rhomolar;
    CoolPropDbl tau = _reducing.T / T;

    // Calculate derivatives
    CoolPropDbl dar_dTau = calc_alphar_deriv_nocache(1, 0, mole_fractions, tau, delta);
    CoolPropDbl da0_dTau = calc_alpha0_deriv_nocache(1, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    CoolPropDbl R_u = gas_constant();

    // Get molar internal energy
    return R_u * T * tau * (da0_dTau + dar_dTau);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_umolar(void) {
    if (isTwoPhase()) {
        if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
        if (std::abs(_Q) < DBL_EPSILON) {
            _umolar = SatL->umolar();
        } else if (std::abs(_Q - 1) < DBL_EPSILON) {
            _umolar = SatV->umolar();
        } else {
            _umolar = _Q * SatV->umolar() + (1 - _Q) * SatL->umolar();
        }
        return static_cast<CoolPropDbl>(_umolar);
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl da0_dTau = dalpha0_dTau();
        CoolPropDbl dar_dTau = dalphar_dTau();
        CoolPropDbl R_u = gas_constant();

        // Get molar internal energy
        _umolar = R_u * _T * _tau.pt() * (da0_dTau + dar_dTau);

        return static_cast<CoolPropDbl>(_umolar);
    } else {
        throw ValueError(format("phase is invalid in calc_umolar"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_cvmolar(void) {
    // Calculate the reducing parameters
    _delta = _rhomolar / _reducing.rhomolar;
    _tau = _reducing.T / _T;

    // Calculate derivatives if needed, or just use cached values
    CoolPropDbl d2ar_dTau2 = d2alphar_dTau2();
    CoolPropDbl d2a0_dTau2 = d2alpha0_dTau2();
    CoolPropDbl R_u = gas_constant();

    // Get cv
    _cvmolar = -R_u * pow(_tau.pt(), 2) * (d2ar_dTau2 + d2a0_dTau2);

    return static_cast<double>(_cvmolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_cpmolar(void) {
    // Calculate the reducing parameters
    _delta = _rhomolar / _reducing.rhomolar;
    _tau = _reducing.T / _T;

    // Calculate derivatives if needed, or just use cached values
    CoolPropDbl d2a0_dTau2 = d2alpha0_dTau2();
    CoolPropDbl dar_dDelta = dalphar_dDelta();
    CoolPropDbl d2ar_dDelta2 = d2alphar_dDelta2();
    CoolPropDbl d2ar_dDelta_dTau = d2alphar_dDelta_dTau();
    CoolPropDbl d2ar_dTau2 = d2alphar_dTau2();
    CoolPropDbl R_u = gas_constant();

    // Get cp
    _cpmolar = R_u
               * (-pow(_tau.pt(), 2) * (d2ar_dTau2 + d2a0_dTau2)
                  + pow(1 + _delta.pt() * dar_dDelta - _delta.pt() * _tau.pt() * d2ar_dDelta_dTau, 2)
                      / (1 + 2 * _delta.pt() * dar_dDelta + pow(_delta.pt(), 2) * d2ar_dDelta2));

    return static_cast<double>(_cpmolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_cpmolar_idealgas(void) {
    // Calculate the reducing parameters
    _delta = _rhomolar / _reducing.rhomolar;
    _tau = _reducing.T / _T;

    // Calculate derivatives if needed, or just use cached values
    CoolPropDbl d2a0_dTau2 = d2alpha0_dTau2();
    CoolPropDbl R_u = gas_constant();

    // Get cp of the ideal gas
    return R_u * (1 + (-pow(_tau.pt(), 2)) * d2a0_dTau2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_speed_sound(void) {
    if (isTwoPhase()) {
        if (std::abs(_Q) < DBL_EPSILON) {
            return SatL->speed_sound();
        } else if (std::abs(_Q - 1) < DBL_EPSILON) {
            return SatV->speed_sound();
        } else {
            throw ValueError(format("Speed of sound is not defined for two-phase states because it depends on the distribution of phases."));
        }
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl d2a0_dTau2 = d2alpha0_dTau2();
        CoolPropDbl dar_dDelta = dalphar_dDelta();
        CoolPropDbl d2ar_dDelta2 = d2alphar_dDelta2();
        CoolPropDbl d2ar_dDelta_dTau = d2alphar_dDelta_dTau();
        CoolPropDbl d2ar_dTau2 = d2alphar_dTau2();
        CoolPropDbl R_u = gas_constant();
        CoolPropDbl mm = molar_mass();

        // Get speed of sound
        _speed_sound = sqrt(
          R_u * _T / mm
          * (1 + 2 * _delta.pt() * dar_dDelta + pow(_delta.pt(), 2) * d2ar_dDelta2
             - pow(1 + _delta.pt() * dar_dDelta - _delta.pt() * _tau.pt() * d2ar_dDelta_dTau, 2) / (pow(_tau.pt(), 2) * (d2ar_dTau2 + d2a0_dTau2))));

        return static_cast<CoolPropDbl>(_speed_sound);
    } else {
        throw ValueError(format("phase is invalid in calc_speed_sound"));
    }
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_gibbsmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    // Calculate the reducing parameters
    CoolPropDbl delta = rhomolar / _reducing.rhomolar;
    CoolPropDbl tau = _reducing.T / T;

    // Calculate derivatives
    CoolPropDbl dar_dDelta = calc_alphar_deriv_nocache(0, 1, mole_fractions, tau, delta);
    CoolPropDbl ar = calc_alphar_deriv_nocache(0, 0, mole_fractions, tau, delta);
    CoolPropDbl a0 = calc_alpha0_deriv_nocache(0, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);

    // Get molar gibbs function
    return gas_constant() * T * (1 + a0 + ar + delta * dar_dDelta);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_gibbsmolar(void) {
    if (isTwoPhase()) {
        if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
        _gibbsmolar = _Q * SatV->gibbsmolar() + (1 - _Q) * SatL->gibbsmolar();
        return static_cast<CoolPropDbl>(_gibbsmolar);
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl ar = alphar();
        CoolPropDbl a0 = alpha0();
        CoolPropDbl dar_dDelta = dalphar_dDelta();
        CoolPropDbl R_u = gas_constant();

        // Get molar gibbs function
        _gibbsmolar = R_u * _T * (1 + a0 + ar + _delta.pt() * dar_dDelta);

        return static_cast<CoolPropDbl>(_gibbsmolar);
    } else {
        throw ValueError(format("phase is invalid in calc_gibbsmolar"));
    }
}
void HelmholtzEOSMixtureBackend::calc_excess_properties(void) {
    _gibbsmolar_excess = this->gibbsmolar(), _smolar_excess = this->smolar(), _hmolar_excess = this->hmolar();
    _umolar_excess = this->umolar();
    _volumemolar_excess = 1 / this->rhomolar();
    for (std::size_t i = 0; i < components.size(); ++i) {
        transient_pure_state.reset(new HelmholtzEOSBackend(components[i].name));
        transient_pure_state->update(PT_INPUTS, p(), T());
        double x_i = mole_fractions[i];
        double R = gas_constant();
        _gibbsmolar_excess = static_cast<double>(_gibbsmolar_excess) - x_i * (transient_pure_state->gibbsmolar() + R * T() * log(x_i));
        _hmolar_excess = static_cast<double>(_hmolar_excess) - x_i * transient_pure_state->hmolar();
        _umolar_excess = static_cast<double>(_umolar_excess) - x_i * transient_pure_state->umolar();
        _smolar_excess = static_cast<double>(_smolar_excess) - x_i * (transient_pure_state->smolar() - R * log(x_i));
        _volumemolar_excess = static_cast<double>(_volumemolar_excess) - x_i / transient_pure_state->rhomolar();
    }
    _helmholtzmolar_excess = static_cast<double>(_umolar_excess) - _T * static_cast<double>(_smolar_excess);
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_helmholtzmolar(void) {
    if (isTwoPhase()) {
        if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
        _helmholtzmolar = _Q * SatV->helmholtzmolar() + (1 - _Q) * SatL->helmholtzmolar();
        return static_cast<CoolPropDbl>(_helmholtzmolar);
    } else if (isHomogeneousPhase()) {
        // Calculate the reducing parameters
        _delta = _rhomolar / _reducing.rhomolar;
        _tau = _reducing.T / _T;

        // Calculate derivatives if needed, or just use cached values
        CoolPropDbl ar = alphar();
        CoolPropDbl a0 = alpha0();
        CoolPropDbl R_u = gas_constant();

        // Get molar Helmholtz energy
        _helmholtzmolar = R_u * _T * (a0 + ar);

        return static_cast<CoolPropDbl>(_helmholtzmolar);
    } else {
        throw ValueError(format("phase is invalid in calc_helmholtzmolar"));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_fugacity_coefficient(std::size_t i) {
    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    return exp(MixtureDerivatives::ln_fugacity_coefficient(*this, i, xN_flag));
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_fugacity(std::size_t i) {
    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    return MixtureDerivatives::fugacity_i(*this, i, xN_flag);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_chemical_potential(std::size_t i) {
    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    double Tci = get_fluid_constant(i, iT_critical);
    double rhoci = get_fluid_constant(i, irhomolar_critical);
    double dnar_dni__const_T_V_nj = MixtureDerivatives::dnalphar_dni__constT_V_nj(*this, i, xN_flag);
    double dna0_dni__const_T_V_nj =
      components[i].EOS().alpha0.base(tau() * (Tci / T_reducing()), delta() / (rhoci / rhomolar_reducing())) + 1 + log(mole_fractions[i]);
    return gas_constant() * T() * (dna0_dni__const_T_V_nj + dnar_dni__const_T_V_nj);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_phase_identification_parameter(void) {
    return 2
           - rhomolar()
               * (second_partial_deriv(iP, iDmolar, iT, iT, iDmolar) / first_partial_deriv(iP, iT, iDmolar)
                  - second_partial_deriv(iP, iDmolar, iT, iDmolar, iT) / first_partial_deriv(iP, iDmolar, iT));
}

SimpleState HelmholtzEOSMixtureBackend::calc_reducing_state_nocache(const std::vector<CoolPropDbl>& mole_fractions) {
    SimpleState reducing;
    if (is_pure_or_pseudopure) {
        reducing = components[0].EOS().reduce;
    } else {
        reducing.T = Reducing->Tr(mole_fractions);
        reducing.rhomolar = Reducing->rhormolar(mole_fractions);
    }
    return reducing;
}
void HelmholtzEOSMixtureBackend::calc_reducing_state(void) {
    if (get_mole_fractions_ref().empty()) {
        throw ValueError("Mole fractions must be set before calling calc_reducing_state");
    }
    /// \todo set critical independently
    _reducing = calc_reducing_state_nocache(mole_fractions);
    _crit = _reducing;
}
void HelmholtzEOSMixtureBackend::calc_all_alphar_deriv_cache(const std::vector<CoolPropDbl>& mole_fractions, const CoolPropDbl& tau,
                                                             const CoolPropDbl& delta) {
    deriv_counter++;
    bool cache_values = true;
    HelmholtzDerivatives derivs = residual_helmholtz->all(*this, get_mole_fractions_ref(), tau, delta, cache_values);
    _alphar = derivs.alphar;
    _dalphar_dDelta = derivs.dalphar_ddelta;
    _dalphar_dTau = derivs.dalphar_dtau;
    _d2alphar_dDelta2 = derivs.d2alphar_ddelta2;
    _d2alphar_dDelta_dTau = derivs.d2alphar_ddelta_dtau;
    _d2alphar_dTau2 = derivs.d2alphar_dtau2;
    _d3alphar_dDelta3 = derivs.d3alphar_ddelta3;
    _d3alphar_dDelta2_dTau = derivs.d3alphar_ddelta2_dtau;
    _d3alphar_dDelta_dTau2 = derivs.d3alphar_ddelta_dtau2;
    _d3alphar_dTau3 = derivs.d3alphar_dtau3;
    _d4alphar_dDelta4 = derivs.d4alphar_ddelta4;
    _d4alphar_dDelta3_dTau = derivs.d4alphar_ddelta3_dtau;
    _d4alphar_dDelta2_dTau2 = derivs.d4alphar_ddelta2_dtau2;
    _d4alphar_dDelta_dTau3 = derivs.d4alphar_ddelta_dtau3;
    _d4alphar_dTau4 = derivs.d4alphar_dtau4;
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl>& mole_fractions,
                                                                  const CoolPropDbl& tau, const CoolPropDbl& delta) {
    bool cache_values = false;
    HelmholtzDerivatives derivs = residual_helmholtz->all(*this, mole_fractions, tau, delta, cache_values);
    return derivs.get(nTau, nDelta);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl>& mole_fractions,
                                                                  const CoolPropDbl& tau, const CoolPropDbl& delta, const CoolPropDbl& Tr,
                                                                  const CoolPropDbl& rhor) {
    CoolPropDbl val;
    if (components.size() == 0) {
        throw ValueError("No alpha0 derivatives are available");
    }
    if (is_pure_or_pseudopure) {
        EquationOfState& E = components[0].EOS();
        // In the case of cubics, we need to use the shifted tau^*=Tc/T and delta^*=rho/rhoc
        // rather than tau=Tr/T and delta=rho/rhor
        // For multiparameter EOS, this changes nothing because Tc/Tr = 1 and rhoc/rhor = 1
        double Tc = get_fluid_constant(0, iT_reducing), rhomolarc = get_fluid_constant(0, irhomolar_reducing);

        // Cache the reducing temperature in some terms that need it (GERG-2004 models)
        E.alpha0.set_Tred(Tc);
        double taustar = Tc / Tr * tau, deltastar = rhor / rhomolarc * delta;
        if (nTau == 0 && nDelta == 0) {
            val = E.base0(taustar, deltastar);
        } else if (nTau == 0 && nDelta == 1) {
            val = E.dalpha0_dDelta(taustar, deltastar);
        } else if (nTau == 1 && nDelta == 0) {
            val = E.dalpha0_dTau(taustar, deltastar);
        } else if (nTau == 0 && nDelta == 2) {
            val = E.d2alpha0_dDelta2(taustar, deltastar);
        } else if (nTau == 1 && nDelta == 1) {
            val = E.d2alpha0_dDelta_dTau(taustar, deltastar);
        } else if (nTau == 2 && nDelta == 0) {
            val = E.d2alpha0_dTau2(taustar, deltastar);
        } else if (nTau == 0 && nDelta == 3) {
            val = E.d3alpha0_dDelta3(taustar, deltastar);
        } else if (nTau == 1 && nDelta == 2) {
            val = E.d3alpha0_dDelta2_dTau(taustar, deltastar);
        } else if (nTau == 2 && nDelta == 1) {
            val = E.d3alpha0_dDelta_dTau2(taustar, deltastar);
        } else if (nTau == 3 && nDelta == 0) {
            val = E.d3alpha0_dTau3(taustar, deltastar);
        } else {
            throw ValueError();
        }
        val *= pow(rhor / rhomolarc, nDelta);
        val /= pow(Tr / Tc, nTau);
        if (!ValidNumber(val)) {
            //calc_alpha0_deriv_nocache(nTau,nDelta,mole_fractions,tau,delta,Tr,rhor);
            throw ValueError(format("calc_alpha0_deriv_nocache returned invalid number with inputs nTau: %d, nDelta: %d, tau: %Lg, delta: %Lg", nTau,
                                    nDelta, tau, delta));
        } else {
            return val;
        }
    } else {
        // See Table B5, GERG 2008 from Kunz Wagner, JCED, 2012
        std::size_t N = mole_fractions.size();
        CoolPropDbl summer = 0;
        CoolPropDbl tau_i, delta_i, rho_ci, T_ci;
        CoolPropDbl Rmix = gas_constant();
        for (unsigned int i = 0; i < N; ++i) {

            rho_ci = get_fluid_constant(i, irhomolar_critical);
            T_ci = get_fluid_constant(i, iT_critical);
            CoolPropDbl Rcomponent = get_fluid_constant(i, igas_constant);
            tau_i = T_ci * tau / Tr;
            delta_i = delta * rhor / rho_ci;
            CoolPropDbl Rratio = Rcomponent / Rmix;

            // Cache the reducing temperature in some terms that need it (GERG-2004 models)
            components[i].EOS().alpha0.set_Tred(Tr);

            if (nTau == 0 && nDelta == 0) {
                double logxi = (std::abs(mole_fractions[i]) > DBL_EPSILON) ? log(mole_fractions[i]) : 0;
                summer += mole_fractions[i] * Rratio * (components[i].EOS().base0(tau_i, delta_i) + logxi);
            } else if (nTau == 0 && nDelta == 1) {
                summer += mole_fractions[i] * Rratio * rhor / rho_ci * components[i].EOS().dalpha0_dDelta(tau_i, delta_i);
            } else if (nTau == 1 && nDelta == 0) {
                summer += mole_fractions[i] * Rratio * T_ci / Tr * components[i].EOS().dalpha0_dTau(tau_i, delta_i);
            } else if (nTau == 0 && nDelta == 2) {
                summer += mole_fractions[i] * Rratio * pow(rhor / rho_ci, 2) * components[i].EOS().d2alpha0_dDelta2(tau_i, delta_i);
            } else if (nTau == 1 && nDelta == 1) {
                summer += mole_fractions[i] * Rratio * rhor / rho_ci * T_ci / Tr * components[i].EOS().d2alpha0_dDelta_dTau(tau_i, delta_i);
            } else if (nTau == 2 && nDelta == 0) {
                summer += mole_fractions[i] * Rratio * pow(T_ci / Tr, 2) * components[i].EOS().d2alpha0_dTau2(tau_i, delta_i);
            } else {
                throw ValueError();
            }
        }
        return summer;
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_alphar(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_alphar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dalphar_dDelta(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_dalphar_dDelta);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dalphar_dTau(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_dalphar_dTau);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alphar_dTau2(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d2alphar_dTau2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alphar_dDelta_dTau(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d2alphar_dDelta_dTau);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alphar_dDelta2(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d2alphar_dDelta2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alphar_dDelta3(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d3alphar_dDelta3);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alphar_dDelta2_dTau(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d3alphar_dDelta2_dTau);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alphar_dDelta_dTau2(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d3alphar_dDelta_dTau2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alphar_dTau3(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d3alphar_dTau3);
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_d4alphar_dDelta4(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d4alphar_dDelta4);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d4alphar_dDelta3_dTau(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d4alphar_dDelta3_dTau);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d4alphar_dDelta2_dTau2(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d4alphar_dDelta2_dTau2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d4alphar_dDelta_dTau3(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d4alphar_dDelta_dTau3);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d4alphar_dTau4(void) {
    calc_all_alphar_deriv_cache(mole_fractions, _tau, _delta);
    return static_cast<CoolPropDbl>(_d4alphar_dTau4);
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_alpha0(void) {
    const int nTau = 0, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dalpha0_dDelta(void) {
    const int nTau = 0, nDelta = 1;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dalpha0_dTau(void) {
    const int nTau = 1, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alpha0_dDelta2(void) {
    const int nTau = 0, nDelta = 2;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alpha0_dDelta_dTau(void) {
    const int nTau = 1, nDelta = 1;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d2alpha0_dTau2(void) {
    const int nTau = 2, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alpha0_dDelta3(void) {
    const int nTau = 0, nDelta = 3;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alpha0_dDelta2_dTau(void) {
    const int nTau = 1, nDelta = 2;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alpha0_dDelta_dTau2(void) {
    const int nTau = 2, nDelta = 1;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_d3alpha0_dTau3(void) {
    const int nTau = 3, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_first_saturation_deriv(parameters Of1, parameters Wrt1, HelmholtzEOSMixtureBackend& SatL,
                                                                    HelmholtzEOSMixtureBackend& SatV) {
    // Derivative of temperature w.r.t. pressure ALONG the saturation curve
    CoolPropDbl dTdP_sat = T() * (1 / SatV.rhomolar() - 1 / SatL.rhomolar()) / (SatV.hmolar() - SatL.hmolar());

    // "Trivial" inputs
    if (Of1 == iT && Wrt1 == iP) {
        return dTdP_sat;
    } else if (Of1 == iP && Wrt1 == iT) {
        return 1 / dTdP_sat;
    }
    // Derivative taken with respect to T
    else if (Wrt1 == iT) {
        return first_partial_deriv(Of1, iT, iP) + first_partial_deriv(Of1, iP, iT) / dTdP_sat;
    }
    // Derivative taken with respect to p
    else if (Wrt1 == iP) {
        return first_partial_deriv(Of1, iP, iT) + first_partial_deriv(Of1, iT, iP) * dTdP_sat;
    } else {
        throw ValueError(
          format("Not possible to take first saturation derivative with respect to %s", get_parameter_information(Wrt1, "short").c_str()));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_first_saturation_deriv(parameters Of1, parameters Wrt1) {
    if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for calc_first_saturation_deriv"));

    // Derivative of temperature w.r.t. pressure ALONG the saturation curve
    CoolPropDbl dTdP_sat = T() * (1 / SatV->rhomolar() - 1 / SatL->rhomolar()) / (SatV->hmolar() - SatL->hmolar());

    // "Trivial" inputs
    if (Of1 == iT && Wrt1 == iP) {
        return dTdP_sat;
    } else if (Of1 == iP && Wrt1 == iT) {
        return 1 / dTdP_sat;
    }
    // Derivative taken with respect to T
    else if (Wrt1 == iT) {
        return first_partial_deriv(Of1, iT, iP) + first_partial_deriv(Of1, iP, iT) / dTdP_sat;
    }
    // Derivative taken with respect to p
    else if (Wrt1 == iP) {
        return first_partial_deriv(Of1, iP, iT) + first_partial_deriv(Of1, iT, iP) * dTdP_sat;
    } else {
        throw ValueError(
          format("Not possible to take first saturation derivative with respect to %s", get_parameter_information(Wrt1, "short").c_str()));
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Wrt2) {
    if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for calc_second_saturation_deriv"));
    if (Wrt1 == iP && Wrt2 == iP) {
        CoolPropDbl dydT_constp = this->first_partial_deriv(Of1, iT, iP);
        CoolPropDbl d2ydTdp = this->second_partial_deriv(Of1, iT, iP, iP, iT);
        CoolPropDbl d2ydp2_constT = this->second_partial_deriv(Of1, iP, iT, iP, iT);
        CoolPropDbl d2ydT2_constp = this->second_partial_deriv(Of1, iT, iP, iT, iP);

        CoolPropDbl dTdp_along_sat = calc_first_saturation_deriv(iT, iP);
        CoolPropDbl dvdrhoL = -1 / POW2(SatL->rhomolar());
        CoolPropDbl dvdrhoV = -1 / POW2(SatV->rhomolar());
        CoolPropDbl DELTAv = 1 / SatV->rhomolar() - 1 / SatL->rhomolar();
        CoolPropDbl dDELTAv_dT_constp = dvdrhoV * SatV->first_partial_deriv(iDmolar, iT, iP) - dvdrhoL * SatL->first_partial_deriv(iDmolar, iT, iP);
        CoolPropDbl dDELTAv_dp_constT = dvdrhoV * SatV->first_partial_deriv(iDmolar, iP, iT) - dvdrhoL * SatL->first_partial_deriv(iDmolar, iP, iT);
        CoolPropDbl DELTAh = SatV->hmolar() - SatL->hmolar();
        CoolPropDbl dDELTAh_dT_constp = SatV->first_partial_deriv(iHmolar, iT, iP) - SatL->first_partial_deriv(iHmolar, iT, iP);
        CoolPropDbl dDELTAh_dp_constT = SatV->first_partial_deriv(iHmolar, iP, iT) - SatL->first_partial_deriv(iHmolar, iP, iT);
        CoolPropDbl ddT_dTdp_along_sat_constp = (DELTAh * (_T * dDELTAv_dT_constp + DELTAv) - _T * DELTAv * dDELTAh_dT_constp) / POW2(DELTAh);
        CoolPropDbl ddp_dTdp_along_sat_constT = (DELTAh * (_T * dDELTAv_dp_constT) - _T * DELTAv * dDELTAh_dp_constT) / POW2(DELTAh);

        double ddp_dydpsigma = d2ydp2_constT + dydT_constp * ddp_dTdp_along_sat_constT + d2ydTdp * dTdp_along_sat;
        double ddT_dydpsigma = d2ydTdp + dydT_constp * ddT_dTdp_along_sat_constp + d2ydT2_constp * dTdp_along_sat;
        return ddp_dydpsigma + ddT_dydpsigma * dTdp_along_sat;
    } else {
        throw ValueError(format("Currently, only possible to take second saturation derivative w.r.t. P (both times)"));
    }
}

CoolPropDbl HelmholtzEOSMixtureBackend::calc_second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2,
                                                                    parameters Constant2) {
    if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for calc_second_two_phase_deriv"));

    if (Of == iDmolar
        && ((Wrt1 == iHmolar && Constant1 == iP && Wrt2 == iP && Constant2 == iHmolar)
            || (Wrt2 == iHmolar && Constant2 == iP && Wrt1 == iP && Constant1 == iHmolar))) {
        parameters h_key = iHmolar, rho_key = iDmolar, p_key = iP;
        // taking the derivative of (drho/dv)*(dv/dh|p) with respect to p with h constant
        CoolPropDbl dv_dh_constp = calc_first_two_phase_deriv(rho_key, h_key, p_key) / (-POW2(rhomolar()));
        CoolPropDbl drhomolar_dp__consth = first_two_phase_deriv(rho_key, p_key, h_key);

        // Calculate the derivative of dvdh|p with respect to p at constant h
        CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(h_key, p_key, *SatL, *SatV);
        CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(h_key, p_key, *SatL, *SatV);
        CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(rho_key, p_key, *SatL, *SatV);
        CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(rho_key, p_key, *SatL, *SatV);
        CoolPropDbl numerator = 1 / SatV->keyed_output(rho_key) - 1 / SatL->keyed_output(rho_key);
        CoolPropDbl denominator = SatV->keyed_output(h_key) - SatL->keyed_output(h_key);
        CoolPropDbl dnumerator = -1 / POW2(SatV->keyed_output(rho_key)) * drhoV_dp_sat + 1 / POW2(SatL->keyed_output(rho_key)) * drhoL_dp_sat;
        CoolPropDbl ddenominator = dhV_dp_sat - dhL_dp_sat;
        CoolPropDbl d_dvdh_dp__consth = (denominator * dnumerator - numerator * ddenominator) / POW2(denominator);
        return -POW2(rhomolar()) * d_dvdh_dp__consth + dv_dh_constp * (-2 * rhomolar()) * drhomolar_dp__consth;
    } else if (Of == iDmass
               && ((Wrt1 == iHmass && Constant1 == iP && Wrt2 == iP && Constant2 == iHmass)
                   || (Wrt2 == iHmass && Constant2 == iP && Wrt1 == iP && Constant1 == iHmass))) {
        parameters h_key = iHmass, rho_key = iDmass, p_key = iP;
        CoolPropDbl rho = keyed_output(rho_key);
        // taking the derivative of (drho/dv)*(dv/dh|p) with respect to p with h constant
        CoolPropDbl dv_dh_constp = calc_first_two_phase_deriv(rho_key, h_key, p_key) / (-POW2(rho));
        CoolPropDbl drho_dp__consth = first_two_phase_deriv(rho_key, p_key, h_key);

        // Calculate the derivative of dvdh|p with respect to p at constant h
        CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(h_key, p_key, *SatL, *SatV);
        CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(h_key, p_key, *SatL, *SatV);
        CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(rho_key, p_key, *SatL, *SatV);
        CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(rho_key, p_key, *SatL, *SatV);
        CoolPropDbl numerator = 1 / SatV->keyed_output(rho_key) - 1 / SatL->keyed_output(rho_key);
        CoolPropDbl denominator = SatV->keyed_output(h_key) - SatL->keyed_output(h_key);
        CoolPropDbl dnumerator = -1 / POW2(SatV->keyed_output(rho_key)) * drhoV_dp_sat + 1 / POW2(SatL->keyed_output(rho_key)) * drhoL_dp_sat;
        CoolPropDbl ddenominator = dhV_dp_sat - dhL_dp_sat;
        CoolPropDbl d_dvdh_dp__consth = (denominator * dnumerator - numerator * ddenominator) / POW2(denominator);
        return -POW2(rho) * d_dvdh_dp__consth + dv_dh_constp * (-2 * rho) * drho_dp__consth;
    } else {
        throw ValueError();
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
    if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for calc_first_two_phase_deriv"));
    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        return -POW2(rhomolar()) * (1 / SatV->rhomolar() - 1 / SatL->rhomolar()) / (SatV->hmolar() - SatL->hmolar());
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return -POW2(rhomass()) * (1 / SatV->rhomass() - 1 / SatL->rhomass()) / (SatV->hmass() - SatL->hmass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        // v = 1/rho; drhodv = -rho^2; dvdrho = -1/rho^2
        CoolPropDbl dvdrhoL = -1 / POW2(SatL->rhomolar());
        CoolPropDbl dvdrhoV = -1 / POW2(SatV->rhomolar());
        CoolPropDbl dvL_dp = dvdrhoL * SatL->calc_first_saturation_deriv(iDmolar, iP, *SatL, *SatV);
        CoolPropDbl dvV_dp = dvdrhoV * SatV->calc_first_saturation_deriv(iDmolar, iP, *SatL, *SatV);
        CoolPropDbl dhL_dp = SatL->calc_first_saturation_deriv(iHmolar, iP, *SatL, *SatV);
        CoolPropDbl dhV_dp = SatV->calc_first_saturation_deriv(iHmolar, iP, *SatL, *SatV);
        CoolPropDbl dxdp_h = (Q() * dhV_dp + (1 - Q()) * dhL_dp) / (SatL->hmolar() - SatV->hmolar());
        CoolPropDbl dvdp_h = dvL_dp + dxdp_h * (1 / SatV->rhomolar() - 1 / SatL->rhomolar()) + Q() * (dvV_dp - dvL_dp);
        return -POW2(rhomolar()) * dvdp_h;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        // v = 1/rho; drhodv = -rho^2; dvdrho = -1/rho^2
        CoolPropDbl dvdrhoL = -1 / POW2(SatL->rhomass());
        CoolPropDbl dvdrhoV = -1 / POW2(SatV->rhomass());
        CoolPropDbl dvL_dp = dvdrhoL * SatL->calc_first_saturation_deriv(iDmass, iP, *SatL, *SatV);
        CoolPropDbl dvV_dp = dvdrhoV * SatV->calc_first_saturation_deriv(iDmass, iP, *SatL, *SatV);
        CoolPropDbl dhL_dp = SatL->calc_first_saturation_deriv(iHmass, iP, *SatL, *SatV);
        CoolPropDbl dhV_dp = SatV->calc_first_saturation_deriv(iHmass, iP, *SatL, *SatV);
        CoolPropDbl dxdp_h = (Q() * dhV_dp + (1 - Q()) * dhL_dp) / (SatL->hmass() - SatV->hmass());
        CoolPropDbl dvdp_h = dvL_dp + dxdp_h * (1 / SatV->rhomass() - 1 / SatL->rhomass()) + Q() * (dvV_dp - dvL_dp);
        return -POW2(rhomass()) * dvdp_h;
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
    }
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) {
    // Note: If you need all three values (drho_dh__p, drho_dp__h and rho_spline),
    // you should calculate drho_dp__h first to avoid duplicate calculations.

    bool drho_dh__p = false;
    bool drho_dp__h = false;
    bool rho_spline = false;

    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        drho_dh__p = true;
        if (_drho_spline_dh__constp) return _drho_spline_dh__constp;
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return first_two_phase_deriv_splined(iDmolar, iHmolar, iP, x_end) * POW2(molar_mass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        drho_dp__h = true;
        if (_drho_spline_dp__consth) return _drho_spline_dp__consth;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        return first_two_phase_deriv_splined(iDmolar, iP, iHmolar, x_end) * molar_mass();
    }
    // Add the special case for the splined density
    else if (Of == iDmolar && Wrt == iDmolar && Constant == iDmolar) {
        rho_spline = true;
        if (_rho_spline) return _rho_spline;
    } else if (Of == iDmass && Wrt == iDmass && Constant == iDmass) {
        return first_two_phase_deriv_splined(iDmolar, iDmolar, iDmolar, x_end) * molar_mass();
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
    }

    if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for calc_first_two_phase_deriv_splined"));
    if (_Q > x_end) {
        throw ValueError(format("Q [%g] is greater than x_end [%Lg]", _Q, x_end).c_str());
    }
    if (_phase != iphase_twophase) {
        throw ValueError(format("state is not two-phase"));
    }

    shared_ptr<HelmholtzEOSMixtureBackend> Liq(new HelmholtzEOSMixtureBackend(this->get_components())),
      End(new HelmholtzEOSMixtureBackend(this->get_components()));

    Liq->specify_phase(iphase_liquid);
    Liq->_Q = -1;
    Liq->update_DmolarT_direct(SatL->rhomolar(), SatL->T());
    End->update(QT_INPUTS, x_end, SatL->T());

    CoolPropDbl Delta = Q() * (SatV->keyed_output(iHmolar) - SatL->keyed_output(iHmolar));
    CoolPropDbl Delta_end = End->keyed_output(iHmolar) - SatL->keyed_output(iHmolar);

    // At the end of the zone to which spline is applied
    CoolPropDbl drho_dh_end = End->calc_first_two_phase_deriv(iDmolar, iHmolar, iP);
    CoolPropDbl rho_end = End->keyed_output(iDmolar);

    // Faking single-phase
    CoolPropDbl rho_liq = Liq->keyed_output(iDmolar);
    CoolPropDbl drho_dh_liq__constp = Liq->first_partial_deriv(iDmolar, iHmolar, iP);

    // Spline coordinates a, b, c, d
    CoolPropDbl Abracket = (2 * rho_liq - 2 * rho_end + Delta_end * (drho_dh_liq__constp + drho_dh_end));
    CoolPropDbl a = 1 / POW3(Delta_end) * Abracket;
    CoolPropDbl b = 3 / POW2(Delta_end) * (-rho_liq + rho_end) - 1 / Delta_end * (drho_dh_end + 2 * drho_dh_liq__constp);
    CoolPropDbl c = drho_dh_liq__constp;
    CoolPropDbl d = rho_liq;

    // Either the spline value or drho/dh|p can be directly evaluated now
    _rho_spline = a * POW3(Delta) + b * POW2(Delta) + c * Delta + d;
    _drho_spline_dh__constp = 3 * a * POW2(Delta) + 2 * b * Delta + c;
    if (rho_spline) return _rho_spline;
    if (drho_dh__p) return _drho_spline_dh__constp;

    // It's drho/dp|h
    // ... calculate some more things

    // Derivatives *along* the saturation curve using the special internal method
    CoolPropDbl dhL_dp_sat = SatL->calc_first_saturation_deriv(iHmolar, iP, *SatL, *SatV);
    CoolPropDbl dhV_dp_sat = SatV->calc_first_saturation_deriv(iHmolar, iP, *SatL, *SatV);
    CoolPropDbl drhoL_dp_sat = SatL->calc_first_saturation_deriv(iDmolar, iP, *SatL, *SatV);
    CoolPropDbl drhoV_dp_sat = SatV->calc_first_saturation_deriv(iDmolar, iP, *SatL, *SatV);
    CoolPropDbl rhoV = SatV->keyed_output(iDmolar);
    CoolPropDbl rhoL = SatL->keyed_output(iDmolar);
    CoolPropDbl drho_dp_end = POW2(End->keyed_output(iDmolar)) * (x_end / POW2(rhoV) * drhoV_dp_sat + (1 - x_end) / POW2(rhoL) * drhoL_dp_sat);

    // Faking single-phase
    //CoolPropDbl drho_dp__consth_liq = Liq->first_partial_deriv(iDmolar, iP, iHmolar);
    CoolPropDbl d2rhodhdp_liq = Liq->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);  // ?

    // Derivatives at the end point
    // CoolPropDbl drho_dp__consth_end = End->calc_first_two_phase_deriv(iDmolar, iP, iHmolar);
    CoolPropDbl d2rhodhdp_end = End->calc_second_two_phase_deriv(iDmolar, iHmolar, iP, iP, iHmolar);

    // Reminder:
    // Delta = Q()*(hV-hL) = h-hL
    // Delta_end = x_end*(hV-hL);
    CoolPropDbl d_Delta_dp__consth = -dhL_dp_sat;
    CoolPropDbl d_Delta_end_dp__consth = x_end * (dhV_dp_sat - dhL_dp_sat);

    // First pressure derivative at constant h of the coefficients a,b,c,d
    // CoolPropDbl Abracket = (2*rho_liq - 2*rho_end + Delta_end * (drho_dh_liq__constp + drho_dh_end));
    CoolPropDbl d_Abracket_dp_consth = (2 * drhoL_dp_sat - 2 * drho_dp_end + Delta_end * (d2rhodhdp_liq + d2rhodhdp_end)
                                        + d_Delta_end_dp__consth * (drho_dh_liq__constp + drho_dh_end));
    CoolPropDbl da_dp = 1 / POW3(Delta_end) * d_Abracket_dp_consth + Abracket * (-3 / POW4(Delta_end) * d_Delta_end_dp__consth);
    CoolPropDbl db_dp = -6 / POW3(Delta_end) * d_Delta_end_dp__consth * (rho_end - rho_liq) + 3 / POW2(Delta_end) * (drho_dp_end - drhoL_dp_sat)
                        + (1 / POW2(Delta_end) * d_Delta_end_dp__consth) * (drho_dh_end + 2 * drho_dh_liq__constp)
                        - (1 / Delta_end) * (d2rhodhdp_end + 2 * d2rhodhdp_liq);
    CoolPropDbl dc_dp = d2rhodhdp_liq;
    CoolPropDbl dd_dp = drhoL_dp_sat;

    _drho_spline_dp__consth =
      (3 * a * POW2(Delta) + 2 * b * Delta + c) * d_Delta_dp__consth + POW3(Delta) * da_dp + POW2(Delta) * db_dp + Delta * dc_dp + dd_dp;
    if (drho_dp__h) return _drho_spline_dp__consth;

    throw ValueError("Something went wrong in HelmholtzEOSMixtureBackend::calc_first_two_phase_deriv_splined");
    return _HUGE;
}

CoolProp::CriticalState HelmholtzEOSMixtureBackend::calc_critical_point(double rho0, double T0) {
    class Resid : public FuncWrapperND
    {
       public:
        HelmholtzEOSMixtureBackend& HEOS;
        double L1, M1;
        Eigen::MatrixXd Lstar, Mstar;
        Resid(HelmholtzEOSMixtureBackend& HEOS) : HEOS(HEOS), L1(_HUGE), M1(_HUGE){};
        std::vector<double> call(const std::vector<double>& tau_delta) {
            double rhomolar = tau_delta[1] * HEOS.rhomolar_reducing();
            double T = HEOS.T_reducing() / tau_delta[0];
            HEOS.update(DmolarT_INPUTS, rhomolar, T);
            Lstar = MixtureDerivatives::Lstar(HEOS, XN_INDEPENDENT);
            Mstar = MixtureDerivatives::Mstar(HEOS, XN_INDEPENDENT, Lstar);
            std::vector<double> o(2);
            o[0] = Lstar.determinant();
            o[1] = Mstar.determinant();
            return o;
        };
        std::vector<std::vector<double>> Jacobian(const std::vector<double>& x) {
            std::size_t N = x.size();
            std::vector<std::vector<double>> J(N, std::vector<double>(N, 0));
            Eigen::MatrixXd adjL = adjugate(Lstar), adjM = adjugate(Mstar), dLdTau = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iTau),
                            dLdDelta = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iDelta),
                            dMdTau = MixtureDerivatives::dMstar_dX(HEOS, XN_INDEPENDENT, iTau, Lstar, dLdTau),
                            dMdDelta = MixtureDerivatives::dMstar_dX(HEOS, XN_INDEPENDENT, iDelta, Lstar, dLdDelta);

            J[0][0] = (adjL * dLdTau).trace();
            J[0][1] = (adjL * dLdDelta).trace();
            J[1][0] = (adjM * dMdTau).trace();
            J[1][1] = (adjM * dMdDelta).trace();
            return J;
        }
        /// Not used, for testing purposes
        std::vector<std::vector<double>> numerical_Jacobian(const std::vector<double>& x) {
            std::size_t N = x.size();
            std::vector<double> rplus, rminus, xp;
            std::vector<std::vector<double>> J(N, std::vector<double>(N, 0));

            double epsilon = 0.0001;

            // Build the Jacobian by column
            for (std::size_t i = 0; i < N; ++i) {
                xp = x;
                xp[i] += epsilon;
                rplus = call(xp);
                xp[i] -= 2 * epsilon;
                rminus = call(xp);

                for (std::size_t j = 0; j < N; ++j) {
                    J[j][i] = (rplus[j] - rminus[j]) / (2 * epsilon);
                }
            }
            std::cout << J[0][0] << " " << J[0][1] << std::endl;
            std::cout << J[1][0] << " " << J[1][1] << std::endl;
            return J;
        };
    };
    Resid resid(*this);
    std::vector<double> x, tau_delta(2);
    tau_delta[0] = T_reducing() / T0;
    tau_delta[1] = rho0 / rhomolar_reducing();
    x = NDNewtonRaphson_Jacobian(&resid, tau_delta, 1e-10, 100);
    _critical.T = T_reducing() / x[0];
    _critical.rhomolar = x[1] * rhomolar_reducing();
    _critical.p = calc_pressure_nocache(_critical.T, _critical.rhomolar);

    CriticalState critical;
    critical.T = _critical.T;
    critical.p = _critical.p;
    critical.rhomolar = _critical.rhomolar;
    if (_critical.p < 0) {
        critical.stable = false;
    } else {
        if (get_config_bool(ASSUME_CRITICAL_POINT_STABLE)) {
            critical.stable = true;
        } else {
            // Otherwise we try to check stability with TPD-based analysis
            StabilityRoutines::StabilityEvaluationClass stability_tester(*this);
            critical.stable = stability_tester.is_stable();
        }
    }
    return critical;
}

/**
 * \brief This class is the objective function for the one-dimensional solver used to find the first intersection with the L1*=0 contour
 */
class OneDimObjective : public FuncWrapper1DWithTwoDerivs
{
   public:
    CoolProp::HelmholtzEOSMixtureBackend& HEOS;
    const double delta;
    double _call, _deriv, _second_deriv;
    OneDimObjective(HelmholtzEOSMixtureBackend& HEOS, double delta0) : HEOS(HEOS), delta(delta0), _call(_HUGE), _deriv(_HUGE), _second_deriv(_HUGE){};
    double call(double tau) {
        double rhomolar = HEOS.rhomolar_reducing() * delta, T = HEOS.T_reducing() / tau;
        HEOS.update_DmolarT_direct(rhomolar, T);
        _call = MixtureDerivatives::Lstar(HEOS, XN_INDEPENDENT).determinant();
        return _call;
    }
    double deriv(double tau) {
        Eigen::MatrixXd adjL = adjugate(MixtureDerivatives::Lstar(HEOS, XN_INDEPENDENT)),
                        dLdTau = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iTau);
        _deriv = (adjL * dLdTau).trace();
        return _deriv;
    };
    double second_deriv(double tau) {
        Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(HEOS, XN_INDEPENDENT),
                        dLstardTau = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iTau),
                        d2LstardTau2 = MixtureDerivatives::d2Lstar_dX2(HEOS, XN_INDEPENDENT, iTau, iTau), adjL = adjugate(Lstar),
                        dadjLstardTau = adjugate_derivative(Lstar, dLstardTau);
        _second_deriv = (dLstardTau * dadjLstardTau + adjL * d2LstardTau2).trace();
        return _second_deriv;
    };
};

/** This class is used to trace the spinodal of the mixture, and is also used to calculate critical points
 */
class L0CurveTracer : public FuncWrapper1DWithDeriv
{
   public:
    CoolProp::HelmholtzEOSMixtureBackend& HEOS;
    double delta, tau,
      M1_last,         ///< The last value that the Mstar determinant had
      theta_last,      ///< The last value that the angle had
      R_tau,           ///< The radius for tau currently being used
      R_delta,         ///< The radius for delta currently being used
      R_tau_tracer,    ///< The radius for tau that should be used in the L1*=0 tracer (user-modifiable after instantiation)
      R_delta_tracer;  ///< The radius for delta that should be used in the L1*=0 tracer (user-modifiable after instantiation)
    std::vector<CoolProp::CriticalState> critical_points;
    int N_critical_points;
    Eigen::MatrixXd Lstar, adjLstar, dLstardTau, d2LstardTau2, dLstardDelta;
    SpinodalData spinodal_values;
    bool
      find_critical_points;  ///< If true, actually calculate the critical points, otherwise, skip evaluation of critical points but still trace the spinodal
    L0CurveTracer(HelmholtzEOSMixtureBackend& HEOS, double tau0, double delta0)
      : HEOS(HEOS), delta(delta0), tau(tau0), M1_last(_HUGE), N_critical_points(0), find_critical_points(true) {
        R_delta_tracer = 0.1;
        R_delta = R_delta_tracer;
        R_tau_tracer = 0.1;
        R_tau = R_tau_tracer;
    };
    /***
     \brief Update values for tau, delta
     @param theta The angle
     @param tau The old value of tau
     @param delta The old value of delta
     @param tau_new The new value of tau
     @param delta_new The new value of delta
     */
    void get_tau_delta(const double theta, const double tau, const double delta, double& tau_new, double& delta_new) {
        tau_new = tau + R_tau * cos(theta);
        delta_new = delta + R_delta * sin(theta);
    };
    /***
     \brief Calculate the value of L1
     @param theta The angle
     */
    double call(double theta) {
        double tau_new, delta_new;
        this->get_tau_delta(theta, tau, delta, tau_new, delta_new);
        double rhomolar = HEOS.rhomolar_reducing() * delta_new, T = HEOS.T_reducing() / tau_new;
        HEOS.update_DmolarT_direct(rhomolar, T);
        Lstar = MixtureDerivatives::Lstar(HEOS, XN_INDEPENDENT);
        adjLstar = adjugate(Lstar);
        dLstardTau = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iTau);
        dLstardDelta = MixtureDerivatives::dLstar_dX(HEOS, XN_INDEPENDENT, iDelta);
        double L1 = Lstar.determinant();
        return L1;
    };
    /***
     \brief Calculate the first partial derivative of L1 with respect to the angle
     @param theta The angle
     */
    double deriv(double theta) {
        double dL1_dtau = (adjLstar * dLstardTau).trace(), dL1_ddelta = (adjLstar * dLstardDelta).trace();
        return -R_tau * sin(theta) * dL1_dtau + R_delta * cos(theta) * dL1_ddelta;
    };

    void trace() {
        bool debug = (get_debug_level() > 0) | false;
        double theta;
        for (int i = 0; i < 300; ++i) {
            if (i == 0) {
                // In the first iteration, search all angles in the positive delta direction using a
                // bounded solver with a very small radius in order to not hit other L1*=0 contours
                // that are in the vicinity
                R_tau = 0.001;
                R_delta = 0.001;
                theta = Brent(this, 0, M_PI, DBL_EPSILON, 1e-10, 100);
            } else {
                // In subsequent iterations, you already have an excellent guess for the direction to
                // be searching in, use Newton's method to refine the solution since we also
                // have an analytic solution for the derivative
                R_tau = R_tau_tracer;
                R_delta = R_delta_tracer;
                try {
                    theta = Newton(this, theta_last, 1e-10, 100);
                } catch (...) {
                    if (N_critical_points > 0 && delta > 1.5 * critical_points[0].rhomolar / HEOS.rhomolar_reducing()) {
                        // Stopping the search - probably we have a kink in the L1*=0 contour
                        // caused by poor low-temperature behavior
                        continue;
                    }
                }

                // If the solver takes a U-turn, going in the opposite direction of travel
                // this is not a good thing, and force a Brent's method solver call to make sure we keep
                // tracing in the same direction
                if (std::abs(angle_difference(theta, theta_last)) > M_PI / 2.0) {
                    // We have found at least one critical point and we are now well above the density
                    // associated with the first critical point
                    if (N_critical_points > 0 && delta > 1.2 * critical_points[0].rhomolar / HEOS.rhomolar_reducing()) {
                        // Stopping the search - probably we have a kink in the L1*=0 contour
                        // caused by poor low-temperature behavior
                        continue;
                    } else {
                        theta = Brent(this, theta_last - M_PI / 3.5, theta_last + M_PI / 3.5, DBL_EPSILON, 1e-10, 100);
                    }
                }
            }

            // Calculate the second criticality condition
            double M1 = MixtureDerivatives::Mstar(HEOS, XN_INDEPENDENT, Lstar).determinant();
            double p_MPa = HEOS.p() / 1e6;

            // Calculate the new tau and delta at the new point
            double tau_new, delta_new;
            this->get_tau_delta(theta, tau, delta, tau_new, delta_new);

            // Stop if bounds are exceeded
            if (p_MPa > 500 || HEOS.get_critical_is_terminated(delta_new, tau_new)) {
                break;
            }

            // If the sign of M1 and the previous value of M1 have different signs, it means that
            // you have bracketed a critical point, run the full critical point solver to
            // find the critical point and store it
            // Only enabled if find_critical_points is true (the default)
            if (i > 0 && M1 * M1_last < 0 && find_critical_points) {
                double rhomolar = HEOS.rhomolar_reducing() * (delta + delta_new) / 2.0, T = HEOS.T_reducing() / ((tau + tau_new) / 2.0);
                CoolProp::CriticalState crit = HEOS.calc_critical_point(rhomolar, T);
                critical_points.push_back(crit);
                N_critical_points++;
                if (debug) {
                    std::cout << HEOS.get_mole_fractions()[0] << " " << crit.rhomolar << " " << crit.T << " " << p_MPa << std::endl;
                }
            }

            // Update the storage values
            this->tau = tau_new;
            this->delta = delta_new;
            this->M1_last = M1;
            this->theta_last = theta;

            this->spinodal_values.tau.push_back(tau_new);
            this->spinodal_values.delta.push_back(delta_new);
            this->spinodal_values.M1.push_back(M1);
        }
    };
};

void HelmholtzEOSMixtureBackend::calc_criticality_contour_values(double& L1star, double& M1star) {
    Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(*this, XN_INDEPENDENT);
    Eigen::MatrixXd Mstar = MixtureDerivatives::Mstar(*this, XN_INDEPENDENT, Lstar);
    L1star = Lstar.determinant();
    M1star = Mstar.determinant();
};

void HelmholtzEOSMixtureBackend::get_critical_point_search_radii(double& R_delta, double& R_tau) {
    R_delta = 0.025;
    R_tau = 0.1;
}
std::vector<CoolProp::CriticalState> HelmholtzEOSMixtureBackend::_calc_all_critical_points(bool find_critical_points) {
    // Populate the temporary class used to calculate the critical point(s)
    add_critical_state();
    if (get_debug_level() > 10) {
        rapidjson::Document doc;
        doc.SetObject();
        rapidjson::Value& val = doc;
        std::vector<std::vector<DepartureFunctionPointer>>& mat = critical_state->residual_helmholtz->Excess.DepartureFunctionMatrix;
        if (mat.size() > 0) {
            mat[0][1]->phi.to_json(val, doc);
            std::cout << cpjson::to_string(doc);
        }
    }
    critical_state->set_mole_fractions(this->get_mole_fractions_ref());

    // Specify state to be something homogeneous to shortcut phase evaluation
    critical_state->specify_phase(iphase_gas);

    double delta0 = _HUGE, tau0 = _HUGE;
    critical_state->get_critical_point_starting_values(delta0, tau0);

    OneDimObjective resid_L0(*critical_state, delta0);

    // If the derivative of L1star with respect to tau is positive,
    // tau needs to be increased such that we sit on the other
    // side of the d(L1star)/dtau = 0 contour
    resid_L0.call(tau0);
    int bump_count = 0;
    while (resid_L0.deriv(tau0) > 0 && bump_count < 3) {
        tau0 *= 1.1;
        bump_count++;
    }
    double tau_L0 = Halley(resid_L0, tau0, 1e-10, 100);
    //double T0 = T_reducing()/tau_L0;
    //double rho0 = delta0*rhomolar_reducing();

    L0CurveTracer tracer(*critical_state, tau_L0, delta0);
    tracer.find_critical_points = find_critical_points;

    double R_delta = 0, R_tau = 0;
    critical_state->get_critical_point_search_radii(R_delta, R_tau);
    tracer.R_delta_tracer = R_delta;
    tracer.R_tau_tracer = R_tau;
    tracer.trace();

    this->spinodal_values = tracer.spinodal_values;

    return tracer.critical_points;
}

double HelmholtzEOSMixtureBackend::calc_tangent_plane_distance(const double T, const double p, const std::vector<double>& w,
                                                               const double rhomolar_guess) {

    const std::vector<CoolPropDbl>& z = this->get_mole_fractions_ref();
    if (w.size() != z.size()) {
        throw ValueError(format("Trial composition vector size [%d] is not the same as bulk composition [%d]", w.size(), z.size()));
    }
    add_TPD_state();
    TPD_state->set_mole_fractions(w);

    CoolPropDbl rho = TPD_state->solver_rho_Tp_global(T, p, 0.9 / TPD_state->SRK_covolume());
    TPD_state->update_DmolarT_direct(rho, T);

    double summer = 0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        summer +=
          w[i] * (log(MixtureDerivatives::fugacity_i(*TPD_state, i, XN_DEPENDENT)) - log(MixtureDerivatives::fugacity_i(*this, i, XN_DEPENDENT)));
    }
    return summer;
}

void HelmholtzEOSMixtureBackend::calc_build_spinodal() {
    // Ok, we are faking a little bit here, hijacking the code for critical points, but skipping evaluation of critical points
    bool find_critical_points = false;
    _calc_all_critical_points(find_critical_points);
}

void HelmholtzEOSMixtureBackend::set_reference_stateS(const std::string& reference_state) {
    for (std::size_t i = 0; i < components.size(); ++i) {
        CoolProp::HelmholtzEOSMixtureBackend HEOS(std::vector<CoolPropFluid>(1, components[i]));
        if (!reference_state.compare("IIR")) {
            if (HEOS.Ttriple() > 273.15) {
                throw ValueError(format("Cannot use IIR reference state; Ttriple [%Lg] is greater than 273.15 K", HEOS.Ttriple()));
            }
            HEOS.update(QT_INPUTS, 0, 273.15);

            // Get current values for the enthalpy and entropy
            double deltah = HEOS.hmass() - 200000;  // offset from 200000 J/kg enthalpy
            double deltas = HEOS.smass() - 1000;    // offset from 1000 J/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(components[i], delta_a1, delta_a2, "IIR");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("ASHRAE")) {
            if (HEOS.Ttriple() > 233.15) {
                throw ValueError(format("Cannot use ASHRAE reference state; Ttriple [%Lg] is greater than than 233.15 K", HEOS.Ttriple()));
            }
            HEOS.update(QT_INPUTS, 0, 233.15);

            // Get current values for the enthalpy and entropy
            double deltah = HEOS.hmass() - 0;  // offset from 0 J/kg enthalpy
            double deltas = HEOS.smass() - 0;  // offset from 0 J/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(components[i], delta_a1, delta_a2, "ASHRAE");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("NBP")) {
            if (HEOS.p_triple() > 101325) {
                throw ValueError(format("Cannot use NBP reference state; p_triple [%Lg Pa] is greater than than 101325 Pa", HEOS.p_triple()));
            }
            // Saturated liquid boiling point at 1 atmosphere
            HEOS.update(PQ_INPUTS, 101325, 0);

            double deltah = HEOS.hmass() - 0;  // offset from 0 kJ/kg enthalpy
            double deltas = HEOS.smass() - 0;  // offset from 0 kJ/kg/K entropy
            double delta_a1 = deltas / (HEOS.gas_constant() / HEOS.molar_mass());
            double delta_a2 = -deltah / (HEOS.gas_constant() / HEOS.molar_mass() * HEOS.get_reducing_state().T);
            // Change the value in the library for the given fluid
            set_fluid_enthalpy_entropy_offset(components[i], delta_a1, delta_a2, "NBP");
            if (get_debug_level() > 0) {
                std::cout << format("set offsets to %0.15g and %0.15g\n", delta_a1, delta_a2);
            }
        } else if (!reference_state.compare("DEF")) {
            set_fluid_enthalpy_entropy_offset(components[i], 0, 0, "DEF");
        } else if (!reference_state.compare("RESET")) {
            set_fluid_enthalpy_entropy_offset(components[i], 0, 0, "RESET");
        } else {
            throw ValueError(format("reference state string is invalid: [%s]", reference_state.c_str()));
        }
    }
}

/// Set the reference state based on a thermodynamic state point specified by temperature and molar density
/// @param T Temperature at reference state [K]
/// @param rhomolar Molar density at reference state [mol/m^3]
/// @param hmolar0 Molar enthalpy at reference state [J/mol]
/// @param smolar0 Molar entropy at reference state [J/mol/K]
void HelmholtzEOSMixtureBackend::set_reference_stateD(double T, double rhomolar, double hmolar0, double smolar0) {
    for (std::size_t i = 0; i < components.size(); ++i) {
        CoolProp::HelmholtzEOSMixtureBackend HEOS(std::vector<CoolPropFluid>(1, components[i]));

        HEOS.update(DmolarT_INPUTS, rhomolar, T);

        // Get current values for the enthalpy and entropy
        double deltah = HEOS.hmolar() - hmolar0;  // offset from specified enthalpy in J/mol
        double deltas = HEOS.smolar() - smolar0;  // offset from specified entropy in J/mol/K
        double delta_a1 = deltas / (HEOS.gas_constant());
        double delta_a2 = -deltah / (HEOS.gas_constant() * HEOS.get_reducing_state().T);
        set_fluid_enthalpy_entropy_offset(components[i], delta_a1, delta_a2, "custom");
    }
}

void HelmholtzEOSMixtureBackend::set_fluid_enthalpy_entropy_offset(CoolPropFluid& component, double delta_a1, double delta_a2,
                                                                   const std::string& ref) {
    component.EOS().alpha0.EnthalpyEntropyOffset.set(delta_a1, delta_a2, ref);

    shared_ptr<CoolProp::HelmholtzEOSBackend> HEOS(new CoolProp::HelmholtzEOSBackend(component));
    HEOS->specify_phase(iphase_gas);  // Something homogeneous;
                                      // Calculate the new enthalpy and entropy values
    HEOS->update(DmolarT_INPUTS, component.EOS().hs_anchor.rhomolar, component.EOS().hs_anchor.T);
    component.EOS().hs_anchor.hmolar = HEOS->hmolar();
    component.EOS().hs_anchor.smolar = HEOS->smolar();

    double f = (HEOS->name() == "Water" || HEOS->name() == "CarbonDioxide") ? 1.00001 : 1.0;

    // Calculate the new enthalpy and entropy values at the reducing state
    HEOS->update(DmolarT_INPUTS, component.EOS().reduce.rhomolar * f, component.EOS().reduce.T * f);
    component.EOS().reduce.hmolar = HEOS->hmolar();
    component.EOS().reduce.smolar = HEOS->smolar();

    // Calculate the new enthalpy and entropy values at the critical state
    HEOS->update(DmolarT_INPUTS, component.crit.rhomolar * f, component.crit.T * f);
    component.crit.hmolar = HEOS->hmolar();
    component.crit.smolar = HEOS->smolar();

    // Calculate the new enthalpy and entropy values
    HEOS->update(DmolarT_INPUTS, component.triple_liquid.rhomolar, component.triple_liquid.T);
    component.triple_liquid.hmolar = HEOS->hmolar();
    component.triple_liquid.smolar = HEOS->smolar();

    // Calculate the new enthalpy and entropy values
    HEOS->update(DmolarT_INPUTS, component.triple_vapor.rhomolar, component.triple_vapor.T);
    component.triple_vapor.hmolar = HEOS->hmolar();
    component.triple_vapor.smolar = HEOS->smolar();

    if (!HEOS->is_pure()) {
        // Calculate the new enthalpy and entropy values
        HEOS->update(DmolarT_INPUTS, component.EOS().max_sat_T.rhomolar, component.EOS().max_sat_T.T);
        component.EOS().max_sat_T.hmolar = HEOS->hmolar();
        component.EOS().max_sat_T.smolar = HEOS->smolar();
        // Calculate the new enthalpy and entropy values
        HEOS->update(DmolarT_INPUTS, component.EOS().max_sat_p.rhomolar, component.EOS().max_sat_p.T);
        component.EOS().max_sat_p.hmolar = HEOS->hmolar();
        component.EOS().max_sat_p.smolar = HEOS->smolar();
    }
}

} /* namespace CoolProp */
