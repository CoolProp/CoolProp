#include "coolpropsolver.h"
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "CPState.h"
#include <iostream>
#include <string>
#include <stdlib.h>

CoolPropSolver::CoolPropSolver(const std::string& mediumName, const std::string& libraryName, const std::string& substanceName)
  : BaseSolver(mediumName, libraryName, substanceName) {

    // Fluid name can be used to pass in other parameters.
    // The string can be composed like "Propane|enable_TTSE=1|calc_transport=0"
    std::vector<std::string> name_options = strsplit(substanceName, '|');

    // Set the defaults
    fluidType = -1;
    enable_TTSE = false;
    debug_level = 0;
    calc_transport = false;
    extend_twophase = false;
    twophase_derivsmoothing_xend = 0;
    rho_smoothing_xend = 0;

    if (name_options.size() > 1) {
        for (unsigned int i = 1; i < name_options.size(); i++) {
            // Split around the equals sign
            std::vector<std::string> param_val = strsplit(name_options[i], '=');
            if (param_val.size() != 2) {
                errorMessage((char*)format("Could not parse the option [%s], must be in the form param=value", name_options[i].c_str()).c_str());
            }

            // Check each of the options in turn
            if (!param_val[0].compare("enable_TTSE")) {
                if (!param_val[1].compare("1") || !param_val[1].compare("true")) {
                    std::cout << "TTSE is on\n";
                    enable_TTSE = true;
                } else if (!param_val[1].compare("0") || !param_val[1].compare("false")) {
                    std::cout << "TTSE is off\n";
                    enable_TTSE = false;
                } else
                    errorMessage((char*)format("I don't know how to handle this option [%s]", name_options[i].c_str()).c_str());
                //throw NotImplementedError((char*)format("I don't know how to handle this option [%s]",name_options[i].c_str()).c_str());
            } else if (!param_val[0].compare("calc_transport")) {
                if (!param_val[1].compare("1") || !param_val[1].compare("true"))
                    calc_transport = true;
                else if (!param_val[1].compare("0") || !param_val[1].compare("false"))
                    calc_transport = false;
                else
                    errorMessage((char*)format("I don't know how to handle this option [%s]", name_options[i].c_str()).c_str());
            } else if (!param_val[0].compare("enable_EXTTP")) {
                if (!param_val[1].compare("1") || !param_val[1].compare("true"))
                    extend_twophase = true;
                else if (!param_val[1].compare("0") || !param_val[1].compare("false"))
                    extend_twophase = false;
                else
                    errorMessage((char*)format("I don't know how to handle this option [%s]", name_options[i].c_str()).c_str());
            } else if (!param_val[0].compare("twophase_derivsmoothing_xend")) {
                twophase_derivsmoothing_xend = strtod(param_val[1].c_str(), NULL);
                if (twophase_derivsmoothing_xend < 0 || twophase_derivsmoothing_xend > 1)
                    errorMessage(
                      (char*)format("I don't know how to handle this twophase_derivsmoothing_xend value [%d]", param_val[0].c_str()).c_str());
            } else if (!param_val[0].compare("rho_smoothing_xend")) {
                rho_smoothing_xend = strtod(param_val[1].c_str(), NULL);
                if (rho_smoothing_xend < 0 || rho_smoothing_xend > 1)
                    errorMessage((char*)format("I don't know how to handle this rho_smoothing_xend value [%d]", param_val[0].c_str()).c_str());
            } else if (!param_val[0].compare("debug")) {
                debug_level = (int)strtol(param_val[1].c_str(), NULL, 0);
                if (debug_level < 0 || debug_level > 1000)
                    errorMessage((char*)format("I don't know how to handle this debug level [%s]", param_val[0].c_str()).c_str());
            } else {
                errorMessage((char*)format("This option [%s] was not understood", name_options[i].c_str()).c_str());
            }

            // Some options were passed in, lets see what we have
            std::cout << param_val[0] << " has the value of " << param_val[1] << std::endl;
        }
    }
    // Handle the name and fill the fluid type
    if (debug_level > 5) std::cout << "Checking fluid " << name_options[0] << " against database." << std::endl;
    fluidType = getFluidType(name_options[0]);  // Throws an error if unknown fluid
    if (debug_level > 5) std::cout << "Check passed, reducing " << substanceName << " to " << name_options[0] << std::endl;
    this->substanceName = name_options[0];
    state = new CoolPropStateClassSI(name_options[0]);
    setFluidConstants();
}

void CoolPropSolver::setFluidConstants() {
    if ((fluidType == FLUID_TYPE_PURE) || (fluidType == FLUID_TYPE_PSEUDOPURE) || (fluidType == FLUID_TYPE_REFPROP)) {
        if (debug_level > 5) std::cout << format("Setting constants for fluid %s \n", substanceName.c_str());
        _fluidConstants.pc = PropsSI((char*)"pcrit", (char*)"T", 0, (char*)"P", 0, (char*)substanceName.c_str());
        _fluidConstants.Tc = PropsSI((char*)"Tcrit", (char*)"T", 0, (char*)"P", 0, (char*)substanceName.c_str());
        _fluidConstants.MM = PropsSI((char*)"molemass", (char*)"T", 0, (char*)"P", 0, (char*)substanceName.c_str());
        _fluidConstants.dc = PropsSI((char*)"rhocrit", (char*)"T", 0, (char*)"P", 0, (char*)substanceName.c_str());
        return;
    }
    if ((fluidType == FLUID_TYPE_INCOMPRESSIBLE_LIQUID) || (fluidType == FLUID_TYPE_INCOMPRESSIBLE_SOLUTION)) {
        if (debug_level > 5) std::cout << format("Setting constants for incompressible fluid %s \n", substanceName.c_str());
        _fluidConstants.pc = -1;
        _fluidConstants.Tc = -1;
        _fluidConstants.MM = -1;
        _fluidConstants.dc = -1;
        return;
    }
}

void CoolPropSolver::preStateChange(void) {
    /// Some common code to avoid pitfalls from incompressibles
    if ((fluidType == FLUID_TYPE_PURE) || (fluidType == FLUID_TYPE_PSEUDOPURE) || (fluidType == FLUID_TYPE_REFPROP)) {
        try {
            if (enable_TTSE)
                state->enable_TTSE_LUT();
            else
                state->disable_TTSE_LUT();

            if (extend_twophase)
                state->enable_EXTTP();
            else
                state->disable_EXTTP();
        } catch (std::exception& e) {
            errorMessage((char*)e.what());
            std::cout << format("Exception from state object: %s \n", (char*)e.what());
        }
    }
}

void CoolPropSolver::postStateChange(ExternalThermodynamicState* const properties) {
    /// Some common code to avoid pitfalls from incompressibles
    switch (fluidType) {
        case FLUID_TYPE_PURE:
        case FLUID_TYPE_PSEUDOPURE:
        case FLUID_TYPE_REFPROP:
            try {
                // Set the values in the output structure
                properties->p = state->p();
                properties->T = state->T();
                properties->d = state->rho();
                properties->h = state->h();
                properties->s = state->s();
                if (state->TwoPhase) {
                    properties->phase = 2;
                } else {
                    properties->phase = 1;
                }
                properties->cp = state->cp();
                properties->cv = state->cv();
                properties->a = state->speed_sound();
                if (state->TwoPhase && state->Q() >= 0 && state->Q() <= twophase_derivsmoothing_xend) {
                    // Use the smoothed derivatives between a quality of 0 and twophase_derivsmoothing_xend
                    properties->ddhp = state->drhodh_constp_smoothed(twophase_derivsmoothing_xend);  // [1/kPa -- > 1/Pa]
                    properties->ddph = state->drhodp_consth_smoothed(twophase_derivsmoothing_xend);  // [1/(kJ/kg) -- > 1/(J/kg)]
                } else if (state->TwoPhase && state->Q() >= 0 && state->Q() <= rho_smoothing_xend) {
                    // Use the smoothed density between a quality of 0 and rho_smoothing_xend
                    double rho_spline;
                    double dsplinedh;
                    double dsplinedp;
                    state->rho_smoothed(rho_smoothing_xend, rho_spline, dsplinedh, dsplinedp);
                    properties->ddhp = dsplinedh;
                    properties->ddph = dsplinedp;
                    properties->d = rho_spline;
                } else {
                    properties->ddhp = state->drhodh_constp();
                    properties->ddph = state->drhodp_consth();
                }
                properties->kappa = state->isothermal_compressibility();
                properties->beta = state->isobaric_expansion_coefficient();

                if (calc_transport) {
                    properties->eta = state->viscosity();
                    properties->lambda = state->conductivity();  //[kW/m/K --> W/m/K]
                } else {
                    properties->eta = -_HUGE;
                    properties->lambda = -_HUGE;
                }
            } catch (std::exception& e) {
                errorMessage((char*)e.what());
            }
            break;
        case FLUID_TYPE_INCOMPRESSIBLE_LIQUID:
        case FLUID_TYPE_INCOMPRESSIBLE_SOLUTION:
            try {
                // Set the values in the output structure
                properties->p = state->p();
                properties->T = state->T();
                properties->d = state->rho();
                properties->h = state->h();
                properties->s = state->s();
                properties->phase = 1;
                properties->cp = state->cp();
                properties->cv = state->cv();
                properties->a = -_HUGE;
                properties->ddhp = state->drhodh_constp();
                properties->ddph = 0.0;  // TODO: Fix this
                properties->kappa = -_HUGE;
                properties->beta = -_HUGE;
                if (calc_transport) {
                    properties->eta = state->viscosity();
                    properties->lambda = state->conductivity();  //[kW/m/K --> W/m/K]
                } else {
                    properties->eta = -_HUGE;
                    properties->lambda = -_HUGE;
                }
            } catch (std::exception& e) {
                errorMessage((char*)e.what());
            }
            break;
        default:
            errorMessage((char*)"Invalid fluid type!");
            break;
    }
}

void CoolPropSolver::setSat_p(double& p, ExternalSaturationProperties* const properties) {

    if (debug_level > 5) std::cout << format("setSat_p(%0.16e)\n", p);

    this->preStateChange();

    try {
        state->update(iP, p, iQ, 0);  // quality only matters for pseudo-pure fluids

        //! Saturation temperature
        properties->Tsat = state->TL();  // Not correct for pseudo-pure fluids
        //! Derivative of Ts wrt pressure
        properties->dTp = state->dTdp_along_sat();
        //! Derivative of dls wrt pressure
        properties->ddldp = state->drhodp_along_sat_liquid();
        //! Derivative of dvs wrt pressure
        properties->ddvdp = state->drhodp_along_sat_vapor();
        //! Derivative of hls wrt pressure
        properties->dhldp = state->dhdp_along_sat_liquid();
        //! Derivative of hvs wrt pressure
        properties->dhvdp = state->dhdp_along_sat_vapor();
        //! Density at bubble line (for pressure ps)
        properties->dl = state->rhoL();
        //! Density at dew line (for pressure ps)
        properties->dv = state->rhoV();
        //! Specific enthalpy at bubble line (for pressure ps)
        properties->hl = state->hL();
        //! Specific enthalpy at dew line (for pressure ps)
        properties->hv = state->hV();
        //! Saturation pressure
        properties->psat = p;
        //! Surface tension
        properties->sigma = state->surface_tension();
        //! Specific entropy at bubble line (for pressure ps)
        properties->sl = state->sL();
        //! Specific entropy at dew line (for pressure ps)
        properties->sv = state->sV();
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

void CoolPropSolver::setSat_T(double& T, ExternalSaturationProperties* const properties) {

    if (debug_level > 5) std::cout << format("setSat_T(%0.16e)\n", T);

    this->preStateChange();

    try {
        state->update(iT, T, iQ, 0);  // Quality only matters for pseudo-pure fluids

        properties->Tsat = T;
        properties->psat = state->p();
        properties->dl = state->rhoL();
        properties->dv = state->rhoV();
        properties->hl = state->hL();
        properties->hv = state->hV();
        properties->dTp = state->dTdp_along_sat();

        properties->ddldp = state->drhodp_along_sat_liquid();
        properties->ddvdp = state->drhodp_along_sat_vapor();
        properties->dhldp = state->dhdp_along_sat_liquid();
        properties->dhvdp = state->dhdp_along_sat_vapor();
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_ph(double& p, double& h, int& phase, ExternalThermodynamicState* const properties) {

    if (debug_level > 5) std::cout << format("setState_ph(p=%0.16e,h=%0.16e)\n", p, h);

    this->preStateChange();

    try {
        // Update the internal variables in the state instance
        state->update(iP, p, iH, h);

        if (!ValidNumber(state->rho()) || !ValidNumber(state->T())) {
            throw ValueError(format("p-h [%g, %g] failed for update", p, h));
        }

        // Set the values in the output structure
        this->postStateChange(properties);
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

void CoolPropSolver::setState_pT(double& p, double& T, ExternalThermodynamicState* const properties) {

    if (debug_level > 5) std::cout << format("setState_pT(p=%0.16e,T=%0.16e)\n", p, T);

    this->preStateChange();

    try {
        // Update the internal variables in the state instance
        state->update(iP, p, iT, T);

        // Set the values in the output structure
        this->postStateChange(properties);
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_dT(double& d, double& T, int& phase, ExternalThermodynamicState* const properties) {

    if (debug_level > 5) std::cout << format("setState_dT(d=%0.16e,T=%0.16e)\n", d, T);

    this->preStateChange();

    try {

        // Update the internal variables in the state instance
        state->update(iD, d, iT, T);

        // Set the values in the output structure
        this->postStateChange(properties);
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_ps(double& p, double& s, int& phase, ExternalThermodynamicState* const properties) {

    if (debug_level > 5) std::cout << format("setState_ps(p=%0.16e,s=%0.16e)\n", p, s);

    this->preStateChange();

    try {
        // Update the internal variables in the state instance
        state->update(iP, p, iS, s);

        // Set the values in the output structure
        this->postStateChange(properties);
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_hs(double& h, double& s, int& phase, ExternalThermodynamicState* const properties) {

    if (debug_level > 5) std::cout << format("setState_hs(h=%0.16e,s=%0.16e)\n", h, s);

    this->preStateChange();

    try {
        // Update the internal variables in the state instance
        state->update(iH, h, iS, s);

        // Set the values in the output structure
        this->postStateChange(properties);
    } catch (std::exception& e) {
        errorMessage((char*)e.what());
    }
}

double CoolPropSolver::Pr(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: Pr() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: Pr() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::T(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: T() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: T() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::a(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: a() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: a() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::beta(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: beta() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: beta() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::cp(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: cp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: cp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::cv(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: cv() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: cv() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::d(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: d() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: d() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::ddhp(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddhp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: ddhp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::ddph(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddph() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: ddph() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::eta(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: eta() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: eta() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::h(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: h() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: h() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::kappa(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: kappa() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: kappa() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::lambda(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: lambda() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: lambda() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::p(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: p() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: p() not implemented in the Solver object");
    return -_HUGE;
}

int CoolPropSolver::phase(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: phase() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: phase() not implemented in the Solver object");
    return -1;
}

double CoolPropSolver::s(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: s() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: s() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::d_der(ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: d_der() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: d_der() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::isentropicEnthalpy(double& p, ExternalThermodynamicState* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: isentropicEnthalpy() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: isentropicEnthalpy() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::dTp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dTp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: dTp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::ddldp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddldp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: ddldp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::ddvdp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: ddvdp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: ddvdp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::dhldp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dhldp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: dhldp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::dhvdp(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dhvdp() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: dhvdp() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::dl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dl() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: dl() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::dv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: dv() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: dv() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::hl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: hl() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: hl() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::hv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: hv() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: hv() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::sigma(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sigma() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: sigma() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::sl(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sl() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: sl() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::sv(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: sv() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: sv() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::psat(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: psat() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: psat() not implemented in the Solver object");
    return -_HUGE;
}

double CoolPropSolver::Tsat(ExternalSaturationProperties* const properties) {
    // Base function returns an error if called - should be redeclared by the solver object
    errorMessage((char*)"Internal error: Tsat() not implemented in the Solver object");
    //throw NotImplementedError((char*)"Internal error: Tsat() not implemented in the Solver object");
    return -_HUGE;
}
