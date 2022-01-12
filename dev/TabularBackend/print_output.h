#include "AbstractState.h"
#include "CoolProp.h"
#include "crossplatform_shared_ptr.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>

using namespace CoolProp;

const int LENGTH_PARAMETER = 20, LENGTH_HEOS = 15, LENGTH_BICUBIC = 15, LENGTH_TTSE = 15;
const int LENGTH_REL_ERROR_BICUBIC = 20, LENGTH_REL_ERROR_TTSE = 20;
const int LINE_LENGTH = LENGTH_PARAMETER + LENGTH_HEOS + LENGTH_BICUBIC + LENGTH_TTSE + LENGTH_REL_ERROR_BICUBIC + LENGTH_REL_ERROR_TTSE - 5;

double relative_error(double value, double expected){
    double rel_error;

    if (expected == 0.0 && value > 0.0){
        return _HUGE;
    }
    else if (expected == 0.0 && value < 0.0){
        return -HUGE;
    }
    else if (expected == 0.0 && value == 0.0){
        return 0.0;
    }
    else {
        rel_error = value / expected - 1.0;
        return rel_error;
    }
};

void print_parameter(std::string parameter, double heos, double bicubic, double ttse){

    double rel_error_bicubic, rel_error_ttse;

    rel_error_bicubic = relative_error(bicubic, heos);
    rel_error_ttse = relative_error(ttse, heos);


    std::cout << std::left
              << std::setw(LENGTH_PARAMETER) << parameter
              << std::setw(LENGTH_HEOS) << heos
              << std::setw(LENGTH_BICUBIC) << bicubic
              << std::setw(LENGTH_TTSE) << ttse
              << std::setw(LENGTH_REL_ERROR_BICUBIC) << rel_error_bicubic
              << std::setw(LENGTH_REL_ERROR_TTSE) << rel_error_ttse
              << std::endl;

};

void print_parameter(std::string parameter, int heos, int bicubic, int ttse){

    bool rel_error_bicubic=false, rel_error_ttse=false;

    if (heos == bicubic){
        rel_error_bicubic = true;
    }
    if (heos == ttse){
        rel_error_ttse = true;
    }


    std::cout << std::left << std::boolalpha
              << std::setw(LENGTH_PARAMETER) << parameter
              << std::setw(LENGTH_HEOS) << heos
              << std::setw(LENGTH_BICUBIC) << bicubic
              << std::setw(LENGTH_TTSE) << ttse
              << std::setw(LENGTH_REL_ERROR_BICUBIC) << rel_error_bicubic
              << std::setw(LENGTH_REL_ERROR_TTSE) << rel_error_ttse
              << std::endl;

};

void print_output(shared_ptr<AbstractState> HEOS, shared_ptr<AbstractState> Bicubic, shared_ptr<AbstractState> TTSE){
    std::cout << std::endl;
    std::cout << std::left
              << std::setw(LENGTH_PARAMETER) << "Parameter"
              << std::setw(LENGTH_HEOS) << "HEOS"
              << std::setw(LENGTH_BICUBIC) << "Bicubic"
              << std::setw(LENGTH_TTSE) << "TTSE"
              << std::setw(LENGTH_REL_ERROR_BICUBIC) << "rel. error Bicubic"
              << std::setw(LENGTH_REL_ERROR_TTSE) << "rel. error TTSE"
              << std::endl;
    std::cout << std::string(LINE_LENGTH, '-') << std::endl;

    print_parameter("umolar", HEOS->umolar(), Bicubic->umolar(), TTSE->umolar());
    print_parameter("rhomolar", HEOS->rhomolar(), Bicubic->rhomolar(), TTSE->rhomolar());
    std::cout << std::string(LINE_LENGTH, '-') << std::endl;

    print_parameter("T", HEOS->T(), Bicubic->T(), TTSE->T());
    print_parameter("p", HEOS->p(), Bicubic->p(), TTSE->p());
    print_parameter("hmolar", HEOS->hmolar(), Bicubic->hmolar(), TTSE->hmolar());
    print_parameter("smolar", HEOS->smolar(), Bicubic->smolar(), TTSE->smolar());
    print_parameter("cpmolar", HEOS->cpmolar(), Bicubic->cpmolar(), TTSE->cpmolar());
    print_parameter("cvmolar", HEOS->cvmolar(), Bicubic->cvmolar(), TTSE->cvmolar());
    print_parameter("viscosity", HEOS->viscosity(), Bicubic->viscosity(), TTSE->viscosity());
    print_parameter("conductivity", HEOS->conductivity(), Bicubic->conductivity(), TTSE->conductivity());
    print_parameter("phase", HEOS->phase(), Bicubic->phase(), TTSE->phase());
    if ( HEOS->phase() == iphase_twophase ){
        try {
            print_parameter("Q", HEOS->Q(), Bicubic->Q(), TTSE->Q());
        }
        catch (ValueError){
            std::cout << "Could not calculate Q for tabular data." << std::endl;
        }
    }
    if ( HEOS->phase() != iphase_twophase ){
        print_parameter("speed_sound", HEOS->speed_sound(), Bicubic->speed_sound(), TTSE->speed_sound());
    }
    else {
        std::cout << "speed_sound not defined for two phase state." << std::endl;
    }

    if ( HEOS->phase() != iphase_twophase ){
        print_parameter("dpdT|rhomolar", HEOS->first_partial_deriv(iP, iT, iDmolar), Bicubic->first_partial_deriv(iP, iT, iDmolar), TTSE->first_partial_deriv(iP, iT, iDmolar));
    }
    else {
        std::cout << "dpdT|rhomolar not defined for two phase state." << std::endl;
    }

    if ( HEOS->phase() != iphase_twophase ){
        print_parameter("d2pdT2|rhomolar", HEOS->second_partial_deriv(iP, iT, iDmolar, iT, iDmolar), Bicubic->second_partial_deriv(iP, iT, iDmolar, iT, iDmolar), TTSE->second_partial_deriv(iP, iT, iDmolar, iT, iDmolar));
        print_parameter("d2pdUmass2|rhomass", HEOS->second_partial_deriv(iP, iUmass, iDmass, iUmass, iDmass), Bicubic->second_partial_deriv(iP, iUmass, iDmass, iUmass, iDmass), TTSE->second_partial_deriv(iP, iUmass, iDmass, iUmass, iDmass));
        print_parameter("dUmassdDmass2|P", HEOS->second_partial_deriv(iUmass, iDmass, iUmass, iDmass, iP), Bicubic->second_partial_deriv(iUmass, iDmass, iUmass, iDmass, iP), TTSE->second_partial_deriv(iUmass, iDmass, iUmass, iDmass, iP));
    }
    else {
        std::cout << "dpdT|rhomolar not defined for two phase state." << std::endl;
    }

    if ( std::abs( HEOS->Q() ) < 1e-6 || std::abs( HEOS->Q() - 1.0 ) < 1e-6 ){
        try {
            print_parameter("dpdT|sat", HEOS->first_saturation_deriv(iP, iT), Bicubic->first_saturation_deriv(iP, iT), TTSE->first_saturation_deriv(iP, iT));
        }
        catch (ValueError) {
            std::cout << "Could not calculate dpdT|sat for tabular data." << std::endl;
        }
    }
    else {
        std::cout << "No outputs for saturated liquid or vapor because state is not saturated." << std::endl;
    }

    if ( HEOS->phase() == iphase_twophase ){
        try {
            print_parameter("drhodh|p, 2p", HEOS->first_two_phase_deriv(iDmolar, iHmolar, iP), Bicubic->first_two_phase_deriv(iP, iT, iDmolar), TTSE->first_two_phase_deriv(iP, iT, iDmolar));
        }
        catch (ValueError) {
            std::cout << "Could not calculate dpdT|p, 2p for tabular data." << std::endl;
        }
    }
    else {
        std::cout << "No outputs for two phase region because state is in single phase." << std::endl;
    }
}
