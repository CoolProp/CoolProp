#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"
#include "Backends/Tabular/TabularBackends.h"
#include "iomanip"
#include <string>
using namespace CoolProp;

const int LENGTH_PARAMETER = 15, LENGTH_HEOS = 15, LENGTH_BICUBIC = 15, LENGTH_TTSE = 15;
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


    std::cout << std::setw(LENGTH_PARAMETER) << std::left << parameter
    		  << std::setw(LENGTH_HEOS) << std::left << heos
			  << std::setw(LENGTH_BICUBIC) << std::left << bicubic
			  << std::setw(LENGTH_TTSE) << std::left << ttse
			  << std::setw(LENGTH_REL_ERROR_BICUBIC) << std::left << rel_error_bicubic
			  << std::setw(LENGTH_REL_ERROR_TTSE) << std::left << rel_error_ttse
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


    std::cout << std::setw(LENGTH_PARAMETER) << std::left << parameter
    		  << std::setw(LENGTH_HEOS) << std::left << heos
			  << std::setw(LENGTH_BICUBIC) << std::left << bicubic
			  << std::setw(LENGTH_TTSE) << std::left << ttse
			  << std::setw(LENGTH_REL_ERROR_BICUBIC) << std::left << rel_error_bicubic
			  << std::setw(LENGTH_REL_ERROR_TTSE) << std::left << rel_error_ttse
			  << std::endl;

};

int main()
{
	set_debug_level(0);

    double rhomolar = 2.7635029493162853e-07, umolar = 13886.488232581964;
    std::string fluid = "helium";

    // Low level interface
    shared_ptr<AbstractState> Bicubic(AbstractState::factory("BICUBIC&HEOS",fluid));
    shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS",fluid));


    Bicubic->update(DmolarUmolar_INPUTS, rhomolar, umolar);
    HEOS->update(DmolarUmolar_INPUTS, rhomolar, umolar);

    std::cout << std::endl;
    std::cout << std::setw(LENGTH_PARAMETER) << std::left << "Parameter"
    		  << std::setw(LENGTH_HEOS) << std::left << "HEOS"
    		  << std::setw(LENGTH_BICUBIC) << std::left << "Bicubic"
			  << std::setw(LENGTH_TTSE) << std::left << "TTSE"
    		  << std::setw(LENGTH_REL_ERROR_BICUBIC) << std::left << "rel. error Bicubic"
			  << std::setw(LENGTH_REL_ERROR_TTSE) << std::left << "rel. error TTSE"
			  << std::endl;
    std::cout << std::string(LINE_LENGTH, '-') << std::endl;

    print_parameter("umolar", HEOS->umolar(), Bicubic->umolar(), 0.0);
    print_parameter("rhomolar", HEOS->rhomolar(), Bicubic->rhomolar(), 0.0);
    std::cout << std::string(LINE_LENGTH, '-') << std::endl;

    print_parameter("T", HEOS->T(), Bicubic->T(), 0.0);
    print_parameter("p", HEOS->p(), Bicubic->p(), 0.0);
    print_parameter("hmolar", HEOS->hmolar(), Bicubic->hmolar(), 0.0);
    print_parameter("smolar", HEOS->smolar(), Bicubic->smolar(), 0.0);
    print_parameter("cpmolar", HEOS->cpmolar(), Bicubic->cpmolar(), 0.0);
    print_parameter("cvmolar", HEOS->cvmolar(), Bicubic->cvmolar(), 0.0);
    print_parameter("viscosity", HEOS->viscosity(), Bicubic->viscosity(), 0.0);
    print_parameter("conductivity", HEOS->conductivity(), Bicubic->conductivity(), 0.0);
    print_parameter("speed_sound", HEOS->speed_sound(), Bicubic->speed_sound(), 0.0);
    print_parameter("dpdT|rhomolar", HEOS->first_partial_deriv(iP, iT, iDmolar), Bicubic->first_partial_deriv(iP, iT, iDmolar), 0.0);
    // print_parameter("dpdT|sat", HEOS->first_saturation_deriv(iP, iT), Bicubic->first_saturation_deriv(iP, iT), 0.0);
    // print_parameter("drhodh|p, 2p", HEOS->first_two_phase_deriv(iDmolar, iHmolar, iP), Bicubic->first_two_phase_deriv(iP, iT, iDmolar), 0.0);
    // print_parameter("drhodh|p, 2p", HEOS->first_two_phase_deriv_splined(iDmolar, iHmolar, iP, 0.0), Bicubic->first_two_phase_deriv_splined(iP, iT, iDmolar, 0.0), 0.0);
    print_parameter("phase", HEOS->phase(), Bicubic->phase(), 0);

    // All done return
    return EXIT_SUCCESS;
}
