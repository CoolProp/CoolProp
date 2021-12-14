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

int main()
{
	set_debug_level(0);

    double rhomolar = 1000.0, umolar = 2500.0;
    std::string fluid = "helium";

    // Low level interface
    shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS",fluid));
    shared_ptr<AbstractState> Bicubic(AbstractState::factory("BICUBIC&HEOS",fluid));
    shared_ptr<AbstractState> TTSE(AbstractState::factory("TTSE&HEOS",fluid));

    HEOS->update(DmolarUmolar_INPUTS, rhomolar, umolar);
    Bicubic->update(DmolarUmolar_INPUTS, rhomolar, umolar);
    TTSE->update(DmolarUmolar_INPUTS, rhomolar, umolar);

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
    print_parameter("speed_sound", HEOS->speed_sound(), Bicubic->speed_sound(), TTSE->speed_sound());
    print_parameter("dpdT|rhomolar", HEOS->first_partial_deriv(iP, iT, iDmolar), Bicubic->first_partial_deriv(iP, iT, iDmolar), TTSE->first_partial_deriv(iP, iT, iDmolar));
    // print_parameter("dpdT|sat", HEOS->first_saturation_deriv(iP, iT), Bicubic->first_saturation_deriv(iP, iT), 0.0);
    // print_parameter("drhodh|p, 2p", HEOS->first_two_phase_deriv(iDmolar, iHmolar, iP), Bicubic->first_two_phase_deriv(iP, iT, iDmolar), 0.0);
    // print_parameter("drhodh|p, 2p", HEOS->first_two_phase_deriv_splined(iDmolar, iHmolar, iP, 0.0), Bicubic->first_two_phase_deriv_splined(iP, iT, iDmolar, 0.0), 0.0);
    print_parameter("phase", HEOS->phase(), Bicubic->phase(), TTSE->phase());

    // All done return
    return EXIT_SUCCESS;
}
