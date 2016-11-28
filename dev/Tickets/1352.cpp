
#include "crossplatform_shared_ptr.h"
#include <AbstractState.h>
#include <CoolProp.h>


int main(int argc, const char* argv[]) {
    shared_ptr<CoolProp::AbstractState> pState;
    pState.reset(CoolProp::AbstractState::factory("HEOS", "Water"));
	double T_test = 25 + 273.15;
	pState->update(CoolProp::QT_INPUTS, 0.3, T_test);
	double rho_test = pState->rhomass();

	pState->specify_phase(CoolProp::iphase_not_imposed);
	pState->update(CoolProp::DmassT_INPUTS, rho_test, T_test);

	pState->specify_phase(CoolProp::iphase_twophase);
	pState->update(CoolProp::DmassT_INPUTS, rho_test, T_test);

}