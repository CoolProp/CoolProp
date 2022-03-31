#include "AbstractState.h"
#include "CoolProp.h"
#include "crossplatform_shared_ptr.h"
#include "print_output.h"
#include <random>
#include <filesystem>
using namespace CoolProp;

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

    print_output(HEOS, Bicubic, TTSE);

    // All done return
    return EXIT_SUCCESS;
}
