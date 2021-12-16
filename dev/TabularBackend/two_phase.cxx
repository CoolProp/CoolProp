#include "AbstractState.h"
#include "CoolProp.h"
#include "crossplatform_shared_ptr.h"
#include "Backends/Tabular/TabularBackends.h"
#include "print_output.h"

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

    HEOS->update(QT_INPUTS, 0.5, HEOS->Tmin() + 0.25*( HEOS->T_critical() - HEOS->Tmin() ) );

    Bicubic->update(DmolarUmolar_INPUTS, HEOS->rhomolar(), HEOS->umolar());
    TTSE->update(DmolarUmolar_INPUTS, HEOS->rhomolar(), HEOS->umolar());

    print_output(HEOS, Bicubic, TTSE);

    // All done return
    return EXIT_SUCCESS;
}
