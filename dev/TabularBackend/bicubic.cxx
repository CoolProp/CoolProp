#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"
#include "Backends/Tabular/TabularBackends.h"
using namespace CoolProp;
int main()
{
	set_debug_level(0);

    // Low level interface
    shared_ptr<AbstractState> Bicubic(AbstractState::factory("BICUBIC&HEOS","n-Heptane"));
    Bicubic->update(DmassT_INPUTS, 4.2692360381933803, 372.49600373915791); // SI units

    shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS","n-Heptane"));
    HEOS->update(DmassT_INPUTS, 4.2692360381933803, 372.49600373915791); // SI units

    std::cout << Bicubic->p() << std::endl;
    std::cout << HEOS->p() << std::endl;

    shared_ptr<AbstractState> AS(AbstractState::factory("HEOS","n-Heptane"));

    LogDUTable single_phase_logdu;

    single_phase_logdu.AS=AS;
    single_phase_logdu.set_limits();
    single_phase_logdu.build(AS);

    std::cout << single_phase_logdu.dTdx[100][100] << std::endl;

    TabularDataSet data_set;
    data_set.build_tables(AS);

    std::cout << "Bicubic" << std::endl;

    std::cout << Bicubic->rhomolar() << std::endl;

    Bicubic->update(DmolarUmolar_INPUTS, 2.7635029493162853e-07, 13886.488232581964);

    std::cout << "rhomolar = " << Bicubic->rhomolar() << ", umolar = " << Bicubic->umolar() << std::endl;


    std::cout << "p = " << Bicubic->p() << ", T = " << Bicubic->T() << std::endl;
    std::cout << "hmolar = " << Bicubic->hmolar() << ", smolar = " << Bicubic->smolar() << std::endl;

    /*

    Bicubic->update(DmassUmass_INPUTS, Bicubic->rhomass(), 2.0*Bicubic->umass());

    std::cout << Bicubic->p() << "\t" << Bicubic->T() << "\t" << Bicubic->hmass() << std::endl;

	*/

    std::cout << "Done" << std::endl;

    // All done return
    return EXIT_SUCCESS;
}
