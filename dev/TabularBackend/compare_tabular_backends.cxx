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

	std::filesystem::path cp_home = get_global_param_string("HOME");
	cp_home /= ".CoolProp";
	cp_home /= "Tables";

	if (std::filesystem::is_directory(cp_home)){
		std::cout << cp_home << std::endl;
		std::filesystem::remove_all(cp_home);
	}

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> _umolar(-50.0, 100.0);
    std::uniform_real_distribution<> _log_rhomolar(-2, 5);

    double rhomolar = 1000.0, umolar = 2500.0;
    std::string fluid = "helium";

    // Low level interface
    shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS",fluid));
    shared_ptr<AbstractState> Bicubic(AbstractState::factory("BICUBIC&HEOS",fluid));
    shared_ptr<AbstractState> TTSE(AbstractState::factory("TTSE&HEOS",fluid));

    for (int i=0; i<=1000; i++){
    	rhomolar = std::pow(10.0, _log_rhomolar(gen));
    	umolar = _umolar(gen);

    	try {
			HEOS->update(DmolarUmolar_INPUTS, rhomolar, umolar);
			Bicubic->update(DmolarUmolar_INPUTS, rhomolar, umolar);
			TTSE->update(DmolarUmolar_INPUTS, rhomolar, umolar);

			std::cout << "sucess: " << i << std::endl;
    	}
    	catch (CoolPropBaseError)
    	{
    	}
    }

    print_output(HEOS, Bicubic, TTSE);

    // All done return
    return EXIT_SUCCESS;
}
