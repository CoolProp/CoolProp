
#include "crossplatform_shared_ptr.h"
#include <AbstractState.h>
#include <CoolProp.h>

int main(int argc, const char* argv[]) {
    shared_ptr<CoolProp::AbstractState> pState;
    pState.reset(CoolProp::AbstractState::factory("TTSE&HEOS", "Water"));
    std::cout << "T crit: " << pState->T_critical() << std::endl;
    pState->update(CoolProp::QT_INPUTS, 0.2, 373.15);
    double res = pState->first_two_phase_deriv_splined(CoolProp::iDmass, CoolProp::iHmass, CoolProp::iP, 0.3);

    /*x, y1 = [], []
		for Q in np.linspace(0, 0.3, steps) :
			AS.update(CoolProp.PQ_INPUTS, 101325, Q)
			x.append(AS.Q())
			try :
			y1.append(AS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3))
			except Exception as e :
	print e
		y1.append(np.NAN)
		break
		plt.plot(x, y1, label = 'Two-phase (splined, tabular)', ls = '--', lw = 3)

		´*/
}