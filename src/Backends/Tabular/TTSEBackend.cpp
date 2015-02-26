#include "TTSEBackend.h"
#include "CoolProp.h"

void CoolProp::TTSEBackend::update(CoolProp::input_pairs input_pair, double val1, double val2)
{
    // Flush the cached indices
    cached_single_phase_i = -1; cached_single_phase_j = -1;
    
    switch(input_pair){
        case HmolarP_INPUTS:{
            _hmolar = val1; _p = val2;
            if (!single_phase_logph.native_inputs_are_in_range(_hmolar, _p)){
                // Use the AbstractState instance
                using_single_phase_table = false;
                if (get_debug_level() > 5){ std::cout << "inputs are not in range"; }
                throw ValueError(format("inputs are not in range, hmolar=%Lg, p=%Lg", static_cast<CoolPropDbl>(_hmolar), _p));
            }
            else{
                using_single_phase_table = true; // Use the table !
                std::size_t iL, iV;
                CoolPropDbl hL = 0, hV = 0;
                if (pure_saturation.is_inside(_p, iHmolar, _hmolar, iL, iV, hL, hV)){
                    using_single_phase_table = false;
                    _Q = (static_cast<double>(_hmolar)-hL)/(hV-hL);
                    if(!is_in_closed_range(0.0,1.0,static_cast<double>(_Q))){
                        throw ValueError("vapor quality is not in (0,1)");
                    }
                    else{
                        cached_saturation_iL = iL;
                        cached_saturation_iV = iV;
                    }
                }
                else{
                    // Find and cache the indices i, j
                    single_phase_logph.find_native_nearest_neighbor(_hmolar, _p, cached_single_phase_i, cached_single_phase_j);
                }
                // Make pointer to the table in use
            }
            break;
        }
        case HmassP_INPUTS:
        {
            // Call again, but this time with molar units
            // H: [J/kg] * [kg/mol] -> [J/mol]
            update(HmolarP_INPUTS, val1 * AS->molar_mass(), val2);
            return;
        }
        case DmolarP_INPUTS:
        case PUmolar_INPUTS:
        case PSmolar_INPUTS:
        case DmolarHmolar_INPUTS:
        case HmolarSmolar_INPUTS:
            throw ValueError("To be implemented as a 1-D iteration using PH table");
        case PT_INPUTS:
            throw ValueError("To be implemented using PT table");
        case SmolarT_INPUTS:
        case DmolarT_INPUTS:
            throw ValueError("To be implemented as a 1-D iteration using PT table");
        default:
            throw ValueError();
    }
}

/// Use the hmolar-log(p) table to evaluate an output
double CoolProp::TTSEBackend::evaluate_single_phase_hmolarp(parameters output, std::size_t i, std::size_t j)
{
    // Define pointers for the matrices to be used; 
    std::vector<std::vector<double> > *y, *dydh, *dydp, *d2ydhdp, *d2ydp2, *d2ydh2;
    
    // Connect the pointers based on the output variable desired
    switch(output){
        case iT:
            y = &single_phase_logph.T; dydh = &single_phase_logph.dTdx; dydp = &single_phase_logph.dTdy;
            d2ydhdp = &single_phase_logph.d2Tdxdy; d2ydh2 = &single_phase_logph.d2Tdx2; d2ydp2 = &single_phase_logph.d2Tdy2;
            break;
        case iDmolar:
            y = &single_phase_logph.rhomolar; dydh = &single_phase_logph.drhomolardx; dydp = &single_phase_logph.drhomolardy;
            d2ydhdp = &single_phase_logph.d2rhomolardxdy; d2ydh2 = &single_phase_logph.d2rhomolardx2; d2ydp2 = &single_phase_logph.d2rhomolardy2;
            break;
        case iSmolar:
            y = &single_phase_logph.smolar; dydh = &single_phase_logph.dsmolardx; dydp = &single_phase_logph.dsmolardy;
            d2ydhdp = &single_phase_logph.d2smolardxdy; d2ydh2 = &single_phase_logph.d2smolardx2; d2ydp2 = &single_phase_logph.d2smolardy2;
            break;
        //case iUmolar:
        default:
            throw ValueError();
    }
    
    // Distances from the node
	double deltah = static_cast<double>(_hmolar) - single_phase_logph.xvec[i];
    double deltap = static_cast<double>(_p) - single_phase_logph.yvec[j];
    
    // Calculate the output value desired
    double val = (*y)[i][j]+deltah*(*dydh)[i][j]+deltap*(*dydp)[i][j]+0.5*deltah*deltah*(*d2ydh2)[i][j]+0.5*deltap*deltap*(*d2ydp2)[i][j]+deltap*deltah*(*d2ydhdp)[i][j];
    
    // Cache the output value calculated
    switch(output){
        case iT:  _T = val; break;
        case iDmolar: _rhomolar = val; break;
        case iSmolar: _smolar = val; break;
        //case iUmolar:
        default: throw ValueError();
    }
    return val;
}
