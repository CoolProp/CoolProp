#ifndef TTSEBACKEND_H
#define TTSEBACKEND_H

#include "TabularBackends.h"

namespace CoolProp
{
    
class BicubicBackend : public TabularBackend
{
    public:
        std::string backend_name(void){return "BicubicBackend";}
        BicubicBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {}
};

class TTSEBackend : public TabularBackend
{
    public:
        std::string backend_name(void){return "TTSEBackend";}
        TTSEBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {}
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
        double evaluate_single_phase_hmolarp(parameters output, std::size_t i, std::size_t j);
        long double calc_T(void){
            if (using_single_phase_table){
                return evaluate_single_phase_hmolarp(iT, cached_single_phase_i, cached_single_phase_j);
            }
            else{
                return pure_saturation.evaluate(iT, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        long double calc_rhomolar(void){
            if (using_single_phase_table){
                return evaluate_single_phase_hmolarp(iDmolar, cached_single_phase_i, cached_single_phase_j);
            }
            else{
                return pure_saturation.evaluate(iDmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
