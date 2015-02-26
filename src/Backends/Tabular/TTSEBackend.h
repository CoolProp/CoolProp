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
        double evaluate_saturation(parameters output)
        {
            std::size_t iL = cached_saturation_iL, iV = cached_saturation_iV;
            double logp = log(_p);
            switch(output){
                case iT:
                {
                    double TV = CubicInterp(pure_saturation.logpV, pure_saturation.TV, iV-2, iV-1, iV, iV+1, logp);
                    double TL = CubicInterp(pure_saturation.logpL, pure_saturation.TL, iL-2, iL-1, iL, iL+1, logp);
                    return _Q*TV + (1-_Q)*TL;
                }
                case iDmolar:
                {
                    double rhoV = exp(CubicInterp(pure_saturation.logpV, pure_saturation.logrhomolarV, iV-2, iV-1, iV, iV+1, logp));
                    double rhoL = exp(CubicInterp(pure_saturation.logpL, pure_saturation.logrhomolarL, iL-2, iL-1, iL, iL+1, logp));
                    if (!ValidNumber(rhoV)){throw ValueError("rhoV is invalid");}
                    if (!ValidNumber(rhoL)){throw ValueError("rhoL is invalid");}
                    return 1/(_Q/rhoV + (1-_Q)/rhoL);
                }
                default:
                    throw ValueError("can't be something other than T or rho");
            }
        }
        long double calc_T(void){
            if (using_single_phase_table){
                return evaluate_single_phase_hmolarp(iT, cached_single_phase_i, cached_single_phase_j);
            }
            else{
                return evaluate_saturation(iT);
            }
        }
        long double calc_rhomolar(void){
            if (using_single_phase_table){
                return evaluate_single_phase_hmolarp(iDmolar, cached_single_phase_i, cached_single_phase_j);
            }
            else{
                return evaluate_saturation(iDmolar);
            }
        }
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
