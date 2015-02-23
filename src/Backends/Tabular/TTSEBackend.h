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
        double evaluate_hmolarp(parameters output, std::size_t i, std::size_t j);
        long double calc_T(void){return evaluate_hmolarp(iT, cached_i, cached_j);};
        long double calc_rhomolar(void){return evaluate_hmolarp(iDmolar, cached_i, cached_j);};
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
