#ifndef TTSEBACKEND_H
#define TTSEBACKEND_H

#include "TabularBackends.h"

namespace CoolProp
{

class TTSEBackend : public TabularBackend
{
    public:
        std::string backend_name(void){return "TTSEBackend";}
        /// Instantiator; base class loads or makes tables
        TTSEBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {}
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
        double evaluate_single_phase(SinglePhaseGriddedTableData &table, parameters output, std::size_t i, std::size_t j);
        double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase(single_phase_logph, output, i, j);
        }
        double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase(single_phase_logpT, output, i, j);
        }
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
