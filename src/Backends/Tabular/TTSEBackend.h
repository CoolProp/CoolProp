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
        enum selected_table_options{SELECTED_NO_TABLE=0, SELECTED_PH_TABLE, SELECTED_PT_TABLE};
        selected_table_options selected_table;
        std::string backend_name(void){return "TTSEBackend";}
        
        TTSEBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {}
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
        double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j);
        double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j);
        CoolPropDbl calc_T(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iT, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return _T;
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
            }
            else{
                return pure_saturation.evaluate(iT, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_rhomolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iDmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iDmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
            }
            else{
                return pure_saturation.evaluate(iDmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_hmolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return _hmolar;
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iHmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
            }
            else{
                return pure_saturation.evaluate(iHmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_smolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iSmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iSmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
            }
            else{
                return pure_saturation.evaluate(iHmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
