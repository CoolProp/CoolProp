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
        TTSEBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {
            // If a pure fluid, don't need to set fractions, go ahead and build
            if (this->AS->get_mole_fractions().size() == 1){
                check_tables();
            }
        }
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
        double evaluate_single_phase(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j);
        double evaluate_single_phase_transport(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j);
        double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase(single_phase_logph, output, _hmolar, _p, i, j);
        }
        double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase(single_phase_logpT, output, _T, _p, i, j);
        }
        double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase_transport(single_phase_logph, output, _hmolar, _p, i, j);
        }
        double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase_transport(single_phase_logpT, output, _T, _p, i, j);
        }
        double invert_single_phase_x(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j);
        double invert_single_phase_y(SinglePhaseGriddedTableData &table, parameters output, double y, double x, std::size_t i, std::size_t j);
        
        /**
         * @brief Evaluate a derivative in terms of the native inputs of the table
         * @param table A reference to the table to be used
         * @param coeffs A reference to the matrix of the coefficients
         * @param output The output variable
         * @param x The 
         * @param y
         * @param i
         * @param j
         * @param Nx The number of derivatives with respect to x with y held constant
         * @param Ny The number of derivatives with respect to y with x held constant
         * @return 
         */
        double evaluate_single_phase_derivative(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);
        double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny){
            return evaluate_single_phase_derivative(single_phase_logph, output, _hmolar, _p, i, j, Nx, Ny);
        };
        double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny){
            return evaluate_single_phase_derivative(single_phase_logpT, output, _T, _p, i, j, Nx, Ny);
        };
};

} // namespace CoolProp

#endif // TTSEBACKEND_H
