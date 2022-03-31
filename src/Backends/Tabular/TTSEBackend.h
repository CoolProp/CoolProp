#ifndef TTSEBACKEND_H
#define TTSEBACKEND_H

#include "TabularBackends.h"
#include "DataStructures.h"

namespace CoolProp {

class TTSEBackend : public TabularBackend
{
   public:
    std::string backend_name(void) {
        return get_backend_string(TTSE_BACKEND);
    }
    /// Instantiator; base class loads or makes tables
    TTSEBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend(AS) {
        imposed_phase_index = iphase_not_imposed;
        // If a pure fluid or a predefined mixture, don't need to set fractions, go ahead and build
        if (!this->AS->get_mole_fractions().empty()) {
            check_tables();
            SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
            SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
            dataset->build_coeffs(single_phase_logph, dataset->coeffs_ph);
            dataset->build_coeffs(single_phase_logpT, dataset->coeffs_pT);
            is_mixture = (this->AS->get_mole_fractions().size() > 1);
        }
    }
    double evaluate_single_phase(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j);
    double evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j);
    double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) {
        SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
        return evaluate_single_phase(single_phase_logph, output, _hmolar, _p, i, j);
    }
    double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) {
        SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
        return evaluate_single_phase(single_phase_logpT, output, _T, _p, i, j);
    }
    double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j) {
        SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
        return evaluate_single_phase_transport(single_phase_logph, output, _hmolar, _p, i, j);
    }
    double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j) {
        SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
        return evaluate_single_phase_transport(single_phase_logpT, output, _T, _p, i, j);
    }
    void invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                               double x, double y, std::size_t i, std::size_t j);
    void invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                               double y, double x, std::size_t i, std::size_t j);

    /// Find the best set of i,j for native inputs.
    virtual void find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, double x,
                                                  double y, std::size_t& i, std::size_t& j) {
        return table.find_native_nearest_good_neighbor(x, y, i, j);
    };
    /// Ask the derived class to find the nearest neighbor (pure virtual)
    virtual void find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                       const parameters variable1, const double value1, const parameters otherkey, const double otherval,
                                       std::size_t& i, std::size_t& j) {
        table.find_nearest_neighbor(variable1, value1, otherkey, otherval, cached_single_phase_i, cached_single_phase_j);
    };

    /**
         * @brief Evaluate a derivative in terms of the native inputs of the table
         * @param table A reference to the table to be used
         * @param output The output variable
         * @param x The
         * @param y
         * @param i
         * @param j
         * @param Nx The number of derivatives with respect to x with y held constant
         * @param Ny The number of derivatives with respect to y with x held constant
         * @return
         */
    double evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j,
                                            std::size_t Nx, std::size_t Ny);
    double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
        return evaluate_single_phase_derivative(single_phase_logph, output, _hmolar, _p, i, j, Nx, Ny);
    };
    double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
        return evaluate_single_phase_derivative(single_phase_logpT, output, _T, _p, i, j, Nx, Ny);
    };
};

}  // namespace CoolProp

#endif  // TTSEBACKEND_H
