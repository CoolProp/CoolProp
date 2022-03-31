#ifndef BICUBICBACKEND_H
#define BICUBICBACKEND_H

#include "TabularBackends.h"
#include "Exceptions.h"
#include "DataStructures.h"
#include "Eigen/Core"

namespace CoolProp {

/** \brief This class implements bicubic interpolation, as very clearly laid out by
 * the page on wikipedia: http://en.wikipedia.org/wiki/Bicubic_interpolation
 *
 * Essentially you have an already-inverted matrix that you need to multiply

\f[
A^{-1} =  \left[ \begin{array}{*{16}c} 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 \\
 -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 \\
 9 & -9 & -9 & 9 & 6 & 3 & -6 & -3 & 6 & -6 & 3 & -3 & 4 & 2 & 2 & 1 \\
 -6 & 6 & 6 & -6 & -3 & -3 & 3 & 3 & -4 & 4 & -2 & 2 & -2 & -2 & -1 & -1 \\
 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 \\
 -6 & 6 & 6 & -6 & -4 & -2 & 4 & 2 & -3 & 3 & -3 & 3 & -2 & -1 & -2 & -1 \\
 4 & -4 & -4 & 4 & 2 & 2 & -2 & -2 & 2 & -2 & 2 & -2 & 1 & 1 & 1 & 1 \end{array} \right]
 \f]

 \f[
 x = \frac{h-h_i}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial x}{\partial h} = \frac{1}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial h}{\partial x} = h_{i+1}-h_i
 \f]

 \f[
 y = \frac{p-p_j}{p_{j+1}-p_j}
 \f]
 \f[
 \frac{\partial y}{\partial p} = \frac{1}{p_{j+1}-p_j}
 \f]
  \f[
 \frac{\partial p}{\partial y} = p_{j+1}-p_j
 \f]

 \f[
 \frac{\partial f}{\partial x} = \frac{\partial f}{\partial h}\cdot\frac{\partial h}{\partial x}
 \f]
 \
*/
typedef std::vector<std::vector<double>> mat;
class BicubicBackend : public TabularBackend
{
   public:
    /// Instantiator; base class loads or makes tables
    BicubicBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend(AS) {
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
    };
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        this->AS->set_mole_fractions(mole_fractions);
        is_mixture = true;
        // Check the tables and build if necessary
        check_tables();
        // For mixtures, the construction of the coefficients is delayed until this
        // function so that the set_mole_fractions function can be called
        SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
        SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
        dataset->build_coeffs(single_phase_logph, dataset->coeffs_ph);
        dataset->build_coeffs(single_phase_logpT, dataset->coeffs_pT);
    };
    std::string backend_name(void) {
        return get_backend_string(BICUBIC_BACKEND);
    }

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
    double evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, parameters output,
                                            double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);
    double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logph, dataset->coeffs_ph, output, _hmolar, _p, i, j, Nx, Ny);
    };
    double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) {
        return evaluate_single_phase_derivative(dataset->single_phase_logpT, dataset->coeffs_pT, output, _T, _p, i, j, Nx, Ny);
    };

    /**
         * @brief
         * @param table A reference to the table that is to be used
         * @param coeffs A reference to the matrix of bicubic coefficients
         * @param output What output is desired
         * @param x The x value for the native inputs
         * @param y
         * @param i
         * @param j
         * @return
         */
    double evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                 const parameters output, const double x, const double y, const std::size_t i, const std::size_t j);
    double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logph, dataset->coeffs_ph, output, _hmolar, _p, i, j);
    };
    double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase(dataset->single_phase_logpT, dataset->coeffs_pT, output, _T, _p, i, j);
    };

    virtual void find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, double x,
                                                  double y, std::size_t& i, std::size_t& j);

    /// Ask the derived class to find the nearest neighbor (pure virtual)
    virtual void find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                       const parameters variable1, const double value1, const parameters otherkey, const double otherval,
                                       std::size_t& i, std::size_t& j);

    /**
         * @brief Evaluate the single-phase transport properties using linear interpolation.  Works well except for near the critical point
         * @param table A reference to the table to be used
         * @param output The output parameter, viscosity or conductivity
         * @param x The
         * @param y
         * @return
         */
    double evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i, std::size_t j);

    double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase_transport(dataset->single_phase_logph, output, _hmolar, _p, i, j);
    };
    double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j) {
        return evaluate_single_phase_transport(dataset->single_phase_logpT, output, _T, _p, i, j);
    };

    /**
         * @brief Use the table to solve for the x variable of the table given the y coordinate of the table and a variable that can yield a unique solution for x
         * @param table The table to be used
         * @param coeffs The matrix of coefficients to be used
         * @param other_key The x variable
         * @param other The value of the x-ish variable to be used to find d
         * @param i The x-coordinate of the cell
         * @param j The y-coordinate of the cell
         */
    void invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters other_key,
                               double other, double y, std::size_t i, std::size_t j);
    void invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, parameters other_key,
                               double other, double x, std::size_t i, std::size_t j);
};

}  // namespace CoolProp

#endif  // BICUBICBACKEND_H
