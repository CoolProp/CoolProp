#ifndef BICUBICBACKEND_H
#define BICUBICBACKEND_H

#include "TabularBackends.h"
#include "Exceptions.h"
#include "Eigen/Core"


namespace CoolProp
{
    
/// This structure holds the coefficients for one cell, the coefficients are stored in matrices
/// and can be obtained by the get() function.  
class CellCoeffs{
    private:
        std::size_t alt_i, alt_j;
        bool _valid, _has_valid_neighbor;
    public:
        double dx_dxhat, dy_dyhat;    
        CellCoeffs(){_valid = false; _has_valid_neighbor = false;} 
        std::vector<double> T, rhomolar, hmolar, p, smolar, umolar;
        /// Return a const reference to the desired matrix
        const std::vector<double> & get(parameters params){
            switch(params){
                case iT: return T;
                case iP: return p;
                case iDmolar: return rhomolar;
                case iHmolar: return hmolar;
                case iSmolar: return smolar;
				case iUmolar: return umolar;
                default: throw KeyError(format("Invalid key to get() function of CellCoeffs"));
            }
        };
        /// Set one of the matrices in this class
        void set(parameters params, const std::vector<double> &mat){
            switch(params){
				case iT: T = mat; break;
                case iP: p = mat; break;
                case iDmolar: rhomolar = mat; break;
                case iHmolar: hmolar = mat; break;
                case iSmolar: smolar = mat; break;
				case iUmolar: umolar = mat; break;
                default: throw KeyError(format("Invalid key to set() function of CellCoeffs"));
            }
        };
        /// Returns true if the cell coefficients seem to have been calculated properly
        bool valid(){return _valid;};
        /// Call this function to set the valid flag to true
        void set_valid(){_valid = true;};
        /// Call this function to set the valid flag to false
        void set_invalid(){_valid = false;};
        /// Set the neighboring (alternate) cell to be used if the cell is invalid
        void set_alternate(std::size_t i, std::size_t j){alt_i = i; alt_j = j; _has_valid_neighbor = true;}
        /// Get neighboring(alternate) cell to be used if this cell is invalid
        void get_alternate(std::size_t &i, std::size_t &j){
            if (_has_valid_neighbor){
                i = alt_i; j = alt_j;
            }
            else{
                throw ValueError("No valid neighbor");
            }
        }
        /// Returns true if cell is invalid and it has valid neighbor
        bool has_valid_neighbor(){
            return _has_valid_neighbor;
        }
};

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
typedef std::vector<std::vector<double> > mat;
class BicubicBackend : public TabularBackend
{
    protected:
        std::vector<std::vector<CellCoeffs> > coeffs_ph, coeffs_pT;
    public:
        /// Instantiator; base class loads or makes tables
		BicubicBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS){
            // If a pure fluid, don't need to set fractions, go ahead and build
            if (this->AS->get_mole_fractions().size() == 1){
                check_tables();
			    build_coeffs(single_phase_logph, coeffs_ph);
			    build_coeffs(single_phase_logpT, coeffs_pT);
            }
		};
        std::string backend_name(void){return "BicubicBackend";}
        /// Build the \f$a_{i,j}\f$ coefficients for bicubic interpolation
        void build_coeffs(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs);
        /** Update the state
         */
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
        
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
        double evaluate_single_phase_derivative(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs, parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny);
		double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny){
            return evaluate_single_phase_derivative(single_phase_logph, coeffs_ph, output, _hmolar, _p, i, j, Nx, Ny);
        };
        double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny){
            return evaluate_single_phase_derivative(single_phase_logpT, coeffs_pT, output, _T, _p, i, j, Nx, Ny);
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
		double evaluate_single_phase(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs, parameters output, double x, double y, std::size_t i, std::size_t j);
        double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j){
			return evaluate_single_phase(single_phase_logph, coeffs_ph, output, _hmolar, _p, i, j);
		};
        double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j){
			return evaluate_single_phase(single_phase_logpT, coeffs_pT, output, _T, _p, i, j);
		};
        
        /**
         * @brief Evaluate the single-phase transport properties using linear interpolation.  Works well except for near the critical point
         * @param table A reference to the table to be used
         * @param output The output parameter, viscosity or conductivity
         * @param x The 
         * @param y
         * @return 
         */
        double evaluate_single_phase_transport(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j);
        
        double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase_transport(single_phase_logph, output, _hmolar, _p, i, j);
        };
		double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j){
            return evaluate_single_phase_transport(single_phase_logpT, output, _T, _p, i, j);
        };
};

}

#endif // BICUBICBACKEND_H
