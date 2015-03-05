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
			build_coeffs(single_phase_logph, coeffs_ph);
			build_coeffs(single_phase_logpT, coeffs_pT);
		};
        std::string backend_name(void){return "BicubicBackend";}
        /// Build the \f$a_{i,j}\f$ coefficients for bicubic interpolation
        void build_coeffs(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs);
        void update(CoolProp::input_pairs input_pair, double val1, double val2);
		double evaluate_single_phase(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs, parameters output, double x, double y, std::size_t i, std::size_t j);
		double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j){
			return evaluate_single_phase(single_phase_logph, coeffs_ph, output, _hmolar, _p, i, j);
		};
        double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j){
			return evaluate_single_phase(single_phase_logpT, coeffs_pT, output, _T, _p, i, j);
		};
};

double do_one();
}

#endif // BICUBICBACKEND_H
