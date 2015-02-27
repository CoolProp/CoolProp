#include "BicubicBackend.h"
#include "MatrixMath.h"

/// The inverse of the A matrix for the bicubic interpolation (http://en.wikipedia.org/wiki/Bicubic_interpolation)
static const double Ainv_data[16*16] = {
     1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,
	-3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,
	 0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,
	 9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1,
	-6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1,
	 2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,
	-6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1,
	 4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1};
static Eigen::Matrix<double, 16, 16> Ainv(Ainv_data);
     
void CoolProp::BicubicBackend::update(CoolProp::input_pairs input_pair, double val1, double val2)
{
}
void CoolProp::BicubicBackend::build_coeffs_ph()
{
	const bool debug = get_debug_level() > 5 || true;
    const int param_count = 5;
    parameters param_list[param_count] = {iT, iDmolar, iSmolar, iHmolar, iP}; //iUmolar
    std::vector<std::vector<double> > *f = NULL, *fx = NULL, *fy = NULL, *fxy = NULL;
    
	SinglePhaseGriddedTableData &table = single_phase_logph;
	std::vector<std::vector<CellCoeffs> > &coeffs = coeffs_ph;
    
	clock_t t1 = clock();

    // Resize the coefficient structures
    coeffs_ph.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));
    int valid_cell_count = 0;
    for (std::size_t k = 0; k < param_count; ++k){
        parameters param = param_list[k];
		if (param == table.xkey || param == table.ykey){continue;} // Skip tables that match either of the input variables
		
		switch(param){
            case iT:
                f = &(table.T); fx = &(table.dTdx); fy = &(table.dTdy); fxy = &(table.d2Tdxdy);
                break;
			case iP:
                f = &(table.p); fx = &(table.dpdx); fy = &(table.dpdy); fxy = &(table.d2pdxdy);
                break;
            case iDmolar:
                f = &(table.rhomolar); fx = &(table.drhomolardx); *fy = (table.drhomolardy); fxy = &(table.d2rhomolardxdy);
                break;
			case iSmolar:
                f = &(table.smolar); fx = &(table.dsmolardx); *fy = (table.dsmolardy); fxy = &(table.d2smolardxdy);
                break;
			case iHmolar:
                f = &(table.hmolar); fx = &(table.dhmolardx); *fy = (table.dhmolardy); fxy = &(table.d2hmolardxdy);
                break;
            default:
                throw ValueError();
		}
        for (std::size_t i = 0; i < table.Nx-1; ++i) // -1 since we have one fewer cells than nodes
        {
            for (std::size_t j = 0; j < table.Ny-1; ++j) // -1 since we have one fewer cells than nodes
            {               
                if (ValidNumber((*f)[i][j]) && ValidNumber((*f)[i+1][j]) && ValidNumber((*f)[i][j+1]) && ValidNumber((*f)[i+1][j+1])){
                    
                    // This will hold the scaled f values for the cell
                    Eigen::Matrix<double, 16, 1> F;
                    // The output values (do not require scaling
                    F(0) = (*f)[i][j]; F(1) = (*f)[i+1][j]; F(2) = (*f)[i][j+1]; F(3) = (*f)[i+1][j+1]; 
                    // Scaling parameter
                    // d(f)/dxhat = df/dx * dx/dxhat, where xhat = (x-x_i)/(x_{i+1}-x_i)
                    coeffs[i][j].dx_dxhat = table.xvec[i+1]-table.xvec[i];
                    double dx_dxhat = coeffs[i][j].dx_dxhat;
                    F(4) = (*fx)[i][j]*dx_dxhat; F(5) = (*fx)[i+1][j]*dx_dxhat; 
                    F(6) = (*fx)[i][j+1]*dx_dxhat; F(7) = (*fx)[i+1][j+1]*dx_dxhat; 
                    // Scaling parameter
                    // d(f)/dyhat = df/dy * dy/dyhat, where yhat = (y-y_j)/(y_{j+1}-y_j)
                    coeffs[i][j].dy_dyhat = table.yvec[j+1]-table.yvec[j];
                    double dy_dyhat = coeffs[i][j].dy_dyhat;
                    F(8) = (*fy)[i][j]*dy_dyhat; F(9) = (*fy)[i+1][j]*dy_dyhat; 
                    F(10) = (*fy)[i][j+1]*dy_dyhat; F(11) = (*fy)[i+1][j+1]*dy_dyhat; 
                    // Cross derivatives are doubly scaled following the examples above
                    F(12) = (*fxy)[i][j]*dy_dyhat*dx_dxhat; F(13) = (*fxy)[i+1][j]*dy_dyhat*dx_dxhat; 
                    F(14) = (*fxy)[i][j+1]*dy_dyhat*dx_dxhat; F(15) = (*fxy)[i+1][j+1]*dy_dyhat*dx_dxhat; 
					// Calculate the alpha coefficients
                    Eigen::MatrixXd alpha = Ainv*F; // 16x1
					std::vector<double> valpha = eigen_to_vec1D(alpha);
                    coeffs[i][j].set(param, valpha);
                    coeffs[i][j].set_valid();
					valid_cell_count++;
                }
                else{
                    coeffs[i][j].set_invalid();
                }
            }
        }
    }
	double elapsed = (clock() - t1)/((double)CLOCKS_PER_SEC);
	if (debug){
		std::cout << format("Calculated bicubic coefficients for %d good cells in %g sec.", valid_cell_count, elapsed);
	}
}
