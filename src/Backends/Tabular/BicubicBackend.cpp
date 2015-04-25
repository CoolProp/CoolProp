#if !defined(NO_TABULAR_BACKENDS)

#include "BicubicBackend.h"
#include "MatrixMath.h"
#include "../Helmholtz/PhaseEnvelopeRoutines.h"

/// The inverse of the A matrix for the bicubic interpolation (http://en.wikipedia.org/wiki/Bicubic_interpolation)
/// NOTE: The matrix is transposed below
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

void CoolProp::BicubicBackend::build_coeffs(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs)
{
    if (!coeffs.empty()){ return; }
	const bool debug = get_debug_level() > 5 || false;
    const int param_count = 6;
    parameters param_list[param_count] = {iDmolar, iT, iSmolar, iHmolar, iP, iUmolar};
    std::vector<std::vector<double> > *f = NULL, *fx = NULL, *fy = NULL, *fxy = NULL;
    
	clock_t t1 = clock();

    // Resize the coefficient structures
    coeffs.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));

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
                f = &(table.rhomolar); fx = &(table.drhomolardx); fy = &(table.drhomolardy); fxy = &(table.d2rhomolardxdy);
                break;
			case iSmolar:
                f = &(table.smolar); fx = &(table.dsmolardx); fy = &(table.dsmolardy); fxy = &(table.d2smolardxdy);
                break;
			case iHmolar:
                f = &(table.hmolar); fx = &(table.dhmolardx); fy = &(table.dhmolardy); fxy = &(table.d2hmolardxdy);
                break;
			case iUmolar:
				f = &(table.umolar); fx = &(table.dumolardx); fy = &(table.dumolardy); fxy = &(table.d2umolardxdy);
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
					Eigen::MatrixXd alpha = Ainv.transpose()*F; // 16x1; Watch out for the transpose!
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
        double elapsed = (clock() - t1)/((double)CLOCKS_PER_SEC);
        if (debug){
            std::cout << format("Calculated bicubic coefficients for %d good cells in %g sec.\n", valid_cell_count, elapsed);
        }
        std::size_t remap_count = 0;
        // Now find invalid cells and give them pointers to a neighboring cell that works
        for (std::size_t i = 0; i < table.Nx-1; ++i) // -1 since we have one fewer cells than nodes
        {
            for (std::size_t j = 0; j < table.Ny-1; ++j) // -1 since we have one fewer cells than nodes
            {
                // Not a valid cell
                if (!coeffs[i][j].valid()){
                    // Offsets that we are going to try in order (left, right, top, bottom, diagonals)
                    int xoffsets[] = {-1,1,0,0,-1,1,1,-1};
                    int yoffsets[] = {0,0,1,-1,-1,-1,1,1};
                    // Length of offset
                    std::size_t N = sizeof(xoffsets)/sizeof(xoffsets[0]);
                    for (std::size_t k = 0; k < N; ++k){
                        std::size_t iplus = i + xoffsets[k];
                        std::size_t jplus = j + yoffsets[k];
                        if (0 < iplus && iplus < table.Nx-1 && 0 < jplus && jplus < table.Ny-1 && coeffs[iplus][jplus].valid()){
                            coeffs[i][j].set_alternate(iplus, jplus);
                            remap_count++;
                            if (debug){std::cout << format("Mapping %d,%d to %d,%d\n",i,j,iplus,jplus);}
                            break;
                        }
                    }
                }
            }
        }
        if (debug){
            std::cout << format("Remapped %d cells\n", remap_count);
        }
    }
}

void CoolProp::BicubicBackend::update(CoolProp::input_pairs input_pair, double val1, double val2)
{
	// Clear cached values
	clear();

    // Check the tables and build if necessary
    check_tables();

    // To start, set quality to value that is for single-phase
    _Q = -1000;

    bool is_mixture = (this->AS->get_mole_fractions().size() >= 2);

    if (is_mixture){
        // For mixtures, the construction of the coefficients is delayed until this 
        // function so that the set_mole_fractions function can be called
        build_coeffs(single_phase_logph, coeffs_ph);
        build_coeffs(single_phase_logpT, coeffs_pT);
    }

	// Flush the cached indices (set to large number)
    cached_single_phase_i = std::numeric_limits<std::size_t>::max(); 
    cached_single_phase_j = std::numeric_limits<std::size_t>::max();
    cached_saturation_iL = std::numeric_limits<std::size_t>::max(); 
    cached_saturation_iV = std::numeric_limits<std::size_t>::max();
    
    switch(input_pair){
        case HmolarP_INPUTS:{
            _hmolar = val1; _p = val2;
            if (!single_phase_logph.native_inputs_are_in_range(_hmolar, _p)){
                // Use the AbstractState instance
                using_single_phase_table = false;
                if (get_debug_level() > 5){ std::cout << "inputs are not in range"; }
                throw ValueError(format("inputs are not in range, hmolar=%Lg, p=%Lg", static_cast<CoolPropDbl>(_hmolar), _p));
            }
            else{
                using_single_phase_table = true; // Use the table !
                std::size_t iL, iV, iclosest = 0;
                CoolPropDbl hL = 0, hV = 0;
                SimpleState closest_state;
                if (
                    (is_mixture && PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, iHmolar, _hmolar, iclosest, closest_state))
                    ||
                    (!is_mixture && pure_saturation.is_inside(iP, _p, iHmolar, _hmolar, iL, iV, hL, hV))
                    )
                {
                    using_single_phase_table = false;
                    _Q = (static_cast<double>(_hmolar)-hL)/(hV-hL);
                    if(!is_in_closed_range(0.0,1.0,static_cast<double>(_Q))){
                        throw ValueError("vapor quality is not in (0,1)");
                    }
                    else{
                        cached_saturation_iL = iL; cached_saturation_iV = iV;
                    }
                }
                else{
                    // Find and cache the indices i, j
                    selected_table = SELECTED_PH_TABLE;
					single_phase_logph.find_native_nearest_good_cell(_hmolar, _p, cached_single_phase_i, cached_single_phase_j);
                    CellCoeffs &cell = coeffs_ph[cached_single_phase_i][cached_single_phase_j];
                    if (!cell.valid()){
                        if (cell.has_valid_neighbor()){
                            // Get new good neighbor
                            cell.get_alternate(cached_single_phase_i, cached_single_phase_j);
                        }
                        else{
                            if (!cell.valid()){throw ValueError(format("Cell is invalid and has no good neighbors for hmolar = %g, p= %g",val1,val2));}
                        }
                    }
                }
            }
            break;
        }
        case HmassP_INPUTS:{
            update(HmolarP_INPUTS, val1 * AS->molar_mass(), val2); // H: [J/kg] * [kg/mol] -> [J/mol]
            return;
        }
		case PT_INPUTS:{
            _p = val1; _T = val2;
            if (!single_phase_logpT.native_inputs_are_in_range(_T, _p)){
                // Use the AbstractState instance
                using_single_phase_table = false;
                if (get_debug_level() > 5){ std::cout << "inputs are not in range"; }
                throw ValueError(format("inputs are not in range, p=%Lg, T=%Lg", _p, _T));
            }
            else{
                using_single_phase_table = true; // Use the table !
                std::size_t iL = 0, iV = 0, iclosest = 0;
                CoolPropDbl TL = 0, TV = 0;
                SimpleState closest_state;
                if (
                    (is_mixture && PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, iT, _T, iclosest, closest_state))
                    ||
                    (!is_mixture && pure_saturation.is_inside(iP, _p, iT, _T, iL, iV, TL, TV))
                    )
                {
                    using_single_phase_table = false;
                    throw ValueError(format("P,T with TTSE cannot be two-phase for now"));
                }
                else{
                    // Find and cache the indices i, j
                    selected_table = SELECTED_PT_TABLE;
					single_phase_logpT.find_native_nearest_good_cell(_T, _p, cached_single_phase_i, cached_single_phase_j);
                }
            }
            break;
        }
		case PQ_INPUTS:{
			std::size_t iL = 0, iV = 0;
			CoolPropDbl hL = 0, hV = 0;
			_p = val1; _Q = val2;
			pure_saturation.is_inside(iP, _p, iQ, _Q, iL, iV, hL, hV);
            using_single_phase_table = false;
            if(!is_in_closed_range(0.0,1.0,static_cast<double>(_Q))){
                throw ValueError("vapor quality is not in (0,1)");
            }
            else{
                cached_saturation_iL = iL; cached_saturation_iV = iV;
            }
			break;
		}
		default:
			throw ValueError("input pair is not currently supported");
	}
}

/** Use the single_phase table to evaluate an output for a transport property
 * 
 * Here we use linear interpolation because we don't have any information about the derivatives with respect to the 
 * independent variables and it is too computationally expensive to build the derivatives numerically
 * 
 * See also http://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
 */
double CoolProp::BicubicBackend::evaluate_single_phase_transport(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j)
{
    // By definition i,i+1,j,j+1 are all in range and valid
    std::vector<std::vector<double> > *f = NULL;
    switch(output){
        case iconductivity:
            f = &table.cond; break;
        case iviscosity:
            f = &table.visc; break;
        default:
            throw ValueError(format("invalid output variable to BicubicBackend::evaluate_single_phase_transport"));
    }
    double x1 = table.xvec[i], x2 = table.xvec[i+1], y1 = table.yvec[j], y2 = table.yvec[j+1];
    double f11 = (*f)[i][j], f12 = (*f)[i][j+1], f21 = (*f)[i+1][j], f22 = (*f)[i+1][j+1];
    double val = 1/((x2-x1)*(y2-y1))*( f11*(x2 - x)*(y2 - y)
                                      +f21*(x - x1)*(y2 - y)
                                      +f12*(x2 - x)*(y - y1)
                                      +f22*(x - x1)*(y - y1));
    
    // Cache the output value calculated
    switch(output){
        case iconductivity: _conductivity = val; break;
        case iviscosity: _viscosity = val; break;
        default: throw ValueError();
    }
    return val;
}
/// Use the single_phase table to evaluate an output
double CoolProp::BicubicBackend::evaluate_single_phase(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs, parameters output, double x, double y, std::size_t i, std::size_t j)
{
    // Get the cell
    CellCoeffs &cell = coeffs[i][j];
    
	// Get the alpha coefficients
    const std::vector<double> &alpha = cell.get(output);
    
    // Normalized value in the range (0, 1)
	double xhat = (x - table.xvec[i])/(table.xvec[i+1] - table.xvec[i]);
    double yhat = (y - table.yvec[j])/(table.yvec[j+1] - table.yvec[j]);
    
    // Calculate the output value desired
	double val = 0;
	for (std::size_t l = 0; l < 4; ++l)
	{
		for(std::size_t m = 0; m < 4; ++m)
		{
			val += alpha[m*4+l]*pow(xhat, static_cast<int>(l))*pow(yhat, static_cast<int>(m));
		}
	}
    
    // Cache the output value calculated
    switch(output){
        case iT:  _T = val; break;
        case iDmolar: _rhomolar = val; break;
        case iSmolar: _smolar = val; break;
		case iHmolar: _hmolar = val; break;
        //case iUmolar:
        default: throw ValueError();
    }
    return val;
}
/// Use the single_phase table to evaluate an output
double CoolProp::BicubicBackend::evaluate_single_phase_derivative(SinglePhaseGriddedTableData &table, std::vector<std::vector<CellCoeffs> > &coeffs, parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny)
{

    // Get the cell
    CellCoeffs &cell = coeffs[i][j];
    
	// Get the alpha coefficients
    const std::vector<double> &alpha = cell.get(output);
    
    // Normalized value in the range (0, 1)
	double xhat = (x - table.xvec[i])/(table.xvec[i+1] - table.xvec[i]);
    double yhat = (y - table.yvec[j])/(table.yvec[j+1] - table.yvec[j]);
    double dxhatdx = 1/(table.xvec[i+1] - table.xvec[i]);
    double dyhatdy = 1/(table.yvec[j+1] - table.yvec[j]);
    
    // Calculate the output value desired
	double val = 0;
    if (Nx == 1 && Ny == 0){
		if (output == table.xkey) { return 1.0; }
		if (output == table.ykey) { return 0.0; }
        for (std::size_t l = 1; l < 4; ++l)
        {
            for(std::size_t m = 0; m < 4; ++m)
            {
                val += alpha[m*4+l]*l*pow(xhat, static_cast<int>(l-1))*pow(yhat, static_cast<int>(m));
            }
        }
        // val is now dz/dxhat|yhat
        return val*dxhatdx;
    }
    else if (Ny == 1 && Nx == 0){
		if (output == table.ykey) { return 1.0; }
		if (output == table.xkey) { return 0.0; }
        for (std::size_t l = 0; l < 4; ++l)
        {
            for(std::size_t m = 1; m < 4; ++m)
            {
                val += alpha[m*4+l]*pow(xhat, static_cast<int>(l))*m*pow(yhat, static_cast<int>(m-1));
            }
        }
        // val is now dz/dyhat|xhat
        return val*dyhatdy;
    }
    else{
        throw ValueError("Invalid input");
    }
}

#endif // !defined(NO_TABULAR_BACKENDS)