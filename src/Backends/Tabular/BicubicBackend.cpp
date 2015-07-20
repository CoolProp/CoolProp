#if !defined(NO_TABULAR_BACKENDS)

#include "BicubicBackend.h"
#include "MatrixMath.h"
#include "Backends/Helmholtz/PhaseEnvelopeRoutines.h"

void CoolProp::BicubicBackend::update(CoolProp::input_pairs input_pair, double val1, double val2)
{
    if (get_debug_level() > 0){ std::cout << format("update(%s,%g,%g)\n", get_input_pair_short_desc(input_pair).c_str(), val1, val2); }
	// Clear cached values
	clear();

    // To start, set quality to value that is for single-phase
    _Q = -1000;

	// Flush the cached indices (set to large number)
    cached_single_phase_i = std::numeric_limits<std::size_t>::max(); 
    cached_single_phase_j = std::numeric_limits<std::size_t>::max();
    cached_saturation_iL = std::numeric_limits<std::size_t>::max(); 
    cached_saturation_iV = std::numeric_limits<std::size_t>::max();
    
    PureFluidSaturationTableData &pure_saturation = dataset->pure_saturation;
    PhaseEnvelopeData & phase_envelope = dataset->phase_envelope;
    SinglePhaseGriddedTableData &single_phase_logph = dataset->single_phase_logph;
    SinglePhaseGriddedTableData &single_phase_logpT = dataset->single_phase_logpT;

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
                    CellCoeffs &cell = dataset->coeffs_ph[cached_single_phase_i][cached_single_phase_j];
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
        case PUmolar_INPUTS:
        case PSmolar_INPUTS:
        case DmolarP_INPUTS:{
            CoolPropDbl otherval; parameters otherkey;
            switch(input_pair){
                case PUmolar_INPUTS: _p = val1; _umolar = val2; otherval = val2; otherkey = iUmolar; break;
                case PSmolar_INPUTS: _p = val1; _smolar = val2; otherval = val2; otherkey = iSmolar; break;
                case DmolarP_INPUTS: _rhomolar = val1; _p = val2; otherval = val1; otherkey = iDmolar; break;
                default: throw ValueError("Bad (impossible) pair");
            }
            
            using_single_phase_table = true; // Use the table (or first guess is that it is single-phase)!
            std::size_t iL, iV;
            CoolPropDbl zL = 0, zV = 0;
            std::size_t iclosest = 0;
            SimpleState closest_state;
            if (
                (is_mixture && PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, otherkey, otherval, iclosest, closest_state))
                ||
                (!is_mixture && pure_saturation.is_inside(iP, _p, otherkey, otherval, iL, iV, zL, zV))
                ){
                using_single_phase_table = false;
                if (otherkey == iDmolar){
                    _Q = (1/otherval - 1/zL)/(1/zV - 1/zL);
                }
                else{
                    _Q = (otherval - zL)/(zV - zL);
                }
                if(!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))){
                    throw ValueError("vapor quality is not in (0,1)");
                }
                else{
                    cached_saturation_iL = iL; cached_saturation_iV = iV;
                }
            }
            else{
                // Find and cache the indices i, j
                selected_table = SELECTED_PH_TABLE;
                single_phase_logph.find_nearest_neighbor(iP, _p, otherkey, otherval, cached_single_phase_i, cached_single_phase_j);
                CellCoeffs &cell = dataset->coeffs_ph[cached_single_phase_i][cached_single_phase_j];
                if (!cell.valid()){
                    if (cell.has_valid_neighbor()){
                        // Get new good neighbor
                        cell.get_alternate(cached_single_phase_i, cached_single_phase_j);
                    }
                    else{
                        if (!cell.valid()){throw ValueError(format("Cell is invalid and has no good neighbors for p = %g Pa, T= %g K",val1,val2));}
                    }
                }
				// Now find hmolar given P, X for X in Hmolar, Smolar, Umolar
                invert_single_phase_x(single_phase_logph, dataset->coeffs_ph, otherkey, otherval, _p, cached_single_phase_i, cached_single_phase_j);
            }
            break;
        }
        case DmassP_INPUTS:{
            // Call again, but this time with molar units; D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmassP_INPUTS, val1 / AS->molar_mass(), val2); return;
        }
        case PUmass_INPUTS:{
            // Call again, but this time with molar units; U: [J/kg] * [kg/mol] -> [J/mol]
            update(PUmolar_INPUTS, val1, val2*AS->molar_mass()); return;
        }
        case PSmass_INPUTS:{
            // Call again, but this time with molar units; S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(PSmolar_INPUTS, val1, val2*AS->molar_mass()); return;
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
                    CellCoeffs &cell = dataset->coeffs_pT[cached_single_phase_i][cached_single_phase_j];
                    if (!cell.valid()){
                        if (cell.has_valid_neighbor()){
                            // Get new good neighbor
                            cell.get_alternate(cached_single_phase_i, cached_single_phase_j);
                        }
                        else{
                            if (!cell.valid()){throw ValueError(format("Cell is invalid and has no good neighbors for p = %g Pa, T= %g K",val1,val2));}
                        }
                    }

                    // If p < pc, you might be getting a liquid solution when you want a vapor solution or vice versa
                    // if you are very close to the saturation curve, so we figure out what the saturation temperature
                    // is for the given pressure
                    if (_p < this->AS->p_critical())
                    {
                        double Ts = pure_saturation.evaluate(iT, _p, _Q, iL, iV);
                        double TL = single_phase_logpT.T[cached_single_phase_i][cached_single_phase_j];
                        double TR = single_phase_logpT.T[cached_single_phase_i+1][cached_single_phase_j];
                        if (TL < Ts && Ts < TR){
                            if (_T < Ts){
                                if (cached_single_phase_i == 0){throw ValueError(format("P, T are near saturation, but cannot move the cell to the left")); }
                                // It's liquid, move the cell to the left
                                cached_single_phase_i--;
                            }else{
                                if (cached_single_phase_i > single_phase_logpT.Nx-2){ throw ValueError(format("P,T are near saturation, but cannot move the cell to the right")); }
                                // It's vapor, move to the right
                                cached_single_phase_i++;
                            }
                        }
                    }   
                }
            }
            break;
        }
        case DmassT_INPUTS:{
            // Call again, but this time with molar units; D: [kg/m^3] / [kg/mol] -> [mol/m^3]
            update(DmolarT_INPUTS, val1 / AS->molar_mass(), val2); return;
        }
        case SmassT_INPUTS:{
            // Call again, but this time with molar units; S: [J/kg/K] * [kg/mol] -> [J/mol/K]
            update(SmolarT_INPUTS, val1*AS->molar_mass(), val2); return;
        }
        case SmolarT_INPUTS:
        case DmolarT_INPUTS:{
            CoolPropDbl otherval; parameters otherkey;
            switch(input_pair){
                case SmolarT_INPUTS: _smolar = val1; _T = val2; otherval = val1; otherkey = iSmolar; break;
                case DmolarT_INPUTS: _rhomolar = val1; _T = val2; otherval = val1; otherkey = iDmolar; break;
                default: throw ValueError("Bad (impossible) pair");
            }
            
            using_single_phase_table = true; // Use the table (or first guess is that it is single-phase)!
            std::size_t iL, iV;
            CoolPropDbl zL = 0, zV = 0;
            std::size_t iclosest = 0;
            SimpleState closest_state;
            if (
                (is_mixture && PhaseEnvelopeRoutines::is_inside(phase_envelope, iT, _T, otherkey, otherval, iclosest, closest_state))
                ||
                (!is_mixture && pure_saturation.is_inside(iT, _T, otherkey, otherval, iL, iV, zL, zV))
                )
            {
                using_single_phase_table = false;
                if (otherkey == iDmolar){
                    _Q = (1/otherval - 1/zL)/(1/zV - 1/zL);
                }
                else{
                    _Q = (otherval - zL)/(zV - zL);
                }
                if(!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))){
                    throw ValueError("vapor quality is not in (0,1)");
                }
                else{
                    cached_saturation_iL = iL; cached_saturation_iV = iV;
                }
                _p = pure_saturation.evaluate(iP, _T, _Q, iL, iV);
            }
            else{
                // Find and cache the indices i, j
                selected_table = SELECTED_PT_TABLE;
                single_phase_logpT.find_nearest_neighbor(iT, _T, otherkey, otherval, cached_single_phase_i, cached_single_phase_j);
                CellCoeffs &cell = dataset->coeffs_pT[cached_single_phase_i][cached_single_phase_j];
                if (!cell.valid()){
                    if (cell.has_valid_neighbor()){
                        // Get new good neighbor
                        cell.get_alternate(cached_single_phase_i, cached_single_phase_j);
                    }
                    else{
                        if (!cell.valid()){throw ValueError(format("Cell is invalid and has no good neighbors for p = %g Pa, T= %g K",val1,val2));}
                    }
                }
				// Now find the y variable (Dmolar or Smolar in this case)
                invert_single_phase_y(single_phase_logpT, dataset->coeffs_pT, otherkey, otherval, _T, cached_single_phase_i, cached_single_phase_j);
            }
            break;
        }
		case PQ_INPUTS:{
			std::size_t iL = 0, iV = 0;
			_p = val1; _Q = val2;
            
            using_single_phase_table = false;
            if(!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))){
                throw ValueError("vapor quality is not in (0,1)");
            }
            else{
                if (is_mixture){
                    std::vector<std::pair<std::size_t, std::size_t> > intersect = PhaseEnvelopeRoutines::find_intersections(phase_envelope, iP, _p);
                    if (intersect.empty()){ throw ValueError(format("p [%g Pa] is not within phase envelope", _p)); }
                    iV = intersect[0].first; iL = intersect[1].first;
                }
                else{
                    CoolPropDbl zL, zV;
                    pure_saturation.is_inside(iP, _p, iQ, _Q, iL, iV, zL, zV);
                }
                cached_saturation_iL = iL; cached_saturation_iV = iV;
            }
			break;
		}
        case QT_INPUTS:{
			std::size_t iL = 0, iV = 0;
			_Q = val1; _T = val2;
            
            using_single_phase_table = false;
            if(!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))){
                throw ValueError("vapor quality is not in (0,1)");
            }
            else{
                if (is_mixture){
                    std::vector<std::pair<std::size_t,std::size_t> > intersect = PhaseEnvelopeRoutines::find_intersections(phase_envelope, iT, _T);
                    if (intersect.empty()){ throw ValueError(format("T [%g K] is not within phase envelope", _T)); }
                    iV = intersect[0].first; iL = intersect[1].first;
                    double pL = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iL);
                    double pV = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iV);
                    _p = _Q*pV + (1-_Q)*pL;
                }
                else{
                    CoolPropDbl zL, zV;
                    pure_saturation.is_inside(iT, _T, iQ, _Q, iL, iV, zL, zV);
                }
                cached_saturation_iL = iL; cached_saturation_iV = iV;
            }
			break;
		}
		default:
			throw ValueError("Sorry, but this set of inputs is not supported for Bicubic backend");
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
        default: throw ValueError("Invalid output variable in evaluate_single_phase_transport");
    }
    return val;
}
// Use the single_phase table to evaluate an output
double CoolProp::BicubicBackend::evaluate_single_phase(const SinglePhaseGriddedTableData &table, const std::vector<std::vector<CellCoeffs> > &coeffs, const parameters output, const double x, const double y, const std::size_t i, const std::size_t j)
{
    // Get the cell
    const CellCoeffs &cell = coeffs[i][j];
    
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
        case iUmolar: _umolar = val; break;
        default: throw ValueError("Invalid output variable in evaluate_single_phase");
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

/// Use the single_phase table to invert for x given a y
void CoolProp::BicubicBackend::invert_single_phase_x(const SinglePhaseGriddedTableData &table, const std::vector<std::vector<CellCoeffs> > &coeffs, parameters other_key, double other, double y, std::size_t i, std::size_t j)
{
    // Get the cell
    const CellCoeffs &cell = coeffs[i][j];
    
	// Get the alpha coefficients
    const std::vector<double> &alpha = cell.get(other_key);
    
    // Normalized value in the range (0, 1)
    double yhat = (y - table.yvec[j])/(table.yvec[j+1] - table.yvec[j]);

    double y_0 = 1, y_1 = yhat, y_2 = yhat*yhat, y_3 = yhat*yhat*yhat;

    double a = alpha[3+0*4]*y_0+alpha[3+1*4]*y_1+alpha[3+2*4]*y_2+alpha[3+3*4]*y_3; // factors of xhat^3
    double b = alpha[2+0*4]*y_0+alpha[2+1*4]*y_1+alpha[2+2*4]*y_2+alpha[2+3*4]*y_3; // factors of xhar^2
    double c = alpha[1+0*4]*y_0+alpha[1+1*4]*y_1+alpha[1+2*4]*y_2+alpha[1+3*4]*y_3; // factors of xhat
    double d = alpha[0+0*4]*y_0+alpha[0+1*4]*y_1+alpha[0+2*4]*y_2+alpha[0+3*4]*y_3 - other; // constant factors
    int N = 0;
    double xhat0, xhat1, xhat2, val, xhat = _HUGE;
    solve_cubic(a, b, c, d, N, xhat0, xhat1, xhat2);
    if (N == 1){
        xhat = xhat0;
    }
    else if (N == 2){
        xhat = std::abs(xhat0) < std::abs(xhat1) ? xhat0 : xhat1;
    }
    else if (N == 3){
        if (std::abs(xhat0) < std::abs(xhat1) && std::abs(xhat0) < std::abs(xhat2)){
            xhat = xhat0;
        }
        // Already know that xhat1 < xhat0 (xhat0 is not the minimum)
        else if (std::abs(xhat1) < std::abs(xhat2)){
            xhat = xhat1;
        }
        else{
            xhat = xhat2;
        }
    }
    else if (N == 0){
        throw ValueError("Could not find a solution in invert_single_phase_x");
    }

    // Unpack xhat into actual value
    // xhat = (x-x_{i})/(x_{i+1}-x_{i})
    val = xhat*(table.xvec[i+1] - table.xvec[i]) + table.xvec[i];
    
    // Cache the output value calculated
    switch(table.xkey){
        case iHmolar: _hmolar = val; break;
        case iT: _T = val; break;
        default: throw ValueError("Invalid output variable in invert_single_phase_x");
    }
}

/// Use the single_phase table to solve for y given an x
void CoolProp::BicubicBackend::invert_single_phase_y(const SinglePhaseGriddedTableData &table, const std::vector<std::vector<CellCoeffs> > &coeffs, parameters other_key, double other, double x, std::size_t i, std::size_t j)
{
    // Get the cell
    const CellCoeffs &cell = coeffs[i][j];
    
	// Get the alpha coefficients
    const std::vector<double> &alpha = cell.get(other_key);
    
    // Normalized value in the range (0, 1)
    double xhat = (x - table.xvec[i])/(table.xvec[i+1] - table.xvec[i]);

    double x_0 = 1, x_1 = xhat, x_2 = xhat*xhat, x_3 = xhat*xhat*xhat;

    double a = alpha[0+3*4]*x_0 + alpha[1+3*4]*x_1 + alpha[2+3*4]*x_2 + alpha[3+3*4]*x_3; // factors of yhat^3 (m= 3)
    double b = alpha[0+2*4]*x_0 + alpha[1+2*4]*x_1 + alpha[2+2*4]*x_2 + alpha[3+2*4]*x_3; // factors of yhat^2
    double c = alpha[0+1*4]*x_0 + alpha[1+1*4]*x_1 + alpha[2+1*4]*x_2 + alpha[3+1*4]*x_3; // factors of yhat
    double d = alpha[0+0*4]*x_0 + alpha[1+0*4]*x_1 + alpha[2+0*4]*x_2 + alpha[3+0*4]*x_3 - other; // constant factors
    int N = 0;
    double yhat0, yhat1, yhat2, val, yhat = _HUGE;
    solve_cubic(a, b, c, d, N, yhat0, yhat1, yhat2);
    if (N == 1){
        yhat = yhat0;
    }
    else if (N == 2){
        yhat = std::abs(yhat0) < std::abs(yhat1) ? yhat0 : yhat1;
    }
    else if (N == 3){
        if (std::abs(yhat0) < std::abs(yhat1) && std::abs(yhat0) < std::abs(yhat2)){
            yhat = yhat0;
        }
        // Already know that yhat1 < yhat0 (yhat0 is not the minimum)
        else if (std::abs(yhat1) < std::abs(yhat2)){
            yhat = yhat1;
        }
        else{
            yhat = yhat2;
        }
    }
    else if (N == 0){
        throw ValueError("Could not find a solution in invert_single_phase_x");
    }

    // Unpack xhat into actual value
    // yhat = (y-y_{j})/(y_{j+1}-y_{j})
    val = yhat*(table.yvec[j+1] - table.yvec[j]) + table.yvec[j];
    
    // Cache the output value calculated
    switch(table.ykey){
        case iP: _p = val; break;
        default: throw ValueError("Invalid output variable in invert_single_phase_x");
    }
}

#endif // !defined(NO_TABULAR_BACKENDS)