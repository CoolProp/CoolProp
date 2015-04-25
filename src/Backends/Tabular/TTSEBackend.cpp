#if !defined(NO_TABULAR_BACKENDS)


#include "TTSEBackend.h"
#include "CoolProp.h"

void CoolProp::TTSEBackend::update(CoolProp::input_pairs input_pair, double val1, double val2)
{
    // Clear cached variables
    clear();

    // Check the tables, build if neccessary
    check_tables();
    
    // Flush the cached indices (set to large number)
    cached_single_phase_i = std::numeric_limits<std::size_t>::max(); 
    cached_single_phase_j = std::numeric_limits<std::size_t>::max();
    cached_saturation_iL = std::numeric_limits<std::size_t>::max(); 
    cached_saturation_iV = std::numeric_limits<std::size_t>::max();

    // To start, set quality to value that is impossible
    _Q = -1000;
    
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
                std::size_t iL, iV;
                CoolPropDbl hL = 0, hV = 0;
                if (pure_saturation.is_inside(iP, _p, iHmolar, _hmolar, iL, iV, hL, hV)){
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
                    single_phase_logph.find_native_nearest_neighbor(_hmolar, _p, cached_single_phase_i, cached_single_phase_j);
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
            if (pure_saturation.is_inside(iP, _p, otherkey, otherval, iL, iV, zL, zV)){
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
				// Now find hmolar
                invert_single_phase_x(single_phase_logph, otherkey, otherval, _p, cached_single_phase_i, cached_single_phase_j);
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
                std::size_t iL = 0, iV = 0;
                CoolPropDbl TL = 0, TV = 0;
                if (pure_saturation.is_inside(iP, _p, iT, _T, iL, iV, TL, TV)){
                    using_single_phase_table = false;
                    throw ValueError(format("P,T with TTSE cannot be two-phase for now"));
                }
                else{
                    // Find and cache the indices i, j
                    selected_table = SELECTED_PT_TABLE;
                    single_phase_logpT.find_native_nearest_neighbor(_T, _p, cached_single_phase_i, cached_single_phase_j);
                }
            }
            break;
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
            if (pure_saturation.is_inside(iT, _T, otherkey, otherval, iL, iV, zL, zV)){
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
                selected_table = SELECTED_PT_TABLE;
                single_phase_logpT.find_nearest_neighbor(iT, _T, otherkey, otherval, cached_single_phase_i, cached_single_phase_j);
				// Now find hmolar
                invert_single_phase_y(single_phase_logpT, otherkey, otherval, _T, cached_single_phase_i, cached_single_phase_j);
            }
            break;
        }
        default:
            throw ValueError("Sorry, but this set of inputs is not supported for TTSE backend");
    }
}
/** Use the single_phase table to evaluate an output for a transport property
 * 
 * Here we use bilinear interpolation because we don't have any information about the derivatives with respect to the 
 * independent variables and it is too computationally expensive to build the derivatives numerically
 * 
 * See also http://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
 */
double CoolProp::TTSEBackend::evaluate_single_phase_transport(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j)
{
    bool in_bounds = (i >=0 && i < table.xvec.size()-1 && j >=0 && j < table.yvec.size()-1);
    if (!in_bounds){
        throw ValueError("Cell to TTSEBackend::evaluate_single_phase_transport is not valid");
    }
    bool is_valid = (ValidNumber(table.smolar[i][j]) && ValidNumber(table.smolar[i+1][j]) && ValidNumber(table.smolar[i][j+1]) && ValidNumber(table.smolar[i+1][j+1]));
    if (!is_valid){
        throw ValueError("Cell to TTSEBackend::evaluate_single_phase_transport must have four valid corners for now");
    }
    const std::vector<std::vector<double> > &f = table.get(output);

    double x1 = table.xvec[i], x2 = table.xvec[i+1], y1 = table.yvec[j], y2 = table.yvec[j+1];
    double f11 = f[i][j], f12 = f[i][j+1], f21 = f[i+1][j], f22 = f[i+1][j+1];
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
/// Solve for deltax
double CoolProp::TTSEBackend::invert_single_phase_x(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j)
{   
    connect_pointers(output, table);
    
    // Distances from the node
    double deltay = y - table.yvec[j];
    
    // Calculate the output value desired
    double a = 0.5*(*d2zdx2)[i][j]; // Term multiplying deltax**2
    double b = (*dzdx)[i][j] + deltay*(*d2zdxdy)[i][j]; // Term multiplying deltax
    double c = (*z)[i][j] - x + deltay*(*dzdy)[i][j] + 0.5*deltay*deltay*(*d2zdy2)[i][j];

    double deltax1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    double deltax2 = (-b - sqrt(b*b - 4*a*c))/(2*a);

    // If only one is less than a multiple of x spacing, thats your solution
    double xspacing, xratio, val;
    if (!table.logx){
        xspacing = table.xvec[1] - table.xvec[0];
        if (std::abs(deltax1) < xspacing && !(std::abs(deltax2) < xspacing) ){
		    val = deltax1 + table.xvec[i];
        }
        else if (std::abs(deltax2) < xspacing && !(std::abs(deltax1) < xspacing) ){
		    val = deltax2 + table.xvec[i];
        }
        else if (std::abs(deltax1) < std::abs(deltax2) && std::abs(deltax1) < 10*xspacing){
            val = deltax1 + table.xvec[i];
        }
        else{
            throw ValueError(format("Cannot find the x solution; xspacing: %g dx1: %g dx2: %g", xspacing, deltax1, deltax2));
        }
    }else{
        xratio = table.xvec[1]/table.xvec[0];
        double xj = table.xvec[j];
        double xratio1 = (xj+deltax1)/xj;
        double xratio2 = (xj+deltax2)/xj;
        if (xratio1 < xratio && xratio1 > 1/xratio ){
		    val = deltax1 + table.xvec[i];
        }
        else if (xratio2 < xratio && xratio2 > 1/xratio ){
		    val = deltax2 + table.xvec[i];
        }
        else if (xratio1 < xratio*5 && xratio1 > 1/xratio/5 ){
		    val = deltax1 + table.xvec[i];
        }
        else{
            throw ValueError(format("Cannot find the x solution; xj: %g xratio: %g xratio1: %g xratio2: %g a: %g b^2-4*a*c %g", xj, xratio, xratio1, xratio2, a, b*b-4*a*c));
        }
    }

    // Cache the output value calculated
    switch(table.xkey){
        case iHmolar: _hmolar = val; break;
        case iT: _T = val; break;
        default: throw ValueError();
    }
    return val;
}
/// Solve for deltay
double CoolProp::TTSEBackend::invert_single_phase_y(SinglePhaseGriddedTableData &table, parameters output, double y, double x, std::size_t i, std::size_t j)
{   
    connect_pointers(output, table);
    
    // Distances from the node
    double deltax = x - table.xvec[i];

    // Calculate the output value desired
    double a = 0.5*(*d2zdy2)[i][j]; // Term multiplying deltay**2
    double b = (*dzdy)[i][j] + deltax*(*d2zdxdy)[i][j]; // Term multiplying deltay
    double c = (*z)[i][j] - y + deltax*(*dzdx)[i][j] + 0.5*deltax*deltax*(*d2zdx2)[i][j];

    double deltay1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    double deltay2 = (-b - sqrt(b*b - 4*a*c))/(2*a);

    // If only one is less than a multiple of x spacing, thats your solution
    double yspacing, yratio, val;
    if (!table.logy){
        yspacing = table.yvec[1] - table.yvec[0];
        if (std::abs(deltay1) < yspacing && !(std::abs(deltay2) < yspacing) ){
		    val = deltay1 + table.yvec[j];
        }
        else if (std::abs(deltay2) < yspacing && !(std::abs(deltay1) < yspacing) ){
		    val = deltay2 + table.yvec[j];
        }
        else if (std::abs(deltay1) < std::abs(deltay2) && std::abs(deltay1) < 10*yspacing){
            val = deltay1 + table.yvec[j];
        }
        else{
            throw ValueError(format("Cannot find the y solution; yspacing: %g dy1: %g dy2: %g", yspacing, deltay1, deltay2));
        }
    }else{
        yratio = table.yvec[1]/table.yvec[0];
        double yj = table.yvec[j];
        double yratio1 = (yj+deltay1)/yj;
        double yratio2 = (yj+deltay2)/yj;
        if (yratio1 < yratio && yratio1 > 1/yratio ){
		    val = deltay1 + table.yvec[j];
        }
        else if (yratio2 < yratio && yratio2 > 1/yratio ){
		    val = deltay2 + table.yvec[j];
        }
        else if (std::abs(yratio1-1) < std::abs(yratio2-1)){
		    val = deltay1 + table.yvec[j];
        }
        else if (std::abs(yratio2-1) < std::abs(yratio1-1)){
		    val = deltay2 + table.yvec[j];
        }
        else{
            throw ValueError(format("Cannot find the y solution; yj: %g yratio: %g yratio1: %g yratio2: %g a: %g b: %g b^2-4ac: %g %d %d", yj, yratio, yratio1, yratio2, a, b, b*b-4*a*c, i, j));
        }
    }

    // Cache the output value calculated
    switch(table.ykey){
        case iHmolar: _hmolar = val; break;
        case iT: _T = val; break;
        case iP: _p = val; break;
        default: throw ValueError();
    }
    return val;
}
/// Use the single-phase table to evaluate an output
double CoolProp::TTSEBackend::evaluate_single_phase(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j)
{
	connect_pointers(output, table);

    // Distances from the node
	double deltax = x - table.xvec[i];
    double deltay = y - table.yvec[j];
    
    // Calculate the output value desired
    double val = (*z)[i][j]+deltax*(*dzdx)[i][j]+deltay*(*dzdy)[i][j]+0.5*deltax*deltax*(*d2zdx2)[i][j]+0.5*deltay*deltay*(*d2zdy2)[i][j]+deltay*deltax*(*d2zdxdy)[i][j];
    
    // Cache the output value calculated
    switch(output){
        case iT:  _T = val; break;
        case iDmolar: _rhomolar = val; break;
        case iSmolar: _smolar = val; break;
        case iHmolar: _hmolar = val; break;
        case iUmolar: _umolar = val; break;
        default: throw ValueError();
    }
    return val;
}
/// Use the single-phase table to evaluate an output
double CoolProp::TTSEBackend::evaluate_single_phase_derivative(SinglePhaseGriddedTableData &table, parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny)
{
	if (Nx == 1 && Ny == 0){
		if (output == table.xkey) { return 1.0; }
		if (output == table.ykey) { return 0.0; }
	}
	else if (Ny == 1 && Nx == 0){
		if (output == table.ykey) { return 1.0; }
		if (output == table.xkey) { return 0.0; }
	}
    
    connect_pointers(output, table);
    
    // Distances from the node
	double deltax = x - table.xvec[i];
    double deltay = y - table.yvec[j];
    double val;
    // Calculate the output value desired
    if (Nx == 1 && Ny == 0){
		if (output == table.xkey) { return 1.0; }
		if (output == table.ykey) { return 0.0; }
		val = (*dzdx)[i][j] + deltax*(*d2zdx2)[i][j] + deltay*(*d2zdxdy)[i][j];
	}
	else if (Ny == 1 && Nx == 0){
		if (output == table.ykey) { return 1.0; }
		if (output == table.xkey) { return 0.0; }
		val = (*dzdy)[i][j] + deltay*(*d2zdy2)[i][j] + deltax*(*d2zdxdy)[i][j];
	}
	else{
		throw NotImplementedError("only first derivatives currently supported");
	}
    return val;
}

#endif // !defined(NO_TABULAR_BACKENDS)