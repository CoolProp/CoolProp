#include "Ancillaries.h"
#include "DataStructures.h"

namespace CoolProp{

SaturationAncillaryFunction::SaturationAncillaryFunction(rapidjson::Value &json_code)
{
	std::string type = cpjson::get_string(json_code,"type");
	if (!type.compare("rational_polynomial"))
	{
		num_coeffs = vec_to_eigen(cpjson::get_double_array(json_code["A"]));
		den_coeffs = vec_to_eigen(cpjson::get_double_array(json_code["B"]));
		max_abs_error = cpjson::get_double(json_code,"max_abs_error");
	}
	else
	{
		n = cpjson::get_double_array(json_code["n"]);
		t = cpjson::get_double_array(json_code["t"]);
		Tmin = cpjson::get_double(json_code,"Tmin");
		Tmax = cpjson::get_double(json_code,"Tmax");
		reducing_value = cpjson::get_double(json_code,"reducing_value");
		using_tau_r = cpjson::get_bool(json_code,"using_tau_r");
		T_r = cpjson::get_double(json_code,"T_r");    
	}   
	
	if (!type.compare("rational_polynomial"))
		this->type = TYPE_RATIONAL_POLYNOMIAL;
	else if (!type.compare("rhoLnoexp"))
		this->type = TYPE_NOT_EXPONENTIAL;
	else
		this->type = TYPE_EXPONENTIAL;
	this->N = n.size();
	s = n;
};
    
double SaturationAncillaryFunction::evaluate(double T)
{
	if (type == TYPE_NOT_SET)
	{
		throw ValueError(format("type not set"));
	}
	else if (type == TYPE_RATIONAL_POLYNOMIAL)
	{
		Polynomial2D poly;
		return poly.evaluate(num_coeffs, T)/poly.evaluate(den_coeffs, T);
	}
	else
	{
		double THETA = 1-T/T_r;

		for (std::size_t i = 0; i < N; ++i)
		{
			s[i] = n[i]*pow(THETA, t[i]);
		}
		double summer = std::accumulate(s.begin(), s.end(), 0.0);

		if (type == TYPE_NOT_EXPONENTIAL)
		{
			return reducing_value*(1+summer);
		}
		else
		{
			double tau_r_value;
			if (using_tau_r)
				tau_r_value = T_r/T;
			else
				tau_r_value = 1.0;
			return reducing_value*exp(tau_r_value*summer);
		}
	}
}
double SaturationAncillaryFunction::invert(double value)
{
	// Invert the ancillary curve to get the temperature as a function of the output variable
	// Define the residual to be driven to zero
	class solver_resid : public FuncWrapper1D
	{
	public:
		int other;
		SaturationAncillaryFunction *anc;
		long double T, value, r, current_value;

		solver_resid(SaturationAncillaryFunction *anc, long double value) : anc(anc), value(value){};

		double call(double T){
			this->T = T;
			current_value = anc->evaluate(T);
			r = current_value - value;
			return r;
		};
	};
	solver_resid resid(this, value);
	std::string errstring;

	try{
		return Brent(resid,Tmin,Tmax,DBL_EPSILON,1e-12,100,errstring);
	}
	catch(std::exception &e){
		return Secant(resid,Tmax, -0.01, 1e-12, 100, errstring);
	}
}

void MeltingLineVariables::set_limits(void)
{
	if (type == MELTING_LINE_SIMON_TYPE){
		MeltingLinePiecewiseSimonSegment &partmin = simon.parts[0];
		MeltingLinePiecewiseSimonSegment &partmax = simon.parts[simon.parts.size()-1];
		Tmin = partmin.T_0;
		Tmax = partmax.T_max;
		pmin = partmin.p_0;
		pmax = evaluate(iP, iT, Tmax);
	}
	else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE){
		MeltingLinePiecewisePolynomialInTrSegment &partmin = polynomial_in_Tr.parts[0];
		MeltingLinePiecewisePolynomialInTrSegment &partmax = polynomial_in_Tr.parts[polynomial_in_Tr.parts.size() - 1];
		Tmin = partmin.T_0;
		Tmax = partmax.T_max;
		pmin = partmin.p_0;
		pmax = evaluate(iP, iT, Tmax);
	}
	else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE){
		MeltingLinePiecewisePolynomialInThetaSegment &partmin = polynomial_in_Theta.parts[0];
		MeltingLinePiecewisePolynomialInThetaSegment &partmax = polynomial_in_Theta.parts[polynomial_in_Theta.parts.size()-1];
		Tmin = partmin.T_0;
		Tmax = partmax.T_max;
		pmin = partmin.p_0;
		pmax = evaluate(iP, iT, Tmax);
	}
	else{
		throw ValueError("only Simon supported now");
	}
}

long double MeltingLineVariables::evaluate(int OF, int GIVEN, long double value)
{
	if (type == MELTING_LINE_NOT_SET){throw ValueError("Melting line curve not set");}
	if (OF == iP && GIVEN == iT){
		long double T = value;
		if (type == MELTING_LINE_SIMON_TYPE){
			// Need to find the right segment
			for (std::size_t i = 0; i < simon.parts.size(); ++i){
				MeltingLinePiecewiseSimonSegment &part = simon.parts[i];
				if (is_in_closed_range(part.T_min, part.T_max, T)){
					return part.p_0 + part.a*(pow(T/part.T_0,part.c)-1);
				}
			}
			throw ValueError("unable to calculate melting line (p,T) for Simon curve");
		}
		else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE){
			// Need to find the right segment
			for (std::size_t i = 0; i < polynomial_in_Tr.parts.size(); ++i){
				MeltingLinePiecewisePolynomialInTrSegment &part = polynomial_in_Tr.parts[i];
				if (is_in_closed_range(part.T_min, part.T_max, T)){
					long double summer = 0;
					for (std::size_t i =0; i < part.a.size(); ++i){
						summer += part.a[i]*(pow(T/part.T_0,part.t[i])-1);
					}
					return part.p_0*(1+summer);
				}
			}
			throw ValueError("unable to calculate melting line (p,T) for polynomial_in_Tr curve");
		}
		else if (type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE){
			// Need to find the right segment
			for (std::size_t i = 0; i < polynomial_in_Theta.parts.size(); ++i){
				MeltingLinePiecewisePolynomialInThetaSegment &part = polynomial_in_Theta.parts[i];
				if (is_in_closed_range(part.T_min, part.T_max, T)){
					long double summer = 0;
					for (std::size_t i =0; i < part.a.size(); ++i){
						summer += part.a[i]*pow(T/part.T_0-1,part.t[i]);
					}
					return part.p_0*(1+summer);
				}
			}
			throw ValueError("unable to calculate melting line (p,T) for polynomial_in_Theta curve");
		}
		else{
			throw ValueError(format("Invalid melting line type [%d]",type));
		}
	}
	else{
		if (type == MELTING_LINE_SIMON_TYPE){
			// Need to find the right segment
			for (std::size_t i = 0; i < simon.parts.size(); ++i){
				MeltingLinePiecewiseSimonSegment &part = simon.parts[i];
				//  p = part.p_0 + part.a*(pow(T/part.T_0,part.c)-1);
				long double T = pow((value-part.p_0)/part.a+1,1/part.c)*part.T_0;
				if (T >= part.T_0 && T <= part.T_max){
					return T;
				}
			}
			throw ValueError("unable to calculate melting line (p,T) for Simon curve");
		}
		else if (type == MELTING_LINE_POLYNOMIAL_IN_TR_TYPE || type == MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE)
		{
			class solver_resid : public FuncWrapper1D
			{
			public:

				MeltingLineVariables *line;
				long double r, given_p, calc_p, T;
				solver_resid(MeltingLineVariables *line, long double p) : line(line), given_p(p){};
				double call(double T){

					this->T = T;

					// Calculate p using melting line
					calc_p = line->evaluate(iP, iT, T);

					// Difference between the two is to be driven to zero
					r = given_p - calc_p;

					return r;
				};
			};
			solver_resid resid(this, value);
			
			double pmin = evaluate(iP, iT, Tmin);
			double pmax = evaluate(iP, iT, Tmax);
			std::string errstr;
			return Brent(resid, Tmin, Tmax, DBL_EPSILON, 1e-12, 100, errstr);
		}
		else{
			throw ValueError(format("Invalid melting line type (T,p) [%d]",type));
		}
	}
}

}; /* namespace CoolProp */
