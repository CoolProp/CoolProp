#include "ReducingFunctions.h"
#include "mixture_reducing_parameters_JSON.h" // Creates the variable mixture_reducing_parameters_JSON
namespace CoolProp{

static MixingParameterLibrary mixturelibrary;

MixingParameterLibrary::MixingParameterLibrary()
{
    rapidjson::Document dd;

    dd.Parse<0>(mixture_reducing_parameters_JSON.c_str());
    if (dd.HasParseError()){throw ValueError();}

    // Iterate over the papers in the listing
    for (rapidjson::Value::ValueIterator itr = dd.Begin(); itr != dd.End(); ++itr)
	{
        // Iterate over the coeffs in the entry
        rapidjson::Value &Coeffs = (*itr)["Coeffs"];
        std::string model = cpjson::get_string(*itr, "Model");
        std::string BibTeX = cpjson::get_string(*itr, "BibTeX");
        for (rapidjson::Value::ValueIterator itr_coeff = Coeffs.Begin(); itr_coeff != Coeffs.End(); ++itr_coeff)
	    {
            // Get the vector of CAS numbers 
            std::vector<std::string> CAS;
            CAS.resize(2);
            CAS[0] = cpjson::get_string(*itr_coeff, "CAS1");
            CAS[1] = cpjson::get_string(*itr_coeff, "CAS2");
            std::string name1 = cpjson::get_string(*itr_coeff, "Name1");
            std::string name2 = cpjson::get_string(*itr_coeff, "Name2");
                
            // Sort the CAS number vector
            std::sort(CAS.begin(), CAS.end());

            // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
            Dictionary d;

            // A sort was carried out, names/CAS were swapped
            bool swapped = CAS[0].compare(cpjson::get_string(*itr_coeff, "CAS1")) != 0; 
            
            if (swapped){
                std::swap(name1,name2);
            }

            // Populate the dictionary with common terms
            d.add_string("model", model);
            d.add_string("name1", name1);
            d.add_string("name2", name2);
            d.add_string("bibtex", BibTeX);

            if (!model.compare("Kunz-JCED-2012"))
            {
                parse_Kunz_JCED_2012(d, *itr_coeff, swapped);
            }
            else if (!model.compare("Lemmon-JPCRD-2004"))
            {
                parse_Lemmon_JPCRD_2004(d, *itr_coeff, swapped);
            }
            else if (!model.compare("Lemmon-JPCRD-2000"))
            {
                parse_Lemmon_JPCRD_2000(d, *itr_coeff, swapped);
            }
            else
            {
                throw ValueError();
            }
                
            // If not in map, add new entry to map with dictionary
            if (reducing_map.find(CAS) == reducing_map.end())
            {
                // One-element vector of the dictionary
                std::vector<Dictionary> vd(1, d);
                // Pair for adding to map
                std::pair<std::vector<std::string>, std::vector<Dictionary> > p(CAS, vd);
                // Add
                reducing_map.insert(p);
            }
            else // If already in map, add entry to the end of the vector
            {
                // Append dictionary to listing
                reducing_map[CAS].push_back(d);
            }
        }
	} 
}

ReducingFunction *ReducingFunction::factory(const std::vector<CoolPropFluid*> &components)
{
    std::string _model;
    std::size_t N = components.size();

    STLMatrix beta_v, gamma_v, beta_T, gamma_T;
    beta_v.resize(N, std::vector<long double>(N, 0));
	gamma_v.resize(N, std::vector<long double>(N, 0));
	beta_T.resize(N, std::vector<long double>(N, 0));
	gamma_T.resize(N, std::vector<long double>(N, 0));

    for (unsigned int i = 0; i < N; ++i)
    {
        for (unsigned int j = 0; j < N; ++j)
        {
            if (i == j){ continue; }

            std::string CAS1 = components[i]->CAS;
            std::vector<std::string> CAS(2,"");
            CAS[0] = components[i]->CAS;
            CAS[1] = components[j]->CAS;
            std::sort(CAS.begin(), CAS.end());

            /// swapped is true if a swap occured.  
            bool swapped = (CAS1.compare(CAS[0]) != 0);

            std::vector<Dictionary> & vd = mixturelibrary.reducing_map[CAS];
            if (vd.size() != 1) { throw NotImplementedError(); }
            // Get a reference to the dictionary itself to save a few dereferences
            Dictionary &d = vd[0];

            std::string model = d.get_string("model");

            if (!model.compare("Kunz-JCED-2012"))
            {
                if (swapped)
                {
                    beta_v[i][j] = 1/d.get_number("betaV");
                    beta_T[i][j] = 1/d.get_number("betaT");
                }
                else
                {
                    beta_v[i][j] = d.get_number("betaV");
                    beta_T[i][j] = d.get_number("betaT");
                }
                gamma_v[i][j] = d.get_number("gammaV");
                gamma_T[i][j] = d.get_number("gammaT");
            }
            else if (!model.compare("Lemmon-JPCRD-2004") || !model.compare("Lemmon-JPCRD-2000"))
            {
                LemmonAirHFCReducingFunction::convert_to_GERG(components,i,j,d,beta_T[i][j],beta_v[i][j],gamma_T[i][j],gamma_v[i][j]);
            }
            else
            {
                throw ValueError();
            }
        }
    }
    return new GERG2008ReducingFunction(components,beta_v, gamma_v, beta_T, gamma_T);
}

long double ReducingFunction::d_ndTrdni_dxj__constxi(const std::vector<long double> &x, int i, int j)
{
	long double s = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		s += x[k]*d2Trdxidxj(x,j,k);
	}
	return d2Trdxidxj(x,i,j)-dTrdxi__constxj(x,j)-s;
}
long double ReducingFunction::d_ndrhorbardni_dxj__constxi(const std::vector<long double> &x, int i, int j)
{
	long double s = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		s += x[k]*d2rhormolardxidxj(x,j,k);
	}
	return d2rhormolardxidxj(x,j,i)-drhormolardxi__constxj(x,j)-s;
}
long double ReducingFunction::ndrhorbardni__constnj(const std::vector<long double> &x, int i)
{
	long double summer_term1 = 0;
	for (unsigned int j = 0; j < N; j++)
	{
		summer_term1 += x[j]*drhormolardxi__constxj(x,j);
	}
	return drhormolardxi__constxj(x,i)-summer_term1;
}
long double ReducingFunction::ndTrdni__constnj(const std::vector<long double> &x, int i)
{
	// GERG Equation 7.54
	long double summer_term1 = 0;
	for (unsigned int j = 0; j < N; j++)
	{
		summer_term1 += x[j]*dTrdxi__constxj(x,j);
	}
	return dTrdxi__constxj(x,i)-summer_term1;
}

long double GERG2008ReducingFunction::Tr(const std::vector<long double> &x)
{
	long double Tr = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = x[i], Tci = pFluids[i]->pEOS->reduce.T;
		Tr += xi*xi*Tci;
		
		// The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
		if (i==N-1){ break; }

		for (unsigned int j = i+1; j < N; j++)
		{
			Tr += c_Y_ij(i, j, beta_T, gamma_T, T_c)*f_Y_ij(x, i, j, beta_T);
		}
	}
	return Tr;
}
long double GERG2008ReducingFunction::dTrdxi__constxj(const std::vector<long double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	long double xi = x[i];
	long double dTr_dxi = 2*xi*pFluids[i]->pEOS->reduce.T;
	for (int k = 0; k < i; k++)
	{
		dTr_dxi += c_Y_ji(k,i,beta_T,gamma_T,T_c)*dfYkidxi__constxk(x,k,i,beta_T);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		dTr_dxi += c_Y_ij(i,k,beta_T,gamma_T,T_c)*dfYikdxi__constxk(x,i,k,beta_T);
	}
	return dTr_dxi;
}
long double GERG2008ReducingFunction::d2Trdxi2__constxj(const std::vector<long double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	long double d2Tr_dxi2 = 2*pFluids[i]->pEOS->reduce.T;
	for (int k = 0; k < i; k++)
	{
		d2Tr_dxi2 += c_Y_ij(k,i,beta_T,gamma_T,T_c)*d2fYkidxi2__constxk(x,k,i,beta_T);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		d2Tr_dxi2 += c_Y_ij(i,k,beta_T,gamma_T,T_c)*d2fYikdxi2__constxk(x,i,k,beta_T);
	}
	return d2Tr_dxi2;
}
long double GERG2008ReducingFunction::d2Trdxidxj(const std::vector<long double> &x, int i, int j)
{
	if (i == j)
	{
		return d2Trdxi2__constxj(x,i);
	}
	else
	{
		// See Table B9 from Kunz Wagner 2012 (GERG 2008)
		return c_Y_ij(i, j, beta_T, gamma_T, T_c)*d2fYijdxidxj(x, i, j, beta_T);
	}
}
long double GERG2008ReducingFunction::dvrmolardxi__constxj(const std::vector<long double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	long double xi = x[i];
	long double dvrbar_dxi = 2*xi/pFluids[i]->pEOS->reduce.rhomolar;

	for (int k = 0; k < i; k++)
	{
		dvrbar_dxi += c_Y_ij(k, i, beta_v, gamma_v, v_c)*dfYkidxi__constxk(x, k, i, beta_v);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		dvrbar_dxi += c_Y_ij(i, k, beta_v, gamma_v, v_c)*dfYikdxi__constxk(x, i, k, beta_v);	
	}
	return dvrbar_dxi;
}
long double GERG2008ReducingFunction::d2vrmolardxidxj(const std::vector<long double> &x, int i, int j)
{
	if (i == j)
	{
		return d2vrmolardxi2__constxj(x, i);
	}
	else
	{
		return c_Y_ij(i, j, beta_v, gamma_v, v_c)*d2fYijdxidxj(x, i, j, beta_v);
	}
}
long double GERG2008ReducingFunction::drhormolardxi__constxj(const std::vector<long double> &x, int i)
{
	return -pow(rhormolar(x),2)*dvrmolardxi__constxj(x,i);
}
long double GERG2008ReducingFunction::d2vrmolardxi2__constxj(const std::vector<long double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double d2vrbardxi2 = 2/pFluids[i]->pEOS->reduce.rhomolar;

	for (int k = 0; k < i; k++)
	{
		d2vrbardxi2 += c_Y_ij(k, i, beta_v, gamma_v, v_c)*d2fYkidxi2__constxk(x, k, i, beta_v);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		d2vrbardxi2 += c_Y_ij(i, k, beta_v, gamma_v, v_c)*d2fYikdxi2__constxk(x, i, k, beta_v);	
	}
	return d2vrbardxi2;
}
long double GERG2008ReducingFunction::d2rhormolardxi2__constxj(const std::vector<long double> &x, int i)
{
	long double rhor = this->rhormolar(x);
	long double dvrbardxi = this->dvrmolardxi__constxj(x,i);
	return 2*pow(rhor,(int)3)*pow(dvrbardxi,(int)2)-pow(rhor,(int)2)*this->d2vrmolardxi2__constxj(x,i);
}
long double GERG2008ReducingFunction::d2rhormolardxidxj(const std::vector<long double> &x, int i, int j)
{
	double rhor = this->rhormolar(x);
	double dvrbardxi = this->dvrmolardxi__constxj(x,i);
	double dvrbardxj = this->dvrmolardxi__constxj(x,j);
	return 2*pow(rhor,(int)3)*dvrbardxi*dvrbardxj-pow(rhor,(int)2)*this->d2vrmolardxidxj(x,i,j);
}

long double GERG2008ReducingFunction::rhormolar(const std::vector<long double> &x)
{
	double vrbar = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = x[i];
		vrbar += xi*xi/pFluids[i]->pEOS->reduce.rhomolar;

		if (i == N-1){ break; }

		for (unsigned int j = i+1; j < N; j++)
		{	
			vrbar += c_Y_ij(i, j, beta_v, gamma_v, v_c)*f_Y_ij(x, i, j, beta_v);
		}
	}
	return 1/vrbar;
}
long double GERG2008ReducingFunction::dfYkidxi__constxk(const std::vector<long double> &x, int k, int i, std::vector< std::vector< long double> > &beta)
{
	double xk = x[k], xi = x[i], beta_Y = beta[k][i];
	return xk*(xk+xi)/(beta_Y*beta_Y*xk+xi)+xk*xi/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi));
}
long double GERG2008ReducingFunction::dfYikdxi__constxk(const std::vector<long double> &x, int i, int k, std::vector< std::vector< long double> > &beta)
{
	double xk = x[k], xi = x[i], beta_Y = beta[i][k];
	return xk*(xi+xk)/(beta_Y*beta_Y*xi+xk)+xi*xk/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk));
}
long double GERG2008ReducingFunction::c_Y_ij(int i, int j, std::vector< std::vector< long double> > &beta, std::vector< std::vector< long double> > &gamma, std::vector< std::vector< long double> > &Y_c)
{
	return 2*beta[i][j]*gamma[i][j]*Y_c[i][j];
}
long double GERG2008ReducingFunction::c_Y_ji(int j, int i, std::vector< std::vector< long double> > &beta, std::vector< std::vector< long double> > &gamma, std::vector< std::vector< long double> > &Y_c)
{
	return 2/beta[i][j]*gamma[i][j]*Y_c[i][j];
}
long double GERG2008ReducingFunction::f_Y_ij(const std::vector<long double> &x, int i, int j, std::vector< std::vector< long double> > &beta)
{
	double xi = x[i], xj = x[j], beta_Y = beta[i][j];
	return xi*xj*(xi+xj)/(beta_Y*beta_Y*xi+xj);
}
long double GERG2008ReducingFunction::d2fYikdxi2__constxk(const std::vector<long double> &x, int i, int k, std::vector< std::vector< long double> > &beta)
{
	double xi = x[i], xk = x[k], beta_Y = beta[i][k];
	return 1/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk))*(2*xk-xi*xk*2*beta_Y*beta_Y/(beta_Y*beta_Y*xi+xk));
}
long double GERG2008ReducingFunction::d2fYkidxi2__constxk(const std::vector<long double> &x, int k, int i, std::vector< std::vector< long double> > &beta)
{
	double xi = x[i], xk = x[k], beta_Y = beta[k][i];
	return 1/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi))*(2*xk-xk*xi*2/(beta_Y*beta_Y*xk+xi));
}
long double GERG2008ReducingFunction::d2fYijdxidxj(const std::vector<long double> &x, int i, int j, std::vector< std::vector< long double> > &beta)
{
	double xi = x[i], xj = x[j], beta_Y = beta[i][j], beta_Y2 = beta_Y*beta_Y;
	return (xi+xj)/(beta_Y2*xi+xj) + xj/(beta_Y2*xi+xj)*(1-(xi+xj)/(beta_Y2*xi+xj))
		+xi/(beta_Y2*xi+xj)*(1-beta_Y2*(xi+xj)/(beta_Y2*xi+xj))
		-xi*xj/pow(beta_Y2*xi+xj,(int)2)*(1+beta_Y2-2*beta_Y2*(xi+xj)/(beta_Y2*xi+xj));
}


} /* namespace CoolProp */
