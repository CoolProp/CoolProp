#include "ExcessHEFunction.h"
#include "mixture_excess_term_JSON.h"

namespace CoolProp{

static MixtureExcessHELibrary mixtureexcesslibrary;

MixtureExcessHELibrary::MixtureExcessHELibrary()
{
    rapidjson::Document dd;

    dd.Parse<0>(mixture_excess_term_JSON.c_str());
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
            std::vector<std::string> CAS1V, Name1V, CAS2V, Name2V;

            if ((itr_coeff->HasMember("Name1") && (*itr_coeff)["Name1"].IsString())
                && (itr_coeff->HasMember("Name2") && (*itr_coeff)["Name2"].IsString())
                && (itr_coeff->HasMember("CAS1") && (*itr_coeff)["CAS1"].IsString())
                && (itr_coeff->HasMember("CAS2") && (*itr_coeff)["CAS2"].IsString()))
            {
                // Turn the strings into one-element vectors
                CAS1V = std::vector<std::string>(1,cpjson::get_string(*itr_coeff,"CAS1"));
                Name1V = std::vector<std::string>(1,cpjson::get_string(*itr_coeff,"Name1"));
                CAS2V = std::vector<std::string>(1,cpjson::get_string(*itr_coeff,"CAS2"));
                Name2V = std::vector<std::string>(1,cpjson::get_string(*itr_coeff,"Name2"));
            }
            else
            {
                CAS1V = cpjson::get_string_array((*itr_coeff)["CAS1"]);
                Name1V = cpjson::get_string_array((*itr_coeff)["Name1"]);
                CAS2V = cpjson::get_string_array((*itr_coeff)["CAS2"]);
                Name2V = cpjson::get_string_array((*itr_coeff)["Name2"]);
            }

            for (unsigned int i = 0; i < CAS1V.size(); ++i)
            {
                // Get the vector of CAS numbers 
                std::vector<std::string> CAS;
                CAS.resize(2);
                CAS[0] = CAS1V[i];
                CAS[1] = CAS2V[i];
                
                // Sort the CAS number vector
                std::sort(CAS.begin(), CAS.end());

                // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
                Dictionary d;

                // Populate the dictionary with common terms
                d.add_string("model", model);
                d.add_string("name1", Name1V[i]);
                d.add_string("name2", Name2V[i]);
                d.add_string("bibtex", BibTeX);

                if (!model.compare("Kunz-JCED-2012"))
                {
                    parse_Kunz_JCED_2012(d, *itr_coeff, i);
                }
                else if (!model.compare("Lemmon-JPCRD-2004"))
                {
                    parse_Lemmon_JPCRD_2004(d, *itr_coeff);
                }
                else if (!model.compare("Lemmon-JPCRD-2000"))
                {
                    parse_Lemmon_JPCRD_2000(d, *itr_coeff);
                }
                else
                {
                    throw ValueError();
                }
                
                // If not in map, add new entry to map with dictionary
                if (excess_map.find(CAS) == excess_map.end())
                {
                    // One-element vector of the dictionary
                    std::vector<Dictionary> vd(1, d);
                    // Pair for adding to map
                    std::pair<std::vector<std::string>, std::vector<Dictionary> > p(CAS, vd);
                    // Add
                    excess_map.insert(p);
                }
                else // If already in map, add entry to the end of the vector
                {
                    // Append dictionary to listing
                    excess_map[CAS].push_back(d);
                }           
            }
        }
	} 
}

void ExcessTerm::construct(const std::vector<CoolPropFluid*> &components)
{
    std::string _model;
    N = components.size();

    F.resize(N, std::vector<double>(N, 0));
    DepartureFunctionMatrix.resize(N, std::vector<DepartureFunction*>(N, NULL));

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

            std::vector<Dictionary> & vd = mixtureexcesslibrary.excess_map[CAS];
            if (vd.size() != 1) { throw NotImplementedError(); }
            // Get a reference to the dictionary itself to save a few dereferences
            Dictionary &dic = vd[0];

            std::string model = dic.get_string("model");

            if (!model.compare("Kunz-JCED-2012"))
            {
                F[i][j] = dic.get_number("F");

                std::vector<double> n = dic.get_double_vector("n");
                std::vector<double> d = dic.get_double_vector("d");
                std::vector<double> t = dic.get_double_vector("t");
                // Terms for the gaussian
                std::vector<double> eta = dic.get_double_vector("eta");
                std::vector<double> epsilon = dic.get_double_vector("epsilon");
                std::vector<double> beta = dic.get_double_vector("beta");
                std::vector<double> gamma = dic.get_double_vector("gamma");
                int Npower = static_cast<int>(dic.get_number("Npower"));
                DepartureFunctionMatrix[i][j] = new GERG2008DepartureFunction(n,d,t,eta,epsilon,beta,gamma,Npower);
            }
            else if (!model.compare("Lemmon-JPCRD-2004") || !model.compare("Lemmon-JPCRD-2000"))
            {
                throw NotImplementedError();
            }
            else
            {
                throw ValueError();
            }
        }
    }
}

ExcessTerm::~ExcessTerm()
{
    std::size_t N = DepartureFunctionMatrix.size();
    if (DepartureFunctionMatrix.empty()){return;}
	for (unsigned int i = 0; i < N; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{	
            delete DepartureFunctionMatrix[i][j]; // It is safe to delete NULL
            DepartureFunctionMatrix[i][j] = NULL;
		}
	}
    DepartureFunctionMatrix.clear();
}
double ExcessTerm::alphar(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dTau(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dDelta(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dDelta(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dDelta2(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dTau2(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dTau2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dDelta_dTau(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dxi(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->alphar(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphardxidxj(double tau, double delta, const std::vector<long double> &x, unsigned int i, unsigned int j)
{
	if (i != j)
	{
		return F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
	}
	else
	{
		return 0;
	}
}
double ExcessTerm::d2alphar_dxi_dTau(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dxi_dDelta(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dDelta(tau,delta);
		}
	}
	return summer;
}

GERG2008DepartureFunction::GERG2008DepartureFunction(const std::vector<double> &n,const std::vector<double> &d,const std::vector<double> &t,
                                                     const std::vector<double> &eta,const std::vector<double> &epsilon,const std::vector<double> &beta,
                                                     const std::vector<double> &gamma, int Npower)
{
    
    /// Break up into power and gaussian terms
    {
        std::vector<long double> _n(n.begin(), n.begin()+Npower);
        std::vector<long double> _d(d.begin(), d.begin()+Npower);
        std::vector<long double> _t(t.begin(), t.begin()+Npower);
        std::vector<long double> _l(Npower, 0.0);
        phi1 = ResidualHelmholtzPower(_n, _d, _t, _l);
    }
    if (n.size() == Npower)
    {
        using_gaussian = false;
    }
    else
    {
        using_gaussian = true;
        std::vector<long double> _n(n.begin()+Npower,                   n.end());
        std::vector<long double> _d(d.begin()+Npower,                   d.end());
        std::vector<long double> _t(t.begin()+Npower,                   t.end());
        std::vector<long double> _eta(eta.begin()+Npower,             eta.end());
        std::vector<long double> _epsilon(epsilon.begin()+Npower, epsilon.end());
        std::vector<long double> _beta(beta.begin()+Npower,          beta.end());
        std::vector<long double> _gamma(gamma.begin()+Npower,       gamma.end());
        phi2 = ResidualHelmholtzGERG2008Gaussian(_n, _d, _t, _eta, _epsilon, _beta, _gamma);
    }
}

double GERG2008DepartureFunction::alphar(double tau, double delta)
{
	if (using_gaussian){
		return phi1.base(tau, delta) + phi2.base(tau, delta);
	}
	else{
		return phi1.base(tau, delta);
	}
}
double GERG2008DepartureFunction::dalphar_dDelta(double tau, double delta)
{
	if (using_gaussian){
		return phi1.dDelta(tau, delta) + phi2.dDelta(tau, delta);
	}
	else{
		return phi1.dDelta(tau, delta);
	}
}
double GERG2008DepartureFunction::d2alphar_dDelta2(double tau, double delta)
{
	if (using_gaussian){
		return phi1.dDelta2(tau, delta) + phi2.dDelta2(tau, delta);
	}
	else{
		return phi1.dDelta2(tau, delta);
	}
}
double GERG2008DepartureFunction::d2alphar_dDelta_dTau(double tau, double delta)
{
	if (using_gaussian){
		return phi1.dDelta_dTau(tau, delta) + phi2.dDelta_dTau(tau, delta);
	}
	else{
		return phi1.dDelta_dTau(tau, delta);
	}
}
double GERG2008DepartureFunction::dalphar_dTau(double tau, double delta)
{
	if (using_gaussian){
		return phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
	}
	else{
		return phi1.dTau(tau, delta);
	}
}
double GERG2008DepartureFunction::d2alphar_dTau2(double tau, double delta)
{
	if (using_gaussian){
		return phi1.dTau2(tau, delta) + phi2.dTau2(tau, delta);
	}
	else{
		return phi1.dTau2(tau, delta);
	}
}

} /* namespace CoolProp */
