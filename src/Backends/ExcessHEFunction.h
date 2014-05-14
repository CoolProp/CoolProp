#ifndef EXCESSHE_FUNCTIONS_H
#define EXCESSHE_FUNCTIONS_H

#include <vector>
#include "../Fluids/CoolPropFluid.h"

namespace CoolProp{

typedef std::vector<std::vector<long double> > STLMatrix;

/// A container for the mixing parameters for CoolProp mixtures
/**

*/
class MixtureExcessHELibrary
{
public:
    /// Map from sorted pair of CAS numbers to excess term dictionary.
    std::map<std::vector<std::string>, std::vector<Dictionary> > excess_map;
    MixtureExcessHELibrary();

    /// Parse a term from GERG 2008
    void parse_Kunz_JCED_2012(Dictionary &d, rapidjson::Value &val, int i)
    {
        assert(val.HasMember("F"));
        if (val["F"].IsDouble())
        {
            d.add_number("F",val["F"].GetDouble());
        }
        else
        {
            std::vector<double> F = cpjson::get_double_array(val["F"]);
            assert(static_cast<std::size_t>(i) < F.size());
            d.add_number("F", F[i]);
        }
        d.add_number("Npower", cpjson::get_double(val,"Npower"));
        
        // Terms for the power
        d.add_double_vector("n", cpjson::get_double_array(val["n"]));
        d.add_double_vector("d", cpjson::get_double_array(val["d"]));
        d.add_double_vector("t", cpjson::get_double_array(val["t"]));
        // Terms for the gaussian
        d.add_double_vector("eta", cpjson::get_double_array(val["eta"]));
        d.add_double_vector("epsilon", cpjson::get_double_array(val["epsilon"]));
        d.add_double_vector("beta", cpjson::get_double_array(val["beta"]));
        d.add_double_vector("gamma", cpjson::get_double_array(val["gamma"]));

    };

    /// Parse a term from HFC mixtures
    void parse_Lemmon_JPCRD_2004(Dictionary &d, rapidjson::Value &val)
    {
    };

    /// Parse a term from Air
    void parse_Lemmon_JPCRD_2000(Dictionary &d, rapidjson::Value &val)
    {
    };
};

/*! 
The abstract base class for departure functions for the excess part of the Helmholtz energy
*/
class DepartureFunction
{
public:
	DepartureFunction(){};
	virtual ~DepartureFunction(){};
	
	/// The excess Helmholtz energy of the binary pair
	/// Pure-virtual function (must be implemented in derived class
	virtual double alphar(double tau, double delta) = 0;
	virtual double dalphar_dDelta(double tau, double delta) = 0;
	virtual double d2alphar_dDelta2(double tau, double delta) = 0;
	virtual double d2alphar_dDelta_dTau(double tau, double delta) = 0;
	virtual double dalphar_dTau(double tau, double delta) = 0;
	virtual double d2alphar_dTau2(double tau, double delta) = 0;
};

class ExcessTerm
{
public:
	unsigned int N;
	std::vector<std::vector<DepartureFunction*> > DepartureFunctionMatrix;
	std::vector<std::vector<double> > F;
	
    ExcessTerm(){};
    void construct(const std::vector<CoolPropFluid*> &components);
	~ExcessTerm();

	double alphar(double tau, double delta, const std::vector<long double> &x);
	double dalphar_dDelta(double tau, double delta, const std::vector<long double> &x);
	double d2alphar_dDelta2(double tau, double delta, const std::vector<long double> &x);
	double d2alphar_dDelta_dTau(double tau, double delta, const std::vector<long double> &x);
	double dalphar_dTau(double tau, double delta, const std::vector<long double> &x);
	double d2alphar_dTau2(double tau, double delta, const std::vector<long double> &x);
	double dalphar_dxi(double tau, double delta, const std::vector<long double> &x, unsigned int i);
	double d2alphardxidxj(double tau, double delta, const std::vector<long double> &x, unsigned int i, unsigned int j);
	double d2alphar_dxi_dTau(double tau, double delta, const std::vector<long double> &x, unsigned int i);
	double d2alphar_dxi_dDelta(double tau, double delta, const std::vector<long double> &x, unsigned int i);
};

class GERG2008DepartureFunction : public DepartureFunction
{
protected:
	bool using_gaussian;
	ResidualHelmholtzPower phi1;
	ResidualHelmholtzGERG2008Gaussian phi2;
public:
	GERG2008DepartureFunction(){};
    GERG2008DepartureFunction(const std::vector<double> &n,const std::vector<double> &d,const std::vector<double> &t,
                              const std::vector<double> &eta,const std::vector<double> &epsilon,const std::vector<double> &beta,
                              const std::vector<double> &gamma, int Npower);
	~GERG2008DepartureFunction(){};
	double alphar(double tau, double delta);
	double dalphar_dDelta(double tau, double delta);
	double d2alphar_dDelta_dTau(double tau, double delta);
	double dalphar_dTau(double tau, double delta);
	double d2alphar_dDelta2(double tau, double delta);
	double d2alphar_dTau2(double tau, double delta);
};

} /* namespace CoolProp */
#endif