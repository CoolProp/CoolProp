#include "MixtureParameters.h"
#include "mixture_departure_functions_JSON.h" // Creates the variable mixture_departure_functions_JSON
#include "mixture_binary_pairs_JSON.h" // Creates the variable mixture_binary_pairs_JSON
#include "predefined_mixtures_JSON.h" // Makes a std::string variable called predefined_mixtures_JSON

namespace CoolProp{

/** \brief A library of predefined mixtures
 *
 * Each entry in the predefined mixture library contains the names and mole fractions for the binary pairs
 */
class PredefinedMixturesLibrary{
    public:
    std::map<std::string, Dictionary> predefined_mixture_map;

    PredefinedMixturesLibrary(){
        rapidjson::Document doc;

        doc.Parse<0>(predefined_mixtures_JSON.c_str());
        if (doc.HasParseError()){throw ValueError();}

        // Iterate over the papers in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr)
        {
            // Instantiate the empty dictionary to be filled
            Dictionary dict;
            // Get the name
            std::string name = cpjson::get_string(*itr, "name")+".mix";
            // Get the fluid names
            dict.add_string_vector("fluids", cpjson::get_string_array(*itr, "fluids"));
            // Get the mole fractions
            dict.add_double_vector("mole_fractions", cpjson::get_double_array(*itr,"mole_fractions"));

            predefined_mixture_map.insert(std::pair<std::string, Dictionary >(name, dict));
        }
    }
};
static PredefinedMixturesLibrary predefined_mixtures_library;

std::string get_csv_predefined_mixtures()
{
    std::vector<std::string> out;
    for (std::map< std::string, Dictionary >::const_iterator it = predefined_mixtures_library.predefined_mixture_map.begin(); it != predefined_mixtures_library.predefined_mixture_map.end(); ++it)
    {
        out.push_back(it->first);
    }
    return strjoin(out, ",");
}

bool is_predefined_mixture(const std::string &name, Dictionary &dict){
    std::map<std::string, Dictionary>::const_iterator iter = predefined_mixtures_library.predefined_mixture_map.find(name);
    if (iter != predefined_mixtures_library.predefined_mixture_map.end()){
        dict = iter->second;
        return true;
    } else { return false; }
}

/** \brief A library of binary pair parameters for the mixture
 *
 * Each entry in the binary pair library includes reducing parameters as well as the name of the reducing function to be used and
 */
class MixtureBinaryPairLibrary{
public:
    /// Map from sorted pair of CAS numbers to reducing parameter map.  The reducing parameter map is a map from key (string) to value (double)
    std::map< std::vector<std::string>, std::vector<Dictionary> > binary_pair_map;

    /** \brief Construct the binary pair library including all the binary pairs that are possible
     *
     * The data structure also includes space for a string that gives the pointer to the departure function to be used for this binary pair.
     */
    MixtureBinaryPairLibrary()
    {
        rapidjson::Document doc;

        doc.Parse<0>(mixture_binary_pairs_JSON.c_str());
        if (doc.HasParseError()){throw ValueError();}

        // Iterate over the papers in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr)
        {
            // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
            Dictionary dict;

            // Get the vector of CAS numbers
            std::vector<std::string> CAS;
            CAS.push_back(cpjson::get_string(*itr, "CAS1"));
            CAS.push_back(cpjson::get_string(*itr, "CAS2"));
            std::string name1 = cpjson::get_string(*itr, "Name1");
            std::string name2 = cpjson::get_string(*itr, "Name2");

            // Sort the CAS number vector
            std::sort(CAS.begin(), CAS.end());

            // A sort was carried out, names/CAS were swapped
            bool swapped = CAS[0].compare(cpjson::get_string(*itr, "CAS1")) != 0;

            if (swapped){ std::swap(name1, name2); }

            // Populate the dictionary with common terms
            dict.add_string("name1", name1);
            dict.add_string("name2", name2);
            dict.add_string("BibTeX", cpjson::get_string(*itr, "BibTeX"));
            dict.add_number("F", cpjson::get_double(*itr, "F"));
            if (std::abs(dict.get_number("F")) > DBL_EPSILON){
                dict.add_string("function", cpjson::get_string(*itr, "function"));
            }

            if (itr->HasMember("xi") && itr->HasMember("zeta")){
                dict.add_string("type","Lemmon-xi-zeta");
                // Air and HFC mixtures from Lemmon - we could also directly do the conversion
                dict.add_number("xi", cpjson::get_double(*itr, "xi"));
                dict.add_number("zeta", cpjson::get_double(*itr, "zeta"));
            }
            else if (itr->HasMember("gammaT") && itr->HasMember("gammaV") && itr->HasMember("betaT") && itr->HasMember("betaV")){
                dict.add_string("type","GERG-2008");
                dict.add_number("gammaV", cpjson::get_double(*itr, "gammaV"));
                dict.add_number("gammaT", cpjson::get_double(*itr, "gammaT"));

                double betaV = cpjson::get_double(*itr, "betaV");
                double betaT = cpjson::get_double(*itr, "betaT");
                if (swapped){
                    dict.add_number("betaV", 1/betaV);
                    dict.add_number("betaT", 1/betaT);
                }
                else{
                    dict.add_number("betaV", betaV);
                    dict.add_number("betaT", betaT);
                }
            }
            else{
                std::cout << "Loading error: binary pair of " << name1 << " & " << name2 << "does not provide either a) xi and zeta b) gammaT, gammaV, betaT, and betaV" << std::endl;
                continue;
            }

            if (binary_pair_map.find(CAS) == binary_pair_map.end()){
                // Add to binary pair map by creating one-element vector
                binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary> >(CAS, std::vector<Dictionary>(1, dict)));
            }
            else
            {
                binary_pair_map[CAS].push_back(dict);
            }
        }
    }
    /// Add a simple mixing rule
    void add_simple_mixing_rule(const std::string &CAS1, const std::string &CAS2, const std::string &rule){
        // Get the empty dictionary to be filled by the appropriate reducing parameter filling function
        Dictionary dict;
        
        // Get the vector of CAS numbers
        std::vector<std::string> CAS;
        CAS.push_back(CAS1);
        CAS.push_back(CAS2);

        // Sort the CAS number vector
        std::sort(CAS.begin(), CAS.end());

        // Get the names of the compounds
        std::vector<std::string> names1(1,CAS[0]);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS1(new CoolProp::HelmholtzEOSMixtureBackend(names1));
        std::string name1 = HEOS1->name();
        std::vector<std::string> names2(1,CAS[1]);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS2(new CoolProp::HelmholtzEOSMixtureBackend(names2));
        std::string name2 = HEOS2->name();
        
        // Populate the dictionary with common terms
        dict.add_string("name1", name1);
        dict.add_string("name2", name2);
        dict.add_string("BibTeX", "N/A - LINEAR MIXING");
        dict.add_number("F", 0);
        dict.add_string("type","GERG-2008");
        
        if (rule == "linear"){
            // Terms for linear mixing
            dict.add_number("gammaT", 0.5*(HEOS1->T_critical()+HEOS2->T_critical())/sqrt(HEOS1->T_critical()*HEOS2->T_critical()));
            double rhoc1 = HEOS1->rhomolar_critical(), rhoc2 = HEOS2->rhomolar_critical();
            dict.add_number("gammaV", 4*(1/rhoc1+1/rhoc2)/pow(1/pow(rhoc1,1.0/3.0)+1/pow(rhoc2,1.0/3.0),3));
            dict.add_number("betaV", 1.0);
            dict.add_number("betaT", 1.0);
        }
        else if (rule == "Lorentz-Berthelot"){
            // Terms for Lorentz-Berthelot quadratic mixing
            dict.add_number("gammaT", 1.0);
            dict.add_number("gammaV", 1.0);
            dict.add_number("betaV", 1.0);
            dict.add_number("betaT", 1.0);
        }
        else{
            throw ValueError(format("Your simple mixing rule [%s] was not understood", rule.c_str()));
        }
        
        if (binary_pair_map.find(CAS) == binary_pair_map.end()){
            // Add to binary pair map by creating one-element vector
            binary_pair_map.insert(std::pair<std::vector<std::string>, std::vector<Dictionary> >(CAS, std::vector<Dictionary>(1, dict)));
        }
        else
        {
            binary_pair_map[CAS].push_back(dict);
        }
    }
};
// The modifiable parameter library
static MixtureBinaryPairLibrary mixturebinarypairlibrary;
// A fixed parameter library containing the default values
static MixtureBinaryPairLibrary mixturebinarypairlibrary_default;

/// Add a simple mixing rule
void apply_simple_mixing_rule(const std::string &CAS1, const std::string &CAS2, const std::string &rule){
    mixturebinarypairlibrary.add_simple_mixing_rule(CAS1, CAS2, rule);
}

std::string get_csv_mixture_binary_pairs()
{
    std::vector<std::string> out;
    for (std::map< std::vector<std::string>, std::vector<Dictionary> >::const_iterator it = mixturebinarypairlibrary.binary_pair_map.begin(); it != mixturebinarypairlibrary.binary_pair_map.end(); ++it)
    {
        out.push_back(strjoin(it->first, "&"));
    }
    return strjoin(out, ",");
}

std::string get_mixture_binary_pair_data(const std::string &CAS1, const std::string &CAS2, const std::string &key)
{
    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    if (mixturebinarypairlibrary.binary_pair_map.find(CAS) != mixturebinarypairlibrary.binary_pair_map.end()){
        std::vector<Dictionary> &v = mixturebinarypairlibrary.binary_pair_map[CAS];
        try{
            if (key == "name1"){ return v[0].get_string("name1"); }
            else if (key == "name2"){ return v[0].get_string("name2"); }
            else if (key == "BibTeX"){ return v[0].get_string("BibTeX"); }
            else if (key == "function"){ return v[0].get_string("function"); }
            else if (key == "type"){ return v[0].get_string("type"); }
            else if (key == "F"){ return format("%0.16g", v[0].get_double("F")); }
            else if (key == "xi"){ return format("%0.16g", v[0].get_double("xi")); }
            else if (key == "zeta"){ return format("%0.16g", v[0].get_double("zeta")); }
            else if (key == "gammaT"){ return format("%0.16g", v[0].get_double("gammaT")); }
            else if (key == "gammaV"){ return format("%0.16g", v[0].get_double("gammaV")); }
            else if (key == "betaT"){ return format("%0.16g", v[0].get_double("betaT")); }
            else if (key == "betaV"){ return format("%0.16g", v[0].get_double("betaV")); }
            else{ }
        }
        catch(...){ }
        throw ValueError(format("Could not match the parameter [%s] for the binary pair [%s,%s] - for now this is an error.", key.c_str(), CAS1.c_str(), CAS2.c_str()));
    }
    else{
        // Sort, see if other order works properly
        std::sort(CAS.begin(), CAS.end());
        if (mixturebinarypairlibrary.binary_pair_map.find(CAS) != mixturebinarypairlibrary.binary_pair_map.end())
        {
            throw ValueError(format("Could not match the binary pair [%s,%s] - order of CAS numbers is backwards; found the swapped CAS numbers.",CAS1.c_str(), CAS2.c_str()));
        }
        else{
            throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.",CAS1.c_str(), CAS2.c_str()));
        }
    }
}
void set_mixture_binary_pair_data(const std::string &CAS1, const std::string &CAS2, const std::string &key, const double value)
{
    // Find pair
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    if (mixturebinarypairlibrary.binary_pair_map.find(CAS) != mixturebinarypairlibrary.binary_pair_map.end()){
        std::vector<Dictionary> &v = mixturebinarypairlibrary.binary_pair_map[CAS];
        if (v[0].has_number(key)){
            v[0].add_number(key, value);
        }
        else{
            throw ValueError(format("Could not set the parameter [%s] for the binary pair [%s,%s] - for now this is an error", 
                                key.c_str(), CAS1.c_str(), CAS2.c_str()));
        }
    }
    else{
        // Sort, see if other order works properly
        std::sort(CAS.begin(), CAS.end());
        if (mixturebinarypairlibrary.binary_pair_map.find(CAS) != mixturebinarypairlibrary.binary_pair_map.end())
        {
            throw ValueError(format("Could not match the binary pair [%s,%s] - order of CAS numbers is backwards; found the swapped CAS numbers.",CAS1.c_str(), CAS2.c_str()));
        }
        else{
            throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.",CAS1.c_str(), CAS2.c_str()));
        }
    }
}

std::string get_reducing_function_name(const std::string &CAS1, const std::string &CAS2)
{
    std::vector<std::string> CAS;
    CAS.push_back(CAS1);
    CAS.push_back(CAS2);

    // Sort the CAS number vector - map is based on sorted CAS codes
    std::sort(CAS.begin(), CAS.end());

    if (mixturebinarypairlibrary.binary_pair_map.find(CAS) != mixturebinarypairlibrary.binary_pair_map.end()){
        return mixturebinarypairlibrary.binary_pair_map[CAS][0].get_string("function");
    }
    else{
        throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.",CAS1.c_str(), CAS2.c_str()));
    }
}

/** \brief A container for the departure functions for CoolProp mixtures
 */
class MixtureDepartureFunctionsLibrary{
public:
    /// Map from sorted pair of CAS numbers to departure term dictionary.
    std::map<std::string, Dictionary> departure_function_map;

    MixtureDepartureFunctionsLibrary()
    {
        rapidjson::Document doc;

        // Load the JSON data for the departure functions
        doc.Parse<0>(mixture_departure_functions_JSON.c_str());
        if (doc.HasParseError()){
            std::cout << mixture_departure_functions_JSON << std::endl ;
            throw ValueError("Unable to parse mixture_departure_functions_JSON.h");
        }

        // Iterate over the departure functions in the listing
        for (rapidjson::Value::ValueIterator itr = doc.Begin(); itr != doc.End(); ++itr)
        {
            // Get the empty dictionary to be filled in
            Dictionary dict;

            // Populate the dictionary with common terms
            std::string Name = cpjson::get_string(*itr, "Name");
            std::string type = cpjson::get_string(*itr, "type");
            dict.add_string("Name", Name);
            dict.add_string("BibTeX", cpjson::get_string(*itr, "BibTeX"));
            dict.add_string_vector("aliases", cpjson::get_string_array(*itr, "aliases"));
            dict.add_string("type", type);

            // Terms for the power (common to both types)
            dict.add_double_vector("n", cpjson::get_double_array(*itr, "n"));
            dict.add_double_vector("d", cpjson::get_double_array(*itr, "d"));
            dict.add_double_vector("t", cpjson::get_double_array(*itr, "t"));

            // Now we need to load additional terms
            if (!type.compare("GERG-2008")){
                // Number of terms that are power terms
                dict.add_number("Npower", cpjson::get_double(*itr, "Npower"));
                // Terms for the gaussian
                dict.add_double_vector("eta", cpjson::get_double_array(*itr, "eta"));
                dict.add_double_vector("epsilon", cpjson::get_double_array(*itr, "epsilon"));
                dict.add_double_vector("beta", cpjson::get_double_array(*itr, "beta"));
                dict.add_double_vector("gamma", cpjson::get_double_array(*itr, "gamma"));
            }
            else if (!type.compare("Exponential")){
                dict.add_double_vector("l", cpjson::get_double_array(*itr, "l"));
            }
            else{
                throw ValueError(format("It was not possible to parse departure function with type [%s]", type.c_str()));
            }

            // Check if this name is already in use
            if (departure_function_map.find(Name) == departure_function_map.end())
            {
                // Not in map, add new entry to map with dictionary as value and Name as key
                departure_function_map.insert(std::pair<std::string, Dictionary>(Name, dict));
            }
            else
            {
                // Error if already in map!
                //
                // Collect all the current names for departure functions for a nicer error message
                std::vector<std::string> names;
                for (std::map<std::string, Dictionary>::const_iterator it = departure_function_map.begin(); it != departure_function_map.end(); ++it)
                {
                    names.push_back(it->first);
                }
                throw ValueError(format("Name of departure function [%s] is already loaded. Current departure function names are: %s", Name.c_str(), strjoin(names,",").c_str() ));
            }
        }
    }
};
static MixtureDepartureFunctionsLibrary mixturedeparturefunctionslibrary;

void MixtureParameters::set_mixture_parameters(HelmholtzEOSMixtureBackend &HEOS)
{
    std::vector<CoolPropFluid> components = HEOS.get_components();

    std::size_t N = components.size();

    STLMatrix beta_v, gamma_v, beta_T, gamma_T;
    beta_v.resize(N, std::vector<CoolPropDbl>(N, 0));
    gamma_v.resize(N, std::vector<CoolPropDbl>(N, 0));
    beta_T.resize(N, std::vector<CoolPropDbl>(N, 0));
    gamma_T.resize(N, std::vector<CoolPropDbl>(N, 0));

    HEOS.residual_helmholtz->Excess.resize(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            if (i == j){ continue; }

            std::string CAS1 = components[i].CAS;
            std::vector<std::string> CAS(2,"");
            CAS[0] = components[i].CAS;
            CAS[1] = components[j].CAS;
            std::sort(CAS.begin(), CAS.end());

            // The variable swapped is true if a swap occured.
            bool swapped = (CAS1.compare(CAS[0]) != 0);

            // ***************************************************
            //         Reducing parameters for binary pair
            // ***************************************************

            if (mixturebinarypairlibrary.binary_pair_map.find(CAS) == mixturebinarypairlibrary.binary_pair_map.end())
            {
                throw ValueError(format("Could not match the binary pair [%s,%s] - for now this is an error.", CAS[0].c_str(), CAS[1].c_str()));
            }

            // Get a reference to the first matching binary pair in the dictionary
            Dictionary &dict_red = mixturebinarypairlibrary.binary_pair_map[CAS][0];

            // Get the name of the type being used, one of GERG-2008, Lemmon-xi-zeta, etc.
            std::string type_red = dict_red.get_string("type");

            if (!type_red.compare("GERG-2008")){
                if (swapped){
                    beta_v[i][j] = 1/dict_red.get_number("betaV");
                    beta_T[i][j] = 1/dict_red.get_number("betaT");
                }
                else{
                    beta_v[i][j] = dict_red.get_number("betaV");
                    beta_T[i][j] = dict_red.get_number("betaT");
                }
                gamma_v[i][j] = dict_red.get_number("gammaV");
                gamma_T[i][j] = dict_red.get_number("gammaT");
            }
            else if (!type_red.compare("Lemmon-xi-zeta")){
                LemmonAirHFCReducingFunction::convert_to_GERG(components,i,j,dict_red,beta_T[i][j],beta_v[i][j],gamma_T[i][j],gamma_v[i][j]);
            }
            else{
                throw ValueError(format("type [%s] for reducing function for pair [%s, %s] is invalid", type_red.c_str(),
                                                                                                        dict_red.get_string("Name1").c_str(),
                                                                                                        dict_red.get_string("Name2").c_str()   ));
            }
            /*
            if (i == 0){
                std::cout << format("betaT %10.9Lg gammaT %10.9Lg betaV %10.9Lg gammaV %10.9Lg %s",  
                                    beta_T[i][j], gamma_T[i][j], beta_v[i][j], gamma_v[i][j], get_mixture_binary_pair_data(CAS[0],CAS[1],"gammaT").c_str()) << std::endl;
            }
            */

            // ***************************************************
            //     Departure functions used in excess term
            // ***************************************************

            // Set the scaling factor F for the excess term
            HEOS.residual_helmholtz->Excess.F[i][j] = dict_red.get_number("F");

            if (std::abs(HEOS.residual_helmholtz->Excess.F[i][j]) < DBL_EPSILON){
                // Empty departure function that will just return 0
                std::vector<double> n(1,0), d(1,1), t(1,1), l(1,0);
                HEOS.residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(new ExponentialDepartureFunction(n, d, t, l));
                continue;
            }

            // Get the name of the departure function to be used for this binary pair
            std::string Name = CoolProp::get_reducing_function_name(components[i].CAS, components[j].CAS);

            // Get the dictionary itself
            Dictionary &dict_dep = mixturedeparturefunctionslibrary.departure_function_map[Name];

            if (dict_dep.is_empty()){throw ValueError(format("Departure function name [%s] seems to be invalid",Name.c_str()));}

            // These terms are common
            std::vector<double> n = dict_dep.get_double_vector("n");
            std::vector<double> d = dict_dep.get_double_vector("d");
            std::vector<double> t = dict_dep.get_double_vector("t");

            std::string type_dep = dict_dep.get_string("type");

            if (!type_dep.compare("GERG-2008")){
                // Number of power terms needed
                int Npower = static_cast<int>(dict_dep.get_number("Npower"));
                // Terms for the gaussian
                std::vector<double> eta = dict_dep.get_double_vector("eta");
                std::vector<double> epsilon = dict_dep.get_double_vector("epsilon");
                std::vector<double> beta = dict_dep.get_double_vector("beta");
                std::vector<double> gamma = dict_dep.get_double_vector("gamma");
                HEOS.residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(new GERG2008DepartureFunction(n, d, t, eta, epsilon, beta, gamma, Npower));
            }
            else if (!type_dep.compare("Exponential"))
            {
                // Powers of the exponents inside the exponential term
                std::vector<double> l = dict_dep.get_double_vector("l");
                HEOS.residual_helmholtz->Excess.DepartureFunctionMatrix[i][j].reset(new ExponentialDepartureFunction(n, d, t, l));
            }
            else
            {
                throw ValueError();
            }
        }
    }
    // We have obtained all the parameters needed for the reducing function, now set the reducing function for the mixture
    HEOS.Reducing = shared_ptr<ReducingFunction>(new GERG2008ReducingFunction(components, beta_v, gamma_v, beta_T, gamma_T));
}

} /* namespace CoolProp */