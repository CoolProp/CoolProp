#ifndef PHASE_ENVELOPE_H
#define PHASE_ENVELOPE_H

#include "Exceptions.h"

namespace CoolProp{
    
class PhaseEnvelopeData
{
public:
    bool built; ///< True if the phase envelope has been constructed
    std::vector< std::vector<long double> > K, lnK, x, y;
    std::vector<long double> T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap;
    
    PhaseEnvelopeData(){ built = false; };
    
    void resize(std::size_t N)
    {
        K.resize(N);
        lnK.resize(N);
        x.resize(N);
        y.resize(N);
    }
    void clear(){
        T.clear(); p.clear(); lnT.clear(); lnp.clear(); rhomolar_liq.clear(); rhomolar_vap.clear(); 
        lnrhomolar_liq.clear(); lnrhomolar_vap.clear(); hmolar_liq.clear(); hmolar_vap.clear(); smolar_liq.clear(); smolar_vap.clear();
        K.clear(); lnK.clear(); x.clear(); y.clear();
    }
    void store_variables(const long double T, 
                         const long double p, 
                         const long double rhomolar_liq, 
                         const long double rhomolar_vap,
                         const long double hmolar_liq, 
                         const long double hmolar_vap,
                         const long double smolar_liq, 
                         const long double smolar_vap,
                         const std::vector<long double> & x, 
                         const std::vector<long double> & y)
    {
        std::size_t N = K.size();
        if (N==0){throw CoolProp::ValueError("Cannot store variables in phase envelope since resize() function has not been called");}
        this->p.push_back(p);
        this->T.push_back(T);
        this->lnT.push_back(log(T));
        this->lnp.push_back(log(p));
        this->rhomolar_liq.push_back(rhomolar_liq);
        this->rhomolar_vap.push_back(rhomolar_vap);
        this->hmolar_liq.push_back(hmolar_liq);
        this->hmolar_vap.push_back(hmolar_vap);
        this->smolar_liq.push_back(smolar_liq);
        this->smolar_vap.push_back(smolar_vap);
        this->lnrhomolar_liq.push_back(log(rhomolar_liq));
        this->lnrhomolar_vap.push_back(log(rhomolar_vap));
        for (unsigned int i = 0; i < N; i++)
        {
            this->K[i].push_back(y[i]/x[i]);
            this->lnK[i].push_back(log(y[i]/x[i]));
            this->x[i].push_back(x[i]);
            this->y[i].push_back(y[i]);
        }
    };
};


} /* namespace CoolProp */

#endif