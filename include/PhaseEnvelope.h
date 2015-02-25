#ifndef PHASE_ENVELOPE_H
#define PHASE_ENVELOPE_H

#include "Exceptions.h"

namespace CoolProp{
    
/** \brief A data structure to hold the data for a phase envelope
 * 
 */
class PhaseEnvelopeData
{
public:
    bool TypeI; ///< True if it is a Type-I mixture that has a phase envelope that looks like a pure fluid more or less
    bool built; ///< True if the phase envelope has been constructed
    std::size_t iTsat_max, ///< The index of the point corresponding to the maximum temperature for Type-I mixtures
                ipsat_max, ///< The index of the point corresponding to the maximum pressure for Type-I mixtures
                icrit; ///< The index of the point corresponding to the critical point
                
    std::vector< std::vector<CoolPropDbl> > K, lnK, x, y;
    std::vector<CoolPropDbl> T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q;
    
    PhaseEnvelopeData(){ built = false; TypeI = false; };
    
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
        K.clear(); lnK.clear(); x.clear(); y.clear(); Q.clear();
    }
    void insert_variables(const CoolPropDbl T, 
                          const CoolPropDbl p, 
                          const CoolPropDbl rhomolar_liq, 
                          const CoolPropDbl rhomolar_vap,
                          const CoolPropDbl hmolar_liq, 
                          const CoolPropDbl hmolar_vap,
                          const CoolPropDbl smolar_liq, 
                          const CoolPropDbl smolar_vap,
                          const std::vector<CoolPropDbl> & x, 
                          const std::vector<CoolPropDbl> & y,
                          std::size_t i)
    {
        std::size_t N = K.size();
        if (N==0){throw CoolProp::ValueError("Cannot insert variables in phase envelope since resize() function has not been called");}
        this->p.insert(this->p.begin() + i, p);
        this->T.insert(this->T.begin() + i, T);
        this->lnT.insert(this->lnT.begin() + i, log(T));
        this->lnp.insert(this->lnp.begin() + i, log(p));
        this->rhomolar_liq.insert(this->rhomolar_liq.begin() + i, rhomolar_liq);
        this->rhomolar_vap.insert(this->rhomolar_vap.begin() + i, rhomolar_vap);
        this->hmolar_liq.insert(this->hmolar_liq.begin() + i, hmolar_liq);
        this->hmolar_vap.insert(this->hmolar_vap.begin() + i, hmolar_vap);
        this->smolar_liq.insert(this->smolar_liq.begin() + i, smolar_liq);
        this->smolar_vap.insert(this->smolar_vap.begin() + i, smolar_vap);
        this->lnrhomolar_liq.insert(this->lnrhomolar_liq.begin() + i, log(rhomolar_liq));
        this->lnrhomolar_vap.insert(this->lnrhomolar_vap.begin() + i, log(rhomolar_vap));
        for (unsigned int j = 0; j < N; j++)
        {
            this->K[j].insert(this->K[j].begin() + i, y[j]/x[j]);
            this->lnK[j].insert(this->lnK[j].begin() + i, log(y[j]/x[j]));
            this->x[j].insert(this->x[j].begin() + i, x[j]);
            this->y[j].insert(this->y[j].begin() + i, y[j]);
        }
        if (rhomolar_liq > rhomolar_vap){
            this->Q.insert(this->Q.begin(), 1);
        }
        else{
            this->Q.insert(this->Q.begin(), 0);
        }
    };
    void store_variables(const CoolPropDbl T, 
                         const CoolPropDbl p, 
                         const CoolPropDbl rhomolar_liq, 
                         const CoolPropDbl rhomolar_vap,
                         const CoolPropDbl hmolar_liq, 
                         const CoolPropDbl hmolar_vap,
                         const CoolPropDbl smolar_liq, 
                         const CoolPropDbl smolar_vap,
                         const std::vector<CoolPropDbl> & x, 
                         const std::vector<CoolPropDbl> & y)
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
        if (rhomolar_liq > rhomolar_vap){
            this->Q.push_back(1);
        }
        else{
            this->Q.push_back(0);
        }
    };
};


} /* namespace CoolProp */

#endif