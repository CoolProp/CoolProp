#ifndef PHASE_ENVELOPE_H
#define PHASE_ENVELOPE_H

class PhaseEnvelopeData
{
public:
    std::vector< std::vector<long double> > K, lnK, x, y;
    std::vector<long double> T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap;
    void resize(std::size_t N)
    {
        K.resize(N);
        lnK.resize(N);
        x.resize(N);
        y.resize(N);
    }
    void store_variables(const long double T, 
                         const long double p, 
                         const long double rhomolar_liq, 
                         const long double rhomolar_vap,
                         const std::vector<long double> & x, 
                         const std::vector<long double> & y)
    {
        std::size_t N = K.size();
        this->p.push_back(p);
        this->T.push_back(T);
        this->lnT.push_back(log(T));
        this->lnp.push_back(log(p));
        this->rhomolar_liq.push_back(rhomolar_liq);
        this->rhomolar_vap.push_back(rhomolar_vap);
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

#endif