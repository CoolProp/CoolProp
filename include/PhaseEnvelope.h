#ifndef PHASE_ENVELOPE_H
#define PHASE_ENVELOPE_H

#include "Exceptions.h"

#define PHASE_ENVELOPE_MATRICES X(K) X(lnK) X(x) X(y)
#define PHASE_ENVELOPE_VECTORS \
    X(T)                       \
    X(p)                       \
    X(lnT)                     \
    X(lnp)                     \
    X(rhomolar_liq)            \
    X(rhomolar_vap)            \
    X(lnrhomolar_liq)          \
    X(lnrhomolar_vap)          \
    X(hmolar_liq)              \
    X(hmolar_vap)              \
    X(smolar_liq)              \
    X(smolar_vap)              \
    X(Q)                       \
    X(cpmolar_liq)             \
    X(cpmolar_vap)             \
    X(cvmolar_liq)             \
    X(cvmolar_vap)             \
    X(viscosity_liq)           \
    X(viscosity_vap)           \
    X(conductivity_liq)        \
    X(conductivity_vap)        \
    X(speed_sound_vap)

namespace CoolProp {

/** \brief A data structure to hold the data for a phase envelope
 *
 */
class PhaseEnvelopeData
{
   public:
    bool TypeI;             ///< True if it is a Type-I mixture that has a phase envelope that looks like a pure fluid more or less
    bool built;             ///< True if the phase envelope has been constructed
    std::size_t iTsat_max,  ///< The index of the point corresponding to the maximum temperature for Type-I mixtures
      ipsat_max,            ///< The index of the point corresponding to the maximum pressure for Type-I mixtures
      icrit;                ///< The index of the point corresponding to the critical point

// Use X macros to auto-generate the variables;
// each will look something like: std::vector<double> T;
#define X(name) std::vector<double> name;
    PHASE_ENVELOPE_VECTORS
#undef X

// Use X macros to auto-generate the variables;
// each will look something like: std::vector<std::vector<double> > K;
#define X(name) std::vector<std::vector<double>> name;
    PHASE_ENVELOPE_MATRICES
#undef X

    PhaseEnvelopeData() : TypeI(false), built(false), iTsat_max(-1), ipsat_max(-1), icrit(-1) {}

    void resize(std::size_t N) {
        K.resize(N);
        lnK.resize(N);
        x.resize(N);
        y.resize(N);
    }
    void clear() {
/* Use X macros to auto-generate the clearing code; each will look something like: T.clear(); */
#define X(name) name.clear();
        PHASE_ENVELOPE_VECTORS
#undef X
#define X(name) name.clear();
        PHASE_ENVELOPE_MATRICES
#undef X
    }
    void insert_variables(const CoolPropDbl T, const CoolPropDbl p, const CoolPropDbl rhomolar_liq, const CoolPropDbl rhomolar_vap,
                          const CoolPropDbl hmolar_liq, const CoolPropDbl hmolar_vap, const CoolPropDbl smolar_liq, const CoolPropDbl smolar_vap,
                          const std::vector<CoolPropDbl>& x, const std::vector<CoolPropDbl>& y, std::size_t i) {
        std::size_t N = K.size();
        if (N == 0) {
            throw CoolProp::ValueError("Cannot insert variables in phase envelope since resize() function has not been called");
        }
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
        for (unsigned int j = 0; j < N; j++) {
            this->K[j].insert(this->K[j].begin() + i, y[j] / x[j]);
            this->lnK[j].insert(this->lnK[j].begin() + i, log(y[j] / x[j]));
            this->x[j].insert(this->x[j].begin() + i, x[j]);
            this->y[j].insert(this->y[j].begin() + i, y[j]);
        }
        if (rhomolar_liq > rhomolar_vap) {
            this->Q.insert(this->Q.begin() + i, 1);
        } else {
            this->Q.insert(this->Q.begin() + i, 0);
        }
    };
    void store_variables(const CoolPropDbl T, const CoolPropDbl p, const CoolPropDbl rhomolar_liq, const CoolPropDbl rhomolar_vap,
                         const CoolPropDbl hmolar_liq, const CoolPropDbl hmolar_vap, const CoolPropDbl smolar_liq, const CoolPropDbl smolar_vap,
                         const std::vector<CoolPropDbl>& x, const std::vector<CoolPropDbl>& y) {
        std::size_t N = K.size();
        if (N == 0) {
            throw CoolProp::ValueError("Cannot store variables in phase envelope since resize() function has not been called");
        }
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
        for (unsigned int i = 0; i < N; i++) {
            this->K[i].push_back(y[i] / x[i]);
            this->lnK[i].push_back(log(y[i] / x[i]));
            this->x[i].push_back(x[i]);
            this->y[i].push_back(y[i]);
        }
        if (rhomolar_liq > rhomolar_vap) {
            this->Q.push_back(1);
        } else {
            this->Q.push_back(0);
        }
    };
};

} /* namespace CoolProp */

#endif
