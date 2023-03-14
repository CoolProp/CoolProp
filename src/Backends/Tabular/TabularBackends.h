#ifndef TABULAR_BACKENDS_H
#define TABULAR_BACKENDS_H

#include "AbstractState.h"
#include "CPmsgpack.h"
#include <msgpack/fbuffer.hpp>
#include "crossplatform_shared_ptr.h"
#include "Exceptions.h"
#include "CoolProp.h"
#include <sstream>
#include "Configuration.h"
#include "Backends/Helmholtz/PhaseEnvelopeRoutines.h"

/** ***MAGIC WARNING***!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_MATRICES \
    X(T)                 \
    X(p)                 \
    X(rhomolar)          \
    X(hmolar)            \
    X(smolar)            \
    X(umolar)            \
    X(dTdx)              \
    X(dTdy)              \
    X(dpdx)              \
    X(dpdy)              \
    X(drhomolardx)       \
    X(drhomolardy)       \
    X(dhmolardx)         \
    X(dhmolardy)         \
    X(dsmolardx)         \
    X(dsmolardy)         \
    X(dumolardx)         \
    X(dumolardy)         \
    X(d2Tdx2)            \
    X(d2Tdxdy)           \
    X(d2Tdy2)            \
    X(d2pdx2)            \
    X(d2pdxdy)           \
    X(d2pdy2)            \
    X(d2rhomolardx2)     \
    X(d2rhomolardxdy)    \
    X(d2rhomolardy2)     \
    X(d2hmolardx2)       \
    X(d2hmolardxdy)      \
    X(d2hmolardy2)       \
    X(d2smolardx2)       \
    X(d2smolardxdy)      \
    X(d2smolardy2)       \
    X(d2umolardx2)       \
    X(d2umolardxdy)      \
    X(d2umolardy2)       \
    X(visc)              \
    X(cond)

/** ***MAGIC WARNING***!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_SATURATION_VECTORS \
    X(TL)                          \
    X(pL)                          \
    X(logpL)                       \
    X(hmolarL)                     \
    X(smolarL)                     \
    X(umolarL)                     \
    X(rhomolarL)                   \
    X(logrhomolarL)                \
    X(viscL)                       \
    X(condL)                       \
    X(logviscL)                    \
    X(TV)                          \
    X(pV)                          \
    X(logpV)                       \
    X(hmolarV)                     \
    X(smolarV)                     \
    X(umolarV)                     \
    X(rhomolarV)                   \
    X(logrhomolarV)                \
    X(viscV)                       \
    X(condV)                       \
    X(logviscV)                    \
    X(cpmolarV)                    \
    X(cpmolarL)                    \
    X(cvmolarV)                    \
    X(cvmolarL)                    \
    X(speed_soundL)                \
    X(speed_soundV)

namespace CoolProp {

class PackablePhaseEnvelopeData : public PhaseEnvelopeData
{

   public:
    int revision;

    PackablePhaseEnvelopeData() : revision(0){};

    void copy_from_nonpackable(const PhaseEnvelopeData& PED) {
/* Use X macros to auto-generate the copying */
#define X(name) name = PED.name;
        PHASE_ENVELOPE_VECTORS
#undef X
/* Use X macros to auto-generate the copying */
#define X(name) name = PED.name;
        PHASE_ENVELOPE_MATRICES
#undef X
    };

    std::map<std::string, std::vector<double>> vectors;
    std::map<std::string, std::vector<std::vector<double>>> matrices;

    MSGPACK_DEFINE(revision, vectors, matrices);  // write the member variables that you want to pack using msgpack

    /// Take all the vectors that are in the class and pack them into the vectors map for easy unpacking using msgpack
    void pack() {
/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::string, std::vector<double> >("T", T)); */
#define X(name) vectors.insert(std::pair<std::string, std::vector<double>>(#name, name));
        PHASE_ENVELOPE_VECTORS
#undef X
/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::string, std::vector<std::vector<CoolPropDbl> > >("T", T)); */
#define X(name) matrices.insert(std::pair<std::string, std::vector<std::vector<double>>>(#name, name));
        PHASE_ENVELOPE_MATRICES
#undef X
    };
    std::map<std::string, std::vector<double>>::iterator get_vector_iterator(const std::string& name) {
        std::map<std::string, std::vector<double>>::iterator it = vectors.find(name);
        if (it == vectors.end()) {
            throw UnableToLoadError(format("could not find vector %s", name.c_str()));
        }
        return it;
    }
    std::map<std::string, std::vector<std::vector<double>>>::iterator get_matrix_iterator(const std::string& name) {
        std::map<std::string, std::vector<std::vector<double>>>::iterator it = matrices.find(name);
        if (it == matrices.end()) {
            throw UnableToLoadError(format("could not find matrix %s", name.c_str()));
        }
        return it;
    }
    /// Take all the vectors that are in the class and unpack them from the vectors map
    void unpack() {
/* Use X macros to auto-generate the unpacking code;
         * each will look something like: T = get_vector_iterator("T")->second
         */
#define X(name) name = get_vector_iterator(#name)->second;
        PHASE_ENVELOPE_VECTORS
#undef X
/* Use X macros to auto-generate the unpacking code;
         * each will look something like: T = get_matrix_iterator("T")->second
         **/
#define X(name) name = get_matrix_iterator(#name)->second;
        PHASE_ENVELOPE_MATRICES
#undef X
        // Find the index of the point with the highest temperature
        iTsat_max = std::distance(T.begin(), std::max_element(T.begin(), T.end()));
        // Find the index of the point with the highest pressure
        ipsat_max = std::distance(p.begin(), std::max_element(p.begin(), p.end()));
    };
    void deserialize(msgpack::object& deserialized) {
        PackablePhaseEnvelopeData temp;
        deserialized.convert(temp);
        temp.unpack();
        if (revision > temp.revision) {
            throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
        }
        std::swap(*this, temp);  // Swap if successful
    };
};

/// Get a conversion factor from mass to molar if needed
inline void mass_to_molar(parameters& param, double& conversion_factor, double molar_mass) {
    conversion_factor = 1.0;
    switch (param) {
        case iDmass:
            conversion_factor = molar_mass;
            param = iDmolar;
            break;
        case iHmass:
            conversion_factor /= molar_mass;
            param = iHmolar;
            break;
        case iSmass:
            conversion_factor /= molar_mass;
            param = iSmolar;
            break;
        case iUmass:
            conversion_factor /= molar_mass;
            param = iUmolar;
            break;
        case iCvmass:
            conversion_factor /= molar_mass;
            param = iCvmolar;
            break;
        case iCpmass:
            conversion_factor /= molar_mass;
            param = iCpmolar;
            break;
        case iDmolar:
        case iHmolar:
        case iSmolar:
        case iUmolar:
        case iCvmolar:
        case iCpmolar:
        case iT:
        case iP:
        case ispeed_sound:
        case iisothermal_compressibility:
        case iisobaric_expansion_coefficient:
        case iviscosity:
        case iconductivity:
            return;
        default:
            throw ValueError("TabularBackends::mass_to_molar - I don't know how to convert this parameter");
    }
}

/** \brief This class holds the data for a two-phase table that is log spaced in p
 *
 * It contains very few members or methods, mostly it just holds the data
 */
class PureFluidSaturationTableData
{
   public:
    std::size_t N;
    shared_ptr<CoolProp::AbstractState> AS;

    PureFluidSaturationTableData() {
        N = 1000;
        revision = 1;
    }

    /// Build this table
    void build(shared_ptr<CoolProp::AbstractState>& AS);

/* Use X macros to auto-generate the variables; each will look something like: std::vector<double> T; */
#define X(name) std::vector<double> name;
    LIST_OF_SATURATION_VECTORS
#undef X

    int revision;
    std::map<std::string, std::vector<double>> vectors;

    MSGPACK_DEFINE(revision, vectors);  // write the member variables that you want to pack

    /***
         * \brief Determine if a set of inputs are single-phase or inside the saturation table
         * @param main The main variable that is being provided (currently T or P)
         * @param mainval The value of the main variable that is being provided
         * @param other The secondary variable
         * @param val The value of the secondary variable
         * @param iL The index associated with the nearest point for the liquid
         * @param iV The index associated with the nearest point for the vapor
         * @param yL The value associated with the nearest point for the liquid (based on interpolation)
         * @param yV The value associated with the nearest point for the vapor (based on interpolation)

         \note If PQ or QT are inputs, yL and yV will correspond to the other main variable: p->T or T->p
         */
    bool is_inside(parameters main, double mainval, parameters other, double val, std::size_t& iL, std::size_t& iV, CoolPropDbl& yL,
                   CoolPropDbl& yV) {
        std::vector<double>*yvecL = NULL, *yvecV = NULL;
        switch (other) {
            case iT:
                yvecL = &TL;
                yvecV = &TV;
                break;
            case iHmolar:
                yvecL = &hmolarL;
                yvecV = &hmolarV;
                break;
            case iQ:
                yvecL = &TL;
                yvecV = &TV;
                break;
            case iSmolar:
                yvecL = &smolarL;
                yvecV = &smolarV;
                break;
            case iUmolar:
                yvecL = &umolarL;
                yvecV = &umolarV;
                break;
            case iDmolar:
                yvecL = &rhomolarL;
                yvecV = &rhomolarV;
                break;
            default:
                throw ValueError("invalid input for other in is_inside");
        }

        // Trivial checks
        if (main == iP) {
            // If p is outside the range (ptriple, pcrit), considered to not be inside
            double pmax = this->pV[pV.size() - 1], pmin = this->pV[0];
            if (mainval > pmax || mainval < pmin) {
                return false;
            }
        } else if (main == iT) {
            // If T is outside the range (Tmin, Tcrit), considered to not be inside
            double Tmax = this->TV[TV.size() - 1], Tmin = this->TV[0];
            if (mainval > Tmax || mainval < Tmin) {
                return false;
            }
        } else {
            throw ValueError("invalid input for other in is_inside");
        }

        // Now check based on a rough analysis using bounding pressure
        std::size_t iLplus, iVplus;
        // Find the indices (iL,iL+1) & (iV,iV+1) that bound the given pressure
        // In general iV and iL will be the same, but if pseudo-pure, they might
        // be different
        if (main == iP) {
            bisect_vector(pV, mainval, iV);
            bisect_vector(pL, mainval, iL);
        } else if (main == iT) {
            bisect_vector(TV, mainval, iV);
            bisect_vector(TL, mainval, iL);
        } else {
            throw ValueError(format("For now, main input in is_inside must be T or p"));
        }

        iVplus = std::min(iV + 1, N - 1);
        iLplus = std::min(iL + 1, N - 1);
        if (other == iQ) {
            // Actually do "saturation" call using cubic interpolation
            if (iVplus < 3) {
                iVplus = 3;
            }
            if (iLplus < 3) {
                iLplus = 3;
            }
            if (main == iP) {
                double logp = log(mainval);
                // Calculate temperature
                yV = CubicInterp(logpV, TV, iVplus - 3, iVplus - 2, iVplus - 1, iVplus, logp);
                yL = CubicInterp(logpL, TL, iLplus - 3, iLplus - 2, iLplus - 1, iLplus, logp);
            } else if (main == iT) {
                // Calculate pressure
                yV = exp(CubicInterp(TV, logpV, iVplus - 3, iVplus - 2, iVplus - 1, iVplus, mainval));
                yL = exp(CubicInterp(TL, logpL, iLplus - 3, iLplus - 2, iLplus - 1, iLplus, mainval));
            }
            return true;
        }
        // Find the bounding values for the other variable
        double ymin = min4((*yvecL)[iL], (*yvecL)[iLplus], (*yvecV)[iV], (*yvecV)[iVplus]);
        double ymax = max4((*yvecL)[iL], (*yvecL)[iLplus], (*yvecV)[iV], (*yvecV)[iVplus]);
        if (val < ymin || val > ymax) {
            return false;
        }
        // Actually do "saturation" call using cubic interpolation
        if (iVplus < 3) {
            iVplus = 3;
        }
        if (iLplus < 3) {
            iLplus = 3;
        }
        if (main == iP) {
            double logp = log(mainval);
            yV = CubicInterp(logpV, *yvecV, iVplus - 3, iVplus - 2, iVplus - 1, iVplus, logp);
            yL = CubicInterp(logpL, *yvecL, iLplus - 3, iLplus - 2, iLplus - 1, iLplus, logp);
        } else if (main == iT) {
            yV = CubicInterp(TV, *yvecV, iVplus - 3, iVplus - 2, iVplus - 1, iVplus, mainval);
            yL = CubicInterp(TL, *yvecL, iLplus - 3, iLplus - 2, iLplus - 1, iLplus, mainval);
        }

        if (!is_in_closed_range(yV, yL, static_cast<CoolPropDbl>(val))) {
            return false;
        } else {
            iL = iLplus - 1;
            iV = iVplus - 1;
            return true;
        }
    }
    /// Resize all the vectors
    void resize(std::size_t N) {
/* Use X macros to auto-generate the code; each will look something like: T.resize(N); std::fill(T.begin(), T.end(), _HUGE); */
#define X(name)     \
    name.resize(N); \
    std::fill(name.begin(), name.end(), _HUGE);
        LIST_OF_SATURATION_VECTORS
#undef X
    };
    /// Take all the vectors that are in the class and pack them into the vectors map for easy unpacking using msgpack
    void pack() {
/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
#define X(name) vectors.insert(std::pair<std::string, std::vector<double>>(#name, name));
        LIST_OF_SATURATION_VECTORS
#undef X
    };
    std::map<std::string, std::vector<double>>::iterator get_vector_iterator(const std::string& name) {
        std::map<std::string, std::vector<double>>::iterator it = vectors.find(name);
        if (it == vectors.end()) {
            throw UnableToLoadError(format("could not find vector %s", name.c_str()));
        }
        return it;
    }
    /// Take all the vectors that are in the class and unpack them from the vectors map
    void unpack() {
/* Use X macros to auto-generate the unpacking code; each will look something like: T = get_vector_iterator("T")->second */
#define X(name) name = get_vector_iterator(#name)->second;
        LIST_OF_SATURATION_VECTORS
#undef X
        N = TL.size();
    };
    void deserialize(msgpack::object& deserialized) {
        PureFluidSaturationTableData temp;
        deserialized.convert(temp);
        temp.unpack();
        if (N != temp.N) {
            throw ValueError(format("old [%d] and new [%d] sizes don't agree", temp.N, N));
        } else if (revision > temp.revision) {
            throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
        }
        std::swap(*this, temp);  // Swap
        this->AS = temp.AS;      // Reconnect the AbstractState pointer
    };
    double evaluate(parameters output, double p_or_T, double Q, std::size_t iL, std::size_t iV) {
        if (iL <= 2) {
            iL = 2;
        } else if (iL + 1 == N) {
            iL = N - 2;
        }
        if (iV <= 2) {
            iV = 2;
        } else if (iV + 1 == N) {
            iV = N - 2;
        }
        double logp = log(p_or_T);
        switch (output) {
            case iP: {
                double _logpV = CubicInterp(this->TV, logpV, iV - 2, iV - 1, iV, iV + 1, p_or_T);
                double _logpL = CubicInterp(this->TL, logpL, iL - 2, iL - 1, iL, iL + 1, p_or_T);
                return Q * exp(_logpV) + (1 - Q) * exp(_logpL);
            }
            case iT: {
                double TV = CubicInterp(logpV, this->TV, iV - 2, iV - 1, iV, iV + 1, logp);
                double TL = CubicInterp(logpL, this->TL, iL - 2, iL - 1, iL, iL + 1, logp);
                return Q * TV + (1 - Q) * TL;
            }
            case iSmolar: {
                double sV = CubicInterp(logpV, smolarV, iV - 2, iV - 1, iV, iV + 1, logp);
                double sL = CubicInterp(logpL, smolarL, iL - 2, iL - 1, iL, iL + 1, logp);
                return Q * sV + (1 - Q) * sL;
            }
            case iHmolar: {
                double hV = CubicInterp(logpV, hmolarV, iV - 2, iV - 1, iV, iV + 1, logp);
                double hL = CubicInterp(logpL, hmolarL, iL - 2, iL - 1, iL, iL + 1, logp);
                return Q * hV + (1 - Q) * hL;
            }
            case iUmolar: {
                double uV = CubicInterp(logpV, umolarV, iV - 2, iV - 1, iV, iV + 1, logp);
                double uL = CubicInterp(logpL, umolarL, iL - 2, iL - 1, iL, iL + 1, logp);
                return Q * uV + (1 - Q) * uL;
            }
            case iDmolar: {
                double rhoV = exp(CubicInterp(logpV, logrhomolarV, iV - 2, iV - 1, iV, iV + 1, logp));
                double rhoL = exp(CubicInterp(logpL, logrhomolarL, iL - 2, iL - 1, iL, iL + 1, logp));
                if (!ValidNumber(rhoV)) {
                    throw ValueError("rhoV is invalid");
                }
                if (!ValidNumber(rhoL)) {
                    throw ValueError("rhoL is invalid");
                }
                return 1 / (Q / rhoV + (1 - Q) / rhoL);
            }
            case iconductivity: {
                double kV = CubicInterp(logpV, condV, iV - 2, iV - 1, iV, iV + 1, logp);
                double kL = CubicInterp(logpL, condL, iL - 2, iL - 1, iL, iL + 1, logp);
                if (!ValidNumber(kV)) {
                    throw ValueError("kV is invalid");
                }
                if (!ValidNumber(kL)) {
                    throw ValueError("kL is invalid");
                }
                return Q * kV + (1 - Q) * kL;
            }
            case iviscosity: {
                double muV = exp(CubicInterp(logpV, logviscV, iV - 2, iV - 1, iV, iV + 1, logp));
                double muL = exp(CubicInterp(logpL, logviscL, iL - 2, iL - 1, iL, iL + 1, logp));
                if (!ValidNumber(muV)) {
                    throw ValueError("muV is invalid");
                }
                if (!ValidNumber(muL)) {
                    throw ValueError("muL is invalid");
                }
                return 1 / (Q / muV + (1 - Q) / muL);
            }
            case iCpmolar: {
                double cpV = CubicInterp(logpV, cpmolarV, iV - 2, iV - 1, iV, iV + 1, logp);
                double cpL = CubicInterp(logpL, cpmolarL, iL - 2, iL - 1, iL, iL + 1, logp);
                if (!ValidNumber(cpV)) {
                    throw ValueError("cpV is invalid");
                }
                if (!ValidNumber(cpL)) {
                    throw ValueError("cpL is invalid");
                }
                return Q * cpV + (1 - Q) * cpL;
            }
            case iCvmolar: {
                double cvV = CubicInterp(logpV, cvmolarV, iV - 2, iV - 1, iV, iV + 1, logp);
                double cvL = CubicInterp(logpL, cvmolarL, iL - 2, iL - 1, iL, iL + 1, logp);
                if (!ValidNumber(cvV)) {
                    throw ValueError("cvV is invalid");
                }
                if (!ValidNumber(cvL)) {
                    throw ValueError("cvL is invalid");
                }
                return Q * cvV + (1 - Q) * cvL;
            }
            case ispeed_sound: {
                double wV = CubicInterp(logpV, speed_soundV, iV - 2, iV - 1, iV, iV + 1, logp);
                double wL = CubicInterp(logpL, speed_soundL, iL - 2, iL - 1, iL, iL + 1, logp);
                if (!ValidNumber(wV)) {
                    throw ValueError("wV is invalid");
                }
                if (!ValidNumber(wL)) {
                    throw ValueError("wL is invalid");
                }
                return Q * wV + (1 - Q) * wL;
            }
            default:
                throw ValueError("Output variable for evaluate is invalid");
        }
    };
    /**
         *  @brief Calculate the first derivative ALONG a saturation curve
         * @param Of1 The parameter that the derivative is to be taken of
         * @param Wrt1 The parameter that the derivative is to be taken with respect to
         * @param Q The vapor quality, 0 or 1
         * @param val The value of the WRT parameter
         * @param i The index in the vectors to be used; must be > 2 and < len-2
         */
    double first_saturation_deriv(parameters Of1, parameters Wrt1, int Q, double val, std::size_t i) {
        if (i < 2 || i > TL.size() - 2) {
            throw ValueError(format("Invalid index (%d) to calc_first_saturation_deriv in TabularBackends", i));
        }
        std::vector<double>*x, *y;
        // Connect pointers for each vector
        switch (Wrt1) {
            case iT:
                x = (Q == 0) ? &TL : &TV;
                break;
            case iP:
                x = (Q == 0) ? &pL : &pV;
                break;
            default:
                throw ValueError(format("Key for Wrt1 is invalid in calc_first_saturation_deriv"));
        }
        CoolPropDbl factor = 1.0;
        switch (Of1) {
            case iT:
                y = (Q == 0) ? &TL : &TV;
                break;
            case iP:
                y = (Q == 0) ? &pL : &pV;
                break;
            case iDmolar:
                y = (Q == 0) ? &rhomolarL : &rhomolarV;
                break;
            case iHmolar:
                y = (Q == 0) ? &hmolarL : &hmolarV;
                break;
            case iSmolar:
                y = (Q == 0) ? &smolarL : &smolarV;
                break;
            case iUmolar:
                y = (Q == 0) ? &umolarL : &umolarV;
                break;
            case iDmass:
                y = (Q == 0) ? &rhomolarL : &rhomolarV;
                factor = AS->molar_mass();
                break;
            case iHmass:
                y = (Q == 0) ? &hmolarL : &hmolarV;
                factor = 1 / AS->molar_mass();
                break;
            case iSmass:
                y = (Q == 0) ? &smolarL : &smolarV;
                factor = 1 / AS->molar_mass();
                break;
            case iUmass:
                y = (Q == 0) ? &umolarL : &umolarV;
                factor = 1 / AS->molar_mass();
                break;
            default:
                throw ValueError(format("Key for Of1 is invalid in calc_first_saturation_deriv"));
        }
        return CubicInterpFirstDeriv((*x)[i - 2], (*x)[i - 1], (*x)[i], (*x)[i + 1], (*y)[i - 2], (*y)[i - 1], (*y)[i], (*y)[i + 1], val) * factor;
    };
    //calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant);
};

/** \brief This class holds the data for a single-phase interpolation table that is regularly spaced
 *
 * It contains very few members or methods, mostly it just holds the data
 */
class SinglePhaseGriddedTableData
{

   public:
    std::size_t Nx, Ny;
    CoolProp::parameters xkey, ykey;
    shared_ptr<CoolProp::AbstractState> AS;
    std::vector<double> xvec, yvec;
    std::vector<std::vector<std::size_t>> nearest_neighbor_i, nearest_neighbor_j;
    bool logx, logy;
    double xmin, ymin, xmax, ymax;

    virtual void set_limits() = 0;

    SinglePhaseGriddedTableData() {
        Nx = 200;
        Ny = 200;
        revision = 0;
        xkey = INVALID_PARAMETER;
        ykey = INVALID_PARAMETER;
        logx = false;
        logy = false;
        xmin = _HUGE;
        xmax = _HUGE;
        ymin = _HUGE;
        ymax = _HUGE;
    }

/* Use X macros to auto-generate the variables; each will look something like: std::vector< std::vector<double> > T; */
#define X(name) std::vector<std::vector<double>> name;
    LIST_OF_MATRICES
#undef X
    int revision;
    std::map<std::string, std::vector<std::vector<double>>> matrices;
    /// Build this table
    void build(shared_ptr<CoolProp::AbstractState>& AS);

    MSGPACK_DEFINE(revision, matrices, xmin, xmax, ymin, ymax);  // write the member variables that you want to pack
    /// Resize all the matrices
    void resize(std::size_t Nx, std::size_t Ny) {
/* Use X macros to auto-generate the code; each will look something like: T.resize(Nx, std::vector<double>(Ny, _HUGE)); */
#define X(name) name.resize(Nx, std::vector<double>(Ny, _HUGE));
        LIST_OF_MATRICES
#undef X
        make_axis_vectors();
    };
    /// Make vectors for the x-axis values and the y-axis values
    void make_axis_vectors(void) {
        if (logx) {
            xvec = logspace(xmin, xmax, Nx);
        } else {
            xvec = linspace(xmin, xmax, Nx);
        }
        if (logy) {
            yvec = logspace(ymin, ymax, Ny);
        } else {
            yvec = linspace(ymin, ymax, Ny);
        }
    };
    /// Make matrices of good neighbors if the current value for i,j corresponds to a bad node
    void make_good_neighbors(void) {
        nearest_neighbor_i.resize(Nx, std::vector<std::size_t>(Ny, std::numeric_limits<std::size_t>::max()));
        nearest_neighbor_j.resize(Nx, std::vector<std::size_t>(Ny, std::numeric_limits<std::size_t>::max()));
        for (std::size_t i = 0; i < xvec.size(); ++i) {
            for (std::size_t j = 0; j < yvec.size(); ++j) {
                nearest_neighbor_i[i][j] = i;
                nearest_neighbor_j[i][j] = j;
                if (!ValidNumber(T[i][j])) {
                    int xoffsets[] = {-1, 1, 0, 0, -1, 1, 1, -1};
                    int yoffsets[] = {0, 0, 1, -1, -1, -1, 1, 1};
                    // Length of offset
                    std::size_t N = sizeof(xoffsets) / sizeof(xoffsets[0]);
                    for (std::size_t k = 0; k < N; ++k) {
                        std::size_t iplus = i + xoffsets[k];
                        std::size_t jplus = j + yoffsets[k];
                        if (0 < iplus && iplus < Nx - 1 && 0 < jplus && jplus < Ny - 1 && ValidNumber(T[iplus][jplus])) {
                            nearest_neighbor_i[i][j] = iplus;
                            nearest_neighbor_j[i][j] = jplus;
                            break;
                        }
                    }
                }
            }
        }
    };
    /// Take all the matrices that are in the class and pack them into the matrices map for easy unpacking using msgpack
    void pack() {
/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
#define X(name) matrices.insert(std::pair<std::string, std::vector<std::vector<double>>>(#name, name));
        LIST_OF_MATRICES
#undef X
    };
    std::map<std::string, std::vector<std::vector<double>>>::iterator get_matrices_iterator(const std::string& name) {
        std::map<std::string, std::vector<std::vector<double>>>::iterator it = matrices.find(name);
        if (it == matrices.end()) {
            throw UnableToLoadError(format("could not find matrix %s", name.c_str()));
        }
        return it;
    }
    /// Take all the matrices that are in the class and pack them into the matrices map for easy unpacking using msgpack
    void unpack() {
/* Use X macros to auto-generate the unpacking code; each will look something like: T = *(matrices.find("T")).second */
#define X(name) name = get_matrices_iterator(#name)->second;
        LIST_OF_MATRICES
#undef X
        Nx = T.size();
        Ny = T[0].size();
        make_axis_vectors();
        make_good_neighbors();
    };
    /// Check that the native inputs (the inputs the table is based on) are in range
    bool native_inputs_are_in_range(double x, double y) {
        double e = 10 * DBL_EPSILON;
        return x >= xmin - e && x <= xmax + e && y >= ymin - e && y <= ymax + e;
    }
    /// @brief Find the nearest neighbor for native inputs (the inputs the table is based on)
    /// Does not check whether this corresponds to a valid node or not
    /// Use bisection since it is faster than calling a logarithm (surprising, but true)
    void find_native_nearest_neighbor(double x, double y, std::size_t& i, std::size_t& j) {
        bisect_vector(xvec, x, i);
        if (i != Nx - 1) {
            if (!logx) {
                if (x > (xvec[i] + xvec[i + 1]) / 2.0) {
                    i++;
                }
            } else {
                if (x > sqrt(xvec[i] * xvec[i + 1])) {
                    i++;
                }
            }
        }
        bisect_vector(yvec, y, j);
        if (j != Ny - 1) {
            if (!logy) {
                if (y > (yvec[j] + yvec[j + 1]) / 2.0) {
                    j++;
                }
            } else {
                if (y > sqrt(yvec[j] * yvec[j + 1])) {
                    j++;
                }
            }
        }
    }
    /// @brief Find the nearest neighbor for one (given) variable native, one variable non-native
    void find_nearest_neighbor(parameters givenkey, double givenval, parameters otherkey, double otherval, std::size_t& i, std::size_t& j) {
        if (givenkey == ykey) {
            bisect_vector(yvec, givenval, j);
            // This one is problematic because we need to make a slice against the grain in the "matrix"
            // which requires a slightly different algorithm
            try {
                bisect_segmented_vector_slice(get(otherkey), j, otherval, i);
            } catch (...) {
                // Now we go for a less intelligent solution, we simply try to find the one that is the closest
                const std::vector<std::vector<double>>& mat = get(otherkey);
                double closest_diff = 1e20;
                std::size_t closest_i = 0;
                for (std::size_t index = 0; index < mat.size(); ++index) {
                    double diff = std::abs(mat[index][j] - otherval);
                    if (diff < closest_diff) {
                        closest_diff = diff;
                        closest_i = index;
                    }
                }
                i = closest_i;
            }
        } else if (givenkey == xkey) {
            bisect_vector(xvec, givenval, i);
            // This one is fine because we now end up with a vector<double> in the other variable
            const std::vector<std::vector<double>>& v = get(otherkey);
            bisect_vector(v[i], otherval, j);
        }
    }
    /// Find the nearest good neighbor node for inputs that are the same as the grid inputs
    /// If the straightforward node (i,j) obtained by bisection is no good, find its nearest good node
    void find_native_nearest_good_neighbor(double x, double y, std::size_t& i, std::size_t& j) {
        // Get the node that is closest
        find_native_nearest_neighbor(x, y, i, j);
        // Check whether found node is good
        if (!ValidNumber(T[i][j])) {
            // If not, find its nearest good neighbor
            // (nearest good neighbors are precalculated and cached)
            std::size_t inew = nearest_neighbor_i[i][j];
            std::size_t jnew = nearest_neighbor_j[i][j];
            i = inew;
            j = jnew;
        }
    }
    /// Find the nearest cell with lower left coordinate (i,j) where (i,j) is a good node, and so are (i+1,j), (i,j+1), (i+1,j+1)
    /// This is needed for bicubic interpolation
    void find_native_nearest_good_cell(double x, double y, std::size_t& i, std::size_t& j) {
        bisect_vector(xvec, x, i);
        bisect_vector(yvec, y, j);
    }
    const std::vector<std::vector<double>>& get(parameters key) {
        switch (key) {
            case iDmolar:
                return rhomolar;
            case iT:
                return T;
            case iUmolar:
                return umolar;
            case iHmolar:
                return hmolar;
            case iSmolar:
                return smolar;
            case iP:
                return p;
            case iviscosity:
                return visc;
            case iconductivity:
                return cond;
            default:
                throw KeyError(format("invalid key"));
        }
    }
};

/// This class holds the single-phase data for a log(p)-h gridded table
class LogPHTable : public SinglePhaseGriddedTableData
{
   public:
    LogPHTable() {
        xkey = iHmolar;
        ykey = iP;
        logy = true;
        logx = false;
    };
    void set_limits() {
        if (this->AS.get() == NULL) {
            throw ValueError("AS is not yet set");
        }
        CoolPropDbl Tmin = std::max(AS->Ttriple(), AS->Tmin());
        // Minimum enthalpy is the saturated liquid enthalpy
        AS->update(QT_INPUTS, 0, Tmin);
        xmin = AS->hmolar();
        ymin = AS->p();

        // Check both the enthalpies at the Tmax isotherm to see whether to use low or high pressure
        AS->update(DmolarT_INPUTS, 1e-10, 1.499 * AS->Tmax());
        CoolPropDbl xmax1 = AS->hmolar();
        AS->update(PT_INPUTS, AS->pmax(), 1.499 * AS->Tmax());
        CoolPropDbl xmax2 = AS->hmolar();
        xmax = std::max(xmax1, xmax2);

        ymax = AS->pmax();
    }
    void deserialize(msgpack::object& deserialized) {
        LogPHTable temp;
        deserialized.convert(temp);
        temp.unpack();
        if (Nx != temp.Nx || Ny != temp.Ny) {
            throw ValueError(format("old [%dx%d] and new [%dx%d] dimensions don't agree", temp.Nx, temp.Ny, Nx, Ny));
        } else if (revision > temp.revision) {
            throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
        } else if ((std::abs(xmin) > 1e-10 && std::abs(xmax) > 1e-10)
                   && (std::abs(temp.xmin - xmin) / xmin > 1e-6 || std::abs(temp.xmax - xmax) / xmax > 1e-6)) {
            throw ValueError(format("Current limits for x [%g,%g] do not agree with loaded limits [%g,%g]", xmin, xmax, temp.xmin, temp.xmax));
        } else if ((std::abs(ymin) > 1e-10 && std::abs(ymax) > 1e-10)
                   && (std::abs(temp.ymin - ymin) / ymin > 1e-6 || std::abs(temp.ymax - ymax) / ymax > 1e-6)) {
            throw ValueError(format("Current limits for y [%g,%g] do not agree with loaded limits [%g,%g]", ymin, ymax, temp.ymin, temp.ymax));
        }
        std::swap(*this, temp);  // Swap
        this->AS = temp.AS;      // Reconnect the AbstractState pointer
    };
};
/// This class holds the single-phase data for a log(p)-T gridded table
class LogPTTable : public SinglePhaseGriddedTableData
{
   public:
    LogPTTable() {
        xkey = iT;
        ykey = iP;
        logy = true;
        logx = false;
        xmin = _HUGE;
        ymin = _HUGE;
        xmax = _HUGE;
        ymax = _HUGE;
    };
    void set_limits() {
        if (this->AS.get() == NULL) {
            throw ValueError("AS is not yet set");
        }
        CoolPropDbl Tmin = std::max(AS->Ttriple(), AS->Tmin());
        AS->update(QT_INPUTS, 0, Tmin);
        xmin = Tmin;
        ymin = AS->p();

        xmax = AS->Tmax() * 1.499;
        ymax = AS->pmax();
    }
    void deserialize(msgpack::object& deserialized) {
        LogPTTable temp;
        deserialized.convert(temp);
        temp.unpack();
        if (Nx != temp.Nx || Ny != temp.Ny) {
            throw ValueError(format("old [%dx%d] and new [%dx%d] dimensions don't agree", temp.Nx, temp.Ny, Nx, Ny));
        } else if (revision > temp.revision) {
            throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
        } else if ((std::abs(xmin) > 1e-10 && std::abs(xmax) > 1e-10)
                   && (std::abs(temp.xmin - xmin) / xmin > 1e-6 || std::abs(temp.xmax - xmax) / xmax > 1e-6)) {
            throw ValueError(format("Current limits for x [%g,%g] do not agree with loaded limits [%g,%g]", xmin, xmax, temp.xmin, temp.xmax));
        } else if ((std::abs(ymin) > 1e-10 && std::abs(ymax) > 1e-10)
                   && (std::abs(temp.ymin - ymin) / ymin > 1e-6 || std::abs(temp.ymax - ymax) / ymax > 1e-6)) {
            throw ValueError(format("Current limits for y [%g,%g] do not agree with loaded limits [%g,%g]", ymin, ymax, temp.ymin, temp.ymax));
        }
        std::swap(*this, temp);  // Swap
        this->AS = temp.AS;      // Reconnect the AbstractState pointer
    };
};

/// This structure holds the coefficients for one cell, the coefficients are stored in matrices
/// and can be obtained by the get() function.
class CellCoeffs
{
   private:
    std::size_t alt_i, alt_j;
    bool _valid, _has_valid_neighbor;

   public:
    double dx_dxhat, dy_dyhat;
    CellCoeffs() {
        _valid = false;
        _has_valid_neighbor = false;
        dx_dxhat = _HUGE;
        dy_dyhat = _HUGE;
        alt_i = 9999999;
        alt_j = 9999999;
    }
    std::vector<double> T, rhomolar, hmolar, p, smolar, umolar;
    /// Return a const reference to the desired matrix
    const std::vector<double>& get(const parameters params) const {
        switch (params) {
            case iT:
                return T;
            case iP:
                return p;
            case iDmolar:
                return rhomolar;
            case iHmolar:
                return hmolar;
            case iSmolar:
                return smolar;
            case iUmolar:
                return umolar;
            default:
                throw KeyError(format("Invalid key to get() function of CellCoeffs"));
        }
    };
    /// Set one of the matrices in this class
    void set(parameters params, const std::vector<double>& mat) {
        switch (params) {
            case iT:
                T = mat;
                break;
            case iP:
                p = mat;
                break;
            case iDmolar:
                rhomolar = mat;
                break;
            case iHmolar:
                hmolar = mat;
                break;
            case iSmolar:
                smolar = mat;
                break;
            case iUmolar:
                umolar = mat;
                break;
            default:
                throw KeyError(format("Invalid key to set() function of CellCoeffs"));
        }
    };
    /// Returns true if the cell coefficients seem to have been calculated properly
    bool valid() const {
        return _valid;
    };
    /// Call this function to set the valid flag to true
    void set_valid() {
        _valid = true;
    };
    /// Call this function to set the valid flag to false
    void set_invalid() {
        _valid = false;
    };
    /// Set the neighboring (alternate) cell to be used if the cell is invalid
    void set_alternate(std::size_t i, std::size_t j) {
        alt_i = i;
        alt_j = j;
        _has_valid_neighbor = true;
    }
    /// Get neighboring(alternate) cell to be used if this cell is invalid
    void get_alternate(std::size_t& i, std::size_t& j) const {
        if (_has_valid_neighbor) {
            i = alt_i;
            j = alt_j;
        } else {
            throw ValueError("No valid neighbor");
        }
    }
    /// Returns true if cell is invalid and it has valid neighbor
    bool has_valid_neighbor() const {
        return _has_valid_neighbor;
    }
};

/// This class contains the data for one set of Tabular data including single-phase and two-phase data
class TabularDataSet
{
   public:
    bool tables_loaded;
    LogPHTable single_phase_logph;
    LogPTTable single_phase_logpT;
    PureFluidSaturationTableData pure_saturation;
    PackablePhaseEnvelopeData phase_envelope;
    std::vector<std::vector<CellCoeffs>> coeffs_ph, coeffs_pT;

    TabularDataSet() {
        tables_loaded = false;
    }
    /// Write the tables to files on the computer
    void write_tables(const std::string& path_to_tables);
    /// Load the tables from file
    void load_tables(const std::string& path_to_tables, shared_ptr<CoolProp::AbstractState>& AS);
    /// Build the tables (single-phase PH, single-phase PT, phase envelope, etc.)
    void build_tables(shared_ptr<CoolProp::AbstractState>& AS);
    /// Build the \f$a_{i,j}\f$ coefficients for bicubic interpolation
    void build_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs);
};

class TabularDataLibrary
{
   private:
    std::map<std::string, TabularDataSet> data;

   public:
    TabularDataLibrary(){};
    std::string path_to_tables(shared_ptr<CoolProp::AbstractState>& AS) {
        std::vector<std::string> fluids = AS->fluid_names();
        std::vector<CoolPropDbl> fractions = AS->get_mole_fractions();
        std::vector<std::string> components;
        for (std::size_t i = 0; i < fluids.size(); ++i) {
            components.push_back(format("%s[%0.10Lf]", fluids[i].c_str(), fractions[i]));
        }
        std::string table_directory = get_home_dir() + "/.CoolProp/Tables/";
        std::string alt_table_directory = get_config_string(ALTERNATIVE_TABLES_DIRECTORY);
        if (!alt_table_directory.empty()) {
            table_directory = alt_table_directory;
        }
        return table_directory + AS->backend_name() + "(" + strjoin(components, "&") + ")";
    }
    /// Return a pointer to the set of tabular datasets
    TabularDataSet* get_set_of_tables(shared_ptr<AbstractState>& AS, bool& loaded);
};

/**
 * @brief This class contains the general code for tabular backends (TTSE, bicubic, etc.)
 *
 * This class layout was used in order to move the general code needed for all backends (building,  writing, loading)
 * into a common base class in order to remove code duplication.  DRY!
 */
class TabularBackend : public AbstractState
{
   protected:
    phases imposed_phase_index;
    bool tables_loaded, using_single_phase_table, is_mixture;
    enum selected_table_options
    {
        SELECTED_NO_TABLE = 0,
        SELECTED_PH_TABLE,
        SELECTED_PT_TABLE
    };
    selected_table_options selected_table;
    std::size_t cached_single_phase_i, cached_single_phase_j;
    std::size_t cached_saturation_iL, cached_saturation_iV;
    std::vector<std::vector<double>> const* z;
    std::vector<std::vector<double>> const* dzdx;
    std::vector<std::vector<double>> const* dzdy;
    std::vector<std::vector<double>> const* d2zdx2;
    std::vector<std::vector<double>> const* d2zdxdy;
    std::vector<std::vector<double>> const* d2zdy2;
    std::vector<CoolPropDbl> mole_fractions;

   public:
    shared_ptr<CoolProp::AbstractState> AS;
    TabularBackend(shared_ptr<CoolProp::AbstractState> AS) : tables_loaded(false), using_single_phase_table(false), is_mixture(false), AS(AS) {
        selected_table = SELECTED_NO_TABLE;
        // Flush the cached indices (set to large number)
        cached_single_phase_i = std::numeric_limits<std::size_t>::max();
        cached_single_phase_j = std::numeric_limits<std::size_t>::max();
        cached_saturation_iL = std::numeric_limits<std::size_t>::max();
        cached_saturation_iV = std::numeric_limits<std::size_t>::max();
        z = NULL;
        dzdx = NULL;
        dzdy = NULL;
        d2zdx2 = NULL;
        d2zdxdy = NULL;
        d2zdy2 = NULL;
        dataset = NULL;
        imposed_phase_index = iphase_not_imposed;
    };

    // None of the tabular methods are available from the high-level interface
    bool available_in_high_level(void) {
        return false;
    }

    std::string calc_name(void) {
        return AS->name();
    }
    std::vector<std::string> calc_fluid_names(void) {
        return AS->fluid_names();
    }

    void connect_pointers(parameters output, const SinglePhaseGriddedTableData& table) {
        // Connect the pointers based on the output variable desired
        switch (output) {
            case iT:
                z = &table.T;
                dzdx = &table.dTdx;
                dzdy = &table.dTdy;
                d2zdxdy = &table.d2Tdxdy;
                d2zdx2 = &table.d2Tdx2;
                d2zdy2 = &table.d2Tdy2;
                break;
            case iDmolar:
                z = &table.rhomolar;
                dzdx = &table.drhomolardx;
                dzdy = &table.drhomolardy;
                d2zdxdy = &table.d2rhomolardxdy;
                d2zdx2 = &table.d2rhomolardx2;
                d2zdy2 = &table.d2rhomolardy2;
                break;
            case iSmolar:
                z = &table.smolar;
                dzdx = &table.dsmolardx;
                dzdy = &table.dsmolardy;
                d2zdxdy = &table.d2smolardxdy;
                d2zdx2 = &table.d2smolardx2;
                d2zdy2 = &table.d2smolardy2;
                break;
            case iHmolar:
                z = &table.hmolar;
                dzdx = &table.dhmolardx;
                dzdy = &table.dhmolardy;
                d2zdxdy = &table.d2hmolardxdy;
                d2zdx2 = &table.d2hmolardx2;
                d2zdy2 = &table.d2hmolardy2;
                break;
            case iUmolar:
                z = &table.umolar;
                dzdx = &table.dumolardx;
                dzdy = &table.dumolardy;
                d2zdxdy = &table.d2umolardxdy;
                d2zdx2 = &table.d2umolardx2;
                d2zdy2 = &table.d2umolardy2;
                break;
            case iviscosity:
                z = &table.visc;
                break;
            case iconductivity:
                z = &table.cond;
                break;
            default:
                throw ValueError();
        }
    }
    TabularDataSet* dataset;

    void recalculate_singlephase_phase() {
        if (p() > p_critical()) {
            if (T() > T_critical()) {
                _phase = iphase_supercritical;
            } else {
                _phase = iphase_supercritical_liquid;
            }
        } else {
            if (T() > T_critical()) {
                _phase = iphase_supercritical_gas;
            } else {
                // Liquid or vapor
                if (rhomolar() > rhomolar_critical()) {
                    _phase = iphase_liquid;
                } else {
                    _phase = iphase_gas;
                }
            }
        }
    }
    /** \brief Specify the phase - this phase will always be used in calculations
        *
        * @param phase_index The index from CoolProp::phases
        */
    void calc_specify_phase(phases phase_index) {
        imposed_phase_index = phase_index;
    };

    /**\brief Unspecify the phase - the phase is no longer imposed, different solvers can do as they like
        */
    void calc_unspecify_phase() {
        imposed_phase_index = iphase_not_imposed;
    };

    virtual double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) = 0;
    virtual double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) = 0;
    virtual double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j) = 0;
    virtual double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j) = 0;
    virtual double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) = 0;
    virtual double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) = 0;

    /// Ask the derived class to find the nearest good set of i,j that it wants to use (pure virtual)
    virtual void find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, double x,
                                                  double y, std::size_t& i, std::size_t& j) = 0;
    /// Ask the derived class to find the nearest neighbor (pure virtual)
    virtual void find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                       const parameters variable1, const double value1, const parameters other, const double otherval, std::size_t& i,
                                       std::size_t& j) = 0;
    ///
    virtual void invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                       parameters output, double x, double y, std::size_t i, std::size_t j) = 0;
    ///
    virtual void invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                       parameters output, double x, double y, std::size_t i, std::size_t j) = 0;

    phases calc_phase(void) {
        return _phase;
    }
    CoolPropDbl calc_T_critical(void) {
        return this->AS->T_critical();
    };
    CoolPropDbl calc_Ttriple(void) {
        return this->AS->Ttriple();
    };
    CoolPropDbl calc_p_triple(void) {
        return this->AS->p_triple();
    };
    CoolPropDbl calc_pmax(void) {
        return this->AS->pmax();
    };
    CoolPropDbl calc_Tmax(void) {
        return this->AS->Tmax();
    };
    CoolPropDbl calc_Tmin(void) {
        return this->AS->Tmin();
    };
    CoolPropDbl calc_p_critical(void) {
        return this->AS->p_critical();
    }
    CoolPropDbl calc_rhomolar_critical(void) {
        return this->AS->rhomolar_critical();
    }
    bool using_mole_fractions(void) {
        return true;
    }
    bool using_mass_fractions(void) {
        return false;
    }
    bool using_volu_fractions(void) {
        return false;
    }
    void update(CoolProp::input_pairs input_pair, double Value1, double Value2);
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        this->AS->set_mole_fractions(mole_fractions);
    };
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
        throw NotImplementedError("set_mass_fractions not implemented for Tabular backends");
    };
    const std::vector<CoolPropDbl>& get_mole_fractions() {
        return AS->get_mole_fractions();
    };
    const std::vector<CoolPropDbl> calc_mass_fractions(void) {
        return AS->get_mass_fractions();
    };

    CoolPropDbl calc_molar_mass(void) {
        return AS->molar_mass();
    };

    CoolPropDbl calc_saturated_liquid_keyed_output(parameters key);
    CoolPropDbl calc_saturated_vapor_keyed_output(parameters key);

    /// Returns the path to the tables that shall be written
    std::string path_to_tables(void);
    /// Load the tables from file; throws UnableToLoadException if there is a problem
    void load_tables();
    void pack_matrices() {
        PackablePhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
        PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
        SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
        SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
        single_phase_logph.pack();
        single_phase_logpT.pack();
        pure_saturation.pack();
        phase_envelope.pack();
    }
    /// Write the tables to file
    void write_tables();

    CoolPropDbl phase_envelope_sat(const PhaseEnvelopeData& env, parameters output, parameters iInput1, double value1) {
        const PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
        CoolPropDbl yL = PhaseEnvelopeRoutines::evaluate(phase_envelope, output, iInput1, value1, cached_saturation_iL);
        CoolPropDbl yV = PhaseEnvelopeRoutines::evaluate(phase_envelope, output, iInput1, value1, cached_saturation_iV);
        return _Q * yV + (1 - _Q) * yL;
    }
    CoolPropDbl calc_cpmolar_idealgas(void) {
        this->AS->set_T(_T);
        return this->AS->cp0molar();
    }
    /// Calculate the surface tension using the wrapped class (fast enough)
    CoolPropDbl calc_surface_tension(void) {
        this->AS->set_T(_T);
        return this->AS->surface_tension();
        this->AS->set_T(_HUGE);
    }
    CoolPropDbl calc_p(void);
    CoolPropDbl calc_T(void);
    CoolPropDbl calc_rhomolar(void);
    CoolPropDbl calc_hmolar(void);
    CoolPropDbl calc_smolar(void);
    CoolPropDbl calc_umolar(void);
    CoolPropDbl calc_cpmolar(void);
    CoolPropDbl calc_cvmolar(void);
    CoolPropDbl calc_viscosity(void);
    CoolPropDbl calc_conductivity(void);
    /// Calculate the speed of sound using a tabular backend [m/s]
    CoolPropDbl calc_speed_sound(void);
    CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);
    /** /brief calculate the derivative along the saturation curve, but only if quality is 0 or 1
        */
    CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1);
    CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant);

    /// If you need all three values (drho_dh__p, drho_dp__h and rho_spline), you should calculate drho_dp__h first to avoid duplicate calculations.
    CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end);

    void check_tables() {
        if (!tables_loaded) {
            try {
                /// Try to load the tables if you can.
                load_tables();
                // Set the flag saying tables have been successfully loaded
                tables_loaded = true;
            } catch (CoolProp::UnableToLoadError& e) {
                if (get_debug_level() > 0) {
                    std::cout << format("Table loading failed with error: %s\n", e.what());
                }
                /// Check directory size
                std::string table_path = path_to_tables();
#if defined(__ISWINDOWS__)
                double directory_size_in_GB = CalculateDirSize(std::wstring(table_path.begin(), table_path.end())) / POW3(1024.0);
#else
                double directory_size_in_GB = CalculateDirSize(table_path) / POW3(1024.0);
#endif
                double allowed_size_in_GB = get_config_double(MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB);
                if (get_debug_level() > 0) {
                    std::cout << "Tabular directory size is " << directory_size_in_GB << " GB\n";
                }
                if (directory_size_in_GB > 1.5 * allowed_size_in_GB) {
                    throw DirectorySizeError(
                      format("Maximum allowed tabular directory size is %g GB, you have exceeded 1.5 times this limit", allowed_size_in_GB));
                } else if (directory_size_in_GB > allowed_size_in_GB) {
                    set_warning_string(format("Maximum allowed tabular directory size is %g GB, you have exceeded this limit", allowed_size_in_GB));
                }
                /// If you cannot load the tables, build them and then write them to file
                dataset->build_tables(this->AS);
                pack_matrices();
                write_tables();
                /// Load the tables back into memory as a consistency check
                load_tables();
                // Set the flag saying tables have been successfully loaded
                tables_loaded = true;
            }
        }
    };
};

} /* namespace CoolProp*/

#endif
