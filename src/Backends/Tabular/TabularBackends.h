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


/** ***MAGIC WARNING***!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_MATRICES X(T) X(p) X(rhomolar) X(hmolar) X(smolar) X(umolar) X(dTdx) X(dTdy) X(dpdx) X(dpdy) X(drhomolardx) X(drhomolardy) X(dhmolardx) X(dhmolardy) X(dsmolardx) X(dsmolardy) X(dumolardx) X(dumolardy) X(d2Tdx2) X(d2Tdxdy) X(d2Tdy2) X(d2pdx2) X(d2pdxdy) X(d2pdy2) X(d2rhomolardx2) X(d2rhomolardxdy) X(d2rhomolardy2) X(d2hmolardx2) X(d2hmolardxdy) X(d2hmolardy2) X(d2smolardx2) X(d2smolardxdy) X(d2smolardy2) X(d2umolardx2) X(d2umolardxdy) X(d2umolardy2) X(visc) X(cond)

/** ***MAGIC WARNING***!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_SATURATION_VECTORS X(TL) X(pL) X(logpL) X(hmolarL) X(smolarL) X(umolarL) X(rhomolarL) X(logrhomolarL) X(viscL) X(condL) X(logviscL) X(TV) X(pV) X(logpV) X(hmolarV) X(smolarV) X(umolarV) X(rhomolarV) X(logrhomolarV) X(viscV) X(condV) X(logviscV)

namespace CoolProp{

/// Get a conversion factor from mass to molar if needed
inline void mass_to_molar(parameters &param, double &conversion_factor, double molar_mass){
    switch(param){
        case iDmass: conversion_factor = molar_mass; param = iDmolar; break;
        case iHmass: conversion_factor = 1/molar_mass; param = iHmolar; break;
        case iSmass: conversion_factor = 1/molar_mass; param = iSmolar; break;
        case iUmass: conversion_factor = 1/molar_mass; param = iUmolar; break;
        case iDmolar:
        case iHmolar:
        case iSmolar:
        case iUmolar:
        case iT:
        case iP:
            return;
        default:
            throw ValueError("I don't know how to convert this parameter");
    }
}

/** \brief This class holds the data for a two-phase table that is log spaced in p
 * 
 * It contains very few members or methods, mostly it just holds the data
 */
class PureFluidSaturationTableData{
	public:
		std::size_t N;
		shared_ptr<CoolProp::AbstractState> AS;
    
		PureFluidSaturationTableData(){N = 1000; revision = 0;}
        
        /// Build this table
        void build(shared_ptr<CoolProp::AbstractState> &AS);
    
		/* Use X macros to auto-generate the variables; each will look something like: std::vector<double> T; */
		#define X(name) std::vector<double> name;
		LIST_OF_SATURATION_VECTORS
		#undef X

		int revision;
		std::map<std::string, std::vector<double> > vectors;
    
		MSGPACK_DEFINE(revision, vectors); // write the member variables that you want to pack

        bool is_inside(parameters main, double mainval, parameters other, double val, std::size_t &iL, std::size_t &iV, CoolPropDbl &yL, CoolPropDbl &yV){
            std::vector<double> *yvecL = NULL, *yvecV = NULL;
            switch(other){
                case iT: yvecL = &TL; yvecV = &TV; break;
                case iHmolar: yvecL = &hmolarL; yvecV = &hmolarV; break;
                case iQ: yvecL = &TL; yvecV = &TV; break;
                case iSmolar: yvecL = &smolarL; yvecV = &smolarV; break;
                case iUmolar: yvecL = &umolarL; yvecV = &umolarV; break;
                case iDmolar: yvecL = &rhomolarL; yvecV = &rhomolarV; break;
                default: throw ValueError("invalid input for other in is_inside");
            }
            
            // Trivial checks
            if (main == iP){
                // If p is outside the range (ptriple, pcrit), considered to not be inside
                double pmax = this->pV[pV.size()-1], pmin = this->pV[0];
                if (mainval > pmax || mainval < pmin){return false;}
            }
            else if (main == iT){
                // If p is outside the range (ptriple, pcrit), considered to not be inside
                double Tmax = this->TV[TV.size()-1], Tmin = this->TV[0];
                if (mainval > Tmax || mainval < Tmin){return false;}
            }
            else{
                throw ValueError("invalid input for other in is_inside");
            }
            
            // Now check based on a rough analysis using bounding pressure
            std::size_t iLplus, iVplus;
            // Find the indices (iL,iL+1) & (iV,iV+1) that bound the given pressure
            // In general iV and iL will be the same, but if pseudo-pure, they might
            // be different
            if (main ==iP){
                bisect_vector(pV, mainval, iV);
                bisect_vector(pL, mainval, iL);
            }
            else if (main == iT){
                bisect_vector(TV, mainval, iV);
                bisect_vector(TL, mainval, iL);
            }
			if (other == iQ){return true;}
            iVplus = std::min(iV+1, N-1);
            iLplus = std::min(iL+1, N-1);
            // Find the bounding values for the other variable
            double ymin = min4((*yvecL)[iL],(*yvecL)[iLplus],(*yvecV)[iV],(*yvecV)[iVplus]);
            double ymax = max4((*yvecL)[iL],(*yvecL)[iLplus],(*yvecV)[iV],(*yvecV)[iVplus]);
            if (val < ymin || val > ymax){ return false;}
            // Actually do "saturation" call using cubic interpolation
            if (iVplus < 3){ iVplus = 3;}
            if (iLplus < 3){ iLplus = 3;}
            if (main==iP){
                double logp = log(mainval);
                yV = CubicInterp(logpV, *yvecV, iVplus-3, iVplus-2, iVplus-1, iVplus, logp);
                yL = CubicInterp(logpL, *yvecL, iLplus-3, iLplus-2, iLplus-1, iLplus, logp);
            }
            else if (main == iT){
                yV = CubicInterp(TV, *yvecV, iVplus-3, iVplus-2, iVplus-1, iVplus, mainval);
                yL = CubicInterp(TL, *yvecL, iLplus-3, iLplus-2, iLplus-1, iLplus, mainval);
            }

            if (!is_in_closed_range(yV, yL, static_cast<CoolPropDbl>(val))){ 
                return false;
            }
            else{
                iL = iLplus-1;
                iV = iVplus-1;
                return true;
            }
        }
		/// Resize all the vectors
		void resize(std::size_t N){
			/* Use X macros to auto-generate the code; each will look something like: T.resize(N); std::fill(T.begin(), T.end(), _HUGE); */
			#define X(name) name.resize(N); std::fill(name.begin(), name.end(), _HUGE);
			LIST_OF_SATURATION_VECTORS
			#undef X
		};
        /// Take all the vectors that are in the class and pack them into the vectors map for easy unpacking using msgpack
		void pack(){
			/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
			#define X(name) vectors.insert(std::pair<std::string, std::vector<double> >(#name, name));
			LIST_OF_SATURATION_VECTORS
			#undef X
		};
        std::map<std::string, std::vector<double> >::iterator get_vector_iterator(const std::string &name){
            std::map<std::string, std::vector<double> >::iterator it = vectors.find(name);
            if (it == vectors.end()){
                throw UnableToLoadError(format("could not find vector %s",name.c_str()));
            }
            return it;
        }
		/// Take all the vectors that are in the class and unpack them from the vectors map
		void unpack(){
			/* Use X macros to auto-generate the unpacking code; each will look something like: T = get_vector_iterator("T")->second */
			#define X(name) name = get_vector_iterator(#name)->second;
			LIST_OF_SATURATION_VECTORS
			#undef X
			N = TL.size();
		};
        void deserialize(msgpack::object &deserialized){       
            PureFluidSaturationTableData temp;
            deserialized.convert(&temp);
            temp.unpack();
            if (N != temp.N)
            {
                throw ValueError(format("old [%d] and new [%d] sizes don't agree", temp.N, N));
            }
            else if (revision > temp.revision)
            {
                throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
            }
            std::swap(*this, temp); // Swap
            this->AS = temp.AS; // Reconnect the AbstractState pointer
        };
        double evaluate(parameters output, double p, double Q, std::size_t iL, std::size_t iV)
        {
            double logp = log(p);
            switch(output){
                case iT:
                {
                    double TV = CubicInterp(logpV, this->TV, iV-2, iV-1, iV, iV+1, logp);
                    double TL = CubicInterp(logpL, this->TL, iL-2, iL-1, iL, iL+1, logp);
                    return Q*TV + (1-Q)*TL;
                }
                case iSmolar:
                {
                    double sV = CubicInterp(logpV, smolarV, iV-2, iV-1, iV, iV+1, logp);
                    double sL = CubicInterp(logpL, smolarL, iL-2, iL-1, iL, iL+1, logp);
                    return Q*sV + (1-Q)*sL;
                }
                case iHmolar:
                {
                    double hV = CubicInterp(logpV, hmolarV, iV-2, iV-1, iV, iV+1, logp);
                    double hL = CubicInterp(logpL, hmolarL, iL-2, iL-1, iL, iL+1, logp);
                    return Q*hV + (1-Q)*hL;
                }
				case iUmolar:
                {
                    double uV = CubicInterp(logpV, umolarV, iV-2, iV-1, iV, iV+1, logp);
                    double uL = CubicInterp(logpL, umolarL, iL-2, iL-1, iL, iL+1, logp);
                    return Q*uV + (1-Q)*uL;
                }
                case iDmolar:
                {
                    double rhoV = exp(CubicInterp(logpV, logrhomolarV, iV-2, iV-1, iV, iV+1, logp));
                    double rhoL = exp(CubicInterp(logpL, logrhomolarL, iL-2, iL-1, iL, iL+1, logp));
                    if (!ValidNumber(rhoV)){throw ValueError("rhoV is invalid");}
                    if (!ValidNumber(rhoL)){throw ValueError("rhoL is invalid");}
                    return 1/(Q/rhoV + (1-Q)/rhoL);
                }
				case iconductivity:
                {
                    double kV = CubicInterp(logpV, condV, iV-2, iV-1, iV, iV+1, logp);
                    double kL = CubicInterp(logpL, condL, iL-2, iL-1, iL, iL+1, logp);
                    if (!ValidNumber(kV)){throw ValueError("kV is invalid");}
                    if (!ValidNumber(kL)){throw ValueError("kL is invalid");}
                    return Q*kV + (1-Q)*kL;
                }
				case iviscosity:
                {
                    double muV = exp(CubicInterp(logpV, logviscV, iV-2, iV-1, iV, iV+1, logp));
                    double muL = exp(CubicInterp(logpL, logviscL, iL-2, iL-1, iL, iL+1, logp));
                    if (!ValidNumber(muV)){throw ValueError("muV is invalid");}
                    if (!ValidNumber(muL)){throw ValueError("muL is invalid");}
                    return 1/(Q/muV + (1-Q)/muL);
                }

                default:
                    throw ValueError("can't be something other than T or rho");
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
        double first_saturation_deriv(parameters Of1, parameters Wrt1, int Q, double val, std::size_t i)
        {
            if (i < 2 || i > TL.size() - 2){throw ValueError("Invalid index to calc_first_saturation_deriv in TabularBackends");}
            std::vector<double> *x, *y;
            // Connect pointers for each vector
            switch(Wrt1){
                case iT: x = (Q == 0) ? &TL : &TV; break;
                case iP: x = (Q == 0) ? &pL : &pV; break;
                default: throw ValueError(format("Key for Wrt1 is invalid in calc_first_saturation_deriv"));
            }
            switch(Of1){
                case iT: y = (Q == 0) ? &TL : &TV; break;
                case iP: y = (Q == 0) ? &pL : &pV; break;
                case iDmolar: y = (Q == 0) ? &rhomolarL : &rhomolarV; break;
                case iHmolar: y = (Q == 0) ? &hmolarL : &hmolarV; break;
                case iSmolar: y = (Q == 0) ? &smolarL : &smolarV; break;
                case iUmolar: y = (Q == 0) ? &umolarL : &umolarV; break;
                default: throw ValueError(format("Key for Of1 is invalid in calc_first_saturation_deriv"));
            }
            return CubicInterpFirstDeriv((*x)[i-2], (*x)[i-1], (*x)[i], (*x)[i+1],
                                         (*y)[i-2], (*y)[i-1], (*y)[i], (*y)[i+1],
                                         val);
        };
        //calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant);
};

/** \brief This class holds the data for a single-phase interpolation table that is regularly spaced
 * 
 * It contains very few members or methods, mostly it just holds the data
 */
class SinglePhaseGriddedTableData{
        
	public:
		std::size_t Nx, Ny;
		CoolProp::parameters xkey, ykey;
		shared_ptr<CoolProp::AbstractState> AS;
		std::vector<double> xvec, yvec;
        std::vector<std::vector<std::size_t> > nearest_neighbor_i, nearest_neighbor_j;
		bool logx, logy;
		double xmin, ymin, xmax, ymax;
        
        virtual void set_limits() = 0;
    
		SinglePhaseGriddedTableData(){Nx = 200; Ny = 200; revision = 0;}
    
		/* Use X macros to auto-generate the variables; each will look something like: std::vector< std::vector<double> > T; */
		#define X(name) std::vector< std::vector<double> > name;
		LIST_OF_MATRICES
		#undef X
		int revision;
		std::map<std::string, std::vector<std::vector<double> > > matrices;
        /// Build this table
        void build(shared_ptr<CoolProp::AbstractState> &AS);
    
		MSGPACK_DEFINE(revision, matrices, xmin, xmax, ymin, ymax); // write the member variables that you want to pack
		/// Resize all the matrices
		void resize(std::size_t Nx, std::size_t Ny){
			/* Use X macros to auto-generate the code; each will look something like: T.resize(Nx, std::vector<double>(Ny, _HUGE)); */
			#define X(name) name.resize(Nx, std::vector<double>(Ny, _HUGE));
			LIST_OF_MATRICES
			#undef X
			make_axis_vectors();
		};
		/// Make vectors for the x-axis values and the y-axis values
		void make_axis_vectors(void){
			if (logx){
				xvec = logspace(xmin, xmax, Nx);
			}
			else{
				xvec = linspace(xmin, xmax, Nx);
			}
			if (logy){
				yvec = logspace(ymin, ymax, Ny);
			}
			else{
				yvec = linspace(ymin, ymax, Ny);
			}
		};
        /// Make matrices of good neighbors if the current value for i,j corresponds to a bad node
		void make_good_neighbors(void){
            nearest_neighbor_i.resize(Nx, std::vector<std::size_t>(Ny, std::numeric_limits<std::size_t>::max()));
            nearest_neighbor_j.resize(Nx, std::vector<std::size_t>(Ny, std::numeric_limits<std::size_t>::max()));
            for (std::size_t i = 0; i < xvec.size(); ++i){
                for (std::size_t j = 0; j < yvec.size(); ++j){
                    nearest_neighbor_i[i][j] = i;
                    nearest_neighbor_j[i][j] = j;
                    if (!ValidNumber(T[i][j])){
                        int xoffsets[] = {-1,1,0,0,-1,1,1,-1};
                        int yoffsets[] = {0,0,1,-1,-1,-1,1,1};
                        // Length of offset
                        std::size_t N = sizeof(xoffsets)/sizeof(xoffsets[0]);
                        for (std::size_t k = 0; k < N; ++k){
                            std::size_t iplus = i + xoffsets[k];
                            std::size_t jplus = j + yoffsets[k];
                            if (0 < iplus && iplus < Nx-1 && 0 < jplus && jplus < Ny-1 && ValidNumber(T[iplus][jplus])){
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
		void pack(){
			/* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
			#define X(name) matrices.insert(std::pair<std::string, std::vector<std::vector<double> > >(#name, name));
			LIST_OF_MATRICES
			#undef X
		};
        std::map<std::string, std::vector<std::vector<double> > >::iterator get_matrices_iterator(const std::string &name){
            std::map<std::string, std::vector<std::vector<double> > >::iterator it = matrices.find(name);
            if (it == matrices.end()){
                throw UnableToLoadError(format("could not find matrix %s",name.c_str()));
            }
            return it;
        }
		/// Take all the matrices that are in the class and pack them into the matrices map for easy unpacking using msgpack
		void unpack(){
			/* Use X macros to auto-generate the unpacking code; each will look something like: T = *(matrices.find("T")).second */
			#define X(name) name = get_matrices_iterator(#name)->second;
			LIST_OF_MATRICES
			#undef X
			Nx = T.size(); Ny = T[0].size();
			make_axis_vectors();
            make_good_neighbors();
		};
		/// Check that the native inputs (the inputs the table is based on) are in range
		bool native_inputs_are_in_range(double x, double y){
			return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
		}
		/// @brief Find the nearest neighbor for native inputs (the inputs the table is based on)
		/// Does not check whether this corresponds to a valid node or not
		/// Use bisection since it is faster than calling a logarithm (surprising, but true)
		void find_native_nearest_neighbor(double x, double y, std::size_t &i, std::size_t &j){
			bisect_vector(xvec, x, i);
			if (i != Nx-1){
				if(!logx){
					if (x > (xvec[i]+xvec[i+1])/2.0){i++;}
				}
				else{
					if (x > sqrt(xvec[i]*xvec[i+1])){i++;}
				}
			}
			bisect_vector(yvec, y, j);
			if (j != Ny-1){
				if(!logy){
					if (y > (yvec[j]+yvec[j+1])/2.0){j++;}
				}
				else{
					if (y > sqrt(yvec[j]*yvec[j+1])){j++;}
				}
			}
		}
        /// @brief Find the nearest neighbor for one (given) variable native, one variable non-native
		void find_nearest_neighbor(parameters givenkey, double givenval, parameters otherkey, double otherval, std::size_t &i, std::size_t &j){
			if (givenkey == ykey){
                bisect_vector(yvec, givenval, j);
                // This one is problematic because we need to make a slice against the grain in the "matrix"
                // which requires a slightly different algorithm
                bisect_segmented_vector_slice(get(otherkey), j, otherval, i);
            }
            else if (givenkey == xkey){
                bisect_vector(xvec, givenval, i);
                // This one is fine because we now end up with a vector<double> in the other variable
                const std::vector<std::vector<double> > & v = get(otherkey);
                bisect_vector(v[i], otherval, j);
                if (j < v[i].size()-1 && std::abs(v[i][j+1] - otherval) < std::abs(v[i][j] - otherval)){
                    j++;
                }
            }
		}
		/// Find the nearest good neighbor node for inputs that are the same as the grid inputs
		/// If the straightforward node (i,j) obtained by bisection is no good, find its nearest good node
		void find_native_nearest_good_neighbor(double x, double y, std::size_t &i, std::size_t &j){
			// Get the node that is closest
			find_native_nearest_neighbor(x,y,i,j);
			// Check whether found node is good
			if (!ValidNumber(T[i][j])){
				// If not, find its nearest good neighbor 
                // (nearest good neighbors are precalculated and cached)
                std::size_t inew = nearest_neighbor_i[i][j];
                std::size_t jnew = nearest_neighbor_j[i][j];
                i = inew; j = jnew;
			}
		}
		/// Find the nearest cell with lower left coordinate (i,j) where (i,j) is a good node, and so are (i+1,j), (i,j+1), (i+1,j+1)
		/// This is needed for bicubic interpolation
		void find_native_nearest_good_cell(double x, double y, std::size_t &i, std::size_t &j){
			bisect_vector(xvec, x, i);
			bisect_vector(yvec, y, j);
		}
        const std::vector<std::vector<double> > get(parameters key){
            switch(key){
                case iDmolar: return rhomolar;
                case iT: return T;
                case iUmolar: return umolar;
                case iHmolar: return hmolar;
                case iSmolar: return smolar;
                case iP: return p;
                case iviscosity: return visc;
                case iconductivity: return cond;
                default: throw KeyError(format("invalid key"));
            }
        }
};

/// This class holds the single-phase data for a log(p)-h gridded table
class LogPHTable : public SinglePhaseGriddedTableData
{
    public:
        LogPHTable(){
            xkey = iHmolar; ykey = iP; logy = true; logx = false;
        };
        void set_limits(){
            if (this->AS.get() == NULL){
                throw ValueError("AS is not yet set");
            }
            // Minimum enthalpy is the saturated liquid enthalpy
            AS->update(QT_INPUTS, 0, AS->Ttriple());
            xmin = AS->hmolar(); ymin = AS->p();
            
            // Check both the enthalpies at the Tmax isotherm to see whether to use low or high pressure
            AS->update(PT_INPUTS, 1e-10, AS->Tmax());
            CoolPropDbl xmax1 = AS->hmolar();
            AS->update(PT_INPUTS, AS->pmax(), AS->Tmax());
            CoolPropDbl xmax2 = AS->hmolar();
            xmax = std::max(xmax1, xmax2);
            
            ymax = AS->pmax();
        }
        void deserialize(msgpack::object &deserialized){       
            LogPHTable temp;
            deserialized.convert(&temp);
            temp.unpack();
            if (Nx != temp.Nx || Ny != temp.Ny)
            {
                throw ValueError(format("old [%dx%d] and new [%dx%d] dimensions don't agree", temp.Nx, temp.Ny, Nx, Ny));
            }
            else if (revision > temp.revision)
            {
                throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
            }
            else if ((std::abs(xmin) > 1e-10 && std::abs(xmax) > 1e-10) && (std::abs(temp.xmin - xmin)/xmin > 1e-6 || std::abs(temp.xmax - xmax)/xmax > 1e-6)){
                throw ValueError(format("Current limits for x [%g,%g] do not agree with loaded limits [%g,%g]", xmin, xmax, temp.xmin, temp.xmax));
            }
            else if ((std::abs(ymin) > 1e-10 && std::abs(ymax) > 1e-10) && (std::abs(temp.ymin - ymin)/ymin > 1e-6 || std::abs(temp.ymax - ymax)/ymax > 1e-6)){
                throw ValueError(format("Current limits for y [%g,%g] do not agree with loaded limits [%g,%g]", ymin, ymax, temp.ymin, temp.ymax));
            }
            std::swap(*this, temp); // Swap
            this->AS = temp.AS; // Reconnect the AbstractState pointer
        };
};
/// This class holds the single-phase data for a log(p)-T gridded table
class LogPTTable : public SinglePhaseGriddedTableData
{
    public:
        LogPTTable(){
            xkey = iT; ykey = iP; logy = true; logx = false;
        };
        void set_limits(){
            if (this->AS.get() == NULL){
                throw ValueError("AS is not yet set");
            }            
            xmin = AS->Ttriple();
            AS->update(QT_INPUTS, 0, AS->Ttriple());
            ymin = AS->p();
            
            xmax = AS->Tmax(); ymax = AS->pmax();
        }
        void deserialize(msgpack::object &deserialized){   
            LogPTTable temp;
            deserialized.convert(&temp);
            temp.unpack();
            if (Nx != temp.Nx || Ny != temp.Ny)
            {
                throw ValueError(format("old [%dx%d] and new [%dx%d] dimensions don't agree",temp.Nx, temp.Ny, Nx, Ny));
            }
            else if (revision > temp.revision)
            {
                throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp.revision, revision));
            }
            else if ((std::abs(xmin) > 1e-10 && std::abs(xmax) > 1e-10) && (std::abs(temp.xmin - xmin)/xmin > 1e-6 || std::abs(temp.xmax - xmax)/xmax > 1e-6)){
                throw ValueError(format("Current limits for x [%g,%g] do not agree with loaded limits [%g,%g]", xmin, xmax, temp.xmin, temp.xmax));
            }
            else if ((std::abs(ymin) > 1e-10 && std::abs(ymax) > 1e-10) && (std::abs(temp.ymin - ymin)/ymin > 1e-6 || std::abs(temp.ymax - ymax)/ymax > 1e-6)){
                throw ValueError(format("Current limits for y [%g,%g] do not agree with loaded limits [%g,%g]", ymin, ymax, temp.ymin, temp.ymax));
            }
            std::swap(*this, temp); // Swap
            this->AS = temp.AS; // Reconnect the AbstractState pointer
        };
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
        bool tables_loaded, using_single_phase_table;
        shared_ptr<CoolProp::AbstractState> AS;
        enum selected_table_options{SELECTED_NO_TABLE=0, SELECTED_PH_TABLE, SELECTED_PT_TABLE};
        selected_table_options selected_table;
        std::size_t cached_single_phase_i, cached_single_phase_j;
        std::size_t cached_saturation_iL, cached_saturation_iV;
        std::vector<std::vector<double> > *z, *dzdx, *dzdy, *d2zdx2, *d2zdxdy, *d2zdy2;
        std::vector<CoolPropDbl> mole_fractions;
    public:

        TabularBackend(shared_ptr<CoolProp::AbstractState> AS) : tables_loaded(false), using_single_phase_table(false), AS(AS) {};

        // None of the tabular methods are available from the high-level interface
        bool available_in_high_level(void){return false;}

        void connect_pointers(parameters output, SinglePhaseGriddedTableData &table)
		{
			// Connect the pointers based on the output variable desired
			switch(output){
				case iT:
					z = &table.T; dzdx = &table.dTdx; dzdy = &table.dTdy;
					d2zdxdy = &table.d2Tdxdy; d2zdx2 = &table.d2Tdx2; d2zdy2 = &table.d2Tdy2;
					break;
				case iDmolar:
					z = &table.rhomolar; dzdx = &table.drhomolardx; dzdy = &table.drhomolardy;
					d2zdxdy = &table.d2rhomolardxdy; d2zdx2 = &table.d2rhomolardx2; d2zdy2 = &table.d2rhomolardy2;
					break;
				case iSmolar:
                    z = &table.smolar; dzdx = &table.dsmolardx; dzdy = &table.dsmolardy;
					d2zdxdy = &table.d2smolardxdy; d2zdx2 = &table.d2smolardx2; d2zdy2 = &table.d2smolardy2;
					break;
				case iHmolar:
					z = &table.hmolar; dzdx = &table.dhmolardx; dzdy = &table.dhmolardy;
					d2zdxdy = &table.d2hmolardxdy; d2zdx2 = &table.d2hmolardx2; d2zdy2 = &table.d2hmolardy2;
					break;
				case iUmolar:
					z = &table.umolar; dzdx = &table.dumolardx; dzdy = &table.dumolardy;
					d2zdxdy = &table.d2umolardxdy; d2zdx2 = &table.d2umolardx2; d2zdy2 = &table.d2umolardy2;
					break;
				case iviscosity:
					z = &table.visc; break;
				case iconductivity:
					z = &table.cond; break;
				default:
					throw ValueError();
			}
		}
        LogPHTable single_phase_logph;
        LogPTTable single_phase_logpT;
        PureFluidSaturationTableData pure_saturation; // This will ultimately be split into pure and mixture backends which derive from this backend
        PhaseEnvelopeData phase_envelope;
        
        bool using_mole_fractions(void){return true;}
        bool using_mass_fractions(void){return false;}
        bool using_volu_fractions(void){return false;}
        void update(CoolProp::input_pairs input_pair, double Value1, double Value2){};
        void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions){this->AS->set_mole_fractions(mole_fractions);};
        void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){};
        const std::vector<long double> & get_mole_fractions(){return AS->get_mole_fractions();};
        CoolPropDbl calc_molar_mass(void){return AS->molar_mass();};
        virtual double evaluate_single_phase_phmolar(parameters output, std::size_t i, std::size_t j) = 0;
        virtual double evaluate_single_phase_pT(parameters output, std::size_t i, std::size_t j) = 0;
        virtual double evaluate_single_phase_phmolar_transport(parameters output, std::size_t i, std::size_t j) = 0;
        virtual double evaluate_single_phase_pT_transport(parameters output, std::size_t i, std::size_t j) = 0;
        virtual double evaluate_single_phase_phmolar_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) = 0;
        virtual double evaluate_single_phase_pT_derivative(parameters output, std::size_t i, std::size_t j, std::size_t Nx, std::size_t Ny) = 0;
        
        /// Returns the path to the tables that shall be written
        std::string path_to_tables(void);
        /// Load the tables from file; throws UnableToLoadException if there is a problem
        void load_tables();
        /// Build the tables
        void build_tables(){
            // Pure or pseudo-pure fluid
            if (AS->get_mole_fractions().size() == 1){
                pure_saturation.build(AS);
            }
            else{
                // Call function to actually construct the phase envelope
                AS->build_phase_envelope("");
                // Copy constructed phase envelope into this class
                phase_envelope = AS->get_phase_envelope_data();
                // Resize so that it will load properly
                pure_saturation.resize(pure_saturation.N);
            }
            single_phase_logph.build(AS); 
            single_phase_logpT.build(AS);
        }
        void pack_matrices(){
            single_phase_logph.pack();
            single_phase_logpT.pack();
            pure_saturation.pack();
            phase_envelope.pack();
        }
        /// Write the tables to file
        void write_tables();        
        
        CoolPropDbl calc_T(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iT, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return _T;
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iT, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_rhomolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iDmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iDmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iDmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_hmolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return _hmolar;
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iHmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iHmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_smolar(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar(iSmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT(iSmolar, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iHmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_cpmolar(void){
            if (using_single_phase_table){
                return calc_first_partial_deriv(iHmolar, iT, iP);
            }
            else{
                throw ValueError("Two-phase not possible for cpmolar currently");
            }
        }
		CoolPropDbl calc_cvmolar(void){
            if (using_single_phase_table){
                return calc_first_partial_deriv(iUmolar, iT, iDmolar);
            }
            else{
                throw ValueError("Two-phase not possible for cvmolar currently");
            }
        }
        
        CoolPropDbl calc_viscosity(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar_transport(iviscosity, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT_transport(iviscosity, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iviscosity, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_conductivity(void){
            if (using_single_phase_table){
                switch(selected_table){
                    case SELECTED_PH_TABLE: return evaluate_single_phase_phmolar_transport(iconductivity, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_PT_TABLE: return evaluate_single_phase_pT_transport(iconductivity, cached_single_phase_i, cached_single_phase_j);
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                return _HUGE; // not needed, will never be hit, just to make compiler happy
            }
            else{
                return pure_saturation.evaluate(iconductivity, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
        CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant){
            if (using_single_phase_table){
                CoolPropDbl dOf_dx, dOf_dy, dWrt_dx, dWrt_dy, dConstant_dx, dConstant_dy;

                // If a mass-based parameter is provided, get a conversion factor and change the key to the molar-based key
                double Of_conversion_factor = 1.0, Wrt_conversion_factor = 1.0, Constant_conversion_factor = 1.0;
                mass_to_molar(Of, Of_conversion_factor, AS->molar_mass());
                mass_to_molar(Wrt, Wrt_conversion_factor, AS->molar_mass());
                mass_to_molar(Constant, Constant_conversion_factor, AS->molar_mass());
                
                switch(selected_table){
                    case SELECTED_PH_TABLE: {
                        dOf_dx = evaluate_single_phase_phmolar_derivative(Of, cached_single_phase_i, cached_single_phase_j,1,0);
                        dOf_dy = evaluate_single_phase_phmolar_derivative(Of, cached_single_phase_i, cached_single_phase_j,0,1);
                        dWrt_dx = evaluate_single_phase_phmolar_derivative(Wrt, cached_single_phase_i, cached_single_phase_j,1,0);
                        dWrt_dy = evaluate_single_phase_phmolar_derivative(Wrt, cached_single_phase_i, cached_single_phase_j,0,1);
                        dConstant_dx = evaluate_single_phase_phmolar_derivative(Constant, cached_single_phase_i, cached_single_phase_j,1,0);
                        dConstant_dy = evaluate_single_phase_phmolar_derivative(Constant, cached_single_phase_i, cached_single_phase_j,0,1);
						break;
                    }
                    case SELECTED_PT_TABLE:{
                        dOf_dx = evaluate_single_phase_pT_derivative(Of, cached_single_phase_i, cached_single_phase_j,1,0);
                        dOf_dy = evaluate_single_phase_pT_derivative(Of, cached_single_phase_i, cached_single_phase_j,0,1);
                        dWrt_dx = evaluate_single_phase_pT_derivative(Wrt, cached_single_phase_i, cached_single_phase_j,1,0);
                        dWrt_dy = evaluate_single_phase_pT_derivative(Wrt, cached_single_phase_i, cached_single_phase_j,0,1);
                        dConstant_dx = evaluate_single_phase_pT_derivative(Constant, cached_single_phase_i, cached_single_phase_j,1,0);
                        dConstant_dy = evaluate_single_phase_pT_derivative(Constant, cached_single_phase_i, cached_single_phase_j,0,1);
						break;
                    }
                    case SELECTED_NO_TABLE: throw ValueError("table not selected");
                }
                double val = (dOf_dx*dConstant_dy-dOf_dy*dConstant_dx)/(dWrt_dx*dConstant_dy-dWrt_dy*dConstant_dx);
                return val*Of_conversion_factor/Wrt_conversion_factor;
            }
            else{
                return pure_saturation.evaluate(iconductivity, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
            
        };

        void check_tables(){
            if (!tables_loaded){
                try{
                    /// Try to load the tables if you can.
                    load_tables();
                    // Set the flag saying tables have been successfully loaded
                    tables_loaded = true;
                }
                catch(CoolProp::UnableToLoadError &){
                    /// Check directory size
                    std::string table_path = get_home_dir() + "/.CoolProp/Tables/";
                    #if defined(__ISWINDOWS__)
                        double directory_size_in_GB = CalculateDirSize(std::wstring(table_path.begin(), table_path.end()))/POW3(1024.0);
                    #else
                        double directory_size_in_GB = CalculateDirSize(table_path)/POW3(1024.0);
                    #endif
                    double allowed_size_in_GB = get_config_double(MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB);
                    if (get_debug_level() > 0){std::cout << "Tabular directory size is " << directory_size_in_GB << " GB\n";}
                    if (directory_size_in_GB > 1.5*allowed_size_in_GB){
                        throw DirectorySizeError(format("Maximum allowed tabular directory size is %g GB, you have exceeded 1.5 times this limit", allowed_size_in_GB));
                    }
                    else if (directory_size_in_GB > allowed_size_in_GB){
                        set_warning_string(format("Maximum allowed tabular directory size is %g GB, you have exceeded this limit", allowed_size_in_GB));
                    }
                    /// If you cannot load the tables, build them and then write them to file
                    build_tables();
                    pack_matrices();
                    write_tables();
                    /// Load the tables back into memory as a consistency check
                    load_tables();
                    // Set the flag saying tables have been successfully loaded
                    tables_loaded = true;
                }
            }
        };
        /** /brief calculate the derivative along the saturation curve, but only if quality is 0 or 1
         */
        CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1){
            if (AS->get_mole_fractions().size() > 1){throw ValueError("calc_first_saturation_deriv not available for mixtures");}
            if (std::abs(_Q) < 1e-6){
                return pure_saturation.first_saturation_deriv(Of1, Wrt1, 0, keyed_output(Wrt1), cached_saturation_iL);
            }
            else if (std::abs(_Q-1) < 1e-6){
                return pure_saturation.first_saturation_deriv(Of1, Wrt1, 1, keyed_output(Wrt1), cached_saturation_iV);
            }
            else{
                throw ValueError(format("Quality [%Lg] must be either 0 or 1 to within 1 ppm", _Q));
            }
        }
        CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant)
        {
            if (Of == iDmolar && Wrt == iHmolar && Constant == iP){
                CoolPropDbl rhoL = pure_saturation.evaluate(iDmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl rhoV = pure_saturation.evaluate(iDmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl hL = pure_saturation.evaluate(iHmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl hV = pure_saturation.evaluate(iHmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
                return -POW2(rhomolar())*(1/rhoV - 1/rhoL)/(hV - hL);
            }
            else if (Of == iDmass && Wrt == iHmass && Constant == iP){
                return first_two_phase_deriv(iDmolar, iHmolar, iP)*POW2(molar_mass());
            }
            else if (Of == iDmolar && Wrt == iP && Constant == iHmolar){
                // v = 1/rho; dvdrho = -rho^2; dvdrho = -1/rho^2
                CoolPropDbl rhoL = pure_saturation.evaluate(iDmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl rhoV = pure_saturation.evaluate(iDmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl hL = pure_saturation.evaluate(iHmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl hV = pure_saturation.evaluate(iHmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
                CoolPropDbl dvdrhoL = -1/POW2(rhoL);
                CoolPropDbl dvdrhoV = -1/POW2(rhoV);
                CoolPropDbl dvL_dp = dvdrhoL*pure_saturation.first_saturation_deriv(iDmolar, iP, 0, _p, cached_saturation_iL);
                CoolPropDbl dvV_dp = dvdrhoV*pure_saturation.first_saturation_deriv(iDmolar, iP, 1, _p, cached_saturation_iV);
                CoolPropDbl dhL_dp = pure_saturation.first_saturation_deriv(iHmolar, iP, 0, _p, cached_saturation_iL);
                CoolPropDbl dhV_dp = pure_saturation.first_saturation_deriv(iHmolar, iP, 1, _p, cached_saturation_iV);
                CoolPropDbl dxdp_h = (Q()*dhV_dp + (1 - Q())*dhL_dp)/(hL - hV);
                CoolPropDbl dvdp_h = dvL_dp + dxdp_h*(1/rhoV - 1/rhoL) + Q()*(dvV_dp - dvL_dp);
                return -POW2(rhomolar())*dvdp_h;
            }
            else if (Of == iDmass && Wrt == iP && Constant == iHmass){
                return first_two_phase_deriv(iDmolar, iP, iHmolar)*molar_mass();
            }
            else{
                throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
            }
        }
};

} /* namespace CoolProp*/

#endif
