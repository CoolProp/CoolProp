#ifndef TABULAR_BACKENDS_H
#define TABULAR_BACKENDS_H

#include "AbstractState.h"
#include "msgpack.hpp"
#include <msgpack/fbuffer.hpp>
#include "crossplatform_shared_ptr.h"
#include "Exceptions.h"

class UnableToLoadException : public CoolProp::CoolPropBaseError{
    public:
        UnableToLoadException(std::string e){err = e;}
        virtual const char* what() const throw(){ return err.c_str(); }
};

/**
 * ***MAGIC WARNING**!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_MATRICES X(T) X(p) X(rhomolar) X(hmolar) X(smolar) X(dTdx) X(dTdy) X(dpdx) X(dpdy) X(drhomolardx) X(drhomolardy) X(dhmolardx) X(dhmolardy) X(dsmolardx) X(dsmolardy) X(d2Tdx2) X(d2Tdxdy) X(d2Tdy2) X(d2pdx2) X(d2pdxdy) X(d2pdy2) X(d2rhomolardx2) X(d2rhomolardxdy) X(d2rhomolardy2) X(d2hmolardx2) X(d2hmolardxdy) X(d2hmolardy2) X(d2smolardx2) X(d2smolardxdy) X(d2smolardy2)

namespace CoolProp{

/** \brief This class holds the data for a single-phase interpolation table that is regularly spaced
 * 
 * It contains very few members or methods, mostly it just holds the data
 */
struct SinglePhaseGriddedTableData{
    std::size_t Nx, Ny;
    CoolProp::parameters xkey, ykey;
    shared_ptr<CoolProp::AbstractState> AS;
    std::vector<double> xvec, yvec;
    bool logx, logy;
    double xmin, ymin, xmax, ymax;
    
    SinglePhaseGriddedTableData(){Nx = 100; Ny = 100; revision = 0;}
    
    /* Use X macros to auto-generate the variables; each will look something like: std::vector< std::vector<double> > T; */
    #define X(name) std::vector< std::vector<double> > name;
    LIST_OF_MATRICES
    #undef X
    int revision;
    std::map<std::string, std::vector<std::vector<double> > > matrices;
    
    MSGPACK_DEFINE(revision, matrices, xmin, xmax, ymin, ymax); // write the member variables that you want to pack
    /// Resize all the matrices
    void resize(std::size_t Nx, std::size_t Ny){
        /* Use X macros to auto-generate the code; each will look something like: T.resize(Nx, std::vector<double>(Ny, _HUGE)); */
        #define X(name) name.resize(Nx, std::vector<double>(Ny, _HUGE));
        LIST_OF_MATRICES
        #undef X
        
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
    /// Take all the matrices that are in the class and pack them into the matrices map for easy unpacking using msgpack
    void pack_matrices(){
        /* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
        #define X(name) matrices.insert(std::pair<std::string, std::vector<std::vector<double> > >(#name, name));
        LIST_OF_MATRICES
        #undef X
    };
    /// Take all the matrices that are in the class and pack them into the matrices map for easy unpacking using msgpack
    void unpack_matrices(){
        /* Use X macros to auto-generate the unpacking code; each will look something like: T = *(matrices.find("T")).second */
        #define X(name) name = (matrices.find(#name))->second;
        LIST_OF_MATRICES
        #undef X
        
        Nx = T.size(); Ny = T[0].size();
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
        bisect_vector(yvec, y, j);
    }
    /// Find the nearest good neighbor node for inputs that are the same as the grid inputs
    /// If the straightforward node (i,j) obtained by bisection is no good, find its nearest good node
    void find_native_nearest_good_neighbor(double x, double y, std::size_t &i, std::size_t &j){
        // Get the node that is closest
        find_native_nearest_neighbor(x,y,i,j);
        // Check whether found node is good
        if (!ValidNumber(T[i][j])){
            // If not, find its nearest good neighbor (nearest good neighbors are precalculated and cached)
        }
    }
    /// Find the nearest cell with lower left coordinate (i,j) where (i,j) is a good node, and so are (i+1,j), (i,j+1), (i+1,j+1)
    /// This is needed for bicubic interpolation
    void find_native_nearest_good_cell(double x, double y, std::size_t &i, std::size_t &j){
    }
    protected:
        /// Build this table; it is protected to remind derived classes to set the necessary variables
        void build(shared_ptr<CoolProp::AbstractState> &AS);
};

/// This class holds the single-phase data for a log(p)-h gridded table
class LogPHTable : public SinglePhaseGriddedTableData
{
    public:
        LogPHTable(){
            xkey = iHmolar; ykey = iP; logy = true; logx = false;
        };
        /// Calculate the limits of the table
        void build(shared_ptr<CoolProp::AbstractState> &AS){
            
            // Minimum enthalpy is the saturated liquid enthalpy
            AS->update(QT_INPUTS, 0, AS->Ttriple());
            xmin = AS->hmolar(); ymin = AS->p();
            
            // Check both the enthalpies at the Tmax isotherm to see whether to use low or high pressure
            AS->update(PT_INPUTS, 1e-10, AS->Tmax());
            long double xmax1 = AS->hmolar();
            AS->update(PT_INPUTS, AS->pmax(), AS->Tmax());
            long double xmax2 = AS->hmolar();
            xmax = std::min(xmax1, xmax2);
            
            ymax = AS->pmax();
            
            // Call the base-class build function since we have set the necessary variables
            SinglePhaseGriddedTableData::build(AS);
        };
};
/// This class holds the single-phase data for a log(p)-T gridded table
class LogPTTable : public SinglePhaseGriddedTableData
{
    public:
        LogPTTable(){
            xkey = iT; ykey = iP; logy = true; logx = false;
        };
        /// Calculate the limits of the table
        void build(shared_ptr<CoolProp::AbstractState> &AS){
            
            xmin = AS->Ttriple(); ymin = AS->p_triple();
            xmax = AS->Tmax(); ymax = AS->pmax();
            
            // Call the base-class build function since we have set the necessary variables
            SinglePhaseGriddedTableData::build(AS);
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
        shared_ptr<CoolProp::AbstractState> AS;
    public:
        LogPHTable single_phase_logph;
        LogPTTable single_phase_logpT;
        
        std::string get_backend_name(void){return "TabularBackend";};
        bool using_mole_fractions(void){return true;}
        bool using_mass_fractions(void){return false;}
        bool using_volu_fractions(void){return false;}
        void update(CoolProp::input_pairs input_pair, double Value1, double Value2){};
        void set_mole_fractions(const std::vector<long double> &mole_fractions){};
        void set_mass_fractions(const std::vector<long double> &mass_fractions){};
        const std::vector<long double> & get_mole_fractions(){throw NotImplementedError("get_mole_fractions not implemented for TTSE");};
        
        /// Returns the path to the tables that shall be written
        std::string path_to_tables(void);
        /// Load the tables from file; throws UnableToLoadException if there is a problem
        void load_tables();
        /// Build the tables
        void build_tables(){
            single_phase_logph.build(AS); 
            //single_phase_logpT.build(AS);
        }
        void pack_matrices(){
            single_phase_logph.pack_matrices();
            //single_phase_logpT.pack_matrices();
        }
        /// Write the tables to file
        void write_tables();
        
        TabularBackend(shared_ptr<CoolProp::AbstractState> AS){
            // Grab onto the pointer to the class
            this->AS = AS;
            
            try{
                /// Try to load the tables if you can.
                load_tables();
            }
            catch(UnableToLoadException &){
                /// If you cannot load the tables, build them and then write them to file
                build_tables();
                pack_matrices();
                write_tables();
                /// Load the tables back into memory as a consistency check
                load_tables();
            }
        }
};

} /* namespace CoolProp*/

#endif
