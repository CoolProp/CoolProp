#ifndef TABULAR_BACKENDS_H
#define TABULAR_BACKENDS_H

#include "AbstractState.h"
#include "msgpack.hpp"
#include <msgpack/fbuffer.hpp>

/**
 * ***MAGIC WARNING**!! X Macros in use
 * See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 */
#define LIST_OF_MATRICES X(T) X(p) X(rhomolar) X(hmolar) X(smolar) X(dTdx) X(dTdy) X(dpdx) X(dpdy) X(drhomolardx) X(drhomolardy) X(dhmolardx) X(dhmolardy) X(dsmolardx) X(dsmolardy) X(d2Tdx2) X(d2Tdxdy) X(d2Tdy2) X(d2pdx2) X(d2pdxdy) X(d2pdy2) X(d2rhomolardx2) X(d2rhomolardxdy) X(d2rhomolardy2) X(d2hmolardx2) X(d2hmolardxdy) X(d2hmolardy2) X(d2smolardx2) X(d2smolardxdy) X(d2smolardy2)

namespace CoolProp{

struct SinglePhaseGriddedTableData{
    CoolProp::parameters xkey, ykey;
    
    /* Use X macros to auto-generate the variables; each will look something like: std::vector< std::vector<double> > T; */
    #define X(name) std::vector< std::vector<double> > name;
    LIST_OF_MATRICES
    #undef X
    int revision;
    std::map<std::string, std::vector<std::vector<double> > > matrices;
    
    MSGPACK_DEFINE(revision, matrices); // write the member variables that you want to pack
    void resize(std::size_t Nx, std::size_t Ny){
        /* Use X macros to auto-generate the variables; each will look something like: T.resize(Nx, std::vector<double>(Ny, _HUGE)); */
        #define X(name) name.resize(Nx, std::vector<double>(Ny, _HUGE));
        LIST_OF_MATRICES
        #undef X
    };
    void pack_matrices(){
        /* Use X macros to auto-generate the packing code; each will look something like: matrices.insert(std::pair<std::vector<std::vector<double> > >("T", T)); */
        #define X(name) matrices.insert(std::pair<std::string, std::vector<std::vector<double> > >(#name, name));
        LIST_OF_MATRICES
        #undef X
    };
};
class GriddedTableBackend : public AbstractState
{
    protected:
    enum tabular_types {LOGPH_TABLE};
    AbstractState *AS;
    public:
    SinglePhaseGriddedTableData single_phase;
    
    bool using_mole_fractions(void){return true;}
    bool using_mass_fractions(void){return false;}
    bool using_volu_fractions(void){return false;}
    
    void update(CoolProp::input_pairs input_pair, double Value1, double Value2){};
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions){};
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){};
    const std::vector<CoolPropDbl> & get_mole_fractions(){throw NotImplementedError("get_mole_fractions not implemented for TTSE");};
    void build_tables(tabular_types type);
    void write_tables(const std::string &directory);
    void load_tables(const std::string &directory);
    void bounding_curves(void);
    
};
    
class TTSEBackend : public GriddedTableBackend
{
public:
    TTSEBackend(AbstractState &AS){
        // Hold a pointer to the abstract state
        this->AS = &AS; 
        build_tables(GriddedTableBackend::LOGPH_TABLE);
        
        std::string path = "here";
        single_phase.pack_matrices();
        write_tables(path);
        load_tables(path);
        int rrgregregstrs = 4;
    };
};

} /* namespace CoolProp*/

#endif
