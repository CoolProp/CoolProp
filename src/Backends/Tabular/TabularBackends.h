#ifndef TABULAR_BACKENDS_H
#define TABULAR_BACKENDS_H

#include "AbstractState.h"

namespace CoolProp{

struct SinglePhaseGriddedTableData{
    CoolProp::parameters xkey, ykey;
    std::vector< std::vector<double> > T, p, rhomolar, hmolar, smolar;
    std::vector< std::vector<double> > dTdx, dTdy, dpdx, dpdy, drhomolardx, drhomolardy, dhmolardx, dhmolardy, dsmolardx, dsmolardy;
    
    void resize(std::size_t Nx, std::size_t Ny){
        // State variables
        T.resize(Nx, std::vector<double>(Ny, _HUGE));
        p.resize(Nx, std::vector<double>(Ny, _HUGE));
        rhomolar.resize(Nx, std::vector<double>(Ny, _HUGE));
        hmolar.resize(Nx, std::vector<double>(Ny, _HUGE));
        smolar.resize(Nx, std::vector<double>(Ny, _HUGE));
        
        // First derivatives of state variables
        dTdx.resize(Nx, std::vector<double>(Ny, _HUGE));
        dpdx.resize(Nx, std::vector<double>(Ny, _HUGE));
        drhomolardx.resize(Nx, std::vector<double>(Ny, _HUGE));
        dhmolardx.resize(Nx, std::vector<double>(Ny, _HUGE));
        dsmolardx.resize(Nx, std::vector<double>(Ny, _HUGE));
        dTdy.resize(Nx, std::vector<double>(Ny, _HUGE));
        dpdy.resize(Nx, std::vector<double>(Ny, _HUGE));
        drhomolardy.resize(Nx, std::vector<double>(Ny, _HUGE));
        dhmolardy.resize(Nx, std::vector<double>(Ny, _HUGE));
        dsmolardy.resize(Nx, std::vector<double>(Ny, _HUGE));
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
    void set_mole_fractions(const std::vector<long double> &mole_fractions){};
    void set_mass_fractions(const std::vector<long double> &mass_fractions){};
    const std::vector<long double> & get_mole_fractions(){throw NotImplementedError("get_mole_fractions not implemented for TTSE");};
    void build_tables(tabular_types type);
    void load_tables(void);
    void bounding_curves(void);
};
    
class TTSEBackend : public GriddedTableBackend
{
public:
    TTSEBackend(AbstractState &AS){
        this->AS = &AS; build_tables(GriddedTableBackend::LOGPH_TABLE);};
};

} /* namespace CoolProp*/

#endif
