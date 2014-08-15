#ifndef TABULAR_BACKENDS_H
#define TABULAR_BACKENDS_H

#include "AbstractState.h"

namespace CoolProp{

class GriddedTableBackend : public AbstractState
{
    protected:
    enum tabular_types {LOGPH_TABLE};
    AbstractState *AS;
    public:
    bool using_mole_fractions(void){return true;}
    bool using_mass_fractions(void){return false;}
    bool using_volu_fractions(void){return false;}
    
    void update(long input_pair, double Value1, double Value2){};
    void set_mole_fractions(const std::vector<long double> &mole_fractions){};
    void set_mass_fractions(const std::vector<long double> &mass_fractions){};
    
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
