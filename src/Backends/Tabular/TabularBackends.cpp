
#include "TabularBackends.h"

namespace CoolProp{
    
void GriddedTableBackend::build_tables(tabular_types type)
{
    std::size_t Nx = 200, Ny = 200;
    long double xmin, xmax, ymin, ymax, x, y;
    bool logy, logx;
    
    switch(type){
        case LOGPH_TABLE:
        {
            parameters xkey, ykey;
            
            // ---------------------------------
            // Calculate the limits of the table
            // ---------------------------------
            
            // Minimum enthalpy is the saturated liquid enthalpy
            AS->update(QT_INPUTS, 0, AS->Ttriple());
            xmin = AS->hmolar();
            ymin = log(AS->p());
            
            // Check both the enthalpies at the Tmax isotherm to see whether to use low or high pressure
            AS->update(PT_INPUTS, 1e-10, AS->Tmax());
            long double xmax1 = AS->hmolar();
            AS->update(PT_INPUTS, AS->pmax(), AS->Tmax());
            long double xmax2 = AS->hmolar();
            xmax = std::max(xmax1, xmax2);
            
            ymax = log(AS->p());
            xkey = iHmolar;
            ykey = iP;
            logy = true;
            logx = false;
            break;
        }
        default:
        {
            throw ValueError(format(""));
        }
    }
    // ------------------------
    // Actually build the table
    // ------------------------
    for (std::size_t i = 0; i < Nx; ++i)
    {
        // Calculate the x value
        x = xmin + (xmax - xmin)/(Nx-1)*i;
        if (logx){
            x = log(x);
        }
        for (std::size_t j = 0; j < Ny; ++j)
        {
            // Calculate the y value
            y = ymin + (ymax - ymin)/(Ny-1)*j;
            if (logy){
                y = log(y);
            }
            
            // Update the state
            
            // Calculate the derivatives
            
        }
    }
}
    
} /* namespace CoolProp */