
#include "TabularBackends.h"
#include "CoolProp.h"

namespace CoolProp{
    
void GriddedTableBackend::build_tables(tabular_types type)
{
    std::size_t Nx = 20, Ny = 20;
    long double xmin, xmax, ymin, ymax, x, y;
    bool logy, logx;
    
    parameters xkey, ykey;
    single_phase.resize(Nx, Ny);
    
    switch(type){
        case LOGPH_TABLE:
        {
            xkey = iHmolar;
            ykey = iP;
            logy = true;
            logx = false;
            
            // ---------------------------------
            // Calculate the limits of the table
            // ---------------------------------
            
            // Minimum enthalpy is the saturated liquid enthalpy
            AS->update(QT_INPUTS, 0, AS->Ttriple());
            xmin = AS->hmolar();
            ymin = AS->p();
            
            // Check both the enthalpies at the Tmax isotherm to see whether to use low or high pressure
            AS->update(PT_INPUTS, 1e-10, AS->Tmax());
            long double xmax1 = AS->hmolar();
            AS->update(PT_INPUTS, AS->pmax(), AS->Tmax());
            long double xmax2 = AS->hmolar();
            xmax = std::min(xmax1, xmax2);
            
            ymax = AS->pmax();
            
            break;
        }
        default:
        {
            throw ValueError(format(""));
        }
    }
    if (get_debug_level() > 5){
        std::cout << format("***********************************************\n");
        std::cout << format(" Single-Phase Table (%s) \n", AS->name().c_str());
        std::cout << format("***********************************************\n");
    }
    // ------------------------
    // Actually build the table
    // ------------------------
    for (std::size_t i = 0; i < Nx; ++i)
    {
        // Calculate the x value
        if (logx){
            // Log spaced
            x = exp(log(xmin) + (log(xmax) - log(xmin))/(Nx-1)*i);
        }
        else{
            // Linearly spaced
            x = xmin + (xmax - xmin)/(Nx-1)*i;
        }
        for (std::size_t j = 0; j < Ny; ++j)
        {
            // Calculate the x value
            if (logy){
                // Log spaced
                y = exp(log(ymin) + (log(ymax/ymin))/(Ny-1)*j);
            }
            else{
                // Linearly spaced
                y = ymin + (ymax - ymin)/(Ny-1)*j;
            }
            
            if (get_debug_level() > 5){std::cout << "x:" << x << " y: " << y;}
            if (xkey == iHmolar && ykey == iP){
                try{
                    // Update the state
                    AS->update(HmolarP_INPUTS, x, y);
                    
                    // Skip two-phase states - they will remain as _HUGE holes in the table
                    if (AS->phase() == iphase_twophase){ 
                        if (get_debug_level() > 5){std::cout << "2Phase" << std::endl;}
                        continue;
                        };
                    
                    // Calculate each of the parameters, and their derivatives with respect to the input variables
                    single_phase.T[i][j] = AS->T();
                    single_phase.p[i][j] = AS->p();
                    single_phase.hmolar[i][j] = AS->hmolar();
                    single_phase.smolar[i][j] = AS->smolar();
                    single_phase.rhomolar[i][j] = AS->rhomolar();
                    if (get_debug_level() > 5){std::cout << "OK" << std::endl;}
                }
                catch(std::exception &e){
                    if (get_debug_level() > 5){std::cout << e.what() << std::endl;}
                    // If there is an error, go to next one
                    continue;
                }
            }
            else{
                throw ValueError("Cannot construct this type of table");
            }
            
            // Calculate the derivatives
        }
    }
}
    
} /* namespace CoolProp */