

#include "TabularBackends.h"
#include "CoolProp.h"
#include <sstream>

namespace CoolProp{
    
void GriddedTableBackend::build_tables(tabular_types type)
{
    std::size_t Nx = 200, Ny = 200;
    CoolPropDbl xmin, xmax, ymin, ymax, x, y;
    bool logy, logx;
    const bool debug = get_debug_level() > 5 || false;
    
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
            CoolPropDbl xmax1 = AS->hmolar();
            AS->update(PT_INPUTS, AS->pmax(), AS->Tmax());
            CoolPropDbl xmax2 = AS->hmolar();
            xmax = std::min(xmax1, xmax2);
            
            ymax = AS->pmax();
            
            break;
        }
        default:
        {
            throw ValueError(format(""));
        }
    }
    if (debug){
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
            
            if (debug){std::cout << "x: " << x << " y: " << y;}
            
            if (xkey == iHmolar && ykey == iP)
            {
                // --------------------
                //   Update the state
                // --------------------
                try{
                    AS->update(HmolarP_INPUTS, x, y);
                }
                catch(std::exception &e){
                    // That failed for some reason, go to the next pair
                    if (debug){std::cout << " " << e.what() << std::endl;}
                    continue;
                }
                
                if (debug){std::cout << " OK" << std::endl;}
            }
            else{
                throw ValueError("Cannot construct this type of table (yet)");
            }
            
            // Skip two-phase states - they will remain as _HUGE holes in the table
            if (AS->phase() == iphase_twophase){ 
                if (debug){std::cout << " 2Phase" << std::endl;}
                continue;
            };
            
            // --------------------
            //   State variables
            // --------------------
            single_phase.T[i][j] = AS->T();
            single_phase.p[i][j] = AS->p();
            single_phase.rhomolar[i][j] = AS->rhomolar();
            single_phase.hmolar[i][j] = AS->hmolar();
            single_phase.smolar[i][j] = AS->smolar();
            
            // ----------------------------------------
            //   First derivatives of state variables
            // ----------------------------------------
            single_phase.dTdx[i][j] = AS->first_partial_deriv(iT, xkey, ykey);
            single_phase.dTdy[i][j] = AS->first_partial_deriv(iT, ykey, xkey);
            single_phase.dpdx[i][j] = AS->first_partial_deriv(iP, xkey, ykey);
            single_phase.dpdy[i][j] = AS->first_partial_deriv(iP, ykey, xkey);
            single_phase.drhomolardx[i][j] = AS->first_partial_deriv(iDmolar, xkey, ykey);
            single_phase.drhomolardy[i][j] = AS->first_partial_deriv(iDmolar, ykey, xkey);
            single_phase.dhmolardx[i][j] = AS->first_partial_deriv(iHmolar, xkey, ykey);
            single_phase.dhmolardy[i][j] = AS->first_partial_deriv(iHmolar, ykey, xkey);
            single_phase.dsmolardx[i][j] = AS->first_partial_deriv(iSmolar, xkey, ykey);
            single_phase.dsmolardy[i][j] = AS->first_partial_deriv(iSmolar, ykey, xkey);
            
            // ----------------------------------------
            //   Second derivatives of state variables
            // ----------------------------------------
            single_phase.d2Tdx2[i][j] = AS->second_partial_deriv(iT, xkey, ykey, xkey, ykey);
            single_phase.d2Tdxdy[i][j] = AS->second_partial_deriv(iT, xkey, ykey, ykey, xkey);
            single_phase.d2Tdy2[i][j] = AS->second_partial_deriv(iT, ykey, xkey, ykey, xkey);
            single_phase.d2pdx2[i][j] = AS->second_partial_deriv(iP, xkey, ykey, xkey, ykey);
            single_phase.d2pdxdy[i][j] = AS->second_partial_deriv(iP, xkey, ykey, ykey, xkey);
            single_phase.d2pdy2[i][j] = AS->second_partial_deriv(iP, ykey, xkey, ykey, xkey);
            single_phase.d2rhomolardx2[i][j] = AS->second_partial_deriv(iDmolar, xkey, ykey, xkey, ykey);
            single_phase.d2rhomolardxdy[i][j] = AS->second_partial_deriv(iDmolar, xkey, ykey, ykey, xkey);
            single_phase.d2rhomolardy2[i][j] = AS->second_partial_deriv(iDmolar, ykey, xkey, ykey, xkey);
            single_phase.d2hmolardx2[i][j] = AS->second_partial_deriv(iHmolar, xkey, ykey, xkey, ykey);
            single_phase.d2hmolardxdy[i][j] = AS->second_partial_deriv(iHmolar, xkey, ykey, ykey, xkey);
            single_phase.d2hmolardy2[i][j] = AS->second_partial_deriv(iHmolar, ykey, xkey, ykey, xkey);
            single_phase.d2smolardx2[i][j] = AS->second_partial_deriv(iSmolar, xkey, ykey, xkey, ykey);
            single_phase.d2smolardxdy[i][j] = AS->second_partial_deriv(iSmolar, xkey, ykey, ykey, xkey);
            single_phase.d2smolardy2[i][j] = AS->second_partial_deriv(iSmolar, ykey, xkey, ykey, xkey);
        }
    }
}

void GriddedTableBackend::write_tables(const std::string &directory){
    std::ofstream ofs("single_phase.bin", std::ofstream::binary);
    msgpack::pack(ofs, single_phase);
    ofs.close();
}
void GriddedTableBackend::load_tables(const std::string &directory){
    
    std::ifstream ifs("single_phase.bin", std::ifstream::binary);
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    msgpack::unpacked upd;
    std::string sbuffer = buffer.str();
    std::size_t N = sbuffer.size();
    msgpack::unpack(upd, sbuffer.c_str(), N);
    msgpack::object deserialized = upd.get();
    deserialized.convert(&single_phase);
}
    
} /* namespace CoolProp */