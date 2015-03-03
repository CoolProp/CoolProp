

#include "TabularBackends.h"
#include "CoolProp.h"
#include <sstream>
#include "time.h"
void CoolProp::PureFluidSaturationTableData::build(shared_ptr<CoolProp::AbstractState> &AS){
    const bool debug = get_debug_level() > 5 || false;
    if (debug){
        std::cout << format("***********************************************\n");
        std::cout << format(" Saturation Table (%s) \n", AS->name().c_str());
        std::cout << format("***********************************************\n");
    }
    resize(N);
    // ------------------------
    // Actually build the table
    // ------------------------
    CoolPropDbl p, pmin = AS->p_triple()*1.001, pmax = 0.9999999*AS->p_critical();
    for (std::size_t i = 0; i < N-1; ++i)
    {
        // Log spaced
        p = exp(log(pmin) + (log(pmax) - log(pmin))/(N-1)*i);
        
        // Saturated liquid
        try{
            AS->update(PQ_INPUTS, p, 0);
            pL[i] = p; TL[i] = AS->T();  rhomolarL[i] = AS->rhomolar(); 
            hmolarL[i] = AS->hmolar(); smolarL[i] = AS->smolar(); umolarL[i] = AS->umolar();
            logpL[i] = log(p); logrhomolarL[i] = log(rhomolarL[i]);
        }
        catch(std::exception &e){
            // That failed for some reason, go to the next pair
            if (debug){std::cout << " " << e.what() << std::endl;}
            continue;
        }
        // Saturated vapor
        try{
            AS->update(PQ_INPUTS, p, 1);
            pV[i] = p; TV[i] = AS->T(); rhomolarV[i] = AS->rhomolar();
            hmolarV[i] = AS->hmolar(); smolarV[i] = AS->smolar(); umolarV[i] = AS->umolar();
            logpV[i] = log(p); logrhomolarV[i] = log(rhomolarV[i]);
        }
        catch(std::exception &e){
            // That failed for some reason, go to the next pair
            if (debug){std::cout << " " << e.what() << std::endl;}
            continue;
        }
    }
    // Last point is at the critical point
    AS->update(PQ_INPUTS, AS->p_critical(), 1);
    std::size_t i = N-1;
    pV[i] = p; TV[i] = AS->T(); rhomolarV[i] = AS->rhomolar();
    hmolarV[i] = AS->hmolar(); smolarV[i] = AS->smolar(); umolarV[i] = AS->umolar();
    pL[i] = p; TL[i] = AS->T();  rhomolarL[i] = AS->rhomolar(); 
    hmolarL[i] = AS->hmolar(); smolarL[i] = AS->smolar(); umolarL[i] = AS->umolar();
}
    
void CoolProp::SinglePhaseGriddedTableData::build(shared_ptr<CoolProp::AbstractState> &AS)
{
    CoolPropDbl x, y;
    const bool debug = get_debug_level() > 5 || true;

    resize(Nx, Ny);
    
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
        xvec[i] = x;
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
            yvec[j] = y;
            
            if (debug){std::cout << "x: " << x << " y: " << y << std::endl;}
            
            // Generate the input pair
            CoolPropDbl v1, v2;
            input_pairs input_pair = generate_update_pair(xkey, x, ykey, y, v1, v2);
            
            // --------------------
            //   Update the state
            // --------------------
            try{
                AS->update(input_pair, v1, v2);
                if (!ValidNumber(AS->rhomolar())){
                    throw ValueError("rhomolar is invalid");
                }
            }
            catch(std::exception &e){
                // That failed for some reason, go to the next pair
                if (debug){std::cout << " " << e.what() << std::endl;}
                continue;
            }
            
            // Skip two-phase states - they will remain as _HUGE holes in the table
            if (AS->phase() == iphase_twophase){ 
                if (debug){std::cout << " 2Phase" << std::endl;}
                continue;
            };
            
            // --------------------
            //   State variables
            // --------------------
            T[i][j] = AS->T();
            p[i][j] = AS->p();
            rhomolar[i][j] = AS->rhomolar();
            hmolar[i][j] = AS->hmolar();
            smolar[i][j] = AS->smolar();
            
            // -------------------------
            //   Transport properties
            // -------------------------
            try{
                visc[i][j] = AS->viscosity();
                cond[i][j] = AS->conductivity();
            }
            catch(std::exception &){}
            
            // ----------------------------------------
            //   First derivatives of state variables
            // ----------------------------------------
            dTdx[i][j] = AS->first_partial_deriv(iT, xkey, ykey);
            dTdy[i][j] = AS->first_partial_deriv(iT, ykey, xkey);
            dpdx[i][j] = AS->first_partial_deriv(iP, xkey, ykey);
            dpdy[i][j] = AS->first_partial_deriv(iP, ykey, xkey);
            drhomolardx[i][j] = AS->first_partial_deriv(iDmolar, xkey, ykey);
            drhomolardy[i][j] = AS->first_partial_deriv(iDmolar, ykey, xkey);
            dhmolardx[i][j] = AS->first_partial_deriv(iHmolar, xkey, ykey);
            dhmolardy[i][j] = AS->first_partial_deriv(iHmolar, ykey, xkey);
            dsmolardx[i][j] = AS->first_partial_deriv(iSmolar, xkey, ykey);
            dsmolardy[i][j] = AS->first_partial_deriv(iSmolar, ykey, xkey);
            
            // ----------------------------------------
            //   Second derivatives of state variables
            // ----------------------------------------
            d2Tdx2[i][j] = AS->second_partial_deriv(iT, xkey, ykey, xkey, ykey);
            d2Tdxdy[i][j] = AS->second_partial_deriv(iT, xkey, ykey, ykey, xkey);
            d2Tdy2[i][j] = AS->second_partial_deriv(iT, ykey, xkey, ykey, xkey);
            d2pdx2[i][j] = AS->second_partial_deriv(iP, xkey, ykey, xkey, ykey);
            d2pdxdy[i][j] = AS->second_partial_deriv(iP, xkey, ykey, ykey, xkey);
            d2pdy2[i][j] = AS->second_partial_deriv(iP, ykey, xkey, ykey, xkey);
            d2rhomolardx2[i][j] = AS->second_partial_deriv(iDmolar, xkey, ykey, xkey, ykey);
            d2rhomolardxdy[i][j] = AS->second_partial_deriv(iDmolar, xkey, ykey, ykey, xkey);
            d2rhomolardy2[i][j] = AS->second_partial_deriv(iDmolar, ykey, xkey, ykey, xkey);
            d2hmolardx2[i][j] = AS->second_partial_deriv(iHmolar, xkey, ykey, xkey, ykey);
            d2hmolardxdy[i][j] = AS->second_partial_deriv(iHmolar, xkey, ykey, ykey, xkey);
            d2hmolardy2[i][j] = AS->second_partial_deriv(iHmolar, ykey, xkey, ykey, xkey);
            d2smolardx2[i][j] = AS->second_partial_deriv(iSmolar, xkey, ykey, xkey, ykey);
            d2smolardxdy[i][j] = AS->second_partial_deriv(iSmolar, xkey, ykey, ykey, xkey);
            d2smolardy2[i][j] = AS->second_partial_deriv(iSmolar, ykey, xkey, ykey, xkey);
        }
    }
}
std::string CoolProp::TabularBackend::path_to_tables(void){
    std::vector<std::string> fluids = AS->fluid_names();
    return get_home_dir() + "/.CoolProp/Tables/" + AS->backend_name() + "(" + strjoin(AS->fluid_names(),"&") + ")";
}
void CoolProp::TabularBackend::write_tables(){
    std::string path = path_to_tables();
    make_dirs(path);
    {
        std::ofstream ofs(std::string(path + "/single_phase_logph.bin").c_str(), std::ofstream::binary);
        msgpack::pack(ofs, single_phase_logph);
    }
    {
        std::ofstream ofs(std::string(path + "/single_phase_logpT.bin").c_str(), std::ofstream::binary);
        msgpack::pack(ofs, single_phase_logpT);
    }
    {
        std::ofstream ofs(std::string(path + "/pure_saturation.bin").c_str(), std::ofstream::binary);
        msgpack::pack(ofs, pure_saturation);
    }
}
void CoolProp::TabularBackend::load_tables(){
    std::string path_to_tables = this->path_to_tables();
    load_table(single_phase_logph, path_to_tables, "single_phase_logph.bin");
    load_table(single_phase_logpT, path_to_tables, "single_phase_logpT.bin");
    load_table(pure_saturation, path_to_tables, "pure_saturation.bin");
}