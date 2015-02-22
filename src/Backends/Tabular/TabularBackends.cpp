

#include "TabularBackends.h"
#include "CoolProp.h"
#include <sstream>

namespace CoolProp{
    
void SinglePhaseGriddedTableData::build(shared_ptr<CoolProp::AbstractState> &AS)
{
    
    long double x, y;
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
            T[i][j] = AS->T();
            p[i][j] = AS->p();
            rhomolar[i][j] = AS->rhomolar();
            hmolar[i][j] = AS->hmolar();
            smolar[i][j] = AS->smolar();
            
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
std::string TabularBackend::path_to_tables(void){
    std::vector<std::string> fluids = AS->fluid_names();
    return get_home_dir() + "/.CoolProp/Tables/" + AS->backend_name() + "(" + strjoin(AS->fluid_names(),"&") + ")";
}
void TabularBackend::write_tables(){
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
}
void TabularBackend::load_tables(){
    std::string path_to_logph = path_to_tables() + "/single_phase_logph.bin";
    std::ifstream ifs(path_to_logph.c_str(), std::ifstream::binary);
    
    if ( (ifs.rdstate() & std::ifstream::failbit ) != 0 ){
        if (get_debug_level() > 0){std::cout << format("Error loading table %s", path_to_logph.c_str()) << std::endl;}
        throw UnableToLoadException(format("Error loading table %s", path_to_logph.c_str()));
    }
    
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    msgpack::unpacked upd;
    std::string sbuffer = buffer.str();
    std::size_t N = sbuffer.size();
    if ( N == 0 ){
        if (get_debug_level() > 0){std::cout << format("No data was read from table %s", path_to_logph.c_str()) << std::endl;}
        throw UnableToLoadException(format("No data was read from table %s", path_to_logph.c_str()));
    }
    try{
        msgpack::unpack(upd, sbuffer.c_str(), N);
        msgpack::object deserialized = upd.get();
        LogPHTable temp_single_phase_logph;
        deserialized.convert(&temp_single_phase_logph);
        temp_single_phase_logph.unpack_matrices();
        if (single_phase_logph.Nx != temp_single_phase_logph.Nx || single_phase_logph.Ny != temp_single_phase_logph.Ny)
        {
            throw ValueError(format("old [%dx%d] and new [%dx%d] dimensions don't agree",temp_single_phase_logph.Nx, temp_single_phase_logph.Ny, single_phase_logph.Nx, single_phase_logph.Ny));
        }
        else if (single_phase_logph.revision > temp_single_phase_logph.revision)
        {
            throw ValueError(format("loaded revision [%d] is older than current revision [%d]", temp_single_phase_logph.revision, single_phase_logph.revision));
        }
        single_phase_logph = temp_single_phase_logph; // Copy
        
        if (get_debug_level() > 0){std::cout << format("Loaded table: %s", path_to_logph.c_str()) << std::endl;}
    }
    catch(std::exception &){
        std::string err = format("Unable to deserialize %s", path_to_logph.c_str());
        if (get_debug_level() > 0){std::cout << err << std::endl;}
        throw UnableToLoadException(err);
    }
}
    
} /* namespace CoolProp */