
#if !defined(NO_TABULAR_BACKENDS)

#include "TabularBackends.h"
#include "CoolProp.h"
#include <sstream>
#include "time.h"
#include "miniz.h"

namespace CoolProp{
/**
 * @brief 
 * @param table
 * @param path_to_tables
 * @param filename
 */
template <typename T> void load_table(T &table, const std::string &path_to_tables, const std::string &filename){
    
    double tic = clock();
    std::string path_to_table = path_to_tables + "/" + filename;
    if (get_debug_level() > -1){std::cout << format("Loading table: %s", path_to_table.c_str()) << std::endl;}
    std::vector<char> raw;
    try{
         raw = get_binary_file_contents(path_to_table.c_str());
    }catch(...){
        std::string err = format("Unable to load file %s", path_to_table.c_str());
        if (get_debug_level() > 0){std::cout << err << std::endl;}
        throw UnableToLoadError(err);
    }
    std::vector<char> newBuffer(raw.size()*5);
    uLong newBufferSize = newBuffer.size();
    int code;
    do{
        code = uncompress((unsigned char *)(&(newBuffer[0])), &newBufferSize, 
                          (unsigned char *)(&(raw[0])), raw.size());
        if (code == Z_BUF_ERROR){ 
            // Output buffer is too small, make it bigger and try again
            newBuffer.resize(newBuffer.size()*5);
            newBufferSize = newBuffer.size();
        }
        else if (code != 0){ // Something else, a big problem
            std::string err = format("Unable to uncompress file %s with miniz code %d", path_to_table.c_str(), code);
            if (get_debug_level() > 0){std::cout << err << std::endl;}
            throw UnableToLoadError(err);
        }
    }while(code != 0);

    try{
        msgpack::unpacked msg;
        msgpack::unpack(&msg, &(newBuffer[0]), newBufferSize);
        msgpack::object deserialized = msg.get();
        
        // Call the class' deserialize function;  if it is an invalid table, it will cause an exception to be thrown
        table.deserialize(deserialized);
        double toc = clock();
        if (get_debug_level() > -1){std::cout << format("Loaded table: %s in %g sec.", path_to_table.c_str(), (toc-tic)/CLOCKS_PER_SEC) << std::endl;}
    }
    catch(std::exception &e){
        std::string err = format("Unable to deserialize %s; err: %s", path_to_table.c_str(), e.what());
        if (get_debug_level() > 0){std::cout << err << std::endl;}
        throw UnableToLoadError(err);
    }
}
template <typename T> void write_table(const T &table, const std::string &path_to_tables, const std::string &name)
{
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, table);
    std::string tabPath = std::string(path_to_tables + "/" + name + ".bin");
    std::string zPath = tabPath + ".z";
    std::vector<char> buffer(sbuf.size());
    uLong outSize = buffer.size();
    compress((unsigned char *)(&(buffer[0])), &outSize, 
             (unsigned char*)(sbuf.data()), sbuf.size());
    std::ofstream ofs2(zPath.c_str(), std::ofstream::binary);
    ofs2.write(&buffer[0], outSize);
    
    if (CoolProp::get_config_bool(SAVE_RAW_TABLES)){
        std::ofstream ofs(tabPath.c_str(), std::ofstream::binary);
        ofs.write(sbuf.data(), sbuf.size());
    }
}

} // namespace CoolProp

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
    AS->update(QT_INPUTS, 0, AS->Ttriple());
    CoolPropDbl p_triple = AS->p();
    // If a mixture, maximum pressure in the phase envelope can be greater than critical pressure
    CoolPropDbl p, pmin = p_triple*1.001, pmax = 0.9999999*AS->p_critical();
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
        // Transport properties - if no transport properties, just keep going
        try{
            viscL[i] = AS->viscosity(); condL[i] = AS->conductivity();
            logviscL[i] = log(viscL[i]);
        }
        catch(std::exception &e){
            if (debug){std::cout << " " << e.what() << std::endl;}
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
        // Transport properties - if no transport properties, just keep going
        try{
            viscV[i] = AS->viscosity(); condV[i] = AS->conductivity();
            logviscV[i] = log(viscV[i]);
        }
        catch(std::exception &e){
            if (debug){std::cout << " " << e.what() << std::endl;}
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
        std::cout << format(" Single-Phase Table (%s) \n", strjoin(AS->fluid_names(), "&").c_str());
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
            if (is_in_closed_range(0.0, 1.0, AS->Q())){ 
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
			umolar[i][j] = AS->umolar();
            
            // -------------------------
            //   Transport properties
            // -------------------------
            try{
                visc[i][j] = AS->viscosity();
                cond[i][j] = AS->conductivity();
            }
            catch(std::exception &){
                // Failures will remain as holes in table
            }
            
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
			dumolardx[i][j] = AS->first_partial_deriv(iUmolar, xkey, ykey);
            dumolardy[i][j] = AS->first_partial_deriv(iUmolar, ykey, xkey);
            
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
			d2umolardx2[i][j] = AS->second_partial_deriv(iUmolar, xkey, ykey, xkey, ykey);
            d2umolardxdy[i][j] = AS->second_partial_deriv(iUmolar, xkey, ykey, ykey, xkey);
            d2umolardy2[i][j] = AS->second_partial_deriv(iUmolar, ykey, xkey, ykey, xkey);
        }
    }
}
std::string CoolProp::TabularBackend::path_to_tables(void){
    std::vector<std::string> fluids = AS->fluid_names();
    std::vector<CoolPropDbl> fractions = AS->get_mole_fractions();
    std::vector<std::string> components;
    for (std::size_t i = 0; i < fluids.size(); ++i){
        components.push_back(format("%s[%0.10Lf]", fluids[i].c_str(), fractions[i]));
    }
    return get_home_dir() + "/.CoolProp/Tables/" + AS->backend_name() + "(" + strjoin(components, "&") + ")";
}

void CoolProp::TabularBackend::write_tables(){
    std::string path_to_tables = this->path_to_tables();
    make_dirs(path_to_tables);
    write_table(single_phase_logph, path_to_tables, "single_phase_logph");
    write_table(single_phase_logpT, path_to_tables, "single_phase_logpT");
    write_table(pure_saturation, path_to_tables, "pure_saturation");
    write_table(phase_envelope, path_to_tables, "phase_envelope");
}
void CoolProp::TabularBackend::load_tables(){
    std::string path_to_tables = this->path_to_tables();
    single_phase_logph.AS = this->AS;
    single_phase_logpT.AS = this->AS;
    single_phase_logph.AS->set_mole_fractions(AS->get_mole_fractions());
    single_phase_logpT.AS->set_mole_fractions(AS->get_mole_fractions());
    single_phase_logph.set_limits();
    single_phase_logpT.set_limits();
    load_table(single_phase_logph, path_to_tables, "single_phase_logph.bin.z");
    load_table(single_phase_logpT, path_to_tables, "single_phase_logpT.bin.z");
    load_table(pure_saturation, path_to_tables, "pure_saturation.bin.z");
    load_table(phase_envelope, path_to_tables, "phase_envelope.bin.z");
}


#endif // !defined(NO_TABULAR_BACKENDS)