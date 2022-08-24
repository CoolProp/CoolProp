
#if !defined(NO_TABULAR_BACKENDS)

#    include "TabularBackends.h"
#    include "CoolProp.h"
#    include <sstream>
#    include "time.h"
#    include "miniz.h"
#    include <fstream>

/// The inverse of the A matrix for the bicubic interpolation (http://en.wikipedia.org/wiki/Bicubic_interpolation)
/// NOTE: The matrix is transposed below
static const double Ainv_data[16 * 16] = {
  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  -3, 3,  0,  0,  -2,
  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2, 0,  0,  1,  1,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  -3, 3,  0,  0,  -2, -1, 0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  -2, 0, 0,  1,  1,  0,  0,  -3, 0,  3,  0,  0,  0, 0,  0,  -2, 0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,
  -3, 0,  3,  0,  0,  0,  0,  0,  -2, 0,  -1, 0, 9,  -9, -9, 9,  6,  3,  -6, -3, 6,  -6, 3, -3, 4,  2,  2, 1,  -6, 6,  6,  -6, -3, -3, 3,  3,  -4,
  4,  -2, 2,  -2, -2, -1, -1, 2,  0,  -2, 0,  0, 0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0, 0,  0,  0,  0, 2,  0,  -2, 0,  0,  0,  0,  0,  1,  0,
  1,  0,  -6, 6,  6,  -6, -4, -2, 4,  2,  -3, 3, -3, 3,  -2, -1, -2, -1, 4,  -4, -4, 4,  2, 2,  -2, -2, 2, -2, 2,  -2, 1,  1,  1,  1};
static Eigen::Matrix<double, 16, 16> Ainv(Ainv_data);

static CoolProp::TabularDataLibrary library;

namespace CoolProp {

/**
 * @brief
 * @param table
 * @param path_to_tables
 * @param filename
 */
template <typename T>
void load_table(T& table, const std::string& path_to_tables, const std::string& filename) {

    double tic = clock();
    std::string path_to_table = path_to_tables + "/" + filename;
    if (get_debug_level() > 0) {
        std::cout << format("Loading table: %s", path_to_table.c_str()) << std::endl;
    }
    std::vector<char> raw;
    try {
        raw = get_binary_file_contents(path_to_table.c_str());
    } catch (...) {
        std::string err = format("Unable to load file %s", path_to_table.c_str());
        if (get_debug_level() > 0) {
            std::cout << "err:" << err << std::endl;
        }
        throw UnableToLoadError(err);
    }
    std::vector<unsigned char> newBuffer(raw.size() * 5);
    uLong newBufferSize = static_cast<uLong>(newBuffer.size());
    mz_ulong rawBufferSize = static_cast<mz_ulong>(raw.size());
    int code;
    do {
        code = uncompress((unsigned char*)(&(newBuffer[0])), &newBufferSize, (unsigned char*)(&(raw[0])), rawBufferSize);
        if (code == Z_BUF_ERROR) {
            // Output buffer is too small, make it bigger and try again
            newBuffer.resize(newBuffer.size() * 5);
            newBufferSize = static_cast<uLong>(newBuffer.size());
        } else if (code != 0) {  // Something else, a big problem
            std::string err = format("Unable to uncompress file %s with miniz code %d", path_to_table.c_str(), code);
            if (get_debug_level() > 0) {
                std::cout << "uncompress err:" << err << std::endl;
            }
            throw UnableToLoadError(err);
        }
    } while (code != 0);
    // Copy the buffer from unsigned char to char (yuck)
    std::vector<char> charbuffer(newBuffer.begin(), newBuffer.begin() + newBufferSize);
    try {
        msgpack::unpacked msg;
        msgpack::unpack(msg, &(charbuffer[0]), charbuffer.size());
        msgpack::object deserialized = msg.get();

        // Call the class' deserialize function;  if it is an invalid table, it will cause an exception to be thrown
        table.deserialize(deserialized);
        double toc = clock();
        if (get_debug_level() > 0) {
            std::cout << format("Loaded table: %s in %g sec.", path_to_table.c_str(), (toc - tic) / CLOCKS_PER_SEC) << std::endl;
        }
    } catch (std::exception& e) {
        std::string err = format("Unable to msgpack deserialize %s; err: %s", path_to_table.c_str(), e.what());
        if (get_debug_level() > 0) {
            std::cout << "err: " << err << std::endl;
        }
        throw UnableToLoadError(err);
    }
}
template <typename T>
void write_table(const T& table, const std::string& path_to_tables, const std::string& name) {
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, table);
    std::string tabPath = std::string(path_to_tables + "/" + name + ".bin");
    std::string zPath = tabPath + ".z";
    std::vector<char> buffer(sbuf.size());
    uLong outSize = static_cast<uLong>(buffer.size());
    compress((unsigned char*)(&(buffer[0])), &outSize, (unsigned char*)(sbuf.data()), static_cast<mz_ulong>(sbuf.size()));
    std::ofstream ofs2(zPath.c_str(), std::ofstream::binary);
    ofs2.write(&buffer[0], outSize);
    ofs2.close();

    if (CoolProp::get_config_bool(SAVE_RAW_TABLES)) {
        std::ofstream ofs(tabPath.c_str(), std::ofstream::binary);
        ofs.write(sbuf.data(), sbuf.size());
    }
}

}  // namespace CoolProp

void CoolProp::PureFluidSaturationTableData::build(shared_ptr<CoolProp::AbstractState>& AS) {
    const bool debug = get_debug_level() > 5 || false;
    if (debug) {
        std::cout << format("***********************************************\n");
        std::cout << format(" Saturation Table (%s) \n", AS->name().c_str());
        std::cout << format("***********************************************\n");
    }
    resize(N);
    // ------------------------
    // Actually build the table
    // ------------------------
    CoolPropDbl Tmin = std::max(AS->Ttriple(), AS->Tmin());
    AS->update(QT_INPUTS, 0, Tmin);
    CoolPropDbl p_triple = AS->p();
    CoolPropDbl p, pmin = p_triple, pmax = 0.9999 * AS->p_critical();
    for (std::size_t i = 0; i < N - 1; ++i) {
        if (i == 0) {
            CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, true);
        }
        // Log spaced
        p = exp(log(pmin) + (log(pmax) - log(pmin)) / (N - 1) * i);

        // Saturated liquid
        try {
            AS->update(PQ_INPUTS, p, 0);
            pL[i] = p;
            TL[i] = AS->T();
            rhomolarL[i] = AS->rhomolar();
            hmolarL[i] = AS->hmolar();
            smolarL[i] = AS->smolar();
            umolarL[i] = AS->umolar();
            logpL[i] = log(p);
            logrhomolarL[i] = log(rhomolarL[i]);
            cpmolarL[i] = AS->cpmolar();
            cvmolarL[i] = AS->cvmolar();
            speed_soundL[i] = AS->speed_sound();
        } catch (std::exception& e) {
            // That failed for some reason, go to the next pair
            if (debug) {
                std::cout << " " << e.what() << std::endl;
            }
            continue;
        }
        // Transport properties - if no transport properties, just keep going
        try {
            viscL[i] = AS->viscosity();
            condL[i] = AS->conductivity();
            logviscL[i] = log(viscL[i]);
        } catch (std::exception& e) {
            if (debug) {
                std::cout << " " << e.what() << std::endl;
            }
        }
        // Saturated vapor
        try {
            AS->update(PQ_INPUTS, p, 1);
            pV[i] = p;
            TV[i] = AS->T();
            rhomolarV[i] = AS->rhomolar();
            hmolarV[i] = AS->hmolar();
            smolarV[i] = AS->smolar();
            umolarV[i] = AS->umolar();
            logpV[i] = log(p);
            logrhomolarV[i] = log(rhomolarV[i]);
            cpmolarV[i] = AS->cpmolar();
            cvmolarV[i] = AS->cvmolar();
            speed_soundV[i] = AS->speed_sound();
        } catch (std::exception& e) {
            // That failed for some reason, go to the next pair
            if (debug) {
                std::cout << " " << e.what() << std::endl;
            }
            continue;
        }
        // Transport properties - if no transport properties, just keep going
        try {
            viscV[i] = AS->viscosity();
            condV[i] = AS->conductivity();
            logviscV[i] = log(viscV[i]);
        } catch (std::exception& e) {
            if (debug) {
                std::cout << " " << e.what() << std::endl;
            }
        }
        if (i == 0) {
            CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, false);
        }
    }
    // Last point is at the critical point
    AS->update(PQ_INPUTS, AS->p_critical(), 1);
    std::size_t i = N - 1;

    pV[i] = AS->p();
    TV[i] = AS->T();
    rhomolarV[i] = AS->rhomolar();
    hmolarV[i] = AS->hmolar();
    smolarV[i] = AS->smolar();
    umolarV[i] = AS->umolar();

    pL[i] = AS->p();
    TL[i] = AS->T();
    rhomolarL[i] = AS->rhomolar();
    hmolarL[i] = AS->hmolar();
    smolarL[i] = AS->smolar();
    umolarL[i] = AS->umolar();

    logpV[i] = log(AS->p());
    logrhomolarV[i] = log(rhomolarV[i]);

    logpL[i] = log(AS->p());
    logrhomolarL[i] = log(rhomolarL[i]);
}

void CoolProp::SinglePhaseGriddedTableData::build(shared_ptr<CoolProp::AbstractState>& AS) {
    CoolPropDbl x, y;
    const bool debug = get_debug_level() > 5 || false;

    resize(Nx, Ny);

    if (debug) {
        std::cout << format("***********************************************\n");
        std::cout << format(" Single-Phase Table (%s) \n", strjoin(AS->fluid_names(), "&").c_str());
        std::cout << format("***********************************************\n");
    }
    // ------------------------
    // Actually build the table
    // ------------------------
    for (std::size_t i = 0; i < Nx; ++i) {
        // Calculate the x value
        if (logx) {
            // Log spaced
            x = exp(log(xmin) + (log(xmax) - log(xmin)) / (Nx - 1) * i);
        } else {
            // Linearly spaced
            x = xmin + (xmax - xmin) / (Nx - 1) * i;
        }
        xvec[i] = x;
        for (std::size_t j = 0; j < Ny; ++j) {
            // Calculate the x value
            if (logy) {
                // Log spaced
                y = exp(log(ymin) + (log(ymax / ymin)) / (Ny - 1) * j);
            } else {
                // Linearly spaced
                y = ymin + (ymax - ymin) / (Ny - 1) * j;
            }
            yvec[j] = y;

            if (debug) {
                std::cout << "x: " << x << " y: " << y << std::endl;
            }

            // Generate the input pair
            CoolPropDbl v1, v2;
            input_pairs input_pair = generate_update_pair(xkey, x, ykey, y, v1, v2);

            // --------------------
            //   Update the state
            // --------------------
            try {
                AS->update(input_pair, v1, v2);
                if (!ValidNumber(AS->rhomolar())) {
                    throw ValueError("rhomolar is invalid");
                }
            } catch (std::exception& e) {
                // That failed for some reason, go to the next pair
                if (debug) {
                    std::cout << " " << e.what() << std::endl;
                }
                continue;
            }

            // Skip two-phase states - they will remain as _HUGE holes in the table
            if (is_in_closed_range(0.0, 1.0, AS->Q())) {
                if (debug) {
                    std::cout << " 2Phase" << std::endl;
                }
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
            try {
                visc[i][j] = AS->viscosity();
                cond[i][j] = AS->conductivity();
            } catch (std::exception&) {
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
std::string CoolProp::TabularBackend::path_to_tables(void) {
    std::vector<std::string> fluids = AS->fluid_names();
    std::vector<CoolPropDbl> fractions = AS->get_mole_fractions();
    std::vector<std::string> components;
    for (std::size_t i = 0; i < fluids.size(); ++i) {
        components.push_back(format("%s[%0.10Lf]", fluids[i].c_str(), fractions[i]));
    }
    std::string table_directory = get_home_dir() + "/.CoolProp/Tables/";
    std::string alt_table_directory = get_config_string(ALTERNATIVE_TABLES_DIRECTORY);
    if (!alt_table_directory.empty()) {
        table_directory = alt_table_directory;
    }
    return table_directory + AS->backend_name() + "(" + strjoin(components, "&") + ")";
}

void CoolProp::TabularBackend::write_tables() {
    std::string path_to_tables = this->path_to_tables();
    make_dirs(path_to_tables);
    bool loaded = false;
    dataset = library.get_set_of_tables(this->AS, loaded);
    PackablePhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
    SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;
    write_table(single_phase_logph, path_to_tables, "single_phase_logph");
    write_table(single_phase_logpT, path_to_tables, "single_phase_logpT");
    write_table(pure_saturation, path_to_tables, "pure_saturation");
    write_table(phase_envelope, path_to_tables, "phase_envelope");
}
void CoolProp::TabularBackend::load_tables() {
    bool loaded = false;
    dataset = library.get_set_of_tables(this->AS, loaded);
    if (loaded == false) {
        throw UnableToLoadError("Could not load tables");
    }
    if (get_debug_level() > 0) {
        std::cout << "Tables loaded" << std::endl;
    }
}

CoolPropDbl CoolProp::TabularBackend::calc_saturated_vapor_keyed_output(parameters key) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    double factor = 1.0;
    mass_to_molar(key, factor, molar_mass());
    if (is_mixture) {
        return phase_envelope_sat(phase_envelope, key, iP, _p) * factor;
    } else {
        return pure_saturation.evaluate(key, _p, 1, cached_saturation_iL, cached_saturation_iV) * factor;
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_saturated_liquid_keyed_output(parameters key) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    double factor = 1.0;
    mass_to_molar(key, factor, molar_mass());
    if (is_mixture) {
        return phase_envelope_sat(phase_envelope, key, iP, _p) * factor;
    } else {
        return pure_saturation.evaluate(key, _p, 0, cached_saturation_iL, cached_saturation_iV) * factor;
    }
};

CoolPropDbl CoolProp::TabularBackend::calc_p(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    if (using_single_phase_table) {
        return _p;
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iP, iT, _T);
        } else {
            return _p;
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_T(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar(iT, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return _T;
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iT, iP, _p);
        } else {
            if (ValidNumber(_T)) {
                return _T;
            } else {
                return pure_saturation.evaluate(iT, _p, _Q, cached_saturation_iL, cached_saturation_iV);
            }
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_rhomolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar(iDmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT(iDmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iDmolar, iP, _p);
        } else {
            return pure_saturation.evaluate(iDmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_hmolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return _hmolar;
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT(iHmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iHmolar, iP, _p);
        } else {
            return pure_saturation.evaluate(iHmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_smolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar(iSmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT(iSmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iSmolar, iP, _p);
        } else {
            return pure_saturation.evaluate(iSmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_umolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar(iUmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT(iUmolar, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            // h = u + pv -> u = h - pv
            CoolPropDbl hmolar = phase_envelope_sat(phase_envelope, iHmolar, iP, _p);
            CoolPropDbl rhomolar = phase_envelope_sat(phase_envelope, iDmolar, iP, _p);
            return hmolar - _p / rhomolar;
        } else {
            return pure_saturation.evaluate(iUmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_cpmolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        return calc_first_partial_deriv(iHmolar, iT, iP);
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iCpmolar, iP, _p);
        } else {
            return pure_saturation.evaluate(iCpmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_cvmolar(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        return calc_first_partial_deriv(iUmolar, iT, iDmolar);
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iCvmolar, iP, _p);
        } else {
            return pure_saturation.evaluate(iCvmolar, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}

CoolPropDbl CoolProp::TabularBackend::calc_viscosity(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar_transport(iviscosity, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT_transport(iviscosity, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iviscosity, iP, _p);
        } else {
            return pure_saturation.evaluate(iviscosity, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_conductivity(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        switch (selected_table) {
            case SELECTED_PH_TABLE:
                return evaluate_single_phase_phmolar_transport(iconductivity, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_PT_TABLE:
                return evaluate_single_phase_pT_transport(iconductivity, cached_single_phase_i, cached_single_phase_j);
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        return _HUGE;  // not needed, will never be hit, just to make compiler happy
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, iconductivity, iP, _p);
        } else {
            return pure_saturation.evaluate(iconductivity, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_speed_sound(void) {
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (using_single_phase_table) {
        return sqrt(1 / molar_mass() * cpmolar() / cvmolar() * first_partial_deriv(iP, iDmolar, iT));
    } else {
        if (is_mixture) {
            return phase_envelope_sat(phase_envelope, ispeed_sound, iP, _p);
        } else {
            return pure_saturation.evaluate(ispeed_sound, _p, _Q, cached_saturation_iL, cached_saturation_iV);
        }
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
    if (using_single_phase_table) {
        CoolPropDbl dOf_dx, dOf_dy, dWrt_dx, dWrt_dy, dConstant_dx, dConstant_dy;

        // If a mass-based parameter is provided, get a conversion factor and change the key to the molar-based key
        double Of_conversion_factor = 1.0, Wrt_conversion_factor = 1.0, Constant_conversion_factor = 1.0, MM = AS->molar_mass();
        mass_to_molar(Of, Of_conversion_factor, MM);
        mass_to_molar(Wrt, Wrt_conversion_factor, MM);
        mass_to_molar(Constant, Constant_conversion_factor, MM);

        switch (selected_table) {
            case SELECTED_PH_TABLE: {
                dOf_dx = evaluate_single_phase_phmolar_derivative(Of, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dOf_dy = evaluate_single_phase_phmolar_derivative(Of, cached_single_phase_i, cached_single_phase_j, 0, 1);
                dWrt_dx = evaluate_single_phase_phmolar_derivative(Wrt, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dWrt_dy = evaluate_single_phase_phmolar_derivative(Wrt, cached_single_phase_i, cached_single_phase_j, 0, 1);
                dConstant_dx = evaluate_single_phase_phmolar_derivative(Constant, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dConstant_dy = evaluate_single_phase_phmolar_derivative(Constant, cached_single_phase_i, cached_single_phase_j, 0, 1);
                break;
            }
            case SELECTED_PT_TABLE: {
                dOf_dx = evaluate_single_phase_pT_derivative(Of, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dOf_dy = evaluate_single_phase_pT_derivative(Of, cached_single_phase_i, cached_single_phase_j, 0, 1);
                dWrt_dx = evaluate_single_phase_pT_derivative(Wrt, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dWrt_dy = evaluate_single_phase_pT_derivative(Wrt, cached_single_phase_i, cached_single_phase_j, 0, 1);
                dConstant_dx = evaluate_single_phase_pT_derivative(Constant, cached_single_phase_i, cached_single_phase_j, 1, 0);
                dConstant_dy = evaluate_single_phase_pT_derivative(Constant, cached_single_phase_i, cached_single_phase_j, 0, 1);
                break;
            }
            case SELECTED_NO_TABLE:
                throw ValueError("table not selected");
        }
        double val = (dOf_dx * dConstant_dy - dOf_dy * dConstant_dx) / (dWrt_dx * dConstant_dy - dWrt_dy * dConstant_dx);
        return val * Of_conversion_factor / Wrt_conversion_factor;
    } else {
        throw ValueError(format("Inputs [rho: %g mol/m3, T: %g K, p: %g Pa] are two-phase; cannot use single-phase derivatives", _rhomolar, _T, _p));
    }
};

CoolPropDbl CoolProp::TabularBackend::calc_first_saturation_deriv(parameters Of1, parameters Wrt1) {
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (AS->get_mole_fractions().size() > 1) {
        throw ValueError("calc_first_saturation_deriv not available for mixtures");
    }
    if (std::abs(_Q) < 1e-6) {
        return pure_saturation.first_saturation_deriv(Of1, Wrt1, 0, keyed_output(Wrt1), cached_saturation_iL);
    } else if (std::abs(_Q - 1) < 1e-6) {
        return pure_saturation.first_saturation_deriv(Of1, Wrt1, 1, keyed_output(Wrt1), cached_saturation_iV);
    } else {
        throw ValueError(format("Quality [%Lg] must be either 0 or 1 to within 1 ppm", _Q));
    }
}
CoolPropDbl CoolProp::TabularBackend::calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        CoolPropDbl rhoL = pure_saturation.evaluate(iDmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl rhoV = pure_saturation.evaluate(iDmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl hL = pure_saturation.evaluate(iHmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl hV = pure_saturation.evaluate(iHmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
        return -POW2(rhomolar()) * (1 / rhoV - 1 / rhoL) / (hV - hL);
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return first_two_phase_deriv(iDmolar, iHmolar, iP) * POW2(molar_mass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        // v = 1/rho; drhodv = -rho^2; dvdrho = -1/rho^2
        CoolPropDbl rhoL = pure_saturation.evaluate(iDmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl rhoV = pure_saturation.evaluate(iDmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl hL = pure_saturation.evaluate(iHmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl hV = pure_saturation.evaluate(iHmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
        CoolPropDbl dvdrhoL = -1 / POW2(rhoL);
        CoolPropDbl dvdrhoV = -1 / POW2(rhoV);
        CoolPropDbl dvL_dp = dvdrhoL * pure_saturation.first_saturation_deriv(iDmolar, iP, 0, _p, cached_saturation_iL);
        CoolPropDbl dvV_dp = dvdrhoV * pure_saturation.first_saturation_deriv(iDmolar, iP, 1, _p, cached_saturation_iV);
        CoolPropDbl dhL_dp = pure_saturation.first_saturation_deriv(iHmolar, iP, 0, _p, cached_saturation_iL);
        CoolPropDbl dhV_dp = pure_saturation.first_saturation_deriv(iHmolar, iP, 1, _p, cached_saturation_iV);
        CoolPropDbl dxdp_h = (Q() * dhV_dp + (1 - Q()) * dhL_dp) / (hL - hV);
        CoolPropDbl dvdp_h = dvL_dp + dxdp_h * (1 / rhoV - 1 / rhoL) + Q() * (dvV_dp - dvL_dp);
        return -POW2(rhomolar()) * dvdp_h;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        return first_two_phase_deriv(iDmolar, iP, iHmolar) * molar_mass();
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
    }
}

CoolPropDbl CoolProp::TabularBackend::calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) {
    // Note: If you need all three values (drho_dh__p, drho_dp__h and rho_spline),
    // you should calculate drho_dp__h first to avoid duplicate calculations.

    bool drho_dh__p = false;
    bool drho_dp__h = false;
    bool rho_spline = false;

    if (Of == iDmolar && Wrt == iHmolar && Constant == iP) {
        drho_dh__p = true;
        if (_drho_spline_dh__constp) return _drho_spline_dh__constp;
    } else if (Of == iDmass && Wrt == iHmass && Constant == iP) {
        return first_two_phase_deriv_splined(iDmolar, iHmolar, iP, x_end) * POW2(molar_mass());
    } else if (Of == iDmolar && Wrt == iP && Constant == iHmolar) {
        drho_dp__h = true;
        if (_drho_spline_dp__consth) return _drho_spline_dp__consth;
    } else if (Of == iDmass && Wrt == iP && Constant == iHmass) {
        return first_two_phase_deriv_splined(iDmolar, iP, iHmolar, x_end) * molar_mass();
    }
    // Add the special case for the splined density
    else if (Of == iDmolar && Wrt == iDmolar && Constant == iDmolar) {
        rho_spline = true;
        if (_rho_spline) return _rho_spline;
    } else if (Of == iDmass && Wrt == iDmass && Constant == iDmass) {
        return first_two_phase_deriv_splined(iDmolar, iDmolar, iDmolar, x_end) * molar_mass();
    } else {
        throw ValueError("These inputs are not supported to calc_first_two_phase_deriv");
    }

    if (_Q > x_end) {
        throw ValueError(format("Q [%g] is greater than x_end [%Lg]", _Q, x_end).c_str());
    }
    if (_phase != iphase_twophase) {
        throw ValueError(format("state is not two-phase"));
    }

    // TODO: replace AS with a cloned instance of TTSE or BICUBIC, add "clone()" as mandatory function in base class
    //shared_ptr<TabularBackend>
    //	Liq(this->clone())),
    //	End(this->clone()));
    //
    //Liq->specify_phase(iphase_liquid);
    //Liq->_Q = -1;
    //Liq->update_DmolarT_direct(SatL->rhomolar(), SatL->T());
    //End->update(QT_INPUTS, x_end, SatL->T());

    // Start with quantities needed for all calculations
    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    CoolPropDbl hL = pure_saturation.evaluate(iHmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
    CoolPropDbl hV = pure_saturation.evaluate(iHmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
    CoolPropDbl hE = pure_saturation.evaluate(iHmolar, _p, x_end, cached_saturation_iL, cached_saturation_iV);

    CoolPropDbl dL = pure_saturation.evaluate(iDmolar, _p, 0, cached_saturation_iL, cached_saturation_iV);
    CoolPropDbl dV = pure_saturation.evaluate(iDmolar, _p, 1, cached_saturation_iL, cached_saturation_iV);
    CoolPropDbl dE = pure_saturation.evaluate(iDmolar, _p, x_end, cached_saturation_iL, cached_saturation_iV);

    CoolPropDbl Delta = Q() * (hV - hL);
    CoolPropDbl Delta_end = hE - hL;

    CoolPropDbl TL = pure_saturation.evaluate(iT, _p, 0, cached_saturation_iL, cached_saturation_iV);
    CoolPropDbl TE = pure_saturation.evaluate(iT, _p, x_end, cached_saturation_iL, cached_saturation_iV);

    // TODO: We cheat here and fake access to a derived class.
    // Liquid state
    AS->specify_phase(iphase_liquid);
    AS->update(DmolarT_INPUTS, dL, TL);
    CoolPropDbl drho_dh_liq__constp = AS->first_partial_deriv(iDmolar, iHmolar, iP);
    CoolPropDbl d2rhodhdp_liq = AS->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);
    // End of spline
    AS->specify_phase(iphase_twophase);
    AS->update(DmolarT_INPUTS, dE, TE);
    CoolPropDbl drho_dh_end__constp = AS->first_partial_deriv(iDmolar, iHmolar, iP);
    CoolPropDbl d2rhodhdp_end = AS->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);

    // Spline coordinates a, b, c, d
    CoolPropDbl Abracket = (2 * dL - 2 * dE + Delta_end * (drho_dh_liq__constp + drho_dh_end__constp));
    CoolPropDbl a = 1 / POW3(Delta_end) * Abracket;
    CoolPropDbl b = 3 / POW2(Delta_end) * (-dL + dE) - 1 / Delta_end * (drho_dh_end__constp + 2 * drho_dh_liq__constp);
    CoolPropDbl c = drho_dh_liq__constp;
    CoolPropDbl d = dL;

    // Either the spline value or drho/dh|p can be directly evaluated now
    _rho_spline = a * POW3(Delta) + b * POW2(Delta) + c * Delta + d;
    _drho_spline_dh__constp = 3 * a * POW2(Delta) + 2 * b * Delta + c;
    if (rho_spline) return _rho_spline;
    if (drho_dh__p) return _drho_spline_dh__constp;

    // It's drho/dp|h
    // ... calculate some more things

    // Derivatives *along* the saturation curve using the special internal method
    CoolPropDbl dhL_dp_sat = pure_saturation.first_saturation_deriv(iHmolar, iP, 0, _p, cached_saturation_iL);
    CoolPropDbl dhV_dp_sat = pure_saturation.first_saturation_deriv(iHmolar, iP, 1, _p, cached_saturation_iV);
    CoolPropDbl drhoL_dp_sat = pure_saturation.first_saturation_deriv(iDmolar, iP, 0, _p, cached_saturation_iL);
    CoolPropDbl drhoV_dp_sat = pure_saturation.first_saturation_deriv(iDmolar, iP, 1, _p, cached_saturation_iV);
    //CoolPropDbl rhoV = SatV->keyed_output(rho_key);
    //CoolPropDbl rhoL = SatL->keyed_output(rho_key);
    CoolPropDbl drho_dp_end = POW2(dE) * (x_end / POW2(dV) * drhoV_dp_sat + (1 - x_end) / POW2(dL) * drhoL_dp_sat);

    // Faking single-phase
    //CoolPropDbl drho_dp__consth_liq = Liq->first_partial_deriv(rho_key, p_key, h_key);
    //CoolPropDbl d2rhodhdp_liq = Liq->second_partial_deriv(rho_key, h_key, p_key, p_key, h_key); // ?

    // Derivatives at the end point
    // CoolPropDbl drho_dp__consth_end = End->calc_first_two_phase_deriv(rho_key, p_key, h_key);
    //CoolPropDbl d2rhodhdp_end = End->calc_second_two_phase_deriv(rho_key, h_key, p_key, p_key, h_key);

    // Reminder:
    // Delta = Q()*(hV-hL) = h-hL
    // Delta_end = x_end*(hV-hL);
    CoolPropDbl d_Delta_dp__consth = -dhL_dp_sat;
    CoolPropDbl d_Delta_end_dp__consth = x_end * (dhV_dp_sat - dhL_dp_sat);

    // First pressure derivative at constant h of the coefficients a,b,c,d
    // CoolPropDbl Abracket = (2*rho_liq - 2*rho_end + Delta_end * (drho_dh_liq__constp + drho_dh_end));
    CoolPropDbl d_Abracket_dp_consth = (2 * drhoL_dp_sat - 2 * drho_dp_end + Delta_end * (d2rhodhdp_liq + d2rhodhdp_end)
                                        + d_Delta_end_dp__consth * (drho_dh_liq__constp + drho_dh_end__constp));
    CoolPropDbl da_dp = 1 / POW3(Delta_end) * d_Abracket_dp_consth + Abracket * (-3 / POW4(Delta_end) * d_Delta_end_dp__consth);
    CoolPropDbl db_dp = -6 / POW3(Delta_end) * d_Delta_end_dp__consth * (dE - dL) + 3 / POW2(Delta_end) * (drho_dp_end - drhoL_dp_sat)
                        + (1 / POW2(Delta_end) * d_Delta_end_dp__consth) * (drho_dh_end__constp + 2 * drho_dh_liq__constp)
                        - (1 / Delta_end) * (d2rhodhdp_end + 2 * d2rhodhdp_liq);
    CoolPropDbl dc_dp = d2rhodhdp_liq;
    CoolPropDbl dd_dp = drhoL_dp_sat;

    _drho_spline_dp__consth =
      (3 * a * POW2(Delta) + 2 * b * Delta + c) * d_Delta_dp__consth + POW3(Delta) * da_dp + POW2(Delta) * db_dp + Delta * dc_dp + dd_dp;
    if (drho_dp__h) return _drho_spline_dp__consth;

    throw ValueError("Something went wrong in TabularBackend::calc_first_two_phase_deriv_splined");
    return _HUGE;
}

void CoolProp::TabularBackend::update(CoolProp::input_pairs input_pair, double val1, double val2) {

    if (get_debug_level() > 0) {
        std::cout << format("update(%s,%g,%g)\n", get_input_pair_short_desc(input_pair).c_str(), val1, val2);
    }

    // Clear cached variables
    clear();

    // Convert to mass-based units if necessary
    CoolPropDbl ld_value1 = val1, ld_value2 = val2;
    mass_to_molar_inputs(input_pair, ld_value1, ld_value2);
    val1 = ld_value1;
    val2 = ld_value2;

    // Check the tables, build if necessary
    check_tables();

    // Flush the cached indices (set to large number)
    cached_single_phase_i = std::numeric_limits<std::size_t>::max();
    cached_single_phase_j = std::numeric_limits<std::size_t>::max();
    cached_saturation_iL = std::numeric_limits<std::size_t>::max();
    cached_saturation_iV = std::numeric_limits<std::size_t>::max();

    // To start, set quality to value that is impossible
    _Q = -1000;

    PureFluidSaturationTableData& pure_saturation = dataset->pure_saturation;
    PhaseEnvelopeData& phase_envelope = dataset->phase_envelope;
    SinglePhaseGriddedTableData& single_phase_logph = dataset->single_phase_logph;
    SinglePhaseGriddedTableData& single_phase_logpT = dataset->single_phase_logpT;

    switch (input_pair) {
        case HmolarP_INPUTS: {
            _hmolar = val1;
            _p = val2;
            if (!single_phase_logph.native_inputs_are_in_range(_hmolar, _p)) {
                // Use the AbstractState instance
                using_single_phase_table = false;
                if (get_debug_level() > 5) {
                    std::cout << "inputs are not in range";
                }
                throw ValueError(format("inputs are not in range, hmolar=%Lg, p=%Lg", static_cast<CoolPropDbl>(_hmolar), _p));
            } else {
                using_single_phase_table = true;  // Use the table !
                std::size_t iL = std::numeric_limits<std::size_t>::max(), iV = std::numeric_limits<std::size_t>::max(), iclosest = 0;
                CoolPropDbl hL = 0, hV = 0;
                SimpleState closest_state;
                bool is_two_phase = false;
                // If phase is imposed, use it, but only if it's single phase:
                //   - Imposed two phase still needs to determine saturation limits by calling pure_saturation.is_inside().
                //   - There's no speed increase to be gained by imposing two phase.
                if ((imposed_phase_index == iphase_not_imposed) || (imposed_phase_index == iphase_twophase)) {
                    if (is_mixture) {
                        is_two_phase = PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, iHmolar, _hmolar, iclosest, closest_state);
                    } else {
                        is_two_phase = pure_saturation.is_inside(iP, _p, iHmolar, _hmolar, iL, iV, hL, hV);
                    }
                }
                // Phase determined or imposed, now interpolate results
                if (is_two_phase) {
                    using_single_phase_table = false;
                    _Q = (static_cast<double>(_hmolar) - hL) / (hV - hL);
                    if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))) {
                        throw ValueError(
                          format("vapor quality is not in (0,1) for hmolar: %g p: %g, hL: %g hV: %g ", static_cast<double>(_hmolar), _p, hL, hV));
                    } else {
                        cached_saturation_iL = iL;
                        cached_saturation_iV = iV;
                        _phase = iphase_twophase;
                    }
                } else {
                    selected_table = SELECTED_PH_TABLE;
                    // Find and cache the indices i, j
                    find_native_nearest_good_indices(single_phase_logph, dataset->coeffs_ph, _hmolar, _p, cached_single_phase_i,
                                                     cached_single_phase_j);
                    // Recalculate the phase
                    recalculate_singlephase_phase();
                }
            }
            break;
        }
        case PT_INPUTS: {
            _p = val1;
            _T = val2;
            if (!single_phase_logpT.native_inputs_are_in_range(_T, _p)) {
                // Use the AbstractState instance
                using_single_phase_table = false;
                if (get_debug_level() > 5) {
                    std::cout << "inputs are not in range";
                }
                throw ValueError(format("inputs are not in range, p=%g Pa, T=%g K", _p, _T));
            } else {
                using_single_phase_table = true;  // Use the table !
                std::size_t iL = 0, iV = 0, iclosest = 0;
                CoolPropDbl TL = 0, TV = 0;
                SimpleState closest_state;
                bool is_two_phase = false;
                // Phase is imposed, use it
                if (imposed_phase_index != iphase_not_imposed) {
                    is_two_phase = (imposed_phase_index == iphase_twophase);
                } else {
                    if (is_mixture) {
                        is_two_phase = PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, iT, _T, iclosest, closest_state);
                    } else {
                        is_two_phase = pure_saturation.is_inside(iP, _p, iT, _T, iL, iV, TL, TV);
                    }
                }
                if (is_two_phase) {
                    using_single_phase_table = false;
                    throw ValueError(format("P,T with TTSE cannot be two-phase for now"));
                } else {
                    selected_table = SELECTED_PT_TABLE;
                    // Find and cache the indices i, j
                    find_native_nearest_good_indices(single_phase_logpT, dataset->coeffs_pT, _T, _p, cached_single_phase_i, cached_single_phase_j);

                    if (imposed_phase_index != iphase_not_imposed) {
                        double rhoc = rhomolar_critical();
                        if (imposed_phase_index == iphase_liquid && cached_single_phase_i > 0) {
                            // We want a liquid solution, but we got a vapor solution
                            if (_p < this->AS->p_critical()) {
                                while (cached_single_phase_i > 0
                                       && single_phase_logpT.rhomolar[cached_single_phase_i + 1][cached_single_phase_j] < rhoc) {
                                    // Bump to lower temperature
                                    cached_single_phase_i--;
                                }
                                double rho = evaluate_single_phase_pT(iDmolar, cached_single_phase_i, cached_single_phase_j);
                                if (rho < rhoc) {
                                    // Didn't work
                                    throw ValueError("Bump unsuccessful");
                                } else {
                                    _rhomolar = rho;
                                }
                            }
                        } else if ((imposed_phase_index == iphase_gas || imposed_phase_index == iphase_supercritical_gas)
                                   && cached_single_phase_i > 0) {

                            // We want a gas solution, but we got a liquid solution
                            if (_p < this->AS->p_critical()) {
                                while (cached_single_phase_i > 0
                                       && single_phase_logpT.rhomolar[cached_single_phase_i][cached_single_phase_j + 1] > rhoc) {
                                    // Bump to lower temperature
                                    cached_single_phase_i++;
                                }
                                double rho = evaluate_single_phase_pT(iDmolar, cached_single_phase_i, cached_single_phase_j);
                                if (rho > rhoc) {
                                    // Didn't work
                                    throw ValueError("Bump unsuccessful");
                                } else {
                                    _rhomolar = rho;
                                }
                            }
                        }
                    } else {

                        // If p < pc, you might be getting a liquid solution when you want a vapor solution or vice versa
                        // if you are very close to the saturation curve, so we figure out what the saturation temperature
                        // is for the given pressure
                        if (_p < this->AS->p_critical()) {
                            double Ts = pure_saturation.evaluate(iT, _p, _Q, iL, iV);
                            double TL = single_phase_logpT.T[cached_single_phase_i][cached_single_phase_j];
                            double TR = single_phase_logpT.T[cached_single_phase_i + 1][cached_single_phase_j];
                            if (TL < Ts && Ts < TR) {
                                if (_T < Ts) {
                                    if (cached_single_phase_i == 0) {
                                        throw ValueError(format("P, T are near saturation, but cannot move the cell to the left"));
                                    }
                                    // It's liquid, move the cell to the left
                                    cached_single_phase_i--;
                                } else {
                                    if (cached_single_phase_i > single_phase_logpT.Nx - 2) {
                                        throw ValueError(format("P,T are near saturation, but cannot move the cell to the right"));
                                    }
                                    // It's vapor, move to the right
                                    cached_single_phase_i++;
                                }
                            }
                        }
                    }
                    // Recalculate the phase
                    recalculate_singlephase_phase();
                }
            }
            break;
        }
        case PUmolar_INPUTS:
        case PSmolar_INPUTS:
        case DmolarP_INPUTS: {
            CoolPropDbl otherval;
            parameters otherkey;
            switch (input_pair) {
                case PUmolar_INPUTS:
                    _p = val1;
                    _umolar = val2;
                    otherval = val2;
                    otherkey = iUmolar;
                    break;
                case PSmolar_INPUTS:
                    _p = val1;
                    _smolar = val2;
                    otherval = val2;
                    otherkey = iSmolar;
                    break;
                case DmolarP_INPUTS:
                    _rhomolar = val1;
                    _p = val2;
                    otherval = val1;
                    otherkey = iDmolar;
                    break;
                default:
                    throw ValueError("Bad (impossible) pair");
            }

            using_single_phase_table = true;  // Use the table (or first guess is that it is single-phase)!
            std::size_t iL = std::numeric_limits<std::size_t>::max(), iV = std::numeric_limits<std::size_t>::max();
            CoolPropDbl zL = 0, zV = 0;
            std::size_t iclosest = 0;
            SimpleState closest_state;
            bool is_two_phase = false;
            // If phase is imposed, use it, but only if it's single phase:
            //   - Imposed two phase still needs to determine saturation limits by calling is_inside().
            //   - There's no speed increase to be gained by imposing two phase.
            if ((imposed_phase_index == iphase_not_imposed) || (imposed_phase_index == iphase_twophase)) {
                if (is_mixture) {
                    is_two_phase = PhaseEnvelopeRoutines::is_inside(phase_envelope, iP, _p, otherkey, otherval, iclosest, closest_state);
                    if (is_two_phase) {
                        std::vector<std::pair<std::size_t, std::size_t>> intersect =
                          PhaseEnvelopeRoutines::find_intersections(phase_envelope, iP, _p);
                        if (intersect.size() < 2) {
                            throw ValueError(format("p [%g Pa] is not within phase envelope", _p));
                        }
                        iV = intersect[0].first;
                        iL = intersect[1].first;
                        zL = PhaseEnvelopeRoutines::evaluate(phase_envelope, otherkey, iP, _p, iL);
                        zV = PhaseEnvelopeRoutines::evaluate(phase_envelope, otherkey, iP, _p, iV);
                    }
                } else {
                    is_two_phase = pure_saturation.is_inside(iP, _p, otherkey, otherval, iL, iV, zL, zV);
                }
            }
            // Phase determined or imposed, now interpolate results
            if (is_two_phase) {
                using_single_phase_table = false;
                if (otherkey == iDmolar) {
                    _Q = (1 / otherval - 1 / zL) / (1 / zV - 1 / zL);
                } else {
                    _Q = (otherval - zL) / (zV - zL);
                }
                if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))) {

                    throw ValueError(format("vapor quality is not in (0,1) for %s: %g p: %g", get_parameter_information(otherkey, "short").c_str(),
                                            otherval, static_cast<double>(_p)));
                } else if (!is_mixture) {
                    cached_saturation_iL = iL;
                    cached_saturation_iV = iV;
                }
                _phase = iphase_twophase;
            } else {
                selected_table = SELECTED_PH_TABLE;
                // Find and cache the indices i, j
                find_nearest_neighbor(single_phase_logph, dataset->coeffs_ph, iP, _p, otherkey, otherval, cached_single_phase_i,
                                      cached_single_phase_j);
                // Now find hmolar given P, X for X in Smolar, Umolar, Dmolar
                invert_single_phase_x(single_phase_logph, dataset->coeffs_ph, otherkey, otherval, _p, cached_single_phase_i, cached_single_phase_j);
                // Recalculate the phase
                recalculate_singlephase_phase();
            }
            break;
        }
        case SmolarT_INPUTS:
        case DmolarT_INPUTS: {
            CoolPropDbl otherval;
            parameters otherkey;
            switch (input_pair) {
                case SmolarT_INPUTS:
                    _smolar = val1;
                    _T = val2;
                    otherval = val1;
                    otherkey = iSmolar;
                    break;
                case DmolarT_INPUTS:
                    _rhomolar = val1;
                    _T = val2;
                    otherval = val1;
                    otherkey = iDmolar;
                    break;
                default:
                    throw ValueError("Bad (impossible) pair");
            }

            using_single_phase_table = true;  // Use the table (or first guess is that it is single-phase)!
            std::size_t iL = std::numeric_limits<std::size_t>::max(), iV = std::numeric_limits<std::size_t>::max();
            CoolPropDbl zL = 0, zV = 0;
            std::size_t iclosest = 0;
            SimpleState closest_state;
            bool is_two_phase = false;
            // If phase is imposed, use it, but only if it's single phase:
            //   - Imposed two phase still needs to determine saturation limits by calling is_inside().
            //   - There's no speed increase to be gained by imposing two phase.
            if ((imposed_phase_index == iphase_not_imposed) || (imposed_phase_index == iphase_twophase)) {
                if (is_mixture) {
                    is_two_phase = PhaseEnvelopeRoutines::is_inside(phase_envelope, iT, _T, otherkey, otherval, iclosest, closest_state);
                    if (is_two_phase) {
                        std::vector<std::pair<std::size_t, std::size_t>> intersect =
                          PhaseEnvelopeRoutines::find_intersections(phase_envelope, iT, _T);
                        if (intersect.size() < 2) {
                            throw ValueError(format("T [%g K] is not within phase envelope", _T));
                        }
                        iV = intersect[0].first;
                        iL = intersect[1].first;
                        zL = PhaseEnvelopeRoutines::evaluate(phase_envelope, otherkey, iT, _T, iL);
                        zV = PhaseEnvelopeRoutines::evaluate(phase_envelope, otherkey, iT, _T, iV);
                    }
                } else {
                    is_two_phase = pure_saturation.is_inside(iT, _T, otherkey, otherval, iL, iV, zL, zV);
                }
            }
            if (is_two_phase) {
                using_single_phase_table = false;
                if (otherkey == iDmolar) {
                    _Q = (1 / otherval - 1 / zL) / (1 / zV - 1 / zL);
                } else {
                    _Q = (otherval - zL) / (zV - zL);
                }
                if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))) {
                    throw ValueError(format("vapor quality is not in (0,1) for %s: %g T: %g", get_parameter_information(otherkey, "short").c_str(),
                                            otherval, static_cast<double>(_T)));
                } else if (!is_mixture) {
                    cached_saturation_iL = iL;
                    cached_saturation_iV = iV;
                    _p = pure_saturation.evaluate(iP, _T, _Q, iL, iV);
                } else {
                    // Mixture
                    std::vector<std::pair<std::size_t, std::size_t>> intersect = PhaseEnvelopeRoutines::find_intersections(phase_envelope, iT, _T);
                    if (intersect.size() < 2) {
                        throw ValueError(format("T [%g K] is not within phase envelope", _T));
                    }
                    iV = intersect[0].first;
                    iL = intersect[1].first;
                    CoolPropDbl pL = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iL);
                    CoolPropDbl pV = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iV);
                    _p = _Q * pV + (1 - _Q) * pL;
                }
            } else {
                selected_table = SELECTED_PT_TABLE;
                // Find and cache the indices i, j
                find_nearest_neighbor(single_phase_logpT, dataset->coeffs_pT, iT, _T, otherkey, otherval, cached_single_phase_i,
                                      cached_single_phase_j);
                // Now find the y variable (Dmolar or Smolar in this case)
                invert_single_phase_y(single_phase_logpT, dataset->coeffs_pT, otherkey, otherval, _T, cached_single_phase_i, cached_single_phase_j);
                // Recalculate the phase
                recalculate_singlephase_phase();
            }
            break;
        }
        case PQ_INPUTS: {
            std::size_t iL = 0, iV = 0;
            _p = val1;
            _Q = val2;
            using_single_phase_table = false;
            if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))) {
                throw ValueError(format("vapor quality [%g] is not in (0,1)", static_cast<double>(val2)));
            } else {
                CoolPropDbl TL = _HUGE, TV = _HUGE;
                if (is_mixture) {
                    std::vector<std::pair<std::size_t, std::size_t>> intersect = PhaseEnvelopeRoutines::find_intersections(phase_envelope, iP, _p);
                    if (intersect.size() < 2) {
                        throw ValueError(format("p [%g Pa] is not within phase envelope", _p));
                    }
                    iV = intersect[0].first;
                    iL = intersect[1].first;
                    TL = PhaseEnvelopeRoutines::evaluate(phase_envelope, iT, iP, _p, iL);
                    TV = PhaseEnvelopeRoutines::evaluate(phase_envelope, iT, iP, _p, iV);
                } else {
                    bool it_is_inside = pure_saturation.is_inside(iP, _p, iQ, _Q, iL, iV, TL, TV);
                    if (!it_is_inside) {
                        throw ValueError("Not possible to determine whether pressure is inside or not");
                    }
                }
                _T = _Q * TV + (1 - _Q) * TL;
                cached_saturation_iL = iL;
                cached_saturation_iV = iV;
                _phase = iphase_twophase;
            }
            break;
        }
        case QT_INPUTS: {
            std::size_t iL = 0, iV = 0;
            _Q = val1;
            _T = val2;

            using_single_phase_table = false;
            if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_Q))) {
                throw ValueError(format("vapor quality [%g] is not in (0,1)", static_cast<double>(val1)));
            } else {
                CoolPropDbl pL = _HUGE, pV = _HUGE;
                if (is_mixture) {
                    std::vector<std::pair<std::size_t, std::size_t>> intersect = PhaseEnvelopeRoutines::find_intersections(phase_envelope, iT, _T);
                    if (intersect.size() < 2) {
                        throw ValueError(format("T [%g K] is not within phase envelope", _T));
                    }
                    iV = intersect[0].first;
                    iL = intersect[1].first;
                    pL = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iL);
                    pV = PhaseEnvelopeRoutines::evaluate(phase_envelope, iP, iT, _T, iV);
                } else {
                    pure_saturation.is_inside(iT, _T, iQ, _Q, iL, iV, pL, pV);
                }
                _p = _Q * pV + (1 - _Q) * pL;
                cached_saturation_iL = iL;
                cached_saturation_iV = iV;
                _phase = iphase_twophase;
            }
            break;
        }
        default:
            throw ValueError("Sorry, but this set of inputs is not supported for Tabular backend");
    }
}

void CoolProp::TabularDataSet::write_tables(const std::string& path_to_tables) {
    make_dirs(path_to_tables);
    write_table(single_phase_logph, path_to_tables, "single_phase_logph");
    write_table(single_phase_logpT, path_to_tables, "single_phase_logpT");
    write_table(pure_saturation, path_to_tables, "pure_saturation");
    write_table(phase_envelope, path_to_tables, "phase_envelope");
}

void CoolProp::TabularDataSet::load_tables(const std::string& path_to_tables, shared_ptr<CoolProp::AbstractState>& AS) {
    single_phase_logph.AS = AS;
    single_phase_logpT.AS = AS;
    pure_saturation.AS = AS;
    single_phase_logph.set_limits();
    single_phase_logpT.set_limits();
    load_table(single_phase_logph, path_to_tables, "single_phase_logph.bin.z");
    load_table(single_phase_logpT, path_to_tables, "single_phase_logpT.bin.z");
    load_table(pure_saturation, path_to_tables, "pure_saturation.bin.z");
    load_table(phase_envelope, path_to_tables, "phase_envelope.bin.z");
    tables_loaded = true;
    if (get_debug_level() > 0) {
        std::cout << "Tables loaded" << std::endl;
    }
};

void CoolProp::TabularDataSet::build_tables(shared_ptr<CoolProp::AbstractState>& AS) {
    // Pure or pseudo-pure fluid
    if (AS->get_mole_fractions().size() == 1) {
        pure_saturation.build(AS);
    } else {
        // Call function to actually construct the phase envelope
        AS->build_phase_envelope("");
        // Copy constructed phase envelope into this class
        PhaseEnvelopeData PED = AS->get_phase_envelope_data();
        // Convert into packable form
        phase_envelope.copy_from_nonpackable(PED);
        // Resize so that it will load properly
        pure_saturation.resize(pure_saturation.N);
    }
    single_phase_logph.build(AS);
    single_phase_logpT.build(AS);
    tables_loaded = true;
}

/// Return the set of tabular datasets
CoolProp::TabularDataSet* CoolProp::TabularDataLibrary::get_set_of_tables(shared_ptr<AbstractState>& AS, bool& loaded) {
    const std::string path = path_to_tables(AS);
    // Try to find tabular set if it is already loaded
    std::map<std::string, TabularDataSet>::iterator it = data.find(path);
    // It is already in the map, return it
    if (it != data.end()) {
        loaded = it->second.tables_loaded;
        return &(it->second);
    }
    // It is not in the map, build it
    else {
        TabularDataSet set;
        data.insert(std::pair<std::string, TabularDataSet>(path, set));
        TabularDataSet& dataset = data[path];
        try {
            if (!dataset.tables_loaded) {
                dataset.load_tables(path, AS);
            }
            loaded = true;
        } catch (std::exception&) {
            loaded = false;
        }
        return &(dataset);
    }
}

void CoolProp::TabularDataSet::build_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs) {
    if (!coeffs.empty()) {
        return;
    }
    const bool debug = get_debug_level() > 5 || false;
    const int param_count = 6;
    parameters param_list[param_count] = {iDmolar, iT, iSmolar, iHmolar, iP, iUmolar};
    std::vector<std::vector<double>>*f = NULL, *fx = NULL, *fy = NULL, *fxy = NULL;

    clock_t t1 = clock();

    // Resize the coefficient structures
    coeffs.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));

    int valid_cell_count = 0;
    for (std::size_t k = 0; k < param_count; ++k) {
        parameters param = param_list[k];
        if (param == table.xkey || param == table.ykey) {
            continue;
        }  // Skip tables that match either of the input variables

        switch (param) {
            case iT:
                f = &(table.T);
                fx = &(table.dTdx);
                fy = &(table.dTdy);
                fxy = &(table.d2Tdxdy);
                break;
            case iP:
                f = &(table.p);
                fx = &(table.dpdx);
                fy = &(table.dpdy);
                fxy = &(table.d2pdxdy);
                break;
            case iDmolar:
                f = &(table.rhomolar);
                fx = &(table.drhomolardx);
                fy = &(table.drhomolardy);
                fxy = &(table.d2rhomolardxdy);
                break;
            case iSmolar:
                f = &(table.smolar);
                fx = &(table.dsmolardx);
                fy = &(table.dsmolardy);
                fxy = &(table.d2smolardxdy);
                break;
            case iHmolar:
                f = &(table.hmolar);
                fx = &(table.dhmolardx);
                fy = &(table.dhmolardy);
                fxy = &(table.d2hmolardxdy);
                break;
            case iUmolar:
                f = &(table.umolar);
                fx = &(table.dumolardx);
                fy = &(table.dumolardy);
                fxy = &(table.d2umolardxdy);
                break;
            default:
                throw ValueError("Invalid variable type to build_coeffs");
        }
        for (std::size_t i = 0; i < table.Nx - 1; ++i)  // -1 since we have one fewer cells than nodes
        {
            for (std::size_t j = 0; j < table.Ny - 1; ++j)  // -1 since we have one fewer cells than nodes
            {
                if (ValidNumber((*f)[i][j]) && ValidNumber((*f)[i + 1][j]) && ValidNumber((*f)[i][j + 1]) && ValidNumber((*f)[i + 1][j + 1])) {

                    // This will hold the scaled f values for the cell
                    Eigen::Matrix<double, 16, 1> F;
                    // The output values (do not require scaling
                    F(0) = (*f)[i][j];
                    F(1) = (*f)[i + 1][j];
                    F(2) = (*f)[i][j + 1];
                    F(3) = (*f)[i + 1][j + 1];
                    // Scaling parameter
                    // d(f)/dxhat = df/dx * dx/dxhat, where xhat = (x-x_i)/(x_{i+1}-x_i)
                    coeffs[i][j].dx_dxhat = table.xvec[i + 1] - table.xvec[i];
                    double dx_dxhat = coeffs[i][j].dx_dxhat;
                    F(4) = (*fx)[i][j] * dx_dxhat;
                    F(5) = (*fx)[i + 1][j] * dx_dxhat;
                    F(6) = (*fx)[i][j + 1] * dx_dxhat;
                    F(7) = (*fx)[i + 1][j + 1] * dx_dxhat;
                    // Scaling parameter
                    // d(f)/dyhat = df/dy * dy/dyhat, where yhat = (y-y_j)/(y_{j+1}-y_j)
                    coeffs[i][j].dy_dyhat = table.yvec[j + 1] - table.yvec[j];
                    double dy_dyhat = coeffs[i][j].dy_dyhat;
                    F(8) = (*fy)[i][j] * dy_dyhat;
                    F(9) = (*fy)[i + 1][j] * dy_dyhat;
                    F(10) = (*fy)[i][j + 1] * dy_dyhat;
                    F(11) = (*fy)[i + 1][j + 1] * dy_dyhat;
                    // Cross derivatives are doubly scaled following the examples above
                    F(12) = (*fxy)[i][j] * dy_dyhat * dx_dxhat;
                    F(13) = (*fxy)[i + 1][j] * dy_dyhat * dx_dxhat;
                    F(14) = (*fxy)[i][j + 1] * dy_dyhat * dx_dxhat;
                    F(15) = (*fxy)[i + 1][j + 1] * dy_dyhat * dx_dxhat;
                    // Calculate the alpha coefficients
                    Eigen::MatrixXd alpha = Ainv.transpose() * F;  // 16x1; Watch out for the transpose!
                    std::vector<double> valpha = eigen_to_vec1D(alpha);
                    coeffs[i][j].set(param, valpha);
                    coeffs[i][j].set_valid();
                    valid_cell_count++;
                } else {
                    coeffs[i][j].set_invalid();
                }
            }
        }
        double elapsed = (clock() - t1) / ((double)CLOCKS_PER_SEC);
        if (debug) {
            std::cout << format("Calculated bicubic coefficients for %d good cells in %g sec.\n", valid_cell_count, elapsed);
        }
        std::size_t remap_count = 0;
        // Now find invalid cells and give them pointers to a neighboring cell that works
        for (std::size_t i = 0; i < table.Nx - 1; ++i)  // -1 since we have one fewer cells than nodes
        {
            for (std::size_t j = 0; j < table.Ny - 1; ++j)  // -1 since we have one fewer cells than nodes
            {
                // Not a valid cell
                if (!coeffs[i][j].valid()) {
                    // Offsets that we are going to try in order (left, right, top, bottom, diagonals)
                    int xoffsets[] = {-1, 1, 0, 0, -1, 1, 1, -1};
                    int yoffsets[] = {0, 0, 1, -1, -1, -1, 1, 1};
                    // Length of offset
                    std::size_t N = sizeof(xoffsets) / sizeof(xoffsets[0]);
                    for (std::size_t k = 0; k < N; ++k) {
                        std::size_t iplus = i + xoffsets[k];
                        std::size_t jplus = j + yoffsets[k];
                        if (0 < iplus && iplus < table.Nx - 1 && 0 < jplus && jplus < table.Ny - 1 && coeffs[iplus][jplus].valid()) {
                            coeffs[i][j].set_alternate(iplus, jplus);
                            remap_count++;
                            if (debug) {
                                std::cout << format("Mapping %d,%d to %d,%d\n", i, j, iplus, jplus);
                            }
                            break;
                        }
                    }
                }
            }
        }
        if (debug) {
            std::cout << format("Remapped %d cells\n", remap_count);
        }
    }
}

#    if defined(ENABLE_CATCH)
#        include <catch2/catch_all.hpp>

// Defined global so we only load once
static shared_ptr<CoolProp::AbstractState> ASHEOS, ASTTSE, ASBICUBIC;

/* Use a fixture so that loading of the tables from memory only happens once in the initializer */
class TabularFixture
{
   public:
    TabularFixture() {}
    void setup() {
        if (ASHEOS.get() == NULL) {
            ASHEOS.reset(CoolProp::AbstractState::factory("HEOS", "Water"));
        }
        if (ASTTSE.get() == NULL) {
            ASTTSE.reset(CoolProp::AbstractState::factory("TTSE&HEOS", "Water"));
        }
        if (ASBICUBIC.get() == NULL) {
            ASBICUBIC.reset(CoolProp::AbstractState::factory("BICUBIC&HEOS", "Water"));
        }
    }
};
TEST_CASE_METHOD(TabularFixture, "Tests for tabular backends with water", "[Tabular]") {
    SECTION("first_saturation_deriv invalid quality") {
        setup();
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 0.5);
        CHECK_THROWS(ASBICUBIC->first_saturation_deriv(CoolProp::iP, CoolProp::iT));
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 0.5);
        CHECK_THROWS(ASTTSE->first_saturation_deriv(CoolProp::iP, CoolProp::iT));
    }

    SECTION("first_saturation_deriv dp/dT") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl expected = ASHEOS->first_saturation_deriv(CoolProp::iP, CoolProp::iT);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_TTSE = ASTTSE->first_saturation_deriv(CoolProp::iP, CoolProp::iT);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_BICUBIC = ASTTSE->first_saturation_deriv(CoolProp::iP, CoolProp::iT);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_saturation_deriv dDmolar/dP") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl expected = ASHEOS->first_saturation_deriv(CoolProp::iDmolar, CoolProp::iP);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_TTSE = ASTTSE->first_saturation_deriv(CoolProp::iDmolar, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_saturation_deriv(CoolProp::iDmolar, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_saturation_deriv dDmass/dP") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl expected = ASHEOS->first_saturation_deriv(CoolProp::iDmass, CoolProp::iP);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_TTSE = ASTTSE->first_saturation_deriv(CoolProp::iDmass, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_saturation_deriv(CoolProp::iDmass, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_saturation_deriv dHmass/dP") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl expected = ASHEOS->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_TTSE = ASTTSE->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_saturation_deriv dHmass/dP w/ QT as inputs") {
        setup();
        ASHEOS->update(CoolProp::QT_INPUTS, 1, 370);
        CoolPropDbl expected = ASHEOS->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        ASTTSE->update(CoolProp::QT_INPUTS, 1, 370);
        CoolPropDbl actual_TTSE = ASTTSE->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        ASBICUBIC->update(CoolProp::QT_INPUTS, 1, 370);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_saturation_deriv(CoolProp::iHmass, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_two_phase_deriv dDmolar/dP|Hmolar") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl expected = ASHEOS->first_two_phase_deriv(CoolProp::iDmolar, CoolProp::iP, CoolProp::iHmolar);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_TTSE = ASTTSE->first_two_phase_deriv(CoolProp::iDmolar, CoolProp::iP, CoolProp::iHmolar);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_two_phase_deriv(CoolProp::iDmolar, CoolProp::iP, CoolProp::iHmolar);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_two_phase_deriv dDmass/dP|Hmass") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl expected = ASHEOS->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iHmass);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_TTSE = ASTTSE->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iHmass);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iHmass);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_two_phase_deriv dDmass/dHmass|P") {
        setup();
        ASHEOS->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl expected = ASHEOS->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iHmass, CoolProp::iP);
        ASTTSE->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_TTSE = ASTTSE->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iHmass, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PQ_INPUTS, 101325, 0.1);
        CoolPropDbl actual_BICUBIC = ASBICUBIC->first_two_phase_deriv(CoolProp::iDmass, CoolProp::iHmass, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-6);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-6);
    }
    SECTION("first_partial_deriv dHmass/dT|P") {
        setup();
        ASHEOS->update(CoolProp::PT_INPUTS, 101325, 300);
        double expected = ASHEOS->cpmass();
        ASTTSE->update(CoolProp::PT_INPUTS, 101325, 300);
        double dhdT_TTSE = ASTTSE->first_partial_deriv(CoolProp::iHmass, CoolProp::iT, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PT_INPUTS, 101325, 300);
        double dhdT_BICUBIC = ASBICUBIC->first_partial_deriv(CoolProp::iHmass, CoolProp::iT, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(dhdT_TTSE);
        CAPTURE(dhdT_BICUBIC);
        CHECK(std::abs((expected - dhdT_TTSE) / expected) < 1e-4);
        CHECK(std::abs((expected - dhdT_BICUBIC) / expected) < 1e-4);
    }
    SECTION("first_partial_deriv dHmolar/dT|P") {
        setup();
        ASHEOS->update(CoolProp::PT_INPUTS, 101325, 300);
        double expected = ASHEOS->cpmolar();
        ASTTSE->update(CoolProp::PT_INPUTS, 101325, 300);
        double dhdT_TTSE = ASTTSE->first_partial_deriv(CoolProp::iHmolar, CoolProp::iT, CoolProp::iP);
        ASBICUBIC->update(CoolProp::PT_INPUTS, 101325, 300);
        double dhdT_BICUBIC = ASBICUBIC->first_partial_deriv(CoolProp::iHmolar, CoolProp::iT, CoolProp::iP);
        CAPTURE(expected);
        CAPTURE(dhdT_TTSE);
        CAPTURE(dhdT_BICUBIC);
        CHECK(std::abs((expected - dhdT_TTSE) / expected) < 1e-4);
        CHECK(std::abs((expected - dhdT_BICUBIC) / expected) < 1e-4);
    }
    SECTION("check isentropic process") {
        setup();
        double T0 = 300;
        double p0 = 1e5;
        double p1 = 1e6;
        ASHEOS->update(CoolProp::PT_INPUTS, p0, 300);
        double s0 = ASHEOS->smolar();
        ASHEOS->update(CoolProp::PSmolar_INPUTS, p1, s0);
        double expected = ASHEOS->T();
        ASTTSE->update(CoolProp::PSmolar_INPUTS, p1, s0);
        double actual_TTSE = ASTTSE->T();
        ASBICUBIC->update(CoolProp::PSmolar_INPUTS, p1, s0);
        double actual_BICUBIC = ASBICUBIC->T();
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-2);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-2);
    }
    SECTION("check D=1 mol/m3, T=500 K inputs") {
        setup();
        double d = 1;
        CAPTURE(d);
        ASHEOS->update(CoolProp::DmolarT_INPUTS, d, 500);
        double expected = ASHEOS->p();
        ASTTSE->update(CoolProp::DmolarT_INPUTS, d, 500);
        double actual_TTSE = ASTTSE->p();
        ASBICUBIC->update(CoolProp::DmolarT_INPUTS, d, 500);
        double actual_BICUBIC = ASBICUBIC->p();
        CAPTURE(expected);
        CAPTURE(actual_TTSE);
        CAPTURE(actual_BICUBIC);
        CHECK(std::abs((expected - actual_TTSE) / expected) < 1e-3);
        CHECK(std::abs((expected - actual_BICUBIC) / expected) < 1e-3);
    }
}
#    endif  // ENABLE_CATCH

#endif  // !defined(NO_TABULAR_BACKENDS)
