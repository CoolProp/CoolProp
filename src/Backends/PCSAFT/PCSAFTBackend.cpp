#include <vector>
#include <string>
#include <cmath>
#include "math.h"
#include <Eigen/Dense>
#include <stdlib.h>

#include "PCSAFTBackend.h"
#include "Backends/Helmholtz/VLERoutines.h"

using std::vector;

/*
References
----------
* J. Gross and G. Sadowski, “Perturbed-Chain SAFT:  An Equation of State
Based on a Perturbation Theory for Chain Molecules,” Ind. Eng. Chem.
Res., vol. 40, no. 4, pp. 1244–1260, Feb. 2001.
* M. Kleiner and G. Sadowski, “Modeling of Polar Systems Using PCP-SAFT: 
An Approach to Account for Induced-Association Interactions,” J. Phys.
Chem. C, vol. 111, no. 43, pp. 15544–15553, Nov. 2007.
* Gross Joachim and Vrabec Jadran, “An equation‐of‐state contribution
for polar components: Dipolar molecules,” AIChE J., vol. 52, no. 3,
pp. 1194–1204, Feb. 2006.
* A. J. de Villiers, C. E. Schwarz, and A. J. Burger, “Improving
vapour–liquid-equilibria predictions for mixtures with non-associating polar
components using sPC-SAFT extended with two dipolar terms,” Fluid Phase
Equilibria, vol. 305, no. 2, pp. 174–184, Jun. 2011.
* S. H. Huang and M. Radosz, “Equation of state for small, large,
polydisperse, and associating molecules,” Ind. Eng. Chem. Res., vol. 29,
no. 11, pp. 2284–2294, Nov. 1990.
* S. H. Huang and M. Radosz, “Equation of state for small, large,
polydisperse, and associating molecules: extension to fluid mixtures,”
Ind. Eng. Chem. Res., vol. 30, no. 8, pp. 1994–2005, Aug. 1991.
* S. H. Huang and M. Radosz, “Equation of state for small, large,
polydisperse, and associating molecules: extension to fluid mixtures.
[Erratum to document cited in CA115(8):79950j],” Ind. Eng. Chem. Res.,
vol. 32, no. 4, pp. 762–762, Apr. 1993.
* J. Gross and G. Sadowski, “Application of the Perturbed-Chain SAFT
Equation of State to Associating Systems,” Ind. Eng. Chem. Res., vol.
41, no. 22, pp. 5510–5515, Oct. 2002.
* L. F. Cameretti, G. Sadowski, and J. M. Mollerup, “Modeling of Aqueous
Electrolyte Solutions with Perturbed-Chain Statistical Associated Fluid
Theory,” Ind. Eng. Chem. Res., vol. 44, no. 9, pp. 3355–3362, Apr. 2005.
* L. F. Cameretti, G. Sadowski, and J. M. Mollerup, “Modeling of Aqueous
Electrolyte Solutions with Perturbed-Chain Statistical Association Fluid
Theory,” Ind. Eng. Chem. Res., vol. 44, no. 23, pp. 8944–8944, Nov. 2005.
* C. Held, L. F. Cameretti, and G. Sadowski, “Modeling aqueous
electrolyte solutions: Part 1. Fully dissociated electrolytes,” Fluid
Phase Equilibria, vol. 270, no. 1, pp. 87–96, Aug. 2008.
* C. Held, T. Reschke, S. Mohammad, A. Luza, and G. Sadowski, “ePC-SAFT
revised,” Chem. Eng. Res. Des., vol. 92, no. 12, pp. 2884–2897, Dec. 2014.
*/

namespace CoolProp {

PCSAFTBackend::PCSAFTBackend(const std::vector<std::string>& component_names, bool generate_SatL_and_SatV) {
    N = component_names.size();
    components.resize(N);
    ion_term = false;
    polar_term = false;
    assoc_term = false;
    water_present = false;
    water_idx = 0;
    for (unsigned int i = 0; i < N; ++i){
        components[i] = PCSAFTLibrary::get_library().get(component_names[i]);
        // Determining which PC-SAFT terms should be used
        if (components[i].getZ() != 0) {
            ion_term = true;
        }
        if (components[i].getDipm() != 0) {
            polar_term = true;
        }
        if (components[i].getVolA() != 0) {
            assoc_term = true;
        }
        if (components[i].getCAS() == "7732-18-5") {
            water_present = true;
            water_idx = i;
        }
    }

    // Set up association scheme
    if (assoc_term) {
        set_assoc_matrix();
    }

    // Set the components and associated flags
    is_pure_or_pseudopure = (N == 1);

    // loading interaction parameters
    std::string kij_string;
    std::string kijT_string;
    if (is_pure_or_pseudopure) {
        this->mole_fractions = std::vector<CoolPropDbl>(1, 1);
    } else {
        k_ij.resize(N * N, 0.0);
        k_ijT.resize(N * N, 0.0);
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < N; ++j) {
                if (i != j) {
                    kij_string = PCSAFTLibrary::get_library().get_binary_interaction_pcsaft(components[i].getCAS(), components[j].getCAS(), "kij");
                    kijT_string = PCSAFTLibrary::get_library().get_binary_interaction_pcsaft(components[i].getCAS(), components[j].getCAS(), "kijT");
                    k_ij[i * N + j] = atof(kij_string.c_str());
                    k_ijT[i * N + j] = atof(kijT_string.c_str());
                }
            }
        }
    }

    if (generate_SatL_and_SatV) {
        bool SatLSatV = false;
        SatL.reset(this->get_copy(SatLSatV));
        SatL->specify_phase(iphase_liquid);
        SatV.reset(this->get_copy(SatLSatV));
        SatV->specify_phase(iphase_gas);
    }

    // Set the phase to default unknown value
    imposed_phase_index = iphase_not_imposed;
    _phase = iphase_unknown;
}

PCSAFTBackend::PCSAFTBackend(const std::vector<PCSAFTFluid>& components_in, bool generate_SatL_and_SatV) {
    components = components_in;
    N = components.size();
    // Determining which PC-SAFT terms should be used
    ion_term = false;
    polar_term = false;
    assoc_term = false;
    water_present = false;
    water_idx = 0;
    for (unsigned int i = 0; i < N; ++i){
        if (components[i].getZ() != 0) {
            ion_term = true;
        }
        if (components[i].getDipm() != 0) {
            polar_term = true;
        }
        if (components[i].getVolA() != 0) {
            assoc_term = true;
        }
        if (components[i].getCAS() == "7732-18-5") {
            water_present = true;
            water_idx = i;
        }
    }

    // Set up association scheme
    if (assoc_term) {
        set_assoc_matrix();
    }

    // Set the components and associated flags
    is_pure_or_pseudopure = (N == 1);

    // loading interaction parameters
    std::string kij_string;
    std::string kijT_string;
    if (is_pure_or_pseudopure) {
        this->mole_fractions = std::vector<CoolPropDbl>(1, 1);
    } else {
        k_ij.resize(N * N, 0.0);
        k_ijT.resize(N * N, 0.0);
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < N; ++j) {
                if (i != j) {
                    kij_string = PCSAFTLibrary::get_library().get_binary_interaction_pcsaft(components[i].getCAS(), components[j].getCAS(), "kij");
                    kijT_string = PCSAFTLibrary::get_library().get_binary_interaction_pcsaft(components[i].getCAS(), components[j].getCAS(), "kijT");
                    k_ij[i * N + j] = atof(kij_string.c_str());
                    k_ijT[i * N + j] = atof(kijT_string.c_str());
                }
            }
        }
    }

    if (generate_SatL_and_SatV) {
        bool SatLSatV = false;
        SatL.reset(this->get_copy(SatLSatV));
        SatL->specify_phase(iphase_liquid);
        SatV.reset(this->get_copy(SatLSatV));
        SatV->specify_phase(iphase_gas);
    }

    // Set the phase to default unknown value
    imposed_phase_index = iphase_not_imposed;
    _phase = iphase_unknown;
}

void PCSAFTBackend::set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
    if (mole_fractions.size() != N) {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]", mole_fractions.size(), N));
    }
    // Copy values without reallocating memory
    this->mole_fractions = mole_fractions;  // Most effective copy
    this->resize(N);                        // No reallocation of this->mole_fractions happens
    // Also store the mole fractions as doubles
    this->mole_fractions_double = std::vector<double>(mole_fractions.begin(), mole_fractions.end());
};

void PCSAFTBackend::set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
    if (mass_fractions.size() != N) {
        throw ValueError(format("size of mass fraction vector [%d] does not equal that of component vector [%d]", mass_fractions.size(), N));
    }
    std::vector<CoolPropDbl> moles;
    CoolPropDbl sum_moles = 0.0;
    CoolPropDbl tmp = 0.0;
    for (unsigned int i = 0; i < components.size(); ++i) {
        tmp = mass_fractions[i] / components[i].molar_mass();
        moles.push_back(tmp);
        sum_moles += tmp;
    }
    std::vector<CoolPropDbl> mole_fractions;
    for (std::vector<CoolPropDbl>::iterator it = moles.begin(); it != moles.end(); ++it) {
        mole_fractions.push_back(*it / sum_moles);
    }
    this->set_mole_fractions(mole_fractions);
};

PCSAFTBackend* PCSAFTBackend::get_copy(bool generate_SatL_and_SatV) {
    // Set up the class with these components
    PCSAFTBackend* ptr = new PCSAFTBackend(components, generate_SatL_and_SatV);
    return ptr;
};

void PCSAFTBackend::resize(std::size_t N) {
    this->mole_fractions.resize(N);
    this->mole_fractions_double.resize(N);
    this->K.resize(N);
    this->lnK.resize(N);
}

CoolPropDbl PCSAFTBackend::update_DmolarT(CoolPropDbl rho) {
    _rhomolar = rho;
    return this->calc_pressure();
}

CoolPropDbl PCSAFTBackend::calc_pressure(void) {
    double den = _rhomolar * N_AV / 1.0e30;

    CoolPropDbl Z = this->calc_compressibility_factor();
    CoolPropDbl P = Z * kb * _T * den * 1.0e30;  // Pa
    return P;
}

CoolPropDbl PCSAFTBackend::calc_alphar(void) {
    int ncomp = N;  // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = components[i].getSigma() * (1 - 0.12 * exp(-3 * components[i].getU() / _T));
    }
    if (ion_term) {
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getZ() != 0) {
                d[i] =
                  components[i].getSigma() * (1 - 0.12);  // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }

    double den = _rhomolar * N_AV / 1.0e30;

    vector<double> zeta(4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j] * components[j].getM() * pow(d[j], i);
        }
        zeta[i] = PI / 6 * den * summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += mole_fractions[i] * components[i].getM();
    }

    vector<double> ghs(ncomp * ncomp, 0);
    vector<double> e_ij(ncomp * ncomp, 0);
    vector<double> s_ij(ncomp * ncomp, 0);
    double m2es3 = 0.;
    double m2e2s3 = 0.;
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            s_ij[idx] = (components[i].getSigma() + components[j].getSigma()) / 2.;
            if (ion_term) {
                if (components[i].getZ() * components[j].getZ()
                    <= 0) {  // for two cations or two anions e_ij is kept at zero to avoid dispersion between like ions (see Held et al. 2014)
                    if (k_ij.empty()) {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                    } else {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                    }
                }
            } else {
                if (k_ij.empty()) {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                } else {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                }
            }
            m2es3 = m2es3 + mole_fractions[i] * mole_fractions[j] * components[i].getM() * components[j].getM() * e_ij[idx] / _T * pow(s_ij[idx], 3);
            m2e2s3 =
              m2e2s3
              + mole_fractions[i] * mole_fractions[j] * components[i].getM() * components[j].getM() * pow(e_ij[idx] / _T, 2) * pow(s_ij[idx], 3);
            ghs[idx] = 1 / (1 - zeta[3]) + (d[i] * d[j] / (d[i] + d[j])) * 3 * zeta[2] / (1 - zeta[3]) / (1 - zeta[3])
                       + pow(d[i] * d[j] / (d[i] + d[j]), 2) * 2 * zeta[2] * zeta[2] / pow(1 - zeta[3], 3);
        }
    }

    double ares_hs = 1 / zeta[0]
                     * (3 * zeta[1] * zeta[2] / (1 - zeta[3]) + pow(zeta[2], 3.) / (zeta[3] * pow(1 - zeta[3], 2))
                        + (pow(zeta[2], 3.) / pow(zeta[3], 2.) - zeta[0]) * log(1 - zeta[3]));

    static double a0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
    static double a1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
    static double a2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
    static double b0[7] = { 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
    static double b1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
    static double b2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };

    vector<double> a(7, 0);
    vector<double> b(7, 0);
    for (int i = 0; i < 7; i++) {
        a[i] = a0[i] + (m_avg - 1.) / m_avg * a1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * a2[i];
        b[i] = b0[i] + (m_avg - 1.) / m_avg * b1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * b2[i];
    }

    double I1 = 0.0;
    double I2 = 0.0;
    for (int i = 0; i < 7; i++) {
        I1 += a[i] * pow(eta, i);
        I2 += b[i] * pow(eta, i);
    }
    double C1 = 1.
                / (1. + m_avg * (8 * eta - 2 * eta * eta) / pow(1 - eta, 4)
                   + (1 - m_avg) * (20 * eta - 27 * eta * eta + 12 * pow(eta, 3) - 2 * pow(eta, 4)) / pow((1 - eta) * (2 - eta), 2.0));

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i] * (components[i].getM() - 1) * log(ghs[i * ncomp + i]);
    }

    double ares_hc = m_avg * ares_hs - summ;
    double ares_disp = -2 * PI * den * I1 * m2es3 - PI * den * m_avg * C1 * I2 * m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double ares_polar = 0.;
    if (polar_term) {
        double A2 = 0.;
        double A3 = 0.;
        vector<double> dipmSQ(ncomp, 0);

        static double a0dip[5] = {0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308};
        static double a1dip[5] = {0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135};
        static double a2dip[5] = {-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575};
        static double b0dip[5] = {0.2187939, -1.1896431, 1.1626889, 0, 0};
        static double b1dip[5] = {-0.5873164, 1.2489132, -0.5085280, 0, 0};
        static double b2dip[5] = {3.4869576, -14.915974, 15.372022, 0, 0};
        static double c0dip[5] = {-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0};
        static double c1dip[5] = {-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0};
        static double c2dip[5] = {-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0};

        const static double conv = 7242.702976750923;  // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(components[i].getDipm(), 2.) / (components[i].getM() * components[i].getU() * pow(components[i].getSigma(), 3.)) * conv;
        }

        vector<double> adip(5, 0);
        vector<double> bdip(5, 0);
        vector<double> cdip(5, 0);
        double J2, J3;
        double m_ij;
        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(components[i].getM() * components[j].getM());
                if (m_ij > 2) {
                    m_ij = 2;
                }
                J2 = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij-1)/m_ij*a1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij-1)/m_ij*b1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*b2dip[l];
                    J2 += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                }
                A2 += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                      * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i]
                      * dipmSQ[j] * J2;

                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((components[i].getM() * components[j].getM() * components[k].getM()), 1 / 3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    J3 = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk - 1) / m_ijk * c1dip[l] + (m_ijk - 1) / m_ijk * (m_ijk - 2) / m_ijk * c2dip[l];
                        J3 += cdip[l] * pow(eta, l);
                    }
                    A3 += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                          * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                          / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum()
                          * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * J3;
                }
            }
        }

        A2 = -PI * den * A2;
        A3 = -4 / 3. * PI * PI * den * den * A3;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            ares_polar = A2/(1-A3/A2);
        }
    }

    // Association term -------------------------------------------------------
    double ares_assoc = 0.;
    if (assoc_term) {
        int num_sites = 0;
        vector<int> iA; //indices of associating compounds
        for(std::vector<int>::iterator it = assoc_num.begin(); it != assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(it - assoc_num.begin());
            }
        }

        vector<double> x_assoc(num_sites); // mole fractions of only the associating compounds
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        // these indices are necessary because we are only using 1D vectors
        vector<double> XA (num_sites, 0);
        vector<double> delta_ij(num_sites * num_sites, 0);
        int idxa = 0;
        int idxi = 0; // index for the ii-th compound
        int idxj = 0; // index for the jj-th compound
        for (int i = 0; i < num_sites; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                idxj = iA[j]*ncomp+iA[j];
                if (assoc_matrix[idxa] != 0) {
                    double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                    double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                        s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                    delta_ij[idxa] = ghs[iA[i]*ncomp+iA[j]]*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                }
                idxa += 1;
            }
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        ares_assoc = 0.;
        for (int i = 0; i < num_sites; i++) {
            ares_assoc += mole_fractions[iA[i]]*(log(XA[i])-0.5*XA[i] + 0.5);
        }
    }

    // Ion term ---------------------------------------------------------------
    double ares_ion = 0.;
    if (ion_term) {
        vector<double> q(ncomp);
        for (int i = 0; i < ncomp; i++) {
            q[i] = components[i].getZ() * E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += components[i].getZ() * components[i].getZ() * mole_fractions[i];
        }
        double kappa =
          sqrt(den * E_CHRG * E_CHRG / kb / _T / (dielc * perm_vac) * summ);  // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            vector<double> chi(ncomp);
            vector<double> sigma_k(ncomp);
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi[i] = 3 / pow(kappa * components[i].getSigma(), 3)
                         * (1.5 + log(1 + kappa * components[i].getSigma()) - 2 * (1 + kappa * components[i].getSigma())
                            + 0.5 * pow(1 + kappa * components[i].getSigma(), 2));
                summ += mole_fractions[i] * q[i] * q[i] * chi[i] * kappa;
            }

            ares_ion = -1 / 12. / PI / kb / _T / (dielc * perm_vac) * summ;
        }
    }

    CoolPropDbl ares = ares_hc + ares_disp + ares_polar + ares_assoc + ares_ion;
    return ares;
}

CoolPropDbl PCSAFTBackend::calc_dadt(void) {
    int ncomp = N;  // number of components
    vector<double> d(ncomp), dd_dt(ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = components[i].getSigma() * (1 - 0.12 * exp(-3 * components[i].getU() / _T));
        dd_dt[i] = components[i].getSigma() * -3 * components[i].getU() / _T / _T * 0.12 * exp(-3 * components[i].getU() / _T);
    }
    if (ion_term) {
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getZ() != 0) {
                d[i] =
                  components[i].getSigma() * (1 - 0.12);  // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
                dd_dt[i] = 0.;
            }
        }
    }

    double den = _rhomolar * N_AV / 1.0e30;

    vector<double> zeta(4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j] * components[j].getM() * pow(d[j], i);
        }
        zeta[i] = PI / 6 * den * summ;
    }

    vector<double> dzeta_dt(4, 0);
    for (int i = 1; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j] * components[j].getM() * i * dd_dt[j] * pow(d[j], (i - 1));
        }
        dzeta_dt[i] = PI / 6 * den * summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += mole_fractions[i] * components[i].getM();
    }

    vector<double> ghs(ncomp * ncomp, 0);
    vector<double> dghs_dt(ncomp * ncomp, 0);
    vector<double> e_ij(ncomp * ncomp, 0);
    vector<double> s_ij(ncomp * ncomp, 0);
    double m2es3 = 0.;
    double m2e2s3 = 0.;
    double ddij_dt;
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            s_ij[idx] = (components[i].getSigma() + components[j].getSigma()) / 2.;
            if (ion_term) {
                if (components[i].getZ() * components[j].getZ()
                    <= 0) {  // for two cations or two anions e_ij is kept at zero to avoid dispersion between like ions (see Held et al. 2014)
                    if (k_ij.empty()) {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                    } else {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                    }
                }
            } else {
                if (k_ij.empty()) {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                } else {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                }
            }
            m2es3 = m2es3 + mole_fractions[i] * mole_fractions[j] * components[i].getM() * components[j].getM() * e_ij[idx] / _T * pow(s_ij[idx], 3);
            m2e2s3 = m2e2s3 + mole_fractions[i] * mole_fractions[j] * components[i].getM() * components[j].getM() * pow(e_ij[idx] / _T, 2) * pow(s_ij[idx], 3);
            ghs[idx] = 1 / (1-zeta[3]) + (d[i]*d[j]/(d[i]+d[j]))*3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                    pow(d[i]*d[j]/(d[i]+d[j]), 2)*2*zeta[2]*zeta[2]/pow(1-zeta[3], 3);
            ddij_dt = (d[i]*d[j]/(d[i]+d[j]))*(dd_dt[i]/d[i]+dd_dt[j]/d[j]-(dd_dt[i]+dd_dt[j])/(d[i]+d[j]));
            dghs_dt[idx] = dzeta_dt[3]/pow(1-zeta[3], 2.)
                + 3*(ddij_dt*zeta[2]+(d[i]*d[j]/(d[i]+d[j]))*dzeta_dt[2])/pow(1-zeta[3], 2.)
                + 4*(d[i]*d[j]/(d[i]+d[j]))*zeta[2]*(1.5*dzeta_dt[3]+ddij_dt*zeta[2]
                + (d[i]*d[j]/(d[i]+d[j]))*dzeta_dt[2])/pow(1-zeta[3], 3.)
                + 6*pow((d[i]*d[j]/(d[i]+d[j]))*zeta[2], 2.)*dzeta_dt[3]/pow(1-zeta[3], 4.);
        }
    }

    double dadt_hs = 1/zeta[0]*(3*(dzeta_dt[1]*zeta[2] + zeta[1]*dzeta_dt[2])/(1-zeta[3])
        + 3*zeta[1]*zeta[2]*dzeta_dt[3]/pow(1-zeta[3], 2.)
        + 3*pow(zeta[2], 2.)*dzeta_dt[2]/zeta[3]/pow(1-zeta[3], 2.)
        + pow(zeta[2],3.)*dzeta_dt[3]*(3*zeta[3]-1)/pow(zeta[3], 2.)/pow(1-zeta[3], 3.)
        + (3*pow(zeta[2], 2.)*dzeta_dt[2]*zeta[3] - 2*pow(zeta[2], 3.)*dzeta_dt[3])/pow(zeta[3], 3.)
        * log(1-zeta[3])
        + (zeta[0]-pow(zeta[2],3)/pow(zeta[3],2.))*dzeta_dt[3]/(1-zeta[3]));

    static double a0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
    static double a1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
    static double a2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
    static double b0[7] = { 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
    static double b1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
    static double b2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };

    vector<double> a(7, 0);
    vector<double> b(7, 0);
    for (int i = 0; i < 7; i++) {
        a[i] = a0[i] + (m_avg - 1.) / m_avg * a1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * a2[i];
        b[i] = b0[i] + (m_avg - 1.) / m_avg * b1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * b2[i];
    }

    double I1 = 0.0;
    double I2 = 0.0;
    double dI1_dt = 0.0, dI2_dt = 0.;
    for (int i = 0; i < 7; i++) {
        I1 += a[i] * pow(eta, i);
        I2 += b[i] * pow(eta, i);
        dI1_dt += a[i] * dzeta_dt[3] * i * pow(eta, i - 1);
        dI2_dt += b[i] * dzeta_dt[3] * i * pow(eta, i - 1);
    }
    double C1 = 1.
                / (1. + m_avg * (8 * eta - 2 * eta * eta) / pow(1 - eta, 4)
                   + (1 - m_avg) * (20 * eta - 27 * eta * eta + 12 * pow(eta, 3) - 2 * pow(eta, 4)) / pow((1 - eta) * (2 - eta), 2.0));
    double C2 = -1 * C1 * C1
                * (m_avg * (-4 * eta * eta + 20 * eta + 8) / pow(1 - eta, 5.)
                   + (1 - m_avg) * (2 * pow(eta, 3) + 12 * eta * eta - 48 * eta + 40) / pow((1 - eta) * (2 - eta), 3));
    double dC1_dt = C2 * dzeta_dt[3];

    summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i] * (components[i].getM() - 1) * dghs_dt[i * ncomp + i] / ghs[i * ncomp + i];
    }

    double dadt_hc = m_avg * dadt_hs - summ;
    double dadt_disp = -2 * PI * den * (dI1_dt - I1 / _T) * m2es3 - PI * den * m_avg * (dC1_dt * I2 + C1 * dI2_dt - 2 * C1 * I2 / _T) * m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double dadt_polar = 0.;
    if (polar_term) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_dt = 0.;
        double dA3_dt = 0.;
        vector<double> dipmSQ(ncomp, 0);

        static double a0dip[5] = {0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308};
        static double a1dip[5] = {0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135};
        static double a2dip[5] = {-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575};
        static double b0dip[5] = {0.2187939, -1.1896431, 1.1626889, 0, 0};
        static double b1dip[5] = {-0.5873164, 1.2489132, -0.5085280, 0, 0};
        static double b2dip[5] = {3.4869576, -14.915974, 15.372022, 0, 0};
        static double c0dip[5] = {-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0};
        static double c1dip[5] = {-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0};
        static double c2dip[5] = {-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0};

        const static double conv = 7242.702976750923;  // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(components[i].getDipm(), 2.) / (components[i].getM() * components[i].getU() * pow(components[i].getSigma(), 3.)) * conv;
        }

        vector<double> adip(5, 0);
        vector<double> bdip(5, 0);
        vector<double> cdip(5, 0);
        double J2, J3, dJ2_dt, dJ3_dt;
        double m_ij;
        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(components[i].getM() * components[j].getM());
                if (m_ij > 2) {
                    m_ij = 2;
                }
                J2 = 0.;
                dJ2_dt = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij - 1) / m_ij * a1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij - 1) / m_ij * b1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * b2dip[l];
                    J2 += (adip[l] + bdip[l] * e_ij[i * ncomp + j] / _T) * pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_dt += adip[l] * l * pow(eta, l - 1) * dzeta_dt[3]
                        + bdip[l] * e_ij[j * ncomp + j] * (1 / _T * l * pow(eta, l - 1) * dzeta_dt[3]
                        - 1 / pow(_T, 2.) * pow(eta, l));
                }
                A2 += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                      * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i]
                      * dipmSQ[j] * J2;
                dA2_dt += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] * e_ij[j * ncomp + j] * pow(s_ij[i * ncomp + i], 3)
                          * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum()
                          * dipmSQ[i] * dipmSQ[j] * (dJ2_dt / pow(_T, 2) - 2 * J2 / pow(_T, 3));

                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((components[i].getM() * components[j].getM() * components[k].getM()), 1 / 3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    J3 = 0.;
                    dJ3_dt = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk - 1) / m_ijk * c1dip[l] + (m_ijk - 1) / m_ijk * (m_ijk - 2) / m_ijk * c2dip[l];
                        J3 += cdip[l] * pow(eta, l);
                        dJ3_dt += cdip[l] * l * pow(eta, l - 1) * dzeta_dt[3];
                    }
                    A3 += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                          * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                          / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum()
                          * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * J3;
                    dA3_dt += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] * e_ij[j * ncomp + j]
                              * e_ij[k * ncomp + k] * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                              / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum()
                              * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k]
                              * (-3 * J3 / pow(_T, 4) + dJ3_dt / pow(_T, 3));
                }
            }
        }

        A2 = -PI * den * A2;
        A3 = -4 / 3. * PI * PI * den * den * A3;
        dA2_dt = -PI * den * dA2_dt;
        dA3_dt = -4 / 3. * PI * PI * den * den * dA3_dt;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            dadt_polar = (dA2_dt - 2 * A3 / A2 * dA2_dt + dA3_dt) / pow(1 - A3 / A2, 2.);
        }
    }

    // Association term -------------------------------------------------------
    double dadt_assoc = 0.;
    if (assoc_term) {
        int num_sites = 0;
        vector<int> iA; //indices of associating compounds
        for(std::vector<int>::iterator it = assoc_num.begin(); it != assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(it - assoc_num.begin());
            }
        }

        vector<double> x_assoc(num_sites); // mole fractions of only the associating compounds
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        // these indices are necessary because we are only using 1D vectors
        vector<double> XA (num_sites, 0);
        vector<double> delta_ij(num_sites * num_sites, 0);
        vector<double> ddelta_dt(num_sites * num_sites, 0);
        int idxa = 0;
        int idxi = 0; // index for the ii-th compound
        int idxj = 0; // index for the jj-th compound
        for (int i = 0; i < num_sites; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                idxj = iA[j]*ncomp+iA[j];
                if (assoc_matrix[idxa] != 0) {
                    double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                    double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                        s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                    delta_ij[idxa] = ghs[iA[i]*ncomp+iA[j]]*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                    ddelta_dt[idxa] = pow(s_ij[idxj],3)*volABij*(-eABij/pow(_T,2)
                        *exp(eABij/_T)*ghs[iA[i]*ncomp+iA[j]] + dghs_dt[iA[i]*ncomp+iA[j]]
                        *(exp(eABij/_T)-1));
                }
                idxa += 1;
            }
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dt(num_sites, 0);
        dXA_dt = dXAdt_find(delta_ij, den, XA, ddelta_dt, x_assoc);

        for (int i = 0; i < num_sites; i++) {
            dadt_assoc += mole_fractions[iA[i]]*(1/XA[i]-0.5)*dXA_dt[i];
        }
    }

    // Ion term ---------------------------------------------------------------
    double dadt_ion = 0.;
    if (ion_term) {
        vector<double> q(ncomp);
        for (int i = 0; i < ncomp; i++) {
            q[i] = components[i].getZ() * E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += components[i].getZ() * components[i].getZ() * mole_fractions[i];
        }
        double kappa =
          sqrt(den * E_CHRG * E_CHRG / kb / _T / (dielc * perm_vac) * summ);  // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        double dkappa_dt;
        if (kappa != 0) {
            vector<double> chi(ncomp);
            vector<double> dchikap_dk(ncomp);
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi[i] = 3 / pow(kappa * components[i].getSigma(), 3)
                         * (1.5 + log(1 + kappa * components[i].getSigma()) - 2 * (1 + kappa * components[i].getSigma())
                            + 0.5 * pow(1 + kappa * components[i].getSigma(), 2));
                dchikap_dk[i] = -2 * chi[i] + 3 / (1 + kappa * components[i].getSigma());
                summ += mole_fractions[i] * components[i].getZ() * components[i].getZ();
            }
            dkappa_dt = -0.5 * den * E_CHRG * E_CHRG / kb / _T / _T / (dielc * perm_vac) * summ / kappa;

            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                summ += mole_fractions[i] * q[i] * q[i] * (dchikap_dk[i] * dkappa_dt / _T - kappa * chi[i] / _T / _T);
            }
            dadt_ion = -1 / 12. / PI / kb / (dielc * perm_vac) * summ;
        }
    }

    double dadt = dadt_hc + dadt_disp + dadt_assoc + dadt_polar + dadt_ion;
    return dadt;
}

CoolPropDbl PCSAFTBackend::calc_hmolar_residual(void) {
    CoolPropDbl Z = calc_compressibility_factor();
    CoolPropDbl dares_dt = calc_dadt();

    CoolPropDbl hres = (-_T * dares_dt + (Z - 1)) * kb * N_AV * _T;  // Equation A.46 from Gross and Sadowski 2001
    return hres;
}

CoolPropDbl PCSAFTBackend::calc_smolar_residual(void) {
    CoolPropDbl dares_dt = calc_dadt();
    CoolPropDbl ares = calc_alphar();

    CoolPropDbl sres = kb * N_AV * (-_T * dares_dt - ares);
    return sres;
}

vector<CoolPropDbl> PCSAFTBackend::calc_fugacity_coefficients(void) {
    int ncomp = N;  // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = components[i].getSigma() * (1 - 0.12 * exp(-3 * components[i].getU() / _T));
    }
    if (ion_term) {
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getZ() != 0) {
                d[i] =
                  components[i].getSigma() * (1 - 0.12);  // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }

    double den = _rhomolar * N_AV / 1.0e30;

    vector<double> zeta(4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j] * components[j].getM() * pow(d[j], i);
        }
        zeta[i] = PI / 6 * den * summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += mole_fractions[i] * components[i].getM();
    }

    vector<double> ghs(ncomp * ncomp, 0);
    vector<double> denghs(ncomp * ncomp, 0);
    vector<double> e_ij(ncomp * ncomp, 0);
    vector<double> s_ij(ncomp * ncomp, 0);
    double m2es3 = 0.;
    double m2e2s3 = 0.;
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
       for (int j = 0; j < ncomp; j++) {
           idx += 1;
           s_ij[idx] = (components[i].getSigma() + components[j].getSigma())/2.;
           if (ion_term) {
               if (components[i].getZ()*components[j].getZ() <= 0) { // for two cations or two anions e_ij is kept at zero to avoid dispersion between like ions (see Held et al. 2014)
                   if (k_ij.empty()) {
                       e_ij[idx] = sqrt(components[i].getU()*components[j].getU());
                   }
                   else {
                       e_ij[idx] = sqrt(components[i].getU()*components[j].getU())*(1 - (k_ij[idx] + k_ijT[idx] * _T));
                   }
               }
           } else {
               if (k_ij.empty()) {
                   e_ij[idx] = sqrt(components[i].getU()*components[j].getU());
               }
               else {
                   e_ij[idx] = sqrt(components[i].getU()*components[j].getU())*(1 - (k_ij[idx] + k_ijT[idx] * _T));
               }
           }
           m2es3 = m2es3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*e_ij[idx]/_T*pow(s_ij[idx], 3);
           m2e2s3 = m2e2s3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*pow(e_ij[idx]/_T,2)*pow(s_ij[idx], 3);
           ghs[idx] = 1/(1-zeta[3]) + (d[i]*d[j]/(d[i]+d[j]))*3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                   pow(d[i]*d[j]/(d[i]+d[j]), 2)*2*zeta[2]*zeta[2]/pow(1-zeta[3], 3);
           denghs[idx] = zeta[3]/(1-zeta[3])/(1-zeta[3]) +
               (d[i]*d[j]/(d[i]+d[j]))*(3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
               6*zeta[2]*zeta[3]/pow(1-zeta[3], 3)) +
               pow(d[i]*d[j]/(d[i]+d[j]), 2)*(4*zeta[2]*zeta[2]/pow(1-zeta[3], 3) +
               6*zeta[2]*zeta[2]*zeta[3]/pow(1-zeta[3], 4));
       }
    }

    double ares_hs = 1/zeta[0]*(3*zeta[1]*zeta[2]/(1-zeta[3]) + pow(zeta[2], 3.)/(zeta[3]*pow(1-zeta[3],2))
           + (pow(zeta[2], 3.)/pow(zeta[3], 2.) - zeta[0])*log(1-zeta[3]));
    double Zhs = zeta[3]/(1-zeta[3]) + 3.*zeta[1]*zeta[2]/zeta[0]/(1.-zeta[3])/(1.-zeta[3]) +
       (3.*pow(zeta[2], 3.) - zeta[3]*pow(zeta[2], 3.))/zeta[0]/pow(1.-zeta[3], 3.);

    static double a0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
    static double a1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
    static double a2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
    static double b0[7] = { 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
    static double b1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
    static double b2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };

    vector<double> a(7, 0);
    vector<double> b(7, 0);
    for (int i = 0; i < 7; i++) {
        a[i] = a0[i] + (m_avg - 1.) / m_avg * a1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * a2[i];
        b[i] = b0[i] + (m_avg - 1.) / m_avg * b1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * b2[i];
    }

    double detI1_det = 0.0;
    double detI2_det = 0.0;
    double I1 = 0.0;
    double I2 = 0.0;
    for (int i = 0; i < 7; i++) {
        detI1_det += a[i] * (i + 1) * pow(eta, i);
        detI2_det += b[i] * (i + 1) * pow(eta, i);
        I2 += b[i] * pow(eta, i);
        I1 += a[i] * pow(eta, i);
    }
    double C1 = 1.
                / (1. + m_avg * (8 * eta - 2 * eta * eta) / pow(1 - eta, 4)
                   + (1 - m_avg) * (20 * eta - 27 * eta * eta + 12 * pow(eta, 3) - 2 * pow(eta, 4)) / pow((1 - eta) * (2 - eta), 2.0));
    double C2 = -1. * C1 * C1
                * (m_avg * (-4 * eta * eta + 20 * eta + 8) / pow(1 - eta, 5)
                   + (1 - m_avg) * (2 * pow(eta, 3) + 12 * eta * eta - 48 * eta + 40) / pow((1 - eta) * (2 - eta), 3.0));

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i] * (components[i].getM() - 1) * log(ghs[i * ncomp + i]);
    }

    double ares_hc = m_avg * ares_hs - summ;
    double ares_disp = -2 * PI * den * I1 * m2es3 - PI * den * m_avg * C1 * I2 * m2e2s3;

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i] * (components[i].getM() - 1) / ghs[i * ncomp + i] * denghs[i * ncomp + i];
    }

    double Zhc = m_avg * Zhs - summ;
    double Zdisp = -2 * PI * den * detI1_det * m2es3 - PI * den * m_avg * (C1 * detI2_det + C2 * eta * I2) * m2e2s3;

    vector<double> dghsii_dx(ncomp * ncomp, 0);
    vector<double> dahs_dx(ncomp, 0);
    vector<double> dzeta_dx(4, 0);
    idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int l = 0; l < 4; l++) {
            dzeta_dx[l] = PI / 6. * den * components[i].getM() * pow(d[i], l);
        }
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            dghsii_dx[idx] =
              dzeta_dx[3] / (1 - zeta[3]) / (1 - zeta[3])
              + (d[j] * d[j] / (d[j] + d[j])) * (3 * dzeta_dx[2] / (1 - zeta[3]) / (1 - zeta[3]) + 6 * zeta[2] * dzeta_dx[3] / pow(1 - zeta[3], 3))
              + pow(d[j] * d[j] / (d[j] + d[j]), 2)
                  * (4 * zeta[2] * dzeta_dx[2] / pow(1 - zeta[3], 3) + 6 * zeta[2] * zeta[2] * dzeta_dx[3] / pow(1 - zeta[3], 4));
        }
        dahs_dx[i] =
          -dzeta_dx[0] / zeta[0] * ares_hs
          + 1 / zeta[0]
              * (3 * (dzeta_dx[1] * zeta[2] + zeta[1] * dzeta_dx[2]) / (1 - zeta[3])
                 + 3 * zeta[1] * zeta[2] * dzeta_dx[3] / (1 - zeta[3]) / (1 - zeta[3])
                 + 3 * zeta[2] * zeta[2] * dzeta_dx[2] / zeta[3] / (1 - zeta[3]) / (1 - zeta[3])
                 + pow(zeta[2], 3) * dzeta_dx[3] * (3 * zeta[3] - 1) / zeta[3] / zeta[3] / pow(1 - zeta[3], 3)
                 + log(1 - zeta[3])
                     * ((3 * zeta[2] * zeta[2] * dzeta_dx[2] * zeta[3] - 2 * pow(zeta[2], 3) * dzeta_dx[3]) / pow(zeta[3], 3) - dzeta_dx[0])
                 + (zeta[0] - pow(zeta[2], 3) / zeta[3] / zeta[3]) * dzeta_dx[3] / (1 - zeta[3]));
    }

    vector<double> dadisp_dx(ncomp, 0);
    vector<double> dahc_dx(ncomp, 0);
    double dzeta3_dx, daa_dx, db_dx, dI1_dx, dI2_dx, dm2es3_dx, dm2e2s3_dx, dC1_dx;
    for (int i = 0; i < ncomp; i++) {
        dzeta3_dx = PI / 6. * den * components[i].getM() * pow(d[i], 3);
        dI1_dx = 0.0;
        dI2_dx = 0.0;
        dm2es3_dx = 0.0;
        dm2e2s3_dx = 0.0;
        for (int l = 0; l < 7; l++) {
            daa_dx = components[i].getM() / m_avg / m_avg * a1[l] + components[i].getM() / m_avg / m_avg * (3 - 4 / m_avg) * a2[l];
            db_dx = components[i].getM() / m_avg / m_avg * b1[l] + components[i].getM() / m_avg / m_avg * (3 - 4 / m_avg) * b2[l];
            dI1_dx += a[l] * l * dzeta3_dx * pow(eta, l - 1) + daa_dx * pow(eta, l);
            dI2_dx += b[l] * l * dzeta3_dx * pow(eta, l - 1) + db_dx * pow(eta, l);
        }
        for (int j = 0; j < ncomp; j++) {
            dm2es3_dx += mole_fractions[j] * components[j].getM() * (e_ij[i * ncomp + j] / _T) * pow(s_ij[i * ncomp + j], 3);
            dm2e2s3_dx += mole_fractions[j] * components[j].getM() * pow(e_ij[i * ncomp + j] / _T, 2) * pow(s_ij[i * ncomp + j], 3);
            dahc_dx[i] += mole_fractions[j] * (components[j].getM() - 1) / ghs[j * ncomp + j] * dghsii_dx[i * ncomp + j];
        }
        dm2es3_dx = dm2es3_dx * 2 * components[i].getM();
        dm2e2s3_dx = dm2e2s3_dx * 2 * components[i].getM();
        dahc_dx[i] = components[i].getM() * ares_hs + m_avg * dahs_dx[i] - dahc_dx[i] - (components[i].getM() - 1) * log(ghs[i * ncomp + i]);
        dC1_dx = C2 * dzeta3_dx
                 - C1 * C1
                     * (components[i].getM() * (8 * eta - 2 * eta * eta) / pow(1 - eta, 4)
                        - components[i].getM() * (20 * eta - 27 * eta * eta + 12 * pow(eta, 3) - 2 * pow(eta, 4)) / pow((1 - eta) * (2 - eta), 2));

        dadisp_dx[i] =
          -2 * PI * den * (dI1_dx * m2es3 + I1 * dm2es3_dx)
          - PI * den * ((components[i].getM() * C1 * I2 + m_avg * dC1_dx * I2 + m_avg * C1 * dI2_dx) * m2e2s3 + m_avg * C1 * I2 * dm2e2s3_dx);
    }

    vector<double> mu_hc(ncomp, 0);
    vector<double> mu_disp(ncomp, 0);
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            mu_hc[i] += mole_fractions[j] * dahc_dx[j];
            mu_disp[i] += mole_fractions[j] * dadisp_dx[j];
        }
        mu_hc[i] = ares_hc + Zhc + dahc_dx[i] - mu_hc[i];
        mu_disp[i] = ares_disp + Zdisp + dadisp_dx[i] - mu_disp[i];
    }

    // Dipole term (Gross and Vrabec term) --------------------------------------
    vector<double> mu_polar(ncomp, 0);
    if (polar_term) {
       double A2 = 0.;
       double A3 = 0.;
       double dA2_det = 0.;
       double dA3_det = 0.;
       vector<double> dA2_dx(ncomp, 0);
       vector<double> dA3_dx(ncomp, 0);

       static double a0dip[5] = { 0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308 };
       static double a1dip[5] = { 0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135 };
       static double a2dip[5] = { -1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575 };
       static double b0dip[5] = { 0.2187939, -1.1896431, 1.1626889, 0, 0 };
       static double b1dip[5] = { -0.5873164, 1.2489132, -0.5085280, 0, 0 };
       static double b2dip[5] = { 3.4869576, -14.915974, 15.372022, 0, 0 };
       static double c0dip[5] = { -0.0646774, 0.1975882, -0.8087562, 0.6902849, 0 };
       static double c1dip[5] = { -0.9520876, 2.9924258, -2.3802636, -0.2701261, 0 };
       static double c2dip[5] = { -0.6260979, 1.2924686, 1.6542783, -3.4396744, 0 };

       const static double conv = 7242.702976750923; // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

       vector<double> dipmSQ (ncomp, 0);
       for (int i = 0; i < ncomp; i++) {
           dipmSQ[i] = pow(components[i].getDipm(), 2.)/(components[i].getM()*components[i].getU()*pow(components[i].getSigma(),3.))*conv;
       }

       vector<double> adip (5, 0);
       vector<double> bdip (5, 0);
       vector<double> cdip (5, 0);
       double J2, dJ2_det, detJ2_det, J3, dJ3_det, detJ3_det;
       double m_ij;
       double m_ijk;
       for (int i = 0; i < ncomp; i++) {
           for (int j = 0; j < ncomp; j++) {
               m_ij = sqrt(components[i].getM()*components[j].getM());
               if (m_ij > 2) {
                   m_ij = 2;
               }
               J2 = 0.;
               dJ2_det = 0.;
               detJ2_det = 0.;
               for (int l = 0; l < 5; l++) {
                   adip[l] = a0dip[l] + (m_ij-1)/m_ij*a1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*a2dip[l];
                   bdip[l] = b0dip[l] + (m_ij-1)/m_ij*b1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*b2dip[l];
                   J2 += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                   dJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*l*pow(eta, l-1);
                   detJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*(l+1)*pow(eta, l);
               }
               A2 += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                   pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*J2;
               dA2_det += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*
                   pow(s_ij[j*ncomp+j],3)/pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*detJ2_det;
               if (i == j) {
                   dA2_dx[i] += e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)
                       /pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*
                       (mole_fractions[i]*mole_fractions[j]*dJ2_det*PI/6.*den*components[i].getM()*pow(d[i],3) + 2*mole_fractions[j]*J2);
               }
               else {
                   dA2_dx[i] += e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)
                       /pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*
                       (mole_fractions[i]*mole_fractions[j]*dJ2_det*PI/6.*den*components[i].getM()*pow(d[i],3) + mole_fractions[j]*J2);
               }

               for (int k = 0; k < ncomp; k++) {
                   m_ijk = pow((components[i].getM()*components[j].getM()*components[k].getM()),1/3.);
                   if (m_ijk > 2) {
                       m_ijk = 2;
                   }
                   J3 = 0.;
                   dJ3_det = 0.;
                   detJ3_det = 0.;
                   for (int l = 0; l < 5; l++) {
                       cdip[l] = c0dip[l] + (m_ijk-1)/m_ijk*c1dip[l] + (m_ijk-1)/m_ijk*(m_ijk-2)/m_ijk*c2dip[l];
                       J3 += cdip[l]*pow(eta, l);
                       dJ3_det += cdip[l]*l*pow(eta, (l-1));
                       detJ3_det += cdip[l]*(l+2)*pow(eta, (l+1));
                   }
                   A3 += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*
                       pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                       s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                       dipmSQ[j]*dipmSQ[k]*J3;
                   dA3_det += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*
                       pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                       s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                       dipmSQ[j]*dipmSQ[k]*detJ3_det;
                   if ((i == j) && (i == k)) {
                       dA3_dx[i] += e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*pow(s_ij[i*ncomp+i],3)
                           *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                           *components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*dipmSQ[j]
                           *dipmSQ[k]*(mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*dJ3_det*PI/6.*den*components[i].getM()*pow(d[i],3)
                           + 3*mole_fractions[j]*mole_fractions[k]*J3);
                   }
                   else if ((i == j) || (i == k)) {
                       dA3_dx[i] += e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*pow(s_ij[i*ncomp+i],3)
                           *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                           *components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*dipmSQ[j]
                           *dipmSQ[k]*(mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*dJ3_det*PI/6.*den*components[i].getM()*pow(d[i],3)
                           + 2*mole_fractions[j]*mole_fractions[k]*J3);
                   }
                   else {
                       dA3_dx[i] += e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*pow(s_ij[i*ncomp+i],3)
                           *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                           *components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*dipmSQ[j]
                           *dipmSQ[k]*(mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*dJ3_det*PI/6.*den*components[i].getM()*pow(d[i],3)
                           + mole_fractions[j]*mole_fractions[k]*J3);
                   }
               }
           }
       }

       A2 = -PI*den*A2;
       A3 = -4/3.*PI*PI*den*den*A3;
       dA2_det = -PI*den/eta*dA2_det;
       dA3_det = -4/3.*PI*PI*den/eta*den/eta*dA3_det;
       for (int i = 0; i < ncomp; i++) {
           dA2_dx[i] = -PI*den*dA2_dx[i];
           dA3_dx[i] = -4/3.*PI*PI*den*den*dA3_dx[i];
       }

       vector<double> dapolar_dx(ncomp);
       for (int i = 0; i < ncomp; i++) {
           dapolar_dx[i] = (dA2_dx[i]*(1-A3/A2) + (dA3_dx[i]*A2 - A3*dA2_dx[i])/A2)/pow(1-A3/A2,2);
       }

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            double ares_polar = A2/(1-A3/A2);
            double Zpolar = eta*((dA2_det*(1-A3/A2)+(dA3_det*A2-A3*dA2_det)/A2)/(1-A3/A2)/(1-A3/A2));
            for (int i = 0; i < ncomp; i++) {
               for (int j = 0; j < ncomp; j++) {
                   mu_polar[i] += mole_fractions[j]*dapolar_dx[j];
               }
               mu_polar[i] = ares_polar + Zpolar + dapolar_dx[i] - mu_polar[i];
            }
        }
    }

    // Association term -------------------------------------------------------
    vector<double> mu_assoc(ncomp, 0);
    if (assoc_term) {
        int num_sites = 0;
        vector<int> iA; //indices of associating compounds
        for(std::vector<int>::iterator it = assoc_num.begin(); it != assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(it - assoc_num.begin());
            }
        }

        vector<double> x_assoc(num_sites); // mole fractions of only the associating compounds
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        // these indices are necessary because we are only using 1D vectors
        vector<double> XA (num_sites, 0);
        vector<double> delta_ij(num_sites * num_sites, 0);
        int idxa = 0;
        int idxi = 0; // index for the ii-th compound
        int idxj = 0; // index for the jj-th compound
        for (int i = 0; i < num_sites; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                idxj = iA[j]*ncomp+iA[j];
                if (assoc_matrix[idxa] != 0) {
                    double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                    double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                        s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                    delta_ij[idxa] = ghs[iA[i]*ncomp+iA[j]]*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                }
                idxa += 1;
            }
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; k++) {
            int idxi = 0; // index for the ii-th compound
            int idxj = 0; // index for the jj-th compound
            idxa = 0;
            for (int i = 0; i < num_sites; i++) {
                idxi = iA[i]*ncomp+iA[i];
                for (int j = 0; j < num_sites; j++) {
                    idxj = iA[j]*ncomp+iA[j];
                    if (assoc_matrix[idxa] != 0) {
                        double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                        double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                            s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                        double dghsd_dx = PI/6.*components[k].getM()*(pow(d[k], 3)/(1-zeta[3])/(1-zeta[3]) + 3*d[iA[i]]*d[iA[j]]/
                            (d[iA[i]]+d[iA[j]])*(d[k]*d[k]/(1-zeta[3])/(1-zeta[3])+2*pow(d[k], 3)*
                            zeta[2]/pow(1-zeta[3], 3)) + 2*pow((d[iA[i]]*d[iA[j]]/(d[iA[i]]+d[iA[j]])), 2)*
                            (2*d[k]*d[k]*zeta[2]/pow(1-zeta[3], 3)+3*(pow(d[k], 3)*zeta[2]*zeta[2]
                            /pow(1-zeta[3], 4))));
                        ddelta_dx[idx_ddelta] = dghsd_dx*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                    }
                    idx_ddelta += 1;
                    idxa += 1;
                }
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = dXAdx_find(assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

        int ij = 0;
        for (int i = 0; i < ncomp; i++) {
           for (int j = 0; j < num_sites; j++) {
               mu_assoc[i] += mole_fractions[iA[j]]*den*dXA_dx[ij]*(1/XA[j]-0.5);
               ij += 1;
           }
        }

        for (int i = 0; i < num_sites; i++) {
           mu_assoc[iA[i]] += log(XA[i]) - 0.5*XA[i] + 0.5;
        }
    }

    // Ion term ---------------------------------------------------------------
    vector<double> mu_ion(ncomp, 0);
    if (ion_term) {
        vector<double> q(ncomp);
        for (int i = 0; i < ncomp; i++) {
            q[i] = components[i].getZ() * E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += components[i].getZ() * components[i].getZ() * mole_fractions[i];
        }
        double kappa =
          sqrt(den * E_CHRG * E_CHRG / kb / _T / (dielc * perm_vac) * summ);  // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            vector<double> chi(ncomp);
            vector<double> sigma_k(ncomp);
            double summ1 = 0.;
            double summ2 = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi[i] = 3 / pow(kappa * components[i].getSigma(), 3)
                         * (1.5 + log(1 + kappa * components[i].getSigma()) - 2 * (1 + kappa * components[i].getSigma())
                            + 0.5 * pow(1 + kappa * components[i].getSigma(), 2));
                sigma_k[i] = -2 * chi[i] + 3 / (1 + kappa * components[i].getSigma());
                summ1 += q[i] * q[i] * mole_fractions[i] * sigma_k[i];
                summ2 += mole_fractions[i] * q[i] * q[i];
            }

            for (int i = 0; i < ncomp; i++) {
                mu_ion[i] = -q[i] * q[i] * kappa / 24. / PI / kb / _T / (dielc * perm_vac) * (2 * chi[i] + summ1 / summ2);
            }
        }
    }

    CoolPropDbl Z = calc_compressibility_factor();

    vector<double> mu(ncomp, 0);
    vector<CoolPropDbl> fugcoef(ncomp, 0);
    for (int i = 0; i < ncomp; i++) {
        mu[i] = mu_hc[i] + mu_disp[i] + mu_polar[i] + mu_assoc[i] + mu_ion[i];
        fugcoef[i] = exp(mu[i] - log(Z));  // the fugacity coefficients
    }

    return fugcoef;
}

CoolPropDbl PCSAFTBackend::calc_gibbsmolar_residual(void) {
    CoolPropDbl ares = calc_alphar();
    CoolPropDbl Z = calc_compressibility_factor();

    CoolPropDbl gres = (ares + (Z - 1) - log(Z)) * kb * N_AV * _T;  // Equation A.50 from Gross and Sadowski 2001
    return gres;
}

CoolPropDbl PCSAFTBackend::calc_compressibility_factor(void) {
    int ncomp = N;  // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = components[i].getSigma() * (1 - 0.12 * exp(-3 * components[i].getU() / _T));
    }
    if (ion_term) {
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getZ() != 0) {
                d[i] =
                  components[i].getSigma() * (1 - 0.12);  // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }

    double den = _rhomolar * N_AV / 1.0e30;

    vector<double> zeta(4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j] * components[j].getM() * pow(d[j], i);
        }
        zeta[i] = PI / 6 * den * summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += mole_fractions[i] * components[i].getM();
    }

    vector<double> ghs (ncomp*ncomp, 0);
    vector<double> denghs (ncomp*ncomp, 0);
    vector<double> e_ij (ncomp*ncomp, 0);
    vector<double> s_ij (ncomp*ncomp, 0);
    double m2es3 = 0.;
    double m2e2s3 = 0.;
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            s_ij[idx] = (components[i].getSigma() + components[j].getSigma()) / 2.;
            if (ion_term) {
                if (components[i].getZ() * components[j].getZ()
                    <= 0) {  // for two cations or two anions e_ij is kept at zero to avoid dispersion between like ions (see Held et al. 2014)
                    if (k_ij.empty()) {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                    } else {
                        e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                    }
                }
            } else {
                if (k_ij.empty()) {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU());
                } else {
                    e_ij[idx] = sqrt(components[i].getU() * components[j].getU()) * (1 - (k_ij[idx] + k_ijT[idx] * _T));
                }
            }
            m2es3 = m2es3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*e_ij[idx]/_T*pow(s_ij[idx], 3);
            m2e2s3 = m2e2s3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*pow(e_ij[idx]/_T,2)*pow(s_ij[idx], 3);
            ghs[idx] = 1/(1-zeta[3]) + (d[i]*d[j]/(d[i]+d[j]))*3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                    pow(d[i]*d[j]/(d[i]+d[j]), 2)*2*zeta[2]*zeta[2]/pow(1-zeta[3], 3);
            denghs[idx] = zeta[3]/(1-zeta[3])/(1-zeta[3]) +
                (d[i]*d[j]/(d[i]+d[j]))*(3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                6*zeta[2]*zeta[3]/pow(1-zeta[3], 3)) +
                pow(d[i]*d[j]/(d[i]+d[j]), 2)*(4*zeta[2]*zeta[2]/pow(1-zeta[3], 3) +
                6*zeta[2]*zeta[2]*zeta[3]/pow(1-zeta[3], 4));
        }
    }

    double Zhs = zeta[3] / (1 - zeta[3]) + 3. * zeta[1] * zeta[2] / zeta[0] / (1. - zeta[3]) / (1. - zeta[3])
                 + (3. * pow(zeta[2], 3.) - zeta[3] * pow(zeta[2], 3.)) / zeta[0] / pow(1. - zeta[3], 3.);

    static double a0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
    static double a1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
    static double a2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
    static double b0[7] = { 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
    static double b1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
    static double b2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };

    vector<double> a(7, 0);
    vector<double> b(7, 0);
    for (int i = 0; i < 7; i++) {
        a[i] = a0[i] + (m_avg - 1.) / m_avg * a1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * a2[i];
        b[i] = b0[i] + (m_avg - 1.) / m_avg * b1[i] + (m_avg - 1.) / m_avg * (m_avg - 2.) / m_avg * b2[i];
    }

    double detI1_det = 0.0;
    double detI2_det = 0.0;
    double I2 = 0.0;
    for (int i = 0; i < 7; i++) {
        detI1_det += a[i] * (i + 1) * pow(eta, i);
        detI2_det += b[i] * (i + 1) * pow(eta, i);
        I2 += b[i] * pow(eta, i);
    }
    double C1 = 1.
                / (1. + m_avg * (8 * eta - 2 * eta * eta) / pow(1 - eta, 4)
                   + (1 - m_avg) * (20 * eta - 27 * eta * eta + 12 * pow(eta, 3) - 2 * pow(eta, 4)) / pow((1 - eta) * (2 - eta), 2.0));
    double C2 = -1. * C1 * C1
                * (m_avg * (-4 * eta * eta + 20 * eta + 8) / pow(1 - eta, 5)
                   + (1 - m_avg) * (2 * pow(eta, 3) + 12 * eta * eta - 48 * eta + 40) / pow((1 - eta) * (2 - eta), 3.0));

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i]*(components[i].getM()-1)/ghs[i*ncomp + i]*denghs[i*ncomp + i];
    }

    double Zid = 1.0;
    double Zhc = m_avg * Zhs - summ;
    double Zdisp = -2 * PI * den * detI1_det * m2es3 - PI * den * m_avg * (C1 * detI2_det + C2 * eta * I2) * m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double Zpolar = 0;
    if (polar_term) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_det = 0.;
        double dA3_det = 0.;
        vector<double> adip(5, 0);
        vector<double> bdip(5, 0);
        vector<double> cdip(5, 0);
        vector<double> dipmSQ(ncomp, 0);
        double J2, detJ2_det, J3, detJ3_det;

        static double a0dip[5] = {0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308};
        static double a1dip[5] = {0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135};
        static double a2dip[5] = {-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575};
        static double b0dip[5] = {0.2187939, -1.1896431, 1.1626889, 0, 0};
        static double b1dip[5] = {-0.5873164, 1.2489132, -0.5085280, 0, 0};
        static double b2dip[5] = {3.4869576, -14.915974, 15.372022, 0, 0};
        static double c0dip[5] = {-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0};
        static double c1dip[5] = {-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0};
        static double c2dip[5] = {-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0};

        const static double conv = 7242.702976750923;  // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(components[i].getDipm(), 2.) / (components[i].getM() * components[i].getU() * pow(components[i].getSigma(), 3.)) * conv;
        }

        double m_ij;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(components[i].getM() * components[j].getM());
                if (m_ij > 2) {
                    m_ij = 2;
                }
                J2 = 0.;
                detJ2_det = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij-1)/m_ij*a1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij-1)/m_ij*b1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*b2dip[l];
                    J2 += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    detJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/_T)*(l+1)*pow(eta, l);
                }
                A2 += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*J2;
                dA2_det += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*pow(s_ij[i*ncomp+i],3)*
                    pow(s_ij[j*ncomp+j],3)/pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*detJ2_det;
            }
        }

        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((components[i].getM() * components[j].getM() * components[k].getM()), 1 / 3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    J3 = 0.;
                    detJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk-1)/m_ijk*c1dip[l] + (m_ijk-1)/m_ijk*(m_ijk-2)/m_ijk*c2dip[l];
                        J3 += cdip[l]*pow(eta, l);
                        detJ3_det += cdip[l]*(l+2)*pow(eta, (l+1));
                    }
                    A3 += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                    dA3_det += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/_T*e_ij[j*ncomp+j]/_T*e_ij[k*ncomp+k]/_T*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*detJ3_det;
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        dA2_det = -PI*den/eta*dA2_det;
        dA3_det = -4/3.*PI*PI*den/eta*den/eta*dA3_det;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
              Zpolar = eta*((dA2_det*(1-A3/A2)+(dA3_det*A2-A3*dA2_det)/A2)/(1-A3/A2)/(1-A3/A2));
        }
    }

    // Association term -------------------------------------------------------
    double Zassoc = 0;
    if (assoc_term) {
        int num_sites = 0;
        vector<int> iA; //indices of associating compounds
        for(std::vector<int>::iterator it = assoc_num.begin(); it != assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(it - assoc_num.begin());
            }
        }

        vector<double> x_assoc(num_sites); // mole fractions of only the associating compounds
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        // these indices are necessary because we are only using 1D vectors
        vector<double> XA (num_sites, 0);
        vector<double> delta_ij(num_sites * num_sites, 0);
        int idxa = 0;
        int idxi = 0; // index for the ii-th compound
        int idxj = 0; // index for the jj-th compound
        for (int i = 0; i < num_sites; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                idxj = iA[j]*ncomp+iA[j];
                if (assoc_matrix[idxa] != 0) {
                    double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                    double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                        s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                    delta_ij[idxa] = ghs[iA[i]*ncomp+iA[j]]*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                }
                idxa += 1;
            }
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; k++) {
            int idxi = 0; // index for the ii-th compound
            int idxj = 0; // index for the jj-th compound
            idxa = 0;
            for (int i = 0; i < num_sites; i++) {
                idxi = iA[i]*ncomp+iA[i];
                for (int j = 0; j < num_sites; j++) {
                    idxj = iA[j]*ncomp+iA[j];
                    if (assoc_matrix[idxa] != 0) {
                        double eABij = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                        double volABij = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                            s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                        double dghsd_dx = PI/6.*components[k].getM()*(pow(d[k], 3)/(1-zeta[3])/(1-zeta[3]) + 3*d[iA[i]]*d[iA[j]]/
                            (d[iA[i]]+d[iA[j]])*(d[k]*d[k]/(1-zeta[3])/(1-zeta[3])+2*pow(d[k], 3)*
                            zeta[2]/pow(1-zeta[3], 3)) + 2*pow((d[iA[i]]*d[iA[j]]/(d[iA[i]]+d[iA[j]])), 2)*
                            (2*d[k]*d[k]*zeta[2]/pow(1-zeta[3], 3)+3*(pow(d[k], 3)*zeta[2]*zeta[2]
                            /pow(1-zeta[3], 4))));
                        ddelta_dx[idx_ddelta] = dghsd_dx*(exp(eABij/_T)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                    }
                    idx_ddelta += 1;
                    idxa += 1;
                }
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-14)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = dXAdx_find(assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

        summ = 0.;
        int ij = 0;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < num_sites; j++) {
                summ += mole_fractions[i]*den*mole_fractions[iA[j]]*(1/XA[j]-0.5)*dXA_dx[ij];
                ij += 1;
            }
        }

        Zassoc = summ;
    }

    // Ion term ---------------------------------------------------------------
    double Zion = 0;
    if (ion_term) {
        vector<double> q(ncomp);
        for (int i = 0; i < ncomp; i++) {
            q[i] = components[i].getZ() * E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += pow(components[i].getZ(), 2.) * mole_fractions[i];
        }

        double kappa =
          sqrt(den * E_CHRG * E_CHRG / kb / _T / (dielc * perm_vac) * summ);  // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            double chi, sigma_k;
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi = 3 / pow(kappa * components[i].getSigma(), 3)
                      * (1.5 + log(1 + kappa * components[i].getSigma()) - 2 * (1 + kappa * components[i].getSigma())
                         + 0.5 * pow(1 + kappa * components[i].getSigma(), 2));
                sigma_k = -2 * chi + 3 / (1 + kappa * components[i].getSigma());
                summ += q[i] * q[i] * mole_fractions[i] * sigma_k;
            }
            Zion = -1 * kappa / 24. / PI / kb / _T / (dielc * perm_vac) * summ;
        }
    }

    double Z = Zid + Zhc + Zdisp + Zpolar + Zassoc + Zion;
    return Z;
}

void PCSAFTBackend::post_update(bool optional_checks) {
    // Check the values that must always be set
    if (!ValidNumber(_p)) {
        throw ValueError("p is not a valid number");
    }
    if (_T < 0) {
        throw ValueError("T is less than zero");
    }
    if (!ValidNumber(_T)) {
        throw ValueError("T is not a valid number");
    }
    if (_rhomolar < 0) {
        throw ValueError("rhomolar is less than zero");
    }
    if (!ValidNumber(_rhomolar)) {
        throw ValueError("rhomolar is not a valid number");
    }

    if (optional_checks) {
        if (!ValidNumber(_Q)) {
            throw ValueError("Q is not a valid number");
        }
        if (_phase == iphase_unknown) {
            throw ValueError("_phase is unknown");
        }
    }
}

void PCSAFTBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    if (get_debug_level() > 10) {
        std::cout << format("%s (%d): update called with (%d: (%s), %g, %g)", __FILE__, __LINE__, input_pair,
                            get_input_pair_short_desc(input_pair).c_str(), value1, value2)
                  << std::endl;
    }

    // Converting input to CoolPropDbl
    CoolPropDbl ld_value1 = value1, ld_value2 = value2;
    value1 = ld_value1;
    value2 = ld_value2;

    // Clear the state
    clear();

    if (is_pure_or_pseudopure == false && mole_fractions.size() == 0) {
        throw ValueError("Mole fractions must be set");
    }

    if (SatL->mole_fractions.empty()) {
        SatL->set_mole_fractions(mole_fractions);
    }
    if (SatV->mole_fractions.empty()) {
        SatV->set_mole_fractions(mole_fractions);
        double summ = 0;
        for (int i = 0; i < N; i++) {
            if (SatV->components[i].getZ() != 0) { // we make the assumption that ions do not appear in the vapor phase
                SatV->mole_fractions[i] = 0;
            }
            else {
                summ += SatV->mole_fractions[i];
            }
        }
        for (int i = 0; i < N; i++) {
            SatV->mole_fractions[i] = SatV->mole_fractions[i] / summ;
        }
    }

    // If the inputs are in mass units, convert them to molar units
    mass_to_molar_inputs(input_pair, value1, value2);

    switch (input_pair) {
        case PT_INPUTS:
            _p = value1;
            _T = value2;
            if (water_present) {
                components[water_idx].calc_water_sigma(_T);
                dielc = dielc_water(
                  _T);  // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
            }

            if (imposed_phase_index != iphase_not_imposed) {
                // Use the imposed phase index
                _phase = imposed_phase_index;
            } else {
                _phase = calc_phase_internal(input_pair);
            }
            _rhomolar = solver_rho_Tp(value2 /*T*/, value1 /*p*/, _phase /*phase*/);
            break;
        case QT_INPUTS:
            _Q = value1;
            _T = value2;
            SatL->_Q = value1;
            SatV->_Q = value1;
            SatL->_T = value2;
            SatV->_T = value2;
            _phase = iphase_twophase;
            if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            if (water_present) {
                components[water_idx].calc_water_sigma(_T);
                SatL->components[water_idx].calc_water_sigma(_T);
                SatV->components[water_idx].calc_water_sigma(_T);
                dielc = dielc_water(
                  _T);  // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
                SatL->dielc = dielc_water(_T);
                SatV->dielc = dielc_water(_T);
            }
            flash_QT(*this);
            break;
        case PQ_INPUTS:
            _p = value1;
            _Q = value2;
            SatL->_p = value1;
            SatV->_p = value1;
            SatL->_Q = value2;
            SatV->_Q = value2;
            _phase = iphase_twophase;
            if ((_Q < 0) || (_Q > 1)) {
                throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
            }
            flash_PQ(*this);
            break;
        case DmolarT_INPUTS:
            _rhomolar = value1; _T = value2;
            SatL->_rhomolar = value1; SatV->_rhomolar = value1;
            SatL->_T = value2; SatV->_T = value2;
            if (water_present) {
                components[water_idx].calc_water_sigma(_T);
                SatL->components[water_idx].calc_water_sigma(_T);
                SatV->components[water_idx].calc_water_sigma(_T);
                dielc = dielc_water(_T); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
                SatL->dielc = dielc_water(_T);
                SatV->dielc = dielc_water(_T);
            }
            _p = update_DmolarT(_rhomolar);

            if (imposed_phase_index != iphase_not_imposed) {
                // Use the imposed phase index
                _phase = imposed_phase_index;
            } else {
                _phase =
                  calc_phase_internal(input_pair);  // if in the two phase region, the pressure is updated by this function to equal the bubble point
            }
            break;
        case SmolarT_INPUTS:
        case DmolarP_INPUTS:
        case DmolarHmolar_INPUTS:
        case DmolarSmolar_INPUTS:
        case DmolarUmolar_INPUTS:
        case HmolarP_INPUTS:
        case PSmolar_INPUTS:
        case PUmolar_INPUTS:
        case HmolarSmolar_INPUTS:
        case QSmolar_INPUTS:
        case HmolarQ_INPUTS:
        case DmolarQ_INPUTS:
        default:
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
    }

    // set Q, if not already set
    if (!ValidNumber(_Q)) {
        if (_phase == iphase_gas) {
            _Q = 1;
        } else if (_phase == iphase_liquid) {
            _Q = 0;
        }
    }

    post_update();
}

phases PCSAFTBackend::calc_phase_internal(CoolProp::input_pairs input_pair) {
    phases phase = iphase_unknown;

    double p_input, rho_input;
    double p_bub, p_dew, p_equil;
    switch(input_pair)
    {
        case PT_INPUTS:
            p_input = _p; rho_input = _rhomolar;
            // first try to estimate without a full flash calculation
            _Q = 0;
            SatL->_Q = _Q; SatV->_Q = _Q;
            SatL->_T = _T; SatV->_T = _T;
            p_equil = estimate_flash_p(*this);
            if (p_input > 1.6 * p_equil) {
                phase = iphase_liquid;
            }
            else if (p_input < 0.5 * p_equil) {
                phase = iphase_gas;
            }
            else {
                // if the pressure is too close to the estimated bubble point, then do a full flash calculation to determine the phase
                _Q = 0;
                SatL->_Q = _Q; SatV->_Q = _Q;
                SatL->_T = _T; SatV->_T = _T;
                try {
                    flash_QT(*this);
                }
                catch (const SolutionError& ex) {
                    phase = iphase_supercritical;
                    break;
                }
                p_bub = _p;
                _p = p_input; _rhomolar = rho_input;
                if (_p > p_bub) {
                    phase = iphase_liquid;
                }
                else if (_p == p_bub) {
                    phase = iphase_twophase;
                }
                else {
                    _Q = 1;
                    SatL->_Q = _Q; SatV->_Q = _Q;
                    flash_QT(*this);
                    p_dew = _p;
                    _p = p_input; _rhomolar = rho_input;
                    if (_p < p_dew) {
                        phase = iphase_gas;
                    }
                    else if ((_p <= p_bub) && (_p >= p_dew)) {
                        phase = iphase_twophase;
                    }
                    else{
                    	phase = iphase_unknown;
                    }
                }
            }
            break;
        case DmolarT_INPUTS:
            double rho_bub, rho_dew;
            p_input = _p; rho_input = _rhomolar;

            _Q = 0;
            SatL->_Q = _Q;
            SatV->_Q = _Q;
            SatL->_T = _T;
            SatV->_T = _T;
            try {
                flash_QT(*this);
            } catch (const SolutionError& ex) {
                phase = iphase_supercritical;
                break;
            }
            rho_bub = _rhomolar;
            p_bub = _p;
            _p = p_input;
            _rhomolar = rho_input;
            if (_rhomolar > rho_bub) {
                phase = iphase_liquid;
            } else if (_rhomolar == rho_bub) {
                phase = iphase_twophase;
                _p = p_bub;
                _Q = 1 - (_rhomolar - SatV->_rhomolar) / (SatL->_rhomolar - SatV->_rhomolar);
            } else {
                _Q = 1;
                SatL->_Q = _Q;
                SatV->_Q = _Q;
                flash_QT(*this);
                rho_dew = _rhomolar;
                _p = p_input;
                _rhomolar = rho_input;
                if (_rhomolar < rho_dew) {
                    phase = iphase_gas;
                } else if ((_rhomolar <= rho_bub) && (_rhomolar >= rho_dew)) {
                    phase = iphase_twophase;
                    _p = p_bub;
                    _Q = 1 - (_rhomolar - SatV->_rhomolar) / (SatL->_rhomolar - SatV->_rhomolar);
                }
            }
            break;
        default:
            throw ValueError(
              format("Phase determination for this pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
    }

    return phase;
}


void PCSAFTBackend::flash_QT(PCSAFTBackend &PCSAFT) {
    bool solution_found = false;
    double p_guess = 0;
    double p = 0;
    try {
        p_guess = estimate_flash_p(PCSAFT);
        p = outerTQ(p_guess, PCSAFT);
        solution_found = true;
    }
    catch (const SolutionError& ex) {}
    catch (const ValueError& ex) {}

    // if solution hasn't been found, try cycling through a range of pressures
    if (!solution_found) {
        double p_lbound = -6; // here we're using log10 of the pressure
        double p_ubound = 9;
        double p_step = 0.1;
        p_guess = p_lbound;
        while (p_guess < p_ubound && !solution_found) {
            try {
                p = outerTQ(pow(10, p_guess), PCSAFT);
                solution_found = true;
            } catch (const SolutionError& ex) {
                p_guess += p_step;
            } catch (const ValueError& ex) {
                p_guess += p_step;
            }
        }
    }

    if (!solution_found) {
        throw SolutionError("solution could not be found for TQ flash");
    }

    // Load the outputs
    PCSAFT._p = p;
    PCSAFT._rhomolar = 1/(PCSAFT._Q/PCSAFT.SatV->_rhomolar + (1 - PCSAFT._Q)/PCSAFT.SatL->_rhomolar);
    PCSAFT._phase = iphase_twophase;
}


void PCSAFTBackend::flash_PQ(PCSAFTBackend &PCSAFT) {
    bool solution_found = false;
    double t_guess = 0;
    double t = 0;
    try {
        t_guess = estimate_flash_t(PCSAFT);
        t = outerPQ(t_guess, PCSAFT);
        solution_found = true;
    }
    catch (const SolutionError& ex) {}
    catch (const ValueError& ex) {}

    // if solution hasn't been found, try calling the flash function directly with a range of initial temperatures
    if (!solution_found) {
        double t_lbound = 1;
        double t_ubound = 800;
        double t_step = 10;
        if (PCSAFT.ion_term) {
            t_lbound = 264;
            t_ubound = 350;
        }
        t_guess = t_ubound;
        while (t_guess > t_lbound && !solution_found) {
            try {
                t = outerPQ(t_guess, PCSAFT);
                solution_found = true;
            } catch (const SolutionError& ex) {
                t_guess -= t_step;
            } catch (const ValueError& ex) {
                t_guess -= t_step;
            }
        }
    }

    if (!solution_found) {
        throw SolutionError("solution could not be found for PQ flash");
    }

    // Load the outputs
    PCSAFT._T = t;
    PCSAFT._rhomolar = 1/(PCSAFT._Q/PCSAFT.SatV->_rhomolar + (1 - PCSAFT._Q)/PCSAFT.SatL->_rhomolar);
    PCSAFT._phase = iphase_twophase;
}


double PCSAFTBackend::outerPQ(double t_guess, PCSAFTBackend &PCSAFT) {
    // Based on the algorithm proposed in H. A. J. Watson, M. Vikse, T. Gundersen, and P. I. Barton, “Reliable Flash Calculations: Part 1. Nonsmooth Inside-Out Algorithms,” Ind. Eng. Chem. Res., vol. 56, no. 4, pp. 960–973, Feb. 2017, doi: 10.1021/acs.iecr.6b03956.
    int ncomp = N; // number of components
    double TOL = 1e-8;
    double MAXITER = 200;

    // Define the residual to be driven to zero
    class SolverInnerResid : public FuncWrapper1D
    {
    public:
        PCSAFTBackend &PCSAFT;
        CoolPropDbl kb0;
        vector<CoolPropDbl> u;

        SolverInnerResid(PCSAFTBackend &PCSAFT, CoolPropDbl kb0, vector<CoolPropDbl> u)
        : PCSAFT(PCSAFT), kb0(kb0), u(u){}
        CoolPropDbl call(CoolPropDbl R){
            int ncomp = PCSAFT.components.size();
            double error = 0;

            vector<double> pp(ncomp, 0);
            double L = 0;
            for (int i = 0; i < ncomp; i++) {
                if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                    pp[i] = PCSAFT.mole_fractions[i] / (1 - R + kb0 * R * exp(u[i]));
                    L += pp[i];
                } else {
                    L += PCSAFT.mole_fractions[i];
                }
            }
            L = (1 - R) * L;

            error = pow((L + PCSAFT._Q - 1), 2.);
            return error;
        };
    };

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (PCSAFT.ion_term && PCSAFT.components[i].getZ() != 0) {
            x_ions += PCSAFT.mole_fractions[i];
        }
    }

    // initialize variables
    vector<double> k(ncomp, 0), u(ncomp, 0), kprime(ncomp, 0), uprime(ncomp, 0);
    double Tref = t_guess - 1;
    double Tprime = t_guess + 1;
    double t = t_guess;

    PCSAFT.SatL->_T = t; // _T must be updated because the density calculation depends on it
    PCSAFT.SatV->_T = t;

    // calculate sigma for water, if it is present
    if (PCSAFT.water_present) {
        PCSAFT.components[water_idx].calc_water_sigma(t);
        PCSAFT.SatL->components[water_idx].calc_water_sigma(t);
        PCSAFT.SatV->components[water_idx].calc_water_sigma(t);
        PCSAFT.dielc = dielc_water(t); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
        PCSAFT.SatL->dielc = dielc_water(t);
        PCSAFT.SatV->dielc = dielc_water(t);
    }

    // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(t, PCSAFT.SatL->_p, iphase_liquid);
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(t, PCSAFT.SatV->_p, iphase_gas);
    if ((PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-4) {
        throw SolutionError("liquid and vapor densities are the same.");
    }
    vector<double> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    vector<double> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();

    double xv_sum = 0;
    double xl_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) { // this if statement sets k to 0 for ionic components
            k[i] = fugcoef_l[i] / fugcoef_v[i];
        } else {
            k[i] = 0;
        }
        PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
        xl_sum += PCSAFT.SatL->mole_fractions[i];
        PCSAFT.SatV->mole_fractions[i] = k[i] * PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
        xv_sum += PCSAFT.SatV->mole_fractions[i];
    }

    if (xv_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            PCSAFT.SatV->mole_fractions[i] = PCSAFT.SatV->mole_fractions[i] / xv_sum;
        }
    }

    if (xl_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            PCSAFT.SatL->mole_fractions[i] = PCSAFT.SatL->mole_fractions[i] / xl_sum;
        }
    }

    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(t, PCSAFT.SatL->_p, iphase_liquid);
    fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(t, PCSAFT.SatV->_p, iphase_gas);
    fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
    for (int i = 0; i < ncomp; i++) {
        k[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    PCSAFT.SatL->_T = Tprime; // _T must be updated because the density calculation depends on it
    PCSAFT.SatV->_T = Tprime;

    if (PCSAFT.water_present) {
        PCSAFT.components[water_idx].calc_water_sigma(Tprime);
        PCSAFT.SatL->components[water_idx].calc_water_sigma(Tprime);
        PCSAFT.SatV->components[water_idx].calc_water_sigma(Tprime);
        PCSAFT.dielc = dielc_water(Tprime); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
        PCSAFT.SatL->dielc = dielc_water(Tprime);
        PCSAFT.SatV->dielc = dielc_water(Tprime);
    }
    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(Tprime, PCSAFT.SatL->_p, iphase_liquid);
    fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(Tprime, PCSAFT.SatV->_p, iphase_gas);
    fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
    for (int i = 0; i < ncomp; i++) {
        kprime[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    vector<double> t_weight(ncomp);
    double t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Tprime - t);
        t_weight[i] = PCSAFT.SatV->mole_fractions[i] * dlnk_dt / (1 + PCSAFT._Q * (k[i] - 1));
        t_sum += t_weight[i];
    }

    double kb = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
            kb += wi * std::log(k[i]);
        }
    }
    kb = std::exp(kb);

    t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Tprime - t);
        t_weight[i] = PCSAFT.SatV->mole_fractions[i] * dlnk_dt / (1 + PCSAFT._Q * (kprime[i] - 1));
        t_sum += t_weight[i];
    }

    double kbprime = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
            kbprime += wi * std::log(kprime[i]);
        }
    }
    kbprime = std::exp(kbprime);
    double kb0 = kbprime;

    for (int i = 0; i < ncomp; i++) {
        u[i] = std::log(k[i] / kb);
        uprime[i] = std::log(kprime[i] / kbprime);
    }

    double B = std::log(kbprime / kb) / (1/Tprime - 1/t);
    double A = std::log(kb) - B * (1/t - 1/Tref);

    // solve
    SolverInnerResid resid(*this, kb0, u);

    vector<double> pp(ncomp, 0);
    double maxdif = 1e10 * TOL;
    int itr = 0;
    double Rmin = 0, Rmax = 1;
    while (maxdif > TOL && itr < MAXITER) {
        // save previous values for calculating the difference at the end of the iteration
        vector<double> u_old = u;
        double A_old = A;

        resid.u = u;
        double R0 = kb * PCSAFT._Q / (kb * PCSAFT._Q + kb0 * (1 - PCSAFT._Q));
        double R = R0;
        if (resid.call(R) > TOL) {
            R = BoundedSecant(resid, R0, Rmin, Rmax, DBL_EPSILON, TOL, MAXITER);
        }

        double pp_sum = 0;
        double eupp_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            pp[i] = PCSAFT.mole_fractions[i] / (1 - R + kb0 * R * std::exp(u[i]));
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                pp_sum += pp[i];
                eupp_sum += std::exp(u[i]) * pp[i];
            }
        }
        kb = pp_sum / eupp_sum;

        t = 1 / (1 / Tref + (std::log(kb) - A) / B);
        for (int i = 0; i < ncomp; i++) {
            if (x_ions == 0) {
                PCSAFT.SatL->mole_fractions[i] = pp[i] / pp_sum;
                PCSAFT.SatV->mole_fractions[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                PCSAFT.SatL->mole_fractions[i] = pp[i] / pp_sum * (1 - x_ions / (1 - PCSAFT._Q));
                PCSAFT.SatV->mole_fractions[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else {
                PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 - PCSAFT._Q);
                PCSAFT.SatV->mole_fractions[i] = 0;
            }
        }

        PCSAFT.SatL->_T = t; // _T must be updated because the density calculation depends on it
        PCSAFT.SatV->_T = t;

        if (PCSAFT.water_present) {
            PCSAFT.components[water_idx].calc_water_sigma(t);
            PCSAFT.SatL->components[water_idx].calc_water_sigma(t);
            PCSAFT.SatV->components[water_idx].calc_water_sigma(t);
            PCSAFT.dielc = dielc_water(t); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
            PCSAFT.SatL->dielc = dielc_water(t);
            PCSAFT.SatV->dielc = dielc_water(t);
        }
        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(t, PCSAFT._p, iphase_liquid);
        vector<CoolPropDbl> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(t, PCSAFT._p, iphase_gas);
        vector<CoolPropDbl> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
        for (int i = 0; i < ncomp; i++) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
            u[i] = std::log(k[i] / kb);
        }

        if (itr == 0) {
            B = std::log(kbprime / kb) / (1/Tprime - 1/t);
            if (B > 0) {
                throw SolutionError("B > 0 in outerPQ");
            }
        }
        A = std::log(kb) - B * (1/t - 1/Tref);

        maxdif = std::abs(A - A_old);
        for (int i = 0; i < ncomp; i++) {
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                double dif = std::abs(u[i] - u_old[i]);
                if (dif > maxdif) {
                    maxdif = dif;
                }
            }
        }

        itr += 1;
    }

    if (!std::isfinite(t) || maxdif > 1e-3 || t < 0) {
        throw SolutionError("outerPQ did not converge to a solution");
    }

    return t;
}


double PCSAFTBackend::outerTQ(double p_guess, PCSAFTBackend &PCSAFT) {
    // Based on the algorithm proposed in H. A. J. Watson, M. Vikse, T. Gundersen, and P. I. Barton, “Reliable Flash Calculations: Part 1. Nonsmooth Inside-Out Algorithms,” Ind. Eng. Chem. Res., vol. 56, no. 4, pp. 960–973, Feb. 2017, doi: 10.1021/acs.iecr.6b03956.
    int ncomp = N; // number of components
    double TOL = 1e-8;
    double MAXITER = 200;

    // Define the residual to be driven to zero
    class SolverInnerResid : public FuncWrapper1D
    {
    public:
        PCSAFTBackend &PCSAFT;
        CoolPropDbl kb0;
        vector<CoolPropDbl> u;

        SolverInnerResid(PCSAFTBackend &PCSAFT, CoolPropDbl kb0, vector<CoolPropDbl> u)
        : PCSAFT(PCSAFT), kb0(kb0), u(u){}
        CoolPropDbl call(CoolPropDbl R){
            int ncomp = PCSAFT.components.size();
            double error = 0;

            vector<double> pp(ncomp, 0);
            double L = 0;

            for (int i = 0; i < ncomp; i++) {
                if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                    pp[i] = PCSAFT.mole_fractions[i] / (1 - R + kb0 * R * exp(u[i]));
                    L += pp[i];
                } else {
                    L += PCSAFT.mole_fractions[i];
                }
            }
            L = (1 - R) * L;

            error = pow((L + PCSAFT._Q - 1), 2.);
            return error;
        };
    };

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (PCSAFT.ion_term && PCSAFT.components[i].getZ() != 0) {
            x_ions += PCSAFT.mole_fractions[i];
        }
    }

    // initialize variables
    vector<double> k(ncomp, 0), u(ncomp, 0), kprime(ncomp, 0), uprime(ncomp, 0);
    double Pref = p_guess - 0.01 * p_guess;
    double Pprime = p_guess + 0.01 * p_guess;
    if (p_guess > 1e6) { // when close to the critical pressure then we need to have Pprime be less than p_guess
        Pprime = p_guess - 0.005 * p_guess;
    }
    double p = p_guess;

    // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT._T, p, iphase_liquid);
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT._T, p, iphase_gas);
    if ((PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-4) {
        throw SolutionError("liquid and vapor densities are the same.");
    }
    vector<double> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    vector<double> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();

    double xv_sum = 0;
    double xl_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) { // this if statement sets k to 0 for ionic components
            k[i] = fugcoef_l[i] / fugcoef_v[i];
        } else {
            k[i] = 0;
        }
        PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
        xl_sum += PCSAFT.SatL->mole_fractions[i];
        PCSAFT.SatV->mole_fractions[i] = k[i] * PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
        xv_sum += PCSAFT.SatV->mole_fractions[i];
    }

    if (xv_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            PCSAFT.SatV->mole_fractions[i] = PCSAFT.SatV->mole_fractions[i] / xv_sum;
        }
    }

    if (xl_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            PCSAFT.SatL->mole_fractions[i] = PCSAFT.SatL->mole_fractions[i] / xl_sum;
        }
    }

    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT._T, p, iphase_liquid);
    fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT._T, p, iphase_gas);
    fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
    for (int i = 0; i < ncomp; i++) {
        k[i] = fugcoef_l[i] / fugcoef_v[i];
        u[i] = std::log(k[i] / kb);
    }

    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT._T, Pprime, iphase_liquid);
    fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT._T, Pprime, iphase_gas);
    fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
    for (int i = 0; i < ncomp; i++) {
        kprime[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    vector<double> t_weight(ncomp);
    double t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = PCSAFT.SatV->mole_fractions[i] * dlnk_dt / (1 + PCSAFT._Q * (k[i] - 1));
        t_sum += t_weight[i];
    }

    double kb = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
            kb += wi * std::log(k[i]);
        }
    }
    kb = std::exp(kb);

    t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = PCSAFT.SatV->mole_fractions[i] * dlnk_dt / (1 + PCSAFT._Q * (kprime[i] - 1));
        t_sum += t_weight[i];
    }

    double kbprime = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
            kbprime += wi * std::log(kprime[i]);
        }
    }
    kbprime = std::exp(kbprime);
    double kb0 = kbprime;

    for (int i = 0; i < ncomp; i++) {
        u[i] = std::log(k[i] / kb);
        uprime[i] = std::log(kprime[i] / kbprime);
    }

    double B = std::log(kbprime / kb) / (1/Pprime - 1/p);
    double A = std::log(kb) - B * (1/p - 1/Pref);

    if (B < 0) {
        throw SolutionError("B < 0 in outerTQ");
    }

    // solve
    SolverInnerResid resid(*this, kb0, u);

    vector<double> pp(ncomp, 0);
    double maxdif = 1e10 * TOL;
    int itr = 0;
    double Rmin = 0, Rmax = 1;
    while (maxdif > TOL && itr < MAXITER) {
        // save previous values for calculating the difference at the end of the iteration
        vector<double> u_old = u;
        double A_old = A;

        double R0 = kb * PCSAFT._Q / (kb * PCSAFT._Q + kb0 * (1 - PCSAFT._Q));
        resid.u = u;
        double R = R0;
        if (resid.call(R) > TOL) {
            R = BoundedSecant(resid, R0, Rmin, Rmax, DBL_EPSILON, TOL, MAXITER);
        }

        double pp_sum = 0;
        double eupp_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            pp[i] = PCSAFT.mole_fractions[i] / (1 - R + kb0 * R * std::exp(u[i]));
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                pp_sum += pp[i];
                eupp_sum += std::exp(u[i]) * pp[i];
            }
        }
        kb = pp_sum / eupp_sum;

        p = 1 / (1 / Pref + (std::log(kb) - A) / B);
        for (int i = 0; i < ncomp; i++) {
            if (x_ions == 0) {
                PCSAFT.SatL->mole_fractions[i] = pp[i] / pp_sum;
                PCSAFT.SatV->mole_fractions[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                PCSAFT.SatL->mole_fractions[i] = pp[i] / pp_sum * (1 - x_ions/(1 - PCSAFT._Q));
                PCSAFT.SatV->mole_fractions[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else {
                PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 - PCSAFT._Q);
                PCSAFT.SatV->mole_fractions[i] = 0;
            }
        }

        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT._T, p, iphase_liquid);
        vector<CoolPropDbl> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT._T, p, iphase_gas);
        vector<CoolPropDbl> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
        for (int i = 0; i < ncomp; i++) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
            u[i] = std::log(k[i] / kb);
        }

        if (itr == 0) {
            B = std::log(kbprime / kb) / (1/Pprime - 1/p);
        }
        A = std::log(kb) - B * (1/p - 1/Pref);

        maxdif = std::abs(A - A_old);
        for (int i = 0; i < ncomp; i++) {
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                double dif = std::abs(u[i] - u_old[i]);
                if (dif > maxdif) {
                    maxdif = dif;
                } else if (!std::isfinite(dif)) {
                    maxdif = dif;
                }
            }
        }
        itr += 1;
    }

    if (!std::isfinite(p) || !std::isfinite(maxdif) || maxdif > 0.1 || p < 0) {
        throw SolutionError("outerTQ did not converge to a solution");
    }

    return p;
}

double PCSAFTBackend::estimate_flash_t(PCSAFTBackend &PCSAFT) {
    /**
    Get a quick estimate of the temperature at which VLE occurs
    */
    double t_guess = _HUGE;
    int ncomp = N; // number of components

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (PCSAFT.ion_term && PCSAFT.components[i].getZ() != 0) {
            x_ions += PCSAFT.mole_fractions[i];
        }
    }

    bool guess_found = false;
    double t_step = 30;
    double t_start = 571;
    double t_lbound = 1;
    if (PCSAFT.ion_term) {
        t_step = 15;
        t_start = 350;
        t_lbound = 264;
    }
    while (!guess_found && t_start > t_lbound) {
        // initialize variables
        double Tprime = t_start - 50;
        double t = t_start;

        PCSAFT.SatL->_T = t; // _T must be updated because the density calculation depends on it
        PCSAFT.SatV->_T = t;

        // calculate sigma for water, if it is present
        if (PCSAFT.water_present) {
            PCSAFT.components[water_idx].calc_water_sigma(t);
            PCSAFT.SatL->components[water_idx].calc_water_sigma(t);
            PCSAFT.SatV->components[water_idx].calc_water_sigma(t);
            PCSAFT.dielc = dielc_water(t); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
            PCSAFT.SatL->dielc = dielc_water(t);
            PCSAFT.SatV->dielc = dielc_water(t);
        }

        try {
            double p1 = estimate_flash_p(PCSAFT);
            PCSAFT.SatL->_T = Tprime;
            PCSAFT.SatV->_T = Tprime;
            double p2 = estimate_flash_p(PCSAFT);
            PCSAFT.SatL->_T = t; // reset to initial value
            PCSAFT.SatV->_T = t;

            double slope = (std::log10(p1) - std::log10(p2)) / (1/t - 1/Tprime);
            double intercept = std::log10(p1) - slope * (1/t);
            t_guess = slope / (std::log10(PCSAFT._p) - intercept);
            guess_found = true;
        } catch (const SolutionError& ex) {
            t_start -= t_step;
        }
    }

    if (!guess_found) {
        throw SolutionError("an estimate for the VLE temperature could not be found");
    }

    return t_guess;
}


double PCSAFTBackend::estimate_flash_p(PCSAFTBackend &PCSAFT) {
    /**
    Get a quick estimate of the pressure at which VLE occurs
    */
    double p_guess = _HUGE;
    int ncomp = N; // number of components

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (PCSAFT.ion_term && PCSAFT.components[i].getZ() != 0) {
            x_ions += PCSAFT.mole_fractions[i];
        }
    }

    bool guess_found = false;
    double p_start = 10000;
    while (!guess_found && p_start < 1e7) {
        // initialize variables
        vector<double> k(ncomp, 0), u(ncomp, 0), kprime(ncomp, 0), uprime(ncomp, 0);
        double Pprime = 0.99 * p_start;
        double p = p_start;

        // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT._T, p, iphase_liquid);
        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT._T, p, iphase_gas);
        if ((PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        vector<double> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
        vector<double> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();


        double xv_sum = 0;
        double xl_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                k[i] = fugcoef_l[i] / fugcoef_v[i];
            } else {
                k[i] = 0; // set k to 0 for ionic components
            }
            PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
            xl_sum += PCSAFT.SatL->mole_fractions[i];
            PCSAFT.SatV->mole_fractions[i] = k[i] * PCSAFT.mole_fractions[i] / (1 + PCSAFT._Q * (k[i] - 1));
            xv_sum += PCSAFT.SatV->mole_fractions[i];
        }

        if (xv_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                PCSAFT.SatV->mole_fractions[i] = PCSAFT.SatV->mole_fractions[i] / xv_sum;
            }
        }

        if (xl_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                PCSAFT.SatL->mole_fractions[i] = PCSAFT.SatL->mole_fractions[i] / xl_sum;
            }
        }

        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT.SatL->_T, p, iphase_liquid);
        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT.SatV->_T, p, iphase_gas);
        if ((PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
        fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
        double numer = 0;
        double denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                numer += PCSAFT.SatL->mole_fractions[i] * fugcoef_l[i];
                denom += PCSAFT.SatV->mole_fractions[i] * fugcoef_v[i];
            }
        }
        double ratio = numer / denom;

        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(PCSAFT.SatL->_T, Pprime, iphase_liquid);
        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(PCSAFT.SatV->_T, Pprime, iphase_gas);
        if ((PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
        fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
        numer = 0;
        denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (!PCSAFT.ion_term || PCSAFT.components[i].getZ() == 0) {
                numer += PCSAFT.SatL->mole_fractions[i] * fugcoef_l[i];
                denom += PCSAFT.SatV->mole_fractions[i] * fugcoef_v[i];
            }
        }
        double ratio_prime = numer / denom;

        double slope = (std::log10(ratio) - std::log10(ratio_prime)) / (std::log10(p) - std::log10(Pprime));
        double intercept = std::log10(ratio) - slope * std::log10(p);
        p_guess = std::pow(10, -intercept / slope);

        guess_found = true;
    }

    if (!guess_found) {
        throw SolutionError("an estimate for the VLE pressure could not be found");
    }

    return p_guess;
}

CoolPropDbl PCSAFTBackend::solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, phases phase) {
    // Define the residual to be driven to zero
    class SolverRhoResid : public FuncWrapper1D
    {
       public:
        PCSAFTBackend& PCSAFT;
        CoolPropDbl T, p;

        SolverRhoResid(PCSAFTBackend& PCSAFT, CoolPropDbl T, CoolPropDbl p) : PCSAFT(PCSAFT), T(T), p(p) {}
        CoolPropDbl call(CoolPropDbl rhomolar) {
            CoolPropDbl peos = PCSAFT.update_DmolarT(rhomolar);
            double cost = (peos - p) / p;
            if (ValidNumber(cost)) {
                return cost;
            } else {
                return 1.0e20;
            }
        };
    };

    SolverRhoResid resid(*this, T, p);

    // split into grid and find bounds for each root
    vector<double> x_lo, x_hi;
    int num_pts = 20;
    double limit_lower = -8; // first use a log scale for the low density region
    double limit_upper = -1;
    double rho_guess = 1e-13;
    double rho_guess_prev = rho_guess;
    double err_prev = (update_DmolarT(reduced_to_molar(rho_guess, T)) - p) / p;
    for (int i = 0; i < num_pts; i++) {
        rho_guess = pow(10, (limit_upper - limit_lower) / (double)num_pts * i + limit_lower);
        double err = (update_DmolarT(reduced_to_molar(rho_guess, T)) - p) / p;
        if (err * err_prev < 0) {
            x_lo.push_back(rho_guess_prev);
            x_hi.push_back(rho_guess);
        }
        err_prev = err;
        rho_guess_prev = rho_guess;
    }

    limit_lower = 0.1; // for the high density region the log scale is not needed
    limit_upper = 0.7405;
    for (int i = 0; i < num_pts; i++) {
        rho_guess = (limit_upper - limit_lower) / (double)num_pts * i + limit_lower;
        double err = (update_DmolarT(reduced_to_molar(rho_guess, T)) - p) / p;
        if (err * err_prev < 0) {
            x_lo.push_back(rho_guess_prev);
            x_hi.push_back(rho_guess);
        }
        err_prev = err;
        rho_guess_prev = rho_guess;
    }

    // solve for appropriate root(s)
    double rho = _HUGE;
    double x_lo_molar = 1e-8, x_hi_molar = 1e7;

    if (x_lo.size() == 1) {
        rho_guess = reduced_to_molar((x_lo[0] + x_hi[0]) / 2., T);
        x_lo_molar = reduced_to_molar(x_lo[0], T);
        x_hi_molar = reduced_to_molar(x_hi[0], T);
        rho = Brent(resid, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
    } else if (x_lo.size() <= 3 && !x_lo.empty()) {
        if ((phase == iphase_liquid) || (phase == iphase_supercritical_liquid)) {
            rho_guess = reduced_to_molar((x_lo.back() + x_hi.back()) / 2., T);
            x_lo_molar = reduced_to_molar(x_lo.back(), T);
            x_hi_molar = reduced_to_molar(x_hi.back(), T);
            rho = Brent(resid, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
        } else if ((phase == iphase_gas) || (phase == iphase_supercritical_gas) || (phase == iphase_supercritical)) {
            rho_guess = reduced_to_molar((x_lo[0] + x_hi[0]) / 40., T);  // starting with a lower guess often provides better results
            x_lo_molar = reduced_to_molar(x_lo[0], T);
            x_hi_molar = reduced_to_molar(x_hi[0], T);
            rho = Brent(resid, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
        }
    } else if (x_lo.size() > 3) {
        // if multiple roots to check, then find the one with the minimum gibbs energy. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
        double g_min = 1e60;
        for (int i = 0; i < x_lo.size(); i++) {
            rho_guess = reduced_to_molar((x_lo[i] + x_hi[i]) / 2., T);
            x_lo_molar = reduced_to_molar(x_lo[i], T);
            x_hi_molar = reduced_to_molar(x_hi[i], T);
            double rho_i = Brent(resid, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
            double rho_original = this->_rhomolar;
            this->_rhomolar = rho_i;
            double g_i = calc_gibbsmolar_residual();
            this->_rhomolar = rho_original;
            if (g_i < g_min) {
                g_min = g_i;
                rho = rho_i;
            }
        }
    } else {
        int num_pts = 25;
        double err_min = 1e40;
        double rho_min;
        for (int i = 0; i < num_pts; i++) {
            double rho_guess = (0.7405 - 1e-8) / (double)num_pts * i + 1e-8;
            double err = (update_DmolarT(reduced_to_molar(rho_guess, T)) - p) / p;
            if (abs(err) < err_min) {
                err_min = abs(err);
                rho_min = reduced_to_molar(rho_guess, T);
            }
        }
        rho = rho_min;
    }

    return rho;
}

CoolPropDbl PCSAFTBackend::reduced_to_molar(CoolPropDbl nu, CoolPropDbl T) {
    vector<CoolPropDbl> d(N);
    CoolPropDbl summ = 0.;
    for (int i = 0; i < N; i++) {
        d[i] = components[i].getSigma() * (1 - 0.12 * exp(-3 * components[i].getU() / T));
        summ += mole_fractions[i] * components[i].getM() * pow(d[i], 3.);
    }
    return 6 / PI * nu / summ * 1.0e30 / N_AV;
}

CoolPropDbl PCSAFTBackend::calc_molar_mass(void) {
    double summer = 0;
    for (unsigned int i = 0; i < N; ++i) {
        summer += mole_fractions[i] * components[i].molar_mass();
    }
    return summer;
}


vector<double> PCSAFTBackend::XA_find(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x) {
    /**Iterate over this function in order to solve for XA*/
    int num_sites = XA_guess.size();
    vector<double> XA = XA_guess;

    int idxij = -1; // index for delta_ij
    for (int i = 0; i < num_sites; i++) {
        double summ = 0.;
        for (int j = 0; j < num_sites; j++) {
            idxij += 1;
            summ += den*x[j]*XA_guess[j]*delta_ij[idxij];
        }
        XA[i] = 1./(1.+summ);
    }

    return XA;
}


vector<double> PCSAFTBackend::dXAdt_find(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x) {
    /**Solve for the derivative of XA with respect to temperature.*/
    int num_sites = XA.size();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    double summ;
    int ij = 0;
    for (int i = 0; i < num_sites; i++) {
        summ = 0;
        for (int j = 0; j < num_sites; j++) {
            B(i) -= x[j]*XA[j]*ddelta_dt[ij];
            A(i,j) = x[j]*delta_ij[ij];
            summ += x[j]*XA[j]*delta_ij[ij];
            ij += 1;
        }
        A(i,i) = pow(1+den*summ, 2.)/den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dt(num_sites);
    for (int i = 0; i < num_sites; i++) {
        dXA_dt[i] = solution(i);
    }
    return dXA_dt;
}


vector<double> PCSAFTBackend::dXAdx_find(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x) {
    /**Solve for the derivative of XA with respect to composition, or actually with respect
    to rho_i (the molar density of component i, which equals x_i * rho).*/
    int num_sites = XA.size();
    int ncomp = assoc_num.size();
    Eigen::MatrixXd B(num_sites*ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites*ncomp, num_sites*ncomp);

    double sum1, sum2;
    int idx1 = 0;
    int ij = 0;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < num_sites; j++) {
            sum1 = 0;
            for (int k = 0; k < num_sites; k++) {
                sum1 = sum1 + den*x[k]*(XA[k]*ddelta_dx[i*num_sites*num_sites + j*num_sites + k]);
                A(ij,i*num_sites+k) = XA[j]*XA[j]*den*x[k]*delta_ij[j*num_sites+k];
            }

            sum2 = 0;
            for (int l = 0; l < assoc_num[i]; l++) {
                sum2 = sum2 + XA[idx1+l]*delta_ij[idx1*num_sites+l*num_sites+j];
            }

            A(ij,ij) = A(ij,ij) + 1;
            B(ij) = -1*XA[j]*XA[j]*(sum1 + sum2);
            ij += 1;
        }
        idx1 += assoc_num[i];
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dx(num_sites*ncomp);
    for (int i = 0; i < num_sites*ncomp; i++) {
        dXA_dx[i] = solution(i);
    }
    return dXA_dx;
}


void PCSAFTBackend::set_assoc_matrix(){
    vector<int> charge; // whether the association site has a partial positive charge (i.e. hydrogen), negative charge, or elements of both (e.g. for acids modelled as type 1)

    for (int i = 0; i < N; i++){
        vector<std::string> assoc_scheme = components[i].getAssocScheme();
        int num_sites = 0;
        int num = assoc_scheme.size();
        for (int j = 0; j < num; j++) {
            switch(get_scheme_index(assoc_scheme[j])) {
                case i1: {
                    charge.push_back(0);
                    num_sites += 1;
                    break;
                }
                case i2a: {
                    vector<int> tmp{0, 0};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 2;
                    break;
                }
                case i2b: {
                    vector<int> tmp{-1, 1};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 2;
                    break;
                }
                case i3a: {
                    vector<int> tmp{0, 0, 0};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 3;
                    break;
                }
                case i3b: {
                    vector<int> tmp{-1, -1, 1};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 3;
                    break;
                }
                case i4a: {
                    vector<int> tmp{0, 0, 0, 0};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 4;
                    break;
                }
                case i4b: {
                    vector<int> tmp{1, 1, 1, -1};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 4;
                    break;
                }
                case i4c: {
                    vector<int> tmp{-1, -1, 1, 1};
                    charge.insert(charge.end(), tmp.begin(), tmp.end());
                    num_sites += 4;
                    break;
                }
                default:
                    throw ValueError(format("%s is not a valid association type.", assoc_scheme[j]));
            }
        }

        assoc_num.push_back(num_sites);
    }

    for (std::vector<int>::iterator i1 = charge.begin(); i1 != charge.end(); i1++) {
        for (std::vector<int>::iterator i2 = charge.begin(); i2 != charge.end(); i2++) {
            if (*i1 == 0 || *i2 == 0) {
                assoc_matrix.push_back(1);
            }
            else if (*i1 == 1 && *i2 == -1) {
                assoc_matrix.push_back(1);
            }
            else if (*i1 == -1 && *i2 == 1) {
                assoc_matrix.push_back(1);
            }
            else {
                assoc_matrix.push_back(0);
            }
        }
    }
}


double PCSAFTBackend::dielc_water(double t) {
    /**
    Return the dielectric constant of water at the given temperature.

    t : double
        Temperature (K)

    This equation was fit to values given in the reference. For temperatures
    from 263.15 to 368.15 K values at 1 bar were used. For temperatures from
    368.15 to 443.15 K values at 10 bar were used.

    Reference:
    D. G. Archer and P. Wang, “The Dielectric Constant of Water and Debye‐Hückel
    Limiting Law Slopes,” J. Phys. Chem. Ref. Data, vol. 19, no. 2, pp. 371–411,
    Mar. 1990.
    */
    double dielc;
    if (t < 263.15) {
        throw ValueError("The current function for the dielectric constant for water is only valid for temperatures above 263.15 K.");
    } else if (t <= 368.15) {
        dielc = 7.6555618295E-04 * _T * _T - 8.1783881423E-01 * _T + 2.5419616803E+02;
    } else if (t <= 443.15) {
        dielc = 0.0005003272124 * _T * _T - 0.6285556029 * _T + 220.4467027;
    } else {
        throw ValueError("The current function for the dielectric constant for water is only valid for temperatures less than 443.15 K.");
    }
    return dielc;
}

} /* namespace CoolProp */
