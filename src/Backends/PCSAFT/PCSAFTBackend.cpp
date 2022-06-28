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
    for (unsigned int i = 0; i < N; ++i) {
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
    for (unsigned int i = 0; i < N; ++i) {
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
                    adip[l] = a0dip[l] + (m_ij - 1) / m_ij * a1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij - 1) / m_ij * b1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * b2dip[l];
                    J2 += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T)
                          * pow(eta, l);  // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
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

        ares_polar = A2 / (1 - A3 / A2);
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    double ares_assoc = 0.;
    if (assoc_term) {
        int a_sites = 2;
        int ncA = 0;     // number of associating compounds
        vector<int> iA;  // indices of associating compounds
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getVolA() != 0) {
                iA.push_back(i);
                ncA += 1;
            }
        }

        vector<double> XA(ncA * a_sites, 0);
        vector<double> eABij(ncA * ncA, 0);
        vector<double> volABij(ncA * ncA, 0);
        vector<double> delta_ij(ncA * ncA, 0);

        // these indices are necessary because we are only using 1D vectors
        int idxa = -1;  // index over only associating compounds
        int idxi = 0;   // index for the ii-th compound
        int idxj = 0;   // index for the jj-th compound
        for (int i = 0; i < ncA; i++) {
            idxi = iA[i] * ncomp + iA[i];
            for (int j = 0; j < ncA; j++) {
                idxa += 1;
                idxj = iA[j] * ncomp + iA[j];
                eABij[idxa] = (components[iA[i]].getUAB() + components[iA[j]].getUAB()) / 2.;
                volABij[idxa] = sqrt(components[iA[i]].getVolA() * components[iA[j]].getVolA())
                                * pow(sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])), 3);
                delta_ij[idxa] = ghs[iA[i] * ncomp + iA[j]] * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
            }
            XA[i * 2] = (-1 + sqrt(1 + 8 * den * delta_ij[i * ncA + i])) / (4 * den * delta_ij[i * ncA + i]);
            if (!ValidNumber(XA[i * 2])) {
                XA[i * 2] = 0.02;
            }
            XA[i * 2 + 1] = XA[i * 2];
        }

        vector<double> x_assoc(ncA);  // mole fractions of only the associating compounds
        for (int i = 0; i < ncA; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 500) && (dif > 1e-9)) {
            ctr += 1;
            XA = XA_find(XA, ncA, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < ncA * 2; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            XA_old = XA;
        }

        ares_assoc = 0.;
        for (int i = 0; i < ncA; i++) {
            for (int k = 0; k < a_sites; k++) {
                ares_assoc += mole_fractions[iA[i]] * (log(XA[i * a_sites + k]) - 0.5 * XA[i * a_sites + k] + 0.5);
            }
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
            m2e2s3 =
              m2e2s3
              + mole_fractions[i] * mole_fractions[j] * components[i].getM() * components[j].getM() * pow(e_ij[idx] / _T, 2) * pow(s_ij[idx], 3);
            ghs[idx] = 1 / (1 - zeta[3]) + (d[i] * d[j] / (d[i] + d[j])) * 3 * zeta[2] / (1 - zeta[3]) / (1 - zeta[3])
                       + pow(d[i] * d[j] / (d[i] + d[j]), 2) * 2 * zeta[2] * zeta[2] / pow(1 - zeta[3], 3);
            ddij_dt = (d[i] * d[j] / (d[i] + d[j])) * (dd_dt[i] / d[i] + dd_dt[j] / d[j] - (dd_dt[i] + dd_dt[j]) / (d[i] + d[j]));
            dghs_dt[idx] = dzeta_dt[3] / pow(1 - zeta[3], 2.)
                           + 3 * (ddij_dt * zeta[2] + (d[i] * d[j] / (d[i] + d[j])) * dzeta_dt[2]) / pow(1 - zeta[3], 2.)
                           + 4 * (d[i] * d[j] / (d[i] + d[j])) * zeta[2]
                               * (1.5 * dzeta_dt[3] + ddij_dt * zeta[2] + (d[i] * d[j] / (d[i] + d[j])) * dzeta_dt[2]) / pow(1 - zeta[3], 3.)
                           + 6 * pow((d[i] * d[j] / (d[i] + d[j])) * zeta[2], 2.) * dzeta_dt[3] / pow(1 - zeta[3], 4.);
        }
    }

    double dadt_hs =
      1 / zeta[0]
      * (3 * (dzeta_dt[1] * zeta[2] + zeta[1] * dzeta_dt[2]) / (1 - zeta[3]) + 3 * zeta[1] * zeta[2] * dzeta_dt[3] / pow(1 - zeta[3], 2.)
         + 3 * pow(zeta[2], 2.) * dzeta_dt[2] / zeta[3] / pow(1 - zeta[3], 2.)
         + pow(zeta[2], 3.) * dzeta_dt[3] * (3 * zeta[3] - 1) / pow(zeta[3], 2.) / pow(1 - zeta[3], 3.)
         + (3 * pow(zeta[2], 2.) * dzeta_dt[2] * zeta[3] - 2 * pow(zeta[2], 3.) * dzeta_dt[3]) / pow(zeta[3], 3.) * log(1 - zeta[3])
         + (zeta[0] - pow(zeta[2], 3) / pow(zeta[3], 2.)) * dzeta_dt[3] / (1 - zeta[3]));

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
                    J2 += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T)
                          * pow(eta, l);  // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_dt += adip[l] * l * pow(eta, l - 1) * dzeta_dt[3]
                              + bdip[l] * e_ij[j * ncomp + j] * (1 / _T * l * pow(eta, l - 1) * dzeta_dt[3] - 1 / pow(_T, 2.) * pow(eta, l));
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

        dadt_polar = (dA2_dt - 2 * A3 / A2 * dA2_dt + dA3_dt) / pow(1 - A3 / A2, 2.);
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    double dadt_assoc = 0.;
    if (assoc_term) {
        int a_sites = 2;
        int ncA = 0;     // number of associating compounds
        vector<int> iA;  // indices of associating compounds
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getVolA() != 0) {
                iA.push_back(i);
                ncA += 1;
            }
        }

        vector<double> XA(ncA * a_sites, 0);
        vector<double> eABij(ncA * ncA, 0);
        vector<double> volABij(ncA * ncA, 0);
        vector<double> delta_ij(ncA * ncA, 0);
        vector<double> ddelta_dt(ncA * ncA, 0);

        // these indices are necessary because we are only using 1D vectors
        int idxa = -1;  // index over only associating compounds
        int idxi = 0;   // index for the ii-th compound
        int idxj = 0;   // index for the jj-th compound
        for (int i = 0; i < ncA; i++) {
            idxi = iA[i] * ncomp + iA[i];
            for (int j = 0; j < ncA; j++) {
                idxa += 1;
                idxj = iA[j] * ncomp + iA[j];
                eABij[idxa] = (components[iA[i]].getUAB() + components[iA[j]].getUAB()) / 2.;
                volABij[idxa] = sqrt(components[iA[i]].getVolA() * components[iA[j]].getVolA())
                                * pow(sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])), 3);
                delta_ij[idxa] = ghs[iA[i] * ncomp + iA[j]] * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
                ddelta_dt[idxa] = pow(s_ij[idxj], 3) * volABij[idxa]
                                  * (-eABij[idxa] / pow(_T, 2) * exp(eABij[idxa] / _T) * ghs[iA[i] * ncomp + iA[j]]
                                     + dghs_dt[iA[i] * ncomp + iA[j]] * (exp(eABij[idxa] / _T) - 1));
            }
            XA[i * 2] = (-1 + sqrt(1 + 8 * den * delta_ij[i * ncA + i])) / (4 * den * delta_ij[i * ncA + i]);
            if (!ValidNumber(XA[i * 2])) {
                XA[i * 2] = 0.02;
            }
            XA[i * 2 + 1] = XA[i * 2];
        }

        vector<double> x_assoc(ncA);  // mole fractions of only the associating compounds
        for (int i = 0; i < ncA; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 500) && (dif > 1e-9)) {
            ctr += 1;
            XA = XA_find(XA, ncA, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < ncA * 2; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            XA_old = XA;
        }

        vector<double> dXA_dt(ncA * a_sites, 0);
        dXA_dt = dXAdt_find(ncA, delta_ij, den, XA, ddelta_dt, x_assoc, a_sites);

        int idx = -1;
        for (int i = 0; i < ncA; i++) {
            for (int j = 0; j < a_sites; j++) {
                idx += 1;
                dadt_assoc += mole_fractions[iA[i]] * (1 / XA[idx] - 0.5) * dXA_dt[idx];
            }
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
            denghs[idx] =
              zeta[3] / (1 - zeta[3]) / (1 - zeta[3])
              + (d[i] * d[j] / (d[i] + d[j])) * (3 * zeta[2] / (1 - zeta[3]) / (1 - zeta[3]) + 6 * zeta[2] * zeta[3] / pow(1 - zeta[3], 3))
              + pow(d[i] * d[j] / (d[i] + d[j]), 2)
                  * (4 * zeta[2] * zeta[2] / pow(1 - zeta[3], 3) + 6 * zeta[2] * zeta[2] * zeta[3] / pow(1 - zeta[3], 4));
        }
    }

    double ares_hs = 1 / zeta[0]
                     * (3 * zeta[1] * zeta[2] / (1 - zeta[3]) + pow(zeta[2], 3.) / (zeta[3] * pow(1 - zeta[3], 2))
                        + (pow(zeta[2], 3.) / pow(zeta[3], 2.) - zeta[0]) * log(1 - zeta[3]));
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

        vector<double> dipmSQ(ncomp, 0);
        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(components[i].getDipm(), 2.) / (components[i].getM() * components[i].getU() * pow(components[i].getSigma(), 3.)) * conv;
        }

        vector<double> adip(5, 0);
        vector<double> bdip(5, 0);
        vector<double> cdip(5, 0);
        double J2, dJ2_det, J3, dJ3_det;
        double m_ij;
        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(components[i].getM() * components[j].getM());
                if (m_ij > 2) {
                    m_ij = 2;
                }
                J2 = 0.;
                dJ2_det = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij - 1) / m_ij * a1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij - 1) / m_ij * b1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * b2dip[l];
                    J2 += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T)
                          * pow(eta, l);  // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_det += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T) * l * pow(eta, l - 1);
                }
                A2 += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                      * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i]
                      * dipmSQ[j] * J2;
                dA2_det += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                           * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum()
                           * dipmSQ[i] * dipmSQ[j] * dJ2_det;
                if (i == j) {
                    dA2_dx[i] += e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3)
                                 / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i] * dipmSQ[j]
                                 * (mole_fractions[i] * mole_fractions[j] * dJ2_det * PI / 6. * den * components[i].getM() * pow(d[i], 3)
                                    + 2 * mole_fractions[j] * J2);
                } else {
                    dA2_dx[i] += e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3)
                                 / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i] * dipmSQ[j]
                                 * (mole_fractions[i] * mole_fractions[j] * dJ2_det * PI / 6. * den * components[i].getM() * pow(d[i], 3)
                                    + mole_fractions[j] * J2);
                }

                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((components[i].getM() * components[j].getM() * components[k].getM()), 1 / 3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    J3 = 0.;
                    dJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk - 1) / m_ijk * c1dip[l] + (m_ijk - 1) / m_ijk * (m_ijk - 2) / m_ijk * c2dip[l];
                        J3 += cdip[l] * pow(eta, l);
                        dJ3_det += cdip[l] * l * pow(eta, (l - 1));
                    }
                    A3 += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                          * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                          / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum()
                          * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * J3;
                    dA3_det += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                               * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                               / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum()
                               * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * dJ3_det;
                    if ((i == j) && (i == k)) {
                        dA3_dx[i] +=
                          e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3)
                          * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3) / s_ij[i * ncomp + j] / s_ij[i * ncomp + k]
                          / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i]
                          * dipmSQ[j] * dipmSQ[k]
                          * (mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * dJ3_det * PI / 6. * den * components[i].getM() * pow(d[i], 3)
                             + 3 * mole_fractions[j] * mole_fractions[k] * J3);
                    } else if ((i == j) || (i == k)) {
                        dA3_dx[i] +=
                          e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3)
                          * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3) / s_ij[i * ncomp + j] / s_ij[i * ncomp + k]
                          / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i]
                          * dipmSQ[j] * dipmSQ[k]
                          * (mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * dJ3_det * PI / 6. * den * components[i].getM() * pow(d[i], 3)
                             + 2 * mole_fractions[j] * mole_fractions[k] * J3);
                    } else {
                        dA3_dx[i] +=
                          e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3)
                          * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3) / s_ij[i * ncomp + j] / s_ij[i * ncomp + k]
                          / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i]
                          * dipmSQ[j] * dipmSQ[k]
                          * (mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * dJ3_det * PI / 6. * den * components[i].getM() * pow(d[i], 3)
                             + mole_fractions[j] * mole_fractions[k] * J3);
                    }
                }
            }
        }

        A2 = -PI * den * A2;
        A3 = -4 / 3. * PI * PI * den * den * A3;
        dA2_det = -PI * den * dA2_det;
        dA3_det = -4 / 3. * PI * PI * den * den * dA3_det;
        for (int i = 0; i < ncomp; i++) {
            dA2_dx[i] = -PI * den * dA2_dx[i];
            dA3_dx[i] = -4 / 3. * PI * PI * den * den * dA3_dx[i];
        }

        vector<double> dapolar_dx(ncomp);
        for (int i = 0; i < ncomp; i++) {
            dapolar_dx[i] = (dA2_dx[i] * (1 - A3 / A2) + (dA3_dx[i] * A2 - A3 * dA2_dx[i]) / A2) / pow(1 - A3 / A2, 2);
        }

        double ares_polar = A2 / (1 - A3 / A2);
        double Zpolar = eta * ((dA2_det * (1 - A3 / A2) + (dA3_det * A2 - A3 * dA2_det) / A2) / (1 - A3 / A2) / (1 - A3 / A2));
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                mu_polar[i] += mole_fractions[j] * dapolar_dx[j];
            }
            mu_polar[i] = ares_polar + Zpolar + dapolar_dx[i] - mu_polar[i];
        }
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    vector<double> mu_assoc(ncomp, 0);
    if (assoc_term) {
        int a_sites = 2;
        int ncA = 0;     // number of associating compounds
        vector<int> iA;  // indices of associating compounds
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getVolA() != 0) {
                iA.push_back(i);
                ncA += 1;
            }
        }

        vector<double> XA(ncA * a_sites, 0);
        vector<double> eABij(ncA * ncA, 0);
        vector<double> volABij(ncA * ncA, 0);
        vector<double> delta_ij(ncA * ncA, 0);
        vector<double> ddelta_dd(ncA * ncA * ncomp, 0);

        // these indices are necessary because we are only using 1D vectors
        int idxa = -1;        // index over only associating compounds
        int idxi = 0;         // index for the ii-th compound
        int idxj = 0;         // index for the jj-th compound
        int idx_ddelta = -1;  // index for ddelta_dd vector
        double dghsd_dd;
        for (int i = 0; i < ncA; i++) {
            idxi = iA[i] * ncomp + iA[i];
            for (int j = 0; j < ncA; j++) {
                idxa += 1;
                idxj = iA[j] * ncomp + iA[j];
                eABij[idxa] = (components[iA[i]].getUAB() + components[iA[j]].getUAB()) / 2.;
                volABij[idxa] = sqrt(components[iA[i]].getVolA() * components[iA[j]].getVolA())
                                * pow(sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])), 3);
                delta_ij[idxa] = ghs[iA[i] * ncomp + iA[j]] * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
                for (int k = 0; k < ncomp; k++) {
                    idx_ddelta += 1;
                    dghsd_dd =
                      PI / 6. * components[k].getM()
                      * (pow(d[k], 3) / (1 - zeta[3]) / (1 - zeta[3])
                         + 3 * d[iA[i]] * d[iA[j]] / (d[iA[i]] + d[iA[j]])
                             * (d[k] * d[k] / (1 - zeta[3]) / (1 - zeta[3]) + 2 * pow(d[k], 3) * zeta[2] / pow(1 - zeta[3], 3))
                         + 2 * pow((d[iA[i]] * d[iA[j]] / (d[iA[i]] + d[iA[j]])), 2)
                             * (2 * d[k] * d[k] * zeta[2] / pow(1 - zeta[3], 3) + 3 * (pow(d[k], 3) * zeta[2] * zeta[2] / pow(1 - zeta[3], 4))));
                    ddelta_dd[idx_ddelta] = dghsd_dd * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
                }
            }
            XA[i * 2] = (-1 + sqrt(1 + 8 * den * delta_ij[i * ncA + i])) / (4 * den * delta_ij[i * ncA + i]);
            if (!ValidNumber(XA[i * 2])) {
                XA[i * 2] = 0.02;
            }
            XA[i * 2 + 1] = XA[i * 2];
        }

        vector<double> x_assoc(ncA);  // mole fractions of only the associating compounds
        for (int i = 0; i < ncA; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 500) && (dif > 1e-9)) {
            ctr += 1;
            XA = XA_find(XA, ncA, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < ncA * 2; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            XA_old = XA;
        }

        vector<double> dXA_dd(ncA * a_sites * ncomp, 0);
        dXA_dd = dXA_find(ncA, ncomp, iA, delta_ij, den, XA, ddelta_dd, x_assoc, a_sites);

        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncA; j++) {
                for (int k = 0; k < a_sites; k++) {
                    mu_assoc[i] += mole_fractions[iA[j]] * den * dXA_dd[i * (ncA * a_sites) + j * a_sites + k] * (1 / XA[j * a_sites + k] - 0.5);
                }
            }
        }

        for (int i = 0; i < ncA; i++) {
            for (int l = 0; l < a_sites; l++) {
                mu_assoc[iA[i]] += log(XA[i * a_sites + l]) - 0.5 * XA[i * a_sites + l];
            }
            mu_assoc[iA[i]] += 0.5 * a_sites;
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

    vector<double> ghs(ncomp, 0);
    vector<double> denghs(ncomp, 0);
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
        }
        ghs[i] = 1 / (1 - zeta[3]) + (d[i] * d[i] / (d[i] + d[i])) * 3 * zeta[2] / (1 - zeta[3]) / (1 - zeta[3])
                 + pow(d[i] * d[i] / (d[i] + d[i]), 2) * 2 * zeta[2] * zeta[2] / pow(1 - zeta[3], 3);
        denghs[i] = zeta[3] / (1 - zeta[3]) / (1 - zeta[3])
                    + (d[i] * d[i] / (d[i] + d[i])) * (3 * zeta[2] / (1 - zeta[3]) / (1 - zeta[3]) + 6 * zeta[2] * zeta[3] / pow(1 - zeta[3], 3))
                    + pow(d[i] * d[i] / (d[i] + d[i]), 2)
                        * (4 * zeta[2] * zeta[2] / pow(1 - zeta[3], 3) + 6 * zeta[2] * zeta[2] * zeta[3] / pow(1 - zeta[3], 4));
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
        summ += mole_fractions[i] * (components[i].getM() - 1) / ghs[i] * denghs[i];
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
        double J2, dJ2_det, J3, dJ3_det;

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
                dJ2_det = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij - 1) / m_ij * a1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij - 1) / m_ij * b1dip[l] + (m_ij - 1) / m_ij * (m_ij - 2) / m_ij * b2dip[l];
                    J2 += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T)
                          * pow(eta, l);  // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_det += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / _T) * l * pow(eta, l - 1);
                }
                A2 += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                      * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum() * dipmSQ[i]
                      * dipmSQ[j] * J2;
                dA2_det += mole_fractions[i] * mole_fractions[j] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T * pow(s_ij[i * ncomp + i], 3)
                           * pow(s_ij[j * ncomp + j], 3) / pow(s_ij[i * ncomp + j], 3) * components[i].getDipnum() * components[j].getDipnum()
                           * dipmSQ[i] * dipmSQ[j] * dJ2_det;
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
                    dJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk - 1) / m_ijk * c1dip[l] + (m_ijk - 1) / m_ijk * (m_ijk - 2) / m_ijk * c2dip[l];
                        J3 += cdip[l] * pow(eta, l);
                        dJ3_det += cdip[l] * l * pow(eta, (l - 1));
                    }
                    A3 += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                          * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                          / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum() * components[j].getDipnum()
                          * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * J3;
                    dA3_det += mole_fractions[i] * mole_fractions[j] * mole_fractions[k] * e_ij[i * ncomp + i] / _T * e_ij[j * ncomp + j] / _T
                               * e_ij[k * ncomp + k] / _T * pow(s_ij[i * ncomp + i], 3) * pow(s_ij[j * ncomp + j], 3) * pow(s_ij[k * ncomp + k], 3)
                               / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k] * components[i].getDipnum()
                               * components[j].getDipnum() * components[k].getDipnum() * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * dJ3_det;
                }
            }
        }

        A2 = -PI * den * A2;
        A3 = -4 / 3. * PI * PI * den * den * A3;
        dA2_det = -PI * den * dA2_det;
        dA3_det = -4 / 3. * PI * PI * den * den * dA3_det;

        Zpolar = eta * ((dA2_det * (1 - A3 / A2) + (dA3_det * A2 - A3 * dA2_det) / A2) / (1 - A3 / A2) / (1 - A3 / A2));
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    double Zassoc = 0;
    if (assoc_term) {
        int a_sites = 2;
        int ncA = 0;     // number of associating compounds
        vector<int> iA;  // indices of associating compounds
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getVolA() != 0) {
                iA.push_back(i);
                ncA += 1;
            }
        }

        vector<double> XA(ncA * a_sites, 0);
        vector<double> eABij(ncA * ncA, 0);
        vector<double> volABij(ncA * ncA, 0);
        vector<double> delta_ij(ncA * ncA, 0);
        vector<double> ddelta_dd(ncA * ncA * ncomp, 0);

        // these indices are necessary because we are only using 1D vectors
        int idxa = -1;        // index over only associating compounds
        int idxi = 0;         // index for the ii-th compound
        int idxj = 0;         // index for the jj-th compound
        int idx_ddelta = -1;  // index for ddelta_dd vector
        double dghsd_dd;
        for (int i = 0; i < ncA; i++) {
            idxi = iA[i] * ncomp + iA[i];
            for (int j = 0; j < ncA; j++) {
                idxa += 1;
                idxj = iA[j] * ncomp + iA[j];
                eABij[idxa] = (components[iA[i]].getUAB() + components[iA[j]].getUAB()) / 2.;
                volABij[idxa] = sqrt(components[iA[i]].getVolA() * components[iA[j]].getVolA())
                                * pow(sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])), 3);
                delta_ij[idxa] = ghs[iA[j]] * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
                for (int k = 0; k < ncomp; k++) {
                    idx_ddelta += 1;
                    dghsd_dd =
                      PI / 6. * components[k].getM()
                      * (pow(d[k], 3) / (1 - zeta[3]) / (1 - zeta[3])
                         + 3 * d[iA[i]] * d[iA[j]] / (d[iA[i]] + d[iA[j]])
                             * (d[k] * d[k] / (1 - zeta[3]) / (1 - zeta[3]) + 2 * pow(d[k], 3) * zeta[2] / pow(1 - zeta[3], 3))
                         + 2 * pow((d[iA[i]] * d[iA[j]] / (d[iA[i]] + d[iA[j]])), 2)
                             * (2 * d[k] * d[k] * zeta[2] / pow(1 - zeta[3], 3) + 3 * (pow(d[k], 3) * zeta[2] * zeta[2] / pow(1 - zeta[3], 4))));
                    ddelta_dd[idx_ddelta] = dghsd_dd * (exp(eABij[idxa] / _T) - 1) * pow(s_ij[iA[i] * ncomp + iA[j]], 3) * volABij[idxa];
                }
            }
            XA[i * 2] = (-1 + sqrt(1 + 8 * den * delta_ij[i * ncA + i])) / (4 * den * delta_ij[i * ncA + i]);
            XA[i * 2 + 1] = XA[i * 2];
        }

        vector<double> x_assoc(ncA);  // mole fractions of only the associating compounds
        for (int i = 0; i < ncA; i++) {
            x_assoc[i] = mole_fractions[iA[i]];
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 500) && (dif > 1e-9)) {
            ctr += 1;
            XA = XA_find(XA, ncA, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < ncA * 2; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            XA_old = XA;
        }

        vector<double> dXA_dd(a_sites * ncA * ncomp, 0);
        dXA_dd = dXA_find(ncA, ncomp, iA, delta_ij, den, XA, ddelta_dd, x_assoc, a_sites);

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncA; j++) {
                for (int k = 0; k < a_sites; k++) {
                    summ += mole_fractions[i] * den * mole_fractions[iA[j]] * (1 / XA[j * a_sites + k] - 0.5)
                            * dXA_dd[i * (ncA * a_sites) + j * (a_sites) + k];
                }
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
    // if (_p < 0){ throw ValueError("p is less than zero");}
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
            if (SatV->components[i].getZ() != 0) {  // we make the assumption that ions do not appear in the vapor phase
                summ -= SatV->mole_fractions[i];
                SatV->mole_fractions[i] = 0;
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
            _rhomolar = value1;
            _T = value2;
            if (water_present) {
                components[water_idx].calc_water_sigma(_T);
                dielc = dielc_water(
                  _T);  // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
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
    phases phase;

    double p_input, rho_input;
    double p_bub, p_dew;
    switch (input_pair) {
        case PT_INPUTS:
            p_input = _p;
            rho_input = _rhomolar;
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
            p_bub = _p;
            _p = p_input;
            _rhomolar = rho_input;
            if (_p > p_bub) {
                phase = iphase_liquid;
            } else if (_p == p_bub) {
                phase = iphase_twophase;
            } else {
                _Q = 1;
                SatL->_Q = _Q;
                SatV->_Q = _Q;
                flash_QT(*this);
                p_dew = _p;
                _p = p_input;
                _rhomolar = rho_input;
                if (_p < p_dew) {
                    phase = iphase_gas;
                } else if ((_p <= p_bub) && (_p >= p_dew)) {
                    phase = iphase_twophase;
                } else {
                    phase = iphase_unknown;
                }
            }
            break;
        case DmolarT_INPUTS:
            double rho_bub, rho_dew;
            p_input = _p;
            rho_input = _rhomolar;
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

void PCSAFTBackend::flash_QT(PCSAFTBackend& PCSAFT) {
    CoolPropDbl T = PCSAFT._T;

    class SolverBubblePResid : public FuncWrapper1D
    {
       public:
        PCSAFTBackend& PCSAFT;
        CoolPropDbl T, p;

        SolverBubblePResid(PCSAFTBackend& PCSAFT, CoolPropDbl T) : PCSAFT(PCSAFT), T(T) {}
        CoolPropDbl call(CoolPropDbl p) {
            double error = 0;
            if (p <= 0) {
                error = 1e20;
            } else {
                if (PCSAFT.is_pure_or_pseudopure) {
                    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(T, p, iphase_liquid);
                    vector<CoolPropDbl> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
                    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(T, p, iphase_gas);
                    vector<CoolPropDbl> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
                    error += 100000 * pow(fugcoef_l[0] - fugcoef_v[0], 2.);
                } else {
                    if (PCSAFT.N > 1) {
                        bool reset_mole_fractions = false;
                        for (int i = 0; i < PCSAFT.N; i++) {
                            if (!ValidNumber(PCSAFT.SatL->mole_fractions[i]) || !ValidNumber(PCSAFT.SatV->mole_fractions[i])) {
                                reset_mole_fractions = true;
                            }
                        }
                        if (reset_mole_fractions) {
                            PCSAFT.SatL->mole_fractions = PCSAFT.mole_fractions;
                            PCSAFT.SatV->mole_fractions = PCSAFT.mole_fractions;
                        }
                    }

                    int itr = 0;
                    double dif = 10000.;
                    vector<CoolPropDbl> fugcoef_l(PCSAFT.N), fugcoef_v(PCSAFT.N);

                    double rhol, rhov, summ;
                    vector<CoolPropDbl> xv_old(PCSAFT.N);
                    double x_ions = 0.;  // overall mole fraction of ions in the system
                    for (int i = 0; i < PCSAFT.N; i++) {
                        if (PCSAFT.components[i].getZ() != 0) {
                            x_ions += PCSAFT.mole_fractions[i];
                        }
                    }
                    while ((dif > 1e-9) && (itr < 100)) {
                        xv_old = PCSAFT.SatV->mole_fractions;
                        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(T, p, iphase_liquid);
                        fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
                        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(T, p, iphase_gas);
                        fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();

                        if (PCSAFT._Q > 0.5) {
                            summ = 0.;
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatL->mole_fractions[i] = fugcoef_v[i] * PCSAFT.SatV->mole_fractions[i] / fugcoef_l[i];
                                    summ += PCSAFT.SatL->mole_fractions[i];
                                }
                            }
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatL->mole_fractions[i] =
                                      PCSAFT.SatL->mole_fractions[i] / summ
                                      * (((1 - PCSAFT._Q) - x_ions) / (1 - PCSAFT._Q));  // ensures that mole fractions add up to 1
                                    PCSAFT.SatV->mole_fractions[i] =
                                      (PCSAFT.mole_fractions[i] - (1 - PCSAFT._Q) * PCSAFT.SatL->mole_fractions[i])
                                      / PCSAFT
                                          ._Q;  // if PCSAFT->_Q is close to zero then this equation behaves poorly, and that is why we use this if statement to switch the equation around
                                } else {
                                    PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 - PCSAFT._Q);
                                    PCSAFT.SatV->mole_fractions[i] = 0.;
                                }
                            }
                        } else {
                            summ = 0.;
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatV->mole_fractions[i] = fugcoef_l[i] * PCSAFT.SatL->mole_fractions[i] / fugcoef_v[i];
                                }
                                summ += PCSAFT.SatV->mole_fractions[i];
                            }
                            for (int i = 0; i < PCSAFT.N; i++) {
                                PCSAFT.SatV->mole_fractions[i] = PCSAFT.SatV->mole_fractions[i] / summ;
                                PCSAFT.SatL->mole_fractions[i] =
                                  (PCSAFT.mole_fractions[i] - (PCSAFT._Q) * PCSAFT.SatV->mole_fractions[i]) / (1 - PCSAFT._Q);
                            }
                        }

                        dif = 0;
                        for (int i = 0; i < PCSAFT.N; i++) {
                            dif += abs(PCSAFT.SatV->mole_fractions[i] - xv_old[i]);
                        }
                        itr += 1;
                    }

                    for (int i = 0; i < PCSAFT.N; i++) {
                        if (PCSAFT.components[i].getZ() == 0) {
                            error += pow(PCSAFT.SatL->mole_fractions[i] * fugcoef_l[i] - PCSAFT.SatV->mole_fractions[i] * fugcoef_v[i], 2.);
                        }
                        error += pow(
                          (PCSAFT.mole_fractions[i] - PCSAFT._Q * PCSAFT.SatV->mole_fractions[i] - (1 - PCSAFT._Q) * PCSAFT.SatL->mole_fractions[i]),
                          2.);
                    }
                }

                if (!ValidNumber(error) || (PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-5) {
                    error = 1e20;
                }
            }
            return error;
        };
    };

    SolverBubblePResid resid(*this, T);

    CoolPropDbl p_guess = _HUGE;
    double x_lo = _HUGE;
    double x_hi = _HUGE;
    double x_lbound = -8;
    double x_ubound = 9;

    try {
        // scanning the range of pressures to find a good initial guess
        int npts = 30;
        double err_min = 1e20;
        int ctr_increasing = 0;  // keeps track of the number of steps where the error is increasing instead of decreasing
        for (int i = 0; i < npts; i++) {
            CoolPropDbl p_i = pow(10, ((x_ubound - x_lbound) / (double)npts * i + x_lbound));
            double err = resid.call(p_i);

            if (err < err_min) {
                err_min = err;
                p_guess = p_i;
                x_lo = pow(10, ((x_ubound - x_lbound) / (double)npts * (i - 1) + x_lbound));
                x_hi = pow(10, ((x_ubound - x_lbound) / (double)npts * (i + 1) + x_lbound));
                ctr_increasing = 0;
            } else if (err_min < 1e20) {
                ctr_increasing += 1;
            }

            if (
              ctr_increasing
              > 2) {  // this is necessary because PC-SAFT often gives a second, erroneous VLE at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
                break;
            }
        }

        if (p_guess == _HUGE) {
            throw SolutionError(format("A suitable initial guess for pressure could not be found for the QT flash."));
        }
    } catch (const SolutionError& ex) {
        // scanning the range of pressures to find a good initial guess
        int npts = 500;
        double err_min = 1e20;
        int ctr_increasing = 0;  // keeps track of the number of steps where the error is increasing instead of decreasing
        for (int i = 0; i < npts; i++) {
            CoolPropDbl p_i = pow(10, ((x_ubound - x_lbound) / (double)npts * i + x_lbound));
            double err = resid.call(p_i);

            if (err < err_min) {
                err_min = err;
                p_guess = p_i;
                x_lo = pow(10, ((x_ubound - x_lbound) / (double)npts * (i - 1) + x_lbound));
                x_hi = pow(10, ((x_ubound - x_lbound) / (double)npts * (i + 1) + x_lbound));
                ctr_increasing = 0;
            } else if (err_min < 1e20) {
                ctr_increasing += 1;
            }

            if (
              ctr_increasing
              > 2) {  // this is necessary because PC-SAFT often gives a second, erroneous VLE at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
                break;
            }
        }

        if (p_guess == _HUGE) {
            throw SolutionError(format("A suitable initial guess for pressure could not be found for the QT flash."));
        }
    }

    CoolPropDbl p;
    try {
        p = BoundedSecant(resid, p_guess, x_lo, x_hi, 0.01 * p_guess, 1e-8, 200);
    } catch (const SolutionError& ex) {
        p = BoundedSecant(resid, p_guess, x_lo, x_hi, 0.01 * p_guess, 0.1, 200);
    }

    // Load the outputs
    PCSAFT._p = p;
    PCSAFT._rhomolar = 1 / (PCSAFT._Q / PCSAFT.SatV->_rhomolar + (1 - PCSAFT._Q) / PCSAFT.SatL->_rhomolar);
    PCSAFT._phase = iphase_twophase;
}

void PCSAFTBackend::flash_PQ(PCSAFTBackend& PCSAFT) {
    CoolPropDbl p = PCSAFT._p;

    class SolverTboilResid : public FuncWrapper1D
    {
       public:
        PCSAFTBackend& PCSAFT;
        CoolPropDbl T, p;

        SolverTboilResid(PCSAFTBackend& PCSAFT, CoolPropDbl p) : PCSAFT(PCSAFT), p(p) {}
        CoolPropDbl call(CoolPropDbl T) {
            double error = 0;
            if (T <= 0) {
                error = 1e20;
            } else {
                PCSAFT.SatL->_T = T;  // _T must be updated because the density calculation depends on it
                PCSAFT.SatV->_T = T;

                if (PCSAFT.water_present) {
                    try {
                        PCSAFT.components[PCSAFT.water_idx].calc_water_sigma(T);
                        PCSAFT.SatL->components[PCSAFT.water_idx].calc_water_sigma(T);
                        PCSAFT.SatV->components[PCSAFT.water_idx].calc_water_sigma(T);
                        PCSAFT.dielc = PCSAFT.dielc_water(
                          T);  // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
                        PCSAFT.SatL->dielc = PCSAFT.dielc_water(T);
                        PCSAFT.SatV->dielc = PCSAFT.dielc_water(T);
                    } catch (const ValueError& ex) {
                        return 1e20;
                    }
                }

                if (PCSAFT.is_pure_or_pseudopure) {
                    PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(T, p, iphase_liquid);
                    vector<CoolPropDbl> fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();

                    PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(T, p, iphase_gas);
                    vector<CoolPropDbl> fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();
                    error += 100000 * pow(fugcoef_l[0] - fugcoef_v[0], 2.);
                } else {
                    if (PCSAFT.N > 1) {
                        bool reset_mole_fractions = false;
                        for (int i = 0; i < PCSAFT.N; i++) {
                            if (!ValidNumber(PCSAFT.SatL->mole_fractions[i]) || !ValidNumber(PCSAFT.SatV->mole_fractions[i])) {
                                reset_mole_fractions = true;
                            }
                        }
                        if (reset_mole_fractions) {
                            PCSAFT.SatL->mole_fractions = PCSAFT.mole_fractions;
                            PCSAFT.SatV->mole_fractions = PCSAFT.mole_fractions;
                        }
                    }

                    int itr = 0;
                    double dif = 10000.;
                    vector<CoolPropDbl> fugcoef_l(PCSAFT.N), fugcoef_v(PCSAFT.N);

                    double rhol, rhov, summ;
                    vector<CoolPropDbl> xv_old(PCSAFT.N);
                    double x_ions = 0.;  // overall mole fraction of ions in the system
                    for (int i = 0; i < PCSAFT.N; i++) {
                        if (PCSAFT.components[i].getZ() != 0) {
                            x_ions += PCSAFT.mole_fractions[i];
                        }
                    }
                    while ((dif > 1e-9) && (itr < 100)) {
                        xv_old = PCSAFT.SatV->mole_fractions;
                        PCSAFT.SatL->_rhomolar = PCSAFT.SatL->solver_rho_Tp(T, p, iphase_liquid);
                        fugcoef_l = PCSAFT.SatL->calc_fugacity_coefficients();
                        PCSAFT.SatV->_rhomolar = PCSAFT.SatV->solver_rho_Tp(T, p, iphase_gas);
                        fugcoef_v = PCSAFT.SatV->calc_fugacity_coefficients();

                        if (PCSAFT._Q > 0.5) {
                            summ = 0.;
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatL->mole_fractions[i] = fugcoef_v[i] * PCSAFT.SatV->mole_fractions[i] / fugcoef_l[i];
                                    summ += PCSAFT.SatL->mole_fractions[i];
                                }
                            }
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatL->mole_fractions[i] =
                                      PCSAFT.SatL->mole_fractions[i] / summ
                                      * (((1 - PCSAFT._Q) - x_ions) / (1 - PCSAFT._Q));  // ensures that mole fractions add up to 1
                                    PCSAFT.SatV->mole_fractions[i] =
                                      (PCSAFT.mole_fractions[i] - (1 - PCSAFT._Q) * PCSAFT.SatL->mole_fractions[i])
                                      / PCSAFT
                                          ._Q;  // if PCSAFT->_Q is close to zero then this equation behaves poorly, and that is why we use this if statement to switch the equation around
                                } else {
                                    PCSAFT.SatL->mole_fractions[i] = PCSAFT.mole_fractions[i] / (1 - PCSAFT._Q);
                                    PCSAFT.SatV->mole_fractions[i] = 0.;
                                }
                            }
                        } else {
                            summ = 0.;
                            for (int i = 0; i < PCSAFT.N; i++) {
                                if (PCSAFT.components[i].getZ() == 0) {
                                    PCSAFT.SatV->mole_fractions[i] = fugcoef_l[i] * PCSAFT.SatL->mole_fractions[i] / fugcoef_v[i];
                                }
                                summ += PCSAFT.SatV->mole_fractions[i];
                            }
                            for (int i = 0; i < PCSAFT.N; i++) {
                                PCSAFT.SatV->mole_fractions[i] = PCSAFT.SatV->mole_fractions[i] / summ;
                                PCSAFT.SatL->mole_fractions[i] =
                                  (PCSAFT.mole_fractions[i] - (PCSAFT._Q) * PCSAFT.SatV->mole_fractions[i]) / (1 - PCSAFT._Q);
                            }
                        }

                        dif = 0;
                        for (int i = 0; i < PCSAFT.N; i++) {
                            dif += abs(PCSAFT.SatV->mole_fractions[i] - xv_old[i]);
                        }
                        itr += 1;
                    }

                    for (int i = 0; i < PCSAFT.N; i++) {
                        if (PCSAFT.components[i].getZ() == 0) {
                            error += pow(PCSAFT.SatL->mole_fractions[i] * fugcoef_l[i] - PCSAFT.SatV->mole_fractions[i] * fugcoef_v[i], 2.);
                        }
                        error += pow(
                          (PCSAFT.mole_fractions[i] - PCSAFT._Q * PCSAFT.SatV->mole_fractions[i] - (1 - PCSAFT._Q) * PCSAFT.SatL->mole_fractions[i]),
                          2.);
                    }
                }

                if (!ValidNumber(error) || (PCSAFT.SatL->_rhomolar - PCSAFT.SatV->_rhomolar) < 1e-5) {
                    error = 1e20;
                }
            }
            return error;
        };
    };

    SolverTboilResid resid(*this, p);

    CoolPropDbl t_guess = _HUGE;
    double x_lo = _HUGE;
    double x_hi = _HUGE;

    double x_lbound = 1;
    double x_ubound = 1000;

    try {
        // scan through the range of temperatures to find a good initial guess
        int npts = 40;
        double err_min = 1e20;
        int ctr_increasing = 0;  // keeps track of the number of steps where the error is increasing instead of decreasing
        for (
          int i = npts; i >= 0;
          i--) {  // here we need to scan in the opposite direction (high T to low T) because a second, erroneous VLE occurs at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
            CoolPropDbl T_i = ((x_ubound - x_lbound) / (double)npts * i + x_lbound);
            double err = resid.call(T_i);

            if (err < err_min) {
                err_min = err;
                t_guess = T_i;
                x_lo = ((x_ubound - x_lbound) / (double)npts * (i - 1) + x_lbound);
                x_hi = ((x_ubound - x_lbound) / (double)npts * (i + 1) + x_lbound);
                ctr_increasing = 0;
            } else if (err_min < 1e20) {
                ctr_increasing += 1;
            }

            if (
              ctr_increasing
              > 2) {  // this is necessary because PC-SAFT often gives a second, erroneous VLE at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
                break;
            }
        }

        if (t_guess == _HUGE) {
            throw SolutionError(format("A suitable initial guess for temperature could not be found for the PQ flash."));
        }
    } catch (const SolutionError& ex) {
        // scan through the range of temperatures to find a good initial guess
        int npts = 1000;
        double err_min = 1e20;
        int ctr_increasing = 0;  // keeps track of the number of steps where the error is increasing instead of decreasing
        for (
          int i = npts; i >= 0;
          i--) {  // here we need to scan in the opposite direction (high T to low T) because a second, erroneous VLE occurs at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
            CoolPropDbl T_i = ((x_ubound - x_lbound) / (double)npts * i + x_lbound);
            double err = resid.call(T_i);

            if (err < err_min) {
                err_min = err;
                t_guess = T_i;
                x_lo = ((x_ubound - x_lbound) / (double)npts * (i - 1) + x_lbound);
                x_hi = ((x_ubound - x_lbound) / (double)npts * (i + 1) + x_lbound);
                ctr_increasing = 0;
            } else if (err_min < 1e20) {
                ctr_increasing += 1;
            }

            if (
              ctr_increasing
              > 2) {  // this is necessary because PC-SAFT often gives a second, erroneous VLE at lower temperatures. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
                break;
            }
        }

        if (t_guess == _HUGE) {
            throw SolutionError(format("A suitable initial guess for temperature could not be found for the PQ flash."));
        }
    }

    CoolPropDbl T;
    try {
        T = BoundedSecant(resid, t_guess, x_lo, x_hi, 0.01 * t_guess, 1e-8, 200);
    } catch (const SolutionError& ex) {
        T = BoundedSecant(resid, t_guess, x_lo, x_hi, 0.01 * t_guess, 0.1, 200);
    }

    // Load the outputs
    PCSAFT._T = T;
    PCSAFT._rhomolar = 1 / (PCSAFT._Q / PCSAFT.SatV->_rhomolar + (1 - PCSAFT._Q) / PCSAFT.SatL->_rhomolar);
    PCSAFT._phase = iphase_twophase;
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
    int num_pts = 25;
    double rho_guess = 1e-13;
    double rho_guess_prev = rho_guess;
    double err_prev = (update_DmolarT(reduced_to_molar(rho_guess, T)) - p) / p;
    for (int i = 0; i < num_pts; i++) {
        rho_guess = 0.7405 / (double)num_pts * i + 6e-3;
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
    double x_lo_molar, x_hi_molar;

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
            double rho_guess = 0.7405 / (double)num_pts * i + 1e-8;
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

vector<double> PCSAFTBackend::XA_find(vector<double> XA_guess, int ncA, vector<double> delta_ij, double den, vector<double> x) {
    /**Iterate over this function in order to solve for XA*/
    int n_sites = XA_guess.size() / ncA;
    double summ2;
    vector<CoolPropDbl> XA = XA_guess;

    for (int i = 0; i < ncA; i++) {
        for (int kout = 0; kout < n_sites; kout++) {
            summ2 = 0.;
            for (int j = 0; j < ncA; j++) {
                for (int kin = 0; kin < n_sites; kin++) {
                    if (kin != kout) {
                        summ2 += den * x[j] * XA_guess[j * n_sites + kin] * delta_ij[i * ncA + j];
                    }
                }
            }
            XA[i * n_sites + kout] = 1. / (1. + summ2);
        }
    }

    return XA;
}

vector<double> PCSAFTBackend::dXA_find(int ncA, int ncomp, vector<int> iA, vector<double> delta_ij, double den, vector<double> XA,
                                       vector<double> ddelta_dd, vector<double> x, int n_sites) {
    /**Solve for the derivative of XA with respect to density.*/
    Eigen::MatrixXd B(n_sites * ncA * ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_sites * ncA * ncomp, n_sites * ncA * ncomp);

    double sum1, sum2;
    int indx1, indx2;
    int indx4 = -1;
    int indx3 = -1;
    for (int i = 0; i < ncomp; i++) {
        indx1 = -1;
        if (find(iA.begin(), iA.end(), i) != iA.end()) {
            indx4 += 1;
        }
        for (int j = 0; j < ncA; j++) {
            for (int h = 0; h < n_sites; h++) {
                indx1 += 1;
                indx3 += 1;
                indx2 = -1;
                sum1 = 0;
                for (int k = 0; k < ncA; k++) {
                    for (int l = 0; l < n_sites; l++) {
                        indx2 += 1;
                        sum1 = sum1
                               + den * x[k]
                                   * (XA[indx2] * ddelta_dd[j * (ncA * ncomp) + k * (ncomp) + i]
                                      * ((indx1 + indx2) % 2));  // (indx1+indx2)%2 ensures that A-A and B-B associations are set to zero
                        A(indx1 + i * n_sites * ncA, indx2 + i * n_sites * ncA) =
                          A(indx1 + i * n_sites * ncA, indx2 + i * n_sites * ncA)
                          + XA[indx1] * XA[indx1] * den * x[k] * delta_ij[j * ncA + k] * ((indx1 + indx2) % 2);
                    }
                }

                sum2 = 0;
                if (find(iA.begin(), iA.end(), i) != iA.end()) {
                    for (int k = 0; k < n_sites; k++) {
                        sum2 = sum2 + XA[n_sites * (indx4) + k] * delta_ij[indx4 * ncA + j] * ((indx1 + k) % 2);
                    }
                }

                A(indx3, indx3) = A(indx3, indx3) + 1;
                B(indx3) = -1 * XA[indx1] * XA[indx1] * (sum1 + sum2);
            }
        }
    }

    Eigen::MatrixXd solution = A.lu().solve(B);  //Solves linear system of equations
    vector<double> dXA_dd(n_sites * ncA * ncomp);
    for (int i = 0; i < n_sites * ncA * ncomp; i++) {
        dXA_dd[i] = solution(i);
    }
    return dXA_dd;
}

vector<double> PCSAFTBackend::dXAdt_find(int ncA, vector<double> delta_ij, double den, vector<double> XA, vector<double> ddelta_dt, vector<double> x,
                                         int n_sites) {
    /**Solve for the derivative of XA with respect to temperature.*/
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n_sites * ncA, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_sites * ncA, n_sites * ncA);

    double summ;
    int i_in, i_out = -1;  // i_out is index of outer iteration loop (follows row of matrices)
    for (int i = 0; i < ncA; i++) {
        for (int ai = 0; ai < n_sites; ai++) {
            i_out += 1;
            i_in = -1;  // index for summation loops
            summ = 0;
            for (int j = 0; j < ncA; j++) {
                for (int bj = 0; bj < n_sites; bj++) {
                    i_in += 1;
                    B(i_out) -= x[j] * XA[i_in] * ddelta_dt[i * ncA + j]
                                * ((i_in + i_out) % 2);  // (i_in+i_out)%2 ensures that A-A and B-B associations are set to zero
                    A(i_out, i_in) = x[j] * delta_ij[i * ncA + j] * ((i_in + i_out) % 2);
                    summ += x[j] * XA[i_in] * delta_ij[i * ncA + j] * ((i_in + i_out) % 2);
                }
            }
            A(i_out, i_out) = A(i_out, i_out) + pow(1 + den * summ, 2.) / den;
        }
    }

    Eigen::MatrixXd solution = A.lu().solve(B);  //Solves linear system of equations
    vector<double> dXA_dt(n_sites * ncA);
    for (int i = 0; i < n_sites * ncA; i++) {
        dXA_dt[i] = solution(i);
    }
    return dXA_dt;
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
