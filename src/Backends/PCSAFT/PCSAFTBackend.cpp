#include <vector>
#include <string>
#include <cmath>
#include "math.h"
#include <Eigen/Dense>
#include <stdlib.h>

#include "PCSAFTBackend.h"
#include "Solvers.h"
#include "Backends/Helmholtz/VLERoutines.h"

using namespace std;

namespace CoolProp {

PCSAFTBackend::PCSAFTBackend(const std::vector<std::string> &component_names) { //, bool generate_SatL_and_SatV
    N = component_names.size();
    std::vector<PCSAFTFluid> components(N);
    ion_term = false;
    polar_term = false;
    assoc_term = false;
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
    }

    // Reset the residual Helmholtz energy class
    // residual_helmholtz.reset(new ResidualHelmholtz());

    // Set the components and associated flags
    is_pure_or_pseudopure = (N == 1);
    // set_components(components, generate_SatL_and_SatV);

    // loading interaction parameters
    std::string kij_string;
    if (!is_pure_or_pseudopure) {
        std::vector<double> k_ij(N*N, 0.0);
        for (unsigned int i = 0; i < N; ++i){
            for (unsigned int j = 0; j < N; ++j){
                if (i != j) {
                    kij_string = PCSAFTLibrary::get_library().get_mixture_binary_pair_pcsaft(components[i].getCAS(), components[j].getCAS(), "kij");
                    k_ij[i*N+j] = atof(kij_string.c_str());
                }
            }
        }
    }

    // Set the phase to default unknown value
    _phase = iphase_unknown;
}

PCSAFTBackend::PCSAFTBackend(const std::vector<PCSAFTFluid> &components) { //bool generate_SatL_and_SatV
    N = components.size();
    // Determining which PC-SAFT terms should be used
    ion_term = false;
    polar_term = false;
    assoc_term = false;
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
    }
    // Reset the residual Helmholtz energy class
    // residual_helmholtz.reset(new ResidualHelmholtz());

    // Set the components and associated flags
    is_pure_or_pseudopure = (N == 1);
    // set_components(components, generate_SatL_and_SatV); //!!! maybe need to implement set_components function

    // loading interaction parameters
    std::string kij_string;
    if (!is_pure_or_pseudopure) {
        std::vector<double> k_ij(N*N, 0.0);
        for (unsigned int i = 0; i < N; ++i){
            for (unsigned int j = 0; j < N; ++j){
                if (i != j) {
                    kij_string = PCSAFTLibrary::get_library().get_mixture_binary_pair_pcsaft(components[i].getCAS(), components[j].getCAS(), "kij");
                    k_ij[i*N+j] = atof(kij_string.c_str());
                }
            }
        }
    }

    // Set the phase to default unknown value
    _phase = iphase_unknown;
}

void PCSAFTBackend::set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions)
{
    if (mole_fractions.size() != N)
    {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]",mole_fractions.size(), N));
    }
    // Copy values without reallocating memory
    this->mole_fractions = mole_fractions; // Most effective copy
    this->resize(N); // No reallocation of this->mole_fractions happens
    // Also store the mole fractions as doubles
    this->mole_fractions_double = std::vector<double>(mole_fractions.begin(), mole_fractions.end());
};

void PCSAFTBackend::set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions)
{
    if (mass_fractions.size() != N)
    {
        throw ValueError(format("size of mass fraction vector [%d] does not equal that of component vector [%d]",mass_fractions.size(), N));
    }
    std::vector<CoolPropDbl> moles;
	CoolPropDbl sum_moles = 0.0;
	CoolPropDbl tmp = 0.0;
    for (unsigned int i = 0; i < components.size(); ++i)
    {
        tmp = mass_fractions[i]/components[i].molar_mass();
        moles.push_back(tmp);
		sum_moles += tmp;
    }
	std::vector<CoolPropDbl> mole_fractions;
	for(std::vector< CoolPropDbl >::iterator it = moles.begin(); it != moles.end(); ++it)
    {
		mole_fractions.push_back(*it/sum_moles);
	}
	this->set_mole_fractions(mole_fractions);
};

void PCSAFTBackend::resize(std::size_t N)
{
    this->mole_fractions.resize(N);
    this->mole_fractions_double.resize(N);
    this->K.resize(N);
    this->lnK.resize(N);
}

void PCSAFTBackend::update_DmolarT() {
    _p = this->calc_pressure_nocache(_T, _rhomolar);
}

CoolPropDbl PCSAFTBackend::calc_pressure_nocache(CoolPropDbl t, CoolPropDbl rho) {
    double den = rho*N_AV/1.0e30;

    double Z = this->calc_compressibility_factor(t, rho);
    double P = Z*kb*t*den*1.0e30; // Pa
    return P;
}

CoolPropDbl PCSAFTBackend::calc_alpha0(void) {

}

CoolPropDbl PCSAFTBackend::calc_alphar(void) {

}

CoolPropDbl PCSAFTBackend::calc_hmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar) {
    // Calculate the reducing parameters
    // CoolPropDbl delta = rhomolar/_reducing.rhomolar;
    // CoolPropDbl tau = _reducing.T/T;
    //
    // // Calculate derivatives if needed, or just use cached values
    // // Calculate derivative if needed
    // CoolPropDbl dar_dDelta = calc_alphar_deriv_nocache(0, 1, mole_fractions, tau, delta);
    // CoolPropDbl dar_dTau = calc_alphar_deriv_nocache(1, 0, mole_fractions, tau, delta);
    // CoolPropDbl da0_dTau = calc_alpha0_deriv_nocache(1, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    // CoolPropDbl R_u = gas_constant();
    //
    // // Get molar enthalpy
    // return R_u*T*(1 + tau*(da0_dTau+dar_dTau) + delta*dar_dDelta);
}

CoolPropDbl PCSAFTBackend::calc_hmolar(void){
	// if (get_debug_level()>=50) std::cout << format("HelmholtzEOSMixtureBackend::calc_hmolar: 2phase: %d T: %g rhomomolar: %g", isTwoPhase(), _T, _rhomolar) << std::endl;
  //   if (isTwoPhase())
  //   {
	// 	if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
  //       if (std::abs(_Q) < DBL_EPSILON){
  //           _hmolar = SatL->hmolar();
  //       }
  //       else if (std::abs(_Q-1) < DBL_EPSILON){
  //           _hmolar = SatV->hmolar();
  //       }
  //       else{
  //           _hmolar = _Q*SatV->hmolar() + (1 - _Q)*SatL->hmolar();
  //       }
  //       return static_cast<CoolPropDbl>(_hmolar);
  //   }
  //   else if (isHomogeneousPhase())
  //   {
  //           // Calculate the reducing parameters
  //       _delta = _rhomolar/_reducing.rhomolar;
  //       _tau = _reducing.T/_T;
  //
  //       // Calculate derivatives if needed, or just use cached values
  //       CoolPropDbl da0_dTau = dalpha0_dTau();
  //       CoolPropDbl dar_dTau = dalphar_dTau();
  //       CoolPropDbl dar_dDelta = dalphar_dDelta();
  //       CoolPropDbl R_u = gas_constant();
  //
  //       // Get molar enthalpy
  //       _hmolar = R_u*_T*(1 + _tau.pt()*(da0_dTau+dar_dTau) + _delta.pt()*dar_dDelta);
  //
  //       return static_cast<CoolPropDbl>(_hmolar);
  //   }
  //   else{
  //       throw ValueError(format("phase is invalid in calc_hmolar"));
  //   }
}

CoolPropDbl PCSAFTBackend::calc_smolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar){
    // Calculate the reducing parameters
    // CoolPropDbl delta = rhomolar/_reducing.rhomolar;
    // CoolPropDbl tau = _reducing.T/T;
    //
    // // Calculate derivatives if needed, or just use cached values
    // // Calculate derivative if needed
    // CoolPropDbl dar_dTau = calc_alphar_deriv_nocache(1, 0, mole_fractions, tau, delta);
    // CoolPropDbl ar = calc_alphar_deriv_nocache(0, 0, mole_fractions, tau, delta);
    // CoolPropDbl da0_dTau = calc_alpha0_deriv_nocache(1, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    // CoolPropDbl a0 = calc_alpha0_deriv_nocache(0, 0, mole_fractions, tau, delta, _reducing.T, _reducing.rhomolar);
    // CoolPropDbl R_u = gas_constant();
    //
    // // Get molar entropy
    // return R_u*(tau*(da0_dTau+dar_dTau) - a0 - ar);
}

CoolPropDbl PCSAFTBackend::calc_smolar(void){
    // if (isTwoPhase())
    // {
		// if (!this->SatL || !this->SatV) throw ValueError(format("The saturation properties are needed for the two-phase properties"));
    //     if (std::abs(_Q) < DBL_EPSILON){
    //         _smolar = SatL->smolar();
    //     }
    //     else if (std::abs(_Q-1) < DBL_EPSILON){
    //         _smolar = SatV->smolar();
    //     }
    //     else{
    //         _smolar = _Q*SatV->smolar() + (1 - _Q)*SatL->smolar();
    //     }
    //     return static_cast<CoolPropDbl>(_smolar);
    // }
    // else if (isHomogeneousPhase())
    // {
    //     // Calculate the reducing parameters
    //     _delta = _rhomolar/_reducing.rhomolar;
    //     _tau = _reducing.T/_T;
    //
    //     // Calculate derivatives if needed, or just use cached values
    //     CoolPropDbl da0_dTau = dalpha0_dTau();
    //     CoolPropDbl ar = alphar();
    //     CoolPropDbl a0 = alpha0();
    //     CoolPropDbl dar_dTau = dalphar_dTau();
    //     CoolPropDbl R_u = gas_constant();
    //
    //     // Get molar entropy
    //     _smolar = R_u*(_tau.pt()*(da0_dTau+dar_dTau) - a0 - ar);
    //
    //     return static_cast<CoolPropDbl>(_smolar);
    // }
    // else{
    //     throw ValueError(format("phase is invalid in calc_smolar"));
    // }
}

CoolPropDbl PCSAFTBackend::calc_fugacity_coefficient(void) {

}

CoolPropDbl PCSAFTBackend::calc_gibbsmolar(void) {

}

CoolPropDbl PCSAFTBackend::calc_cpmolar(void) {

}

CoolPropDbl PCSAFTBackend::calc_compressibility_factor(double t, double rho){
    int ncomp = N; // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = components[i].getSigma()*(1-0.12*exp(-3*components[i].getU()/t));
    }
    if (ion_term) {
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getZ() != 0) {
                d[i] = components[i].getSigma()*(1-0.12); // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }

    double den = rho*N_AV/1.0e30;

    vector<double> zeta (4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += mole_fractions[j]*components[j].getM()*pow(d[j], i);
        }
        zeta[i] = PI/6*den*summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += mole_fractions[i]*components[i].getM();
    }

    vector<double> ghs (ncomp, 0);
    vector<double> denghs (ncomp, 0);
    vector<double> e_ij (ncomp*ncomp, 0);
    vector<double> s_ij (ncomp*ncomp, 0);
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
                        e_ij[idx] = sqrt(components[i].getU()*components[j].getU())*(1-k_ij[idx]);
                    }
                }
            } else {
                if (k_ij.empty()) {
                    e_ij[idx] = sqrt(components[i].getU()*components[j].getU());
                }
                else {
                    e_ij[idx] = sqrt(components[i].getU()*components[j].getU())*(1-k_ij[idx]);
                }
            }
            m2es3 = m2es3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*e_ij[idx]/t*pow(s_ij[idx], 3);
            m2e2s3 = m2e2s3 + mole_fractions[i]*mole_fractions[j]*components[i].getM()*components[j].getM()*pow(e_ij[idx]/t,2)*pow(s_ij[idx], 3);
        }
        ghs[i] = 1/(1-zeta[3]) + (d[i]*d[i]/(d[i]+d[i]))*3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
            pow(d[i]*d[i]/(d[i]+d[i]), 2)*2*zeta[2]*zeta[2]/pow(1-zeta[3], 3);
        denghs[i] = zeta[3]/(1-zeta[3])/(1-zeta[3]) +
            (d[i]*d[i]/(d[i]+d[i]))*(3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
            6*zeta[2]*zeta[3]/pow(1-zeta[3], 3)) +
            pow(d[i]*d[i]/(d[i]+d[i]), 2)*(4*zeta[2]*zeta[2]/pow(1-zeta[3], 3) +
            6*zeta[2]*zeta[2]*zeta[3]/pow(1-zeta[3], 4));
    }

    double Zhs = zeta[3]/(1-zeta[3]) + 3.*zeta[1]*zeta[2]/zeta[0]/(1.-zeta[3])/(1.-zeta[3]) +
        (3.*pow(zeta[2], 3.) - zeta[3]*pow(zeta[2], 3.))/zeta[0]/pow(1.-zeta[3], 3.);

    static double a0[7] = { 0.910563145, 0.636128145, 2.686134789, -26.54736249, 97.75920878, -159.5915409, 91.29777408 };
    static double a1[7] = { -0.308401692, 0.186053116, -2.503004726, 21.41979363, -65.25588533, 83.31868048, -33.74692293 };
    static double a2[7] = { -0.090614835, 0.452784281, 0.596270073, -1.724182913, -4.130211253, 13.77663187, -8.672847037 };
    static double b0[7] = { 0.724094694, 2.238279186, -4.002584949, -21.00357682, 26.85564136, 206.5513384, -355.6023561 };
    static double b1[7] = { -0.575549808, 0.699509552, 3.892567339, -17.21547165, 192.6722645, -161.8264617, -165.2076935 };
    static double b2[7] = { 0.097688312, -0.255757498, -9.155856153, 20.64207597, -38.80443005, 93.62677408, -29.66690559 };

    vector<double> a(7, 0);
    vector<double> b(7, 0);
    for (int i = 0; i < 7; i++) {
        a[i] = a0[i] + (m_avg-1.)/m_avg*a1[i] + (m_avg-1.)/m_avg*(m_avg-2.)/m_avg*a2[i];
        b[i] = b0[i] + (m_avg-1.)/m_avg*b1[i] + (m_avg-1.)/m_avg*(m_avg-2.)/m_avg*b2[i];
    }

    double detI1_det = 0.0;
    double detI2_det = 0.0;
    double I2 = 0.0;
    for (int i = 0; i < 7; i++) {
        detI1_det += a[i]*(i+1)*pow(eta, i);
        detI2_det += b[i]*(i+1)*pow(eta, i);
        I2 += b[i]*pow(eta, i);
    }
    double C1 = 1./(1. + m_avg*(8*eta-2*eta*eta)/pow(1-eta, 4) + (1-m_avg)*(20*eta-27*eta*eta+12*pow(eta, 3)-2*pow(eta, 4))/pow((1-eta)*(2-eta), 2.0));
    double C2 = -1.*C1*C1*(m_avg*(-4*eta*eta+20*eta+8)/pow(1-eta, 5) + (1-m_avg)*(2*pow(eta, 3)+12*eta*eta-48*eta+40)/pow((1-eta)*(2-eta), 3.0));

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += mole_fractions[i]*(components[i].getM()-1)/ghs[i]*denghs[i];
    }

    double Zid = 1.0;
    double Zhc = m_avg*Zhs - summ;
    double Zdisp = -2*PI*den*detI1_det*m2es3 - PI*den*m_avg*(C1*detI2_det + C2*eta*I2)*m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double Zpolar = 0;
    if (polar_term) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_det = 0.;
        double dA3_det = 0.;
        vector<double> adip (5, 0);
        vector<double> bdip (5, 0);
        vector<double> cdip (5, 0);
        vector<double> dipmSQ (ncomp, 0);
        double J2, dJ2_det, J3, dJ3_det;

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

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(components[i].getDipm(), 2.)/(components[i].getM()*components[i].getU()*pow(components[i].getSigma(),3.))*conv;
        }

        double m_ij;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(components[i].getM()*components[j].getM());
                if (m_ij > 2) {
                    m_ij = 2;
                }
                J2 = 0.;
                dJ2_det = 0.;
                for (int l = 0; l < 5; l++) {
                    adip[l] = a0dip[l] + (m_ij-1)/m_ij*a1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*a2dip[l];
                    bdip[l] = b0dip[l] + (m_ij-1)/m_ij*b1dip[l] + (m_ij-1)/m_ij*(m_ij-2)/m_ij*b2dip[l];
                    J2 += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*pow(eta, l); // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_det += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*l*pow(eta, l-1);
                }
                A2 += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*J2;
                dA2_det += mole_fractions[i]*mole_fractions[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*
                    pow(s_ij[j*ncomp+j],3)/pow(s_ij[i*ncomp+j],3)*components[i].getDipnum()*components[j].getDipnum()*dipmSQ[i]*dipmSQ[j]*dJ2_det;
            }
        }

        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((components[i].getM()*components[j].getM()*components[k].getM()),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    J3 = 0.;
                    dJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        cdip[l] = c0dip[l] + (m_ijk-1)/m_ijk*c1dip[l] + (m_ijk-1)/m_ijk*(m_ijk-2)/m_ijk*c2dip[l];
                        J3 += cdip[l]*pow(eta, l);
                        dJ3_det += cdip[l]*l*pow(eta, (l-1));
                    }
                    A3 += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                    dA3_det += mole_fractions[i]*mole_fractions[j]*mole_fractions[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*components[i].getDipnum()*components[j].getDipnum()*components[k].getDipnum()*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*dJ3_det;
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        dA2_det = -PI*den*dA2_det;
        dA3_det = -4/3.*PI*PI*den*den*dA3_det;

        Zpolar = eta*((dA2_det*(1-A3/A2)+(dA3_det*A2-A3*dA2_det)/A2)/(1-A3/A2)/(1-A3/A2));
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    double Zassoc = 0;
    if (assoc_term) {
        int a_sites = 2;
        int ncA = 0; // number of associating compounds
        vector<int> iA; // indices of associating compounds
        for (int i = 0; i < ncomp; i++) {
            if (components[i].getVolA() != 0) {
                iA.push_back(i);
                ncA += 1;
            }
        }

        vector<double> XA (ncA*a_sites, 0);
        vector<double> eABij (ncA*ncA, 0);
        vector<double> volABij (ncA*ncA, 0);
        vector<double> delta_ij (ncA*ncA, 0);
        vector<double> ddelta_dd (ncA*ncA*ncomp, 0);

        // these indices are necessary because we are only using 1D vectors
        int idxa = -1; // index over only associating compounds
        int idxi = 0; // index for the ii-th compound
        int idxj = 0; // index for the jj-th compound
        int idx_ddelta = -1; // index for ddelta_dd vector
        double dghsd_dd;
        for (int i = 0; i < ncA; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < ncA; j++) {
                idxa += 1;
                idxj = iA[j]*ncomp+iA[j];
                eABij[idxa] = (components[iA[i]].getUAB()+components[iA[j]].getUAB())/2.;
                volABij[idxa] = sqrt(components[iA[i]].getVolA()*components[iA[j]].getVolA())*pow(sqrt(s_ij[idxi]*
                    s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                delta_ij[idxa] = ghs[iA[j]]*(exp(eABij[idxa]/t)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij[idxa];
                for (int k = 0; k < ncomp; k++) {
                    idx_ddelta += 1;
                    dghsd_dd = PI/6.*components[k].getM()*(pow(d[k], 3)/(1-zeta[3])/(1-zeta[3]) + 3*d[iA[i]]*d[iA[j]]/
                        (d[iA[i]]+d[iA[j]])*(d[k]*d[k]/(1-zeta[3])/(1-zeta[3])+2*pow(d[k], 3)*
                        zeta[2]/pow(1-zeta[3], 3)) + 2*pow((d[iA[i]]*d[iA[j]]/(d[iA[i]]+d[iA[j]])), 2)*
                        (2*d[k]*d[k]*zeta[2]/pow(1-zeta[3], 3)+3*(pow(d[k], 3)*zeta[2]*zeta[2]
                        /pow(1-zeta[3], 4))));
                    ddelta_dd[idx_ddelta] = dghsd_dd*(exp(eABij[idxa]/t)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij[idxa];
                }
            }
            XA[i*2] = (-1 + sqrt(1+8*den*delta_ij[i*ncA+i]))/(4*den*delta_ij[i*ncA+i]);
            XA[i*2+1] = XA[i*2];
        }

        vector<double> x_assoc(ncA); // mole fractions of only the associating compounds
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
            for (int i = 0; i < ncA*2; i++) {
                dif += abs(XA[i] - XA_old[i]);
            }
            XA_old = XA;
        }

        vector<double> dXA_dd(a_sites*ncA*ncomp, 0);
        dXA_dd = dXA_find(ncA, ncomp, iA, delta_ij, den, XA, ddelta_dd, x_assoc, a_sites);

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncA; j++) {
                for (int k = 0; k < a_sites; k++) {
                    summ += mole_fractions[i]*den*mole_fractions[iA[j]]*(1/XA[j*a_sites+k]-0.5)*dXA_dd[i*(ncA*a_sites)+j*(a_sites)+k];
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
            q[i] = components[i].getZ()*E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += pow(components[i].getZ(),2.)*mole_fractions[i];
        }

        double kappa = sqrt(den*E_CHRG*E_CHRG/kb/t/(dielc*perm_vac)*summ); // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            double chi, sigma_k;
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi = 3/pow(kappa*components[i].getSigma(), 3)*(1.5 +
                    log(1+kappa*components[i].getSigma()) - 2*(1+kappa*components[i].getSigma()) +
                    0.5*pow(1+kappa*components[i].getSigma(), 2));
                sigma_k = -2*chi+3/(1+kappa*components[i].getSigma());
                summ += q[i]*q[i]*mole_fractions[i]*sigma_k;
            }
            Zion = -1*kappa/24./PI/kb/t/(dielc*perm_vac)*summ;
        }
    }

    double Z = Zid + Zhc + Zdisp + Zpolar + Zassoc + Zion;
    return Z;
}

void PCSAFTBackend::post_update(bool optional_checks){
    // Check the values that must always be set
    //if (_p < 0){ throw ValueError("p is less than zero");}
    if (!ValidNumber(_p)){
        throw ValueError("p is not a valid number");}
    //if (_T < 0){ throw ValueError("T is less than zero");}
    if (!ValidNumber(_T)){ throw ValueError("T is not a valid number");}
    if (_rhomolar < 0){ throw ValueError("rhomolar is less than zero");}
    if (!ValidNumber(_rhomolar)){ throw ValueError("rhomolar is not a valid number");}

    if (optional_checks){
        if (!ValidNumber(_Q)){ throw ValueError("Q is not a valid number");}
        if (_phase == iphase_unknown){
                throw ValueError("_phase is unknown");
        }
    }
}

void PCSAFTBackend::update(CoolProp::input_pairs input_pair, double value1, double value2){
    if (get_debug_level() > 10){std::cout << format("%s (%d): update called with (%d: (%s), %g, %g)",__FILE__,__LINE__, input_pair, get_input_pair_short_desc(input_pair).c_str(), value1, value2) << std::endl;}

    // Converting input to CoolPropDbl
    CoolPropDbl ld_value1 = value1, ld_value2 = value2;
    value1 = ld_value1; value2 = ld_value2;

    // Clear the state
    clear();

    if (is_pure_or_pseudopure == false && mole_fractions.size() == 0) {
        throw ValueError("Mole fractions must be set");
    }

    // If the inputs are in mass units, convert them to molar units
    mass_to_molar_inputs(input_pair, value1, value2);

    switch(input_pair)
    {
        case PT_INPUTS:
            dielc = dielc_water(value2); // Right now only aqueous mixtures are supported. Other solvents could be modeled by replacing the dielc_water function.
            _p = value1; _T = value2; _rhomolar = solver_rho_Tp(value2/*T*/, value1/*p*/, -1.0/*rho_guess*/); break;
        case QT_INPUTS:
            dielc = dielc_water(value2);
            _Q = value1; _T = value2; flash_QT(*this); break;
        case PQ_INPUTS:
            _p = value1; _Q = value2; flash_PQ(*this); break;
        case DmolarT_INPUTS:
            dielc = dielc_water(value2);
            _rhomolar = value1; _T = value2; update_DmolarT(); break;
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

    post_update();
}

/// Set binary mixture floating point parameter for this instance
// void PCSAFTBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter, const double value){
//     if (parameter == "Fij"){
//         residual_helmholtz->Excess.F[i][j] = value;
//         residual_helmholtz->Excess.F[j][i] = value;
//     }
//     else{
//         Reducing->set_binary_interaction_double(i,j,parameter,value);
//     }
//     /// Also set the parameters in the managed pointers for other states
//     for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend> >::iterator it = linked_states.begin(); it != linked_states.end(); ++it){
//         it->get()->set_binary_interaction_double(i, j, parameter, value);
//     }
// }
//
// /// Get binary mixture floating point parameter for this instance
// double PCSAFTBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter){
//     if (parameter == "Fij"){
//         return residual_helmholtz->Excess.F[i][j];
//     }
//     else{
//         return Reducing->get_binary_interaction_double(i,j,parameter);
//     }
// }

double PCSAFTBackend::calc_SatLiquid(bool is_temperature_input, CoolPropDbl value1) {
}

double PCSAFTBackend::calc_SatVapor(bool is_temperature_input, CoolPropDbl value1) {
}

void PCSAFTBackend::flash_QT(PCSAFTBackend &PCSAFT){
}

void PCSAFTBackend::flash_PQ(PCSAFTBackend &PCSAFT){
    //!!! remember to update dielc constant
}

CoolPropDbl PCSAFTBackend::solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_guess){
    phases phase;

    // SolverTPResid resid(this,T,p);
    //
    // // Check if the phase is imposed
    // if (imposed_phase_index != iphase_not_imposed)
    //     // Use the imposed phase index
    //     phase = imposed_phase_index;
    // else
    //     // Use the phase index in the class
    //     phase = _phase;
    //
    // if (rhomolar_guess < 0) // Not provided
    // {
    //     // A gas-like phase, ideal gas might not be the perfect model, but probably good enough
    //     if (phase == iphase_gas || phase == iphase_supercritical_gas || phase == iphase_supercritical)
    //     {
    //         if (rhomolar_guess < 0 || !ValidNumber(rhomolar_guess)) // If the guess is bad, probably high temperature, use ideal gas
    //         {
    //             rhomolar_guess = p/(gas_constant()*T);
    //         }
    //     }
    //     else if (phase == iphase_liquid)
    //     {
    //         double rhomolar;
    //         if (is_pure_or_pseudopure){
    //             mw = calc_molar_mass();
    //             rhomolar_guess = 1100000.0/mw; // initial guess for liquid chosen to be 1100000 g m^-3
    //             try{
    //                 // First we try with Halley's method starting at saturated liquid
    //                 rhomolar = Halley(resid, rhomolar_guess, 1e-8, 100);
    //                 if (!ValidNumber(rhomolar) || first_partial_deriv(iP, iDmolar, iT) < 0 || second_partial_deriv(iP, iDmolar, iT, iDmolar, iT) < 0){
    //                     throw ValueError("Liquid density is invalid");
    //                 }
    //             }
    //             catch(std::exception &){
    //                 // Next we try with a Brent method bounded solver since the function is 1-1
    //                 rhomolar = Brent(resid, rhomolar_guess*0.5, rhomolar_guess*1.4, DBL_EPSILON,1e-8,100);
    //                 if (!ValidNumber(rhomolar)){throw ValueError();}
    //             }
    //         }
    //         else{
    //             // Try with 4th order Householder method starting at a very high density
    //             rhomolar = Householder4(&resid, 3*rhomolar_guess, 1e-8, 100);
    //         }
    //         return rhomolar;
    //     }
    //     else if (phase == iphase_supercritical_liquid){
    //         mw = calc_molar_mass();
    //         rhomolar_guess = 1100000.0/mw; // initial guess chosen to be 1100000 g m^-3
    //
    //         double rhomolar = Brent(resid, rhomolar_guess*0.5, rhomolar_guess*2, DBL_EPSILON,1e-8,100);
    //         if (!ValidNumber(rhomolar)){throw ValueError();}
    //         return rhomolar;
    //     }
    // }
    //
    // try{
    //     double rhomolar = rhomolar = Householder4(resid, rhomolar_guess, 1e-8, 20);
    //     if (!ValidNumber(rhomolar) || rhomolar < 0) {
    //         throw ValueError();
    //     }
    //     if (phase == iphase_liquid){
    //         double dpdrho = first_partial_deriv(iP, iDmolar, iT);
    //         double d2pdrho2 = second_partial_deriv(iP, iDmolar, iT, iDmolar, iT);
    //         if(dpdrho < 0 || d2pdrho2 < 0){
    //             // Try again with a larger density in order to end up at the right solution
    //             rhomolar = Householder4(resid, 3*rhomolar_guess, 1e-8, 100);
    //             return rhomolar;
    //         }
    //     }
    //     else if (phase == iphase_gas){
    //         double dpdrho = first_partial_deriv(iP, iDmolar, iT);
    //         double d2pdrho2 = second_partial_deriv(iP, iDmolar, iT, iDmolar, iT);
    //         if(dpdrho < 0 || d2pdrho2 > 0){
    //             // Try again with a tiny density in order to end up at the right solution
    //             rhomolar = Householder4(resid, 1e-6, 1e-8, 100);
    //             return rhomolar;
    //         }
    //     }
    //     return rhomolar;
    // }
    // catch(std::exception &e)
    // {
    //     if (phase == iphase_supercritical || phase == iphase_supercritical_gas){
    //         mw = calc_molar_mass();
    //         rhomolar_guess = 350000.0/mw; // initial guess chosen to be 350000 g m^-3
    //         double rhomolar = Brent(resid, 1e-10, 3*rhomolar_guess, DBL_EPSILON, 1e-8, 100);
    //         return rhomolar;
    //     }
    //     throw ValueError(format("solver_rho_Tp was unable to find a solution for T=%10Lg, p=%10Lg, with guess value %10Lg with error: %s",T,p,rhomolar_guess, e.what()));
    // }
    return 1100000.0;
}

CoolPropDbl PCSAFTBackend::calc_molar_mass(void){
    double summer = 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        summer += mole_fractions[i] * components[i].molar_mass();
    }
    return summer;
}

vector<double> PCSAFTBackend::XA_find(vector<double> XA_guess, int ncomp, vector<double> delta_ij, double den,
    vector<double> x) {
    /**Iterate over this function in order to solve for XA*/
    int n_sites = XA_guess.size()/ncomp;
    double summ2;
    vector<CoolPropDbl> XA = XA_guess;

    for (int i = 0; i < ncomp; i++) {
        for (int kout = 0; kout < n_sites; kout++) {
            summ2 = 0.;
            for (int j = 0; j < ncomp; j++) {
                for (int kin = 0; kin < n_sites; kin++) {
                    if (kin != kout) {
                        summ2 += den*mole_fractions[j]*XA_guess[j*ncomp+kin]*delta_ij[i*ncomp+j];
                    }
                }
            }
            XA[i*ncomp+kout] = 1./(1.+summ2);
        }
    }

    return XA;
    }

vector<double> PCSAFTBackend::dXA_find(int ncA, int ncomp, vector<int> iA, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dd, vector<double> x, int n_sites) {
    /**Solve for the derivative of XA with respect to density.*/
    Eigen::MatrixXd B(n_sites*ncA*ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_sites*ncA*ncomp, n_sites*ncA*ncomp);

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
                        sum1 = sum1 + den*mole_fractions[k]*(XA[indx2]*ddelta_dd[j*(ncA*ncomp)+k*(ncomp)+i]*((indx1+indx2)%2)); // (indx1+indx2)%2 ensures that A-A and B-B associations are set to zero
                        A(indx1+i*n_sites*ncA,indx2+i*n_sites*ncA) =
                        A(indx1+i*n_sites*ncA,indx2+i*n_sites*ncA) +
                        XA[indx1]*XA[indx1]*den*mole_fractions[k]*delta_ij[j*ncA+k]*((indx1+indx2)%2);
                    }
                }

                sum2 = 0;
                if (find(iA.begin(), iA.end(), i) != iA.end()) {
                    for (int k = 0; k < n_sites; k++) {
                        sum2 = sum2 + XA[n_sites*(indx4)+k]*delta_ij[indx4*ncA+j]*((indx1+k)%2);
                    }
                }

                A(indx3,indx3) = A(indx3,indx3) + 1;
                B(indx3) = -1*XA[indx1]*XA[indx1]*(sum1 + sum2);
            }
        }
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dd(n_sites*ncA*ncomp);
    for (int i = 0; i < n_sites*ncA*ncomp; i++) {
        dXA_dd[i] = solution(i);
    }
    return dXA_dd;
}


vector<double> PCSAFTBackend::dXAdt_find(int ncA, vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x, int n_sites) {
    /**Solve for the derivative of XA with respect to temperature.*/
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n_sites*ncA, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_sites*ncA, n_sites*ncA);

    double summ;
    int i_in, i_out = -1; // i_out is index of outer iteration loop (follows row of matrices)
    for (int i = 0; i < ncA; i++) {
        for (int ai = 0; ai < n_sites; ai++) {
            i_out += 1;
            i_in = -1; // index for summation loops
            summ = 0;
            for (int j = 0; j < ncA; j++) {
                for (int bj = 0; bj < n_sites; bj++) {
                    i_in += 1;
                    B(i_out) -= mole_fractions[j]*XA[i_in]*ddelta_dt[i*ncA+j]*((i_in+i_out)%2); // (i_in+i_out)%2 ensures that A-A and B-B associations are set to zero
                    A(i_out,i_in) = mole_fractions[j]*delta_ij[i*ncA+j]*((i_in+i_out)%2);
                    summ += mole_fractions[j]*XA[i_in]*delta_ij[i*ncA+j]*((i_in+i_out)%2);
                }
            }
            A(i_out,i_out) = A(i_out,i_out) + pow(1+den*summ, 2.)/den;
        }
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dt(n_sites*ncA);
    for (int i = 0; i < n_sites*ncA; i++) {
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
    }
    else if (t <= 368.15) {
        dielc = 7.6555618295E-04*t*t - 8.1783881423E-01*t + 2.5419616803E+02;
    }
    else if (t <= 443.15) {
        dielc = 0.0005003272124*t*t - 0.6285556029*t + 220.4467027;
    }
    else {
        throw ValueError("The current function for the dielectric constant for water is only valid for temperatures less than 443.15 K.");
    }
    return dielc;
}

} /* namespace CoolProp */
