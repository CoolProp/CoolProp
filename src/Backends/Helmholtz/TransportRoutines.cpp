
#include "TransportRoutines.h"
#include "CoolPropFluid.h"

namespace CoolProp {

CoolPropDbl TransportRoutines::viscosity_dilute_kinetic_theory(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        CoolPropDbl Tstar = HEOS.T() / HEOS.components[0].transport.epsilon_over_k;
        CoolPropDbl sigma_nm = HEOS.components[0].transport.sigma_eta * 1e9;  // 1e9 to convert from m to nm
        CoolPropDbl molar_mass_kgkmol = HEOS.molar_mass() * 1000;             // 1000 to convert from kg/mol to kg/kmol

        // The nondimensional empirical collision integral from Neufeld
        // Neufeld, P. D.; Janzen, A. R.; Aziz, R. A. Empirical Equations to Calculate 16 of the Transport Collision Integrals (l,s)*
        // for the Lennard-Jones (12-6) Potential. J. Chem. Phys. 1972, 57, 1100-1102
        CoolPropDbl OMEGA22 =
          1.16145 * pow(Tstar, static_cast<CoolPropDbl>(-0.14874)) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar);

        // The dilute gas component -
        return 26.692e-9 * sqrt(molar_mass_kgkmol * HEOS.T()) / (pow(sigma_nm, 2) * OMEGA22);  // Pa-s
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_dilute_kinetic_theory is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_dilute_collision_integral(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteGasCollisionIntegralData& data = HEOS.components[0].transport.viscosity_dilute.collision_integral;
        const std::vector<CoolPropDbl>&a = data.a, &t = data.t;
        const CoolPropDbl C = data.C, molar_mass = data.molar_mass;

        CoolPropDbl S;
        // Unit conversions and variable definitions
        const CoolPropDbl Tstar = HEOS.T() / HEOS.components[0].transport.epsilon_over_k;
        const CoolPropDbl sigma_nm = HEOS.components[0].transport.sigma_eta * 1e9;  // 1e9 to convert from m to nm
        const CoolPropDbl molar_mass_kgkmol = molar_mass * 1000;                    // 1000 to convert from kg/mol to kg/kmol

        /// Both the collision integral \f$\mathfrak{S}^*\f$ and effective cross section \f$\Omega^{(2,2)}\f$ have the same form,
        /// in general we don't care which is used.  The are related through \f$\Omega^{(2,2)} = (5/4)\mathfrak{S}^*\f$
        /// see Vesovic(JPCRD, 1990) for CO\f$_2\f$ for further information
        CoolPropDbl summer = 0, lnTstar = log(Tstar);
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * pow(lnTstar, t[i]);
        }
        S = exp(summer);

        // The dilute gas component
        return C * sqrt(molar_mass_kgkmol * HEOS.T()) / (pow(sigma_nm, 2) * S);  // Pa-s
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_dilute_collision_integral is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_dilute_powers_of_T(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteGasPowersOfT& data = HEOS.components[0].transport.viscosity_dilute.powers_of_T;
        const std::vector<CoolPropDbl>&a = data.a, &t = data.t;

        CoolPropDbl summer = 0, T = HEOS.T();
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * pow(T, t[i]);
        }
        return summer;
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_dilute_powers_of_T is only for pure and pseudo-pure");
    }
}
CoolPropDbl TransportRoutines::viscosity_dilute_powers_of_Tr(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteGasPowersOfTr& data = HEOS.components[0].transport.viscosity_dilute.powers_of_Tr;
        const std::vector<CoolPropDbl>&a = data.a, &t = data.t;
        CoolPropDbl summer = 0, Tr = HEOS.T() / data.T_reducing;
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * pow(Tr, t[i]);
        }
        return summer;
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_dilute_powers_of_Tr is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_dilute_collision_integral_powers_of_T(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteCollisionIntegralPowersOfTstarData& data =
          HEOS.components[0].transport.viscosity_dilute.collision_integral_powers_of_Tstar;
        const std::vector<CoolPropDbl>&a = data.a, &t = data.t;

        CoolPropDbl summer = 0, Tstar = HEOS.T() / data.T_reducing;
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * pow(Tstar, t[i]);
        }
        return data.C * sqrt(HEOS.T()) / summer;
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_dilute_collision_integral_powers_of_T is only for pure and pseudo-pure");
    }
}
CoolPropDbl TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        CoolProp::ViscosityModifiedBatschinskiHildebrandData& HO =
          HEOS.components[0].transport.viscosity_higher_order.modified_Batschinski_Hildebrand;

        CoolPropDbl delta = HEOS.rhomolar() / HO.rhomolar_reduce, tau = HO.T_reduce / HEOS.T();

        // The first term that is formed of powers of tau (Tc/T) and delta (rho/rhoc)
        CoolPropDbl S = 0;
        for (unsigned int i = 0; i < HO.a.size(); ++i) {
            S += HO.a[i] * pow(delta, HO.d1[i]) * pow(tau, HO.t1[i]) * exp(HO.gamma[i] * pow(delta, HO.l[i]));
        }

        // For the terms that multiplies the bracketed term with delta and delta0
        CoolPropDbl F = 0;
        for (unsigned int i = 0; i < HO.f.size(); ++i) {
            F += HO.f[i] * pow(delta, HO.d2[i]) * pow(tau, HO.t2[i]);
        }

        // for delta_0
        CoolPropDbl summer_numer = 0;
        for (unsigned int i = 0; i < HO.g.size(); ++i) {
            summer_numer += HO.g[i] * pow(tau, HO.h[i]);
        }
        CoolPropDbl summer_denom = 0;
        for (unsigned int i = 0; i < HO.p.size(); ++i) {
            summer_denom += HO.p[i] * pow(tau, HO.q[i]);
        }
        CoolPropDbl delta0 = summer_numer / summer_denom;

        // The higher-order-term component
        return S + F * (1 / (delta0 - delta) - 1 / delta0);  // Pa-s
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityRainWaterFriendData& data = HEOS.components[0].transport.viscosity_initial.rainwater_friend;
        const std::vector<CoolPropDbl>&b = data.b, &t = data.t;

        CoolPropDbl B_eta, B_eta_star;
        CoolPropDbl Tstar = HEOS.T() / HEOS.components[0].transport.epsilon_over_k;  // [no units]
        CoolPropDbl sigma = HEOS.components[0].transport.sigma_eta;                  // [m]

        CoolPropDbl summer = 0;
        for (unsigned int i = 0; i < b.size(); ++i) {
            summer += b[i] * pow(Tstar, t[i]);
        }
        B_eta_star = summer;                                 // [no units]
        B_eta = 6.02214129e23 * pow(sigma, 3) * B_eta_star;  // [m^3/mol]
        return B_eta;                                        // [m^3/mol]
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_initial_density_dependence_empirical(HelmholtzEOSMixtureBackend& HEOS) {
    // Inspired by the form from Tariq, JPCRD, 2014
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ViscosityInitialDensityEmpiricalData& data = HEOS.components[0].transport.viscosity_initial.empirical;
        const std::vector<CoolPropDbl>&n = data.n, &d = data.d, &t = data.t;

        CoolPropDbl tau = data.T_reducing / HEOS.T();                  // [no units]
        CoolPropDbl delta = HEOS.rhomolar() / data.rhomolar_reducing;  // [no units]

        CoolPropDbl summer = 0;
        for (unsigned int i = 0; i < n.size(); ++i) {
            summer += n[i] * pow(delta, d[i]) * pow(tau, t[i]);
        }
        return summer;  // [Pa-s]
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_initial_density_dependence_empirical is only for pure and pseudo-pure");
    }
}

static void visc_Helper(double Tbar, double rhobar, double* mubar_0, double* mubar_1) {
    std::vector<std::vector<CoolPropDbl>> H(6, std::vector<CoolPropDbl>(7, 0));
    double sum;
    int i, j;

    // Dilute-gas component
    *mubar_0 = 100.0 * sqrt(Tbar) / (1.67752 + 2.20462 / Tbar + 0.6366564 / powInt(Tbar, 2) - 0.241605 / powInt(Tbar, 3));

    //Fill in zeros in H
    for (i = 0; i <= 5; i++) {
        for (j = 0; j <= 6; j++) {
            H[i][j] = 0;
        }
    }

    //Set non-zero parameters of H
    H[0][0] = 5.20094e-1;
    H[1][0] = 8.50895e-2;
    H[2][0] = -1.08374;
    H[3][0] = -2.89555e-1;

    H[0][1] = 2.22531e-1;
    H[1][1] = 9.99115e-1;
    H[2][1] = 1.88797;
    H[3][1] = 1.26613;
    H[5][1] = 1.20573e-1;

    H[0][2] = -2.81378e-1;
    H[1][2] = -9.06851e-1;
    H[2][2] = -7.72479e-1;
    H[3][2] = -4.89837e-1;
    H[4][2] = -2.57040e-1;

    H[0][3] = 1.61913e-1;
    H[1][3] = 2.57399e-1;

    H[0][4] = -3.25372e-2;
    H[3][4] = 6.98452e-2;

    H[4][5] = 8.72102e-3;

    H[3][6] = -4.35673e-3;
    H[5][6] = -5.93264e-4;

    // Finite density component
    sum = 0;
    for (i = 0; i <= 5; i++) {
        for (j = 0; j <= 6; j++) {
            sum += powInt(1 / Tbar - 1, i) * (H[i][j] * powInt(rhobar - 1, j));
        }
    }
    *mubar_1 = exp(rhobar * sum);
}
CoolPropDbl TransportRoutines::viscosity_heavywater_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    double Tbar = HEOS.T() / 643.847, rhobar = HEOS.rhomass() / 358;
    double A[] = {1.000000, 0.940695, 0.578377, -0.202044};
    int I[] = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3, 5, 0, 1, 5, 3};
    int J[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
    double Bij[] = {0.4864192,  -0.2448372,  -0.8702035, 0.8716056,    -1.051126,   0.3458395,    0.3509007,   1.315436,    1.297752,
                    1.353448,   -0.2847572,  -1.037026,  -1.287846,    -0.02148229, 0.07013759,   0.4660127,   0.2292075,   -0.4857462,
                    0.01641220, -0.02884911, 0.1607171,  -0.009603846, -0.01163815, -0.008239587, 0.004559914, -0.003886659};
    double mu0 = sqrt(Tbar) / (A[0] + A[1] / Tbar + A[2] / POW2(Tbar) + A[3] / POW3(Tbar));
    double summer = 0;
    for (int i = 0; i < 26; ++i) {
        summer += Bij[i] * pow(1 / Tbar - 1, I[i]) * pow(rhobar - 1, J[i]);
    }
    double mu1 = exp(rhobar * summer);
    double mubar = mu0 * mu1;
    return 55.2651e-6 * mubar;
}
CoolPropDbl TransportRoutines::viscosity_water_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    double x_mu = 0.068, qc = 1 / 1.9, qd = 1 / 1.1, nu = 0.630, gamma = 1.239, zeta_0 = 0.13, LAMBDA_0 = 0.06, Tbar_R = 1.5, pstar, Tstar, rhostar;
    double delta, tau, mubar_0, mubar_1, mubar_2, drhodp, drhodp_R, DeltaChibar, zeta, w, L, Y, psi_D, Tbar, rhobar;
    double drhobar_dpbar, drhobar_dpbar_R, R_Water;

    pstar = 22.064e6;  // [Pa]
    Tstar = 647.096;   // [K]
    rhostar = 322;     // [kg/m^3]
    Tbar = HEOS.T() / Tstar;
    rhobar = HEOS.rhomass() / rhostar;
    R_Water = HEOS.gas_constant() / HEOS.molar_mass();  // [J/kg/K]

    // Dilute and finite gas portions
    visc_Helper(Tbar, rhobar, &mubar_0, &mubar_1);

    // **********************************************************************
    // ************************ Critical Enhancement ************************
    // **********************************************************************
    delta = rhobar;
    // "Normal" calculation
    drhodp = 1 / (R_Water * HEOS.T() * (1 + 2 * delta * HEOS.dalphar_dDelta() + delta * delta * HEOS.d2alphar_dDelta2()));
    drhobar_dpbar = pstar / rhostar * drhodp;
    // "Reducing" calculation
    tau = 1 / Tbar_R;
    drhodp_R = 1
               / (R_Water * Tbar_R * Tstar
                  * (1 + 2 * rhobar * HEOS.calc_alphar_deriv_nocache(0, 1, HEOS.mole_fractions, tau, delta)
                     + delta * delta * HEOS.calc_alphar_deriv_nocache(0, 2, HEOS.mole_fractions, tau, delta)));
    drhobar_dpbar_R = pstar / rhostar * drhodp_R;

    DeltaChibar = rhobar * (drhobar_dpbar - drhobar_dpbar_R * Tbar_R / Tbar);
    if (DeltaChibar < 0) DeltaChibar = 0;
    zeta = zeta_0 * pow(DeltaChibar / LAMBDA_0, nu / gamma);
    if (zeta < 0.3817016416) {
        Y = 1.0 / 5.0 * qc * zeta * powInt(qd * zeta, 5) * (1 - qc * zeta + powInt(qc * zeta, 2) - 765.0 / 504.0 * powInt(qd * zeta, 2));
    } else {
        psi_D = acos(pow(1 + powInt(qd * zeta, 2), -1.0 / 2.0));
        w = sqrt(std::abs((qc * zeta - 1) / (qc * zeta + 1))) * tan(psi_D / 2.0);
        if (qc * zeta > 1) {
            L = log((1 + w) / (1 - w));
        } else {
            L = 2 * atan(std::abs(w));
        }
        Y = 1.0 / 12.0 * sin(3 * psi_D) - 1 / (4 * qc * zeta) * sin(2 * psi_D)
            + 1.0 / powInt(qc * zeta, 2) * (1 - 5.0 / 4.0 * powInt(qc * zeta, 2)) * sin(psi_D)
            - 1.0 / powInt(qc * zeta, 3) * ((1 - 3.0 / 2.0 * powInt(qc * zeta, 2)) * psi_D - pow(std::abs(powInt(qc * zeta, 2) - 1), 3.0 / 2.0) * L);
    }
    mubar_2 = exp(x_mu * Y);

    return (mubar_0 * mubar_1 * mubar_2) / 1e6;
}
CoolPropDbl TransportRoutines::viscosity_toluene_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl Tr = HEOS.T() / 591.75, rhor = HEOS.keyed_output(CoolProp::iDmass) / 291.987;
    CoolPropDbl c[] = {19.919216, -2.6557905, -135.904211, -7.9962719, -11.014795, -10.113817};
    return 1e-6 * pow(static_cast<double>(rhor), 2.0 / 3.0) * sqrt(Tr)
           * ((c[0] * rhor + c[1] * pow(rhor, 4)) / Tr + c[2] * rhor * rhor * rhor / (rhor * rhor + c[3] + c[4] * Tr) + c[5] * rhor);
}

CoolPropDbl TransportRoutines::viscosity_hydrogen_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl Tr = HEOS.T() / 33.145, rhor = HEOS.keyed_output(CoolProp::iDmass) * 0.011;
    CoolPropDbl c[] = {0, 6.43449673e-6, 4.56334068e-2, 2.32797868e-1, 9.58326120e-1, 1.27941189e-1, 3.63576595e-1};
    return c[1] * pow(rhor, 2) * exp(c[2] * Tr + c[3] / Tr + c[4] * pow(rhor, 2) / (c[5] + Tr) + c[6] * pow(rhor, 6));
}
CoolPropDbl TransportRoutines::viscosity_benzene_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl Tr = HEOS.T() / 562.02, rhor = HEOS.rhomass() / 304.792;
    CoolPropDbl c[] = {-9.98945, 86.06260, 2.74872, 1.11130, -1.0, -134.1330, -352.473, 6.60989, 88.4174};
    return 1e-6 * pow(rhor, static_cast<CoolPropDbl>(2.0 / 3.0)) * sqrt(Tr)
           * (c[0] * pow(rhor, 2) + c[1] * rhor / (c[2] + c[3] * Tr + c[4] * rhor)
              + (c[5] * rhor + c[6] * pow(rhor, 2)) / (c[7] + c[8] * pow(rhor, 2)));
}
CoolPropDbl TransportRoutines::viscosity_hexane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {

    CoolPropDbl Tr = HEOS.T() / 507.82, rhor = HEOS.keyed_output(CoolProp::iDmass) / 233.182;

    // Output is in Pa-s
    double c[] = {2.53402335 / 1e6, -9.724061002 / 1e6, 0.469437316, 158.5571631, 72.42916856 / 1e6,
                  10.60751253,      8.628373915,        -6.61346441, -2.212724566};
    return pow(rhor, static_cast<CoolPropDbl>(2.0 / 3.0)) * sqrt(Tr)
           * (c[0] / Tr + c[1] / (c[2] + Tr + c[3] * rhor * rhor)
              + c[4] * (1 + rhor) / (c[5] + c[6] * Tr + c[7] * rhor + rhor * rhor + c[8] * rhor * Tr));
}

CoolPropDbl TransportRoutines::viscosity_heptane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    /// From Michailidou-JPCRD-2014-Heptane
    CoolPropDbl Tr = HEOS.T() / 540.13, rhor = HEOS.rhomass() / 232;

    // Output is in Pa-s
    double c[] = {0, 22.15000 / 1e6, -15.00870 / 1e6, 3.71791 / 1e6, 77.72818 / 1e6, 9.73449, 9.51900, -6.34076, -2.51909};
    return pow(rhor, static_cast<CoolPropDbl>(2.0L / 3.0L)) * sqrt(Tr)
           * (c[1] * rhor + c[2] * pow(rhor, 2) + c[3] * pow(rhor, 3)
              + c[4] * rhor / (c[5] + c[6] * Tr + c[7] * rhor + rhor * rhor + c[8] * rhor * Tr));
}

CoolPropDbl TransportRoutines::viscosity_higher_order_friction_theory(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        CoolProp::ViscosityFrictionTheoryData& F = HEOS.components[0].transport.viscosity_higher_order.friction_theory;

        CoolPropDbl tau = F.T_reduce / HEOS.T(), kii = 0, krrr = 0, kaaa = 0, krr, kdrdr;

        double psi1 = exp(tau) - F.c1;
        double psi2 = exp(pow(tau, 2)) - F.c2;

        double ki = (F.Ai[0] + F.Ai[1] * psi1 + F.Ai[2] * psi2) * tau;

        double ka = (F.Aa[0] + F.Aa[1] * psi1 + F.Aa[2] * psi2) * pow(tau, F.Na);
        double kr = (F.Ar[0] + F.Ar[1] * psi1 + F.Ar[2] * psi2) * pow(tau, F.Nr);
        double kaa = (F.Aaa[0] + F.Aaa[1] * psi1 + F.Aaa[2] * psi2) * pow(tau, F.Naa);
        if (F.Arr.empty()) {
            krr = 0;
            kdrdr = (F.Adrdr[0] + F.Adrdr[1] * psi1 + F.Adrdr[2] * psi2) * pow(tau, F.Nrr);
        } else {
            krr = (F.Arr[0] + F.Arr[1] * psi1 + F.Arr[2] * psi2) * pow(tau, F.Nrr);
            kdrdr = 0;
        }
        if (!F.Aii.empty()) {
            kii = (F.Aii[0] + F.Aii[1] * psi1 + F.Aii[2] * psi2) * pow(tau, F.Nii);
        }
        if (!F.Arrr.empty() && !F.Aaaa.empty()) {
            krrr = (F.Arrr[0] + F.Arrr[1] * psi1 + F.Arrr[2] * psi2) * pow(tau, F.Nrrr);
            kaaa = (F.Aaaa[0] + F.Aaaa[1] * psi1 + F.Aaaa[2] * psi2) * pow(tau, F.Naaa);
        }

        double p = HEOS.p() / 1e5;  // [bar]; 1e5 for conversion from Pa -> bar
        double pr =
          HEOS.T() * HEOS.first_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar) / 1e5;  // [bar/K]; 1e5 for conversion from Pa -> bar
        double pa = p - pr;                                                                          //[bar]
        double pid = HEOS.rhomolar() * HEOS.gas_constant() * HEOS.T() / 1e5;                         // [bar]; 1e5 for conversion from Pa -> bar
        double deltapr = pr - pid;

        double eta_f = ka * pa + kr * deltapr + ki * pid + kaa * pa * pa + kdrdr * deltapr * deltapr + krr * pr * pr + kii * pid * pid
                       + krrr * pr * pr * pr + kaaa * pa * pa * pa;

        return eta_f;  //[Pa-s]
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_higher_order_friction_theory is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::viscosity_helium_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    double eta_0, eta_0_slash, eta_E_slash, B, C, D, ln_eta, x;
    //
    // Arp, V.D., McCarty, R.D., and Friend, D.G.,
    // "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa",
    // NIST Technical Note 1334 (revised), 1998.
    //
    // Using Arp NIST report
    // Report is not clear on viscosity, referring to REFPROP source code for clarity

    // Correlation wants density in g/cm^3; kg/m^3 --> g/cm^3, divide by 1000
    CoolPropDbl rho = HEOS.keyed_output(CoolProp::iDmass) / 1000.0, T = HEOS.T();

    if (T <= 300) {
        x = log(T);
    } else {
        x = log(300.0);
    }
    // Evaluate the terms B,C,D
    B = -47.5295259 / x + 87.6799309 - 42.0741589 * x + 8.33128289 * x * x - 0.589252385 * x * x * x;
    C = 547.309267 / x - 904.870586 + 431.404928 * x - 81.4504854 * x * x + 5.37008433 * x * x * x;
    D = -1684.39324 / x + 3331.08630 - 1632.19172 * x + 308.804413 * x * x - 20.2936367 * x * x * x;
    eta_0_slash = -0.135311743 / x + 1.00347841 + 1.20654649 * x - 0.149564551 * x * x + 0.012520841 * x * x * x;
    eta_E_slash = rho * B + rho * rho * C + rho * rho * rho * D;

    if (T <= 100) {
        ln_eta = eta_0_slash + eta_E_slash;
        // Correlation yields viscosity in micro g/(cm-s); to get Pa-s, divide by 10 to get micro Pa-s, then another 1e6 to get Pa-s
        return exp(ln_eta) / 10.0 / 1e6;
    } else {
        ln_eta = eta_0_slash + eta_E_slash;
        eta_0 = 196 * pow(T, static_cast<CoolPropDbl>(0.71938)) * exp(12.451 / T - 295.67 / T / T - 4.1249);
        // Correlation yields viscosity in micro g/(cm-s); to get Pa-s, divide by 10 to get micro Pa-s, then another 1e6 to get Pa-s
        return (exp(ln_eta) + eta_0 - exp(eta_0_slash)) / 10.0 / 1e6;
    }
}

CoolPropDbl TransportRoutines::viscosity_methanol_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl B_eta, C_eta, epsilon_over_k = 577.87, /* [K]*/
      sigma0 = 0.3408e-9,                              /* [m] */
      delta = 0.4575,                                  /* NOT the reduced density, that is rhor here*/
      N_A = 6.02214129e23, M = 32.04216,               /* kg/kmol */
      T = HEOS.T();
    CoolPropDbl rhomolar = HEOS.rhomolar();

    CoolPropDbl B_eta_star, C_eta_star;
    CoolPropDbl Tstar = T / epsilon_over_k;  // [no units]
    CoolPropDbl rhor = HEOS.rhomass() / 273;
    CoolPropDbl Tr = T / 512.6;

    // Rainwater-Friend initial density terms
    {  // Scoped here so that we can re-use the b variable
        CoolPropDbl b[9] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
        CoolPropDbl t[9] = {0, -0.25, -0.5, -0.75, -1.0, -1.25, -1.5, -2.5, -5.5};
        CoolPropDbl summer = 0;
        for (unsigned int i = 0; i < 9; ++i) {
            summer += b[i] * pow(Tstar, t[i]);
        }
        B_eta_star = summer;                        // [no units]
        B_eta = N_A * pow(sigma0, 3) * B_eta_star;  // [m^3/mol]

        CoolPropDbl c[2] = {1.86222085e-3, 9.990338};
        C_eta_star = c[0] * pow(Tstar, 3) * exp(c[1] * pow(Tstar, static_cast<CoolPropDbl>(-0.5)));  // [no units]
        C_eta = pow(N_A * pow(sigma0, 3), 2) * C_eta_star;                                           // [m^6/mol^2]
    }

    CoolPropDbl eta_g = 1 + B_eta * rhomolar + C_eta * rhomolar * rhomolar;
    CoolPropDbl a[13] = {1.16145, -0.14874, 0.52487, -0.77320, 2.16178,  -2.43787, 0.95976e-3,
                         0.10225, -0.97346, 0.10657, -0.34528, -0.44557, -2.58055};
    CoolPropDbl d[7] = {-1.181909, 0.5031030, -0.6268461, 0.5169312, -0.2351349, 5.3980235e-2, -4.9069617e-3};
    CoolPropDbl e[10] = {0, 4.018368, -4.239180, 2.245110, -0.5750698, 2.3021026e-2, 2.5696775e-2, -6.8372749e-3, 7.2707189e-4, -2.9255711e-5};

    CoolPropDbl OMEGA_22_star_LJ = a[0] * pow(Tstar, a[1]) + a[2] * exp(a[3] * Tstar) + a[4] * exp(a[5] * Tstar);
    CoolPropDbl OMEGA_22_star_delta = a[7] * pow(Tstar, a[8]) + a[9] * exp(a[10] * Tstar) + a[11] * exp(a[12] * Tstar);
    CoolPropDbl OMEGA_22_star_SM = OMEGA_22_star_LJ * (1 + pow(delta, 2) / (1 + a[6] * pow(delta, 6)) * OMEGA_22_star_delta);
    CoolPropDbl eta_0 = 2.66957e-26 * sqrt(M * T) / (pow(sigma0, 2) * OMEGA_22_star_SM);

    CoolPropDbl summerd = 0;
    for (int i = 0; i < 7; ++i) {
        summerd += d[i] / pow(Tr, i);
    }
    for (int j = 1; j < 10; ++j) {
        summerd += e[j] * pow(rhor, j);
    }
    CoolPropDbl sigmac = 0.7193422e-9;                                                                    // [m]
    CoolPropDbl sigma_HS = summerd * sigmac;                                                              // [m]
    CoolPropDbl b = 2 * M_PI * N_A * pow(sigma_HS, 3) / 3;                                                // [m^3/mol]
    CoolPropDbl zeta = b * rhomolar / 4;                                                                  // [-]
    CoolPropDbl g_sigma_HS = (1 - 0.5 * zeta) / pow(1 - zeta, 3);                                         // [-]
    CoolPropDbl eta_E = 1 / g_sigma_HS + 0.8 * b * rhomolar + 0.761 * g_sigma_HS * pow(b * rhomolar, 2);  // [-]

    CoolPropDbl f = 1 / (1 + exp(5 * (rhor - 1)));
    return eta_0 * (f * eta_g + (1 - f) * eta_E);
}

CoolPropDbl TransportRoutines::viscosity_R23_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    double C1 = 1.3163,  //
      C2 = 0.1832, DeltaGstar = 771.23, rhoL = 32.174, rhocbar = 7.5114, Tc = 299.2793, DELTAeta_max = 3.967, Ru = 8.31451, molar_mass = 70.014;

    double a[] = {0.4425728, -0.5138403, 0.1547566, -0.02821844, 0.001578286};
    double e_k = 243.91, sigma = 0.4278;
    double Tstar = HEOS.T() / e_k;
    double logTstar = log(Tstar);
    double Omega = exp(a[0] + a[1] * logTstar + a[2] * pow(logTstar, 2) + a[3] * pow(logTstar, 3) + a[4] * pow(logTstar, 4));
    double eta_DG = 1.25 * 0.021357 * sqrt(molar_mass * HEOS.T()) / (sigma * sigma * Omega);  // uPa-s

    double rhobar = HEOS.rhomolar() / 1000;  // [mol/L]
    double eta_L = C2 * (rhoL * rhoL) / (rhoL - rhobar) * sqrt(HEOS.T()) * exp(rhobar / (rhoL - rhobar) * DeltaGstar / (Ru * HEOS.T()));

    double chi = rhobar - rhocbar;
    double tau = HEOS.T() - Tc;

    double DELTAeta_c = 4 * DELTAeta_max / ((exp(chi) + exp(-chi)) * (exp(tau) + exp(-tau)));

    return (pow((rhoL - rhobar) / rhoL, C1) * eta_DG + pow(rhobar / rhoL, C1) * eta_L + DELTAeta_c) / 1e6;
}

CoolPropDbl TransportRoutines::viscosity_o_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    // From CAO, JPCRD, 2016
    double D[] = {-2.05581e-3, 2.38762, 0, 10.4497, 15.9587};
    double n[] = {10.3, 3.3, 25, 0.7, 0.4};
    double E[] = {2.65651e-3, 0, 1.77616e-12, -18.2446, 0};
    double k[] = {0.8, 0, 4.4};
    double Tr = HEOS.T() / 630.259, rhor = HEOS.rhomolar() / 1000.0 / 2.6845;

    double A0 = -1.4933, B0 = 473.2, C0 = -57033, T = HEOS.T();
    double ln_Seta = A0 + B0 / T + C0 / (T * T);
    double eta0 = 0.22225 * sqrt(T) / exp(ln_Seta);  // [uPa-s]

    double A1 = 13.2814, B1 = -10862.4, C1 = 1664060, rho_molL = HEOS.rhomolar() / 1000.0;
    double eta1 = (A1 + B1 / T + C1 / (T * T)) * rho_molL;  // [uPa-s]

    double f = (D[0] + E[0] * pow(Tr, -k[0])) * pow(rhor, n[0]) + D[1] * pow(rhor, n[1]) + E[2] * pow(rhor, n[2]) / pow(Tr, k[2])
               + (D[3] * rhor + E[3] * Tr) * pow(rhor, n[3]) + D[4] * pow(rhor, n[4]);
    double DELTAeta = pow(rhor, 2.0 / 3.0) * sqrt(Tr) * f;  // [uPa-s]

    return (eta0 + eta1 + DELTAeta) / 1e6;
}
CoolPropDbl TransportRoutines::viscosity_m_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    // From CAO, JPCRD, 2016
    double D[] = {-0.268950, -0.0290018, 0, 14.7728, 17.1128};
    double n[] = {6.8, 3.3, 22.0, 0.6, 0.4};
    double E[] = {0.320971, 0, 1.72866e-10, -18.9852, 0};
    double k[] = {0.3, 0, 3.2};
    double Tr = HEOS.T() / 616.89, rhor = HEOS.rhomolar() / 1000.0 / 2.665;

    double A0 = -1.4933, B0 = 473.2, C0 = -57033, T = HEOS.T();
    double ln_Seta = A0 + B0 / T + C0 / (T * T);
    double eta0 = 0.22115 * sqrt(T) / exp(ln_Seta);  // [uPa-s]

    double A1 = 13.2814, B1 = -10862.4, C1 = 1664060, rho_molL = HEOS.rhomolar() / 1000.0;
    double eta1 = (A1 + B1 / T + C1 / (T * T)) * rho_molL;  // [uPa-s]

    double f = (D[0] + E[0] * pow(Tr, -k[0])) * pow(rhor, n[0]) + D[1] * pow(rhor, n[1]) + E[2] * pow(rhor, n[2]) / pow(Tr, k[2])
               + (D[3] * rhor + E[3] * Tr) * pow(rhor, n[3]) + D[4] * pow(rhor, n[4]);
    double DELTAeta = pow(rhor, 2.0 / 3.0) * sqrt(Tr) * f;  // [uPa-s]

    return (eta0 + eta1 + DELTAeta) / 1e6;  // [Pa-s]
}
CoolPropDbl TransportRoutines::viscosity_p_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    // From Balogun, JPCRD, 2016
    double Tr = HEOS.T() / 616.168, rhor = HEOS.rhomolar() / 1000.0 / 2.69392;

    double A0 = -1.4933, B0 = 473.2, C0 = -57033, T = HEOS.T();
    double ln_Seta = A0 + B0 / T + C0 / (T * T);
    double eta0 = 0.22005 * sqrt(T) / exp(ln_Seta);  // [uPa-s]

    double A1 = 13.2814, B1 = -10862.4, C1 = 1664060, rho_molL = HEOS.rhomolar() / 1000.0;
    double eta1 = (A1 + B1 / T + C1 / (T * T)) * rho_molL;  // [uPa-s]

    double sum1 = 122.919 * pow(rhor, 1.5) - 282.329 * pow(rhor, 2) + 279.348 * pow(rhor, 3) - 146.776 * pow(rhor, 4) + 28.361 * pow(rhor, 5)
                  - 0.004585 * pow(rhor, 11);
    double sum2 = 15.337 * pow(rhor, 1.5) - 0.0004382 * pow(rhor, 11) + 0.00002307 * pow(rhor, 15);
    double DELTAeta = pow(rhor, 2.0 / 3.0) * (sum1 + 1 / sqrt(Tr) * sum2);

    return (eta0 + eta1 + DELTAeta) / 1e6;  // [Pa-s]
}

CoolPropDbl TransportRoutines::viscosity_dilute_ethane(HelmholtzEOSMixtureBackend& HEOS) {
    double C[] = {
      0, -3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390, -9.6570437074e-3};
    double OMEGA_2_2 = 0, e_k = 245, Tstar;

    Tstar = HEOS.T() / e_k;
    for (int i = 1; i <= 9; i++) {
        OMEGA_2_2 += C[i] * pow(Tstar, (i - 1) / 3.0 - 1);
    }

    return 12.0085 * sqrt(Tstar) * OMEGA_2_2 / 1e6;  //[Pa-s]
}
CoolPropDbl TransportRoutines::viscosity_dilute_cyclohexane(HelmholtzEOSMixtureBackend& HEOS) {
    // From Tariq, JPCRD, 2014
    CoolPropDbl T = HEOS.T();
    CoolPropDbl S_eta = exp(-1.5093 + 364.87 / T - 39537 / pow(T, 2));  //[nm^2]
    return 0.19592 * sqrt(T) / S_eta / 1e6;                             //[Pa-s]
}
CoolPropDbl TransportRoutines::viscosity_ethane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS) {
    double r[] = {0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 1, 1};
    double s[] = {0, 0, 1, 0, 1, 1.5, 0, 2, 0, 1, 0, 1};
    double g[] = {0,           0.47177003,  -0.23950311, 0.39808301,   -0.27343335, 0.35192260,
                  -0.21101308, -0.00478579, 0.07378129,  -0.030435255, -0.30435286, 0.001215675};

    double sum1 = 0, sum2 = 0, tau = 305.33 / HEOS.T(), delta = HEOS.rhomolar() / 6870;

    for (int i = 1; i <= 9; ++i) {
        sum1 += g[i] * pow(delta, r[i]) * pow(tau, s[i]);
    }
    for (int i = 10; i <= 11; ++i) {
        sum2 += g[i] * pow(delta, r[i]) * pow(tau, s[i]);
    }
    return 15.977 * sum1 / (1 + sum2) / 1e6;
}
CoolPropDbl TransportRoutines::viscosity_Chung(HelmholtzEOSMixtureBackend& HEOS) {
    // Retrieve values from the state class
    CoolProp::ViscosityChungData& data = HEOS.components[0].transport.viscosity_Chung;

    double a0[] = {0, 6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540, -1.89992, 24.27450, 0.79716, -0.23816, 0.68629e-1};
    double a1[] = {0, 50.41190, -0.11536e-2, 254.20900, 38.09570, 7.63034, -12.53670, 3.44945, 1.11764, 0.67695e-1, 0.34793};
    double a2[] = {0, -51.68010, -0.62571e-2, -168.48100, -8.46414, -14.35440, 4.98529, -11.29130, 0.12348e-1, -0.81630, 0.59256};
    double a3[] = {0, 1189.02000, 0.37283e-1, 3898.27000, 31.41780, 31.52670, -18.15070, 69.34660, -4.11661, 4.02528, -0.72663};
    double A[11];

    if (HEOS.is_pure_or_pseudopure) {
        double Vc_cm3mol = 1 / (data.rhomolar_critical / 1e6);  // [cm^3/mol]
        double acentric = data.acentric;                        // [-]
        double M_gmol = data.molar_mass * 1000.0;               // [g/mol]
        double Tc = data.T_critical;                            // [K]
        double mu_D = data.dipole_moment_D;                     // [D]
        double kappa = 0;

        double mu_r = 131.3 * mu_D / sqrt(Vc_cm3mol * Tc);  // [-]

        for (int i = 1; i <= 10; ++i) {
            A[i] = a0[i] + a1[i] * acentric + a2[i] * pow(mu_r, 4) + a3[i] * kappa;
        }
        double F_c = 1 - 0.2756 * acentric + 0.059035 * pow(mu_r, 4) + kappa;  // [-]
        double epsilon_over_k = Tc / 1.2593;                                   // [K]

        double rho_molcm3 = HEOS.rhomolar() / 1e6;
        double T = HEOS.T();
        double Tstar = T / epsilon_over_k;
        double Omega_2_2 = 1.16145 * pow(Tstar, -0.14874) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar)
                           - 6.435e-4 * pow(Tstar, 0.14874) * sin(18.0323 * pow(Tstar, -0.76830) - 7.27371);  // [-]
        double eta0_P = 4.0785e-5 * sqrt(M_gmol * T) / (pow(Vc_cm3mol, 2.0 / 3.0) * Omega_2_2) * F_c;         // [P]

        double Y = rho_molcm3 * Vc_cm3mol / 6.0;
        double G_1 = (1.0 - 0.5 * Y) / pow(1 - Y, 3);
        double G_2 = (A[1] * (1 - exp(-A[4] * Y)) / Y + A[2] * G_1 * exp(A[5] * Y) + A[3] * G_1) / (A[1] * A[4] + A[2] + A[3]);
        double eta_k_P = eta0_P * (1 / G_2 + A[6] * Y);  // [P]

        double eta_p_P = (36.344e-6 * sqrt(M_gmol * Tc) / pow(Vc_cm3mol, 2.0 / 3.0)) * A[7] * pow(Y, 2) * G_2
                         * exp(A[8] + A[9] / Tstar + A[10] / pow(Tstar, 2));  // [P]

        return (eta_k_P + eta_p_P) / 10.0;  // [P] -> [Pa*s]
    } else {
        throw NotImplementedError("TransportRoutines::viscosity_Chung is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::conductivity_dilute_ratio_polynomials(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ConductivityDiluteRatioPolynomialsData& data = HEOS.components[0].transport.conductivity_dilute.ratio_polynomials;

        CoolPropDbl summer1 = 0, summer2 = 0, Tr = HEOS.T() / data.T_reducing;
        for (std::size_t i = 0; i < data.A.size(); ++i) {
            summer1 += data.A[i] * pow(Tr, data.n[i]);
        }
        for (std::size_t i = 0; i < data.B.size(); ++i) {
            summer2 += data.B[i] * pow(Tr, data.m[i]);
        }

        return summer1 / summer2;
    } else {
        throw NotImplementedError("TransportRoutines::conductivity_dilute_ratio_polynomials is only for pure and pseudo-pure");
    }
};

CoolPropDbl TransportRoutines::conductivity_residual_polynomial(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ConductivityResidualPolynomialData& data = HEOS.components[0].transport.conductivity_residual.polynomials;

        CoolPropDbl summer = 0, tau = data.T_reducing / HEOS.T(), delta = HEOS.keyed_output(CoolProp::iDmass) / data.rhomass_reducing;
        for (std::size_t i = 0; i < data.B.size(); ++i) {
            summer += data.B[i] * pow(tau, data.t[i]) * pow(delta, data.d[i]);
        }
        return summer;
    } else {
        throw NotImplementedError("TransportRoutines::conductivity_residual_polynomial is only for pure and pseudo-pure");
    }
};

CoolPropDbl TransportRoutines::conductivity_residual_polynomial_and_exponential(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Retrieve values from the state class
        CoolProp::ConductivityResidualPolynomialAndExponentialData& data =
          HEOS.components[0].transport.conductivity_residual.polynomial_and_exponential;

        CoolPropDbl summer = 0, tau = HEOS.tau(), delta = HEOS.delta();
        for (std::size_t i = 0; i < data.A.size(); ++i) {
            summer += data.A[i] * pow(tau, data.t[i]) * pow(delta, data.d[i]) * exp(-data.gamma[i] * pow(delta, data.l[i]));
        }
        return summer;
    } else {
        throw NotImplementedError("TransportRoutines::conductivity_residual_polynomial_and_exponential is only for pure and pseudo-pure");
    }
};

CoolPropDbl TransportRoutines::conductivity_critical_simplified_Olchowy_Sengers(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // Olchowy and Sengers cross-over term

        // Retrieve values from the state class
        CoolProp::ConductivityCriticalSimplifiedOlchowySengersData& data = HEOS.components[0].transport.conductivity_critical.Olchowy_Sengers;

        double k = data.k, R0 = data.R0, nu = data.nu, gamma = data.gamma, GAMMA = data.GAMMA, zeta0 = data.zeta0, qD = data.qD,
               Tc = HEOS.get_reducing_state().T,      // [K]
          rhoc = HEOS.get_reducing_state().rhomolar,  // [mol/m^3]
          Pcrit = HEOS.get_reducing_state().p,        // [Pa]
          Tref,                                       // [K]
          cp, cv, delta, num, zeta, mu, pi = M_PI, OMEGA_tilde, OMEGA_tilde0;

        if (ValidNumber(data.T_ref))
            Tref = data.T_ref;
        else
            Tref = 1.5 * Tc;

        delta = HEOS.delta();

        double dp_drho = HEOS.gas_constant() * HEOS.T() * (1 + 2 * delta * HEOS.dalphar_dDelta() + delta * delta * HEOS.d2alphar_dDelta2());
        double X = Pcrit / pow(rhoc, 2) * HEOS.rhomolar() / dp_drho;

        double tau_ref = Tc / Tref;
        double dp_drho_ref = HEOS.gas_constant() * Tref
                             * (1 + 2 * delta * HEOS.calc_alphar_deriv_nocache(0, 1, HEOS.mole_fractions, tau_ref, delta)
                                + delta * delta * HEOS.calc_alphar_deriv_nocache(0, 2, HEOS.mole_fractions, tau_ref, delta));
        double Xref = Pcrit / pow(rhoc, 2) * HEOS.rhomolar() / dp_drho_ref * Tref / HEOS.T();
        num = X - Xref;

        // No critical enhancement if numerator is negative, zero, or just a tiny bit positive due to roundoff
        // See also Lemmon, IJT, 2004, page 27
        if (num < DBL_EPSILON * 10)
            return 0.0;
        else
            zeta = zeta0 * pow(num / GAMMA, nu / gamma);  //[m]

        cp = HEOS.cpmolar();    //[J/mol/K]
        cv = HEOS.cvmolar();    //[J/mol/K]
        mu = HEOS.viscosity();  //[Pa-s]

        OMEGA_tilde = 2.0 / pi * ((cp - cv) / cp * atan(zeta * qD) + cv / cp * (zeta * qD));                                        //[-]
        OMEGA_tilde0 = 2.0 / pi * (1.0 - exp(-1.0 / (1.0 / (qD * zeta) + 1.0 / 3.0 * (zeta * qD) * (zeta * qD) / delta / delta)));  //[-]

        double lambda = HEOS.rhomolar() * cp * R0 * k * HEOS.T() / (6 * pi * mu * zeta) * (OMEGA_tilde - OMEGA_tilde0);  //[W/m/K]
        return lambda;                                                                                                   //[W/m/K]
    } else {
        throw NotImplementedError("TransportRoutines::conductivity_critical_simplified_Olchowy_Sengers is only for pure and pseudo-pure");
    }
};

CoolPropDbl TransportRoutines::conductivity_critical_hardcoded_R123(HelmholtzEOSMixtureBackend& HEOS) {
    double a13 = 0.486742e-2, a14 = -100, a15 = -7.08535;
    return a13 * exp(a14 * pow(HEOS.tau() - 1, 4) + a15 * pow(HEOS.delta() - 1, 2));
};

CoolPropDbl TransportRoutines::conductivity_critical_hardcoded_CO2_ScalabrinJPCRD2006(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl nc = 0.775547504e-3 * 4.81384, Tr = HEOS.T() / 304.1282, alpha, rhor = HEOS.keyed_output(iDmass) / 467.6;
    static CoolPropDbl a[] = {0.0, 3.0, 6.70697, 0.94604, 0.30, 0.30, 0.39751, 0.33791, 0.77963, 0.79857, 0.90, 0.02, 0.20};

    // Equation 6 from Scalabrin
    alpha = 1 - a[10] * acosh(1 + a[11] * pow(pow(1 - Tr, 2), a[12]));

    // Equation 5 from Scalabrin
    CoolPropDbl numer = rhor * exp(-pow(rhor, a[1]) / a[1] - pow(a[2] * (Tr - 1), 2) - pow(a[3] * (rhor - 1), 2));
    CoolPropDbl braced = (1 - 1 / Tr) + a[4] * pow(pow(rhor - 1, 2), 0.5 / a[5]);
    CoolPropDbl denom = pow(pow(pow(braced, 2), a[6]) + pow(pow(a[7] * (rhor - alpha), 2), a[8]), a[9]);
    return nc * numer / denom;
}

CoolPropDbl TransportRoutines::conductivity_dilute_hardcoded_CO2(HelmholtzEOSMixtureBackend& HEOS) {

    double e_k = 251.196, Tstar;
    double b[] = {0.4226159, 0.6280115, -0.5387661, 0.6735941, 0, 0, -0.4362677, 0.2255388};
    double c[] = {0, 2.387869e-2, 4.350794, -10.33404, 7.981590, -1.940558};

    //Vesovic Eq. 31 [no units]
    double summer = 0;
    for (int i = 1; i <= 5; i++)
        summer += c[i] * pow(HEOS.T() / 100.0, 2 - i);
    double cint_k = 1.0 + exp(-183.5 / HEOS.T()) * summer;

    //Vesovic Eq. 12 [no units]
    double r = sqrt(2.0 / 5.0 * cint_k);

    // According to REFPROP, 1+r^2 = cp-2.5R.  This is unclear to me but seems to suggest that cint/k is the difference
    // between the ideal gas specific heat and a monatomic specific heat of 5/2*R. Using the form of cint/k from Vesovic
    // does not yield exactly the correct values

    Tstar = HEOS.T() / e_k;
    //Vesovic Eq. 30 [no units]
    summer = 0;
    for (int i = 0; i <= 7; i++)
        summer += b[i] / pow(Tstar, i);
    double Gstar_lambda = summer;

    //Vesovic Eq. 29 [W/m/K]
    double lambda_0 = 0.475598e-3 * sqrt(HEOS.T()) * (1 + r * r) / Gstar_lambda;

    return lambda_0;
}

CoolPropDbl TransportRoutines::conductivity_dilute_hardcoded_ethane(HelmholtzEOSMixtureBackend& HEOS) {

    double e_k = 245.0;
    double tau = 305.33 / HEOS.T(), Tstar = HEOS.T() / e_k;
    double fint = 1.7104147 - 0.6936482 / Tstar;
    double lambda_0 = 0.276505e-3 * (HEOS.calc_viscosity_dilute() * 1e6) * (3.75 - fint * (tau * tau * HEOS.d2alpha0_dTau2() + 1.5));  //[W/m/K]

    return lambda_0;
}

CoolPropDbl TransportRoutines::conductivity_dilute_eta0_and_poly(HelmholtzEOSMixtureBackend& HEOS) {

    if (HEOS.is_pure_or_pseudopure) {
        CoolProp::ConductivityDiluteEta0AndPolyData& E = HEOS.components[0].transport.conductivity_dilute.eta0_and_poly;

        double eta0_uPas = HEOS.calc_viscosity_dilute() * 1e6;  // [uPa-s]
        double summer = E.A[0] * eta0_uPas;
        for (std::size_t i = 1; i < E.A.size(); ++i)
            summer += E.A[i] * pow(static_cast<CoolPropDbl>(HEOS.tau()), E.t[i]);
        return summer;
    } else {
        throw NotImplementedError("TransportRoutines::conductivity_dilute_eta0_and_poly is only for pure and pseudo-pure");
    }
}

CoolPropDbl TransportRoutines::conductivity_hardcoded_heavywater(HelmholtzEOSMixtureBackend& HEOS) {
    double Tbar = HEOS.T() / 643.847, rhobar = HEOS.rhomass() / 358;
    double A[] = {1.00000, 37.3223, 22.5485, 13.0465, 0, -2.60735};
    double lambda0 = A[0] + A[1] * Tbar + A[2] * POW2(Tbar) + A[3] * POW3(Tbar) + A[4] * POW4(Tbar) + A[5] * POW5(Tbar);
    double Be = -2.506, B[] = {-167.310, 483.656, -191.039, 73.0358, -7.57467};
    double DELTAlambda = B[0] * (1 - exp(Be * rhobar)) + B[1] * rhobar + B[2] * POW2(rhobar) + B[3] * POW3(rhobar) + B[4] * POW4(rhobar);
    double f_1 = exp(0.144847 * Tbar + -5.64493 * POW2(Tbar));
    double f_2 = exp(-2.80000 * POW2(rhobar - 1)) - 0.080738543 * exp(-17.9430 * POW2(rhobar - 0.125698));
    double tau = Tbar / (std::abs(Tbar - 1.1) + 1.1);
    double f_3 = 1 + exp(60 * (tau - 1) + 20);
    double f_4 = 1 + exp(100 * (tau - 1) + 15);
    double DELTAlambda_c = 35429.6 * f_1 * f_2 * (1 + POW2(f_2) * (5000.0e6 * POW4(f_1) / f_3 + 3.5 * f_2 / f_4));
    double DELTAlambda_L = -741.112 * pow(f_1, 1.2) * (1 - exp(-pow(rhobar / 2.5, 10)));
    double lambdabar = lambda0 + DELTAlambda + DELTAlambda_c + DELTAlambda_L;
    return lambdabar * 0.742128e-3;
}

CoolPropDbl TransportRoutines::conductivity_hardcoded_water(HelmholtzEOSMixtureBackend& HEOS) {

    double L[5][6] = {{1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258},
                      {2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245},
                      {2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816},
                      {-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0, 0},
                      {-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842}};

    double lambdabar_0, lambdabar_1, lambdabar_2, rhobar, Tbar, sum;
    double Tstar = 647.096, rhostar = 322, pstar = 22064000, lambdastar = 1e-3, mustar = 1e-6;
    double xi;
    int i, j;
    double R = 461.51805;  //[J/kg/K]

    Tbar = HEOS.T() / Tstar;
    rhobar = HEOS.keyed_output(CoolProp::iDmass) / rhostar;

    // Dilute gas contribution
    lambdabar_0 =
      sqrt(Tbar) / (2.443221e-3 + 1.323095e-2 / Tbar + 6.770357e-3 / pow(Tbar, 2) - 3.454586e-3 / pow(Tbar, 3) + 4.096266e-4 / pow(Tbar, 4));

    sum = 0;
    for (i = 0; i <= 4; i++) {
        for (j = 0; j <= 5; j++) {
            sum += L[i][j] * powInt(1.0 / Tbar - 1.0, i) * powInt(rhobar - 1, j);
        }
    }
    // Finite density contribution
    lambdabar_1 = exp(rhobar * sum);

    double nu = 0.630, GAMMA = 177.8514, gamma = 1.239, xi_0 = 0.13, Lambda_0 = 0.06, Tr_bar = 1.5, qd_bar = 1 / 0.4, pi = 3.141592654,
           delta = HEOS.delta();

    double drhodp = 1 / (R * HEOS.T() * (1 + 2 * rhobar * HEOS.dalphar_dDelta() + rhobar * rhobar * HEOS.d2alphar_dDelta2()));
    double drhobar_dpbar = pstar / rhostar * drhodp;
    double drhodp_Trbar = 1
                          / (R * Tr_bar * Tstar
                             * (1 + 2 * rhobar * HEOS.calc_alphar_deriv_nocache(0, 1, HEOS.mole_fractions, 1 / Tr_bar, delta)
                                + delta * delta * HEOS.calc_alphar_deriv_nocache(0, 2, HEOS.mole_fractions, 1 / Tr_bar, delta)));
    double drhobar_dpbar_Trbar = pstar / rhostar * drhodp_Trbar;
    double cp = HEOS.cpmass();  // [J/kg/K]
    double cv = HEOS.cvmass();  // [J/kg/K]
    double cpbar = cp / R;      //[-]
    double mubar = HEOS.viscosity() / mustar;
    double DELTAchibar_T = rhobar * (drhobar_dpbar - drhobar_dpbar_Trbar * Tr_bar / Tbar);
    if (DELTAchibar_T < 0)
        xi = 0;
    else
        xi = xi_0 * pow(DELTAchibar_T / Lambda_0, nu / gamma);
    double y = qd_bar * xi;

    double Z;
    double kappa = cp / cv;
    if (y < 1.2e-7)
        Z = 0;
    else
        Z = 2 / (pi * y) * (((1 - 1 / kappa) * atan(y) + y / kappa) - (1 - exp(-1 / (1 / y + y * y / 3 / rhobar / rhobar))));

    lambdabar_2 = GAMMA * rhobar * cpbar * Tbar / mubar * Z;

    return (lambdabar_0 * lambdabar_1 + lambdabar_2) * lambdastar;
}

CoolPropDbl TransportRoutines::conductivity_hardcoded_R23(HelmholtzEOSMixtureBackend& HEOS) {

    double B1 = -2.5370,     // [mW/m/K]
      B2 = 0.05366,          // [mW/m/K^2]
      C1 = 0.94215,          // [-]
      C2 = 0.14914,          // [mW/m/K^2]
      DeltaGstar = 2508.58,  //[J/mol]
      rhoL = 68.345,         // [mol/dm^3] = [mol/L]
      rhocbar = 7.5114,      // [mol/dm^3]
      DELTAlambda_max = 25,  //[mW/m/K]
      Ru = 8.31451,          // [J/mol/K]
      Tc = 299.2793,         //[K]
      T = HEOS.T();          //[K]

    double lambda_DG = B1 + B2 * T;

    double rhobar = HEOS.rhomolar() / 1000;  // [mol/L]
    double lambda_L = C2 * (rhoL * rhoL) / (rhoL - rhobar) * sqrt(T) * exp(rhobar / (rhoL - rhobar) * DeltaGstar / (Ru * T));

    double chi = rhobar - rhocbar;
    double tau = T - Tc;

    double DELTAlambda_c = 4 * DELTAlambda_max / ((exp(chi) + exp(-chi)) * (exp(tau) + exp(-tau)));

    return (pow((rhoL - rhobar) / rhoL, C1) * lambda_DG + pow(rhobar / rhoL, C1) * lambda_L + DELTAlambda_c) / 1e3;
}

CoolPropDbl TransportRoutines::conductivity_critical_hardcoded_ammonia(HelmholtzEOSMixtureBackend& HEOS) {

    /*
    From "Thermal Conductivity of Ammonia in a Large
    Temperature and Pressure Range Including the Critical Region"
    by R. Tufeu, D.Y. Ivanov, Y. Garrabos, B. Le Neindre,
    Bereicht der Bunsengesellschaft Phys. Chem. 88 (1984) 422-427
    */

    double T = HEOS.T(), Tc = 405.4, rhoc = 235, rho;
    double LAMBDA = 1.2, nu = 0.63, gamma = 1.24, DELTA = 0.50, t, zeta_0_plus = 1.34e-10, a_zeta = 1, GAMMA_0_plus = 0.423e-8;
    double pi = 3.141592654, a_chi, k_B = 1.3806504e-23, X_T, DELTA_lambda, dPdT, eta_B, DELTA_lambda_id, DELTA_lambda_i;

    rho = HEOS.keyed_output(CoolProp::iDmass);
    t = std::abs((T - Tc) / Tc);
    a_chi = a_zeta / 0.7;
    eta_B = (2.60 + 1.6 * t) * 1e-5;
    dPdT = (2.18 - 0.12 / exp(17.8 * t)) * 1e5;  // [Pa-K]
    X_T = 0.61 * rhoc + 16.5 * log(t);
    // Along the critical isochore (only a function of temperature) (Eq. 9)
    DELTA_lambda_i = LAMBDA * (k_B * T * T) / (6 * pi * eta_B * (zeta_0_plus * pow(t, -nu) * (1 + a_zeta * pow(t, DELTA)))) * dPdT * dPdT
                     * GAMMA_0_plus * pow(t, -gamma) * (1 + a_chi * pow(t, DELTA));
    DELTA_lambda_id = DELTA_lambda_i * exp(-36 * t * t);
    if (rho < 0.6 * rhoc) {
        DELTA_lambda = DELTA_lambda_id * (X_T * X_T) / (X_T * X_T + powInt(0.6 * rhoc - 0.96 * rhoc, 2)) * powInt(rho, 2) / powInt(0.6 * rhoc, 2);
    } else {
        DELTA_lambda = DELTA_lambda_id * (X_T * X_T) / (X_T * X_T + powInt(rho - 0.96 * rhoc, 2));
    }

    return DELTA_lambda;
}

CoolPropDbl TransportRoutines::conductivity_hardcoded_helium(HelmholtzEOSMixtureBackend& HEOS) {
    /*
    What an incredibly annoying formulation!  Implied coefficients?? Not cool.
    */
    double rhoc = 68.0, lambda_e, lambda_c, T = HEOS.T(), rho = HEOS.rhomass();
    double summer = 3.739232544 / T - 2.620316969e1 / T / T + 5.982252246e1 / T / T / T - 4.926397634e1 / T / T / T / T;
    double lambda_0 = 2.7870034e-3 * pow(T, 7.034007057e-1) * exp(summer);
    double c[] = {1.862970530e-4,  -7.275964435e-7, -1.427549651e-4, 3.290833592e-5, -5.213335363e-8, 4.492659933e-8,
                  -5.924416513e-9, 7.087321137e-6,  -6.013335678e-6, 8.067145814e-7, 3.995125013e-7};
    // Equation 17
    lambda_e = (c[0] + c[1] * T + c[2] * pow(T, 1 / 3.0) + c[3] * pow(T, 2.0 / 3.0)) * rho
               + (c[4] + c[5] * pow(T, 1.0 / 3.0) + c[6] * pow(T, 2.0 / 3.0)) * rho * rho * rho
               + (c[7] + c[8] * pow(T, 1.0 / 3.0) + c[9] * pow(T, 2.0 / 3.0) + c[10] / T) * rho * rho * log(rho / rhoc);

    // Critical component
    lambda_c = 0.0;

    if (3.5 < T && T < 12) {
        double x0 = 0.392, E1 = 2.8461, E2 = 0.27156, beta = 0.3554, gamma = 1.1743, delta = 4.304, rhoc_crit = 69.158, Tc = 5.18992, pc = 2.2746e5;

        double DeltaT = std::abs(1 - T / Tc), DeltaRho = std::abs(1 - rho / rhoc_crit);
        double eta = HEOS.viscosity();  // [Pa-s]
        double K_T = HEOS.isothermal_compressibility(), K_Tprime, K_Tbar;
        double dpdT = HEOS.first_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar);

        double W = pow(DeltaT / 0.2, 2) + pow(DeltaRho / 0.25, 2);

        if (W > 1) {
            K_Tbar = K_T;
        } else {
            double x = pow(DeltaT / DeltaRho, 1 / beta);
            double h = E1 * (1 + x / x0) * pow(1 + E2 * pow(1 + x / x0, 2 / beta), (gamma - 1) / (2 * beta));

            /**
            dh/dx derived using sympy:

                E1,x,x0,E2,beta,gamma = symbols('E1,x,x0,E2,beta,gamma')
                h = E1*(1 + x/x0)*pow(1 + E2*pow(1 + x/x0, 2/beta), (gamma-1)/(2*beta))
                ccode(simplify(diff(h,x)))
            */
            double dhdx =
              E1
              * (E2 * pow((x + x0) / x0, 2 / beta) * (gamma - 1) * pow(E2 * pow((x + x0) / x0, 2 / beta) + 1, (1.0 / 2.0) * (gamma - 1) / beta)
                 + pow(beta, 2) * pow(E2 * pow((x + x0) / x0, 2 / beta) + 1, (1.0 / 2.0) * (2 * beta + gamma - 1) / beta))
              / (pow(beta, 2) * x0 * (E2 * pow((x + x0) / x0, 2 / beta) + 1));
            // Right-hand-side of Equation 9
            double RHS = pow(DeltaRho, delta - 1) * (delta * h - x / beta * dhdx);
            K_Tprime = 1 / (RHS * pow(rho / rhoc_crit, 2) * pc);
            K_Tbar = W * K_T + (1 - W) * K_Tprime;
        }

        // 3.4685233d-17 and 3.726229668d0 are "magical" coefficients that are present in the REFPROP source to yield the right values. Not clear why these values are needed.
        // Also, the form of the critical term in REFPROP does not agree with Hands paper. EL and MH from NIST are not sure where these coefficients come from.
        lambda_c =
          3.4685233e-17 * 3.726229668 * sqrt(K_Tbar) * pow(T, 2) / rho / eta * pow(dpdT, 2) * exp(-18.66 * pow(DeltaT, 2) - 4.25 * pow(DeltaRho, 4));
    }
    return lambda_0 + lambda_e + lambda_c;
}
CoolPropDbl TransportRoutines::conductivity_hardcoded_methane(HelmholtzEOSMixtureBackend& HEOS) {

    double delta = HEOS.rhomolar() / 10139.0, tau = 190.55 / HEOS.T();
    double lambda_dilute, lambda_residual, lambda_critical;

    // Viscosity formulation from Friend, JPCRD, 1989
    // Dilute
    double C[] = {
      0, -3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390, -9.6570437074e-3};
    double OMEGA22_summer = 0;
    double t = HEOS.T() / 174.0;
    for (int i = 1; i <= 9; ++i) {
        OMEGA22_summer += C[i] * pow(t, (i - 1.0) / 3.0 - 1.0);
    }
    double eta_dilute = 10.50 * sqrt(t) * OMEGA22_summer;
    double re[] = {0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 1, 1};
    double se[] = {0, 0, 1, 0, 1, 1.5, 0, 2, 0, 1, 0, 1};
    double ge[] = {0,           0.41250137, -0.14390912, 0.10366993,  0.40287464,  -0.24903524,
                   -0.12953131, 0.06575776, 0.02566628,  -0.03716526, -0.38798341, 0.03533815};
    double summer1 = 0, summer2 = 0;
    for (int i = 1; i <= 9; ++i) {
        summer1 += ge[i] * pow(delta, re[i]) * pow(tau, se[i]);
    }
    for (int i = 10; i <= 11; ++i) {
        summer2 += ge[i] * pow(delta, re[i]) * pow(tau, se[i]);
    }
    double eta_residual = 12.149 * summer1 / (1 + summer2);
    double eta = eta_residual + eta_dilute;

    // Dilute
    double f_int = 1.458850 - 0.4377162 / t;
    lambda_dilute = 0.51828 * eta_dilute * (3.75 - f_int * (POW2(HEOS.tau()) * HEOS.d2alpha0_dTau2() + 1.5));  // [mW/m/K]
    // Residual
    double rl[] = {0, 1, 3, 4, 4, 5, 5, 2};
    double sl[] = {0, 0, 0, 0, 1, 0, 1, 0};
    double jl[] = {0, 2.4149207, 0.55166331, -0.52837734, 0.073809553, 0.24465507, -0.047613626, 1.5554612};
    double summer = 0;
    for (int i = 1; i <= 6; ++i) {
        summer += jl[i] * pow(delta, rl[i]) * pow(tau, sl[i]);
    }
    double delta_sigma_star = 1.0;  // Looks like a typo in Friend - should be 1 instead of 11
    if (HEOS.T() < HEOS.T_critical() && HEOS.rhomolar() < HEOS.rhomolar_critical()) {
        delta_sigma_star = HEOS.saturation_ancillary(iDmolar, 1, iT, HEOS.T()) / HEOS.keyed_output(CoolProp::irhomolar_critical);
    }
    lambda_residual = 6.29638 * (summer + jl[7] * POW2(delta) / delta_sigma_star);  // [mW/m/K]
    // Critical region
    double Tstar = 1 - 1 / tau;
    double rhostar = 1 - delta;
    double F_T = 2.646, F_rho = 2.678, F_A = -0.637;
    double F = exp(-F_T * sqrt(std::abs(Tstar)) - F_rho * POW2(rhostar) - F_A * rhostar);
    double CHI_T_star;
    if (std::abs(Tstar) < 0.03) {
        if (std::abs(rhostar) < 1e-16) {
            // Equation 26
            const double LAMBDA = 0.0801, gamma = 1.190;
            CHI_T_star = LAMBDA * pow(std::abs(Tstar), -gamma);
        } else if (std::abs(rhostar) < 0.03) {
            // Equation 23
            const double beta = 0.355, W = -1.401, S = -6.098, E = 0.287, a = 3.352, b = 0.732, R = 0.535, Q = 0.1133;
            double OMEGA = W * Tstar * pow(std::abs(rhostar), -1 / beta);
            double theta = 1;
            if (Tstar < -pow(std::abs(rhostar), -1 / beta) / S) {
                theta = 1 + E * pow(1 + S * Tstar * pow(std::abs(rhostar), -1 / beta), 2 * beta);
            }
            CHI_T_star = Q * pow(std::abs(rhostar), -a) * pow(theta, b) / (theta + OMEGA * (theta + R));
        } else {
            // Equation 19a
            CHI_T_star = 0.28631 * delta * tau / (1 + 2 * delta * HEOS.dalphar_dDelta() + POW2(delta) * HEOS.d2alphar_dDelta2());
        }
    } else {
        // Equation 19a
        CHI_T_star = 0.28631 * delta * tau / (1 + 2 * delta * HEOS.dalphar_dDelta() + POW2(delta) * HEOS.d2alphar_dDelta2());
    }

    lambda_critical = 91.855 / (eta * POW2(tau)) * POW2(1 + delta * HEOS.dalphar_dDelta() - delta * tau * HEOS.d2alphar_dDelta_dTau())
                      * pow(CHI_T_star, 0.4681) * F;  //[mW/m/K]
    return (lambda_dilute + lambda_residual + lambda_critical) * 0.001;
}

void TransportRoutines::conformal_state_solver(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference, CoolPropDbl& T0,
                                               CoolPropDbl& rhomolar0) {
    int iter = 0;
    double resid = 9e30, resid_old = 9e30;
    CoolPropDbl alphar = HEOS.alphar();
    CoolPropDbl Z = HEOS.keyed_output(iZ);

    Eigen::Vector2d r;
    Eigen::Matrix2d J;
    HEOS_Reference.specify_phase(iphase_gas);  // Something homogeneous, not checked
    // Update the reference fluid with the conformal state
    HEOS_Reference.update_DmolarT_direct(rhomolar0, T0);
    do {
        CoolPropDbl dtau_dT = -HEOS_Reference.T_critical() / (T0 * T0);
        CoolPropDbl ddelta_drho = 1 / HEOS_Reference.rhomolar_critical();
        // Independent variables are T0 and rhomolar0, residuals are matching alphar and Z
        r(0) = HEOS_Reference.alphar() - alphar;
        r(1) = HEOS_Reference.keyed_output(iZ) - Z;
        J(0, 0) = HEOS_Reference.dalphar_dTau() * dtau_dT;
        J(0, 1) = HEOS_Reference.dalphar_dDelta() * ddelta_drho;
        // Z = 1+delta*dalphar_ddelta(tau,delta)
        // dZ_dT
        J(1, 0) = HEOS_Reference.delta() * HEOS_Reference.d2alphar_dDelta_dTau() * dtau_dT;
        // dZ_drho
        J(1, 1) = (HEOS_Reference.delta() * HEOS_Reference.d2alphar_dDelta2() + HEOS_Reference.dalphar_dDelta()) * ddelta_drho;
        // Step in v obtained from Jv = -r
        Eigen::Vector2d v = J.colPivHouseholderQr().solve(-r);
        bool good_solution = false;
        double T0_init = HEOS_Reference.T(), rhomolar0_init = HEOS_Reference.rhomolar();
        // Calculate the old residual after the last step
        resid_old = sqrt(POW2(r(0)) + POW2(r(1)));
        for (double frac = 1.0; frac > 0.001; frac /= 2) {
            try {
                // Calculate new values
                double T_new = T0_init + frac * v(0);
                double rhomolar_new = rhomolar0_init + frac * v(1);
                // Update state with step
                HEOS_Reference.update_DmolarT_direct(rhomolar_new, T_new);
                resid = sqrt(POW2(HEOS_Reference.alphar() - alphar) + POW2(HEOS_Reference.keyed_output(iZ) - Z));
                if (resid > resid_old) {
                    continue;
                }
                good_solution = true;
                T0 = T_new;
                rhomolar0 = rhomolar_new;
                break;
            } catch (...) {
                continue;
            }
        }
        if (!good_solution) {
            throw ValueError(format("Not able to get a solution"));
        }
        iter++;
        if (iter > 50) {
            throw ValueError(format("conformal_state_solver took too many iterations; residual is %g; prior was %g", resid, resid_old));
        }
    } while (std::abs(resid) > 1e-9);
}

CoolPropDbl TransportRoutines::viscosity_ECS(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference) {
    // Collect some parameters
    CoolPropDbl M = HEOS.molar_mass(), M0 = HEOS_Reference.molar_mass(), Tc = HEOS.T_critical(), Tc0 = HEOS_Reference.T_critical(),
                rhocmolar = HEOS.rhomolar_critical(), rhocmolar0 = HEOS_Reference.rhomolar_critical();

    // Get a reference to the ECS data
    CoolProp::ViscosityECSVariables& ECS = HEOS.components[0].transport.viscosity_ecs;

    // The correction polynomial psi_eta
    double psi = 0;
    for (std::size_t i = 0; i < ECS.psi_a.size(); i++) {
        psi += ECS.psi_a[i] * pow(HEOS.rhomolar() / ECS.psi_rhomolar_reducing, ECS.psi_t[i]);
    }

    // The dilute gas portion for the fluid of interest [Pa-s]
    CoolPropDbl eta_dilute = viscosity_dilute_kinetic_theory(HEOS);

    // ************************************
    // Start with a guess for theta and phi
    // ************************************
    CoolPropDbl theta = 1;
    CoolPropDbl phi = 1;

    // The equivalent substance reducing ratios
    CoolPropDbl f = Tc / Tc0 * theta;
    CoolPropDbl h = rhocmolar0 / rhocmolar * phi;  // Must be the ratio of MOLAR densities!!

    // To be solved for
    CoolPropDbl T0 = HEOS.T() / f;
    CoolPropDbl rhomolar0 = HEOS.rhomolar() * h;

    // **************************
    // Solver for conformal state
    // **************************

    //
    HEOS_Reference.specify_phase(iphase_gas);  // something homogeneous

    conformal_state_solver(HEOS, HEOS_Reference, T0, rhomolar0);

    // Update the reference fluid with the updated conformal state
    HEOS_Reference.update_DmolarT_direct(rhomolar0 * psi, T0);

    // Recalculate ESRR
    f = HEOS.T() / T0;
    h = rhomolar0 / HEOS.rhomolar();  // Must be the ratio of MOLAR densities!!

    // **********************
    // Remaining calculations
    // **********************

    // The reference fluid's contribution to the viscosity [Pa-s]
    CoolPropDbl eta_resid = HEOS_Reference.calc_viscosity_background();

    // The F factor
    CoolPropDbl F_eta = sqrt(f) * pow(h, -static_cast<CoolPropDbl>(2.0L / 3.0L)) * sqrt(M / M0);

    // The total viscosity considering the contributions of the fluid of interest and the reference fluid [Pa-s]
    CoolPropDbl eta = eta_dilute + eta_resid * F_eta;

    return eta;
}

CoolPropDbl TransportRoutines::viscosity_rhosr(HelmholtzEOSMixtureBackend& HEOS) {

    // Get a reference to the data
    const CoolProp::ViscosityRhoSrVariables& data = HEOS.components[0].transport.viscosity_rhosr;

    // The dilute gas portion for the fluid of interest [Pa-s]
    CoolPropDbl eta_dilute = viscosity_dilute_kinetic_theory(HEOS);

    // Calculate x
    double x = HEOS.rhomolar() * HEOS.gas_constant() * (HEOS.tau() * HEOS.dalphar_dTau() - HEOS.alphar()) / data.rhosr_critical;

    // Crossover variable
    double psi_liq = 1 / (1 + exp(-100.0 * (x - 2)));

    // Evaluated using Horner's method
    const std::vector<double>&cL = data.c_liq, cV = data.c_vap;
    double f_liq = cL[0] + x * (cL[1] + x * (cL[2] + x * (cL[3])));
    double f_vap = cV[0] + x * (cV[1] + x * (cV[2] + x * (cV[3])));

    // Evaluate the reference fluid
    double etastar_ref = exp(psi_liq * f_liq + (1 - psi_liq) * f_vap);

    // Get the non-dimensionalized viscosity
    double etastar_fluid = 1 + data.C * (etastar_ref - 1);

    return etastar_fluid * eta_dilute;
}

CoolPropDbl TransportRoutines::conductivity_ECS(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference) {
    // Collect some parameters
    CoolPropDbl M = HEOS.molar_mass(), M_kmol = M * 1000, M0 = HEOS_Reference.molar_mass(), Tc = HEOS.T_critical(), Tc0 = HEOS_Reference.T_critical(),
                rhocmolar = HEOS.rhomolar_critical(), rhocmolar0 = HEOS_Reference.rhomolar_critical(), R_u = HEOS.gas_constant(),
                R = HEOS.gas_constant() / HEOS.molar_mass(),  //[J/kg/K]
      R_kJkgK = R_u / M_kmol;

    // Get a reference to the ECS data
    CoolProp::ConductivityECSVariables& ECS = HEOS.components[0].transport.conductivity_ecs;

    // The correction polynomial psi_eta in rho/rho_red
    double psi = 0;
    for (std::size_t i = 0; i < ECS.psi_a.size(); ++i) {
        psi += ECS.psi_a[i] * pow(HEOS.rhomolar() / ECS.psi_rhomolar_reducing, ECS.psi_t[i]);
    }

    // The correction polynomial f_int in T/T_red
    double fint = 0;
    for (std::size_t i = 0; i < ECS.f_int_a.size(); ++i) {
        fint += ECS.f_int_a[i] * pow(HEOS.T() / ECS.f_int_T_reducing, ECS.f_int_t[i]);
    }

    // The dilute gas density for the fluid of interest [uPa-s]
    CoolPropDbl eta_dilute = viscosity_dilute_kinetic_theory(HEOS) * 1e6;

    // The mass specific ideal gas constant-pressure specific heat [J/kg/K]
    CoolPropDbl cp0 = HEOS.calc_cpmolar_idealgas() / HEOS.molar_mass();

    // The internal contribution to the thermal conductivity [W/m/K]
    CoolPropDbl lambda_int = fint * eta_dilute * (cp0 - 2.5 * R) / 1e3;

    // The dilute gas contribution to the thermal conductivity [W/m/K]
    CoolPropDbl lambda_dilute = 15.0e-3 / 4.0 * R_kJkgK * eta_dilute;

    // ************************************
    // Start with a guess for theta and phi
    // ************************************

    CoolPropDbl theta = 1;
    CoolPropDbl phi = 1;

    // The equivalent substance reducing ratios
    CoolPropDbl f = Tc / Tc0 * theta;
    CoolPropDbl h = rhocmolar0 / rhocmolar * phi;  // Must be the ratio of MOLAR densities!!

    // Initial values for the conformal state
    CoolPropDbl T0 = HEOS.T() / f;
    CoolPropDbl rhomolar0 = HEOS.rhomolar() * h;

    // **************************
    // Solver for conformal state
    // **************************

    try {
        conformal_state_solver(HEOS, HEOS_Reference, T0, rhomolar0);
    } catch (std::exception& e) {
        throw ValueError(format("Conformal state solver failed; error was %s", e.what()));
    }

    // Update the reference fluid with the conformal state
    HEOS_Reference.update(DmolarT_INPUTS, rhomolar0 * psi, T0);

    // Recalculate ESRR
    f = HEOS.T() / T0;
    h = rhomolar0 / HEOS.rhomolar();  // Must be the ratio of MOLAR densities!!

    // The reference fluid's contribution to the conductivity [W/m/K]
    CoolPropDbl lambda_resid = HEOS_Reference.calc_conductivity_background();

    // The F factor
    CoolPropDbl F_lambda = sqrt(f) * pow(h, static_cast<CoolPropDbl>(-2.0 / 3.0)) * sqrt(M0 / M);

    // The critical contribution from the fluid of interest [W/m/K]
    CoolPropDbl lambda_critical = conductivity_critical_simplified_Olchowy_Sengers(HEOS);

    // The total thermal conductivity considering the contributions of the fluid of interest and the reference fluid [Pa-s]
    CoolPropDbl lambda = lambda_int + lambda_dilute + lambda_resid * F_lambda + lambda_critical;

    return lambda;
}

}; /* namespace CoolProp */
