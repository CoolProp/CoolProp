"""Reference implementation of the RES thermal-conductivity model of Li et al. (2024).

    Li, Duan, Yang,
    "Linking Thermal Conductivity to Equations of State Using the Residual Entropy Scaling Theory",
    Ind. Eng. Chem. Res. 63, 18160-18175 (2024).  https://doi.org/10.1021/acs.iecr.4c02946

Redistributed and modified under the Creative Commons Attribution 4.0 International licence
(CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), which the article and its supporting
information are published under.  The accompanying data files in this directory
(Dilute_gas_TC.txt, Fluid_Constants.txt, RES_Parameter.txt, Samples_*.txt, Table_S*_SI.txt) are
the authors' own and are reproduced VERBATIM.

CC BY 4.0 requires that changes be indicated.  This file is `code_SI/code_SI.py` from the
supporting information with the following modifications, each marked `# ADAPTED:` at the point of
change:

  1. The equation of state is a parameter (`eos=`) rather than the hardcoded "REFPROP", which the
     original names both as `TheEOS` and in 23 separate 'REFPROP::' string literals.
  2. The dilute-gas term is selectable (`dilute_source=`) rather than always the backend's native
     lambda0 for mixtures.
  3. The viscosity consumed by the Olchowy-Sengers enhancement is selectable
     (`enhancement_viscosity=`).  The original takes it from the backend's own transport model;
     CoolProp uses the RES viscosity instead, which is the only option on backends that have no
     transport model at all.  "res" routes the call into the Martinek 2025 module next door.
  4. Data files are located relative to THIS FILE rather than the process working directory.
  5. `delim_whitespace=True` -> `sep=r'\\s+'`; the former was removed in pandas 2.2.
  6. `TC_RES` additionally returns the dilute-gas term it used, so a caller can weight a
     comparison by how much of the conductivity is actually residual.
  7. The `if __name__ == '__main__'` sample driver is split into `write_sample_table()` and
     `write_binary_sample_table()`, which do the same thing on demand.
  8. That driver formats the conductivity through float() first.  TC_RES returns a 1-element
     numpy array for a pure fluid, and `"%f" % array` raises TypeError from numpy 2.0 onwards
     where it used to work; the published Table_S5_SI.txt predates that change.  Presentation
     only -- the value written is the same one.
  9. One added comment in TC_RES() recording the coefficient order.  No code change.

NOT changed, deliberately: the Olchowy functions use kB = 1.38064852e-23 where the rest of this
file uses 1.380649e-23 (CoolProp uses the latter throughout).  That is a 1.6e-7 relative
difference on the enhancement term alone, and rewriting the authors' numerics would make this
less of a reference.

Every default reproduces the published behaviour, so running with defaults regenerates both
Table_S5_SI.txt and Table_S6_SI.txt byte-for-byte -- which is how dev/RES_reference_run.py
--self-test checks that these modifications did not change the model.
"""

import os
import sys

import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd

# ADAPTED (3): used by _enhancement_viscosity() below when enhancement_viscosity="res".
from ..martinek2025_viscosity.viscosity_RES import viscosity_RES as _martinek_viscosity_RES

# ADAPTED (4): the originals read these by bare filename from the working directory.
_HERE = os.path.dirname(os.path.realpath(__file__))


def geometric_mean(fraction, aa):
    if abs(np.sum(fraction) - 1) > 1e-6:
        print('sum(fraction)~=1 in geometric_mean')
        sys.exit()

    ncomp = len(aa)
    aiaj = np.zeros((ncomp, ncomp))
    aa_mix = 0
    for ii in range(ncomp):
        for jj in range(ncomp):
            aiaj[ii, jj] = np.sqrt(aa[ii] * aa[jj])
            aa_mix += fraction[ii] * fraction[jj] * aiaj[ii, jj]
    return aa_mix


def get_paramters(AllMaterial):

    # ADAPTED (4, 5)
    Fluid_Constants = pd.read_csv(os.path.join(_HERE, "Fluid_Constants.txt"), sep=r'\s+', header=0, skipinitialspace=True)
    RES_Parameter = pd.read_csv(os.path.join(_HERE, "RES_Parameter.txt"), sep=r'\s+', header=0, skipinitialspace=True)
    Dilute_Parameter = pd.read_csv(os.path.join(_HERE, "Dilute_gas_TC.txt"), skipinitialspace=True, header=1, delimiter=';')
    Len_FC = len(Fluid_Constants['Material'])

    if isinstance(AllMaterial, str):
        nfluid = 1
    else:
        nfluid = len(AllMaterial)

    parameter_N = np.zeros((nfluid, 4))
    epsilon_kB_K = np.zeros(nfluid)
    sigma_nm = np.zeros(nfluid)
    Tc_K = np.zeros(nfluid)
    cp0_k0 = np.zeros(nfluid)
    cp0_k1 = np.zeros(nfluid)
    R_D = np.zeros(nfluid)
    gamma = np.zeros(nfluid)
    xi0 = np.zeros(nfluid)
    Gamma = np.zeros(nfluid)
    qDinv = np.zeros(nfluid)
    Tref = np.zeros(nfluid)
    group_n = np.zeros(nfluid)
    xita = np.zeros(nfluid)
    Tc0_N = np.zeros((nfluid, 5))
    for ifluid in range(nfluid):
        if isinstance(AllMaterial, str):
            material = AllMaterial
        else:
            material = AllMaterial[ifluid]
        for ii in range(Len_FC):
            if material.upper() == Fluid_Constants['Material'][ii].upper():
                if Fluid_Constants['Material'][ii].upper() != Fluid_Constants['Material'][ii].upper():
                    print(material + ' are not the same in Fluid_Constants.txt and RES_Parameter.txt')
                    sys.exit()
                epsilon_kB_K[ifluid] = Fluid_Constants['epsilon_kB_K'][ii]
                sigma_nm[ifluid] = Fluid_Constants['sigma_nm'][ii]
                Tc_K[ifluid] = Fluid_Constants['Tc_K'][ii]
                cp0_k0[ifluid] = Fluid_Constants['k0_cp0_JmolK'][ii]
                cp0_k1[ifluid] = Fluid_Constants['k1_cp0_JmolK'][ii]
                R_D[ifluid] = 1.02
                gamma[ifluid] = Fluid_Constants['gamma_uni'][ii]
                xi0[ifluid] = Fluid_Constants['xi0'][ii]
                Gamma[ifluid] = Fluid_Constants['Gamma'][ii]
                qDinv[ifluid] = Fluid_Constants['qDinv'][ii]
                Tref[ifluid] = Fluid_Constants['Tref'][ii]
                group_n[ifluid] = RES_Parameter['GroupN'][ii]
                xita[ifluid] = RES_Parameter['xita'][ii]
                Tc0_N[ifluid, :] = [
                    Dilute_Parameter['n0'][ii],
                    Dilute_Parameter['n1'][ii],
                    Dilute_Parameter['n2'][ii],
                    Dilute_Parameter['n3'][ii],
                    Dilute_Parameter['n4'][ii],
                ]
                if RES_Parameter['ind_fit'][ii] == 1:
                    parameter_N[ifluid, :] = [
                        RES_Parameter['n1_ind'][ii],
                        RES_Parameter['n2_ind'][ii],
                        RES_Parameter['n3_ind'][ii],
                        RES_Parameter['n4_ind'][ii],
                    ]
                    xita[ifluid] = 1
                else:
                    parameter_N[ifluid, :] = [
                        RES_Parameter['n1_glb'][ii],
                        RES_Parameter['n2_glb'][ii],
                        RES_Parameter['n3_glb'][ii],
                        RES_Parameter['n4_glb'][ii],
                    ]
                break
            if ii == Len_FC - 1:
                print(material + ' not found in data/Fluid_Constants.txt')
                sys.exit()
    return (
        epsilon_kB_K, sigma_nm, parameter_N, xita, group_n, Tc0_N, Tc_K,
        cp0_k0, cp0_k1, R_D, gamma, xi0, Gamma, qDinv, Tref,
    )


def lambda_dilute_gas_multi_wilke_fast(lambda0, mass_kg, MoleFrac):
    ncomp = len(MoleFrac)
    fai = np.zeros((ncomp, ncomp))
    summ = np.zeros(ncomp)
    lambda_mix = 0
    for icomp in range(ncomp):
        for jcomp in range(ncomp):
            fai[icomp][jcomp] = (1 + np.sqrt(lambda0[icomp] / lambda0[jcomp]) * (mass_kg[jcomp] / mass_kg[icomp]) ** (1 / 4)) ** 2 / np.sqrt(
                8 * (1 + mass_kg[icomp] / mass_kg[jcomp])
            )
            summ[icomp] += MoleFrac[jcomp] * fai[icomp][jcomp]
        lambda_mix += MoleFrac[icomp] * lambda0[icomp] / summ[icomp]
    return lambda_mix


def _enhancement_viscosity(mode, eos, spec, MoleFrac, T, rho, p_Pa):
    """ADAPTED (3): which viscosity the enhancement term consumes.

    "native" is the paper: the backend's own transport model at (T, rho).
    "res"    is what CoolProp does: the Martinek 2025 RES viscosity at the SAME state, reached
             through pressure because that is the entry point that module offers.  `spec` is the
             fluid spelling that module expects -- a bare name for a pure fluid, a list for a
             mixture -- and `p_Pa` the pressure of the state being evaluated.
    """
    if mode == "native":
        return CP.PropsSI('V', 'T', T, 'Dmass', rho, eos + '::' + spec if isinstance(spec, str) else spec)
    val = _martinek_viscosity_RES(spec, MoleFrac, p_Pa, T, eos=eos, dilute_source="polynomial")[0]
    return float(np.ravel(val)[0])


def Olchowy_critical_enhancement(Material, T, rho, R_D, gamma, xi0, Gamma, qDinv, Tref,
                                 eos="REFPROP", enhancement_viscosity="native", MoleFrac=None, p_Pa=None):
    # K  kg/m3 none none  m    none    m      K
    if rho == 0.0:
        return 0.0

    Tc = CP.PropsSI('Tcrit', 'T', 0, 'p', 0, eos + '::' + Material)
    rhoc = CP.PropsSI('rhocrit', 'T', 0, 'p', 0, eos + '::' + Material)
    if rho / rhoc >= 2 or T / Tc > 1.4:
        return 0.0

    kB = 1.38064852e-23
    pi = 3.1415926535

    # R_D = 1.02
    nu = 0.63
    # gamma = 1.239

    # ADAPTED (3): was CP.PropsSI('V', 'T', T, 'Dmass', rho, 'REFPROP::' + Material)
    eta = _enhancement_viscosity(enhancement_viscosity, eos, Material, MoleFrac, T, rho, p_Pa)  # Pa s
    molemass = CP.PropsSI('molemass', 'T', 0, 'p', 0, eos + '::' + Material)  # kg/mol
    cp = CP.PropsSI('Cpmolar', 'T', T, 'Dmass', rho, eos + '::' + Material)  # J/mol/K
    cv = CP.PropsSI('Cvmolar', 'T', T, 'Dmass', rho, eos + '::' + Material)  # J/mol/K

    pc = CP.PropsSI('Pcrit', 'T', 0, 'p', 0, eos + '::' + Material)  # Pa

    dpdrho_T = CP.PropsSI('d(P)/d(Dmass)|T', 'T', T, 'Dmass', rho, eos + '::' + Material)
    dpdrho_Tref = CP.PropsSI('d(P)/d(Dmass)|T', 'T', Tref, 'Dmass', rho, eos + '::' + Material)

    kappa = cp / cv
    delta = rho / rhoc
    arg = 1 / dpdrho_T - Tref / T / dpdrho_Tref

    if arg < 0:
        return 0.0
    xi = xi0 * ((pc * rho) / (Gamma * rhoc**2)) ** (nu / gamma) * (arg) ** (nu / gamma)
    y = xi / qDinv
    Omega = 2 / pi * ((1 - 1 / kappa) * np.arctan(y) + y / kappa)
    Omega0 = 2 / pi * (1 - np.exp(-1 / (1 / y + (y / delta) ** 2 / 3)))
    val = kB * R_D * rho * cp * T / (eta * xi) / (6 * pi) * (Omega - Omega0) / molemass

    return val


def Olchowy_critical_enhancement_mix(fluid_ref, MoleFrac, T, rho, R_D, gamma, xi0, Gamma, qDinv, Tref,
                                     eos="REFPROP", enhancement_viscosity="native", components=None, p_Pa=None):
    # K  kg/m3 none none  m    none    m      K
    if rho == 0.0:
        return 0.0
    try:

        Ncomp = len(MoleFrac)
        RefMixName = fluid_ref

        Tc = CP.PropsSI('Tcrit', 'T', 0, 'p', 0, eos + '::' + RefMixName)
        rhoc = CP.PropsSI('rhocrit', 'T', 0, 'p', 0, eos + '::' + RefMixName)
        if rho / rhoc >= 2 or T / Tc > 1.4:
            return 0.0

        kB = 1.38064852e-23
        pi = 3.1415926535

        # R_D = 1.02
        nu = 0.63
        # gamma = 1.239

        # ADAPTED (3): was CP.PropsSI('V', 'T', T, 'Dmass', rho, 'REFPROP::' + RefMixName)
        eta = _enhancement_viscosity(
            enhancement_viscosity, eos, RefMixName if enhancement_viscosity == "native" else components,
            MoleFrac, T, rho, p_Pa,
        )  # Pa s
        molemass = CP.PropsSI('molemass', 'T', 0, 'p', 0, eos + '::' + RefMixName)  # kg/mol
        cp = CP.PropsSI('Cpmolar', 'T', T, 'Dmass', rho, eos + '::' + RefMixName)  # J/mol/K
        cv = CP.PropsSI('Cvmolar', 'T', T, 'Dmass', rho, eos + '::' + RefMixName)  # J/mol/K
        pc = CP.PropsSI('Pcrit', 'T', 0, 'p', 0, eos + '::' + RefMixName)  # Pa

        Gamma_mix, R_D_mix, gamma_mix, xi0_mix, qDinv_mix = [0, 0, 0, 0, 0]
        for icomp in range(Ncomp):
            R_D_mix += MoleFrac[icomp] * R_D[icomp]
            gamma_mix += MoleFrac[icomp] * gamma[icomp]
            xi0_mix += MoleFrac[icomp] * xi0[icomp]
            Gamma_mix += MoleFrac[icomp] * Gamma[icomp]
            qDinv_mix += MoleFrac[icomp] * qDinv[icomp]

        Tref_Mix = 1.5 * Tc

        dpdrho_T = CP.PropsSI('d(P)/d(Dmass)|T', 'T', T, 'Dmass', rho, eos + '::' + RefMixName)
        dpdrho_Tref = CP.PropsSI('d(P)/d(Dmass)|T', 'T', Tref_Mix, 'Dmass', rho, eos + '::' + RefMixName)

        kappa = cp / cv
        delta = rho / rhoc
        arg = 1 / dpdrho_T - Tref_Mix / T / dpdrho_Tref

        if arg < 0:
            return 0.0
        xi = xi0_mix * ((pc * rho) / (Gamma_mix * rhoc**2)) ** (nu / gamma_mix) * (arg) ** (nu / gamma_mix)
        y = xi / qDinv_mix
        Omega = 2 / pi * ((1 - 1 / kappa) * np.arctan(y) + y / kappa)
        Omega0 = 2 / pi * (1 - np.exp(-1 / (1 / y + (y / delta) ** 2 / 3)))
        val = kB * R_D_mix * rho * cp * T / (eta * xi) / (6 * pi) * (Omega - Omega0) / molemass
        if val < 0:
            return 0.0
        return val
    except Exception:
        return 0.0


def s_plus_lambda_plus(x, n1, n2, n3, n4, c):
    return n1 * (x / c) + n2 * (x / c) ** 1.5 + n3 * (x / c) ** 2 + n4 * (x / c) ** 2.5


def fluid_mix_name_cp(material, molefrac):
    mix_name = material[0] + '[' + str(molefrac[0]) + ']'
    ncomp = len(material)
    for ii in range(1, ncomp):
        mix_name = mix_name + '&' + material[ii] + '[' + str(molefrac[ii]) + ']'
    return mix_name


def TC_RES(fluid, MoleFrac, p_Pa, T_K, eos="REFPROP", dilute_source="native", enhancement_viscosity="native"):
    """RES thermal conductivity in W/m/K.

    ADAPTED (1, 2, 3, 6).  All three options default to what the paper does, so the published
    sample table is reproduced with no arguments beyond the state.

    Returns (lambda_res, rho, S_resi, Group_N, lambda0) -- the trailing lambda0 is the addition.
    """
    xitapower = np.array([1, 1.5, 2, 2.5])
    R, kB = [8.314462618, 1.380649e-23]
    N_A = R / kB
    TheEOS = eos  # ADAPTED (1): was TheEOS = "REFPROP"
    (epsilon_kB_K, sigma_nm, parameter_N, xita, Group_N, Tc0_N, Tc_K,
     cp0_k0, cp0_k1, R_D, gamma, xi0, Gamma, qDinv, Tref) = get_paramters(fluid)

    if isinstance(fluid, str):
        Purefluid = True
        fluid_ref = fluid
    else:
        Purefluid = False
        Ncomp = len(fluid)
        fluid_ref = fluid[0]
        for ii in range(1, Ncomp):
            fluid_ref = fluid_ref + '&' + fluid[ii]

    AS = CP.AbstractState(TheEOS, fluid_ref)
    if not Purefluid:
        AS.set_mole_fractions(MoleFrac)
    AS.update(CP.PT_INPUTS, p_Pa, T_K)
    MoleMass = AS.molar_mass()
    S_resi = AS.smolar_residual()
    S_plus = -S_resi / R
    rho = AS.rhomass()
    rhoN = rho / MoleMass * N_A
    mass_kg = MoleMass / N_A

    if Purefluid:

        # ADAPTED (9): comment only, no code change.  Index 0 multiplies T^4 here, matching
        # Dilute_gas_TC.txt in this directory.  CoolProp stores the same numbers ascending and
        # reverses them once on import; both agree.  See the same note in the viscosity module.
        lambda0 = Tc0_N[0, 0] * T_K**4 + Tc0_N[0, 1] * T_K**3 + Tc0_N[0, 2] * T_K**2 + Tc0_N[0, 3] * T_K + Tc0_N[0, 4]

        if p_Pa == 0:
            return lambda0, rho, S_resi, Group_N, lambda0
        lambda_plus_calc = s_plus_lambda_plus(
            S_plus, parameter_N[0, 0], parameter_N[0, 1], parameter_N[0, 2], parameter_N[0, 3], xita
        )
        try:
            lambda_c = Olchowy_critical_enhancement(
                fluid_ref, T_K, rho, R_D, gamma, xi0, Gamma, qDinv, Tref,
                eos=TheEOS, enhancement_viscosity=enhancement_viscosity, MoleFrac=MoleFrac, p_Pa=p_Pa,
            )
        except Exception:
            lambda_c = 0
    else:
        MoleMass_pure = np.zeros(Ncomp)
        mass_kg_pure = np.zeros(Ncomp)
        lambd0_pure = np.zeros(Ncomp)
        para_n = np.zeros(4)
        SumMass = 0
        for icomp in range(Ncomp):
            lambd0_pure[icomp] = (
                Tc0_N[icomp][0] * T_K**4
                + Tc0_N[icomp][1] * T_K**3
                + Tc0_N[icomp][2] * T_K**2
                + Tc0_N[icomp][3] * T_K
                + Tc0_N[icomp][4]
            )
            MoleMass_pure[icomp] = CP.PropsSI('molemass', 'T', 0, 'p', 0, TheEOS + '::' + fluid[icomp])
            mass_kg_pure[icomp] = MoleMass_pure[icomp] / N_A
            SumMass += MoleFrac[icomp] * MoleMass_pure[icomp]
            for ii in range(4):
                para_n[ii] += MoleFrac[icomp] * parameter_N[icomp, ii] / (xita[icomp] ** xitapower[ii])

        MassFrac = np.zeros(Ncomp)
        for icomp in range(Ncomp):
            MassFrac[icomp] = MoleFrac[icomp] * MoleMass_pure[icomp] / SumMass
        mass_kg = geometric_mean(MassFrac, MoleMass_pure) / N_A
        RefMixName = fluid_mix_name_cp(fluid, MoleFrac)
        # ADAPTED (2): the original always attempted the native lookup and fell back to
        # polynomial + Wilke only when it raised.  "polynomial" now skips the attempt outright.
        if dilute_source == "polynomial":
            lambda0 = lambda_dilute_gas_multi_wilke_fast(lambd0_pure, mass_kg_pure, MoleFrac)
        else:
            try:
                lambda0 = CP.PropsSI('L', 'T', T_K, 'P', 0.1, TheEOS + '::' + RefMixName)
            except Exception:
                lambda0 = lambda_dilute_gas_multi_wilke_fast(lambd0_pure, mass_kg_pure, MoleFrac)

        if p_Pa == 0:
            return lambda0, rho, S_resi, Group_N, lambda0
        lambda_plus_calc = s_plus_lambda_plus(S_plus, para_n[0], para_n[1], para_n[2], para_n[3], 1)
        try:
            lambda_c = Olchowy_critical_enhancement_mix(
                RefMixName, MoleFrac, T_K, rho, R_D, gamma, xi0, Gamma, qDinv, Tref,
                eos=TheEOS, enhancement_viscosity=enhancement_viscosity, components=fluid, p_Pa=p_Pa,
            )
        except Exception:
            lambda_c = 0

    lambda_res = (lambda_plus_calc / S_plus ** (2 / 3)) * (rhoN ** (2 / 3) * kB * (kB * T_K / mass_kg) ** (1 / 2)) + lambda0 + lambda_c

    return lambda_res, rho, S_resi, Group_N, lambda0


def write_sample_table(out_path):
    """ADAPTED (7): the `if __name__ == '__main__'` driver, callable.

    Regenerates Table_S5_SI.txt from Samples_pure_fluids.txt exactly as the published script
    does.  Used by dev/RES_reference_run.py --self-test to confirm that the modifications above
    left the model itself untouched.
    """
    FID_pure = open(out_path, "w")
    FID_pure.write("           Material        T/K      p/kPa    den/kg/m3  s_resi/JmolK  TC_EXP/Wm-1K-1 TC_RES TC_REF \n")
    PureData = pd.read_csv(os.path.join(_HERE, "Samples_pure_fluids.txt"), sep=r'\s+', header=0, skipinitialspace=True)
    for ifluid in range(len(PureData['Material'])):
        p_Pa = PureData['p/kPa'][ifluid] * 1000
        T_K = PureData['T/K'][ifluid]
        try:
            TheEOS = "REFPROP"
            AS = CP.AbstractState(TheEOS, PureData['Material'][ifluid])
            AS.update(CP.PT_INPUTS, p_Pa, T_K)
            lambda_REF = AS.conductivity()
        except Exception:
            lambda_REF = 0
        lambda_res, rho, S_resi, Group_N, _lambda0 = TC_RES(PureData['Material'][ifluid], 1, p_Pa, T_K)
        FID_pure.write(
            "%21s  %8.3f %11.6e %10.3f %10.3f %10.5f %10.5f %10.5f\n"
            % (
                PureData['Material'][ifluid],
                T_K,
                p_Pa * 1e-3,
                rho,
                S_resi,
                PureData['TC_exp/Wm-1K-1'][ifluid],
                # ADAPTED (8): TC_RES returns a 1-element array for a pure fluid, and
                # "%f" % array raises from numpy 2.0 on where it used to work.
                float(np.ravel(lambda_res)[0]),
                lambda_REF,
            )
        )
    FID_pure.close()


def write_binary_sample_table(out_path):
    """ADAPTED (7): the binary-mixture half of the `__main__` driver, callable.

    Regenerates Table_S6_SI.txt.  Worth having alongside the pure-fluid table because it is the
    only thing that exercises the Wilke rule for the dilute term, the mole-fraction averaging of
    the residual coefficients, and the mixture branch of the critical enhancement.
    """
    with open(out_path, "w") as FID_bin:
        FID_bin.write(
            "            Comp1          Comp2    GroupN1 N2    MF1       MF2      T/K    p/kPa  den/kg/m3  "
            "s_resi/JmolK  TC_EXP/Wm-1K-1 TC_RES TC_REF \n"
        )
        BinaryData = pd.read_csv(os.path.join(_HERE, "Samples_binaries.txt"), sep=r'\s+', header=0, skipinitialspace=True)
        for imix in range(len(BinaryData['Comp1'])):
            Mixture = [BinaryData['Comp1'][imix], BinaryData['Comp2'][imix]]
            MoleFrac = [BinaryData['MF1'][imix], BinaryData['MF2'][imix]]
            p_Pa = BinaryData['p/kPa'][imix] * 1e3
            T_K = BinaryData['T/K'][imix]
            try:
                TheEOS = "REFPROP"
                AS = CP.AbstractState(TheEOS, Mixture[0] + '&' + Mixture[1])
                AS.set_mole_fractions(MoleFrac)
                AS.update(CP.PT_INPUTS, p_Pa, T_K)
                lambda_REF = AS.conductivity()
            except Exception:
                lambda_REF = 0
            lambda_res, rho, S_resi, Group_N, _lambda0 = TC_RES(Mixture, MoleFrac, p_Pa, T_K)
            FID_bin.write(
                "%17s %17s %4.0f %4.0f %8.6f %8.6f %8.3f %10.6f %8.3f %8.3f %14.3f %10.5f %10.5f \n"
                % (
                    BinaryData['Comp1'][imix],
                    BinaryData['Comp2'][imix],
                    Group_N[0],
                    Group_N[1],
                    BinaryData['MF1'][imix],
                    BinaryData['MF2'][imix],
                    BinaryData['T/K'][imix],
                    BinaryData['p/kPa'][imix],
                    rho,
                    S_resi,
                    BinaryData['TC_EXP/Wm-1K-1'][imix],
                    lambda_res,
                    lambda_REF,
                )
            )
