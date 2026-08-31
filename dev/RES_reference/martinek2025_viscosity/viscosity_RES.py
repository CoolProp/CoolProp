"""Reference implementation of the RES viscosity model of Martinek et al. (2025).

    Martinek, Bell, Herzog, Richter, Yang,
    "Entropy Scaling of Viscosity IV -- Application to 124 Industrially Important Fluids",
    J. Chem. Eng. Data 70, 727-742 (2025).  https://doi.org/10.1021/acs.jced.4c00451

Redistributed and modified under the Creative Commons Attribution 4.0 International licence
(CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), which the article and its supporting
information are published under.  The accompanying data files in this directory
(Dilute_gas_viscosity.txt, RES_Parameter.txt, Samples_*.txt) are the authors' own and are
reproduced VERBATIM.

CC BY 4.0 requires that changes be indicated.  This file is `code_SI/main.py` from the supporting
information with the following modifications, each marked `# ADAPTED:` at the point of change:

  1. The equation of state is a parameter (`eos=`) rather than the hardcoded "REFPROP".
  2. The dilute-gas term is selectable (`dilute_source=`) rather than always REFPROP's native
     eta0 for pure fluids.
  3. Data files are located relative to THIS FILE rather than the process working directory.
  4. `viscosity_RES` additionally returns the dilute-gas term it used, so a caller can weight a
     comparison by how much of the viscosity is actually residual.
  5. The `if __name__ == '__main__'` sample driver is split into `write_sample_table()` and
     `write_binary_sample_table()`, which do the same thing on demand.
  6. This docstring.  `import sys` is kept: get_paramters() still calls sys.exit().
  7. One added comment in vis_dilute_gas_fast() recording the coefficient order.  No code change.

Every default reproduces the published behaviour, so running with defaults regenerates both
Samples_pure_fluids_output.txt and Samples_binaries_output.txt byte-for-byte -- which is how
dev/RES_reference_run.py --self-test checks that these modifications did not change the model.
"""

import os
import sys

import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd

# ADAPTED (3): the originals read these by bare filename from the working directory.
_HERE = os.path.dirname(os.path.realpath(__file__))


def get_paramters(AllMaterial):

    Dilute_gas_viscosity = pd.read_csv(
        os.path.join(_HERE, "Dilute_gas_viscosity.txt"), skipinitialspace=True, header=1, delimiter=';'
    )
    RES_Parameter = pd.read_csv(os.path.join(_HERE, "RES_Parameter.txt"), sep=r'\s+', header=0, skipinitialspace=True)
    Len_all = len(Dilute_gas_viscosity['Material'])

    if isinstance(AllMaterial, str):
        nfluid = 1
    else:
        nfluid = len(AllMaterial)

    parameter_N = np.zeros((nfluid, 4))
    vis0_N = np.zeros((nfluid, 5))
    group_n = np.zeros(nfluid)
    xita = np.zeros(nfluid)
    for ifluid in range(nfluid):
        if isinstance(AllMaterial, str):
            material = AllMaterial
        else:
            material = AllMaterial[ifluid]
        for ii in range(Len_all):
            if material.upper() == Dilute_gas_viscosity['Material'][ii].upper():
                if Dilute_gas_viscosity['Material'][ii].upper() != RES_Parameter['Material'][ii].upper():
                    print(material + ' are not the same in Fluid_Constants.txt and RES_Parameter.txt')
                    sys.exit()

                vis0_N[ifluid, :] = [
                    Dilute_gas_viscosity['n0'][ii],
                    Dilute_gas_viscosity['n1'][ii],
                    Dilute_gas_viscosity['n2'][ii],
                    Dilute_gas_viscosity['n3'][ii],
                    Dilute_gas_viscosity['n4'][ii],
                ]
                group_n[ifluid] = RES_Parameter['GroupN'][ii]
                xita[ifluid] = RES_Parameter['xita'][ii]
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
            if ii == Len_all - 1:
                print(material + ' not found in data/Fluid_Constants.txt')
                sys.exit()
    return vis0_N, parameter_N, xita, group_n


def vis_dilute_gas_fast(vis0_n, T_K):
    # ADAPTED (7): comment only, no code change.  Dilute_gas_viscosity.txt in this directory heads
    # its columns n0..n4 with n0 the T^4 coefficient, so index 0 multiplies T^4 here.  CoolProp
    # stores the same five numbers ascending and evaluates them that way, reversing them once on
    # import in dev/convert_RES_csv_to_json.py::dilute_ascending().  Both are self-consistent and
    # agree exactly; only a reader who mixes the two conventions gets a result off by ~1e9.
    ncomp = np.size(vis0_n, 0)
    vis0 = np.zeros(ncomp)
    for icomp in range(ncomp):
        vis0[icomp] = (
            vis0_n[icomp][0] * T_K**4
            + vis0_n[icomp][1] * T_K**3
            + vis0_n[icomp][2] * T_K**2
            + vis0_n[icomp][3] * T_K
            + vis0_n[icomp][4]
        ) * 1e-6
    return vis0


def vis_dilute_gas_bin_wilke_fast(vis0, mass, mf):
    fai00 = 1
    fai11 = 1
    fai01 = (1 + np.sqrt(vis0[0] / vis0[1]) * (mass[1] / mass[0]) ** (1 / 4)) ** 2 / np.sqrt(8 * (1 + mass[0] / mass[1]))
    fai10 = (1 + np.sqrt(vis0[1] / vis0[0]) * (mass[0] / mass[1]) ** (1 / 4)) ** 2 / np.sqrt(8 * (1 + mass[1] / mass[0]))
    sum0 = mf[0] * fai00 + mf[1] * fai01
    sum1 = mf[0] * fai10 + mf[1] * fai11
    vis_mix = mf[0] * vis0[0] / sum0 + mf[1] * vis0[1] / sum1
    return vis_mix


def s_plus_eta_plus(x, n1, n2, n3, c):
    return n1 * (x / c) ** 1.8 + n2 * (x / c) ** 2.4 + n3 * (x / c) ** 2.8


def viscosity_RES(fluid, MoleFrac, p_Pa, T_K, eos="REFPROP", dilute_source="native"):
    """RES viscosity in Pa.s.

    ADAPTED (1, 2, 4).  `eos` and `dilute_source` both default to what the paper does, so the
    published sample table is reproduced with no arguments beyond the state.

      eos            -- CoolProp backend name to evaluate the equation of state on.
      dilute_source  -- "native" uses the backend's own eta0 for a PURE fluid, as the paper does;
                        "polynomial" always uses the fitted polynomial in this directory.  A
                        mixture uses the polynomial + Wilke either way, which is also what the
                        paper does.

    Returns (vis_res, rho, S_resi, Group_N, vis0) -- the trailing vis0 is the addition.
    """
    xitapower = np.array([1.8, 2.4, 2.8])
    R, kB = [8.314462618, 1.380649e-23]
    N_A = R / kB
    TheEOS = eos  # ADAPTED (1): was TheEOS = "REFPROP"
    vis0_N, Parameter_N, Xita, Group_N = get_paramters(fluid)

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
        # ADAPTED (2): the original always attempted the native lookup and fell back to the
        # polynomial only when it raised.  "polynomial" now skips the attempt outright.
        if dilute_source == "polynomial":
            vis0 = vis_dilute_gas_fast(vis0_N, T_K)[0]
        else:
            try:
                vis0 = CP.PropsSI('V', 'T', T_K, 'P', 1e-9, TheEOS + '::' + fluid)
            except Exception:
                vis0 = vis_dilute_gas_fast(vis0_N, T_K)
        if p_Pa == 0:
            return vis0, rho, S_resi, Group_N, vis0
        vis_plus_calc = np.exp(s_plus_eta_plus(S_plus, Parameter_N[0, 0], Parameter_N[0, 1], Parameter_N[0, 2], Xita)) - 1
    else:
        Ncomp = len(fluid)
        MoleMass_pure = np.zeros(Ncomp)
        mass_kg_pure = np.zeros(Ncomp)
        para_n = np.zeros(3)
        for icomp in range(Ncomp):
            MoleMass_pure[icomp] = CP.PropsSI('molemass', 'T', 0, 'p', 0, TheEOS + '::' + fluid[icomp])
            mass_kg_pure[icomp] = MoleMass_pure[icomp] / N_A
            for ii in range(3):
                para_n[ii] += MoleFrac[icomp] * Parameter_N[icomp, ii] / (Xita[icomp] ** xitapower[ii])
        vis_pure0 = vis_dilute_gas_fast(vis0_N, T_K)
        vis0 = vis_dilute_gas_bin_wilke_fast(vis_pure0, mass_kg_pure, MoleFrac)
        if p_Pa == 0:
            return vis0, rho, S_resi, Group_N, vis0
        vis_plus_calc = np.exp(s_plus_eta_plus(S_plus, para_n[0], para_n[1], para_n[2], 1)) - 1
    vis_res = (vis_plus_calc / S_plus ** (2 / 3)) * (rhoN ** (2 / 3) * (mass_kg * kB * T_K) ** (1 / 2)) + vis0
    return vis_res, rho, S_resi, Group_N, vis0


def write_sample_table(out_path):
    """ADAPTED (5): the `if __name__ == '__main__'` driver, callable.

    Regenerates Samples_pure_fluids_output.txt from Samples_pure_fluids_input.txt exactly as the
    published script does.  Used by dev/RES_reference_run.py --self-test to confirm that the
    modifications above left the model itself untouched.
    """
    FID_pure = open(out_path, "w")
    FID_pure.write("           Material        T/K      p/MPa    den/kg/m3  s_resi/JmolK  vis_exp/uPas vis_res vis_REF \n")
    PureData = pd.read_csv(os.path.join(_HERE, "Samples_pure_fluids_input.txt"), sep=r'\s+', header=0, skipinitialspace=True)
    for ifluid in range(len(PureData['Material'])):
        p_Pa = PureData['p/MPa'][ifluid] * 1000000
        T_K = PureData['T/K'][ifluid]
        try:
            TheEOS = "REFPROP"
            AS = CP.AbstractState(TheEOS, PureData['Material'][ifluid])
            AS.update(CP.PT_INPUTS, p_Pa, T_K)
            vis_REF = AS.viscosity()
        except Exception:
            vis_REF = 0
        vis_res, rho, S_resi, Group_N, _vis0 = viscosity_RES(PureData['Material'][ifluid], 1, p_Pa, T_K)
        FID_pure.write(
            "%21s  %8.3f %11.6e %10.3f %10.3f %10.3f %10.3f %10.3f \n"
            % (
                PureData['Material'][ifluid],
                T_K,
                p_Pa * 1e-6,
                rho,
                S_resi,
                PureData['vis_exp/uPas'][ifluid],
                vis_res[0] * 1e6,
                vis_REF * 1e6,
            )
        )
    FID_pure.close()


def write_binary_sample_table(out_path):
    """ADAPTED (5): the binary-mixture half of the `__main__` driver, callable.

    Regenerates Samples_binaries_output.txt.  Worth having alongside the pure-fluid table because
    it is the only thing that exercises the Wilke rule for the dilute term and the mole-fraction
    averaging of the residual coefficients.
    """
    with open(out_path, "w") as FID_bin:
        FID_bin.write(
            "            Comp1          Comp2    GroupN1 N2    MF1       MF2      T/K    p/MPa  den/kg/m3  "
            "s_resi/JmolK  vis_exp/uPas vis_res   vis_REF \n"
        )
        BinaryData = pd.read_csv(os.path.join(_HERE, "Samples_binaries_input.txt"), sep=r'\s+', header=0, skipinitialspace=True)
        for imix in range(len(BinaryData['Comp1'])):
            Mixture = [BinaryData['Comp1'][imix], BinaryData['Comp2'][imix]]
            MoleFrac = [BinaryData['MF1'][imix], BinaryData['MF2'][imix]]
            p_Pa = BinaryData['p/MPa'][imix] * 1e6
            T_K = BinaryData['T/K'][imix]
            try:
                TheEOS = "REFPROP"
                AS = CP.AbstractState(TheEOS, Mixture[0] + '&' + Mixture[1])
                AS.set_mole_fractions(MoleFrac)
                AS.update(CP.PT_INPUTS, p_Pa, T_K)
                vis_REF = AS.viscosity()
            except Exception:
                vis_REF = 0
            vis_res, rho, S_resi, Group_N, _vis0 = viscosity_RES(Mixture, MoleFrac, p_Pa, T_K)
            FID_bin.write(
                "%17s %17s %4.0f %4.0f %8.6f %8.6f %8.3f %10.6f %8.3f %8.3f %14.3f %10.3f %10.3f  \n"
                % (
                    BinaryData['Comp1'][imix],
                    BinaryData['Comp2'][imix],
                    Group_N[0],
                    Group_N[1],
                    BinaryData['MF1'][imix],
                    BinaryData['MF2'][imix],
                    BinaryData['T/K'][imix],
                    BinaryData['p/MPa'][imix],
                    rho,
                    S_resi,
                    BinaryData['vis_exp/uPas'][imix],
                    vis_res * 1e6,
                    vis_REF * 1e6,
                )
            )
