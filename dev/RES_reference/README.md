# Vendored RES reference implementations

The three papers CoolProp's RES transport model comes from publish their reference code and
parameter tables as supporting information.  All of it is copied here — the two code
implementations adapted as described below, every data table verbatim — and used as the
independent truth the CoolProp implementation is checked against.

| Directory | Source | DOI |
|---|---|---|
| `martinek2025_viscosity/` | Martinek, Bell, Herzog, Richter, Yang, *Entropy Scaling of Viscosity IV — Application to 124 Industrially Important Fluids*, J. Chem. Eng. Data **70**, 727–742 (2025) | [10.1021/acs.jced.4c00451](https://doi.org/10.1021/acs.jced.4c00451) |
| `li2024_conductivity/` | Li, Duan, Yang, *Linking Thermal Conductivity to Equations of State Using the Residual Entropy Scaling Theory*, Ind. Eng. Chem. Res. **63**, 18160–18175 (2024) | [10.1021/acs.iecr.4c02946](https://doi.org/10.1021/acs.iecr.4c02946) |
| `yang2025/` | Yang, *Viscosity and Thermal Conductivity Models of 151 Common Fluids Based on Residual Entropy Scaling and Cubic Equations of State*, ACS Omega (2025) | [10.1021/acsomega.4c10815](https://doi.org/10.1021/acsomega.4c10815) |

## Licence

All articles and their supporting information are open access under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/),
which permits redistribution and modification for any purpose, including commercially, provided
the source is attributed and **changes are indicated**.

Attribution is in the table above and repeated in each module's docstring.  Changes are indicated
in three places: the module docstring lists them, every changed line carries an `# ADAPTED:`
marker, and the section below summarises them.

All `.txt` data files are reproduced **verbatim**.  Only the two `.py` files are modified;
`yang2025/` is data only, and carries no code to modify.

`yang2025/` holds the four coefficient sets CoolProp's cubic backends can use, taken from
*Supporting Information/Parameters/* in the article's `si_002` archive.  The same archive also
ships PTV and YFR fits, which are omitted because CoolProp has no backend for either.

## What was changed, and why

Both modules hardcode a single configuration: REFPROP as the equation of state, and each paper's
own choice of where the dilute-gas term and the enhancement viscosity come from.  CoolProp needs
the same model evaluated on a different equation of state and with the dilute term sourced the
way CoolProp sources it, otherwise a comparison measures those choices rather than the thing
under test.

Every change is a new keyword argument whose **default reproduces the published behaviour**:

| Option | Default (= the paper) | Used by the harness |
|---|---|---|
| `eos=` | `"REFPROP"` | `"HEOS"` for the CoolProp-side comparison |
| `dilute_source=` | `"native"` — backend's own `eta0` / `lambda0` | `"polynomial"` — the fitted polynomial, as CoolProp uses |
| `enhancement_viscosity=` (conductivity only) | `"native"` — backend's own viscosity | `"res"` — the RES viscosity, as CoolProp uses |
| `zero_ind_fit=` (conductivity only) | `"published"` — trust the `ind_fit` flag | `"global"` — ignore an all-zero individual row |

`zero_ind_fit=` is the one option that corrects a source file rather than re-aiming it.  The
tables carry two coefficient sets per fluid, an individual fit and a group fit, chosen by an
`ind_fit` flag.  Exactly one row sets that flag over four all-zero individual coefficients:
NEOPENTN, in `li2024_conductivity/RES_Parameter.txt` and in both of Yang 2025's conductivity
tables — and Li's published paper lists that fluid with `ind_fit = 0`.  The flag is a
transcription error in the supporting information; followed literally it deletes the residual
term, which is **-89%** against REFPROP's own conductivity in the liquid.  CoolProp ships the
group row the paper intends and the harness pins this option to match.  The default still follows
the flag, which is why the sample tables below still regenerate byte-for-byte.

Plus five mechanical changes that cannot affect a result: data files resolved relative to the
module instead of the working directory; `delim_whitespace=` replaced by `sep=r'\s+'` (removed in
pandas 2.2); the `__main__` sample driver split into callable `write_sample_table()` and
`write_binary_sample_table()`; both entry points additionally return the dilute-gas term they
used, so the report can weight each comparison by how much of the property is actually residual;
and Li's sample driver formats through `float()`, because `"%f" % one_element_array` raises from
numpy 2.0 onwards where it used to work.

Deliberately **not** changed: Li's Olchowy functions use `kB = 1.38064852e-23` where the rest of
the same file — and CoolProp — use `1.380649e-23`.  That is 1.6e-7 relative on the enhancement
term alone.  Correcting the authors' numerics would make this less of a reference.

## Checking the copy is faithful

Because the defaults are the published behaviour, the vendored code must still reproduce the
papers' own published sample tables:

```bash
python dev/RES_reference_run.py --self-test
```

This regenerates all four published sample tables -- pure fluids and binary mixtures, for both
properties -- and compares them byte-for-byte against the copies shipped in the supporting
information.  The binary tables are the ones that matter most here: they are the only thing that
exercises the Wilke rule, the mole-fraction averaging of the residual coefficients, and the
mixture branch of the critical enhancement.  Needs REFPROP, since that is what the papers ran on.
Run it after touching either module.
