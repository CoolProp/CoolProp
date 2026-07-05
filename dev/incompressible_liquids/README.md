# Incompressible fluids: fitting pipeline

This directory turns tabulated fluid property data into the coefficient
files that CoolProp's `INCOMP` backend ships. You do **not** need to be a
programmer to add a fluid: you copy an existing fluid definition, paste in
your data, and run one script.

## Layout

| Path | What it is |
|---|---|
| `CPIncomp/PureFluids.py`, `SolutionFluids.py`, ... | The fluid definitions (one Python class per fluid) |
| `CPIncomp/ExampleObjects.py` | Commented templates for every fluid kind |
| `CPIncomp/data/` | Tabulated source data for the file-backed fluids |
| `json/` | **The output.** One coefficient file per fluid; these are baked into the C++ library |
| `all_incompressibles.py` | The driver script that fits everything and writes `json/` |
| `test_json_sanity.py` | Checks that no unfitted placeholder coefficients are shipped |
| `test_fitting_regression.py` | Re-fits three fluids and compares against the committed JSON |

## Prerequisites

Python 3 with `numpy` and `scipy` (see `requirements.txt`):

```bash
pip install -r requirements.txt
```

That is enough to fit fluids and write JSON. `matplotlib` and a built
CoolProp Python package are only needed if you also want the fitting
reports, plots and documentation tables.

## Adding a new pure fluid, step by step

1. Open `CPIncomp/PureFluids.py` and copy an existing class (or start from
   `PureExample` in `CPIncomp/ExampleObjects.py`). Rename the class and give
   it a unique `self.name`.

2. Paste your data as numpy arrays. All arrays must have the same length as
   `self.temperature.data`, and everything is in SI units:

   | Attribute | Unit |
   |---|---|
   | `self.temperature.data` | K |
   | `self.density.data` | kg/m³ |
   | `self.specific_heat.data` | J/(kg·K) |
   | `self.conductivity.data` | W/(m·K) |
   | `self.viscosity.data` | Pa·s |
   | `self.saturation_pressure.data` | Pa |

   Use `np.nan` for points you do not have. If you have no data at all for a
   property, do not set it -- the fluid will then cleanly report that the
   property is not available instead of returning made-up numbers.

3. For each property you filled in, mark the data source, e.g.
   `self.density.source = self.density.SOURCE_DATA`, and set `self.Tmin`,
   `self.Tmax`, `self.TminPsat` (minimum temperature at which the vapour
   pressure fit is valid), `self.description` and `self.reference`
   (a citation for where the data came from; a BibTeX key from
   `CoolPropBibTeXLibrary.bib` if possible). Finish with `self.reshapeAll()`.

4. Run the fitter from this directory (the flags skip the optional reports,
   plots, tables and statistics):

   ```bash
   python all_incompressibles.py -nr -ns -nt -nst
   ```

   Watch the output for your fluid. A message like
   `MyFluid: could not fit viscosity, marking it as not defined` means the
   fit did not succeed (usually: no data was set) and the property was
   left out on purpose.

5. Check the result: `json/<name>.json` should now exist and contain
   coefficients for each property you provided. Run the sanity tests:

   ```bash
   python -m pytest test_json_sanity.py test_fitting_regression.py
   ```

6. Rebuild CoolProp. The build regenerates the embedded fluid library
   (`dev/generate_headers.py` bundles everything in `json/` into a C++
   header) so your fluid is available as `INCOMP::<name>`:

   ```python
   PropsSI("D", "T", 300, "P", 101325, "INCOMP::MyFluid")
   ```

7. Commit your edited fluid module and the new `json/<name>.json` file.

Binary mixtures (water solutions and the like) work the same way but live in
`CPIncomp/SolutionFluids.py` and need a concentration vector as the second
data axis -- see `SolutionExample` in `CPIncomp/ExampleObjects.py` and the
"Adding New Fluids" section of the online documentation
(`Web/fluid_properties/Incompressibles.rst`).

## Chebyshev caloric entries

Every fluid's JSON also carries optional `density_cheb` /
`specific_heat_cheb` entries (`type: "chebyshev"`): Chebyshev-in-T
coefficients on the explicit `Trange`, one column per `(x - xbase)` power.
They are produced automatically by the writer — fitted from your raw data
when the fluid has any, otherwise an exact basis conversion of the
polynomial fit — and were introduced because enthalpy/entropy integrals of
a Chebyshev cp fit are exact and singularity-free (see
`NOTES_thermodynamic_consistency.md` and `CPIncomp/ChebyshevFits.py`).
`add_chebyshev_entries.py` is the one-shot migration tool that added them
to the committed files without refitting the polynomial entries;
`test_chebyshev_entries.py` guards their quality (positivity across the
whole domain, agreement with data and with the polynomial fits).

## Notes

- The fit is centred around `Tbase`/`xbase` (defaults: the midpoint of your
  data range). The enthalpy/entropy the backend later reports are derived
  exactly from the fitted density and heat-capacity polynomials.
- Four fluids (Air, Acetone, Ethanol, Hexane) are sampled from CoolProp's
  own equations of state, so *re-fitting those four* needs the CoolProp
  Python package installed; without it they are skipped and their committed
  JSON stays as-is. All other fluids fit from the data in this directory.
- Fit output can drift slightly across numpy/scipy versions; that is what
  `test_fitting_regression.py` guards. If it fails after a fluid addition
  that did not touch the three reference fluids, suspect your environment,
  not your fluid.
- Solutions carry no enthalpy/entropy of mixing (documented limitation, see
  `NOTES_mixing_models.md`).
