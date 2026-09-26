# AD spike: autodiff vs num-dual-style duals vs mcx / complex step / FD

This is a throwaway spike, not built by CoolProp's CMake. It asks one question: if CoolProp grew
AD-based models, would it inherit teqp's compile-time and binary-size scaling from
`autodiff`? Or does a num-dual (FeOs) style dual library avoid it?

## Setup

- **Model:** PC-SAFT (hard chain + dispersion, Gross & Sadowski 2001), written once in
  `pcsaft_generic.hpp` and generic over the number type. Every method runs identical model
  code. It is dependency-free (`std::array`, so no heap in the timed path). It matches teqp's
  `get_Ar00/Ar01/Ar11` to 1 ulp.
- **Quantities** (teqp conventions): `Ar0n` = {A00..A03} in ρ; `Arn0` = {A00, A10, A20} in 1/T;
  `Ar11`; `gradx` = ∂αr/∂x_i with independent compositions.
- **Methods,** one TU each in `methods/`:

| method | what |
|---|---|
| `double` | value only (floor for compile time and run time) |
| `ad_real` | autodiff as teqp uses it: `Real<N>` for the pure ρ and 1/T derivatives, `dual2nd` for A11, one seeded `dual` pass per component |
| `ad_dual` | autodiff with nested expression-template duals only (`dual3rd`, `dual2nd`) |
| `numdual` | `numdual.hpp`, ~200 lines imitating Rust num-dual: eager structs `Dual3`, `Dual2`, `HyperDual`, `DualVec<N>`, no expression templates |
| `mcx` | multicomplex step (teqp's alternative backend) |
| `cstep` | `std::complex` complex step, first derivatives only |
| `fd` | central finite differences, h ~ ε^(1/(n+2)) |

- **Reference:** an independent pure-Python PC-SAFT evaluated at 60 digits in mpmath and
  differentiated with `mp.diff`.
- **Scaling:** the same TU compiled with `KMODELS` = 1, 4, 16 distinct model types (`PCSAFTModel<Tag>`).
  This mimics adding models to a library.

Run it with `python3 run_spike.py` (it needs the autodiff and mcx headers; the defaults point into a
teqp checkout).

## Results (Apple clang 21, M-series, -O2)

### Accuracy: max relative error vs 60-digit reference

| method | A01 | A02 | A03 | A10 | A20 | A11 | ∂/∂x |
|---|---|---|---|---|---|---|---|
| all AD (ad_real, ad_dual, numdual, mcx) | ≤1e-15 | ≤1e-14 | ≤2e-14 (dense) | ≤1e-15 | ≤2e-15 | ≤4e-16 | ≤1e-15 |
| cstep | ≤4e-15 | — | — | ≤1e-15 | — | — | ≤1e-15 |
| fd | 4e-11 | 2e-7 | 8e-6 … 1e-1 | 2e-12 | 1e-7 | 9e-9 | 2e-10 |

All the exact methods are indistinguishable at the level of double roundoff. At the dilute gas
state, A03 ≈ −1.6e-6 carries ~1e-11 relative error for **every** exact method, which comes from
the conditioning of αr in double and not from the AD method. FD loses 5–10 digits there.

### Speed: ns per call (3-component mixture, dense)

| method | Ar0n (to 3rd) | Arn0 (to 2nd) | Ar11 | gradx (N=3) |
|---|---|---|---|---|
| double (1 eval) | ~100 | | | |
| ad_real | 209 | 218 | 299 | 351 |
| ad_dual | 606 | 298 | 300 | 352 |
| **numdual** | **150** | **153** | **187** | **154** |
| mcx | 28 000 | 27 000 | 27 000 | 47 000 |
| cstep | 325 (1st only) | 521 (1st only) | — | 1135 |
| fd | 733 | 442 | 340 | 495 |

- The pure-component states show the same ordering. There, numdual is about 1.3–1.6× faster than
  teqp-style autodiff.
- Most of the gain comes from purpose-built storage. `Dual3` holds 4 numbers, whereas `dual3rd` holds
  2³ = 8, so it is 4× slower. `HyperDual` holds 4 coefficients, while `dual2nd` holds 4 plus
  expression-template machinery.
- `DualVec<N>` computes the whole composition gradient in one pass instead of N passes.
- `std::complex` is slow here because its division carries C99 inf/NaN recovery.
  `-fcx-limited-range` would help.

### Compile time and code size per TU

| method | K=1 | K=4 | K=16 | **s / extra model** | **kB text / extra model** | -O0 K=1 text |
|---|---|---|---|---|---|---|
| double | 0.26 s | 0.30 | 0.50 | 0.02 | 3.7 | 3 kB |
| ad_real | 0.92 s | 1.44 | 3.46 | **0.17** | **21** | 89 kB |
| ad_dual | 0.92 s | 1.41 | 3.49 | **0.17** | **18** | 75 kB |
| numdual | 0.43 s | 0.94 | 3.02 | **0.17** | **22** | 30 kB |
| mcx | 1.24 s | 2.77 | 9.11 | 0.52 | 71 | 83 kB |
| cstep | 0.56 s | 0.89 | 2.23 | 0.11 | 16 | 17 kB |
| fd | 0.26 s | 0.32 | 0.57 | 0.02 | 4.8 | 5 kB |

Header parsing alone costs: autodiff `dual.hpp` 0.47 s, `real.hpp` 0.46 s, `numdual.hpp` 0.24 s,
`<complex>` 0.45 s. The empty-TU floor is about 0.23 s.

## Conclusions

1. **The dual library changes the intercept, not the slope.** At K=1, autodiff costs about 0.5 s more
   per TU than numdual. Roughly half of that is header parsing and half is expression-template
   instantiation. But the marginal cost per additional model is identical: 0.17 s and about 20 kB
   of optimized code per model, for this set of 4 derivative kinds. The cost scales as
   (models × derivative kinds × number types) × (model size), and every forward-mode
   operator-overloading approach pays that. Rust monomorphization pays it too, which is why
   FeOs builds are also slow.
2. **What num-dual really buys is speed, plus debug builds that aren't bloated.** It runs 1.3–2.3×
   faster at runtime, with the largest gain for composition gradients. `-O0` text is 3× smaller
   because there are no unevaluated expression trees. Accuracy is identical.
3. **The levers for compile time are architectural, not the choice of dual type:**
   - Fix a small, closed set of number types, e.g. `double`, `Dual3` (ρ), `Dual2` (1/T),
     `HyperDual` (1/T, ρ), `DualVec<N>` (x). Don't let a template-int family like
     `HigherOrderDual<iT+iD+iX>` mint new types per derivative request.
   - Put each model's instantiations in its own `.cpp` via explicit instantiation. The cost then
     becomes linear, is paid once, and builds in parallel, instead of being paid in every TU that
     includes the model.
   - Keep the AD types behind CoolProp's existing runtime boundary (`AbstractState`), so no
     public header ever sees them.
4. **mcx** is a fine reference oracle, but it is about 100× slower and about 3× the compile cost per
   model. **FD** is not acceptable past the first derivative.
