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

---

# Experiment 2: the whole half-matrix A_ij (i + j ≤ N) in one call

This models CoolProp's approach of filling every derivative up to 4th order at once
(`HelmholtzDerivatives`). The half-matrix has 15 entries for N=4. Run it with
`python3 run_all.py`; the code is in `methods_all/`.

| method | how |
|---|---|
| `taylor2` | `taylor2.hpp`: a bivariate truncated Taylor type carrying all (N+1)(N+2)/2 coefficients in **one pass**, in scaled variables (1/T₀)(1+u), ρ₀(1+v), so A_ij = i!j!·c_ij |
| `polar` | polarization: **N+1 passes** of autodiff's univariate `Real<N>` along directions (1, s_m) at Chebyshev nodes s_m, then a Vandermonde solve per degree recovers the mixed terms |
| `teqp` | teqp's scheme: `Real<N>` in ρ, `Real<N>` in 1/T, plus **one `HigherOrderDual<i+j>` pass per mixed (i,j)**; that is 6 extra passes at N=4 |

### Speed: ns per call for the whole triangle

| | N=1 (3 values) | N=2 (6) | N=3 (10) | N=4 (15) |
|---|---|---|---|---|
| **pure propane** (1 plain evaluation = 46 ns) | | | | |
| taylor2 | 84 | 167 | 358 | **789** |
| polar | 197 | 367 | 772 | 1474 |
| teqp | 169 | 438 | 1764 | 7004 |
| **C1/C2/C3 mixture** (1 plain evaluation = 88 ns) | | | | |
| taylor2 | 134 | 269 | 708 | **1656** |
| polar | 276 | 624 | 1311 | 2484 |
| teqp | 254 | 669 | 2988 | 12616 |

For comparison, requesting only what you need with the purpose-built types costs 112–153 ns for
A00..A03 (`Dual3`) and 126–189 ns for A11 (`HyperDual`).

### Accuracy: max relative error vs 60-digit mpmath, by total order k

| state | method | k=0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|
| dense (liq, mixture) | all three | ≤1e-15 | ≤2e-15 | ≤3e-15 | ≤1e-14 | ≤2e-14 |
| dilute gas | taylor2 / teqp | 1e-15 | 2e-16 | ≤4e-14 | ≤3e-11 | 1e-9 |
| dilute gas | polar | 1e-15 | 4e-16 | 2e-13 | 5e-11 | 3e-9 |

The dilute-gas losses are in A03 ≈ −1.6e-6 and A04 ≈ −1.5e-7, which are tiny numbers built
from cancellation. All the exact methods lose those digits equally. Polarization's Vandermonde
solve costs a further factor of about 2–5.

### Compile cost: one TU, only order N instantiated

| method | N | K=1 | K=16 | s / extra model | kB / extra model |
|---|---|---|---|---|---|
| taylor2 | 2 | 0.74 s / 8 kB | 1.48 s / 126 kB | 0.05 | 8 |
| taylor2 | 4 | 1.28 s / 28 kB | 5.69 s / 424 kB | 0.29 | 26 |
| polar | 2 | 0.56 s / 8 kB | 1.38 s / 113 kB | 0.05 | 7 |
| polar | 4 | 0.64 s / 14 kB | 2.36 s / 213 kB | **0.11** | **13** |
| teqp | 2 | 0.78 s / 21 kB | 3.01 s / 272 kB | 0.15 | 17 |
| teqp | 4 | 1.33 s / 63 kB | 6.57 s / 636 kB | 0.35 | 38 |

### Findings

1. **For all-at-once to 4th order, one-shot is the fastest method: about 17× a plain evaluation.**
   That is 9× faster than teqp's per-derivative scheme and 2× faster than polarization. teqp's
   cost explodes because the mixed terms go through nested duals, where `HigherOrderDual<4>`
   carries 16 components per pass and needs 6 separate passes.
2. **Polarization is the compile-time winner.** Every model instantiates exactly one number type,
   `Real<N>`. It has no expression templates and its body is small, so it costs 0.11 s and 13 kB
   per model at N=4, a third of teqp's cost. It runs at about 2× the one-shot cost, and it reuses a
   type teqp already has.
3. **One-shot costs compile time and code size at high order.** The product is C(N+4,4) = 70
   FMAs per multiply at N=4, and it must be fully unrolled to be fast. Without
   `#pragma clang loop unroll(full)`, clang keeps an indexed, latency-bound loop from N=3 upward,
   which is about 4–10× slower. The unrolling is what grows the code. It also makes the result
   depend on the compiler: MSVC needs its own unroll strategy or hand-generated kernels.
4. **The order is the real knob.** The cost roughly doubles per order: 167 → 358 → 789 ns for the
   one-shot on propane. At N=2 every approach is cheap to build (≤0.15 s, ≤17 kB per model), and
   the one-shot runs at 3.6× a plain evaluation. Most properties (p, s, h, cv, cp, w, dp/dT…) need
   only i+j ≤ 2. Third order is needed for things like the fundamental derivative and d²p/dρ², and
   4th order mostly for critical-point conditions.
5. **All-at-once vs on-demand.** If a flash iteration only needs p and ∂p/∂ρ, filling the N=4
   triangle costs 7–11× more than the targeted `Dual2`/`Dual3` request (789 vs ~110 ns).
   Filling to N=2 costs about 1.5×. A tiered cache sidesteps the choice: fill to N=2 by default
   and upgrade to N=4 lazily when a 3rd/4th-order quantity is requested.
