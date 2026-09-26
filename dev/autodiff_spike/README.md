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

---

# Experiment 3: exploiting structure — hand derivatives vs structure-aware AD

This asks whether the ~17× penalty is intrinsic to AD or comes from ignoring the model's
structure. The experiment has two halves:

- **Multiparameter (CoolProp, hand-coded).** `coolprop_bench.cpp` times
  `ResidualHelmholtzGeneralizedExponential::all()`, which fills the whole i+j ≤ 4 triangle, and
  `all_deltaonly()`. It compares them against a value-only loop over the same term arrays, a naive
  `Taylor2<4>` pass over every term, and a **separable AD** version. The separable version runs a
  single-variable `Taylor1<4>` in δ and in τ per term, then forms the outer product n·F_j(δ)·G_i(τ).
  The bench checks that the GenExp block is the fluid's entire residual.
- **PC-SAFT, structure-aware AD** (`methods_all/a_structured.cpp`). The model is written in
  (η, T) with η = ρ·q(T). At fixed composition, αr = Σ cₖ(T)·φₖ(η) + Σᵢ ln gᵢᵢ(η; T). All the η
  arithmetic (I₁, I₂, C₁, the hard-sphere rationals) runs in single-variable Taylor arithmetic,
  and so does all the T arithmetic. The only bivariate work is composing onto
  η(u, v) = ρ₀(1+v)·q(u): shared powers of H = η − η₀, one axpy per φₖ, one
  single-variable × bivariate product per coefficient, and one bivariate log per component.

Build the multiparameter bench against a Release static CoolProp:

```bash
cmake -B build_rel -S ../.. -DCMAKE_BUILD_TYPE=Release -DCOOLPROP_STATIC_LIBRARY=ON && cmake --build build_rel -j8
# then compile coolprop_bench.cpp with the CXX_INCLUDES from
# build_rel/CMakeFiles/CoolProp.dir/flags.make, plus -I. and build_rel/libCoolProp.a
```

## Multiparameter: ns per call, whole N=4 triangle (15 values)

| fluid, state | terms | value only | **hand `all()`** | hand δ-only | naive AD | separable AD |
|---|---|---|---|---|---|---|
| propane, 300 K, 11 000 mol/m³ | 18 | 111 | **315** (2.8×) | 209 (1.9×) | 1670 (15×) | 576 (5.2×) |
| propane, 300 K, 100 mol/m³ | 18 | 105 | **317** (3.0×) | 209 | 1628 (16×) | 564 (5.4×) |
| nitrogen, 100 K, 25 000 mol/m³ | 36 | 235 | **678** (2.9×) | 471 | 4037 (17×) | 1324 (5.6×) |
| nitrogen, 300 K, 400 mol/m³ | 36 | 246 | **707** (2.9×) | 495 | 3996 (16×) | 1350 (5.5×) |
| R1234yf, 300 K, 10 000 mol/m³ | 17 | 106 | **298** (2.8×) | 196 | 1450 (14×) | 551 (5.2×) |
| n-decane, 400 K, 4500 mol/m³ | 12 | 80 | **208** (2.6×) | 137 | 1075 (13×) | 315 (3.9×) |

Agreement with the hand values is ≤2e-14 at the dense states. At the dilute states, the naive and
separable AD differ from the hand code by up to 4e-11 and 5e-12 respectively, concentrated in the
small high-order δ-terms. There is no independent reference here, so which side is right is
unresolved.

## PC-SAFT: ns per call, whole triangle (same run and machine state as the table below)

| | value only | N=2 | N=4 | N=4 / value |
|---|---|---|---|---|
| propane (liquid) | 55 | | | |
| naive one-shot (`taylor2`) | | 200 | 1025 | 19× |
| teqp scheme | | 545 | 8671 | 158× |
| **structure-aware** | | **117** | **320** | **5.8×** |
| C1/C2/C3 mixture | 114 | | | |
| naive one-shot | | 323 | 1982 | 17× |
| **structure-aware** | | **188** | **518** | **4.5×** |

Compile cost at N=4 is 0.16 s and 12 kB per model, about the same as polarization and a third
of the naive one-shot or teqp scheme.

**Accuracy: this corrects Experiments 1–2.** The structure-aware form stays at ≤9e-15 through
4th order *at the dilute gas state*, where every generic formulation lost 5–6 digits (1e-9).
The loss attributed above to "evaluating αr in double" is really cancellation in the generic
formulation's arithmetic: the 1/ζ₀ division and ζ₂³/ζ₃² − ζ₀ in a_hs, where the terms scale like
1/ρ and cancel. It is a property of that computational graph. Writing the model in η removes it.

## Summary: cost of the full 4th-order triangle relative to one value evaluation

| model | hand | structure-aware AD | naive one-shot AD |
|---|---|---|---|
| multiparameter (GenExp) | 2.6–3.0× | 3.9–5.6× | 13–17× |
| PC-SAFT | (not hand-coded) | 4.5–5.8× | 17–19× |

- **Hand derivatives are cheap, but not quite "a couple of flops per order".** The recursion adds
  about 70 flops per term against 1–3 transcendentals in the value. Measured, the full triangle is
  about 2.8× a value evaluation, and the δ-only path is about 1.9×.
- **Most of the naive AD penalty is structure, not AD.** Generic bivariate AD costs about 15× on
  the multiparameter form too. With structure respected, AD lands within 1.5–2× of hand code on
  the multiparameter form. PC-SAFT then costs about the same per value evaluation as CoolProp's
  hand-coded multiparameter derivatives.
- **The remaining multiparameter gap to hand code is the exponentials.** Separable AD pays for
  2–4 single-variable Taylor `exp`s per term, while the hand code needs one scalar `exp` plus
  recursions on the argument.

---

# Experiment 4: Chebyshev density rootfinding for PC-SAFT (`pcsaft_cheb.py`)

This asks whether a Bell & Alpert (2018) style all-roots density solver can work for PC-SAFT,
including mixtures.

**Reformulation.** With η = ζ₃ = ρ·q(T), fixed (T, x) makes the pressure equation univariate
in η:

  F(η) := ηZ = η + η²·dαr/dη = p·q(T)/(RT)

The target only shifts the constant term. η ∈ [0, 0.74] is a natural bounded interval: close
packing is π/(3√2) ≈ 0.7405, and the η = 1 pole lies outside it. So the Chebyshev domain is
*fixed and bounded for every fluid, mixture and temperature*, with no compactification. T and x
enter only through scalar weights (the same decomposition as Experiment 3). That decomposition
reproduces ηZ from the direct model to ≤7e-12.

**Singularities.** The distance of the nearest complex singularity from the interval sets the
degree per piece.
- η = 1: a pole and a log branch point.
- Zeros of g_hs: η ≥ 1.5 or ≤ −1 for a ∈ [1, 3]. Note b_i = 2a_i²/9 identically, so the
  chain term is a one-parameter family in a_i = 3·D_i·r₂.
- **Poles of C1:** 1/C1 = P(η, m̄) / ((1−η)⁴(2−η)²), with P a degree-6 polynomial that is
  *affine in m̄*. For long chains P has a real zero just left of η = 0: η ≈ −0.036 at m̄ = 10 and
  −0.011 at m̄ = 30. That forces 11–14 pieces at n=12.
  **Fix:** root-find G = P²·(F − target) instead. That is (P²F) − target·P², so the target still
  enters linearly. The C1 poles are gone and the piece count at n=12 no longer depends on m̄.
- **Gas root precision:** fit Q = P²Z, which is O(1), and apply the factor η exactly with a
  Chebyshev multiply-by-x. The absolute error near η = 0 then scales with η.

**Results.** The cases are propane, C1/C2/C3, 70/30 methane/n-decane and n-decane, at
T = 40–800 K and p = 1e-4 to 1000 MPa (288 cases).
- **Root counts:** the Chebyshev root count matched a 200 001-point scan with bisection in every
  case.
- **Spurious roots:** it finds Privat-type roots, with four roots at low T. The fourth, at
  η ≈ 0.62–0.72, is past the physical liquid root. Those must be rejected by Gibbs energy, not by
  mechanical stability.
- **Accuracy:** the median relative root error is 1.4e-11. With the singularity-based 5 pieces
  the worst is 4e-6, in n-decane at 40 K, where Z spans orders of magnitude inside the first
  piece. The singularity bound ignores that amplitude variation. **About 10 pieces at n=12
  give ≤2e-9 everywhere tested (≤1e-11 at 400 K), and 20 pieces give ≤1e-10.**

**Towards universal tables (not built here).** After multiplying by P²:
- every term of P²F is a scalar weight w_k(T, x) times a *fluid-independent* function of η;
- P² is quadratic in m̄, and I₁ and I₂ are affine in (c₁ₘ, c₂ₘ);
- the one exception is the per-component chain term, a smooth one-parameter family in a_i,
  i.e. a 2-D table on the (η, a) rectangle.

So the Bell–Alpert "stacked per-term representation" carries over: universal Chebyshev tables
on fixed rectangles, with runtime assembly as weighted sums of coefficient vectors. Assembly
cost is small next to the colleague-matrix eigensolves on the 1–3 pieces that can hold a root.
Neither cost has been measured here.

**Caveats.** Association (the site-fraction solve) and polar terms (Padé forms with their own
poles, like C1) are not included. Association is smooth in η at fixed (T, x), so a per-(T, x) fit
from node evaluations works, though it isn't universal. The polar terms' Padé poles need the
same singularity analysis.

### Where the spurious roots live (`spurious_scan.py`)

Brute-force root counts over T = 20–400 K and p = 1e-4 to 1e3 MPa.

| fluid | GS2001: highest T with >3 roots | extra root η | triple point | Liang 2012 / 2014 constants |
|---|---|---|---|---|
| methane | 28 K | 0.71–0.74 | 90.7 K | none |
| propane | 84 K | 0.63–0.74 | 85.5 K | none |
| n-decane | 140 K | 0.56–0.73 | 243.5 K | none |
| C1/C2/C3 | 54 K | 0.66–0.74 | — | none |
| C1/nC10 70/30 | 88 K | 0.63–0.74 | — | none |

- **Gross & Sadowski constants:** the extra roots occur only at or below the pure-fluid triple
  point. For the mixtures, they appear only at temperatures far below where a liquid mixture
  exists.
- **Liang constants:** these remove the extra roots entirely over this range. Caveat: here they
  are paired with the GS2001 pure-component parameters, not the refitted parameters they were
  published with, so this is qualitative only.

---

# Experiment 5: speed of the universal-table Chebyshev solver (`chebsolve.cpp`)

Build it with
`c++ -std=c++20 -O3 -DNDEBUG -mcpu=native -I<eigen> chebsolve.cpp`.

## What is computed when

| tier | when | what | cost |
|---|---|---|---|
| universal | offline, once for all fluids | 8 fixed η-pieces on [0, 0.74]; degree-12 Chebyshev coefficients of 27 basis functions (×η, applied exactly); Π₀, Π₁, Π₂ of P² = Π₀ + m̄Π₁ + m̄²Π₂; three 2-D (η, a) chain tables, degree 10 in a ∈ [0.75, 2.5] | 55 kB; 0.4 ms to build |
| components | once per fluid set | m_i, σ_i, ε_i, k_ij. Nothing else: no table is component-specific | — |
| state (T, x) | every call | N `exp`s for d_i(T); s₀..s₃, r_n, m̄, c₁, c₂, E₁, E₂, A₁..A₃, B₁, B₂, a_i, w_i, q (about 30 scalars); then G = Σ W_j·U_j plus the chain contraction, per piece | 0.5–0.8 µs |
| pressure | every call | t = p·q/(RT); subtract t·P²; certified subdivision on the pieces not proven root-free; ρ = η/q | 0.7–0.9 µs |

P depends on composition only through m̄, and only affinely. So "rebuilding P and Q" at a new
(T, x) is a linear combination of stored vectors with about 30 scalar weights. It is never a refit.
q depends on (T, x), not on ρ, so t is known before rootfinding starts.

## Results

Cases: 7 temperatures (150–500 K) × 7 pressures (1 kPa – 100 MPa), 147 states, 227 roots. The
reference is a 200 000-point scan of the explicit Z with 80-step bisection.

| case | assemble | certified roots | total | max rel. root error | counts | colleague eigensolve | grid + Illinois | direct node fit (instead of assemble) |
|---|---|---|---|---|---|---|---|---|
| propane | 0.55 µs | 0.90 µs | **1.5 µs** | 1.5e-11 | all ok | 38 µs, 7e-5 | 1.5 µs, 5e-5 | 17 µs |
| C1/nC10 70/30 | 0.52 | 0.91 | **1.4** | 2.2e-11 | all ok | 36 µs, 3e-7 | 1.5 µs, 1e-6 | 17 |
| C1/C2/C3 | 0.79 | 0.68 | **1.5** | 2.9e-11 | all ok | 33 µs, 4e-5 | 1.3 µs, 4e-6 | 19 |

- **Assembly accuracy:** G agrees with η·P²·Z from the model to ≤1e-9, scaled by
  ηP²(1+|Z|). The worst point is at the close-packing end, η = 0.74.
- **Certified subdivision.** On each piece: exclude if |c₀| > Σ|c_k|. If the derivative series
  passes the same test, the piece is monotone, so bracket by endpoint signs and refine with
  Illinois. Otherwise re-expand onto the two halves and recurse. The half-interval re-expansions
  are two precomputed 14×14 matrices. This beats the eigensolve by 40× and is more accurate on
  the tiny gas roots, because the recursion zooms in. A root can only be missed if it is a
  tangency (double root) finer than the depth cap of 12.
- **Direct fitting** at the nodes per (T, x) costs about 17 µs. That is why the universal tables
  matter.

**Not covered:**
- **Temperature range:** T < 150 K is not timed. Deeply subcooled states need more pieces (Experiment 4).
- **Table rectangle:** states with a_i outside [0.75, 2.5] are rejected (no fallback
  implemented); none occurred here.
- **Warm cache:** timings are warm, with the 55 kB of tables resident.
- **Association and polar terms** are not included.

---

# Experiment 6: GERG-2008 multi-fluid Chebyshev density solver (`gerg_cheb.cpp`)

Build it against the Release static library. `-force_load` is needed so that the GERG backend's
static registration survives linking:

```bash
c++ -std=c++20 -O3 -DNDEBUG -mcpu=native -DNQ_DEG=20 <CoolProp CXX_INCLUDES> -I. gerg_cheb.cpp \
    -Wl,-force_load,build_rel/libCoolProp.a -o gerg_cheb
./gerg_cheb [delta_max=4] [tol=1e-12] [fast]
```

`gen_cheb2bern.py` emits `cheb2bern.hpp`: exact Chebyshev-to-Bernstein matrices, computed in
rational arithmetic and then rounded.

## Structure

Every GERG-2008 residual term, pure-fluid and departure, is n·τᵗ·e^{u_τ(τ)} · δᵈ·e^{u_δ(δ)}. This
includes water and CO₂, because the GERG backend builds its own pure-fluid EOS from power terms
only; the GERG departure β(δ−γ) terms are stored as δ-linear. With ρ = δ·ρ_r(x):

  G(δ) = δ·Z(δ) − p/(ρ_r·R·T),  Z = 1 + Σ_g W_g(T, x)·χ_g(δ)

where χ_g = δ·d/dδ[δᵈe^{u_δ}] and W_g = Σ over the group's terms of X_k·κ_k(τ), with X_k equal to
x_i or x_i·x_j·F_ij.

| tier | what | cost |
|---|---|---|
| component set | collect terms; **group terms with identical δ-parameters** (91→32, 164→38 and 92→30 δ-functions for C1/C2/C3, 5-component natural gas and humid air); adaptive shared δ-pieces on [0, 4]; degree-20 Chebyshev tables | 5–12 ms; 90–260 kB |
| per (T, x) | τ, ρ_r from CoolProp's reducing function; one `exp` per term for κ; G per piece as Σ_g W_g·C_g | 2–8 µs |
| per p | subtract t; per piece: coefficient exclusion, then Bernstein conversion (exact matrix, rigorous roundoff bound), Descartes + de Casteljau subdivision, safeguarded Newton on the Chebyshev series | 2–6 µs |

## What had to be fixed along the way

- **Piece acceptance.** A test relative to each term's own size never converges for the
  exp(−δ⁶) terms, which are about 1e-24 by δ ≈ 2 (this produced 1 100–1 900 pieces). The error
  has to be budgeted in units of Z, weighted by max|κ(τ)| over the τ range of use. Weights reach
  n·τ³⁰ at low T.
- **Budget relative to the smallest local magnitude, not the piece maximum.** Otherwise the gas
  root at δ ~ 1e-4 (1 kPa, 100 K) came out at 5e-6. There is no roundoff floor on the piece
  touching δ = 0.
- **Rootfinding.** Chebyshev-coefficient exclusion and monotonicity tests are loose on strongly
  varying pieces. Bernstein/Descartes is tighter. The bigger fix was the refinement stopping test:
  a 1e-16 absolute bracket width never triggered, which cost 50–77 iterations per root. Safeguarded
  Newton with a ulp-relative stop takes 10–14.
- **Memory layout.** Piece-major tables let the assembly stream contiguously.

## Results

Degree 20, table tolerance 1e-12. The states are 9 temperatures × 8 pressures (1 kPa – 100 MPa)
per mixture. The reference is a 100 000-point scan of the direct term sum on δ ∈ (0, 6] with
bisection. The term sum reproduces CoolProp's GERG2008 pressure to ≤2e-15.

| mixture | δ-groups | pieces | roots | max rel. error | counts | roots with δ > 4 | assemble | roots | total | CoolProp `solver_rho_Tp` (one root) |
|---|---|---|---|---|---|---|---|---|---|---|
| C1/C2 50/50 | 25 | 32 | 184 | 1.5e-13 | all ok | 0 | 3.9 µs | 5.3 µs | 9.2 µs | 18 µs |
| C1/C2/C3 50/30/20 | 32 | 38 | 180 | 4.5e-14 | all ok | 0 | 5.7 | 5.5 | 11.2 | 32 |
| natural gas (C1, N₂, CO₂, C2, C3) | 38 | 41 | 158 | 3.0e-14 | all ok | 0 | 7.6 | 4.9 | 12.5 | 54 |
| C1/H₂S 50/50 (type III) | 20 | 28 | 166 | 4.1e-13 | all ok | 0 | 2.1 | 4.5 | 6.6 | 8 |
| humid air (N₂, O₂, Ar, CO₂, H₂O) | 30 | 37 | 72 | 3.4e-14 | all ok | 0 | 5.3 | 2.2 | 7.5 | 17 |

**Caveat: timings were taken on a heavily loaded machine.** The load average was 12–34 from
other sessions' test suites; CoolProp's own solver timed 39 µs, then 106 µs, for the same
natural-gas calls. Ratios against `solver_rho_Tp` in the same run (0.23–0.83) are more
meaningful than the absolute numbers. Note that `solver_rho_Tp` returns one root from a guess;
this solver returns all roots in [0, δ_max].

Low-temperature states have up to 5 real roots of the GERG equation, all of them found.

**Not covered:**
- **Range:** δ_max = 4 is fixed, and the scan found no roots in (4, 6].
- **Non-analytic terms:** none are handled. Reference-EOS water and CO₂ (IAPWS-95, Span–Wagner)
  would need a direct add-in.
- **Tangencies:** the depth cap of 16 means a root pair closer than 2⁻¹⁶ of a piece width
  (tangency at a spinodal) can be missed.
- **Wider sampling:** this is a grid of states, not random sampling.

### Loose tables plus a Newton polish on the true equation

`POLISH=1 ./gerg_cheb 4 <tol>` runs one Newton step or more on the *true* grouped term sum,
starting from each Chebyshev root. Each step reuses the assembly's group weights, so it costs one
pow and one exp per δ-group. A step that moves the root more than 1e-3 relative is rejected, and
the Chebyshev root is kept.

| tables | pieces | Newton steps / root | rejected | max rel. error after polish | root counts |
|---|---|---|---|---|---|
| degree 12, 1e-6 | 33–49 | 1.36–1.49 | 0 of 760 | ≤3.3e-15 | all ok |
| degree 12, 1e-4 | 22–42 | 1.38–1.65 | 0 | ≤3.3e-15 | all ok |
| degree 12, 1e-3 | 18–38 | 1.41–1.69 | 0 | ≤3.3e-15 | all ok |
| degree 16, 1e-6 | 22–38 | 1.31–1.43 | 0 | ≤3.5e-15 | all ok |
| degree 16, 1e-4 | 14–33 | 1.36–1.57 | 0 | ≤3.5e-15 | all ok |
| degree 16, 1e-3 | 13–33 | 1.36–1.61 | 0 | ≤3.5e-15 | all ok |

Without the polish, 1e-6 tables already give roots within ≤3.3e-8 (degree 16), because the
tolerance is in units of Z and applies at every point. Loosening from 1e-12 to 1e-6 cuts the
pieces from 43–60 to 22–38. Below about 1e-4 the count stops falling: it is then set by the
absolute budget on the piece touching δ = 0 and by the roundoff floor.

**Risk.** Loosening affects the *root count*, not Newton convergence. A fit error ε in Z can merge
or split a near-tangent root pair, i.e. at pressures within about ε (relative) of a spinodal
pressure. Near a spinodal the Newton basin shrinks (radius ~ |G′/G″|), so a start within 1 ppm can
converge onto the partner root. The 1e-3 guard does not catch that. The fix is a polish that is
bracketed: the subdivision's bracket ends have certified signs for the *true* G wherever
|G_cheb| > ε there.
