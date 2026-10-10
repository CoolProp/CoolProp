# Chebyshev mixture density rootfinder — spec (v2, corrected)

**Issue:** [CoolProp-3gt]
**Status:** spec
**Depends on:** [CoolProp-aor] (closed)

## Architecture (corrected v2)

> The v1 of this spec proposed building 2D Chebyshev expansions of `P(τ, δ; z)`,
> either per-z (cached) or per-call. **Both ideas were wrong.** The correct
> architecture is to expand each individual EOS term in **δ only**, with the
> τ-dependence carried as a closed-form scalar prefactor evaluated at runtime.

### Why δ-only per-term works

A residual Helmholtz EOS is a sum of terms. The dominant family,
`ResidualHelmholtzGeneralizedExponentialElement`
(`include/Helmholtz.h:300`), has the form

```
αᵢʳ(τ, δ) = nᵢ · δ^{dᵢ} · τ^{tᵢ}
          · exp(-cᵢ·δ^{lᵢ}
                - ωᵢ·τ^{mᵢ}
                - η₁,ᵢ(δ-ε₁,ᵢ) - η₂,ᵢ(δ-ε₂,ᵢ)²
                - β₁,ᵢ(τ-γ₁,ᵢ) - β₂,ᵢ(τ-γ₂,ᵢ)²)
```

This factors exactly into

```
αᵢʳ(τ, δ) = κᵢ(τ) · φᵢ(δ)
```

where

```
κᵢ(τ) = nᵢ · τ^{tᵢ}
       · exp(-ωᵢ·τ^{mᵢ} - β₁,ᵢ(τ-γ₁,ᵢ) - β₂,ᵢ(τ-γ₂,ᵢ)²)

φᵢ(δ) = δ^{dᵢ} · exp(-cᵢ·δ^{lᵢ}
                     - η₁,ᵢ(δ-ε₁,ᵢ) - η₂,ᵢ(δ-ε₂,ᵢ)²)
```

`φᵢ(δ)` depends only on the term's **structural** parameters
(`dᵢ, cᵢ, lᵢ, η, ε`) — not on T, ρ, or z. It is therefore a
**static, per-term object**: build one Chebyshev-in-δ expansion per term, once,
at EOS load time, and reuse forever.

`κᵢ(τ)` is a single scalar per term per τ-evaluation, evaluable in O(1) closed form.

### Mixture extension is free

The corresponding-states multifluid model is

```
αʳ_mix(τ, δ; z) = Σᵢ zᵢ · αᵢʳ(τ, δ)
                + Σᵢ Σⱼ zᵢ zⱼ Fᵢⱼ · αᵢⱼ^{dep}(τ, δ)
```

where each `αᵢʳ` and each departure function `αᵢⱼ^{dep}` is itself a sum of
generalized-exponential terms. Therefore at any (τ, z),

```
αʳ_mix(δ; τ_eval, z) = Σ_{terms k} wₖ(τ_eval, z) · φₖ(δ)
```

with weights

```
wₖ(τ_eval, z) = (linear in z·zⱼFᵢⱼ) · κₖ(τ_eval)
```

The right-hand side is a **linear combination of Chebyshev polynomials in δ**,
which is itself a Chebyshev polynomial in δ. So:

> Building the mixture `P(δ)` Chebyshev at any (τ, z) costs **one pass over
> the term list**: evaluate `κₖ(τ_eval)`, scale `φₖ` (its precomputed coef
> vector), accumulate. No `alphar` evaluations. No rootfinding inside the build.

Differentiation in δ is a Chebyshev coefficient operation (linear), so we get
`∂P/∂ρ` and `∂²P/∂ρ²` for free.

### Memory and build budget

- **Memory:** O(N_terms_total) per EOS, independent of mixture composition.
  A typical pure-fluid EOS has ~10–25 terms; binary departure functions add
  ~5–15 each. At n=8 Chebyshev nodes per term per piece (say 4 pieces over
  the δ rectangle), each `φₖ` is ~256 doubles = 2 KB. Per fluid: ~50 KB.
  Per binary pair: ~30 KB. AGA-8 Amarillo (10 fluids, up to 45 binary pairs):
  ~1.8 MB total. Loaded once per process.

- **Build:** one-shot at EOS-load time, in C++ (or precomputed and shipped
  as JSON, like the existing pure-fluid superanc). For each term: sample
  `φₖ(δ)` at Chebyshev nodes, run `dyadic_splitting` until residuals are
  below tolerance.

- **Runtime per density rootfind:**
  - Pre-build mixture `P(δ)` Chebyshev at this (τ, z): O(N_terms · n_coefs) ≈ 2–10 μs.
  - Rootfind: companion-matrix eigenvalues on a degree-n polynomial, or
    bracketed Brent on the polynomial, ≈ 1–5 μs.
  - **Total ≈ 5–20 μs per call**, independent of how many calls came before.

This beats both options in v1 of the spec, and beats the baseline
`solver_rho_Tp` p95 (33–1376 μs across the four systems measured in `aor`).

## Term-type coverage

| Term class                                | Separable? | Plan |
|---|---|---|
| `ResidualHelmholtzGeneralizedExponential` | yes (proof above) | precompute `φₖ(δ)` |
| `ResidualHelmholtzGaoB`                   | likely yes (Gao 2020 form is separable) | verify, precompute |
| `ResidualHelmholtzXiangDeiters`           | composes GenExp internally | recurse |
| `ResidualHelmholtzGeneralizedCubic`       | **no** (cubic α^r mixes τ, δ) | direct eval, add to baseline of mixture-`P(δ)` Chebyshev |
| `ResidualHelmholtzNonAnalytic`            | **no** (Span-Wagner Δ^b·δ·ψ) | direct eval, sample to add a δ-correction Chebyshev at runtime if needed |
| `ResidualHelmholtzSAFTAssociating`        | **no** (Δ̄(τ,δ) cross) | direct eval |

For non-separable terms we have two choices:

- **(A) Direct add-in.** Evaluate the non-separable term's `α^r` and its
  δ-derivatives directly at the Chebyshev nodes of the rectangle (n+1 points),
  build a single δ-Chebyshev for *that term* at the current τ, and add into
  the mixture polynomial. Cost: a handful of `alphar` calls per rootfind for
  the non-separable terms only. Most fluids (incl. AGA-8 components) have
  no non-analytic terms; CO2 has a few.

- **(B) Pre-built 2D fallback** for non-separable terms only — small enough
  in count to be tractable. Defer to a follow-up issue.

Default to **(A)**.

## API sketch

```cpp
// New: include/superancillary/term_density_chebyshev.h
namespace CoolProp::superancillary::density {

/// Static, per-term δ-Chebyshev expansion of φ(δ) = α^r_term(τ, δ) / κ(τ).
/// Built once at EOS load. Owns dyadic-split pieces of degree n.
struct TermPhi {
    std::vector<ChebyshevExpansion<std::vector<double>>> pieces;
    // ... metadata: term index in EOS, structural params used to build, ...
};

/// Per-EOS-term cache populated lazily on first density solve.
class EosTermCache {
public:
    explicit EosTermCache(const ResidualHelmholtzContainer& alphar,
                          double delta_min, double delta_max,
                          int Ndegree = 8);

    /// Number of separable terms (those that have φₖ).
    std::size_t separable_count() const;

    /// Get φₖ for term k (separable terms only).
    const TermPhi& phi(std::size_t k) const;

    /// Evaluate κₖ(τ) for term k (closed-form, ~O(1)).
    double kappa(std::size_t k, double tau) const;

    /// True if term k is non-separable (caller must direct-eval).
    bool is_nonseparable(std::size_t k) const;
};

/// Composes per-fluid + per-binary-departure caches into a mixture rootfinder.
class MixtureDensityChebyshev {
public:
    /// Build the mixture P(δ) Chebyshev at this (τ, z). One linear pass.
    /// Returns coefs in the same dyadic-split structure used by TermPhi.
    std::vector<ChebyshevExpansion<std::vector<double>>>
    build_P_of_delta(double tau, const std::vector<double>& z) const;

    /// Solve for δ such that P(δ) = p_target. Returns roots in [δ_min, δ_max],
    /// sorted ascending. Empty if none in domain.
    std::vector<double>
    solve(double tau, const std::vector<double>& z, double p_target) const;
};

} // namespace
```

## Pseudocode — entry point

```
function solver_rho_Tp_chebyshev(T, p, z):
    tau   = T_reducing(z) / T
    delta_min, delta_max = domain_for_z(z)   // hard rectangle, no extrapolation
    P_cheb = mixture_cheb.build_P_of_delta(tau, z)
        // ≈ O(N_terms · n_coefs), 2–10 μs
        // - linear scaling of precomputed φₖ by κₖ(tau) and (zᵢ, zⱼ Fᵢⱼ)
        // - direct eval for non-separable terms, add into accumulator

    roots = P_cheb.real_roots_in([delta_min, delta_max])
    if len(roots) == 0:
        throw "No density solution in domain"

    if len(roots) == 1:
        return roots[0] · ρ_reducing(z)

    // Multiple roots — vapor-like + liquid-like + (rare) supercritical bumps.
    // Pick by lower Gibbs energy. g(δ) needs μ(δ); reuse same per-term machinery
    // (we already have α^r as a Chebyshev, ∂α^r/∂τ via ∂κₖ/∂τ).
    return argmin_g(roots, tau, z) · ρ_reducing(z)
```

## Decisions locked

- **Domain rectangle in δ:** single common reduced-δ rectangle `[1e-6, 6]`
  (revised from `[1e-10, 6]` after Phase-A empirical work — `δ ≤ 1e-6`
  forces dyadic splits that buy no useful accuracy because the function
  values themselves are already below double-precision noise)
  applied to every per-term `φₖ(δ)` across all fluids. δ here is mixture-
  reduced (δ = ρ / ρ_reducing(z)), which degenerates to native δ_critical for
  pure-fluid mode. Per-EOS validity bounds become advisory warnings, not
  rectangle constraints. Queries outside the rectangle throw `OutOfDomainError`;
  caller (flash, `solver_rho_Tp`) falls back to existing
  `solver_rho_Tp_global`. "Good extrapolation" (asymptotic-aware Chebyshev
  extension) is out of scope; file as a future enhancement if needed.

- **Where `φₖ(δ)` is built:** JSON-preferred, C++ fallback. Same pattern as
  saturation superanc. Production fluids ship a JSON blob (extension of
  `dev/fluids/<fluid>.json`) emitted by an external Python tool analogous
  to `fastchebpure`. Runtime deserializes when present; if absent (new EOS
  in development, or a fluid not yet processed by the offline tool), C++
  builds in-memory at EOS load via `dyadic_splitting` — same engine the
  saturation superanc would use, just sampling `φₖ` at Chebyshev nodes
  rather than `rhoL/rhoV/p` at saturation.

- **Build parameters:** degree n = 12, tolerance 1e-8, max_refine_passes = 20
  (revised from n=8/tol=1e-12 after CoolProp-42i tuning sweep on methane —
  higher degree dominates tighter tol; recommended config gives ~520 KB per
  fluid for a 40-term EOS with ~10⁻¹² accuracy in α^r reconstruction).

- **Performance target (revised):** 50–150 μs per density-solve at 5–10 component
  mixtures with full τ flattening, matching Bell-Alpert §5 Fig. 11 on modern
  hardware. The earlier "5–20 μs" target in this spec was aspirational and not
  paper-grounded. With τ fixed (table walk along an isotherm), the paper
  achieves ~30 μs and that's a reasonable stretch goal for cached calls.

- **Non-separable terms (`NonAnalytic`, `SAFTAssociating`, `GenCubic`): not
  modeled by the Chebyshev rootfinder.** Rationale:
    - `ResidualHelmholtzNonAnalytic` (Span-Wagner critical-region patches)
      is only relevant in the immediate critical neighborhood; the Chebyshev
      rootfinder is explicitly a "good-enough" path elsewhere, and callers
      that need critical-region accuracy fall back to `solver_rho_Tp_global`.
    - `ResidualHelmholtzSAFTAssociating` is acceptable to drop in this scope.
    - `ResidualHelmholtzGeneralizedCubic` lives in cubic-EOS backends, which
      already have explicit/analytic density solutions — they don't need or
      want a Chebyshev rootfinder.
  Implementation: silently skip non-separable terms when summing the mixture
  `P(δ)` polynomial. EOS that contain only separable terms get full Chebyshev
  treatment; EOS with non-separable terms get a Chebyshev that's biased by
  exactly the missing contribution. The bias is small away from the critical
  patches the non-analytic terms were added to fix.

- **Polynomial rootfind: iterative Chebyshev (paper §2.4.2 + §4.4),
  *not* companion matrix.** Per dyadic piece: oversample the residue
  `r̃(δ) = p̃(δ) - p_target` at `2N+1` Chebyshev-Lobatto nodes; if all
  oversampled values are "far from zero" relative to a tolerance, skip the
  piece (rootfinding-avoidance heuristic). Otherwise, locate sign-change
  brackets among the oversampled nodes, fit a local quadratic, polish with
  modified secant (Pegasus/King). Companion-matrix eigenvalue solve is
  reserved for verification — at N=50 it's ~700 μs per piece (paper §2.4.1),
  too slow for the hot path. Existing `superancillary::detail::balance_matrix`
  (James-Langou-Lowery, already cited) is for the verification path.

- **Solution polish on the FULL EOS, not the Chebyshev** (paper §4.4.1).
  After all candidate δ roots are collected from the Chebyshev rootfinder,
  each one gets 1–2 Halley steps using `α^r_δ` and `α^r_δδ` from the
  production EOS (Eq 60–63). Two iterations reach <10⁻¹⁰% pressure error.
  Phase pruning / Gibbs comparison among polished roots also uses full-EOS
  values. So: Chebyshev = reliable enumeration; full EOS = accuracy + phase
  selection. This sidesteps the bias from skipping non-analytic terms in
  the Chebyshev — those terms are still in the polish step and the Gibbs
  comparison.

- **Near-critical fallback: geometric, in (τ, δ) space.** The reducing point
  (τ=1, δ=1) is the natural origin of the critical neighborhood, regardless
  of mixture composition (the reducing function defines τ, δ). Two-phase
  check:
    - **Pre-rootfind:** if `|τ - 1| > τ_far` (e.g. 0.20), the query is far
      enough from critical that the Chebyshev path is faithful; proceed.
    - **Post-rootfind:** after solving for δ, if `|τ - 1| < τ_near` (e.g.
      0.10) AND `|δ - 1| < δ_near` (e.g. 1.0), the query landed inside the
      critical neighborhood where dropped non-analytic terms bias the answer
      — redo via `solver_rho_Tp_global` and return that.
    - Between τ_near and τ_far: trust the Chebyshev answer.
  Defaults conservative; tunable per-EOS via metadata (some fluids have
  wider non-analytic decay than others). For pure-fluid mode the same
  thresholds apply directly. For mixtures this is a "reducing-point
  neighborhood," not the true mixture critical line — the true critical
  line is more expensive to characterize and not needed here.

## Optimization layers (Q5 follow-up)

Layered, so we ship the simple version first and add as benchmarks justify.

1. **Colleague matrix instead of standard companion.** The natural eigenvalue
   form for a Chebyshev expansion is the *colleague matrix* (Boyd 2002,
   "Computing zeros on a real interval through Chebyshev expansion and
   polynomial rootfinding"; Chebfun uses this). Tridiagonal-plus-rank-1
   structure → cheaper than a dense O(n³) general eigenvalue solve. Existing
   `superancillary::detail::balance_matrix` already operates on a similar
   matrix; reuse with the colleague form.

2. **Interval-bound pruning over the dyadic tree.** Each `φₖ` piece has a
   precomputable `[P_min_piece, P_max_piece]` envelope (sum of |coefs|
   gives a quick upper bound on |P|; tighter bounds via piecewise extrema
   at build time). At runtime, prune any piece whose envelope does not
   straddle `p_target`. In typical queries this eliminates most pieces
   without any eigenvalue work — log(n_pieces) work instead of n_pieces.

3. **Closed-form roots for low-degree pieces (n ≤ 3).** Quadratic or cubic
   formula instead of eigenvalues. Most dyadic splits converge to small n
   in smooth regions of P(δ).

4. **Newton-Raphson polish after eigenvalue extraction.** Eigenvalue roots
   are accurate to ~n·εᵐᵃᶜʰ. One or two Newton steps on the polynomial
   bring them to εᵐᵃᶜʰ. Cost negligible (one polynomial evaluation per
   step); critical when downstream code (e.g. Gibbs comparison between
   coalesced roots) is sensitive.

5. **Single P(δ) build feeds all derivative needs.** Once the mixture
   polynomial `P(δ)` is assembled at (τ, z), `∂P/∂ρ` and `∂²P/∂ρ²` are
   linear coefficient transforms. So the Newton polish in (4), the spinodal
   check (`∂P/∂ρ = 0`), and the second-derivative phase classification
   (`∂²P/∂ρ²`) all come from the same coefficient vector. No additional
   per-call sampling.

6. **Group τ-prefactor evaluations by exponent.** Many terms share the
   same `t` (τ-power) or `m` (exp(-ω·τ^m) exponent). Hoist common
   subexpressions in the κ(τ) sweep — measurable on Amarillo where 10
   fluids' worth of GenExp terms have a lot of structural overlap.

Defer until measured worth: term-set SIMD across fluids, JSON-cached
per-piece interval envelopes, custom small-matrix QR for n=4–8.

## Open questions

All resolved — see "Decisions locked" and "Optimization layers (Q5 follow-up)" above.

[CoolProp-3gt]: ../.beads/  "bd show CoolProp-3gt"
[CoolProp-aor]: ../.beads/  "bd show CoolProp-aor"
