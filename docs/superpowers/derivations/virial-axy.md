# Zero-Density-Safe `Axy` Derivatives and Exact Virial Coefficients

**Beads:** CoolProp-6gf (prep) → CoolProp-izv (impl) · **GitHub:** #2991 (child of #1676)
**Source script:** `dev/derivations/virial_axy_derivations.py` (sympy; reproduce with `python3 -u dev/derivations/virial_axy_derivations.py`)

This is the front-loaded derivation that turns the implementation (CoolProp-izv) into transcription. Every closed form below was produced symbolically by sympy and the finiteness/virial claims are machine-checked; the core identity is additionally cross-checked numerically against CoolProp's own derivative engine (see §7).

---

## 1. Notation and the two things we need

Reduced Helmholtz `a = a0(τ,δ) + ar(τ,δ)`. Define the dimensionless group

```
Axy = tau^x * delta^y * d^{x+y}(a) / (dtau^x ddelta^y)
```

`get_Ar(x,y)` returns the residual `Axy`; `get_Aig(x,y)` the ideal-gas one (teqp-style naming, user-selected).

**Why `Axy` and not the raw derivative:** the raw `d(ar)/d(delta)` is coded in CoolProp via a log-derivative form multiplied by `1/delta` (`Helmholtz.cpp:267`, `:164`), GaoB has an explicit `/delta` (`:735`), and the ideal lead is `1/delta` (`:1148`). All are `Inf`/`NaN` at exactly `delta=0`. The `delta^y`-scaled `Axy` cancels the `1/delta` analytically, so it is finite at `delta=0`.

**(a) Pressure.** `p = rho*R*T*(1 + Ar(0,1))`. At `rho=0`, `delta=0`, `Ar(0,1)=0` (finite) ⇒ `p=0` *exactly*. The current `rho*R*T*(1 + delta*dalphar_dDelta)` evaluates `0*Inf = NaN` (verified: `update(DmolarT,0,T)` throws "p is not a valid number" for every HEOS fluid).

**(b) Virial coefficients** are the `delta`-Taylor coefficients of `ar` at `delta=0` (since `ar = c1*delta + c2*delta^2 + ...`):

```
c1 = d(ar)/d(delta)|_0           B(T)  = c1 / rho_r
c2 = (1/2) d^2(ar)/d(delta^2)|_0  C(T)  = 2*c2 / rho_r^2
dB/dT = (dc1/dtau) * dtau/dT / rho_r          dtau/dT = -T_r / T^2
dC/dT = 2*(dc2/dtau) * dtau/dT / rho_r^2
```

`c1`/`c2` get contributions **only from low `delta`-power terms** (`d=1` → c1; `d∈{1,2}` → c2). This is the "separate path for virial coefficients at zero density" #2991 asks for; `get_Ar`/`get_Aig` is the "general (non-virial) derivative" path.

**Key principle for every term:** accumulate `Axy` directly in the `delta^y`-scaled form so each `1/delta^k` cancels against an explicit `delta^k` *in source* — never compute `delta * (something/delta)` at runtime.

> sympy caveat: `limit(delta**d * ..., delta, 0)` with *symbolic* positive `d` wrongly returns `oo`. Finiteness is therefore checked with concrete `d` (`d=1` is the critical lowest power). The product forms themselves are correct for symbolic `d`.

---

## 2. Ideal-gas lead — `a0_lead = log(delta) + a1 + a2*tau`

The canonical zero-density case (`d(a0)/d(delta) = 1/delta` is `Inf`, but `Aig(0,1)=1`):

| Axy | value | ccode | finite@0 |
|-----|-------|-------|----------|
| `Aig(0,1)` | `1`        | `1`       | ✓ (`=1`) |
| `Aig(0,2)` | `-1`       | `-1`      | ✓ |
| `Aig(1,0)` | `a2*tau`   | `a2*tau`  | ✓ |
| `Aig(1,1)` | `0`        | `0`       | ✓ |
| `Aig(1,2)` | `0`        | `0`       | ✓ |

**All other ideal-gas terms are functions of `tau` only** → their `delta`-derivatives are exactly 0; they contribute only to `A[itau][0] = tau^itau * d^itau(a0)/dtau^itau` (mechanical `tau`-scaling of existing `tau`-derivatives, already finite). `A00` contains `log(delta)` (non-finite at 0) and is **never read** by pressure/virials.

---

## 3. Generalized exponential (workhorse — most fluids)

```
ar_elem = n*tau^t*delta^d*exp(u),
u = -c*delta^l - omega*tau^m - eta1*(delta-eps1) - eta2*(delta-eps2)^2
    - beta1*(tau-gam1) - beta2*(tau-gam2)^2
```

Let `Bdelta = d - c*l*delta^l - delta*(eta1 + 2*eta2*(delta-eps2))`  (= `delta * dln(ar_elem)/ddelta`, the `1/delta`-free bracket)
and `Btau   = t - omega*m*tau^m - tau*(beta1 + 2*beta2*(tau-gam2))`   (note sympy printed `-tau*(beta1 - 2*beta2*(gamma2-tau))`).

**Axy product forms** (all verified finite=0 at `delta=0` for `d=1,2`):

```c
// A01 = delta*d(ar)/ddelta
A01 = pow(delta,d)*n*pow(tau,t)*( -c*pow(delta,l)*l + d - delta*(eta1 + 2*eta2*(delta-epsilon2)) )*exp(u);

// A02 = delta^2 * d2(ar)/ddelta2
A02 = -pow(delta,d)*n*pow(tau,t)*( c*pow(delta,l)*pow(l,2) - c*pow(delta,l)*l - d*(d-1)
        + 2*d*(c*pow(delta,l)*l + delta*(eta1+2*eta2*(delta-epsilon2))) + 2*pow(delta,2)*eta2
        - pow(c*pow(delta,l)*l + delta*(eta1+2*eta2*(delta-epsilon2)), 2) )*exp(u);

// A10 = tau*d(ar)/dtau
A10 = pow(delta,d)*n*pow(tau,t)*( -m*omega*pow(tau,m) + t - tau*(beta1 - 2*beta2*(gamma2-tau)) )*exp(u);
```

(`A11`, `A12` likewise — full ccode in the script output; both finite at `delta=0`.)

**Virial coefficients** (with `u0 := u(tau, delta=0) = beta1*gam1 - beta1*tau - beta2*gam2^2 + 2*beta2*gam2*tau - beta2*tau^2 + eps1*eta1 - eps2^2*eta2 - omega*tau^m`):

| case | `c1` (→ B) | `c2` (→ C) |
|------|-----------|-----------|
| `d=1, l=1` | `n*tau^t*exp(u0)` | `n*tau^t*(-c + 2*eps2*eta2 - eta1)*exp(u0)` |
| `d=1, l=2` | `n*tau^t*exp(u0)` | `n*tau^t*(2*eps2*eta2 - eta1)*exp(u0)` |
| `d=2` (any l) | `0` | `n*tau^t*exp(u0)` |
| `d>=3` | `0` | `0` |

Note the `-c` cross-term enters `c2` **only when `l=1`** (then `delta^{d+l}=delta^2`); for `l>=2` it does not. `dc1/dtau`, `dc2/dtau` are in the script output (used for dB/dT, dC/dT).

**XiangDeiters** = `phi0 + acentric*phi1 + theta*phi2`, each `phi_k` a GenExp sub-term → apply the above per sub-term and combine linearly. No new derivation.

---

## 4. Non-analytic (IAPWS-95-style; Water, CO2, …)

```
ar_elem = delta*n*DELTA^b*PSI
theta = (1-tau) + A*((delta-1)^2)^(1/(2*beta));  PSI = exp(-C*(delta-1)^2 - D*(tau-1)^2)
DELTA = theta^2 + B*((delta-1)^2)^a
```

No `1/delta` in source ⇒ raw derivatives already finite at `delta=0`; the implementation just scales the existing coded contributions by `delta^y * tau^x`. **Numerically verified** (params `n=0.7,a=3.5,b=0.85,beta=0.3,A=0.32,B=0.2,C=28,D=700,tau=1.1`): `A01,A02,A10,A11,A12` all → 0 and finite at `delta→0`.

**Virials:** the explicit `delta` factor gives
```
c1 = (n*DELTA^b*PSI)|_{delta=0} = n * ((1-tau+A)^2 + B)^b * exp(-C - D*(tau-1)^2)
c2 = (1/2) d^2/ddelta^2|_0   (use existing coded d2alphar_ddelta2 contribution at delta=0)
```
For the test params `c1 ≈ 1.35e-16` — exponentially tiny because `PSI(0)=exp(-C)`, i.e. non-analytic terms contribute negligibly to B (physically expected; they are localized near the critical point). Finite and correct.

---

## 5. Gao-B

```
ar_elem = n*tau^t*exp(1/(b+beta*(tau-gam)^2)) * delta^d*exp(eta*(delta-eps)^2)
```

The source `/delta` (`Helmholtz.cpp:735`) is a coding choice; functionally `ar_elem` is analytic at `delta=0`. Axy product forms (finite at 0 for `d=1,2`):

```c
// A01
A01 = n*pow(tau,t)*( d*pow(delta,d) + 2*pow(delta,d+1)*eta*(delta-epsilon) )
      * exp( (eta*(b+beta*pow(gamma-tau,2))*pow(delta-epsilon,2) + 1)/(b+beta*pow(gamma-tau,2)) );

// A02
A02 = pow(delta,d)*n*pow(tau,t)*( 4*d*delta*eta*(delta-epsilon) + d*(d-1)
        + 2*pow(delta,2)*eta*(2*eta*pow(delta-epsilon,2)+1) )
      * exp( eta*pow(delta-epsilon,2) + 1.0/(b+beta*pow(gamma-tau,2)) );
```

(`A10`,`A11`,`A12` in script output; all finite at 0.)

**Virials** (with `w := eta*eps^2 + 1/(b+beta*(tau-gam)^2)`):

| case | `c1` | `c2` |
|------|------|------|
| `d=1` | `n*tau^t*exp(w)` | `-2*eps*eta*n*tau^t*exp(w)` |
| `d=2` | `0` | `n*tau^t*exp(w)` |
| `d>=3` | `0` | `0` |

`dc1/dtau`, `dc2/dtau` in script output.

---

## 6. SAFT association

```
ar = m*a*(log(X) - X/2 + 1/2),  X = 2/(sqrt(1+4*Deltabar*delta)+1)
Deltabar = g(vbarn*delta)*(exp(epsbar*tau)-1)*kappabar,  g(eta)=0.5*(2-eta)/(1-eta)^3
```

`X(delta=0) = 1` (verified) ⇒ `ar(delta=0)=0`. Denominators are powers of `(2*Deltabar*X*delta+1) → 1`, so all coded `X`-derivatives are finite at `delta=0`; scale by `delta^y*tau^x` for `Axy`. **Numerically verified** (`m=a=1, vbarn=0.45, epsbar=5, kappabar=0.002, tau=1.2`): `A01,A02,A10,A11,A12` all finite, → 0 as `delta→0`; `c1=-0.402`, `c2=-0.129` (finite, nonzero — SAFT does contribute to B, C). Implementation reuses the existing closed-form `X`-derivatives (`Helmholtz.cpp:861-1110`).

---

## 7. Generalized cubic

Delegates to `m_abstractcubic->alphar(tau,delta,z,itau,idelta)` (raw mixed partial). Cubic EOS are analytic at `rho=0` (`alphar ~ ln`-of-linear-in-`delta`), so `Axy = tau^x*delta^y * alphar(...,x,y)` is finite at `delta=0`. No `1/delta` to cancel; just apply the `delta^y*tau^x` scaling.

---

## 8. Numeric cross-check of the foundation

`B*rho_r = c1 = lim_{delta→0} d(ar)/d(delta)`, checked against CoolProp's own `dalphar_dDelta` (independent of the internal `1e-12` hack), at `T=300 K`:

| fluid | `B_ref` (m³/mol) | `c1=B*rho_r` | `dar_dDelta@δ=1e-9` |
|-------|------------------|--------------|---------------------|
| Nitrogen | −4.5537e−6 | −5.09280e−2 | −5.09280e−2 |
| CO2 | −1.2127e−4 | −1.288523 | −1.288523 |
| Methane | −4.2210e−5 | −0.4279731 | −0.4279731 |
| Water | −1.2013e−3 | −21.47170 | −21.47170 |

Convergence to ~6 sig figs confirms the Taylor-coefficient identity the virial path relies on. (Water differs at δ=1e-3 — the non-analytic terms — but converges by δ=1e-6.)

---

## 9. Implementation mapping (for CoolProp-izv)

| Term | `all()` in `src/Helmholtz.cpp` | Work |
|------|-------------------------------|------|
| Ideal lead | `:1143` | add `A[0][1]=1, A[0][2]=-1, A[1][0]=a2*tau`; rest 0 |
| Other ideal | `:1154–1430` | `A[itau][0] = tau^itau * (existing dtau-deriv)`; `delta`-orders 0 |
| GenExp | `:138` | accumulate §3 `A01..A12` (no `one_over_delta`); virial `c1/c2` per `d,l` table |
| NonAnalytic | `:386` | scale existing raw derivs by `delta^y*tau^x`; `c1=§4` |
| GaoB | `:665` | accumulate §5 forms (drop `/delta`); virial table |
| SAFT | `:1077` | scale existing `X`-deriv-based contributions by `delta^y*tau^x` |
| cubic | `:615` | `A[x][y]=tau^x*delta^y*abstractcubic->alphar(...,x,y)` |
| XiangDeiters | `:800` | linear combine GenExp sub-terms |

`HelmholtzDerivatives` stores `A[itau][idelta]` (`itau∈{0,1}, idelta∈{0,1,2}`) + `getA(itau,idelta)`; backend `calc_Ar`/`calc_Aig`; public `get_Ar`/`get_Aig`; `calc_pressure = rho*R*T*(1+Ar(0,1))` with `rho==0` short-circuit; virials from `c1/c2`. Per-term gate: hand-derived `Axy` matches `delta^y*tau^x*(raw deriv)` to ~1e-9 away from 0 (fast-path oracle) and `one_mcx` autodiff (absolute oracle).
