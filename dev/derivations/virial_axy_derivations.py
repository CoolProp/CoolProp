#!/usr/bin/env python3
"""
Zero-density-safe Axy derivations for the CoolProp Helmholtz terms.

Beads: CoolProp-6gf (prep for CoolProp-izv / GH #2991).  Run unbuffered:
    python3 -u dev/derivations/virial_axy_derivations.py

For each term we derive:
  * Axy = tau^x * delta^y * d^{x+y}(a)/(dtau^x ddelta^y), 1/delta-free -> finite at delta=0.
    Feeds get_Ar/get_Aig and p = rho*R*T*(1 + Ar(0,1)) (so p(rho=0)=0 exactly).
  * Virial Taylor coefficients c1 = da/ddelta|_0 (=B*rho_r), c2 = (1/2)d2a/ddelta2|_0
    (=> C = 2*c2/rho_r^2), and tau-derivatives (dB/dT, dC/dT).

Symbolic+ccode for terms being rewritten in source (ideal lead, GenExp, GaoB).
Numeric finiteness+virials for NonAnalytic/SAFT (fractional/sqrt powers blow up
symbolic simplify; their source derivatives are already 1/delta-free, so we only
need to confirm finiteness and the virial values).  Follows src/Helmholtz.cpp:643.

NOTE on sympy: limit(delta**d * ..., delta, 0) with *symbolic* positive d wrongly
returns oo, so finiteness is checked with concrete d (d=1 is the critical lowest power).
"""
import sympy as sy

tau, delta = sy.symbols('tau delta', positive=True)
AXY = [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]


def log(m): print(m, flush=True)


def axy_symbolic(name, alpha, concrete=None):
    """Full simplified Axy + ccode; finiteness via concrete substitution dict `concrete`."""
    log('=' * 78); log(name); log('=' * 78)
    log(f'  alpha = {alpha}')
    for (x, y) in AXY:
        log(f'  ... A{x}{y}')
        a = sy.simplify(tau**x * delta**y * sy.diff(alpha, tau, x, delta, y))
        log(f'  A{x}{y} = {a}')
        log(f'        ccode: {sy.ccode(a)}')
        if concrete:
            checks = []
            for cs in concrete:
                lim = sy.limit(a.subs(cs), delta, 0)
                checks.append(f'{cs}->{sy.nsimplify(lim) if lim.is_number else lim} (finite={lim.is_finite})')
            log('        finiteness@delta=0: ' + ' ; '.join(checks))


def virials_symbolic(name, alpha, d_sym, l_sym=None):
    log('-' * 78); log(f'{name}: virial coefficients'); log('-' * 78)
    da, dda = sy.diff(alpha, delta, 1), sy.diff(alpha, delta, 2)
    combos = []
    for dv in (1, 2, 3):
        if l_sym is None:
            combos.append(({d_sym: dv}, f'd={dv}'))
        else:
            for lv in (1, 2):
                combos.append(({d_sym: dv, l_sym: lv}, f'd={dv},l={lv}'))
    for s, label in combos:
        c1 = sy.simplify(sy.limit(da.subs(s), delta, 0))
        c2 = sy.simplify(sy.limit(dda.subs(s), delta, 0) / 2)
        log(f'  [{label}] c1={c1}')
        log(f'  [{label}] c2={c2}')
        if c1 != 0:
            log(f'  [{label}] dc1/dtau={sy.simplify(sy.diff(c1, tau))}')
        if c2 != 0:
            log(f'  [{label}] dc2/dtau={sy.simplify(sy.diff(c2, tau))}')


def numeric_check(name, alpha, subs):
    """Numeric finiteness of Axy at delta->0 and virial c1,c2 for concrete params."""
    log('=' * 78); log(name + '  (numeric, concrete params)'); log('=' * 78)
    log(f'  params: {subs}')
    for (x, y) in AXY:
        a = (tau**x * delta**y * sy.diff(alpha, tau, x, delta, y)).subs(subs)
        vals = [complex(a.subs(delta, dd)) for dd in (1e-4, 1e-7, 1e-10)]
        finite = all(abs(v) < 1e30 for v in vals)
        log(f'  A{x}{y}: vals@(1e-4,1e-7,1e-10)={[round(v.real,8) for v in vals]} finite={finite}')
    da = sy.diff(alpha, delta, 1).subs(subs)
    dda = sy.diff(alpha, delta, 2).subs(subs)
    c1 = complex(da.subs(delta, 1e-9)).real
    c2 = complex(dda.subs(delta, 1e-9)).real / 2
    log(f'  c1(~B*rho_r) = {c1:.10g}   c2 = {c2:.10g}')


# 1) IDEAL LEAD
a1, a2 = sy.symbols('a1 a2')
axy_symbolic('IDEAL LEAD', sy.log(delta) + a1 + a2 * tau)
log('  NOTE: A00 = log(delta)+a1+a2*tau; A00 not read by p/virials.\n')

# 2) GENEXP element: n*tau^t*delta^d*exp(u)
n, t, d, c, l = sy.symbols('n t d c l', positive=True)
omega, m = sy.symbols('omega m', positive=True)
eta1, eps1, eta2, eps2 = sy.symbols('eta1 epsilon1 eta2 epsilon2')
beta1, gam1, beta2, gam2 = sy.symbols('beta1 gamma1 beta2 gamma2')
u = (-c * delta**l - omega * tau**m - eta1 * (delta - eps1) - eta2 * (delta - eps2)**2
     - beta1 * (tau - gam1) - beta2 * (tau - gam2)**2)
alpha_genexp = n * tau**t * delta**d * sy.exp(u)
axy_symbolic('GENEXP element', alpha_genexp, concrete=[{d: 1, l: 1}, {d: 2, l: 1}])
virials_symbolic('GENEXP element', alpha_genexp, d_sym=d, l_sym=l)

# 3) GAOB element
ng, tg, dg = sy.symbols('n t d', positive=True)
bg, betag, gamg, etag, epsg = sy.symbols('b beta gamma eta epsilon')
alpha_gaob = (ng * tau**tg * sy.exp(1 / (bg + betag * (tau - gamg)**2))
              * delta**dg * sy.exp(etag * (delta - epsg)**2))
axy_symbolic('GAOB element', alpha_gaob, concrete=[{dg: 1}, {dg: 2}])
virials_symbolic('GAOB element', alpha_gaob, d_sym=dg)

# 4) NONANALYTIC element: delta*n*DELTA^b*PSI  (numeric)
nn, a_, b_, beta_ = sy.symbols('n a b beta', positive=True)
A_, B_, C_, D_ = sy.symbols('A B C D', positive=True)
theta = (1 - tau) + A_ * (((delta - 1)**2)**(1 / (2 * beta_)))
PSI = sy.exp(-C_ * (delta - 1)**2 - D_ * (tau - 1)**2)
DELTA = theta**2 + B_ * (((delta - 1)**2)**a_)
alpha_na = delta * nn * DELTA**b_ * PSI
numeric_check('NONANALYTIC element', alpha_na,
              {nn: 0.7, a_: 3.5, b_: 0.85, beta_: 0.3, A_: 0.32, B_: 0.2, C_: 28.0, D_: 700.0, tau: 1.1})

# 5) SAFT association (numeric)
ms, as_ = sy.symbols('m a', positive=True)
vbarn, epsbar, kappabar = sy.symbols('vbarn epsilonbar kappabar', positive=True)
eta_s = vbarn * delta
g_s = sy.Rational(1, 2) * (2 - eta_s) / (1 - eta_s)**3
Deltabar = g_s * (sy.exp(epsbar * tau) - 1) * kappabar
X = 2 / (sy.sqrt(1 + 4 * Deltabar * delta) + 1)
alpha_saft = ms * as_ * (sy.log(X) - X / 2 + sy.Rational(1, 2))
log('=' * 78); log('SAFT association'); log('=' * 78)
log(f'  X(delta->0) = {sy.limit(X, delta, 0)}  (expect 1)')
numeric_check('SAFT association', alpha_saft,
              {ms: 1.0, as_: 1.0, vbarn: 0.45, epsbar: 5.0, kappabar: 0.002, tau: 1.2})

log('\nDONE.')
