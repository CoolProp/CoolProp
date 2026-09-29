"""Corrected per-state verdicts (v2) for the ratio maps.  Changes vs v1 (tools/verdict.cpp):
 - rp_wrong is re-judged with REFPROP's own FGCTY2/PRESS (no CoolProp code) and split into
   rp_wrong (split off equilibrium by >= 1e-3 in ln f or p, and density differs from CoolProp)
   and rp_minor (density matches CoolProp within 1e-4 - phase label only; or off equilibrium by
   only 1e-5..1e-3; or not reproduced on a re-call from 10-digit-rounded CSV inputs);
 - adds cp_wrong for states BOTH codes call single phase where a REFPROP-only brute-force
   tangent-plane test finds a split (tm < -1e-8) - v1 only judged disagreements;
 - rp_scope kept only where the REFPROP-only stability test confirms instability (see stab_scope_*)."""
import csv, os
pairs=[('v_dense_%d.csv'%k,'dense_%d.csv'%k,'stab_dense_%d.out'%k) for k in (1,2,4,5,6)]+[('v_dense_3_fixed.csv','dense_3_fixed.csv','stab_dense_3_fixed.out')]+[('v_hard_%d_k1.csv'%k,'hard_%d_k1.csv'%k,'stab_hard_%d_k1.out'%k) for k in (1,2,3,4)]
selfj={}
for l in open('rp_self_perstate.tsv'):
    m,i,w=l.rstrip('\n').split('\t'); selfj[(m,i)]=w
scope_ok={}
if os.path.exists('stab_scope.tsv'):
    for l in open('stab_scope.tsv'):
        m,i,tm=l.rstrip('\n').split('\t'); scope_ok[(m,i)]=float(tm)<-1e-8
from collections import Counter
tot=Counter()
for v,d,st in pairs:
    D={(r['mixture'],r['i']):r for r in csv.DictReader(open(d))}
    mix=next(iter(D))[0]
    rows=[]
    for r in csv.DictReader(open(v)):
        key=(r['mixture'],r['i']); verdict=r['verdict']; det=r['detail']
        if verdict=='rp_wrong':
            s=D[key]; rc,rr=float(s['rho_cp']),float(s['rho_rp']); Qc=float(s['Q_cp']); w=selfj[key]
            if w=='not2ph': verdict,det='rp_minor','TPFLSH two-phase not reproduced on re-call from 10-digit CSV inputs (same knife-edge artifact as cp_history)'
            elif not (0<Qc<1) and abs(rr-rc)<=1e-4*abs(rc): verdict,det='rp_minor','TPFLSH density matches the single-phase answer; only the two-phase label is wrong (REFPROP self-check %s)'%w
            elif float(w)<1e-3: verdict,det='rp_minor','TPFLSH split off equilibrium by only %s in its own ln f / p (loose convergence?)'%w
            else: det='TPFLSH split off equilibrium by %s in REFPROP\'s own ln f / p'%w
        elif verdict=='rp_scope' and scope_ok:
            if not scope_ok.get(key,False): verdict,det='unresolved','CoolProp split not confirmed by the REFPROP-only stability test'
        rows.append((r['mixture'],r['i'],verdict,det))
    for l in (open(st) if os.path.exists(st) else []):
        f=l.split()
        if len(f)==4 and f[3]!='NOREF' and float(f[3])<-1e-8:
            rows.append((mix,f[0],'cp_wrong','both codes single phase; REFPROP-only stability test finds a split (tm=%s)'%f[3]))
    with open('v2_'+d,'w',newline='') as o:
        w=csv.writer(o); w.writerow(['mixture','i','verdict','detail']); w.writerows(rows)
    c=Counter(x[2] for x in rows); tot+=c
    print('%-28s'%mix[:28], dict(c))
print('TOTAL', dict(tot))
