import csv, collections, sys
CATS=['ok','loose 1e-4..1e-3','wrong density >1e-3','missed split','false split','error','no reference']
def score(phase_code, rho_code, err, truth):
    if err: return 'error'
    if truth is None: return 'no reference'
    tph, tr = truth
    if tph=='2ph' and not phase_code: return 'missed split'
    if tph=='1ph' and phase_code: return 'false split'
    if tr is None: return 'no reference'
    d=abs(rho_code/tr-1)
    return 'ok' if d<=1e-4 else 'loose 1e-4..1e-3' if d<=1e-3 else 'wrong density >1e-3'
allv={}
for t,name in (('99','CO2/H2O 99/1 (Gernert)'),('50','CO2/H2O 50/50 (Gernert)')):
    T={}
    for l in open(f'gw{t}_truth.txt'):
        f=l.split()
        if f[1]=='NOREF': continue
        rho=float(f[2]); T[f[0]]=(f[1], rho if rho>0 else None)
    K1={r['i']:r for r in csv.DictReader(open(f'gw{t}_k1.csv'))}; K0={r['i']:r for r in csv.DictReader(open(f'gw{t}_k0.csv'))}
    res={'kernel on':collections.Counter(),'kernel off':collections.Counter(),'TPFLSH':collections.Counter()}
    for i in K1:
        tr=T.get(i)
        for lab,R in (('kernel on',K1),('kernel off',K0)):
            r=R[i]; c=score(0<float(r['Q_cp'])<1, float(r['rho_cp']), r['cp_fail']=='1', tr); res[lab][c]+=1
            allv[(f"{name}, {lab}", i)]=c
        r=K1[i]; q=float(r['q_rp']); c=score(0<=q<=1, float(r['rho_rp']), int(r['ierr_rp'])>0, tr); res['TPFLSH'][c]+=1
        allv[(f"{name}, TPFLSH", i)]=c
    print(f"\n{name}")
    print('%-12s'%'', ''.join('%22s'%c for c in CATS))
    for lab,c in res.items(): print('%-12s'%lab, ''.join('%22d'%c[k] for k in CATS))
with open('co2w_scores.csv','w',newline='') as o:
    w=csv.writer(o); w.writerow(['mixture','i','verdict','detail'])
    for (m,i),c in allv.items(): w.writerow([m,i,c,''])
