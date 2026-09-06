#!/usr/bin/env python3
"""sens_medium_kcat.py - two sensitivity analyses on the 4-species predicted optima:
(1) MEDIUM: restricted vs rich amino-acid availability (AA_POOL exchange bound).
(2) kcat: MEDIAN-over-substrates kcat instead of the MAX ("best") used in the base run.
Report predicted-optimum spread under each; the compression should persist.
"""
import build_etcgem_tpc as B, cobra, json, statistics, numpy as np, pandas as pd, warnings
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
warnings.filterwarnings('ignore')
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
grid=np.arange(22,44.1,1.0)

def median_kcat_map(sp):
    raw=pd.read_csv(f'{TABLES}/kcat_{sp}.tsv',sep='\t')   # raw per-substrate DLKcat output, from 11_aggregate_kcat.py
    raw['gene']=raw['Substrate Name'].str.split('|').str[0]
    gmed=raw.groupby('gene')['Kcat value (1/s)'].median().to_dict()
    kc=pd.read_csv(f'{TABLES}/kcat_reaction_{sp}.csv')
    return {r['reaction']:(gmed.get(r['best_gene'],float(r['kcat_per_s'])),r['best_gene']) for _,r in kc.iterrows()}

def opt(sp,aa_bound=None,kmap_override=None):
    tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
    if aa_bound is not None:
        med=pd.read_csv(f'{INPUTS}/'+({'parapsilosis':'medium_iDC1003_parapsilosis.csv'}.get(sp,'medium_iRV973_auris.csv')))
        for _,row in med[med.setting=='AA_POOL'].iterrows():
            rid=str(row['exchange_id'])
            if rid in EX: m.reactions.get_by_id(rid).lower_bound=aa_bound
    if kmap_override: kmap=kmap_override
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        base=kmap.get(r.id); ks,bg=base if base else (13.7,None)
        terms.append((r.forward_variable,r.reverse_variable,statistics.mean(mws)/(ks*3600.0),tov.get(bg,med_to),tmv.get(bg,med_tm)))
    con=m.problem.Constraint(0.0,ub=P0,name='pp'); m.add_cons_vars(con); m.solver.update()
    best=(-1,None)
    for t in grid:
        co={}
        for fv,rv,bc,To,Tm in terms:
            a=max(B.act(t,To,Tm,SIG,W),1e-3); c=bc/a; co[fv]=c; co[rv]=c
        con.set_linear_coefficients(co); g=(m.slim_optimize() or 0)
        if g>best[0]: best=(g,t)
    return best[1]

print("(1) MEDIUM sensitivity (AA_POOL bound): optima per species")
print(f"{'species':<16}{'restricted(-1)':>15}{'base(-3)':>10}{'rich(-1000)':>12}")
for tag,rows in [('',None)]:
    R=[];Bs=[];Ri=[]
    for sp in SPP:
        r=opt(sp,-1); b=opt(sp,-3); ri=opt(sp,-1000)
        R.append(r);Bs.append(b);Ri.append(ri)
        print(f"{sp:<16}{r:>14.0f}C{b:>9.0f}C{ri:>11.0f}C")
    print(f"  optima spread: restricted={max(R)-min(R):.0f}C  base={max(Bs)-min(Bs):.0f}C  rich={max(Ri)-min(Ri):.0f}C")

print("\n(2) kcat MEDIAN vs MAX: optima per species")
print(f"{'species':<16}{'kcat=MAX(base)':>15}{'kcat=MEDIAN':>13}")
Mx=[];Md=[]
for sp in SPP:
    mx=opt(sp); md=opt(sp,kmap_override=median_kcat_map(sp))
    Mx.append(mx);Md.append(md); print(f"{sp:<16}{mx:>14.0f}C{md:>12.0f}C")
print(f"  optima spread: MAX={max(Mx)-min(Mx):.0f}C  MEDIAN={max(Md)-min(Md):.0f}C")
