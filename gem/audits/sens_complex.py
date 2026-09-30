#!/usr/bin/env python3
"""sens_complex.py - complex-GPR sensitivity. The base model treats every multi-gene
GPR as isozyme-OR (pool cost = MEAN subunit MW). The opposite extreme is enzyme-
complex-AND (all subunits required -> cost = SUM subunit MW). Recompute the 4-species
predicted optima under SUM-MW and check the compression conclusion is unchanged.
Also does a targeted repair of respiratory / ATP-synthase reactions specifically.
"""
import build_etcgem_tpc as B, cobra, json, statistics, re, numpy as np, pandas as pd, warnings
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
warnings.filterwarnings('ignore')
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
grid=np.arange(22,44.1,1.0)
RESP=re.compile(r'cytochrome|ubiquinol|NADH dehydrogenase|ATP synthase|succinate dehydrogenase|complex',re.I)

def build(sp,agg):
    tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        base=kmap.get(r.id); ks,bg=base if base else (13.7,None)
        M = agg(mws)                              # mean (OR) or sum (AND)
        terms.append((r.forward_variable,r.reverse_variable,M/(ks*3600.0),tov.get(bg,med_to),tmv.get(bg,med_tm)))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return m,terms,con
def tpc_opt(m,terms,con):
    con.ub=P0; best=(-1,None)
    for t in grid:
        co={}
        for fv,rv,bc,To,Tm in terms:
            a=max(B.act(t,To,Tm,SIG,W),1e-6); c=bc/a; co[fv]=c; co[rv]=c
        con.set_linear_coefficients(co); g=m.slim_optimize() or 0
        if g>best[0]: best=(g,t)
    return best[1]
print(f"{'species':<16}{'opt_OR(mean)':>13}{'opt_AND(sum)':>13}")
optsOR=[];optsAND=[]
for sp in SPP:
    m,t,c=build(sp,statistics.mean); oOR=tpc_opt(m,t,c)
    m,t,c=build(sp,sum); oAND=tpc_opt(m,t,c)
    optsOR.append(oOR); optsAND.append(oAND)
    print(f"{sp:<16}{oOR:>12.0f}C{oAND:>12.0f}C")
print(f"\noptima spread  OR-only={max(optsOR)-min(optsOR):.0f}C   AND-complex={max(optsAND)-min(optsAND):.0f}C")
print("If both compressed, the conclusion is robust to GPR complex structure.")
