#!/usr/bin/env python3
"""common_network.py - COMMON-NETWORK CONTROL. Run each species' enzyme profile
(Topt/Tm/kcat/MW) on the IDENTICAL auris reaction scaffold, mapping auris genes to
each species' orthologs. This removes any network/reconstruction difference, so any
remaining TPC difference is purely the enzyme-SEQUENCE effect. If optima still
compress here, network similarity was not what drove the compression.
"""
import cobra, json, statistics, numpy as np, pandas as pd, warnings
import build_etcgem_tpc as B
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
warnings.filterwarnings('ignore')
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
grid=np.arange(22,44.1,1.0)
# auris scaffold + medium (reuse setup)
def load_species_data(sp):
    to=pd.read_csv(f'{TABLES}/thermal_topt.csv'); to=to[to.id.str.startswith(sp+'|')]
    to={r.id.split('|')[1]:r.pred_topt for r in to.itertuples()}
    tm=pd.read_csv(f'{TABLES}/thermal_tm.csv'); tm=tm[tm.id.str.startswith(sp+'|')]
    tm={r.id.split('|')[1]:r.pred_tm for r in tm.itertuples()}
    mw={r['gene']:float(r['MW_kDa']) for _,r in pd.read_csv(f'{TABLES}/enzyme_mw_{sp}.csv').iterrows()}
    kc={r['reaction']:float(r['kcat_per_s']) for _,r in pd.read_csv(f'{TABLES}/kcat_reaction_{sp}.csv').iterrows()}
    return to,tm,mw,kc
# auris gene -> species ortholog (RBH)
def ortho_map(sp):
    if sp=='auris': return None
    f={'haemulonii':'pid_haemulonii','duobushaemulonii':'pid_duobushaemulonii','parapsilosis':'pid_parapsilosis'}[sp]
    d=pd.read_csv(f'{TABLES}/rbh/{f}.tsv',sep='\t',names=['a','b','pid'])   # written by 07_transfer_orthologs.py --rbh-out
    return {r.a:r.b for r in d.itertuples()}

# build auris scaffold once
m0,EX,mw_au,kmap_au=B.setup_pool('auris')
auto,autm,aumw,aukc=load_species_data('auris')

print(f"{'species(on auris net)':<24}{'opt':>6}{'peak_mu':>9}")
opts=[]
for sp in SPP:
    to,tm,mw,kc=load_species_data(sp); omap=ortho_map(sp)
    m=m0  # same scaffold object; rebuild constraint each species
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        base=kmap_au.get(r.id); aug = base[1] if base else None   # auris best gene
        # map to species ortholog
        sg = aug if sp=='auris' else (omap.get(aug) if omap else None)
        Topt = to.get(sg, statistics.median(list(to.values())))
        Tm   = tm.get(sg, statistics.median(list(tm.values())))
        MWv  = mw.get(sg, statistics.median(list(mw.values())))
        ks   = kc.get(r.id, base[0] if base else 13.7)
        terms.append((r.forward_variable,r.reverse_variable,MWv/(ks*3600.0),Topt,Tm))
    con=m.problem.Constraint(0.0,ub=P0,name='pp'); m.add_cons_vars(con); m.solver.update()
    best=(-1,None)
    for t in grid:
        co={}
        for fv,rv,bc,To,Tm in terms:
            a=max(B.act(t,To,Tm,SIG,W),1e-3); c=bc/a; co[fv]=c; co[rv]=c
        con.set_linear_coefficients(co); g=(m.slim_optimize() or 0)
        if g>best[0]: best=(g,t)
    opts.append(best[1]); print(f"{sp:<24}{best[1]:>5.0f}C{best[0]*SC:>9.3f}")
    m.remove_cons_vars([con]); m.solver.update()
print(f"\noptima spread on COMMON auris network = {max(opts)-min(opts):.0f} C")
print("If still ~1-2C, compression is a SEQUENCE effect, not a network-reconstruction artifact.")
