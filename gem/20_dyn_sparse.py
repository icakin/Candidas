#!/usr/bin/env python3
"""Strengthened sparse counterfactual: DYNAMIC best-first greedy (marginals recomputed
after every addition, so bottlenecks exposed by rerouting are caught) + BEAM search over
pairs and triples among the strongest candidates. Requires BOTH constraints at the end:
mu(40) < detection AND mu(34) >= permissive floor. Also documents where flux reroutes.
Run:  python3 20_dyn_sparse.py Tm   /   python3 20_dyn_sparse.py Topt
"""
import os, sys, cobra, pandas as pd, numpy as np, json, statistics, itertools
from pathlib import Path
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
PARAM=sys.argv[1] if len(sys.argv)>1 else 'Tm'
DMAX=float(sys.argv[2]) if len(sys.argv)>2 else 8.0
G = GEM
SP={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
    'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
    'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
    'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
RELS=['haemulonii','duobushaemulonii','parapsilosis']
calib=json.load(open(TABLES / 'etcgem_calib.json')); SIG,W,P0,SCALE=calib['sig'],calib['w'],calib['P'],calib['SCALE']
topt=pd.read_csv(TABLES / 'thermal_topt.csv'); tm=pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt,tm): df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]
def therm_maps(sp):
    tov={r.g:r.pred_topt for r in topt[topt.sp==sp].itertuples()}; tmv={r.g:r.pred_tm for r in tm[tm.sp==sp].itertuples()}
    return tov,tmv,(statistics.median(tov.values()) if tov else 35.0),(statistics.median(tmv.values()) if tmv else 54.0)
def setup_pool(sp):
    xml,medf=SP[sp]; m=cobra.io.read_sbml_model(str(MODELS / xml)); med=pd.read_csv(INPUTS / medf)
    EX={r.id for r in m.reactions if r.id.startswith(("EX_","Drain")) or r.boundary}
    for r in m.reactions:
        if r.id in EX: r.lower_bound=0.0
    for _,row in med.iterrows():
        rid=str(row['exchange_id'])
        if rid in EX:
            if row['setting']=='OPEN': m.reactions.get_by_id(rid).lower_bound=-1000
            elif row['setting']=='FIT': m.reactions.get_by_id(rid).lower_bound=-10
    for _,row in med[med.setting=='AA_POOL'].iterrows():
        if str(row['exchange_id']) in EX: m.reactions.get_by_id(str(row['exchange_id'])).lower_bound=-3
    # ATP maintenance irreversible and forced in every species. iDC1003 (C. parapsilosis)
    # ships it REVERSIBLE at lb=-3.9 and the solver runs it backwards, making ATP from
    # ADP + Pi. Must match 19_etcgem_counterfactual.py or the two disagree about baselines.
    for _rid in ('ATP_maintenance_NGAM__cyto','ATP_Maintenance__cyto'):
        if _rid in m.reactions:
            _r=m.reactions.get_by_id(_rid); _v=abs(_r.upper_bound); _r.bounds=(_v,_v)
    mw={r['gene']:float(r['MW_kDa']) for _,r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc=pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv'); kmap={r['reaction']:(float(r['kcat_per_s']),r['best_gene']) for _,r in kc.iterrows()}
    return m,EX,mw,kmap
def precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm):
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        kcat_s,bg=kmap.get(r.id,(13.7,None))
        terms.append((r.id,r.forward_variable,r.reverse_variable,statistics.mean(mws)/(kcat_s*3600.0),tov.get(bg,med_to),tmv.get(bg,med_tm)))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return terms,con
def act(T,Topt,Tm): return np.exp(-((T-Topt)**2)/(2.0*SIG*SIG))*(1.0/(1.0+np.exp((T-Tm)/W)))
def mu(St,T,sel):
    coeffs={}
    for rid,fv,rv,bc,Topt,Tm in St['terms']:
        dt=sel.get(rid,0.0)
        a=max(act(T,Topt+(dt if PARAM=='Topt' else 0.0),Tm+(dt if PARAM=='Tm' else 0.0)),1e-6)
        c=bc/a; coeffs[fv]=c; coeffs[rv]=c
    St['con'].ub=P0; St['con'].set_linear_coefficients(coeffs)
    try: g=St['m'].slim_optimize()
    except Exception: g=float('nan')
    return (g if g==g else 0.0)*SCALE
THR=0.05; PERM_FRAC=0.7; POOL=80; MAXSTEPS=60; BEAM=12
_only=os.environ.get('ONLY','').strip()
_dest=TABLES / f'dyn_sparse_{PARAM}.json'
if _only and _dest.exists():
    out=json.load(open(_dest)); RELS=[r for r in RELS if r==_only]
    print(f'merging: recomputing only {_only}, keeping '
          f'{[r for r in out["relatives"] if r!=_only]}',flush=True)
else:
    out={'param':PARAM,'Dmax':DMAX,'thr':THR,'relatives':{}}
for rel in RELS:
    m,EX,mw,kmap=setup_pool(rel); tov,tmv,mto,mtm=therm_maps(rel)
    terms,con=precompute(m,EX,mw,kmap,tov,tmv,mto,mtm)
    St=dict(m=m,con=con,terms=terms)
    base40=mu(St,40,{}); base34=mu(St,34,{}); floor=PERM_FRAC*base34
    sol=m.optimize()
    contrib=[]
    for rid,fv,rv,bc,Topt,Tm in terms:
        v=abs(sol.fluxes.get(rid,0.0))
        if v<=0: continue
        contrib.append((rid,(bc/max(act(40,Topt,Tm),1e-6))*v))
    contrib.sort(key=lambda x:-x[1]); pool=[c[0] for c in contrib[:POOL]]
    # DYNAMIC greedy
    sel={}; traj=[]; achieved=False
    for step in range(MAXSTEPS):
        best=None
        for rid in pool:
            if rid in sel: continue
            trial=dict(sel); trial[rid]=-DMAX
            g40=mu(St,40,trial)
            if best is None or g40<best[1]: best=(rid,g40)
        if best is None: break
        rid,g40=best; sel[rid]=-DMAX
        g34=mu(St,34,sel); traj.append(dict(n=len(sel),rid=rid,mu40=round(g40,4),mu34=round(g34,4)))
        if len(sel)>=2 and g40>=traj[-2]['mu40']-1e-4: pass  # no strong improvement, but keep going a bit
        if g40<THR:
            achieved=(g34>=floor); break
    greedy=dict(achieved_boundary=bool(mu(St,40,sel)<THR),permissive_preserved=bool(mu(St,34,sel)>=floor),
                n_enzymes=len(sel),final_mu40=round(mu(St,40,sel),4),final_mu34=round(mu(St,34,sel),4),
                trajectory=traj[:60])
    # BEAM pairs/triples among top-BEAM by single marginal
    singles=[]
    for rid in pool[:max(BEAM,15)]:
        singles.append((rid,mu(St,40,{rid:-DMAX})))
    singles.sort(key=lambda x:x[1]); top=[s[0] for s in singles[:BEAM]]
    beam_best={'pair':None,'triple':None}
    bp=None
    for a,b in itertools.combinations(top,2):
        g=mu(St,40,{a:-DMAX,b:-DMAX})
        if bp is None or g<bp[1]: bp=((a,b),g)
    if bp: beam_best['pair']=dict(rids=list(bp[0]),mu40=round(bp[1],4),reaches_detection=bool(bp[1]<THR))
    bt=None
    for a,b,c in itertools.combinations(top[:10],3):
        g=mu(St,40,{a:-DMAX,b:-DMAX,c:-DMAX})
        if bt is None or g<bt[1]: bt=((a,b,c),g)
    if bt: beam_best['triple']=dict(rids=list(bt[0]),mu40=round(bt[1],4),reaches_detection=bool(bt[1]<THR))
    # rerouting: shift single strongest enzyme, compare flux
    topr=singles[0][0]
    coeffs={}
    for rid,fv,rv,bc,Topt,Tm in terms:
        dt=-DMAX if rid==topr else 0.0
        a=max(act(40,Topt+(dt if PARAM=='Topt' else 0.0),Tm+(dt if PARAM=='Tm' else 0.0)),1e-6); coeffs[fv]=bc/a; coeffs[rv]=bc/a
    con.ub=P0; con.set_linear_coefficients(coeffs); sol2=m.optimize()
    dflux=[]
    for rid,*_ in terms:
        d=sol2.fluxes.get(rid,0)-sol.fluxes.get(rid,0)
        if abs(d)>1e-6: dflux.append((rid,round(sol.fluxes.get(rid,0),3),round(sol2.fluxes.get(rid,0),3),round(d,3)))
    dflux.sort(key=lambda x:-abs(x[3])); reroute=[dict(rid=r,before=b,after=a,delta=d) for r,b,a,d in dflux[:8]]
    out['relatives'][rel]=dict(base40=round(base40,3),base34=round(base34,3),floor=round(floor,3),
        n_flux_carrying=len(contrib),greedy=greedy,beam=beam_best,top_enzyme=topr,reroute_on_top_shift=reroute)
    json.dump(out,open(TABLES / f'dyn_sparse_{PARAM}.json','w'),indent=2)
    print(f"{rel}: base40={base40:.3f} | greedy n={greedy['n_enzymes']} final_mu40={greedy['final_mu40']} "
          f"reached_det={greedy['achieved_boundary']} perm_ok={greedy['permissive_preserved']} | "
          f"beam pair mu40={beam_best['pair']['mu40']} triple mu40={beam_best['triple']['mu40']}",flush=True)
print(f"wrote dyn_sparse_{PARAM}.json")
