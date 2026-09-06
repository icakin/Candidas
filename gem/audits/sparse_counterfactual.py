#!/usr/bin/env python3
"""Sparse bottleneck counterfactual: can a modest Tm/Topt shift in ONE or a FEW
essential/bottleneck enzymes create the 40 C boundary (mu40 < detection) while keeping
permissive 34 C growth above floor? If yes with small shifts, the global 32 C claim
does not generalize. If no, the falsification is genuinely strong.

Approach (per relative):
  1. Solve at 40 C, rank costed reactions by pool contribution (coef40 * |flux|).
  2. Candidate set = top-K pool contributors carrying flux.
  3. Rank each candidate by its marginal effect: mu40 when ONLY it is shifted by -Dmax.
  4. Add candidates best-first (cumulative shift -Dmax each) until mu40 < thr; record N.
  5. Check permissive: mu34 with all selected shifted vs floor.
Repeat for Dmax in {3,5,8} C, param in {Tm,Topt}, for each relative.
"""
import cobra, pandas as pd, numpy as np, json, statistics
from pathlib import Path
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
SP={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
    'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
    'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
    'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
RELS=['haemulonii','duobushaemulonii','parapsilosis']
calib=json.load(open(TABLES / 'etcgem_calib.json')); SIG,W,P0,SCALE=calib['sig'],calib['w'],calib['P'],calib['SCALE']
topt=pd.read_csv(TABLES / 'thermal_topt.csv'); tm=pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt,tm):
    df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]
def therm_maps(sp):
    tov={r.g:r.pred_topt for r in topt[topt.sp==sp].itertuples()}
    tmv={r.g:r.pred_tm for r in tm[tm.sp==sp].itertuples()}
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
    mw={r['gene']:float(r['MW_kDa']) for _,r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc=pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv'); kmap={r['reaction']:(float(r['kcat_per_s']),r['best_gene']) for _,r in kc.iterrows()}
    return m,EX,mw,kmap
def precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm):
    terms=[]  # (rid, fv, rv, base_coef, Topt, Tm, gene)
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        kcat_s,bg=kmap.get(r.id,(13.7,None))
        base_coef=statistics.mean(mws)/(kcat_s*3600.0)
        terms.append((r.id,r.forward_variable,r.reverse_variable,base_coef,tov.get(bg,med_to),tmv.get(bg,med_tm),bg))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return terms,con
def act(T,Topt,Tm): return np.exp(-((T-Topt)**2)/(2.0*SIG*SIG))*(1.0/(1.0+np.exp((T-Tm)/W)))
def mu(St,T,dTm=None,dTopt=None):
    dTm=dTm or {}; dTopt=dTopt or {}
    coeffs={}
    for rid,fv,rv,bc,Topt,Tm in St['terms_min']:
        a=max(act(T,Topt+dTopt.get(rid,0.0),Tm+dTm.get(rid,0.0)),1e-6); c=bc/a; coeffs[fv]=c; coeffs[rv]=c
    St['con'].ub=P0; St['con'].set_linear_coefficients(coeffs)
    try: g=St['m'].slim_optimize()
    except Exception: g=float('nan')
    return (g if g==g else 0.0)*SCALE

STATE={}
for sp in ['duobushaemulonii','haemulonii','parapsilosis']:
    tov,tmv,mto,mtm=therm_maps(sp); m,EX,mw,kmap=setup_pool(sp)
    terms,con=precompute(m,EX,mw,kmap,tov,tmv,mto,mtm)
    terms_min=[(t[0],t[1],t[2],t[3],t[4],t[5]) for t in terms]
    STATE[sp]=dict(m=m,con=con,terms=terms,terms_min=terms_min)
    print(f"setup {sp}: {len(terms)} costed reactions",flush=True)

THR=0.05; PERM_FRAC=0.7
out={'thr':THR,'perm_frac':PERM_FRAC,'relatives':{}}
for rel in RELS:
    St=STATE[rel]
    base40=mu(St,40); base34=mu(St,34); floor=PERM_FRAC*base34
    # rank pool contributors at 40C
    sol=St['m'].optimize()
    contrib=[]
    for rid,fv,rv,bc,Topt,Tm in St['terms_min']:
        v=abs(sol.fluxes.get(rid,0.0));
        if v<=0: continue
        a=max(act(40,Topt,Tm),1e-6); contrib.append((rid,(bc/a)*v,Tm,Topt))
    contrib.sort(key=lambda x:-x[1]); cand=[c[0] for c in contrib[:60]]
    relres={'base40':round(base40,3),'base34':round(base34,3),'floor':round(floor,3),
            'n_flux_carrying':len(contrib),'sparse':{}}
    for param in ('Tm','Topt'):
        relres['sparse'][param]={}
        for Dmax in (3.0,5.0,8.0):
            # rank candidates by marginal single-enzyme effect on mu40
            marg=[]
            for rid in cand:
                shift={rid:-Dmax}
                m40=mu(St,40, dTm=shift if param=='Tm' else None, dTopt=shift if param=='Topt' else None)
                marg.append((rid,m40))
            marg.sort(key=lambda x:x[1])  # lowest mu40 first = strongest
            # add best-first cumulatively
            sel={}; achieved=False; nN=None
            for i,(rid,_) in enumerate(marg,1):
                sel[rid]=-Dmax
                m40=mu(St,40, dTm=sel if param=='Tm' else None, dTopt=sel if param=='Topt' else None)
                if m40<THR:
                    achieved=True; nN=i; break
            m34=mu(St,34, dTm=sel if param=='Tm' else None, dTopt=sel if param=='Topt' else None) if sel else base34
            relres['sparse'][param][f'Dmax={Dmax:.0f}']=dict(
                achieved_boundary=bool(achieved), n_enzymes=nN if achieved else f">{len(marg)}",
                mu34_after=round(m34,3), permissive_preserved=bool(m34>=floor),
                single_best_mu40=round(marg[0][1],3))
    # progression curve for the figure: cumulative mu40 as best-first enzymes are shifted -5C
    prog=[]
    for param in ('Tm','Topt'):
        Dmax=5.0
        marg=[]
        for rid in cand:
            shift={rid:-Dmax}
            m40=mu(St,40, dTm=shift if param=='Tm' else None, dTopt=shift if param=='Topt' else None)
            marg.append((rid,m40))
        marg.sort(key=lambda x:x[1])
        sel={}
        for i,(rid,_) in enumerate(marg,1):
            sel[rid]=-Dmax
            m40=mu(St,40, dTm=sel if param=='Tm' else None, dTopt=sel if param=='Topt' else None)
            prog.append(dict(relative=rel,param=param,n=i,mu40=round(m40,4)))
    relres['progression']=prog
    out['relatives'][rel]=relres
    print(f"\n{rel}: base mu40={base40:.3f} mu34={base34:.3f} floor={floor:.3f}")
    for param in ('Tm','Topt'):
        for Dk,v in relres['sparse'][param].items():
            print(f"  {param} {Dk}: boundary={v['achieved_boundary']} n_enzymes={v['n_enzymes']} "
                  f"| single-best mu40={v['single_best_mu40']} | mu34_after={v['mu34_after']} perm_ok={v['permissive_preserved']}")
json.dump(out,open(TABLES / 'sparse_counterfactual_results.json','w'),indent=2)
prog_all=[p for rel in RELS for p in out['relatives'][rel]['progression']]
pd.DataFrame(prog_all).to_csv(G/'sparse_progress.csv',index=False)
# exact detection counts at 44 C from the isolate matrix
M=pd.read_csv(G/'fig3_rate_matrix.csv')
cnt44={}
for sp in ['auris','haemulonii','duobushaemulonii','parapsilosis']:
    b=M[M.sp==sp]; cnt44[sp]=(int(b['p44'].notna().sum()),len(b))
json.dump(cnt44,open(TABLES / 'counts_at_44.json','w'),indent=2)
print("detected at 44 C (n/N):",cnt44)
print("\nwrote sparse_counterfactual_results.json, sparse_progress.csv, counts_at_44.json")
