#!/usr/bin/env python3
"""18_build_etcgem_tpc.py - assemble sequence-driven growth TPCs from the enzyme- and
temperature-constrained GEMs, calibrate on C. auris, and predict the other species.

Temperature model (transparent, uses the directly-predicted per-enzyme parameters):
  Each reaction r carries a DLKcat kcat_r (its optimal kcat) tied to its best_gene g,
  for which we predicted Topt_g and Tm_g from sequence. Relative activity:

     act_g(T) = peak_g(T) * death_g(T)
     peak_g(T)  = exp( -(T - Topt_g)^2 / (2 sig^2) )     # Gaussian peak at Topt_g
     death_g(T) = 1 / (1 + exp( (T - Tm_g)/w ))          # denaturation cutoff at Tm_g

  Effective kcat_r(T) = kcat_r * act_g(T)  (peak ~1 at the optimum). In the sMOMENT
  pool constraint  sum_r (MW_r / kcat_r(T)) v_r <= P, moving away from each enzyme's
  optimum tightens the pool and lowers growth -> a growth-vs-T curve.

  ONLY {Topt_g, Tm_g} differ between species; the shape hyperparameters {sig, w, P}
  are shared, calibrated ONCE on C. auris, then FROZEN to predict the other three.
  So any between-species TPC difference comes only from the predicted thermal params.
"""
import cobra, pandas as pd, numpy as np, json, statistics
from pathlib import Path
from scipy.optimize import minimize
# The gem/ directory, resolved from this file's own location. This used to be a
# hard-coded absolute path from the machine the script was written on, which meant the
# module raised FileNotFoundError at IMPORT time on any other machine, before any caller
# could override it. Set CANDIDAS_ROOT to point elsewhere.
import os as _os
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
SP={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
    'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
    'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
    'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
MEAS={'auris':'auris','haemulonii':'Hae','duobushaemulonii':'Duo','parapsilosis':'para'}

topt=pd.read_csv(TABLES / 'thermal_topt.csv'); tm=pd.read_csv(TABLES / 'thermal_tm.csv')
for df,col in ((topt,'pred_topt'),(tm,'pred_tm')):
    df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]

def therm_maps(sp):
    tov={r.g:r.pred_topt for r in topt[topt.sp==sp].itertuples()}
    tmv={r.g:r.pred_tm  for r in tm[tm.sp==sp].itertuples()}
    med_to=statistics.median(tov.values()) if tov else 35.0
    med_tm=statistics.median(tmv.values()) if tmv else 54.0
    return tov,tmv,med_to,med_tm

def setup_pool(sp):
    xml,medf=SP[sp]
    m=cobra.io.read_sbml_model(str(MODELS / xml)); med=pd.read_csv(INPUTS / medf)
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
    kc=pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv')
    kmap={r['reaction']:(float(r['kcat_per_s']),r['best_gene']) for _,r in kc.iterrows()}
    return m,EX,mw,kmap

def act(T,Topt,Tm,sig,w):
    peak=np.exp(-((T-Topt)**2)/(2.0*sig*sig))
    death=1.0/(1.0+np.exp((T-Tm)/w))
    return peak*death

def precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm):
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        base=kmap.get(r.id)
        kcat_s,bg=base if base else (13.7,None)
        base_coef=statistics.mean(mws)/(kcat_s*3600.0)
        Topt=tov.get(bg,med_to); Tm=tmv.get(bg,med_tm)
        terms.append((r.forward_variable,r.reverse_variable,base_coef,Topt,Tm))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return terms,con

def growth_at_T(m,terms,con,T_C,sig,w,P):
    coeffs={}
    for fv,rv,bc,Topt,Tm in terms:
        a=max(act(T_C,Topt,Tm,sig,w),1e-6); c=bc/a
        coeffs[fv]=c; coeffs[rv]=c
    con.ub=P; con.set_linear_coefficients(coeffs)
    try: g=m.slim_optimize()
    except Exception: g=float('nan')
    return g if g==g else 0.0

def main():
    # honest measured TPCs: dead wells counted as zero growth (thermal-tolerance
    # curve), not survivor-only medians. Columns species,T,mu_honest.
    hon=pd.read_csv(TABLES / 'measured_tpc_honest.csv')
    meas=hon.rename(columns={'mu_honest':'median'}).copy()
    meas['Species']=meas['species'].map({'auris':'auris','haemulonii':'Hae',
        'duobushaemulonii':'Duo','parapsilosis':'para'})
    S={}
    for sp in SP:
        tov,tmv,med_to,med_tm=therm_maps(sp); m,EX,mw,kmap=setup_pool(sp)
        terms,con=precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm)
        S[sp]=dict(m=m,terms=terms,con=con)
        print(f"setup {sp}: {len(terms)} enzyme-costed reactions",flush=True)
    au=meas[meas.Species=='auris'].dropna(subset=['median']); Ts=au['T'].values; mu=au['median'].values
    d=S['auris']; neval=[0]
    def loss(x):
        sig,w,P=x
        if sig<=1 or w<=0.3 or P<=0.02: return 1e3
        pred=np.array([growth_at_T(d['m'],d['terms'],d['con'],t,sig,w,P) for t in Ts])
        sc=mu.max()/max(pred.max(),1e-6); neval[0]+=1
        e=float(np.mean((pred*sc-mu)**2))
        if neval[0]%25==0: print(f"  eval {neval[0]} sig={sig:.1f} w={w:.1f} P={P:.3f} mse={e:.4f}",flush=True)
        return e
    best=None
    for x0 in ([8,3,0.20],[12,5,0.25],[6,2,0.15],[15,8,0.30]):
        r=minimize(loss,x0,method='Nelder-Mead',options=dict(maxiter=120,xatol=.3,fatol=1e-4))
        if best is None or r.fun<best.fun: best=r
    sig,w,P=best.x
    predau=np.array([growth_at_T(d['m'],d['terms'],d['con'],t,sig,w,P) for t in Ts])
    SCALE=mu.max()/max(predau.max(),1e-6)
    json.dump(dict(sig=sig,w=w,P=P,SCALE=SCALE,auris_fit_mse=best.fun),open(TABLES / 'etcgem_calib.json','w'),indent=2)
    grid=sorted(meas['T'].unique()); rows=[]
    for sp in SP:
        d=S[sp]
        for t in grid:
            g=growth_at_T(d['m'],d['terms'],d['con'],t,sig,w,P)*SCALE
            rows.append(dict(species=sp,T=t,pred_mu=round(g,4)))
    pred=pd.DataFrame(rows); pred.to_csv(TABLES / 'etcgem_tpc_pred.csv',index=False)
    print("calibration:",json.load(open(TABLES / 'etcgem_calib.json')))
    piv=pred.pivot_table(index='T',columns='species',values='pred_mu')
    print(piv.round(2).to_string())
    print("\npredicted peak T per species:")
    for sp in SP: print(f"  {sp}: peak {piv[sp].max():.2f} at {piv[sp].idxmax():.0f}C")
    print("\nwrote etcgem_tpc_pred.csv, etcgem_calib.json")

if __name__=='__main__': main()
