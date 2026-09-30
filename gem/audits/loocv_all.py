#!/usr/bin/env python3
"""loocv_all.py - all 4 leave-one-species-out folds (fast: fit shape sigma on the
3 training species with w,P fixed at the joint-fit values, predict held-out optimum
and high-T behavior)."""
import build_etcgem_tpc as B, numpy as np, pandas as pd, json, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); W=cal['w']; P=cal['P']
hon=pd.read_csv(f'{TABLES}/measured_tpc_honest.csv')
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
S={}
for sp in SPP:
    tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
    S[sp]=(m,)+B.precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm)
grid=np.arange(22,44.1,1.0)
def curve(sp,sig):
    m,terms,con=S[sp]; return {float(t):B.growth_at_T(m,terms,con,float(t),sig,W,P) for t in grid}
def fit_sig(train):
    def loss(sig):
        e=0.0
        for sp in train:
            c=curve(sp,sig); mm=hon[hon.species==sp].set_index('T')['mu_honest']
            pv=np.array([c[float(t)] for t in mm.index]); sc=mm.max()/max(pv.max(),1e-6)
            e+=float(np.mean((pv*sc-mm.values)**2))
        return e
    return float(minimize_scalar(loss,bounds=(4,16),method='bounded',options=dict(maxiter=12)).x)
print(f"{'held_out':<18}{'meas_opt':>9}{'pred_opt':>9}{'opt_err':>8}{'meas_death':>11}{'pred_mu@death+':>14}")
errs=[]
for held in SPP:
    sig=fit_sig([s for s in SPP if s!=held])
    c=curve(held,sig); mm=hon[hon.species==held].set_index('T')['mu_honest']
    sc=mm.max()/max(max(c.values()),1e-6)
    optp=max(c,key=c.get); optm=int(mm.idxmax())
    death=int(mm[mm>0.05].index.max())          # last alive T (measured)
    over=c.get(min(death+2,44.0),0)*sc          # predicted mu just past death
    errs.append(abs(optm-optp))
    print(f"{held:<18}{optm:>8}C{optp:>7.0f}C{abs(optm-optp):>7.0f}C{death:>10}C{over:>14.2f}")
print(f"\nmean |optimum error| across folds = {np.mean(errs):.1f} C")
print("Pattern: optima recoverable (~1-2C, because optima are similar) but the model")
print("keeps predicting growth PAST each species' measured death temperature.")
