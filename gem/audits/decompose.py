#!/usr/bin/env python3
"""decompose.py - unified OUT-OF-FOLD performance decomposition.
Uses the proven persistent-constraint pattern (B.precompute + B.growth_at_T).
Two passes: genome-specific and identical-null; plus a trivial median-TPC baseline.
Metrics per held-out species: raw RMSE, max-normalised RMSE, optimum, predicted
growth at first observed no-growth temperature.
"""
import build_etcgem_tpc as B, numpy as np, pandas as pd, json, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); W=cal['w']; P=cal['P']
hon=pd.read_csv(f'{TABLES}/measured_tpc_honest.csv')
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
GTo=float(pd.read_csv(f'{TABLES}/thermal_topt.csv')['pred_topt'].median())
GTm=float(pd.read_csv(f'{TABLES}/thermal_tm.csv')['pred_tm'].median())
grid=np.arange(22,44.1,1.0)

def make(null):
    """persistent precompute per species; null=True uses global-median Topt/Tm."""
    S={}
    for sp in SPP:
        tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
        if null: tov,tmv,med_to,med_tm={},{},GTo,GTm
        terms,con=B.precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm)
        S[sp]=(m,terms,con)
    return S
def curve(S,sp,sig):
    m,terms,con=S[sp]; return {float(t):B.growth_at_T(m,terms,con,float(t),sig,W,P) for t in grid}
def fit_sig(S,train):
    def loss(sig):
        e=0.0
        for sp in train:
            c=curve(S,sp,sig); mm=hon[hon.species==sp].set_index('T')['mu_honest']
            pv=np.array([c[float(t)] for t in mm.index]); sc=mm.max()/max(pv.max(),1e-6)
            e+=float(np.mean((pv*sc-mm.values)**2))
        return e
    return float(minimize_scalar(loss,bounds=(4,16),method='bounded',options=dict(maxiter=12)).x)

def oof_curves(S):
    res={}
    for held in SPP:
        sig=fit_sig(S,[s for s in SPP if s!=held])
        c=curve(S,held,sig); mm=hon[hon.species==held].set_index('T')['mu_honest']
        sc=mm.max()/max(max(c.values()),1e-6)
        res[held]={float(t):c[float(t)]*sc for t in grid}
    return res

print("genome pass...",flush=True); Sg=make(False); GEN=oof_curves(Sg)
print("null pass...",flush=True);   Sn=make(True);  NUL=oof_curves(Sn)

rows=[]; oof={}
for held in SPP:
    mm=hon[hon.species==held].set_index('T')['mu_honest']; Ts=mm.index.values; obs=mm.values
    cg=GEN[held]; cn=NUL[held]
    pg=np.array([cg[float(t)] for t in Ts]); pn=np.array([cn[float(t)] for t in Ts])
    # trivial median-TPC baseline = mean of other 3 measured curves
    tri=[]
    for t in Ts:
        v=[hon[(hon.species==s)&(hon.T==t)]['mu_honest'].values for s in SPP if s!=held]
        v=[x[0] for x in v if len(x)]; tri.append(np.mean(v) if v else np.nan)
    tri=np.array(tri)
    def rmse(a,b): mk=~np.isnan(b); return float(np.sqrt(np.mean((np.asarray(a)[mk]-b[mk])**2)))
    def nrmse(a,b):
        mk=~np.isnan(b); an=np.asarray(a)/max(np.asarray(a).max(),1e-6); bn=b/max(np.nanmax(b),1e-6)
        return float(np.sqrt(np.mean((an[mk]-bn[mk])**2)))
    nz=mm[mm<=0.05]; fT=float(nz.index.min()) if len(nz) else None
    def at(cv,T): return cv[float(T)] if T is not None else np.nan
    opt_tri = int(Ts[np.nanargmax(tri)]) if not np.all(np.isnan(tri)) else None
    rows.append(dict(species=held, opt_meas=int(mm.idxmax()),
        opt_gen=int(Ts[np.argmax(pg)]), opt_null=int(Ts[np.argmax(pn)]), opt_tri=opt_tri,
        rmse_gen=round(rmse(pg,obs),3), rmse_null=round(rmse(pn,obs),3), rmse_tri=round(rmse(tri,obs),3),
        nrmse_gen=round(nrmse(pg,obs),3), nrmse_null=round(nrmse(pn,obs),3),
        firstNoGrowT=fT, mu_nogrow_gen=round(at(cg,fT),3) if fT else None,
        mu_nogrow_null=round(at(cn,fT),3) if fT else None))
    oof[held]=dict(T=[float(t) for t in grid],gen=[cg[float(t)] for t in grid],null=[cn[float(t)] for t in grid])
df=pd.DataFrame(rows); df.to_csv(f'{TABLES}/oof_decomposition.csv',index=False); json.dump(oof,open(f'{TABLES}/oof_curves.json','w'))
print(df.to_string(index=False))
print("\n== mean over folds ==")
tri_ok=df.dropna(subset=['opt_tri'])
print(f"optimum MAE:  genome {np.mean(np.abs(df.opt_gen-df.opt_meas)):.1f}  null {np.mean(np.abs(df.opt_null-df.opt_meas)):.1f}  trivial {np.mean(np.abs(tri_ok.opt_tri-tri_ok.opt_meas)):.1f} (n={len(tri_ok)})")
print(f"raw RMSE:     genome {df.rmse_gen.mean():.3f}  null {df.rmse_null.mean():.3f}  trivial {df.rmse_tri.mean():.3f}")
print(f"norm RMSE:    genome {df.nrmse_gen.mean():.3f}  null {df.nrmse_null.mean():.3f}")
gg=df.dropna(subset=['mu_nogrow_gen'])
print(f"mu at first no-growth T: genome {gg.mu_nogrow_gen.mean():.2f}  null {gg.mu_nogrow_null.mean():.2f}  (observed 0)")
print("wrote oof_decomposition.csv, oof_curves.json")
