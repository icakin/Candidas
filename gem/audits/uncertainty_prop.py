#!/usr/bin/env python3
"""uncertainty_prop.py - propagate Seq2Topt/Seq2Tm predictor error into the etcGEM
TPCs (vectorised). Independent per-enzyme jitter: sigma_topt=12 C (Seq2Topt RMSE),
sigma_tm=7 C (Seq2Tm, R2~0.76). A systematic per-species bias is common-mode and
cancels in the interspecies comparison, so independent jitter is the correct model
for the DIFFERENCE/compression claim. Reports per-T 5/50/95% bands + optimum CI.
"""
import build_etcgem_tpc as B, numpy as np, pandas as pd, json, warnings, statistics
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
warnings.filterwarnings('ignore')
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
S_TOPT,S_TM=12.0,7.0; NREP=20; CLIP=20.0
grid=np.arange(22,44.1,2.0); rng=np.random.default_rng(7)

def skeleton(sp):
    tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
    fvs=[];rvs=[];bc=[];To=[];Tm=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        base=kmap.get(r.id); ks,bg=base if base else (13.7,None)
        fvs.append(r.forward_variable); rvs.append(r.reverse_variable)
        bc.append(statistics.mean(mws)/(ks*3600.0)); To.append(tov.get(bg,med_to)); Tm.append(tmv.get(bg,med_tm))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return m,fvs,rvs,np.array(bc),np.array(To),np.array(Tm),con

def act_vec(T,To,Tm,sig,w):
    return np.exp(-((T-To)**2)/(2*sig*sig))*(1.0/(1.0+np.exp((T-Tm)/w)))

out={}
for sp in SPP:
    m,fvs,rvs,bc,To0,Tm0,con=skeleton(sp)
    curves=np.zeros((NREP,len(grid))); opts=np.zeros(NREP); con.ub=P0
    for k in range(NREP):
        To=To0+np.clip(rng.normal(0,S_TOPT,len(To0)),-CLIP,CLIP); Tm=Tm0+np.clip(rng.normal(0,S_TM,len(Tm0)),-CLIP,CLIP)
        row=[]
        for t in grid:
            a=np.maximum(act_vec(t,To,Tm,SIG,W),1e-3); c=bc/a
            co={}
            for i in range(len(fvs)): co[fvs[i]]=c[i]; co[rvs[i]]=c[i]
            con.set_linear_coefficients(co); g=m.slim_optimize(); row.append((g if g==g else 0.0)*SC)
        curves[k]=row; opts[k]=grid[int(np.argmax(row))]
        if (k+1)%15==0: print(f"  {sp} {k+1}/{NREP}",flush=True)
    out[sp]=dict(T=grid.tolist(),lo=np.percentile(curves,5,0).tolist(),
                 md=np.percentile(curves,50,0).tolist(),hi=np.percentile(curves,95,0).tolist(),
                 opt_med=float(np.median(opts)),opt_lo=float(np.percentile(opts,5)),opt_hi=float(np.percentile(opts,95)))
    print(f"{sp}: optimum {out[sp]['opt_med']:.0f}C [{out[sp]['opt_lo']:.0f}-{out[sp]['opt_hi']:.0f}]",flush=True)
json.dump(out,open(f'{TABLES}/uncertainty_bands.json','w')); print("wrote uncertainty_bands.json")
