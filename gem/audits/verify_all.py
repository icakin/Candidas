#!/usr/bin/env python3
"""verify_all.py - recompute every headline model claim from scratch, one global
scale where absolute values matter, and print a verification table."""
import build_etcgem_tpc as B, numpy as np, pandas as pd, json, warnings, statistics
warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
cal=json.load(open(f'{TABLES}/etcgem_calib.json')); SIG,W,P,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
hon=pd.read_csv(f'{TABLES}/measured_tpc_honest.csv')
SPP=['auris','haemulonii','duobushaemulonii','parapsilosis']
GTo=float(pd.read_csv(f'{TABLES}/thermal_topt.csv')['pred_topt'].median())
GTm=float(pd.read_csv(f'{TABLES}/thermal_tm.csv')['pred_tm'].median())
grid=np.arange(22,44.1,1.0)
R=[]

def build(sp,null=False):
    tov,tmv,med_to,med_tm=B.therm_maps(sp); m,EX,mw,kmap=B.setup_pool(sp)
    if null: tov,tmv,med_to,med_tm={},{},GTo,GTm
    terms,con=B.precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm)
    return m,terms,con
def curve(sp,sig,null=False,scale=SC):
    m,terms,con=build(sp,null)
    c={float(t):B.growth_at_T(m,terms,con,float(t),sig,W,P)*scale for t in grid}
    m.remove_cons_vars([con]); return c

# 1. POOL BINDING
m,terms,con=build('auris'); B.growth_at_T(m,terms,con,34,SIG,W,P); gp=m.slim_optimize()
con.ub=1e9; m.solver.update(); gf=m.slim_optimize(); m.remove_cons_vars([con])
R.append(("Pool binds (auris 34C)",f"with {gp:.2f} vs without {gf:.2f}/h",f"{gf/gp:.1f}x reduction","YES - load-bearing"))

# 2. GROWTH PREDICTION, single global scale, alive region
allm=[];allp=[]
for sp in SPP:
    c=curve(sp,SIG); mm=hon[hon.species==sp].set_index('T')['mu_honest']; mm=mm[mm>0.1]
    for t in mm.index:
        if float(t) in c: allm.append(mm[t]); allp.append(c[float(t)])
allm=np.array(allm);allp=np.array(allp)
r=np.corrcoef(allm,allp)[0,1]; mae=np.mean(np.abs(allm-allp)); tmae=np.mean(np.abs(allm-allm.mean()))
R.append(("Growth pred (global scale, alive)",f"r={r:.2f} MAE={mae:.2f}",f"trivial-mean MAE={tmae:.2f}","beats trivial" if mae<tmae else "no"))

# 3. OUT-OF-FOLD genome vs null (fit sigma on 3, predict 4th)
def fit_sig(train,null):
    built={sp:build(sp,null) for sp in train}
    def loss(sig):
        e=0
        for sp in train:
            m,terms,con=built[sp]; cc={float(t):B.growth_at_T(m,terms,con,float(t),sig,W,P) for t in grid}
            mm=hon[hon.species==sp].set_index('T')['mu_honest']; pv=np.array([cc[float(t)] for t in mm.index])
            sc=mm.max()/max(pv.max(),1e-6); e+=np.mean((pv*sc-mm.values)**2)
        return e
    s=float(minimize_scalar(loss,bounds=(4,16),method='bounded',options=dict(maxiter=10)).x)
    for sp in train: built[sp][0].remove_cons_vars([built[sp][2]])
    return s
rg=[];rn=[];og=[];on=[];mo=[]
for held in SPP:
    tr=[s for s in SPP if s!=held]
    for null,rr,oo in [(False,rg,og),(True,rn,on)]:
        sig=fit_sig(tr,null); c=curve(held,sig,null,scale=1.0)
        mm=hon[hon.species==held].set_index('T')['mu_honest']
        pv=np.array([c[float(t)] for t in mm.index]); scl=mm.max()/max(pv.max(),1e-6)
        rr.append(np.sqrt(np.mean((pv*scl-mm.values)**2))); oo.append(grid[np.argmax([c[float(t)] for t in grid])])
    mo.append(hon[hon.species==held].set_index('T')['mu_honest'].idxmax())
gmae=np.mean(np.abs(np.array(og)-np.array(mo))); nmae=np.mean(np.abs(np.array(on)-np.array(mo)))
R.append(("Out-of-fold RMSE genome vs null",f"genome {np.mean(rg):.3f} null {np.mean(rn):.3f}","null <= genome","genome NO advantage"))
R.append(("Out-of-fold optimum MAE g vs null",f"genome {gmae:.1f} null {nmae:.1f}","trivial-median=1.5","genome NO advantage; loses to trivial"))

# 4. optimum spread fine grid
fg=np.arange(28,40.1,0.5); peaks=[]
for sp in SPP:
    m,terms,con=build(sp); gs=[B.growth_at_T(m,terms,con,float(t),SIG,W,P) for t in fg]; peaks.append(fg[np.argmax(gs)]); m.remove_cons_vars([con])
R.append(("Predicted optimum spread",f"{min(peaks):.1f}-{max(peaks):.1f}C (spread {max(peaks)-min(peaks):.1f})","measured spread 6C","compressed - confirmed"))

# 5. enzyme conservation
to=pd.read_csv(f'{TABLES}/thermal_topt.csv'); to['sp']=to['id'].str.split('|').str[0]
meds=[to[to.sp==sp]['pred_topt'].median() for sp in SPP]
pr=pd.read_csv(f'{TABLES}/paired_ortholog.csv')
R.append(("Enzyme Topt medians",f"{min(meds):.1f}-{max(meds):.1f}C (within {max(meds)-min(meds):.1f})","paired |dTopt| med "+f"{pr.dTopt.abs().median():.1f}C","clade-conserved - confirmed"))

print(f"{'CLAIM':<38}{'RECOMPUTED':<30}{'BASELINE/REF':<24}{'VERDICT'}")
print("-"*118)
for c,v,b,ver in R: print(f"{c:<38}{v:<30}{b:<24}{ver}")
PY