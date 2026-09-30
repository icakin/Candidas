"""Step 2, parallel + resumable. Frozen windows; classifications read, never recomputed.
r_se = bootstrap SD of r (AR(1)-matched). Downstream quantity only."""
import pandas as pd, numpy as np, sys, os, time, hashlib
from multiprocessing import Pool
sys.path.insert(0,'.')
import classifier as C
BUDGET=float(sys.argv[1]) if len(sys.argv)>1 else 150.0
OUT="automatic_continuous_partial.csv"
st=pd.read_csv("unblinded_results.csv"); tr=pd.read_csv("blind_traces.csv.gz")
G={w:(g.Time.values,g.Oxygen.values) for w,g in tr.groupby("WELLID")}
def one(w):
    t,y=G[w]; m,_=C.usable_interval(t,y)
    rec=dict(WELLID=w,win_start=np.nan,win_end=np.nan,dur_h_auto=np.nan,n_pts=np.nan,
             r_auto=np.nan,K_auto=np.nan,r_se=np.nan,resid_sd=np.nan,O2_0_auto=np.nan,dd_auto=np.nan)
    if m is None: return rec
    tt,yy=C.decimate(t[m],y[m]); tr2=tt-tt[0]
    rec.update(win_start=float(tt[0]),win_end=float(tt[-1]),dur_h_auto=(tt[-1]-tt[0])/60.0,
               n_pts=len(tt),dd_auto=float(yy[0]-yy.min()))
    f2=C.fit_m2(tr2,yy)
    if f2 is None: return rec
    rho=C.lag1(f2.fun); sd=float(np.std(f2.fun))
    rec.update(O2_0_auto=float(f2.x[0]),K_auto=float(f2.x[1]),r_auto=float(f2.x[2]),resid_sd=sd)
    rr=np.random.default_rng(int(hashlib.md5(w.encode()).hexdigest()[:8],16))
    fit=C._m2(f2.x,tr2); boot=[]
    for _ in range(C.N_BOOT):
        e=np.empty(len(tr2)); e[0]=rr.normal(0,sd)
        nz=rr.normal(0,sd*np.sqrt(max(1-rho**2,1e-6)),len(tr2))
        for q in range(1,len(tr2)): e[q]=rho*e[q-1]+nz[q]
        b=C.fit_m2(tr2,fit+e,start=f2.x)
        if b is not None: boot.append(b.x[2])
    if len(boot)>=50: rec["r_se"]=float(np.std(boot))
    return rec
done=set(pd.read_csv(OUT).WELLID) if os.path.exists(OUT) else set()
todo=[w for w in st.WELLID if w not in done]
print(f"{len(done)}/{len(st)} done, {len(todo)} to go",flush=True)
t0=time.time(); rows=[]
with Pool(4) as p:
    for r in p.imap_unordered(one,todo,chunksize=4):
        rows.append(r)
        if time.time()-t0>BUDGET: break
    p.terminate()
if rows: pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
print(f"now {len(pd.read_csv(OUT))}/{len(st)}",flush=True)
