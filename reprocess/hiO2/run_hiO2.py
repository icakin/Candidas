"""High-oxygen refit: rerun the per-well likelihood model on each trace truncated where dissolved
oxygen first falls below FRAC of its value at the window start (90 min). Prespecified FRAC = 0.5
(top half of the drawdown). Wells whose truncated window is shorter than MIN_H hours or has
fewer than MIN_N points are reported as unidentifiable. Output: lik_hiO2_<frac>.csv"""
import pandas as pd, numpy as np, sys, os, hashlib, time
from multiprocessing import Pool
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import lik_model as L
FRAC=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
MIN_H=2.0; MIN_N=60; NB=199
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),f"lik_hiO2_{FRAC:.2f}.csv")
HERE=os.path.dirname(os.path.abspath(__file__)); tr=pd.read_csv(os.path.join(HERE,"..","blind_traces.csv.gz")); G={w:(g.Time.values,g.Oxygen.values) for w,g in tr.groupby("WELLID")}
def one(w):
    t,y=G[w]; ok=np.isfinite(t)&np.isfinite(y); t,y=t[ok],y[ok]
    m=(t>=L.START_MIN)&(t<=L.START_MIN+L.MAX_DUR_H*60); t,y=t[m],y[m]
    y0=np.median(y[:5]); cut=np.argmax(y<FRAC*y0) if (y<FRAC*y0).any() else len(y)
    tt,yy=t[:cut],y[:cut]
    base=dict(WELLID=w,frac=FRAC,n_kept=len(tt),dur_kept_h=(tt[-1]-tt[0])/60 if len(tt) else 0.0,identifiable=False)
    if len(tt)<MIN_N or base["dur_kept_h"]<MIN_H: return base
    o=L.analyse(tt,yy,seed=int(hashlib.md5(w.encode()).hexdigest()[:8],16),n_boot=NB)
    base.update(o); base["identifiable"]=bool(o.get("ok",False)); return base
if __name__=="__main__":
    done=set(pd.read_csv(OUT).WELLID) if os.path.exists(OUT) else set()
    todo=[w for w in sorted(G) if w not in done]; print(f"{len(done)} done, {len(todo)} to go",flush=True)
    t0=time.time(); rows=[]
    with Pool(2) as p:
        for k,res in enumerate(p.imap_unordered(one,todo,chunksize=4)):
            rows.append(res)
            if len(rows)%50==0:
                pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False); rows=[]
                print(f"{k+1}/{len(todo)} {time.time()-t0:.0f}s",flush=True)
    if rows: pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
    print("COMPLETE",time.time()-t0,flush=True)
