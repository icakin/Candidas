"""Resumable runner for lik_model.analyse over blind_traces.csv.gz (4 workers, time budget)."""
import pandas as pd, numpy as np, sys, os, hashlib, time; sys.path.insert(0,'.')
from multiprocessing import Pool
import lik_model as L
BUDGET=float(sys.argv[1]) if len(sys.argv)>1 else 140.0
OUT="lik_partial.csv"
tr=pd.read_csv("blind_traces.csv.gz")
G={w:(g.Time.values,g.Oxygen.values) for w,g in tr.groupby("WELLID")}
def one(w):
    t,y=G[w]; o=L.analyse(t,y,seed=int(hashlib.md5(w.encode()).hexdigest()[:8],16),n_boot=499)
    return dict(WELLID=w,**o)
if __name__=="__main__":
    done=set(pd.read_csv(OUT).WELLID) if os.path.exists(OUT) else set()
    todo=[w for w in sorted(G) if w not in done]
    print(f"{len(done)}/{len(G)} done, {len(todo)} to go",flush=True)
    t0=time.time(); rows=[]
    with Pool(4) as p:
        for res in p.imap_unordered(one,todo,chunksize=2):
            rows.append(res)
            if time.time()-t0>BUDGET: break
        p.terminate()
    if rows: pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
    n=len(pd.read_csv(OUT)); print(f"now {n}/{len(G)} ({time.time()-t0:.0f}s)",flush=True)
    if n==len(G): pd.read_csv(OUT).to_csv("lik_results.csv",index=False); print("COMPLETE")
