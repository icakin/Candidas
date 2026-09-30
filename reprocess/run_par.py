"""Parallel resumable runner (4 workers). Identical frozen classifier; parallelism
is an implementation detail and changes no threshold."""
import pandas as pd, numpy as np, sys, os, time, hashlib
from multiprocessing import Pool
sys.path.insert(0,'.')
import classifier as C
MODE=sys.argv[1]; BUDGET=float(sys.argv[2]) if len(sys.argv)>2 else 150.0
R_GRID=[0.0,0.005,0.01,0.02,0.05,0.1,0.2,0.4,0.8]; N_PER=200
def sim_one(k):
    rr=SIM_RR; NOISE=SIM_NOISE
    r_true,_=[(r,i) for r in R_GRID for i in range(N_PER)][k] if False else (R_GRID[k//N_PER],k%N_PER)
    s=np.random.default_rng(1000003*k+7)
    j=s.integers(len(rr)); O0=rr.O0.iloc[j]; dur_h=rr.dur_h.iloc[j]; dd=rr.dd.iloc[j]
    rho=rr.rho.iloc[j]; sig=NOISE[s.integers(len(NOISE))]
    t=np.arange(0,dur_h*60+1.0,1.0); rm=r_true/60.0
    K=(dd*rm/np.expm1(rm*t[-1])) if rm>1e-12 else dd/t[-1]
    clean=C._m2([O0,K,rm],t)
    e=np.empty(len(t)); e[0]=s.normal(0,sig)
    nz=s.normal(0,sig*np.sqrt(max(1-rho**2,1e-6)),len(t))
    for q in range(1,len(t)): e[q]=rho*e[q-1]+nz[q]
    floor=max(0.05, rr.O0.iloc[j]-dd)   # v4: realistic plateau floor
    y=np.clip(clean+e,floor,None)
    o=C.classify(np.arange(len(t),dtype=float)+95.0,y,seed=int(s.integers(1e9)))
    return dict(k=k,r_true=r_true,rt_true=r_true*dur_h,dd_true=dd,dur_h=dur_h,sigma=sig,
                rho_true=rho,**{kk:o[kk] for kk in ("state","reason","drawdown","rt","r_hat",
                "r_boot_lo","dlpd","dlpd_se","G","eff_n")})
def real_one(w):
    g=REAL_G[w]
    o=C.classify(g[0],g[1],seed=int(hashlib.md5(w.encode()).hexdigest()[:8],16))
    return dict(WELLID=w,**o)
if MODE=="sim":
    SIM_RR=pd.read_csv("blind_trace_ranges.csv"); SIM_NOISE=pd.read_csv("blind_noise.csv").sigma.values
    OUT="sim_partial.csv"; N=len(R_GRID)*N_PER
    done=set(pd.read_csv(OUT).k) if os.path.exists(OUT) else set()
    todo=[k for k in range(N) if k not in done]; fn=sim_one
else:
    tr=pd.read_csv("blind_traces.csv.gz")
    REAL_G={w:(g.Time.values,g.Oxygen.values) for w,g in tr.groupby("WELLID")}
    OUT="blind_states_partial.csv"; N=len(REAL_G)
    done=set(pd.read_csv(OUT).WELLID) if os.path.exists(OUT) else set()
    todo=[w for w in sorted(REAL_G) if w not in done]; fn=real_one
print(f"{MODE}: {N-len(todo)}/{N} done, {len(todo)} to go",flush=True)
t0=time.time(); rows=[]
with Pool(4) as p:
    for res in p.imap_unordered(fn,todo,chunksize=4):
        rows.append(res)
        if time.time()-t0>BUDGET: break
    p.terminate()
if rows: pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
print(f"{MODE}: now {len(pd.read_csv(OUT))}/{N}",flush=True)
