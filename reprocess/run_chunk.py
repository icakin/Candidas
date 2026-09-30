"""Resumable runner. Processes until the time budget, appends, exits. Re-run to continue."""
import pandas as pd, numpy as np, sys, os, time, hashlib; sys.path.insert(0,'.')
import classifier as C
MODE=sys.argv[1]; BUDGET=float(sys.argv[2]) if len(sys.argv)>2 else 150.0
t0=time.time()
if MODE=="real":
    OUT="blind_states_partial.csv"
    tr=pd.read_csv("blind_traces.csv.gz")
    allw=sorted(tr.WELLID.unique())
    done=set(pd.read_csv(OUT).WELLID) if os.path.exists(OUT) else set()
    todo=[w for w in allw if w not in done]
    print(f"real: {len(done)} done, {len(todo)} remaining",flush=True)
    grp=dict(list(tr.groupby("WELLID"))); rows=[]
    for w in todo:
        g=grp[w]
        o=C.classify(g.Time.values,g.Oxygen.values,seed=int(hashlib.md5(w.encode()).hexdigest()[:8],16))
        rows.append(dict(WELLID=w,**o))
        if time.time()-t0>BUDGET: break
    if rows:
        pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
    n=len(pd.read_csv(OUT)); print(f"real: now {n}/{len(allw)}",flush=True)
else:
    OUT="sim_partial.csv"
    rng=np.random.default_rng(20260921)
    rr=pd.read_csv("blind_trace_ranges.csv"); NOISE=pd.read_csv("blind_noise.csv").sigma.values
    R_GRID=[0.0,0.005,0.01,0.02,0.05,0.1,0.2,0.4,0.8]; N_PER=200
    jobs=[(r,i) for r in R_GRID for i in range(N_PER)]
    done=len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0
    print(f"sim: {done}/{len(jobs)} done",flush=True)
    rows=[]
    for k in range(done,len(jobs)):
        r_true,i=jobs[k]
        s=np.random.default_rng(1000003*k+7)
        j=s.integers(len(rr)); O0=rr.O0.iloc[j]; dur_h=rr.dur_h.iloc[j]; dd=rr.dd.iloc[j]
        rho=rr.rho.iloc[j]; sig=NOISE[s.integers(len(NOISE))]
        t=np.arange(0,dur_h*60+1.0,1.0); rm=r_true/60.0
        K=(dd*rm/np.expm1(rm*t[-1])) if rm>1e-12 else dd/t[-1]
        clean=C._m2([O0,K,rm],t)
        e=np.empty(len(t)); e[0]=s.normal(0,sig)
        nz=s.normal(0,sig*np.sqrt(max(1-rho**2,1e-6)),len(t))
        for q in range(1,len(t)): e[q]=rho*e[q-1]+nz[q]
        y=np.clip(clean+e,0.02,None)
        o=C.classify(np.arange(len(t),dtype=float)+95.0,y,seed=int(s.integers(1e9)))
        rows.append(dict(r_true=r_true,rt_true=r_true*dur_h,dd_true=dd,dur_h=dur_h,sigma=sig,
                         rho_true=rho,**{kk:o[kk] for kk in ("state","reason","drawdown","rt",
                         "r_hat","r_boot_lo","dlpd","dlpd_se","G","eff_n")}))
        if time.time()-t0>BUDGET: break
    if rows: pd.DataFrame(rows).to_csv(OUT,mode="a",header=not os.path.exists(OUT),index=False)
    n=len(pd.read_csv(OUT)); print(f"sim: now {n}/{len(jobs)}",flush=True)
