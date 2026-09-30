"""PROTOCOL v1 sec.6 / v1.1 acceptance criteria A1-A6.
Simulated traces use ONLY blind trace-shape distributions (no labels)."""
import numpy as np, pandas as pd, sys; sys.path.insert(0,'.')
import classifier as C
rng=np.random.default_rng(20260921)
rr=pd.read_csv("blind_trace_ranges.csv")
# measurement noise from second differences (robust to smooth trend), estimated blind
NOISE=pd.read_csv("blind_noise.csv").sigma.values
R_GRID=[0.0,0.005,0.01,0.02,0.05,0.1,0.2,0.4,0.8]   # h^-1
N_PER=200
rows=[]
for r_true in R_GRID:
    for i in range(N_PER):
        j=rng.integers(len(rr)); O0=rr.O0.iloc[j]; dur_h=rr.dur_h.iloc[j]; dd=rr.dd.iloc[j]
        rho=rr.rho.iloc[j]; sig=NOISE[rng.integers(len(NOISE))]
        dt=1.0; t=np.arange(0,dur_h*60+dt,dt)
        rm=r_true/60.0
        K=(dd*rm/np.expm1(rm*t[-1])) if rm>1e-12 else dd/t[-1]
        clean=C._m2([O0,K,rm],t)
        e=np.empty(len(t)); e[0]=rng.normal(0,sig)
        nz=rng.normal(0,sig*np.sqrt(max(1-rho**2,1e-6)),len(t))
        for k in range(1,len(t)): e[k]=rho*e[k-1]+nz[k]
        y=np.clip(clean+e,0.02,None)
        tfull=np.arange(0,len(t)*dt,dt)+95.0     # place after equilibration
        o=C.classify(tfull,y,seed=int(rng.integers(1e9)))
        rows.append(dict(r_true=r_true,rt_true=r_true*dur_h,dd_true=dd,dur_h=dur_h,
                         sigma=sig,rho_true=rho,**{k:o[k] for k in
                         ("state","reason","drawdown","rt","r_hat","r_boot_lo","dlpd","dlpd_se","G","eff_n")}))
    print(f"  r={r_true}: done",flush=True)
d=pd.DataFrame(rows); d.to_csv("sim_results.csv",index=False)
print("\n=== state by true r ===")
print(pd.crosstab(d.r_true,d.state,normalize="index").round(3).to_string())
print("\n=== ACCEPTANCE (v1.1) ===")
z=d[d.r_true==0]; a1=(z.state==1).mean()
hi=d[d.rt_true>=1.0]; a2=(hi.state==1).mean(); a4=(hi.state==2).mean()
nb=d[(d.rt_true>=0.4)&(d.rt_true<=0.6)]; a3=(nb.state==2).mean() if len(nb) else np.nan
a5=(d.state==3).mean(); big=d[d.dd_true>2.0]; a6=(big.state==4).mean()
for nm,v,op,lim in (("A1 false growth at r=0",a1,"<=",0.05),("A2 power at rt>=1",a2,">=",0.80),
                    ("A3 STATE2 at rt~0.5",a3,"<=",0.50),("A4 STATE2 at rt>=1",a4,"<=",0.10),
                    ("A5 ambiguous overall",a5,"<=",0.25),("A6 STATE4 when dd>2",a6,"<=",0.05)):
    ok = (v<=lim) if op=="<=" else (v>=lim)
    print(f"  {nm:26s} {v:.3f} {op} {lim}   {'PASS' if ok else 'FAIL'}")
