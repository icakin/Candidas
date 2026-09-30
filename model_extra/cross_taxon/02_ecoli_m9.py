"""E. coli in R2A, 12 temperatures 15-50 C (data/Bacteria/Ecoli). Per-well fit of the same
oxygen-dynamics model used for the Candida: O2(t) = O2_0 - (K/r)(exp(rt)-1), then
Sharpe-Schoolfield on r and Arrhenius on K, with a stratified well bootstrap."""
import glob, re, os, numpy as np, pandas as pd
from scipy.optimize import least_squares
k=8.617e-5; TREF=293.15; NBOOT=600; MIN_DD=0.5
rng=np.random.default_rng(11)
D=os.path.expanduser("~/mnt/Candidas/data/Bacteria/Ecoli")

def fit_well(t, y):
    """window: from the post-equilibration peak to where the trace flattens."""
    i0=int(np.argmax(y[:max(5,len(y)//10)]))          # end of the warming rise
    t=t[i0:]; y=y[i0:]
    if len(t)<30: return None
    ymin=y.min(); thr=ymin+0.02*(y[0]-ymin)
    below=np.where(y<=thr)[0]
    i1=below[0] if len(below) else len(y)-1
    t=t[:i1+1]-t[0]; y=y[:i1+1]
    if len(t)<30 or (y[0]-y[-1])<MIN_DD: return None
    def f(p):
        O,K,r=p
        return (O-(K/r)*(np.expm1(r*t)) if abs(r)>1e-9 else O-K*t)-y
    best=None
    for r0 in (0.0005,0.003,0.01):
        try:
            res=least_squares(f,[y[0],max((y[0]-y[-1])/t[-1],1e-6),r0],
                              bounds=([-5,1e-9,1e-6],[25,5,0.5]))
            if best is None or res.cost<best.cost: best=res
        except Exception: pass
    if best is None: return None
    O,K,r=best.x
    ss=np.sum(best.fun**2); r2=1-ss/np.sum((y-y.mean())**2)
    if r2<0.90 or r>0.49 or K>4.9: return None
    return dict(r=r*60, K=K*60, r2=r2, n=len(t))   # per hour

rows=[]
for f in sorted(glob.glob(f"{D}/*.csv"), key=lambda x:int(re.search(r"Ecoli_(\d+)_",x).group(1))):
    T=int(re.search(r"Ecoli_(\d+)_",f).group(1)); d=pd.read_csv(f)
    t=d.Time.values.astype(float)
    for c in [c for c in d.columns if c.startswith("OTU")]:
        y=pd.to_numeric(d[c],errors="coerce").values
        m=np.isfinite(y)&np.isfinite(t)
        if m.sum()<40: continue
        out=fit_well(t[m],y[m])
        if out: rows.append(dict(T=T,well=c,strain=c.split("_")[0],**out))
w=pd.DataFrame(rows)
print("wells fitted per temperature and strain:")
print(w.groupby(["strain","T"]).size().unstack(fill_value=0).to_string()); print()
print(w.groupby(["strain","T"]).agg(r_h=("r","median"),K_h=("K","median")).round(4).to_string()); print()

def ss_topt(E,Eh,Th):
    return 1/(1/Th-np.log(E/(Eh-E))/(Eh/k))-273.15 if 0<E<Eh else np.nan
def fit_curves(T,r,K):
    x=1/(k*TREF)-1/(k*(T+273.15)); Tk=T+273.15
    g=least_squares(lambda p: p[0]+p[1]*x-np.log1p(np.exp(np.clip(p[2]*(1/(k*p[3])-1/(k*Tk)),-50,50)))-np.log(r),
                    [np.log(np.median(r)),0.8,3.0,np.median(Tk)+5],
                    bounds=([-40,0.05,0.3,285],[30,4,25,345]))
    lnB0,E,Eh,Th=g.x
    a=least_squares(lambda p:(p[0]+p[1]*x)-np.log(K),[np.log(np.median(K)),0.5],bounds=([-80,0.0],[30,3]))
    return E,a.x[1],ss_topt(E,Eh,Th),ss_topt(E-a.x[1],Eh,Th)

for strain,s in list(w.groupby("strain"))+[("POOLED",w)]:
    E,ER,Topt,Tc=fit_curves(s["T"].values.astype(float),s.r.values,s.K.values)
    B=[]
    for _ in range(NBOOT):
        b=s.groupby("T",group_keys=False).apply(lambda u:u.sample(len(u),replace=True,random_state=rng.integers(1e9)))
        try: B.append(fit_curves(b["T"].values.astype(float),b.r.values,b.K.values))
        except Exception: pass
    B=np.array(B); tc=B[:,3]; tc=tc[np.isfinite(tc)]
    print(f"{strain:8s} n={len(s):3d}  E={E:.2f}  E_R={ER:.2f}  Topt={Topt:.1f} "
          f"[{np.nanpercentile(B[:,2],2.5):.1f},{np.nanpercentile(B[:,2],97.5):.1f}]  "
          f"T_carbon={Tc:.1f} [{np.percentile(tc,2.5):.1f},{np.percentile(tc,97.5):.1f}]  "
          f"P(T_carbon<37)={np.mean(tc<37):.3f}")
w.to_csv(os.path.expanduser("~/mnt/Candidas/model_extra/cross_taxon/ecoli_r2a_wells.csv"),index=False)
