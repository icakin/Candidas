"""R2A panel: E. coli (mammal-associated) and four environmental isolates, same medium,
same instrument, same model. Per-well oxygen fit, then Sharpe-Schoolfield on r and
Arrhenius on K, with a temperature-cluster bootstrap."""
import glob, re, os, numpy as np, pandas as pd
from scipy.optimize import least_squares
k=8.617e-5; TREF=293.15; rng=np.random.default_rng(5); NB=600
B=os.path.expanduser("~/mnt/Candidas/data/Bacteria")

def fit_well(t,y,min_dd=0.4):
    i0=int(np.argmax(y[:max(5,len(y)//10)])); t=t[i0:]; y=y[i0:]
    if len(t)<30: return None
    ymin=y.min(); thr=ymin+0.02*(y[0]-ymin); b=np.where(y<=thr)[0]
    i1=b[0] if len(b) else len(y)-1
    t=t[:i1+1]-t[0]; y=y[:i1+1]
    if len(t)<30 or (y[0]-y[-1])<min_dd: return None
    best=None
    for r0 in (0.0005,0.003,0.01):
        try:
            r=least_squares(lambda p:(p[0]-(p[1]/p[2])*np.expm1(p[2]*t))-y,[y[0],max((y[0]-y[-1])/t[-1],1e-6),r0],
                            bounds=([-5,1e-9,1e-6],[25,5,0.5]))
            if best is None or r.cost<best.cost: best=r
        except Exception: pass
    if best is None: return None
    O,K,r=best.x; r2=1-np.sum(best.fun**2)/np.sum((y-y.mean())**2)
    if r2<0.90 or r>0.49 or K>4.9: return None
    return r*60,K*60

def load(pattern,rx,strain_prefix=""):
    rows=[]
    for f in sorted(glob.glob(pattern),key=lambda x:int(re.search(rx,x).group(1))):
        T=int(re.search(rx,f).group(1)); d=pd.read_csv(f); t=d.Time.values.astype(float)
        for c in [c for c in d.columns if c.startswith("OTU")]:
            y=pd.to_numeric(d[c],errors="coerce").values; m=np.isfinite(y)&np.isfinite(t)
            if m.sum()<40: continue
            o=fit_well(t[m],y[m])
            if o: rows.append(dict(T=T,strain=strain_prefix+c.split("_")[0],r=o[0],K=o[1]))
    return pd.DataFrame(rows)

envs=[]
for folder in ("7_13_14_15","2_6_9_11","1_4_8_20"):
    ids=folder.split("_")
    e=load(f"{B}/{folder}/*.csv", r"R2A_(\d+)")
    e["strain"]=e.strain.map({f"OTU{i+1}":f"OTU{ids[i]}" for i in range(len(ids))})
    envs.append(e)
env=pd.concat(envs)
ec=load(f"{B}/Ecoli/*.csv", r"Ecoli_(\d+)_")
ec["strain"]=ec.strain.map({"OTU1":"E. coli A","OTU2":"E. coli B"})
w=pd.concat([ec,env]); w["growing"]=w.r>1e-3
print("wells fitted (growing / total):")
print(w.groupby(["strain","T"]).apply(lambda s: f"{int(s.growing.sum())}/{len(s)}").unstack(fill_value="-").to_string()); print()

def topt(E,Eh,Th): return 1/(1/Th-np.log(E/(Eh-E))/(Eh/k))-273.15 if 0<E<Eh else np.nan
def fit(g,a):
    T=g["T"].values.astype(float); x=1/(k*TREF)-1/(k*(T+273.15)); Tk=T+273.15
    r=least_squares(lambda p: p[0]+p[1]*x-np.log1p(np.exp(np.clip(p[2]*(1/(k*p[3])-1/(k*Tk)),-50,50)))-np.log(g.r.values),
        [np.log(np.median(g.r)),0.9,5.0,np.median(Tk)+5],bounds=([-40,0.05,0.3,285],[30,4,60,350]))
    E,Eh,Th=r.x[1],r.x[2],r.x[3]
    xa=1/(k*TREF)-1/(k*(a["T"].values.astype(float)+273.15))
    b=least_squares(lambda p:(p[0]+p[1]*xa)-np.log(a.K.values),[np.log(np.median(a.K)),0.5],bounds=([-80,0.0],[30,3]))
    return E,b.x[1],topt(E,Eh,Th),topt(E-b.x[1],Eh,Th)

out=[]
for st,s in w.groupby("strain"):
    g=s[s.growing]; 
    if g["T"].nunique()<5: continue
    base=fit(g,s); Ts=np.array(sorted(g["T"].unique())); BB=[]
    for _ in range(NB):
        pick=rng.choice(Ts,len(Ts),replace=True)
        gg=pd.concat([g[g["T"]==t].sample(len(g[g["T"]==t]),replace=True,random_state=rng.integers(1e9)) for t in pick])
        aa=pd.concat([s[s["T"]==t] for t in np.unique(pick)])
        if gg["T"].nunique()<4: continue
        try: BB.append(fit(gg,aa))
        except Exception: pass
    BB=np.array(BB)
    q=lambda i,p: np.nanpercentile(BB[:,i],p)
    to=BB[:,2]; to=to[np.isfinite(to)]
    out.append(dict(strain=st,n=len(s),nT=g["T"].nunique(),E=base[0],ER=base[1],
        Topt=base[2],To_lo=q(2,2.5),To_hi=q(2,97.5),Tcarbon=base[3],Tc_lo=q(3,2.5),Tc_hi=q(3,97.5),
        P_Topt_below_37=float(np.mean(to<37))))
o=pd.DataFrame(out); pd.set_option("display.width",220)
print(o.round(2).to_string(index=False))
o.to_csv(os.path.expanduser("~/mnt/Candidas/model_extra/cross_taxon/r2a_panel.csv"),index=False)
