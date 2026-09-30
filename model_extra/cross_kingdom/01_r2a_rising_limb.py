"""R2A soil OTUs and E. coli, same instrument, same per-well model as the manuscript.
Per strain: Arrhenius activation energy of r (E_G) and of K (E_R) on the temperatures
where the strain grows, residual bootstrap; and the state of every well (growing /
respiring with growth below detection / inert), same criteria as manuscript Fig. 5."""
import glob, re, os, numpy as np, pandas as pd
from scipy.optimize import least_squares
k=8.617e-5; TREF=293.15; rng=np.random.default_rng(7); NB=1000
B=os.path.expanduser("~/mnt/Candidas/data/Bacteria"); OUT=os.path.expanduser("~/mnt/Candidas/model_extra/cross_kingdom")
R_GROW=0.05  # h-1: wells fall either below 0.02 or above 0.08 h-1, the threshold sits in that gap
def fit_well(t,y):
    i0=int(np.argmax(y[:max(5,len(y)//10)])); t=t[i0:]; y=y[i0:]
    if len(t)<30: return None
    ymin=y.min(); thr=ymin+0.02*(y[0]-ymin); b=np.where(y<=thr)[0]; i1=b[0] if len(b) else len(y)-1
    t=t[:i1+1]-t[0]; y=y[:i1+1]; dd=y[0]-y[-1]
    if len(t)<30: return dict(r=np.nan,K=np.nan,r2=np.nan,dd=dd)
    best=None
    for r0 in (0.0005,0.003,0.01):
        try:
            r=least_squares(lambda p:(p[0]-(p[1]/p[2])*np.expm1(p[2]*t))-y,[y[0],max(dd/t[-1],1e-6),r0],bounds=([-5,1e-9,1e-6],[25,5,0.5]))
            if best is None or r.cost<best.cost: best=r
        except Exception: pass
    if best is None: return dict(r=np.nan,K=np.nan,r2=np.nan,dd=dd)
    O,K,r=best.x; r2=1-np.sum(best.fun**2)/np.sum((y-y.mean())**2)
    # linear slope over first 5 h for the non-growing state (mg/L/h)
    m5=t<=300; slope=-np.polyfit(t[m5],y[m5],1)[0]*60 if m5.sum()>10 else np.nan
    m1=(t>=120)&(t<=180); s1=-np.polyfit(t[m1],y[m1],1)[0]*60 if m1.sum()>10 else np.nan   # observed slope over hours 2-3 (after equilibration), mg/L/h
    rh=r*60; pred1=K*60*np.exp(2*rh)*np.expm1(rh)/rh if rh>1e-6 else K*60                  # model's own mean slope over hours 2-3
    return dict(r=rh,K=K*60,r2=r2,dd=dd,slope=slope,lag_ratio=s1/pred1 if pred1>0 else np.nan)
def load(pattern,rx,names):
    rows=[]
    for f in sorted(glob.glob(pattern),key=lambda x:int(re.search(rx,x).group(1))):
        T=int(re.search(rx,f).group(1)); d=pd.read_csv(f); t=d.Time.values.astype(float)
        for c in [c for c in d.columns if c.startswith("OTU")]:
            y=pd.to_numeric(d[c],errors="coerce").values; m=np.isfinite(y)&np.isfinite(t)
            if m.sum()<40: continue
            o=fit_well(t[m],y[m])
            if o: rows.append(dict(T=T,strain=names[c.split("_")[0]],well=c,**o))
    return pd.DataFrame(rows)
parts=[]
for folder in ("7_13_14_15","2_6_9_11","1_4_8_20"):
    ids=folder.split("_")
    parts.append(load(f"{B}/{folder}/*.csv", r"R2A_(\d+)", {f"OTU{i+1}":f"OTU{ids[i]}" for i in range(4)}))
parts.append(load(f"{B}/Ecoli/*.csv", r"Ecoli_(\d+)_", {"OTU1":"E. coli A","OTU2":"E. coli B"}))
w=pd.concat(parts,ignore_index=True)
# a well whose observed first-hour consumption is >3x what its own fit implies has a lag followed by
# late growth: the exponential model does not describe its start, and it is counted as respiring
# (growth below detection, late recovery), as the manuscript does for late-recovering wells.
w["late_recovery"]=(w.r>R_GROW)&(w.lag_ratio>3)
w["state"]=np.where((w.r>R_GROW)&(w.r2>=0.9)&(~w.late_recovery),"growing",np.where(w.dd>=2.0,"respiring","inert"))
print("late-recovery wells:"); print(w[w.late_recovery][["strain","T","well","r","K","lag_ratio"]].round(3).to_string())
w.to_csv(f"{OUT}/r2a_wells.csv",index=False)
print(w.groupby(["strain","T"]).state.apply(lambda s:f"{(s=='growing').sum()}g/{(s=='respiring').sum()}r/{(s=='inert').sum()}i").unstack(fill_value="-").to_string()); print()
def arr(T,y):
    x=1/(k*TREF)-1/(k*(T+273.15)); A=np.c_[np.ones_like(x),x]; b,res,_,_=np.linalg.lstsq(A,np.log(y),rcond=None); return b,A
out=[]
for st,s in w.groupby("strain"):
    g=s[s.state=="growing"]; nt=g.groupby("T").size(); Tg=nt[nt>=3].index.values
    g=g[g["T"].isin(Tg)]
    if len(Tg)<3: out.append(dict(strain=st,nT=len(Tg),note="fewer than 3 growing temperatures")); continue
    (bG,A)=arr(g["T"].values,g.r.values); (bR,_)=arr(g["T"].values,g.K.values)
    resG=np.log(g.r.values)-A@bG; resR=np.log(g.K.values)-A@bR; EG=[];ER=[]
    for _ in range(NB):
        i=rng.integers(0,len(resG),len(resG)); j=rng.integers(0,len(resR),len(resR))
        EG.append(np.linalg.lstsq(A,A@bG+resG[i],rcond=None)[0][1]); ER.append(np.linalg.lstsq(A,A@bR+resR[j],rcond=None)[0][1])
    EG=np.array(EG); ER=np.array(ER)
    above=s[s["T"]>Tg.max()]
    out.append(dict(strain=st,nT=len(Tg),Tmin=Tg.min(),Tmax=Tg.max(),n=len(g),
        E_G=bG[1],E_G_lo=np.quantile(EG,.025),E_G_hi=np.quantile(EG,.975),
        E_R=bR[1],E_R_lo=np.quantile(ER,.025),E_R_hi=np.quantile(ER,.975),
        P_EG_gt_ER=np.mean(EG-ER>0),
        wells_above_limit=len(above),respiring_above_limit=int((above.state=="respiring").sum()),
        inert_above_limit=int((above.state=="inert").sum()),
        slope_respiring_above=above[above.state=="respiring"].slope.median() if (above.state=="respiring").any() else np.nan,
        slope_growing_Tmax=g[g["T"]==Tg.max()].K.median()))
o=pd.DataFrame(out); o.to_csv(f"{OUT}/r2a_strain_fits.csv",index=False)
print(o.round(3).to_string())
