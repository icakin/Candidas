"""Scale-free comparison of the primary manual-window fits with their high-oxygen refits
(84_manual_window_hiO2_refit.py at cutoffs 0.5, 0.7, 0.8): model-free CUE optimum
(y = ln K - ln r - r t_s, isolate-first resampling, as in 75_modelfree_topt_n0.py),
Arrhenius activation energy of K over 22-38 C, and the 38->40 C change in K/(rN)."""
import os, sys, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); C=os.path.abspath(os.path.join(HERE,"..",".."))
names=pd.read_csv(f"{C}/results/tables/otu_names.csv")[["otu_name","group"]]
k=8.617e-5; rng=np.random.default_rng(3)
def topt_grid(df,col):
    m=df.groupby("T")[col].mean(); Ts=m.index.values.astype(float); v=m.values
    if len(Ts)<3: return np.nan
    i=int(np.argmin(v))
    if i==0 or i==len(v)-1: return Ts[i]
    x=Ts[i-1:i+2]; y=v[i-1:i+2]; a,b,c=np.polyfit(x,y,2); return -b/(2*a) if a>0 else Ts[i]
def boot(df,fn,nb=1000):
    isos=df.otu_name.unique(); bs=[]
    for _ in range(nb):
        pick=rng.choice(isos,len(isos),replace=True); parts=[]
        for iso in pick:
            si=df[df.otu_name==iso]
            for T,st in si.groupby("T"): parts.append(st.sample(len(st),replace=True))
        bs.append(fn(pd.concat(parts)))
    return np.array(bs,float)
rows=[]
for frac in (0.5,0.7,0.8):
    R=pd.read_csv(f"{HERE}/manual_refit_arms_{frac:.2f}.csv").merge(names,on="otu_name")
    R=R[R.identifiable&R.has_curvature&(R.r>0)&(R.K>0)]
    d=pd.read_csv(f"{C}/results/tables/derived_N0_R_results_with_carbon.csv")[["T","OTU","Replicate","fit_start_time"]]
    R=R.merge(d,on=["T","OTU","Replicate"]); R["y"]=np.log(R.K)-np.log(R.r)-R.r*R.fit_start_time
    for arm in ("full","hiO2"):
        if arm=="full" and frac!=0.5: continue
        a=R[R.arm==arm]
        for g,s in a.groupby("group"):
            if g not in ("Clade1","Clade2","Clade3","Clade4","para","Hae","Duo"): continue
            e=topt_grid(s,"y"); bs=boot(s,lambda x: topt_grid(x,"y"),300); bs=bs[np.isfinite(bs)]
            sub=s[s["T"]<=38]; x=1/(k*(sub["T"]+273.15)); EK=-np.polyfit(x,np.log(sub.K),1)[0]
            bE=boot(sub,lambda z: -np.polyfit(1/(k*(z["T"]+273.15)),np.log(z.K),1)[0],150)
            m=s.groupby("T").y.mean(); step=np.exp(m.get(40,np.nan)-m.get(38,np.nan))
            rows.append(dict(cutoff=frac,arm=arm,group=g,n=len(s),Topt_CUE=e,lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5),P37=np.mean(bs<37),E_K=EK,E_K_lo=np.percentile(bE,2.5),E_K_hi=np.percentile(bE,97.5),fold_38_40=step))
out=pd.DataFrame(rows); pd.set_option("display.width",220); print(out.round(2).to_string(index=False)); out.to_csv(f"{HERE}/manual_arms_comparison.csv",index=False)
