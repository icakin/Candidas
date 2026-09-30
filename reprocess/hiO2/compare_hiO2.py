"""Compare the high-oxygen refit (traces truncated at 50% of the window-start oxygen) with the
full-window likelihood model: growth calls, per-well r and K, the model-free CUE optimum
(y = ln K - ln r - r t_s, isolate-first resampling) and the 38->40 C change in K/(r N)."""
import pandas as pd, numpy as np, sys
FRAC=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
import os; HERE=os.path.dirname(os.path.abspath(__file__)); full=pd.read_csv(os.path.join(HERE,"..","lik_results_labelled.csv")); hi=pd.read_csv(os.path.join(HERE,f"lik_hiO2_{FRAC:.2f}.csv"))
full=full[full.group.isin(["Clade1","Clade2","Clade3","Clade4","para","Hae","Duo"])]
d=full.merge(hi,on="WELLID",suffixes=("","_hi")); print("wells",len(d))
d["grow"]=(d.p_boot<0.01)&(d.drawdown>=1); d["grow_hi"]=d.identifiable&(d.p_boot_hi<0.01)&(d.drawdown>=1)
print("identifiable after truncation:",d.identifiable.sum(),"of",len(d)); 
ct=pd.crosstab(d.grow,d.grow_hi); print(ct)
u=d[~d.identifiable]; print("unidentifiable by full-window state: growth",int(u.grow.sum()),"non-growth",int((~u.grow).sum()))
print("unidentifiable by T:",u.groupby("T").size().to_dict())
# isolate-level counts at 40
for T in (40,42):
    a=d[(d["T"]==T)]; print(f"T={T}: full growth {int(a.grow.sum())}/{len(a)}  hiO2 growth {int(a.grow_hi.sum())}/{len(a)}  (hiO2 among identifiable {int(a[a.identifiable].grow_hi.sum())}/{int(a.identifiable.sum())})")
    for sp,lab in ((True,"auris"),(False,"relatives")):
        b=a[a.auris==sp] if "auris" in a else a[a.group.str.startswith("Clade")==sp]
        print(f"   {lab}: full {int(b.grow.sum())}/{len(b)}, hiO2 {int(b.grow_hi.sum())}/{len(b)}")
# per-well r and K agreement among wells growing in both
g=d[d.grow&d.grow_hi]
print("r: median ratio hi/full",np.median(g.r_hat_hi/g.r_hat).round(3),"IQR",np.percentile(g.r_hat_hi/g.r_hat,[25,75]).round(3),"Spearman",g[["r_hat","r_hat_hi"]].corr("spearman").iloc[0,1].round(3))
print("K: median ratio hi/full",np.median(g.K_hat_hi/g.K_hat).round(3),"IQR",np.percentile(g.K_hat_hi/g.K_hat,[25,75]).round(3),"Spearman",g[["K_hat","K_hat_hi"]].corr("spearman").iloc[0,1].round(3))
# model-free CUE optimum both arms
rng=np.random.default_rng(3); ts=1.5
def topt_grid(df,col):
    m=df.groupby("T")[col].mean(); Ts=m.index.values.astype(float); v=m.values
    if len(Ts)<3: return np.nan
    k=int(np.argmin(v))
    if k==0 or k==len(v)-1: return Ts[k]
    x=Ts[k-1:k+2]; y=v[k-1:k+2]; a,b,c=np.polyfit(x,y,2); return -b/(2*a) if a>0 else Ts[k]
def est(df,col,nb=1000):
    isos=df.otu_name.unique(); e=topt_grid(df,col); bs=[]
    for _ in range(nb):
        pick=rng.choice(isos,len(isos),replace=True); parts=[]
        for iso in pick:
            si=df[df.otu_name==iso]
            for T,st in si.groupby("T"): parts.append(st.sample(len(st),replace=True))
        bs.append(topt_grid(pd.concat(parts),col))
    bs=np.array(bs); bs=bs[np.isfinite(bs)]; return e,np.percentile(bs,2.5),np.percentile(bs,97.5),np.mean(bs<37)
rows=[]
for grp,s in d.groupby("group"):
    a=s[s.grow].copy(); a["y"]=np.log(a.K_hat)-np.log(a.r_hat)-a.r_hat*ts
    b=s[s.grow_hi].copy(); b["y"]=np.log(b.K_hat_hi)-np.log(b.r_hat_hi)-b.r_hat_hi*ts
    ef=est(a,"y"); eh=est(b,"y")
    # 38->40 change in ln(K/(rN))
    def step(x):
        m=x.groupby("T").y.mean(); return (m.get(40,np.nan)-m.get(38,np.nan))
    rows.append(dict(group=grp,n_full=len(a),n_hi=len(b),Topt_full=ef[0],lo_f=ef[1],hi_f=ef[2],P37_full=ef[3],Topt_hi=eh[0],lo_h=eh[1],hi_h=eh[2],P37_hi=eh[3],
                     fold38_40_full=np.exp(step(a)),fold38_40_hi=np.exp(step(b))))
R=pd.DataFrame(rows); pd.set_option("display.width",220); print(R.round(2).to_string(index=False))
R.to_csv(os.path.join(HERE,f"hiO2_comparison_{FRAC:.2f}.csv"),index=False)
d.to_csv(os.path.join(HERE,f"lik_hiO2_{FRAC:.2f}_labelled.csv"),index=False)
