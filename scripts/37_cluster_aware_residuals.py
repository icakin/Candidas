"""Task 2: isolate-clustered re-test of taxon x temperature residual structure,
after correcting the well-inclusion definition found in the denominator audit."""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
rng=np.random.default_rng(17); NPERM=5000
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv"); on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
ws=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv")[["T","OTU","Replicate","state"]]
d=cd.merge(on,on="OTU").merge(ws,on=["T","OTU","Replicate"],how="left")
d=d[(d.keep)&(d.fit_valid)&(d.K>0)&(d.r>0)].copy()
d["y"]=np.log(d.K/d.r); d["TK"]=d["T"]+273.15
OLD=d.copy(); NEW=d[d.state=="growing"].copy()
print(f"previous inclusion (fit-valid)      : {len(OLD)} wells")
print(f"corrected inclusion (state=growing) : {len(NEW)} wells   [drops {len(OLD)-len(NEW)}]")
# common support: temperatures where EVERY taxon keeps >=3 growing wells
g=NEW.groupby(["group","T"]).size().unstack(fill_value=0)
ok=[T for T in g.columns if (g[T]>=3).all()]
print(f"\ngrowing wells per taxon x T:\n{g.to_string()}")
print(f"\ncommon support (all taxa >=3 growing wells): {min(ok)}-{max(ok)} C")
def build(sub):
    TAX=sorted(sub.group.unique()); ti={t:i for i,t in enumerate(TAX)}; n=len(TAX)
    def pr(p,df,m):
        k=df.group.map(ti).values; TK=df.TK.values; a=p[:n]
        if m==1: dE=np.full(n,p[n]); Eh=np.full(n,p[n+1]); Th=p[n+2:2*n+2]
        else:    dE=p[n:2*n];        Eh=p[2*n:3*n];        Th=p[3*n:4*n]
        u=np.clip(Eh[k]*CINV*(1/Th[k]-1/TK),-50,50)
        return a[k]+dE[k]*CINV*(1/TREF-1/TK)+np.log1p(np.exp(u))
    a0=[np.log(sub[sub.group==t].K.median()/sub[sub.group==t].r.median()) for t in TAX]
    st={1:np.array(a0+[-.4,3.]+[309.]*n),3:np.array(a0+[-.4]*n+[3.]*n+[309.]*n)}
    bd={1:([-20]*n+[-4,.3]+[295]*n,[20]*n+[4,40]+[330]*n),
        3:([-20]*n+[-4]*n+[.3]*n+[295]*n,[20]*n+[4]*n+[40]*n+[330]*n)}
    return pr,lambda df,m: least_squares(lambda p:pr(p,df,m)-df.y.values,st[m],bounds=bd[m],max_nfev=20000)
def cluster_test(sub,label):
    pr,fit=build(sub); f=fit(sub,1); sub=sub.copy(); sub["res"]=sub.y.values-pr(f.x,sub,1)
    cell=sub.groupby(["otu_name","group","T"]).res.mean().reset_index()   # isolate x T unit
    cell["res_c"]=cell.res-cell.groupby("T").res.transform("mean")        # remove T main effect
    def stat(lab):
        m=cell.otu_name.map(lab)
        return cell.assign(g=m).groupby(["g","T"]).res_c.mean().pow(2).sum()
    true_lab=dict(zip(cell.otu_name,cell.group)); obs=stat(true_lab)
    isos=list(true_lab); labs=[true_lab[i] for i in isos]
    null=[]
    for _ in range(NPERM):
        p=rng.permutation(labs); null.append(stat(dict(zip(isos,p))))
    null=np.array(null); pval=(np.sum(null>=obs)+1)/(NPERM+1)
    print(f"\n  [{label}] isolate-clustered permutation of taxon labels ({len(isos)} isolates)")
    print(f"    observed taxon x T structure = {obs:.4f}; null median {np.median(null):.4f} "
          f"(95th pct {np.percentile(null,95):.4f})   p = {pval:.4f}")
    return cell,pval
print("\n=== (2a) previous inclusion, full range ===");  c_old,_=cluster_test(OLD,"fit-valid, 22-44")
print("\n=== (2b) corrected inclusion, full range ==="); c_new,_=cluster_test(NEW,"growing, 22-44")
print("\n=== (2c) corrected inclusion, common support only ===")
c_cs,_=cluster_test(NEW[NEW["T"]<=max(ok)],f"growing, {min(ok)}-{max(ok)}")
print("\n=== per-isolate residual trajectory, 36-44 C (corrected inclusion, M1) ===")
piv=c_new[c_new["T"]>=36].pivot_table(index=["group","otu_name"],columns="T",values="res")
print(piv.round(2).to_string())
c_new.to_csv(f"{C}/results/tables/m1_isolate_residuals.csv",index=False)
