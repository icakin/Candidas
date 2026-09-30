"""Per-isolate thermal limits from the per-well likelihood model (reprocess/lik_model.py).
A well is scored as growing when the bootstrap p-value of the exponential-vs-linear
likelihood ratio is below P_GROW. Isolate-level loss temperature is the lowest measured
temperature at which fewer than half of the wells grow and this persists at every higher
temperature (the rule of 42_phase2a_transition.py); censored at >44 if never.
The grid interval [last T with >=50% growing, T_loss] is the resolution of the assay.
Also: log-rank comparison of C. auris vs relatives with an exact permutation null."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P_GROW=0.01
d=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv")
d["grow"]=d.p_boot<P_GROW
def first_persistent_below(s,thresh=0.5):
    ts=sorted(s.index)
    for i,T in enumerate(ts):
        if s[T]<thresh and all(s[u]<thresh for u in ts[i:]): return T
    return np.nan
rows=[]
for (g,iso),s in d.groupby(["group","otu_name"]):
    pg=s.groupby("T").grow.mean(); tl=first_persistent_below(pg)
    last_ok=max([T for T in pg.index if T<tl]) if np.isfinite(tl) else 44.0
    rows.append(dict(group=g,isolate=iso,T_growthloss=tl,growth_censored=np.isnan(tl),
                     T_last_growth=last_ok,p_growth_40=pg.get(40,np.nan),p_growth_42=pg.get(42,np.nan),
                     p_growth_44=pg.get(44,np.nan),n_wells=len(s)))
t=pd.DataFrame(rows); t["clade"]=np.where(t.group.str.startswith("Clade"),"C. auris","relatives")
t.to_csv(f"{C}/results/tables/lik_transition.csv",index=False)
pd.set_option("display.width",200); print(t.to_string(index=False))
def logrank(time,event,grp):
    ts=np.sort(np.unique(time[event==1])); O_E=0.0; V=0.0
    for tt in ts:
        at=time>=tt; n=at.sum(); n1=(at&(grp==1)).sum()
        dd=((time==tt)&(event==1)).sum(); d1=((time==tt)&(event==1)&(grp==1)).sum()
        if n<2: continue
        O_E+=d1-dd*n1/n; V+=dd*(n1/n)*(1-n1/n)*((n-dd)/(n-1))
    return (O_E**2/V) if V>0 else np.nan
rng=np.random.default_rng(2026)
for excl,tag in ((None,"all isolates"),("Hae_1724","excluding Hae_1724")):
    x=t if excl is None else t[t.isolate!=excl]
    time=np.where(x.growth_censored,44.0,x.T_growthloss.values); event=(~x.growth_censored).astype(int).values
    grp=(x.clade=="C. auris").astype(int).values; obs=logrank(time,event,grp)
    null=np.array([logrank(time,event,rng.permutation(grp)) for _ in range(20000)]); null=null[np.isfinite(null)]
    p=(np.sum(null>=obs)+1)/(len(null)+1)
    for cl in ("C. auris","relatives"):
        y=x[x.clade==cl]; o=y[~y.growth_censored]
        print(f"  {tag}: {cl:10s} n={len(y)} observed loss n={len(o)} median {o.T_growthloss.median() if len(o) else float('nan')} censored {int(y.growth_censored.sum())}")
    print(f"  {tag}: log-rank chi2 = {obs:.2f}, exact permutation p = {p:.4f}\n")
