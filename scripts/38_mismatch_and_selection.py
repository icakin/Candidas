"""Task 3: are the 125 definition-mismatch wells' r values real or detection-floor?
Plus: which isolates actually carry the high-temperature growing wells that the
Bayesian fits condition on (selection / survivor structure)."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv"); on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
ws=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv")
d=cd.merge(on,on="OTU").merge(ws[["T","OTU","Replicate","state","drawdown","window_h","slope_early"]],
                              on=["T","OTU","Replicate"],how="left")
d=d[(d.keep)&(d.fit_valid)&(d.K>0)&(d.r>0)].copy()
d["rh"]=d.r*60
d["rt"]=d.r*60*d.window_h   # e-folds over the fitting window
mis=d[d.state!="growing"]; gro=d[d.state=="growing"]
print(f"mismatch wells: {len(mis)} | growing: {len(gro)}\n")
print("=== (3) mismatch wells: r (h^-1), rt, drawdown, window ===")
print(mis.groupby(["group","T"]).agg(n=("rh","size"),r_med=("rh","median"),r_max=("rh","max"),
      rt_med=("rt","median"),dd_med=("drawdown","median"),win_med=("window_h","median")).round(3).to_string())
print("\n=== same cells, ACCEPTED growing wells, for comparison ===")
cmp=gro[gro.set_index(["group","T"]).index.isin(mis.set_index(["group","T"]).index)]
print(cmp.groupby(["group","T"]).agg(n=("rh","size"),r_med=("rh","median"),rt_med=("rt","median"),
      dd_med=("drawdown","median")).round(3).to_string())
print("\n=== rt distribution (rt = r * window; <1 means less than one e-fold over the whole window) ===")
for nm,s in (("mismatch",mis),("growing",gro)):
    q=np.nanpercentile(s.rt,[5,25,50,75,95])
    print(f"  {nm:9s} rt  5/25/50/75/95 pct: "+"  ".join(f"{v:.2f}" for v in q)+
          f"   frac rt<0.5: {np.mean(s.rt<0.5):.2f}  frac rt<1: {np.mean(s.rt<1):.2f}")
print("\n  -> earlier pipeline rule (rt > 0.5): mismatch wells passing = "
      f"{int((mis.rt>0.5).sum())}/{len(mis)}; growing wells passing = {int((gro.rt>0.5).sum())}/{len(gro)}")
print("\n=== r relative to the same taxon's r at its permissive optimum ===")
for g_,s in mis.groupby("group"):
    ref=gro[(gro.group==g_)&(gro["T"].between(30,34))].rh.median()
    print(f"  {g_:7s} mismatch r median {s.rh.median():.4f} h-1  vs permissive {ref:.3f} h-1  = {100*s.rh.median()/ref:.1f}%")
print("\n=== (5) SELECTION: which isolates carry the growing wells at 38-44 C ===")
hi=gro[gro["T"]>=38]
tab=hi.groupby(["group","T","otu_name"]).size().unstack(fill_value=0)
print(tab.to_string())
print("\n  taxa whose high-T growing wells come from a SINGLE isolate:")
for (g_,T),s in hi.groupby(["group","T"]):
    if s.otu_name.nunique()==1:
        print(f"    {g_:7s} {T} C : all {len(s)} wells from {s.otu_name.iloc[0]}")
mis[["T","OTU","Replicate","otu_name","group","r","K","rt","drawdown","window_h","state"]].to_csv(
    f"{C}/results/tables/mismatch_wells.csv",index=False)
