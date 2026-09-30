"""Task 1: reconcile the dataset denominator. Exclusion flow, per taxon and temperature."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv")
on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
ws=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv")[["T","OTU","Replicate","state"]]
cd=cd.merge(on,on="OTU",how="left")
print("=== raw table sizes ===")
print(f"fit_coefficients_wide rows : {len(cd)}   (23 isolates x 12 T x 5 wells = 1380 expected)")
print(f"fig4_well_states rows      : {len(ws)}")
ng=cd; ngw=ws.merge(on,on="OTU",how="left")
print(f"\nfit table {len(ng)} | well-states table {len(ngw)}")
iso=sorted(ng.otu_name.unique()); print(f"isolates: {len(iso)}  -> expected 20 x 12 x 5 = 1200")
# per isolate presence
pres=ng.groupby(["group","otu_name"]).size().rename("present")
print("\nwells present per isolate (expect 60):"); print(pres.to_string())
d=ng.merge(ws,on=["T","OTU","Replicate"],how="left")
flow=[]
for (g,T),s in d.groupby(["group","T"]):
    niso=s.otu_name.nunique()
    flow.append(dict(taxon=g,T=T,expected=niso*5,present=len(s),
        keep=int(s.keep.sum()),fit_valid=int(s.fit_valid.sum()),
        K_pos=int((s.K>0).sum()),r_pos=int((s.r>0).sum()),
        growing=int((s.state=="growing").sum()),resp_only=int((s.state=="respiring").sum()),
        inert=int((s.state=="inert").sum()),
        in_logKr=int((s.keep&s.fit_valid&(s.K>0)&(s.r>0)).sum())))
f=pd.DataFrame(flow)
tot=f[["expected","present","keep","fit_valid","K_pos","r_pos","growing","resp_only","inert","in_logKr"]].sum()
print("\n=== TOTALS ==="); print(tot.to_string())
print("\n=== per taxon ==="); print(f.groupby("taxon")[["expected","present","growing","resp_only","inert","in_logKr"]].sum().to_string())
print("\n=== exclusion flow by taxon x temperature (in_logKr) ===")
print(f.pivot(index="taxon",columns="T",values="in_logKr").to_string())
print("\n=== growing-well count by taxon x temperature (state-based) ===")
print(f.pivot(index="taxon",columns="T",values="growing").to_string())
print("\n=== RECONCILIATION ===")
print(f"  expected (20 x 12 x 5)            : 1200")
print(f"  physically present                : {tot.present}")
print(f"  manuscript states                 : 1155")
print(f"  used in log(K/r) analysis         : {tot.in_logKr}")
print(f"  state=growing (Fig 5 definition)  : {tot.growing}")
miss=d[d.state.isna()]
print(f"  rows in fit table with no state   : {len(miss)}")
bad=d[(d.keep&d.fit_valid&(d.K>0)&(d.r>0))&(d.state!="growing")]
print(f"  in log(K/r) but NOT state=growing : {len(bad)}  <- definition mismatch")
if len(bad): print(bad.groupby(["group","T"]).size().to_string())
f.to_csv(f"{C}/results/tables/exclusion_flow.csv",index=False)
