"""PHASE 2B - is respiration without confirmed growth metabolically substantial?
PRESPECIFIED window, identical for every well, chosen BEFORE looking at species:
  start = 90 min   (the manuscript's 1.5 h thermal-equilibration exclusion)
  end   = 90 + 300 min = 390 min  (EARLY_H = 5.0 h, from scripts/20_fig4.py)
No window is selected per well, per taxon or per state.
CENSORING: if the trace reaches its oxygen floor inside the window, consumption is
RIGHT-CENSORED (the well could have consumed more) and flagged, never treated as low
metabolism. Floor = ymin + 0.02*(y_start - ymin), the legacy 2% depletion rule.
"""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W0,W1=90.0,390.0
tr=pd.read_csv(f"{C}/reprocess/blind_traces.csv.gz")
st=pd.read_csv(f"{C}/reprocess/unblinded_results.csv"); ISO=set(pd.read_csv(f"{C}/results/tables/otu_names.csv").OTU); st=st[st.OTU.isin(ISO)].copy()  # study isolates only
rows=[]
for w,g in tr.groupby("WELLID"):
    t=g.Time.values; y=g.Oxygen.values
    m=(t>=W0)&(t<=W1)
    if m.sum()<30: rows.append(dict(WELLID=w,cons=np.nan,rate=np.nan,censored=True,n=int(m.sum()))); continue
    tt,yy=t[m],y[m]
    post=t>=W0; ymin=float(np.min(y[post])); ystart=float(yy[0])
    floor=ymin+0.02*(ystart-ymin)
    cens=bool(np.any(yy<=floor))
    cons=float(yy[0]-yy.min())                      # mg/L consumed in the fixed window
    rate=cons/((tt[-1]-tt[0])/60.0)                 # mg/L/h
    rows.append(dict(WELLID=w,cons=cons,rate=rate,censored=cens,n=int(m.sum())))
cons=pd.DataFrame(rows).merge(st,on="WELLID")
NM={1:"growth",2:"resp_only",3:"ambig",4:"no_resp"}; cons["S"]=cons.state.map(NM)
cons["clade"]=np.where(cons.group.str.startswith("Clade"),"C. auris","relatives")
cons.to_csv(f"{C}/results/tables/phase2b_consumption.csv",index=False)
print(f"wells: {len(cons)} | censored (hit floor inside the fixed window): "
      f"{cons.censored.sum()} ({100*cons.censored.mean():.1f}%)")
print("\n=== censoring by state ==="); print(cons.groupby("S").censored.agg(['sum','size','mean']).round(3).to_string())
print("\n=== O2 consumption in the fixed 1.5-6.5 h window (mg/L), by state and temperature ===")
print(cons.groupby(["S","T"]).cons.median().unstack().round(2).to_string())
print("\n\n=== THE 2B COMPARISONS, at 40 °C ===")
r40=cons[(cons["T"]==40)&(cons.S=="resp_only")]
print(f"\nrespiration-only wells at 40 °C: n={len(r40)}  median consumption {r40.cons.median():.2f} mg/L "
      f"({100*r40.censored.mean():.0f}% censored)")
print(f"  by taxon:"); print(r40.groupby("group").cons.agg(['size','median']).round(2).to_string())
print("\n(i) vs confirmed-growing wells of the SAME ISOLATE at 40 °C:")
any_=False
for iso,s in r40.groupby("otu_name"):
    gw=cons[(cons["T"]==40)&(cons.S=="growth")&(cons.otu_name==iso)]
    if len(gw): any_=True; print(f"    {iso:12s} resp-only {s.cons.median():.2f}  vs growing {gw.cons.median():.2f}  "
                                 f"ratio {s.cons.median()/gw.cons.median():.2f}")
if not any_: print("    none - no isolate has both states at 40 °C")
print("\n(ii) vs the SAME ISOLATE's confirmed-growing wells at 38 °C:")
for iso,s in r40.groupby("otu_name"):
    g38=cons[(cons["T"]==38)&(cons.S=="growth")&(cons.otu_name==iso)]
    if len(g38):
        print(f"    {iso:12s} resp-only@40 {s.cons.median():.2f}  vs growing@38 {g38.cons.median():.2f}  "
              f"ratio {s.cons.median()/g38.cons.median():.2f}")
    else:
        print(f"    {iso:12s} resp-only@40 {s.cons.median():.2f}  vs growing@38  (none)")
a40=cons[(cons["T"]==40)&(cons.S=="growth")&(cons.clade=="C. auris")]
print(f"\n(iii) vs C. auris confirmed-growing wells at 40 °C: median {a40.cons.median():.2f} mg/L (n={len(a40)})")
print(f"      relatives' respiration-only wells reach "
      f"{100*r40.cons.median()/a40.cons.median():.0f}% of that")
