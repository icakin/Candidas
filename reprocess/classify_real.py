"""Apply the frozen classifier to the blinded traces. No labels are read."""
import pandas as pd, numpy as np, sys, hashlib; sys.path.insert(0,'.')
import classifier as C
tr=pd.read_csv("blind_traces.csv.gz")
rows=[]
for i,(w,g) in enumerate(tr.groupby("WELLID")):
    o=C.classify(g.Time.values,g.Oxygen.values,seed=int(hashlib.md5(w.encode()).hexdigest()[:8],16))
    rows.append(dict(WELLID=w,**o))
    if (i+1)%200==0: print(f"  {i+1} wells",flush=True)
d=pd.DataFrame(rows)
d.to_csv("blind_states.csv",index=False)
print(f"\nwells classified: {len(d)}")
print("\n=== blind state counts ===")
NAMES={1:"1 growth+respiration",2:"2 respiration, no detectable growth",
       3:"3 ambiguous",4:"4 no detectable respiration"}
vc=d.state.value_counts().sort_index()
for s,n in vc.items(): print(f"  {NAMES[s]:40s} {n:5d}  ({100*n/len(d):5.1f}%)")
print("\n=== reason codes ==="); print(d.reason.value_counts().to_string())
print("\n=== quality summary by state ===")
print(d.groupby("state")[["drawdown","dur_h","eff_n","rho","rt","dlpd"]].median().round(3).to_string())
