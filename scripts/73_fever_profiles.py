"""Idea 4: net growth under host temperature profiles, from the fitted isolate curves.
r(T) in h^-1 = exp(lnB0 + lnG(T)) / quota (as fig_common.R). Net daily gain
ln(N24/N0) = integral of r(T(t)) dt over 24 h, assuming instantaneous thermal response and
no death; expressed as doublings per day. Profiles: constant 37; constant 40; 12 h at 40 and
12 h at 37 (sustained fever with defervescence); 4 h at 40 in 24 h (spiking fever); constant
38.5. Posterior medians and 95% CrI over draws, per isolate; species contrast on the
difference between profiles."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv"); pp=pp[~pp.Isolate.str.startswith("GROUP_")]
q=pd.read_csv(f"{C}/results/tables/group_quota.csv").set_index("Group").quota_fgC
def lnG(TC,E,Eh,Th):
    TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK); return E*CINV*(1/TREF-1/TK)-np.log1p(np.exp(np.minimum(u,700)))
def r_h(TC,g): return np.exp(g.lnB0.values+lnG(TC,g.E.values,g.Eh.values,g.Th.values))/q[g.Group.iloc[0]]
PROFILES={"37 constant":[(37,24)],"38.5 constant":[(38.5,24)],"40 constant":[(40,24)],
          "40 for 12 h, 37 for 12 h":[(40,12),(37,12)],"40 for 4 h, 37 for 20 h":[(40,4),(37,20)]}
rows=[]
for iso,g in pp.groupby("Isolate"):
    row=dict(isolate=iso,group=g.Group.iloc[0],auris=g.Group.iloc[0].startswith("Clade"))
    for name,segs in PROFILES.items():
        gain=sum(r_h(T,g)*h for T,h in segs)/np.log(2)   # doublings per day
        row[name]=np.median(gain); row[name+"_lo"]=np.percentile(gain,2.5); row[name+"_hi"]=np.percentile(gain,97.5)
    rows.append(row)
R=pd.DataFrame(rows); R.to_csv(f"{C}/results/tables/fever_profiles.csv",index=False)
pd.set_option("display.width",220)
print("doublings per day (posterior median):"); print(R[["isolate"]+list(PROFILES)].round(1).to_string(index=False))
print("\nby species group (median of isolate medians):"); print(R.groupby("auris")[list(PROFILES)].median().round(1).to_string())
R["fever_penalty_40"]=R["40 constant"]/R["37 constant"]; R["fever_penalty_12h"]=R["40 for 12 h, 37 for 12 h"]/R["37 constant"]
print("\ngrowth retained under fever relative to 37 C (median across isolates): 40 constant auris %.2f rel %.2f | 12h/12h auris %.2f rel %.2f"%(
    R[R.auris].fever_penalty_40.median(),R[~R.auris].fever_penalty_40.median(),R[R.auris].fever_penalty_12h.median(),R[~R.auris].fever_penalty_12h.median()))
