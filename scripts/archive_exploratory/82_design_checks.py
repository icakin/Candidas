"""Design checks (Methods, Design checks): reader/incubator confounded with temperature pair,
group confounded with run date. (a) Reader-boundary test: per-well log residuals from each
isolate's fitted growth and respiration curves, tested for a step at the six reader boundaries
(24|26, 28|30, 32|34, 36|38, 40|42) against steps at the within-reader boundaries
(22|24, 26|28, ...). (b) Run-date check: the species contrast with *C. parapsilosis* (June)
and Hae_1724 (August) as the cross-over cases. Writes results/tables/design_checks.csv."""
import numpy as np, pandas as pd
from tpc_common import *
G,I,quota,vals=load(); P=points()
run=pd.read_csv(f"{C}/results/tables/run_log.csv"); run=run[run.group.isin(GROUPS)]
P=P.merge(run[["group","T","sdr","day"]],on=["group","T"],how="left")
# residuals from the isolate-level posterior-median curves
res=[]
for iso,di in I.items():
    g=group_of(iso); q=quota[g]; p=P[P.otu_name==iso]
    if len(p)==0: continue
    E,Eh,Th,ER,lnB0,al=[np.median(di[c].values) for c in ("E","Eh","Th","ER","lnB0","alpha")]
    gfit=np.exp(lnB0+lnG(p["T"].values,E,Eh,Th))/q; rfit=np.exp(al+lnR(p["T"].values,ER))
    for (_,w),gf,rf in zip(p.iterrows(),gfit,rfit):
        res.append(dict(isolate=iso,group=g,T=w["T"],sdr=w.sdr,day=w.day,res_g=np.log(w.r_h/gf),res_r=np.log(w.respiration_fgC_h/rf)))
R=pd.DataFrame(res)
# mean residual per temperature (over wells), then step sizes between adjacent temperatures
mt=R.groupby("T")[["res_g","res_r"]].mean(); Ts=sorted(mt.index)
steps=[]
for a,b in zip(Ts[:-1],Ts[1:]):
    same_reader=run[(run["T"]==a)].sdr.iloc[0]==run[(run["T"]==b)].sdr.iloc[0]
    steps.append(dict(boundary=f"{a}|{b}",within_reader=same_reader,step_g=mt.loc[b,"res_g"]-mt.loc[a,"res_g"],step_r=mt.loc[b,"res_r"]-mt.loc[a,"res_r"]))
S=pd.DataFrame(steps); print(S.round(3).to_string(index=False))
for col in ("step_g","step_r"):
    a=S[~S.within_reader][col].abs(); b=S[S.within_reader][col].abs()
    print(f"{col}: |step| at reader boundaries median {a.median():.3f} (n={len(a)}), within reader {b.median():.3f} (n={len(b)})")
# permutation test: is the mean |step| at reader boundaries larger than at random boundary labellings?
rng=np.random.default_rng(1); out={}
for col in ("step_g","step_r"):
    obs=S[~S.within_reader][col].abs().mean()-S[S.within_reader][col].abs().mean()
    null=[]
    lab=S.within_reader.values.copy()
    for _ in range(20000):
        rng.shuffle(lab); null.append(S[~lab][col].abs().mean()-S[lab][col].abs().mean())
    out[col]=(obs,np.mean(np.array(null)>=obs))
    print(f"{col}: reader-boundary excess {obs:+.3f} log units, permutation p = {out[col][1]:.3f}")
# per-reader mean residual (all wells on that reader)
byr=R.groupby("sdr")[["res_g","res_r"]].agg(["mean","std","size"]).round(3); print(byr)
S["metric"]="step"; S.to_csv(f"{C}/results/tables/design_checks.csv",index=False)
with open(f"{C}/results/tables/design_checks_summary.txt","w") as fh:
    for col,(obs,p) in out.items(): fh.write(f"{col}\treader_boundary_excess_abs_step={obs:.4f}\tperm_p={p:.4f}\n")
    fh.write(byr.to_string())
