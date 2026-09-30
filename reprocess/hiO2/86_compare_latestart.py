import os, numpy as np, pandas as pd, sys
HERE=os.path.dirname(os.path.abspath(__file__)); C=os.path.abspath(os.path.join(HERE,"..",".."))
exec(open(f"{HERE}/85_compare_manual_arms.py").read().split("rows=[]")[0])   # reuse helpers
names=pd.read_csv(f"{C}/results/tables/otu_names.csv")[["otu_name","group"]]
d=pd.read_csv(f"{C}/results/tables/derived_N0_R_results_with_carbon.csv")[["T","OTU","Replicate","fit_start_time"]]
R=pd.read_csv(f"{HERE}/manual_refit_arms_0.70_latestart.csv").merge(names,on="otu_name").merge(d,on=["T","OTU","Replicate"])
R=R[R.identifiable&R.has_curvature&(R.r>0)&(R.K>0)]; R["y"]=np.log(R.K)-np.log(R.r)-R.r*R.fit_start_time
rows=[]
for arm in ("full_late","hiO2_late"):
    a=R[R.arm==arm]
    for g,s in a.groupby("group"):
        if g not in ("Clade1","Clade2","Clade3","Clade4","para","Hae","Duo"): continue
        e=topt_grid(s,"y"); bs=boot(s,lambda x: topt_grid(x,"y"),200); bs=bs[np.isfinite(bs)]
        sub=s[s["T"]<=38]; EK=-np.polyfit(1/(k*(sub["T"]+273.15)),np.log(sub.K),1)[0]
        rows.append(dict(arm=arm,group=g,n=len(s),Topt_CUE=e,lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5),P37=np.mean(bs<37),E_K=EK))
out=pd.DataFrame(rows); pd.set_option("display.width",200); print(out.round(2).to_string(index=False)); out.to_csv(f"{HERE}/manual_arms_latestart_comparison.csv",index=False)
