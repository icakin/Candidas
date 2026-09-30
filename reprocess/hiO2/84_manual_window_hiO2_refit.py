"""Apples-to-apples high-oxygen arm for the PRIMARY (manual-window) pipeline.
For every well that enters the primary fits, the oxygen-dynamics model
    O2(t) = O2_0 + (K/r)(1 - exp(r t)),  t from the manual window start,
is refitted (a) on the full manual window, reproducing the pipeline, and (b) on the window
truncated at the first point where measured oxygen falls below 50% of its measured value at
the window start. Per-cell growth and respiration are then derived with the pipeline's own
formulas (N0 = N_inoc exp(r delta); R per cell = K/N0), giving two datasets in exactly the
form 09_bayesian_models.R consumes. 41_refit_hiO2_brms.R fits the hierarchical models to both.
Output: reprocess/hiO2/manual_refit_arms.csv"""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
HERE=os.path.dirname(os.path.abspath(__file__)); C=os.path.abspath(os.path.join(HERE,"..",".."))
import sys
FRAC=float(sys.argv[1]) if len(sys.argv)>1 else 0.5; MIN_PTS=20; MIN_MIN=60.0
MG_TO_FG=1e12; O2_TO_C=12.011/31.998; RQ=1.0; MIN_TO_H=60.0
ox=pd.read_csv(f"{C}/results/tables/Oxygen_Data_Filtered.csv")
d=pd.read_csv(f"{C}/results/tables/derived_N0_R_results_with_carbon.csv")
names=pd.read_csv(f"{C}/results/tables/otu_names.csv")[["OTU","otu_name","group"]] if "OTU" in pd.read_csv(f"{C}/results/tables/otu_names.csv",nrows=1).columns else None
d=d[d.keep&d.fit_valid].copy()
G={k:g for k,g in ox.groupby(["T","OTU","Replicate"])}
def model(p,t): O0,K,r=p; return O0+(K/r)*(1-np.exp(r*t))
def fit(t,y):
    ys=np.ptp(y); O0s=np.median(y[:3]); Ks=max(1e-6,-np.polyfit(t[:max(5,len(t)//3)],y[:max(5,len(t)//3)],1)[0])
    best=None
    for r0 in (1e-4,1e-3,5e-3,2e-2):
        try:
            f=least_squares(lambda p: model(p,t)-y,[O0s,Ks,r0],bounds=([y.min()-0.5*ys,1e-10,1e-6],[y.max()+0.5*ys,max(1.0,Ks*50),0.15]),method="trf",max_nfev=2000)
        except Exception: continue
        if best is None or f.cost<best.cost: best=f
    if best is None: return None
    J=best.jac; res=best.fun; dof=max(1,len(y)-3); s2=np.sum(res**2)/dof
    try: cov=np.linalg.inv(J.T@J)*s2; se=np.sqrt(np.clip(np.diag(cov),0,None))
    except np.linalg.LinAlgError: se=np.full(3,np.nan)
    return dict(O2_0=best.x[0],K=best.x[1],r=best.x[2],K_se=se[1],r_se=se[2],rss=np.sum(res**2))
rows=[]
for _,w in d.iterrows():
    key=(w["T"],w.OTU,w.Replicate)
    if key not in G: continue
    g=G[key].sort_values("Time"); m=(g.Time>=w.fit_start_time)&(g.Time<=w.fit_end_time); g=g[m]
    if len(g)<MIN_PTS: continue
    t=g.Time.values-w.fit_start_time; y=g.Oxygen.values; O2start=np.median(y[:3])
    cut=np.argmax(y<FRAC*O2start) if (y<FRAC*O2start).any() else len(y)
    for arm,(tt,yy) in (("full",(t,y)),("hiO2",(t[:cut],y[:cut]))):
        ok=len(tt)>=MIN_PTS and (tt[-1]-tt[0])>=MIN_MIN
        f=fit(tt,yy) if ok else None
        row=dict(arm=arm,T=w["T"],OTU=w.OTU,Replicate=w.Replicate,otu_name=w.otu_name,cell_carbon_fg=w.cell_carbon_fg,
                 has_curvature=bool(w.has_curvature),O2_start_meas=O2start,n_pts=len(tt),dur_min=(tt[-1]-tt[0]) if len(tt) else 0.0,identifiable=f is not None,
                 pipeline_r=w.r,pipeline_K=w.K)
        if f:
            N0=w.N_inoculation_cells_per_L*np.exp(f["r"]*w.delta_Ninoc_to_N0_min)
            row.update(f, N0_cells_per_L=N0, growth_fgC_h=f["r"]*w.cell_carbon_fg*MIN_TO_H,
                       respiration_fgC_h=(f["K"]/N0)*MG_TO_FG*O2_TO_C*RQ*MIN_TO_H,
                       se_log_r=f["r_se"]/f["r"] if f["r"]>0 else np.nan,
                       se_log_resp=np.sqrt((f["K_se"]/f["K"])**2+(w.delta_Ninoc_to_N0_min*f["r_se"])**2) if f["K"]>0 else np.nan)
        rows.append(row)
R=pd.DataFrame(rows); R.to_csv(f"{HERE}/manual_refit_arms_{FRAC:.2f}.csv",index=False)
full=R[R.arm=="full"]; hi=R[R.arm=="hiO2"]
print("wells",len(full),"| full identifiable",full.identifiable.sum(),"| hiO2 identifiable",hi.identifiable.sum())
print("full vs pipeline r: median ratio",np.nanmedian(full.r/full.pipeline_r).round(3),"K:",np.nanmedian(full.K/full.pipeline_K).round(3))
mm=full.merge(hi,on=["T","OTU","Replicate"],suffixes=("_f","_h")); mm=mm[mm.identifiable_f&mm.identifiable_h]
print("hiO2 vs full: r ratio median",np.nanmedian(mm.r_h/mm.r_f).round(3),"IQR",np.nanpercentile(mm.r_h/mm.r_f,[25,75]).round(3),"| K ratio",np.nanmedian(mm.K_h/mm.K_f).round(3),np.nanpercentile(mm.K_h/mm.K_f,[25,75]).round(3))
print("measured O2 at window start: ",full.O2_start_meas.describe().round(2).to_dict()); print(full.groupby("T").O2_start_meas.median().round(2).to_dict())
print("truncated duration (min):",hi.dur_min.describe().round(0).to_dict())
