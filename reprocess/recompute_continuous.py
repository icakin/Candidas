"""Step 2: regenerate r, K and their uncertainties from the RAW traces using the
FROZEN automatic windows. Classifications are NOT recomputed or changed - they are
read from blind_states.csv. r_se is the bootstrap SD of r (AR(1)-matched residuals),
a downstream quantity only; it enters no classification decision."""
import pandas as pd, numpy as np, sys, hashlib; sys.path.insert(0,'.')
import classifier as C
st=pd.read_csv("unblinded_results.csv"); st=st[st.OTU.isin(set(pd.read_csv("../results/tables/otu_names.csv").OTU))].copy(); tr=pd.read_csv("blind_traces.csv.gz")
g=dict(list(tr.groupby("WELLID")))
rows=[]
for i,w in enumerate(st.WELLID):
    d=g[w]; t=d.Time.values; y=d.Oxygen.values
    m,reason=C.usable_interval(t,y)
    rec=dict(WELLID=w,win_start=np.nan,win_end=np.nan,dur_h=np.nan,n_pts=np.nan,
             r_auto=np.nan,K_auto=np.nan,r_se=np.nan,resid_sd=np.nan,O2_0_auto=np.nan,drawdown=np.nan)
    if m is not None:
        tt,yy=C.decimate(t[m],y[m]); t0=tt[0]; tr2=tt-t0
        rec.update(win_start=float(tt[0]),win_end=float(tt[-1]),dur_h=(tt[-1]-tt[0])/60.0,
                   n_pts=len(tt),drawdown=float(yy[0]-yy.min()))
        f2=C.fit_m2(tr2,yy)
        if f2 is not None:
            rho=C.lag1(f2.fun); sd=np.std(f2.fun)
            rec.update(O2_0_auto=float(f2.x[0]),K_auto=float(f2.x[1]),r_auto=float(f2.x[2]),resid_sd=float(sd))
            rr=np.random.default_rng(int(hashlib.md5(w.encode()).hexdigest()[:8],16))
            fit=C._m2(f2.x,tr2); boot=[]
            for _ in range(C.N_BOOT):
                e=np.empty(len(tr2)); e[0]=rr.normal(0,sd)
                nz=rr.normal(0,sd*np.sqrt(max(1-rho**2,1e-6)),len(tr2))
                for q in range(1,len(tr2)): e[q]=rho*e[q-1]+nz[q]
                b=C.fit_m2(tr2,fit+e,start=f2.x)
                if b is not None: boot.append(b.x[2])
            if len(boot)>=50: rec["r_se"]=float(np.std(boot))
    rows.append(rec)
    if (i+1)%300==0: print(f"  {i+1}",flush=True)
out=pd.DataFrame(rows).merge(st,on="WELLID")
out["r_auto_h"]=out.r_auto*60; out["K_auto_h"]=out.K_auto*60; out["r_se_h"]=out.r_se*60
out.to_csv("automatic_continuous.csv",index=False)
NM={1:"growth",2:"resp_only",3:"ambig",4:"no_resp"}
out["S"]=out.state.map(NM)
for s in ("growth","ambig","resp_only"):
    sub=out[(out.S==s)]
    sub.to_csv(f"dataset_{s}.csv",index=False)
    print(f"dataset_{s}.csv : {len(sub)} wells")
print(f"\nfull audit table automatic_continuous.csv : {len(out)} wells")
print(out.groupby("S")[["r_auto_h","K_auto_h","r_se_h","dur_h"]].median().round(4).to_string())
