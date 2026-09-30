"""PHASE 3 - predict the interval-censored upper thermal limit from 22-38 °C predictors.
LOIO, complete isolate withheld. Censored isolates (>44) enter training as Tobit-style
lower bounds and are scored separately as 'correct if predicted >= 44'."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D=pd.read_csv(f"{C}/results/tables/phase3_design.csv")
D["lim"]=np.where(D.cen_R==1,44.0,D.limit_R); D["cen"]=(D.cen_R==1).astype(int)
SETS={"M1_growth":["r38","topt"],"M2_metab":["q38","qslope"],"M3_combined":None}
def zs(a,b):
    mu=a.mean(0); sd=a.std(0); sd=np.where(sd<1e-9,1,sd); return (a-mu)/sd,(b-mu)/sd
def ols(X,y):
    Xd=np.c_[np.ones(len(X)),X]
    return np.linalg.lstsq(Xd,y,rcond=None)[0]
rows=[]
for name,cols in SETS.items():
    preds=[]
    for i in range(len(D)):
        trn=np.ones(len(D),bool); trn[i]=False
        if cols is None:
            bg=max(["r38","topt"],key=lambda c: abs(np.corrcoef(D[c].values[trn],D.lim.values[trn])[0,1]))
            bm=max(["q38","qslope"],key=lambda c: abs(np.corrcoef(D[c].values[trn],D.lim.values[trn])[0,1]))
            cc=[bg,bm]
        else: cc=cols
        Xtr,Xte=zs(D.loc[trn,cc].values,D.loc[[i],cc].values)
        b=ols(Xtr,D.lim.values[trn]); preds.append(float(np.r_[1,Xte[0]]@b))
    D[f"pred_{name}"]=preds
    obs=D[D.cen==0]
    err=np.abs(obs[f"pred_{name}"]-obs.lim)
    cen=D[D.cen==1]
    ok=(cen[f"pred_{name}"]>=44).sum() if len(cen) else np.nan
    rows.append(dict(model=name,MAE_uncensored=err.mean(),median_AE=err.median(),
                     n_uncensored=len(obs),censored_correct=f"{ok}/{len(cen)}"))
r=pd.DataFrame(rows); r.to_csv(f"{C}/results/tables/phase3_limit_results.csv",index=False)
print("=== upper thermal limit, leave-one-isolate-out (bound R) ===")
print(r.round(3).to_string(index=False))
print(f"\n  M2-M1 MAE difference: {r.set_index('model').MAE_uncensored['M2_metab']-r.set_index('model').MAE_uncensored['M1_growth']:+.3f} °C "
      f"(positive = metabolism WORSE)")
print(f"  M3-M1 MAE difference: {r.set_index('model').MAE_uncensored['M3_combined']-r.set_index('model').MAE_uncensored['M1_growth']:+.3f} °C")
print("\n=== every isolate: observed vs predicted limit ===")
D["observed"]=np.where(D.cen==1,">44",D.lim.round(0).astype(str))
print(D[["isolate","observed","pred_M1_growth","pred_M2_metab","pred_M3_combined",
         "r38","topt","q38","qslope"]].round(2).to_string(index=False))
D.to_csv(f"{C}/results/tables/phase3_limit_per_isolate.csv",index=False)
print("\n=== the three isolates the plan names specifically ===")
for iso in ("Hae_1724","para_2052","Clade2_2073"):
    s=D[D.isolate==iso].iloc[0]
    print(f"  {iso:12s} observed {s.observed:>4s}  growth-model {s.pred_M1_growth:.1f}  "
          f"metabolism-model {s.pred_M2_metab:.1f}  combined {s.pred_M3_combined:.1f}")
