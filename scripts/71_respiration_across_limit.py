"""Idea 6: does volumetric respiration keep rising through the growth limit?
Per well, the likelihood model gives K (initial O2 consumption rate, mg/L/min -> /h) for
growing and non-growing wells alike. Within each isolate, K is placed on an Arrhenius axis
and wells are labelled by their position relative to that isolate's growth limit (below,
at/above). Test: within isolate, does ln K above the limit continue the sub-limit Arrhenius
line (no break), fall below it (respiration collapses with growth), or exceed it?
Model: ln K = a_iso + E_K * (1/kT_ref - 1/kT) + delta * above + eps, fitted by OLS with
isolate fixed effects on wells with detectable respiration; delta is the log-offset of
above-limit wells from the isolate's own line; also a slope change (E_K x above)."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
d=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv"); d=d[d.drawdown>=1.0].copy()
d["grow"]=d.p_boot<0.01; d["K_h"]=d.K_hat  # K_hat already in mg/L/h (x60 in lik_model)
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv").set_index("isolate")
d["Tloss"]=d.otu_name.map(t.T_growthloss); d["cens"]=d.otu_name.map(t.growth_censored)
d["above"]=(~d.cens)&(d["T"]>=d.Tloss)
d["x"]=CINV*(1/TREF-1/(d["T"]+273.15)); d["lnK"]=np.log(d.K_h.clip(lower=1e-6))
d["auris"]=d.group.str.startswith("Clade")
import statsmodels.api as sm
def fit(df,slope_change=False):
    X=pd.get_dummies(df.otu_name,drop_first=False).astype(float); X["x"]=df.x; X["above"]=df.above.astype(float)
    if slope_change: X["x_above"]=df.x*df.above
    m=sm.OLS(df.lnK.values,X.values).fit(cov_type="cluster",cov_kwds={"groups":pd.factorize(df.otu_name+"_"+df["T"].astype(str))[0]})
    names=list(X.columns); return m,names
for lab,sub in (("all isolates",d),("relatives only",d[~d.auris]),("C. auris only",d[d.auris])):
    m,names=fit(sub); i=names.index("x"); j=names.index("above")
    print(f"{lab}: n wells {len(sub)}, above-limit wells {int(sub.above.sum())} | E_K = {m.params[i]:.2f} eV (SE {m.bse[i]:.2f}) | above-limit offset delta = {m.params[j]:+.2f} log units = x{np.exp(m.params[j]):.2f} (95% CI x{np.exp(m.params[j]-1.96*m.bse[j]):.2f} to x{np.exp(m.params[j]+1.96*m.bse[j]):.2f}), p = {m.pvalues[j]:.4f}")
    m2,n2=fit(sub,True); k=n2.index("x_above"); print(f"    slope change above the limit: {m2.params[k]:+.2f} eV (SE {m2.bse[k]:.2f}), p = {m2.pvalues[k]:.3f}")
# descriptive: median K by isolate at last-growth T vs first-loss T
rows=[]
for iso,g in d.groupby("otu_name"):
    if g.cens.iloc[0]: continue
    Tl=g.Tloss.iloc[0]; before=g[(g["T"]==Tl-2)]; at=g[g["T"]==Tl]; after=g[g["T"]==Tl+2]
    rows.append(dict(isolate=iso,T_loss=Tl,K_before=before.K_h.median(),K_at=at.K_h.median(),K_after=after.K_h.median() if len(after) else np.nan,
                     n_grow_at=int(at.grow.sum()),n_at=len(at)))
R=pd.DataFrame(rows); R["ratio_at_before"]=R.K_at/R.K_before; R["ratio_after_before"]=R.K_after/R.K_before
pd.set_option("display.width",200); print("\nK (mg O2 L-1 h-1) around the limit, per isolate (median of wells):"); print(R.round(2).to_string(index=False))
print(f"\nmedian K ratio at limit / 2 C below: {R.ratio_at_before.median():.2f}; 2 C above / 2 C below: {R.ratio_after_before.median():.2f}")
R.to_csv(f"{C}/results/tables/respiration_across_limit.csv",index=False)
