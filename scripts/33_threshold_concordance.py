"""Concordance (NOT independent validation) between the cost-multiplier thresholds of Fig. 3
and the temperature at which growth is actually lost.

r enters both quantities: the cost ratio is K/r, and the growth-state call depends on the
pipeline accepting a growth fit, which depends on r and trace curvature. So agreement gives
the threshold biological grounding; it does not prove it independently.

Observed growth loss per taxon: hierarchical-free logistic of P(growing) on temperature,
fitted per taxon with isolate-clustered bootstrap, T_loss = temperature where P = 0.5.
Right-censoring: taxa still growing at 44 C have T_loss above the assay ceiling; these are
reported as censored and excluded from the concordance error, not extrapolated.
"""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); rng=np.random.default_rng(3); NB=2000
w=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv")
w["grow"]=(w.state=="growing").astype(int)
dr=pd.read_csv(f"{C}/results/tables/thermal_headroom_draws.csv")
MULT={"1.5x":"T_1.5xcost","2x":"T_2xcost","3x":"T_3xcost"}
def t50(df):
    """logistic P(grow) ~ T ; return T at P=0.5, or nan if never crosses in range"""
    T=df["T"].values.astype(float); y=df.grow.values
    if y.mean()>0.99 or y.mean()<0.01: return np.nan
    def res(p): 
        z=np.clip(p[0]*(T-p[1]),-30,30); return 1/(1+np.exp(-z))-y
    try:
        r=least_squares(res,[-1.0,np.median(T)],bounds=([-20,15],[-0.01,80]))
        return r.x[1]
    except Exception: return np.nan
rows=[]
for g,s in w.groupby("group"):
    obs=t50(s)
    # isolate-clustered bootstrap
    isos=s.OTU.unique(); bs=[]
    for _ in range(NB):
        pick=rng.choice(isos,len(isos),replace=True)
        bs.append(t50(pd.concat([s[s.OTU==i] for i in pick])))
    bs=np.array(bs); bs=bs[np.isfinite(bs)]
    cens = (s[s["T"]==s["T"].max()].grow.mean()>0.5)   # still growing at ceiling
    row=dict(taxon=g,T_growth_loss=obs,lo=np.percentile(bs,2.5) if len(bs)>50 else np.nan,
             hi=np.percentile(bs,97.5) if len(bs)>50 else np.nan,censored=cens,
             frac_grow_44=s[s["T"]==44].grow.mean())
    for k,col in MULT.items(): row[k]=dr[dr.Group==g][col].median()
    rows.append(row)
d=pd.DataFrame(rows).sort_values("T_growth_loss")
pd.set_option("display.width",200)
print(d.round(2).to_string(index=False)); print()
un=d[~d.censored]
print(f"Uncensored taxa (growth loss observed within 22-44 C): {list(un.taxon)}\n")
print("Concordance with observed growth loss, uncensored taxa only:")
best=None
for k in MULT:
    e=(un[k]-un.T_growth_loss); mae=e.abs().mean(); bias=e.mean()
    print(f"  {k:>4} : mean signed error {bias:+.2f} C   mean |error| {mae:.2f} C")
    if best is None or mae<best[1]: best=(k,mae)
print(f"\nClosest multiplier: {best[0]} (mean |error| {best[1]:.2f} C)")
d.to_csv(f"{C}/results/tables/threshold_concordance.csv",index=False)
