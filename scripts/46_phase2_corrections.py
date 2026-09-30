"""PHASE 2 CORRECTIONS 1-3. No state is changed; no threshold is touched."""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
t=pd.read_csv(f"{C}/results/tables/phase2a_transition.csv")
cons=pd.read_csv(f"{C}/results/tables/phase2b_consumption.csv")
rng=np.random.default_rng(2026)
# ============ CORRECTION 1: right-censoring-aware comparison =================
def logrank(time,event,grp):
    """standard log-rank statistic; event=1 observed, 0 right-censored"""
    ts=np.sort(np.unique(time[event==1])); O_E=0.0; V=0.0
    for tt in ts:
        at=time>=tt
        n=at.sum(); n1=(at&(grp==1)).sum()
        d=((time==tt)&(event==1)).sum(); d1=((time==tt)&(event==1)&(grp==1)).sum()
        if n<2: continue
        E1=d*n1/n; O_E+=d1-E1
        V+=d*(n1/n)*(1-n1/n)*((n-d)/(n-1))
    return (O_E**2/V) if V>0 else np.nan, O_E, V
print("=== CORRECTION 1: censoring-aware comparison of confirmed-growth loss ===")
print("   Right-censored isolates are those still growing at 44 °C: reported as >44,")
print("   never imputed as 46. Log-rank statistic with EXACT permutation (n=20).\n")
for b,lab in (("R","bound R (ambiguous = respiration-only)"),("P","bound P (ambiguous = productive)")):
    for excl,tag in ((None,"all isolates"),("Hae_1724","excluding Hae_1724")):
        x=t[t.bound==b].copy()
        if excl: x=x[x.isolate!=excl]
        time=np.where(x.growth_censored,44.0,x.T_growthloss.values)
        event=(~x.growth_censored).astype(int).values
        grp=(x.clade=="C. auris").astype(int).values
        obs,OE,V=logrank(time,event,grp)
        null=[]
        for _ in range(20000):
            g=rng.permutation(grp); s,_,_=logrank(time,event,g)
            if np.isfinite(s): null.append(s)
        p=(np.sum(np.array(null)>=obs)+1)/(len(null)+1)
        na=int(grp.sum()); nr=int((1-grp).sum())
        cen_a=int(x[x.clade=="C. auris"].growth_censored.sum()); cen_r=int(x[x.clade=="relatives"].growth_censored.sum())
        print(f"  {lab}, {tag}")
        print(f"    C. auris  n={na} ({cen_a} censored >44)   relatives n={nr} ({cen_r} censored >44)")
        print(f"    log-rank chi2 = {obs:.2f}   exact permutation p = {p:.4f}")
    print()
print("=== descriptive bounds, no imputation ===")
for b in ("R","P"):
    x=t[t.bound==b]
    for cl in ("C. auris","relatives"):
        y=x[x.clade==cl]; obs=y[~y.growth_censored]
        s=", ".join(sorted(f"{v:.0f}" for v in obs.T_growthloss)) or "none"
        print(f"  bound {b} {cl:10s}: observed losses at [{s}] °C; {int(y.growth_censored.sum())} isolate(s) >44")
print("\n=== respiration-loss temperature ===")
print(f"  detectable respiration never fell below 50% in any isolate at any temperature.")
print(f"  T_resploss > 44 °C for all {t.isolate.nunique()} isolates; the respiration-without-")
print( "  confirmed-growth interval is LOWER-BOUNDED (>= 44 - T_growthloss), never estimated.")
