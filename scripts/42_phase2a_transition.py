"""PHASE 2A - thermal transition, isolate-level. Frozen states, never recomputed.
Ambiguity handled ONLY as two bounds:
  BOUND P ("potentially productive") : ambiguous counted as confirmed growth
  BOUND R ("potentially resp-only")  : ambiguous counted as respiration without growth
T_growthloss = lowest temperature at which P(confirmed growth) first falls below 0.5
               AND stays below at every higher measured temperature (no recovery).
               Right-censored (>44) if it never falls below 0.5.
T_resploss   = lowest temperature at which P(detectable respiration) falls below 0.5,
               same rule; censored at >44 otherwise.
Window = T_resploss - T_growthloss : the respiration-without-growth interval.
"""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d=pd.read_csv(f"{C}/reprocess/unblinded_results.csv"); ISO=set(pd.read_csv(f"{C}/results/tables/otu_names.csv").OTU); d=d[d.OTU.isin(ISO)].copy()  # study isolates only
NM={1:"growth",2:"resp_only",3:"ambig",4:"no_resp"}; d["S"]=d.state.map(NM)
Ts=sorted(d["T"].unique())
def first_below(series_by_T, thresh=0.5):
    """lowest T where value<thresh and it stays <thresh at all higher measured T"""
    ts=sorted(series_by_T.index)
    for i,T in enumerate(ts):
        if series_by_T[T] < thresh and all(series_by_T[u] < thresh for u in ts[i:]):
            return T
    return np.nan            # censored: never persistently below
rows=[]
for (g,iso),s in d.groupby(["group","otu_name"]):
    for bound,growset in (("P",{"growth","ambig"}),("R",{"growth"})):
        pg=s.groupby("T").apply(lambda x: x.S.isin(growset).mean(), include_groups=False)
        # detectable respiration = anything that is not "no_resp";
        # under bound R an ambiguous well still shows respiration (it passed the R gate)
        pr=s.groupby("T").apply(lambda x: (x.S!="no_resp").mean(), include_groups=False)
        tg=first_below(pg); tr_=first_below(pr)
        rows.append(dict(group=g,isolate=iso,bound=bound,
            T_growthloss=tg, growth_censored=np.isnan(tg),
            T_resploss=tr_, resp_censored=np.isnan(tr_),
            window=(tr_-tg) if (np.isfinite(tg) and np.isfinite(tr_)) else np.nan,
            p_growth_38=pg.get(38,np.nan), p_growth_40=pg.get(40,np.nan),
            p_growth_42=pg.get(42,np.nan), p_growth_44=pg.get(44,np.nan),
            p_responly_40=(s[s["T"]==40].S.eq("resp_only").mean() if 40 in s["T"].values else np.nan)))
t=pd.DataFrame(rows)
t["clade"]=np.where(t.group.str.startswith("Clade"),"C. auris","relatives")
t.to_csv(f"{C}/results/tables/phase2a_transition.csv",index=False)
pd.set_option("display.width",220)
for b,lab in (("R","BOUND R - ambiguous = respiration without confirmed growth"),
              ("P","BOUND P - ambiguous = potentially productive (counted as growth)")):
    x=t[t.bound==b]
    print(f"\n=== {lab} ===")
    print(x[["clade","group","isolate","T_growthloss","growth_censored","T_resploss",
             "resp_censored","window","p_growth_40"]].to_string(index=False))
    print("\n  isolate-level summary (censored at >44 shown separately):")
    for cl,y in x.groupby("clade"):
        obs=y[~y.growth_censored]; cen=y[y.growth_censored]
        med=obs.T_growthloss.median() if len(obs) else np.nan
        print(f"    {cl:10s} n={len(y)}  growth-loss observed in {len(obs)} "
              f"(median {med if np.isfinite(med) else float('nan'):.0f} °C), censored >44 in {len(cen)}")
print("\n=== respiration-without-growth interval (BOUND R, observed only) ===")
w=t[(t.bound=="R")&t.window.notna()]
print(w[["clade","isolate","T_growthloss","T_resploss","window"]].to_string(index=False))
print(f"\n  isolates whose respiration outlasts growth by >=2 °C: {(w.window>=2).sum()} of {len(w)}")
