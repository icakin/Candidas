"""PHASE 3 - can metabolism at 22-38 °C predict frozen 40-44 °C outcomes better than growth?
PREDICTORS: original HAND-TRIMMED fits, restricted to 22-38 °C. Taxon identity is never used.
OUTCOMES  : frozen v4 states at 40-44 °C, and the interval-censored upper thermal limit.
Ambiguous outcomes are handled as INTERVALS (both bounds), never reallocated.
Validation: leave-one-ISOLATE-out, the complete isolate withheld.
"""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
leg=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv")
on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
leg=leg.merge(on,on="OTU"); leg=leg[(leg.keep)&(leg.fit_valid)&(leg.r>0)&(leg.K>0)]
tr=leg[leg["T"]<=38].copy()                      # <-- ONLY 22-38 °C
tr["rh"]=tr.r*60; tr["q"]=np.log((tr.K*60)/tr.rh)
# ---- PREDICTORS, per isolate -------------------------------------------------
rows=[]
for iso,s in tr.groupby("otu_name"):
    r38=s[s["T"]==38].rh.median()
    m=s.groupby("T").rh.median()
    Tv=m.index.values.astype(float); lv=np.log(m.values)
    topt=np.nan
    if len(Tv)>=4:
        cf=np.polyfit(Tv,lv,2)
        if cf[0]<0: topt=-cf[1]/(2*cf[0])
    q38=s[s["T"]==38].q.median()
    w=s[s["T"].isin([34,36,38])].groupby("T").q.median()
    qslope=np.polyfit(w.index.values.astype(float),w.values,1)[0] if len(w)>=3 else np.nan
    rows.append(dict(isolate=iso,group=s.group.iloc[0],r38=r38,topt=topt,q38=q38,qslope=qslope))
X=pd.DataFrame(rows)
# ---- OUTCOMES from the frozen states -----------------------------------------
st=pd.read_csv(f"{C}/reprocess/unblinded_results.csv"); ISO=set(pd.read_csv(f"{C}/results/tables/otu_names.csv").OTU); st=st[st.OTU.isin(ISO)].copy()  # study isolates only
NM={1:"growth",2:"resp_only",3:"ambig",4:"no_resp"}; st["S"]=st.state.map(NM)
out=[]
for iso,s in st.groupby("otu_name"):
    d={"isolate":iso}
    for T in (40,42,44):
        w=s[s["T"]==T]
        d[f"g{T}_lo"]=int((w.S=="growth").sum()>=3)                       # ambiguous NOT growth
        d[f"g{T}_hi"]=int((w.S.isin(["growth","ambig"])).sum()>=3)        # ambiguous IS growth
        d[f"resp{T}"]=int((w.S=="resp_only").sum()>=3)
    out.append(d)
Y=pd.DataFrame(out)
t2=pd.read_csv(f"{C}/results/tables/phase2a_transition.csv")
lim=t2.pivot_table(index="isolate",columns="bound",values="T_growthloss",dropna=False)
cen=t2.pivot_table(index="isolate",columns="bound",values="growth_censored",dropna=False)
lim.columns=["limit_P","limit_R"]; cen.columns=["cen_P","cen_R"]
D=X.merge(Y,on="isolate").merge(lim.reset_index(),on="isolate").merge(cen.reset_index(),on="isolate")
D.to_csv(f"{C}/results/tables/phase3_design.csv",index=False)
pd.set_option("display.width",220)
print("=== predictors from 22-38 °C only, and frozen 40-44 °C outcomes ===")
print(D[["isolate","r38","topt","q38","qslope","g40_lo","g40_hi","g42_lo","g44_lo","limit_R","cen_R"]].round(3).to_string(index=False))
print(f"\n  predictor completeness: r38 {D.r38.notna().sum()}/20, topt {D.topt.notna().sum()}/20, "
      f"q38 {D.q38.notna().sum()}/20, qslope {D.qslope.notna().sum()}/20")
print("\n  outcome balance:")
for k in ["g40_lo","g40_hi","g42_lo","g42_hi","g44_lo","g44_hi"]:
    print(f"    {k}: {int(D[k].sum())} positive of {len(D)}")
