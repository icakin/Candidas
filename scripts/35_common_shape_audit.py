"""Audit of 34_common_shape.py: residual structure, robustness, per-isolate deltas."""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy import stats
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv"); on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
d=cd.merge(on,on="OTU"); d=d[(d.keep)&(d.fit_valid)&(d.r>0)&(d.K>0)].copy()
d["y"]=np.log(d.K/d.r); d["TK"]=d["T"]+273.15
def build(sub):
    TAX=sorted(sub.group.unique()); ti={t:i for i,t in enumerate(TAX)}; n=len(TAX)
    def predict(p,df,m):
        k=df.group.map(ti).values; TK=df.TK.values; a=p[:n]
        if m==1:   dE=np.full(n,p[n]); Eh=np.full(n,p[n+1]); Th=p[n+2:2*n+2]
        elif m==2: dE=np.full(n,p[n]); Eh=p[n+1:2*n+1];      Th=p[2*n+1:3*n+1]
        else:      dE=p[n:2*n];        Eh=p[2*n:3*n];        Th=p[3*n:4*n]
        u=np.clip(Eh[k]*CINV*(1/Th[k]-1/TK),-50,50)
        return a[k]+dE[k]*CINV*(1/TREF-1/TK)+np.log1p(np.exp(u))
    a0=[np.log(sub[sub.group==t].K.median()/sub[sub.group==t].r.median()) for t in TAX]
    st={1:np.array(a0+[-0.4,3.0]+[309.]*n),2:np.array(a0+[-0.4]+[3.]*n+[309.]*n),3:np.array(a0+[-0.4]*n+[3.]*n+[309.]*n)}
    bd={1:([-20]*n+[-4,0.3]+[295]*n,[20]*n+[4,40]+[330]*n),
        2:([-20]*n+[-4]+[0.3]*n+[295]*n,[20]*n+[4]+[40]*n+[330]*n),
        3:([-20]*n+[-4]*n+[0.3]*n+[295]*n,[20]*n+[4]*n+[40]*n+[330]*n)}
    def fit(df,m): return least_squares(lambda p:predict(p,df,m)-df.y.values,st[m],bounds=bd[m],max_nfev=20000)
    return predict,fit
def loio(sub,label):
    predict,fit=build(sub); rows=[]
    for iso in sorted(sub.otu_name.unique()):
        tr=sub[sub.otu_name!=iso]; te=sub[sub.otu_name==iso]
        if tr[tr.group==te.group.iloc[0]].empty: continue
        r={}
        for m in (1,2,3):
            f=fit(tr,m); sd=np.sqrt(np.mean((predict(f.x,tr,m)-tr.y.values)**2))
            e=predict(f.x,te,m)-te.y.values
            r[f"lpd{m}"]=np.mean(-0.5*np.log(2*np.pi*sd**2)-e**2/(2*sd**2))
        rows.append(dict(isolate=iso,taxon=te.group.iloc[0],**r))
    cv=pd.DataFrame(rows); dl=cv.lpd3-cv.lpd1
    print(f"  {label:38s} n={len(cv):2d}  dLPD(M3-M1) {dl.mean():+.4f} +/- {dl.std(ddof=1)/np.sqrt(len(dl)):.4f}   favour M3: {(dl>0).sum()}/{len(dl)}")
    return cv
# ---- (5) per-isolate deltas, full data -------------------------------------
print("=== (5) per-isolate M3 - M1 predictive difference, full data ===")
cv=loio(d,"ALL (reference)")
cv["d31"]=cv.lpd3-cv.lpd1; cv["d21"]=cv.lpd2-cv.lpd1
print(cv[["isolate","taxon","d31","d21"]].round(4).to_string(index=False))
t=stats.ttest_1samp(cv.d31,0); w=stats.wilcoxon(cv.d31)
print(f"\n  paired t on dLPD(M3-M1): t={t.statistic:+.2f} p={t.pvalue:.3f} | Wilcoxon p={w.pvalue:.3f}")
# ---- (6) residual structure under M1 ---------------------------------------
print("\n=== (6) residual structure under M1 (shared shape) ===")
predict,fit=build(d); f1=fit(d,1); d["res1"]=d.y.values-predict(f1.x,d,1)
piv=d.groupby(["group","T"]).res1.mean().unstack()
print("  mean M1 residual by taxon x temperature:"); print(piv.round(2).to_string())
# does taxon x temperature interaction explain residual variance beyond noise?
import itertools
grand=d.res1.var(ddof=1)
cell=d.groupby(["group","T"]).res1.agg(["mean","count"])
ss_between=(cell["count"]*cell["mean"]**2).sum(); df_b=len(cell)-1
ss_within=sum(((g.res1-g.res1.mean())**2).sum() for _,g in d.groupby(["group","T"])); df_w=len(d)-len(cell)
F=(ss_between/df_b)/(ss_within/df_w); p=1-stats.f.cdf(F,df_b,df_w)
print(f"\n  cell-mean F-test on M1 residuals: F({df_b},{df_w})={F:.2f}  p={p:.3g}")
print(f"  -> systematic taxon x temperature structure REMAINS under the shared shape" if p<0.05 else "  -> no detectable structure")
f3=fit(d,3); rss1=np.sum((predict(f1.x,d,1)-d.y.values)**2); rss3=np.sum((predict(f3.x,d,3)-d.y.values)**2)
k1,k3=len(f1.x),len(f3.x); n=len(d); Fm=((rss1-rss3)/(k3-k1))/(rss3/(n-k3)); pm=1-stats.f.cdf(Fm,k3-k1,n-k3)
print(f"  in-sample M1 vs M3 nested F({k3-k1},{n-k3})={Fm:.2f}  p={pm:.3g}")
# ---- (7) robustness ---------------------------------------------------------
print("\n=== (7) robustness ===")
loio(d[d.otu_name!="Hae_1724"],"excluding Hae_1724")
for t_ in sorted(d.group.unique()): loio(d[d.group!=t_],f"excluding taxon {t_}")
