"""Idea 3: kinetics of metabolic arrest in wells respiring without growth.
For each non-growing well (likelihood model, p >= 0.01, drawdown >= 1 mg/L), fit
O2(t) = O0 - (K/delta)(1 - exp(-delta t)) + plateau floor, i.e. a consumption rate that
decays exponentially at rate delta (h^-1) from K, on the same window as lik_model
(90 min to 15 h), against the constant-rate alternative (delta = 0). delta is the rate at
which respiration winds down; a half-life ln2/delta. Compared by AR(1) likelihood ratio as in
lik_model; delta is reported where the decaying model is preferred (LR > 10). Then: does
delta rise with temperature (Arrhenius), and does it differ between taxa?"""
import os, sys, numpy as np, pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"reprocess")); import lik_model as L
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
d=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv"); d=d[(d.p_boot>=0.01)&(d.drawdown>=1.0)].copy()
tr=pd.read_csv(f"{C}/reprocess/blind_traces.csv.gz"); tr=tr[tr.WELLID.isin(d.WELLID)]
def mean_decay(p,t):
    O0,K,dl,F=p; a=O0-(K*t if dl<1e-9 else K*(1-np.exp(-dl*t))/dl); return L._smoothmax(a,F)
def fit_decay(t,y):
    y0,ymin=float(y[0]),float(y.min()); K0=max((y0-ymin)/max(t[-1],1.0),1e-6); best=None
    for d0 in (0.0,0.002,0.01):
        try:
            r=least_squares(lambda p: mean_decay(p,t)-y,[y0,K0,d0,ymin],bounds=([-5,0,0,-1],[30,5,0.2,15]),max_nfev=3000)
            if best is None or r.cost<best.cost: best=r
        except Exception: pass
    return best
rows=[]
for w,g in tr.groupby("WELLID"):
    t,y=L._window(g.Time.values.astype(float),g.Oxygen.values.astype(float))
    if len(t)<60: continue
    fl=L._fit(t,y,False); fd=fit_decay(t,y)
    if fl is None or fd is None: continue
    ll_l,_,_=L._ar1_loglik(y-L._mean(fl.x,t,False)); ll_d,_,_=L._ar1_loglik(y-mean_decay(fd.x,t))
    rows.append(dict(WELLID=w,LR_decay=2*(ll_d-ll_l),K0_h=fd.x[1]*60,delta_h=fd.x[2]*60,halflife_h=np.log(2)/(fd.x[2]*60) if fd.x[2]>0 else np.inf))
R=pd.DataFrame(rows).merge(d[["WELLID","otu_name","group","T","K_hat","cons_v4" if "cons_v4" in d else "drawdown"]],on="WELLID")
R["decaying"]=R.LR_decay>10; R["auris"]=R.group.str.startswith("Clade")
R.to_csv(f"{C}/results/tables/arrest_kinetics.csv",index=False)
pd.set_option("display.width",200)
print(f"non-growing respiring wells: {len(R)}; decaying-rate model preferred (LR > 10) in {int(R.decaying.sum())}")
print("\nby group: n, n decaying, median delta (h^-1) and half-life (h) among decaying wells")
print(R.groupby("group").apply(lambda x: pd.Series(dict(n=len(x),n_decay=int(x.decaying.sum()),delta_med=x[x.decaying].delta_h.median(),halflife_med=x[x.decaying].halflife_h.median())),include_groups=False).round(3).to_string())
print("\nby temperature (all groups): n, n decaying, median delta, median half-life")
print(R.groupby("T").apply(lambda x: pd.Series(dict(n=len(x),n_decay=int(x.decaying.sum()),delta_med=x[x.decaying].delta_h.median(),halflife_med=x[x.decaying].halflife_h.median())),include_groups=False).round(3).to_string())
# Arrhenius of delta among decaying wells
z=R[R.decaying&(R.delta_h>0)].copy(); z["x"]=CINV*(1/TREF-1/(z["T"]+273.15)); z["ln_delta"]=np.log(z.delta_h)
import statsmodels.api as sm
X=pd.get_dummies(z.group,drop_first=True).astype(float); X["x"]=z.x; X=sm.add_constant(X)
m=sm.OLS(z.ln_delta.values,X.values).fit(cov_type="cluster",cov_kwds={"groups":pd.factorize(z.otu_name)[0]})
i=list(X.columns).index("x"); print(f"\nArrhenius slope of delta (with group intercepts): E_delta = {m.params[i]:.2f} eV (SE {m.bse[i]:.2f}), p = {m.pvalues[i]:.4f}, n = {len(z)}")
