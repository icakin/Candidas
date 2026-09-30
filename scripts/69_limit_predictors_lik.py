"""Idea 2: which fitted quantity predicts the isolate's growth limit best?
Candidates (per-isolate posterior medians): the 1.5x/2x/3x cost crossings, the growth
optimum T_opt, the carbon-economy optimum T_min(CUE), the deactivation midpoint Th, and the
temperature at which growth falls to half its peak (T50). Each is used alone as the
predictor in the interval-censored normal model of the limit (loss in (T_last, T_loss],
censored >44). Compared by maximised log-likelihood (same n, same df) and by
leave-one-isolate-out predictive log density."""
import os, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv")
cs=pd.read_csv(f"{C}/results/tables/isolate_crossings_summary.csv").rename(columns={"Isolate":"isolate"})
dr=pd.read_csv(f"{C}/results/tables/isolate_crossings_draws.csv").rename(columns={"Isolate":"isolate"})
pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv").rename(columns={"Isolate":"isolate"})
pp=pp[~pp.isolate.str.startswith("GROUP_")]
def lnG(TC,E,Eh,Th):
    TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK); return E*CINV*(1/TREF-1/TK)-np.log1p(np.exp(np.minimum(u,700)))
def topt(Es,Eh,Th):
    s=Es/Eh; s=np.where((s<=0)|(s>=1),np.nan,s); return 1/(1/Th-np.log(s/(1-s))/(Eh*CINV))-273.15
# T50 per draw: growth falls to half its peak on the warm limb (bisection on a grid)
rows=[]
for iso,g in pp.groupby("isolate"):
    To=topt(g.E.values,g.Eh.values,g.Th.values); T50=np.full(len(g),np.nan)
    grid=np.linspace(30,70,801)
    for i in range(len(g)):
        if not np.isfinite(To[i]): continue
        gg=lnG(grid,g.E.values[i],g.Eh.values[i],g.Th.values[i]); g0=lnG(To[i],g.E.values[i],g.Eh.values[i],g.Th.values[i])
        ok=(grid>To[i])&(gg<g0+np.log(0.5)); T50[i]=grid[ok][0] if ok.any() else np.nan
    rows.append(dict(isolate=iso,Th_C=np.median(g.Th.values)-273.15,T50=np.nanmedian(T50),Topt_draws=np.nanmedian(To)))
X=pd.DataFrame(rows).merge(cs[["isolate","T_opt_med","T_1.5x_med","T_2x_med","T_3x_med"]],on="isolate")
tm=dr.groupby("isolate").T_min_cost.median().rename("Tmin_CUE").reset_index(); X=X.merge(tm,on="isolate")
d=t.merge(X,on="isolate"); lo=d.T_last_growth.values.astype(float); hi=np.where(d.growth_censored,np.inf,d.T_growthloss.values.astype(float))
PRED={"3x crossing":"T_3x_med","2x crossing":"T_2x_med","1.5x crossing":"T_1.5x_med","growth optimum T_opt":"T_opt_med",
      "CUE optimum T_min":"Tmin_CUE","deactivation midpoint Th":"Th_C","growth half-peak T50":"T50"}
def nll(p,x,lo,hi):
    a,b,ls=p; m=a+b*x; sd=np.exp(ls); F=lambda v: np.where(np.isinf(v),1.0,norm.cdf((v-m)/sd))
    return -np.sum(np.log(np.clip(F(hi)-F(lo),1e-300,None)))
def fit(x,lo,hi):
    best=None
    for a0 in (0,20,40):
        for b0 in (0.5,1.0,1.5):
            r=minimize(nll,[a0,b0,0.3],args=(x,lo,hi),method="Nelder-Mead",options=dict(xatol=1e-7,fatol=1e-9,maxiter=20000))
            if best is None or r.fun<best.fun: best=r
    return best
null=fit(np.zeros(len(d)),lo,hi)
res=[]
for lab,col in PRED.items():
    x=d[col].values.astype(float); b=fit(x,lo,hi)
    # LOO predictive log density
    lpd=[]
    for i in range(len(d)):
        m=np.ones(len(d),bool); m[i]=False; bi=fit(x[m],lo[m],hi[m])
        lpd.append(-nll(bi.x,x[[i]],lo[[i]],hi[[i]]))
    res.append(dict(predictor=lab,logLik=-b.fun,dLL_vs_null=null.fun-b.fun,slope=b.x[1],intercept=b.x[0],sigma=np.exp(b.x[2]),LOO_lpd=np.sum(lpd),
                    spearman=pd.Series(x).corr(pd.Series(np.where(d.growth_censored,46.0,d.T_growthloss)),method="spearman")))
R=pd.DataFrame(res).sort_values("logLik",ascending=False); pd.set_option("display.width",200)
print(R.round(3).to_string(index=False)); print(f"\nnull (no predictor) logLik = {-null.fun:.2f}")
R.to_csv(f"{C}/results/tables/limit_predictors_lik.csv",index=False); d.to_csv(f"{C}/results/tables/limit_predictors_isolates.csv",index=False)
print("\nper-isolate predictors:"); print(d[["isolate","T_last_growth","T_growthloss","T_3x_med","T_opt_med","Tmin_CUE","Th_C","T50"]].round(1).to_string(index=False))
