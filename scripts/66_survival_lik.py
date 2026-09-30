"""Interval-censored model for the isolate thermal limit (likelihood-model limits).
Each isolate's loss temperature lies in (T_last_growth, T_growthloss], the 2 C grid step over
which growth was lost, or in (44, inf) if still growing at the assay ceiling. Two parametric
forms, both with a C. auris indicator:
  location-shift (normal):   T_loss ~ N(mu + beta * auris, sigma)  -> beta = shift in C
  Weibull proportional hazards on (T - 30): hazard ratio exp(-b) for C. auris vs relatives
Profile-likelihood 95% CIs and LR tests; a sensitivity fit without Hae_1724."""
import os, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm, chi2
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv")
def data(df):
    lo=df.T_last_growth.values.astype(float); hi=np.where(df.growth_censored,np.inf,df.T_growthloss.values.astype(float))
    x=(df.group.str.startswith("Clade")).astype(float).values; return lo,hi,x
def nll_norm(p,lo,hi,x):
    mu,b,ls=p; s=np.exp(ls); m=mu+b*x
    F=lambda v: np.where(np.isinf(v),1.0,norm.cdf((v-m)/s))
    return -np.sum(np.log(np.clip(F(hi)-F(lo),1e-300,None)))
def nll_weib(p,lo,hi,x,T0=30.0):
    la,b,lk=p; k=np.exp(lk); lam=np.exp(la+b*x)   # scale per group
    S=lambda v: np.where(np.isinf(v),0.0,np.exp(-((np.maximum(v-T0,1e-9))/lam)**k))
    return -np.sum(np.log(np.clip(S(lo)-S(hi),1e-300,None)))
def fit(nll,p0,args):
    best=None
    n=len(p0); starts=[np.zeros(n)]+[np.eye(n)[i]*d for i in range(n) for d in (1.5,-1.5)]
    for s in starts:
        r=minimize(nll,np.array(p0,float)+s,args=args,method="Nelder-Mead",options=dict(xatol=1e-7,fatol=1e-9,maxiter=20000))
        if best is None or r.fun<best.fun: best=r
    return best
def profile_ci(nll,best,idx,args,grid):
    out=[]
    for v in grid:
        def f(q):
            p=np.insert(q,idx,v); return nll(p,*args)
        q0=np.delete(best.x,idx); r=minimize(f,q0,method="Nelder-Mead",options=dict(xatol=1e-7,fatol=1e-9,maxiter=20000))
        out.append(r.fun)
    out=np.array(out); ok=grid[out<=best.fun+chi2.ppf(0.95,1)/2]
    return (ok.min(),ok.max()) if len(ok) else (np.nan,np.nan)
for lab,df in (("all 20 isolates",t),("excluding Hae_1724",t[t.isolate!="Hae_1724"])):
    lo,hi,x=data(df); print(f"=== {lab} ===")
    b=fit(nll_norm,[41,3,np.log(2)],(lo,hi,x)); b0=fit(lambda p,lo,hi,x: nll_norm([p[0],0.0,p[1]],lo,hi,x),[41,np.log(2)],(lo,hi,x))
    LR=2*(b0.fun-b.fun); ci=profile_ci(nll_norm,b,1,(lo,hi,x),np.linspace(-2,14,161))
    print(f"  location shift: relatives mu = {b.x[0]:.2f} C, C. auris shift = {b.x[1]:+.2f} C (95% profile CI {ci[0]:+.2f} to {ci[1]:+.2f}), sigma = {np.exp(b.x[2]):.2f}; LR = {LR:.2f}, p = {chi2.sf(LR,1):.4f}")
    w=fit(nll_weib,[np.log(10),0.3,np.log(6)],(lo,hi,x)); w0=fit(lambda p,lo,hi,x: nll_weib([p[0],0.0,p[1]],lo,hi,x),[np.log(10),np.log(6)],(lo,hi,x))
    LRw=2*(w0.fun-w.fun); ciw=profile_ci(nll_weib,w,1,(lo,hi,x),np.linspace(-0.5,2.0,126))
    k=np.exp(w.x[2]); hr=np.exp(-k*w.x[1]); hr_ci=(np.exp(-k*ciw[1]),np.exp(-k*ciw[0]))
    print(f"  Weibull PH: shape k = {k:.2f}; hazard of growth loss, C. auris vs relatives = {hr:.3f} (95% profile CI {hr_ci[0]:.3f} to {hr_ci[1]:.3f}); LR = {LRw:.2f}, p = {chi2.sf(LRw,1):.4f}")
    med_rel=w.x[0]; print(f"  Weibull median loss temperature: relatives {30+np.exp(w.x[0])*np.log(2)**(1/k):.1f} C, C. auris {30+np.exp(w.x[0]+w.x[1])*np.log(2)**(1/k):.1f} C")
