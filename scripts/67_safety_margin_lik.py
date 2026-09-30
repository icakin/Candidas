"""Does the carbon-economy collapse predict the growth limit, isolate by isolate?
Per isolate: the temperature at which R/G reaches M x its own minimum (posterior from
65_isolate_crossings.R) against the observed loss of growth under the likelihood model
(lik_transition.csv: loss in (T_last, T_loss], or >44 if censored).
  - distance of the crossing (posterior median) to the observed interval (0 if inside)
  - posterior probability that the crossing lies inside the interval (from the draws)
  - interval-censored regression of T_loss on the crossing (slope 1, intercept 0 = prediction)
  - Spearman rank correlation (censored limits ranked as ties above 44)"""
import os, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm, chi2, spearmanr
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
s=pd.read_csv(f"{C}/results/tables/isolate_crossings_summary.csv").rename(columns={"Isolate":"isolate"})
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv")
d=t.merge(s,on="isolate")
d["lo"]=d.T_last_growth; d["hi"]=np.where(d.growth_censored,np.inf,d.T_growthloss)
dr=pd.read_csv(f"{C}/results/tables/isolate_crossings_draws.csv").rename(columns={"Isolate":"isolate"})
rows=[]
for M,col in (("1.5x","T_1.5x_med"),("2x","T_2x_med"),("3x","T_3x_med")):
    x=d[col].values; lo=d.lo.values; hi=d.hi.values
    dist=np.where(x<lo,x-lo,np.where(x>hi,x-hi,0.0))          # signed distance to interval
    unc=~d.growth_censored.values
    # posterior P(crossing in interval)
    dcol={"1.5x":"T_1.5xcost","2x":"T_2xcost","3x":"T_3xcost"}[M]
    pin=[]
    for iso,g in dr.groupby("isolate"):
        v=g[dcol].dropna().values; r=d[d.isolate==iso].iloc[0]
        pin.append(np.mean((v>r.lo)&(v<=r.hi)))
    pin=np.array(pin)
    # interval-censored regression T_loss ~ N(a + b x, s)
    def nll(p):
        a,b,ls=p; m=a+b*x; sd=np.exp(ls)
        F=lambda v: np.where(np.isinf(v),1.0,norm.cdf((v-m)/sd))
        return -np.sum(np.log(np.clip(F(hi)-F(lo),1e-300,None)))
    best=min((minimize(nll,[a0,b0,0.5],method="Nelder-Mead",options=dict(xatol=1e-7,fatol=1e-9,maxiter=20000)) for a0 in (0,20,40) for b0 in (0.5,1.0)),key=lambda r:r.fun)
    null=minimize(lambda p: nll([p[0],0.0,p[1]]),[41,0.5],method="Nelder-Mead").fun
    a,b,sd=best.x[0],best.x[1],np.exp(best.x[2]); LR=2*(null-best.fun)
    # fixed slope 1, intercept 0: is the crossing an unbiased predictor?
    nll10=nll([0.0,1.0,np.log(sd)]); 
    rho=spearmanr(x,np.where(d.growth_censored,46.0,d.T_growthloss)).correlation
    rows.append(dict(M=M,MAE_uncensored=np.abs(dist[unc]).mean(),bias_uncensored=dist[unc].mean(),
                     n_inside=int((dist==0).sum()),n_below=int((dist<0).sum()),n_above=int((dist>0).sum()),
                     mean_P_inside=pin.mean(),censored_consistent=int(((x>=44)&~unc).sum()),n_censored=int((~unc).sum()),
                     slope=b,intercept=a,sigma=sd,LR_vs_null=LR,p=chi2.sf(LR,1),spearman=rho))
    d[f"dist_{M}"]=dist; d[f"Pin_{M}"]=pin
pd.set_option("display.width",220)
print(pd.DataFrame(rows).round(3).to_string(index=False)); print()
print(d[["isolate","T_last_growth","T_growthloss","growth_censored","T_1.5x_med","T_2x_med","T_3x_med","dist_3x","Pin_3x"]].round(2).to_string(index=False))
d.to_csv(f"{C}/results/tables/safety_margin_lik.csv",index=False)

# ---- identity test at 3x and a supplementary figure ----------------------------
import sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
x=d.T_3x_med.values; lo=d.lo.values; hi=d.hi.values
def nll_id(ls):
    sd=np.exp(ls); F=lambda v: np.where(np.isinf(v),1.0,norm.cdf((v-x)/sd))
    return -np.sum(np.log(np.clip(F(hi)-F(lo),1e-300,None)))
def nll_free(p):
    a,b,ls=p; m=a+b*x; sd=np.exp(ls); F=lambda v: np.where(np.isinf(v),1.0,norm.cdf((v-m)/sd))
    return -np.sum(np.log(np.clip(F(hi)-F(lo),1e-300,None)))
fi=minimize(nll_id,[0.0],method="Nelder-Mead"); ff=min((minimize(nll_free,[a0,b0,0.3],method="Nelder-Mead",options=dict(maxiter=20000,xatol=1e-7,fatol=1e-9)) for a0 in (0,7) for b0 in (0.85,1.0)),key=lambda r:r.fun)
LRid=2*(fi.fun-ff.fun); print(f"\n3x crossing as identity predictor (T_loss = T_3x + e): sigma = {np.exp(fi.x[0]):.2f} C; LR vs free line = {LRid:.2f} on 2 df, p = {chi2.sf(LRid,2):.3f}")
fig,ax=plt.subplots(figsize=(3.6,3.4))
ax.plot([34,48],[34,48],color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1)
for _,r in d.iterrows():
    c=OI[r.group]
    ax.plot([r.T_3x_lo,r.T_3x_hi],[ (r.lo+r.hi)/2 if np.isfinite(r.hi) else 45.0]*2,color=c,lw=0.8,alpha=0.6,zorder=2)
    if np.isfinite(r.hi):
        ax.plot([r.T_3x_med]*2,[r.lo,r.hi],color=c,lw=2.4,alpha=0.35,solid_capstyle="butt",zorder=2)
        ax.scatter(r.T_3x_med,(r.lo+r.hi)/2,s=22,color=c,edgecolor="white",linewidth=0.5,zorder=3)
    else:
        ax.scatter(r.T_3x_med,45.0,s=26,marker="^",color=c,zorder=3)
ax.axhspan(44,46.3,color=AMB,alpha=0.25,lw=0); ax.text(34.6,45.1,"still growing at 44 °C (censored)",fontsize=5.8,color=MUTED,va="center")
ax.set_xlim(34,48); ax.set_ylim(34,46.3); ax.set_xticks([36,40,44,48]); ax.set_yticks([36,38,40,42,44])
ax.set_xlabel("temperature at which R/G reaches 3× its minimum (°C)\nposterior median, 95% CrI",fontsize=7); ax.set_ylabel("observed loss of growth (°C)\ngrid interval",fontsize=7)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.text(0.97,0.04,f"Spearman ρ = {rows[2]['spearman']:.2f} (n = 20)\nMAE {rows[2]['MAE_uncensored']:.2f} °C (12 uncensored)",transform=ax.transAxes,fontsize=6.2,va="bottom",ha="right",color=INK)
group_legend(ax,loc="upper left",bbox_to_anchor=(1.01,1.0),fontsize=5.8,handletextpad=0.3,labelspacing=0.25,frameon=False)
save(fig,"FIG_SUPP_safety_margin")
