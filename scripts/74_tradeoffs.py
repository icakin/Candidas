"""Idea 5: trade-offs across the 20 isolates. Per isolate and posterior draw: growth limit
(likelihood-model limit; censored isolates imputed uniformly on (44,48]), peak growth rate
r_max (h^-1), growth optimum T_opt, CUE optimum T_min, maximum CUE proxy = min(K/r)
(scale-free, lower = cheaper), and the 37->40 C fever cost. Draw-wise Spearman correlations
(posterior median and 95% CrI) between the limit and each trait, over all 20 isolates and
within C. auris (12). Phylogenetic caveat applies: four clades of one species."""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15; rng=np.random.default_rng(7)
pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv"); pp=pp[~pp.Isolate.str.startswith("GROUP_")]
q=pd.read_csv(f"{C}/results/tables/group_quota.csv").set_index("Group").quota_fgC
t=pd.read_csv(f"{C}/results/tables/lik_transition.csv").set_index("isolate")
def lnG(TC,E,Eh,Th):
    TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK); return E*CINV*(1/TREF-1/TK)-np.log1p(np.exp(np.minimum(u,700)))
def ln_cost(TC,E,ER,Eh,Th):
    TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK); return (ER-E)*CINV*(1/TREF-1/TK)+np.where(u>30,u,np.log1p(np.exp(np.minimum(u,30))))
def topt(Es,Eh,Th):
    s=Es/Eh
    with np.errstate(invalid="ignore",divide="ignore"): return np.where((s>0)&(s<1),1/(1/Th-np.log(s/(1-s))/(Eh*CINV))-273.15,np.nan)
isos=sorted(pp.Isolate.unique()); ND=pp.groupby("Isolate").size().min()
TR={k:np.full((len(isos),ND),np.nan) for k in ["limit","r_max","T_opt","T_min","cost37_40"]}
for i,iso in enumerate(isos):
    g=pp[pp.Isolate==iso].iloc[:ND]; E,Eh,Th,ER,lnB0=[g[c].values for c in ("E","Eh","Th","ER","lnB0")]
    To=topt(E,Eh,Th); Tm=topt(E-ER,Eh,Th)
    TR["T_opt"][i]=To; TR["T_min"][i]=Tm
    TR["r_max"][i]=np.exp(lnB0+lnG(To,E,Eh,Th))/q[g.Group.iloc[0]]
    TR["cost37_40"][i]=np.exp(ln_cost(40,E,ER,Eh,Th)-ln_cost(37,E,ER,Eh,Th))
    r=t.loc[iso]; TR["limit"][i]=np.where(r.growth_censored,rng.uniform(44,48,ND),r.T_growthloss)
auris=np.array([iso.startswith("Clade") for iso in isos])
def corr(a,b,mask):
    out=[]
    for j in range(ND):
        x,y=a[mask,j],b[mask,j]
        if np.all(np.isfinite(x)) and np.all(np.isfinite(y)): out.append(spearmanr(x,y).correlation)
    out=np.array(out); return np.median(out),np.percentile(out,2.5),np.percentile(out,97.5),np.mean(out>0)
rows=[]
for trait in ["r_max","T_opt","T_min","cost37_40"]:
    for lab,mask in (("all 20",np.ones(len(isos),bool)),("C. auris 12",auris),("relatives 8",~auris)):
        m,lo,hi,pp_=corr(TR["limit"],TR[trait],mask); rows.append(dict(trait=trait,set=lab,rho_med=m,lo=lo,hi=hi,P_pos=pp_))
R=pd.DataFrame(rows); pd.set_option("display.width",200); print("Spearman rho between the growth limit and each trait (posterior median, 95% CrI, P(rho>0)):"); print(R.round(2).to_string(index=False))
med=pd.DataFrame({k:np.nanmedian(v,axis=1) for k,v in TR.items()},index=isos); med["auris"]=auris
print("\nper-isolate medians:"); print(med.round(2).to_string())
R.to_csv(f"{C}/results/tables/tradeoffs.csv",index=False); med.to_csv(f"{C}/results/tables/tradeoffs_isolates.csv")
