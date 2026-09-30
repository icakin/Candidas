"""Shared maths and data for the Python versions of Figs 4-5 and Supp Figs 1-4.
Same functions as fig_common.R (lnG, lnR, topt, tax, gpct), same constants."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
CINV=11604.51812; TREF=293.15; T_BODY,T_FEVER,T_MIN,T_MAX=37.0,40.0,22.0,44.0
TG=np.arange(T_MIN,T_MAX+1e-9,0.25)
def lnG(TC,E,Eh,Th):
    TC=np.asarray(TC,float); TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK)
    return E*CINV*(1/TREF-1/TK)-np.log1p(np.exp(np.minimum(u,700)))
def lnR(TC,ER): TK=np.asarray(TC,float)+273.15; return ER*CINV*(1/TREF-1/TK)
def ln_cost(TC,E,ER,Eh,Th):
    TK=np.asarray(TC,float)+273.15; u=Eh*CINV*(1/Th-1/TK)
    return (ER-E)*CINV*(1/TREF-1/TK)+np.where(u>30,u,np.log1p(np.exp(np.minimum(u,30))))
def topt(Es,Eh,Th):
    s=Es/Eh
    with np.errstate(invalid="ignore",divide="ignore"):
        return np.where((s>0)&(s<1),1/(1/Th-np.log(s/(1-s))/(Eh*CINV))-273.15,np.nan)
def tax(TC,E,ER,Eh,Th):
    t0=topt(E-ER,Eh,Th); t0=np.where(np.isfinite(t0),t0,T_MIN); t0=np.clip(t0,T_MIN,T_MAX)
    return np.exp(ln_cost(TC,E,ER,Eh,Th)-ln_cost(t0,E,ER,Eh,Th))
def gpct(TC,E,Eh,Th):
    tp=topt(E,Eh,Th); out=np.full(len(E),np.nan); k=np.isfinite(tp)
    out[k]=100*np.exp(lnG(TC,E[k],Eh[k],Th[k])-lnG(tp[k],E[k],Eh[k],Th[k])); return out
def load():
    pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv")
    G={g:pp[pp.Isolate==f"GROUP_{g}"].reset_index(drop=True) for g in GROUPS}
    I={iso:pp[pp.Isolate==iso].reset_index(drop=True) for iso in pp.Isolate.unique() if not iso.startswith("GROUP_")}
    quota=pd.read_csv(f"{C}/results/tables/group_quota.csv").set_index("Group").quota_fgC
    vals=pd.read_csv(f"{C}/results/tables/fig_values.csv").set_index("Group")
    return G,I,quota,vals
def curve(draws,fn,grid=TG):
    """apply fn(T, draws) over grid -> (median, lo, hi) arrays"""
    M=np.stack([fn(t,draws) for t in grid],1)
    return np.nanmedian(M,0),np.nanpercentile(M,2.5,0),np.nanpercentile(M,97.5,0)
def points():
    """per-well growth r (h^-1) and per-cell respiration, the wells that feed the fits"""
    d=pd.read_csv(f"{C}/results/tables/derived_N0_R_results_with_carbon.csv")
    d=d.merge(pd.read_csv(f"{C}/results/tables/otu_names.csv")[["otu_name","group"]],on="otu_name")
    d=d[d.group.isin(GROUPS)&(d.keep==True)&(d.fit_valid==True)&(d.has_curvature==True)].copy()
    d["r_h"]=d.r*60.0; return d
