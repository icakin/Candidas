"""Supplementary Figs 8 and 9 redrawn in the manuscript style: per-cell respiration (Arrhenius)
and apparent CUE per group, credible band plus the temperature-equilibration envelope
(same quadrature widening as 11_bayesian_plots.R). Curves from the posterior draws dumped by
68_isolate_params_dump.R; CUE curves from results/tables/bayes_cue_by_clade.csv; points from
results/tables/derived_N0_R_results_with_carbon.csv."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
CINV=11604.51812; TREF=293.15
teq=pd.read_csv(f"{C}/results/tables/temperature_equilibration_sensitivity.csv")
pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv"); G=pp[pp.Isolate.str.startswith("GROUP_")].copy(); G["Group"]=G.Isolate.str.replace("GROUP_","")
pts=pd.read_csv(f"{C}/results/tables/derived_N0_R_results_with_carbon.csv").merge(pd.read_csv(f"{C}/results/tables/otu_names.csv")[["otu_name","group"]],on="otu_name")
pts=pts[pts.group.isin(GROUPS)&np.isfinite(pts.respiration_fgC_h)&(pts.respiration_fgC_h>0)]
TG=np.linspace(22,44,89)
_nw=pts.groupby(["group","T"]).size()
Tsup={g:max(t for (gg,t),k in _nw.items() if gg==g and k>=0.5*_nw[g].max()) for g in GROUPS}   # last T with >= half the group's wells growing
def rib(ax,T,lo,hi,lo_w,hi_w,med,c):
    ax.fill_between(T,lo_w,hi_w,color=c,alpha=0.08,lw=0); ax.plot(T,lo_w,color=c,lw=0.5,ls=(0,(2,2))); ax.plot(T,hi_w,color=c,lw=0.5,ls=(0,(2,2)))
    ax.fill_between(T,lo,hi,color=c,alpha=0.22,lw=0); ax.plot(T,med,color=c,lw=1.4)
# ---- Supp 8: respiration ----
fig,axes=plt.subplots(2,4,figsize=(7.2,4.0),sharex=True,sharey=True,gridspec_kw=dict(wspace=0.12,hspace=0.35)); axes=axes.ravel()
for ax,g in zip(axes,GROUPS):
    x=G[G.Group==g]; curves=np.exp(x.alpha.values[:,None]+x.ER.values[:,None]*CINV*(1/TREF-1/(TG[None,:]+273.15)))
    med=np.median(curves,0); lo=np.percentile(curves,2.5,0); hi=np.percentile(curves,97.5,0)
    rf=np.interp(TG,teq["T"],np.abs(teq.resp_swing_med)/100); hw=np.sqrt(((np.log(hi)-np.log(lo))/2)**2+np.log1p(rf)**2)
    p=pts[pts.group==g]; ax.scatter(p["T"]+np.random.RandomState(1).uniform(-.3,.3,len(p)),p.respiration_fgC_h,s=5,color=OI[g],alpha=0.35,lw=0,zorder=2)
    rib(ax,TG,lo,hi,np.exp(np.log(med)-hw),np.exp(np.log(med)+hw),med,OI[g])
    if Tsup[g]<44: ax.axvspan(Tsup[g],44,color="white",alpha=0.55,lw=0,zorder=3); ax.axvline(Tsup[g],color=MUTED,lw=0.6,ls=(0,(2,2)),zorder=4)
    ax.set_yscale("log"); ax.set_title(SHORT[g],fontsize=7.2,style="italic",loc="left",pad=3); ax.set_xticks([22,30,38,44])
    ax.axvspan(36.8,37.2,color=INK,alpha=0.08,lw=0); ax.axvspan(39.8,40.2,color=INK,alpha=0.05,lw=0)
axes[-1].axis("off")
axes[-1].text(0.02,0.85,"band: 95% credible interval\ndashed envelope: widened by the\ntemperature-equilibration component\npoints: vials\nfaded: beyond the last temperature at which\nat least half the group's vials grew",transform=axes[-1].transAxes,fontsize=6.4,va="top",color=MUTED,linespacing=1.4)
fig.supxlabel("temperature (°C)",fontsize=7.5,y=0.02); fig.supylabel("per-cell respiration (fg C cell⁻¹ h⁻¹)",fontsize=7.5,x=0.04)
save(fig,"FIG_SUPP_resp_envelope")
# ---- Supp 9: CUE ----
cue=pd.read_csv(f"{C}/results/tables/bayes_cue_by_clade.csv")
fig,axes=plt.subplots(2,4,figsize=(7.2,4.0),sharex=True,sharey=True,gridspec_kw=dict(wspace=0.12,hspace=0.35)); axes=axes.ravel()
for ax,g in zip(axes,GROUPS):
    x=cue[cue.Group==g].sort_values("T"); T=x["T"].values; med=x.CUE.values; lo=x.lo.values; hi=x.hi.values
    cf=np.interp(T,teq["T"],np.abs(teq.cue_swing_med)/100); hw=np.sqrt(((hi-lo)/2)**2+(cf*med)**2)
    rib(ax,T,lo,hi,np.clip(med-hw,0,1),np.clip(med+hw,0,1),med,OI[g])
    if Tsup[g]<44: ax.axvspan(Tsup[g],44,color="white",alpha=0.55,lw=0,zorder=3); ax.axvline(Tsup[g],color=MUTED,lw=0.6,ls=(0,(2,2)),zorder=4)
    ax.set_ylim(0,1); ax.set_title(SHORT[g],fontsize=7.2,style="italic",loc="left",pad=3); ax.set_xticks([22,30,38,44])
    ax.axvspan(36.8,37.2,color=INK,alpha=0.08,lw=0); ax.axvspan(39.8,40.2,color=INK,alpha=0.05,lw=0)
axes[-1].axis("off")
axes[-1].text(0.02,0.85,"band: 95% credible interval\ndashed envelope: widened by the\ntemperature-equilibration component\nfaded: beyond the last temperature at which\nat least half the group's vials grew",transform=axes[-1].transAxes,fontsize=6.4,va="top",color=MUTED,linespacing=1.4)
fig.supxlabel("temperature (°C)",fontsize=7.5,y=0.02); fig.supylabel("apparent CUE, G/(G+R)",fontsize=7.5,x=0.04)
save(fig,"FIG_SUPP_cue_envelope")
