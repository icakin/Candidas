"""Supplementary Fig. 7 redrawn in the manuscript style from
results/tables/temperature_equilibration_sensitivity.csv (written by 08_...R)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
d=pd.read_csv(f"{C}/results/tables/temperature_equilibration_sensitivity.csv")
fig,axes=plt.subplots(1,3,figsize=(7.2,2.4),sharey=True,gridspec_kw=dict(wspace=0.18))
panels=[("growth per cell",None,None,None,MUTED),("respiration per cell","resp_swing_med","resp_swing_lo","resp_swing_p75",RED),
        ("apparent CUE","cue_swing_med","cue_swing_p25","cue_swing_hi",GROW)]
for ax,(lab,m,lo,hi,c),L in zip(axes,panels,"abc"):
    ax.axhline(0,color=AMB,lw=0.8,zorder=1)
    ax.axvspan(36.8,37.2,color=INK,alpha=0.08,lw=0); ax.axvspan(39.8,40.2,color=INK,alpha=0.05,lw=0)
    if m is None: ax.plot(d["T"],d.growth_swing,color=c,lw=1.6)
    else:
        ax.fill_between(d["T"],d[lo],d[hi],color=c,alpha=0.18,lw=0); ax.plot(d["T"],d[m],color=c,lw=1.6); ax.scatter(d["T"],d[m],s=10,color=c,zorder=3)
    ax.set_title(lab,fontsize=7.5,loc="left",pad=10); ax.set_xlabel("setpoint (°C)"); ax.set_xticks([22,26,30,34,38,42])
    ax.set_ylim(-35,30); panel_label(ax,L,x=-0.12 if L=="a" else -0.06,y=1.14)
axes[0].set_ylabel("% shift from the pipeline value")
save(fig,"FIG_SUPP_equilibration")
