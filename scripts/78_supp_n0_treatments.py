"""Supplementary Fig. 10 redrawn in the manuscript style from
results/tables/bayes_resp_arr_E_three_treatments.csv (15_n0_treatment_panel.R; the five
groups of the original five-group refit)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
d=pd.read_csv(f"{C}/results/tables/bayes_resp_arr_E_three_treatments.csv")
TR=["With term (published)","Equilibration-corrected","Term-free"]; TL=["with the back-projection term","equilibration-corrected","term-free (N0 = N_inoc)"]
order=["Clade4","para","Clade3","Clade2","Clade1"]
fig,axes=plt.subplots(1,3,figsize=(7.2,2.3),sharey=True,sharex=True,gridspec_kw=dict(wspace=0.15))
for ax,t,tl,L in zip(axes,TR,TL,"abc"):
    x=d[d.treatment==t].set_index("clade")
    for i,g in enumerate(order):
        r=x.loc[g]; ax.plot([r.lwr,r.upr],[i,i],color=OI[g],lw=1.4,solid_capstyle="round"); ax.scatter(r.E,i,s=24,color=OI[g],edgecolor="white",linewidth=0.5,zorder=3)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([SHORT[g] for g in order],fontsize=6.6)
    ax.set_title(tl,fontsize=7.2,loc="left",pad=10); ax.set_xlabel("respiration activation energy E (eV)"); ax.set_xlim(0.15,0.75)
    ax.tick_params(axis="y",length=0); ax.spines["left"].set_visible(False); panel_label(ax,L,x=-0.3 if L=="a" else -0.06,y=1.14)
save(fig,"FIG_SUPP_n0_treatments")
