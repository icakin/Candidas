"""Fever's cost to the pathogen against fever's cost to the host.
Host: resting metabolic rate rises ~10-13% per degree C, so 37->40 C costs the host
1.33-1.44 fold (compounded). Pathogen: fever_cost from results/tables/fig_values.csv,
the rise in respiration per unit growth over the same 3 C step.
ILLUSTRATIVE COMPARISON, not a like-for-like measurement: the host figure is a rise in
absolute energy expenditure, the pathogen figure a rise in carbon burned per unit biomass.
Both are 'extra resource per unit of function maintained' but they are not the same currency.
"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
INK,MUTED,FAINT,ORANGE,AQUA,VIOLET="#0b0b0b","#52514e","#e8e8e6","#eb6834","#1baf7a","#4a3aa7"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
HOST_LO, HOST_HI = 1.10**3, 1.13**3          # 10-13% per degree, 3 degrees
d=pd.read_csv(f"{C}/results/tables/fig_values.csv")
lab={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
     "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
col={"Clade1":"#9ec4ef","Clade2":"#5c9ae0","Clade3":"#2a78d6","Clade4":"#12508f",
     "para":ORANGE,"Hae":VIOLET,"Duo":AQUA}
d=d[d.Group.isin(lab)].sort_values("fever_cost").reset_index(drop=True)
fig,ax=plt.subplots(figsize=(6.6,3.3))
ax.axvspan(HOST_LO,HOST_HI,color=INK,alpha=0.11,lw=0,zorder=1)
ax.axvline(1.0,color=FAINT,lw=1.0,zorder=1)
for i,r in d.iterrows():
    c=col[r.Group]
    ax.plot([r.fc_lo,r.fc_hi],[i,i],color=c,lw=1.3,solid_capstyle="round",zorder=3)
    ax.scatter(r.fever_cost,i,s=34,color=c,edgecolor="white",linewidth=0.7,zorder=4)
    ax.text(r.fc_hi*1.06,i,f"{r.fever_cost:.2f}×",va="center",fontsize=6.3,color=MUTED)
ax.set_yticks(range(len(d))); ax.set_yticklabels([lab[g] for g in d.Group],fontsize=7)
ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xticks([1,1.5,2,3,5,10]); ax.set_xticklabels(["1","1.5","2","3","5","10"])
ax.set_xlim(0.95,16)
ax.set_xlabel("Cost of a 3 °C fever (37 → 40 °C): fold increase in resource burned per unit of function maintained",labelpad=7)
ax.text(np.sqrt(HOST_LO*HOST_HI),-0.95,"what fever\ncosts the host",ha="center",va="bottom",
        fontsize=6.5,color=INK,linespacing=1.2)
ax.annotate("",xy=(HOST_LO,-0.1),xytext=(HOST_HI,-0.1),arrowprops=dict(arrowstyle="-",color=MUTED,lw=0.6))
ax.text(0.985,0.5,"← cheaper for the pathogen than for the host",transform=ax.transAxes,
        fontsize=6.2,color=MUTED,va="center",ha="right",style="italic")
ax.spines[["top","right","left"]].set_visible(False); ax.tick_params(axis="y",length=0); ax.tick_params(length=2.5)
fig.savefig(f"{FIG}/FIG_fever_asymmetry.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG_fever_asymmetry.pdf",bbox_inches="tight")
print(f"host band {HOST_LO:.2f}-{HOST_HI:.2f}")
print(d[["Group","fever_cost","fc_lo","fc_hi"]].round(2).to_string())
print("\ncheaper for pathogen than host (fever_cost < host_lo):",
      [lab[g] for g in d[d.fever_cost<HOST_LO].Group])
print("within the host band:",[lab[g] for g in d[(d.fever_cost>=HOST_LO)&(d.fever_cost<=HOST_HI)].Group])
print("costlier than the host:",[lab[g] for g in d[d.fever_cost>HOST_HI].Group])
