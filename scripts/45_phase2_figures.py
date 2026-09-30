"""PHASE 2 figures: transition (2A) and productive-respiration fraction (2C)."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
BLUE,ORANGE,INK,MUTED,FAINT="#2a78d6","#eb6834","#0b0b0b","#52514e","#e8e8e6"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
t=pd.read_csv(f"{C}/results/tables/phase2a_transition.csv")
p=pd.read_csv(f"{C}/results/tables/phase2c_productive.csv")
cons=pd.read_csv(f"{C}/results/tables/phase2b_consumption.csv")
ORD=[i for i in ["Clade1_2068","Clade1_2069","Clade1_2070","Clade2_2071","Clade2_2072","Clade2_2073",
 "Clade3_2074","Clade3_2075","Clade3_2076","Clade4_2077","Clade4_2078","Clade4_2079",
 "para_2051","para_2052","para_2053","Hae_1724","Hae_1768","Hae_1769","Duo_1770","Duo_1771"]]
col=lambda i:(BLUE if i.startswith("Clade") else ORANGE)
# ================= FIGURE: TRANSITION =========================================
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.2,4.6),gridspec_kw=dict(width_ratios=[1,1.02],wspace=0.42))
for ax,b,lab in ((axA,"R","ambiguous counted as\nrespiration without growth"),
                 (axB,"P","ambiguous counted as\npotentially productive")):
    x=t[t.bound==b].set_index("isolate")
    for yi,iso in enumerate(ORD):
        r=x.loc[iso]; c=col(iso)
        if r.growth_censored:
            ax.plot([44,46.2],[yi,yi],color=c,lw=1.6,solid_capstyle="butt",zorder=2)
            ax.scatter(44,yi,s=26,marker=">",color=c,zorder=3)
        else:
            ax.scatter(r.T_growthloss,yi,s=34,color=c,edgecolor="white",linewidth=0.6,zorder=3)
            ax.plot([r.T_growthloss,46.2],[yi,yi],color=c,lw=1.0,alpha=0.30,zorder=2)
    ax.axvspan(37,40,color=INK,alpha=0.08,lw=0,zorder=1)
    ax.axvline(44,color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1)
    ax.set_yticks(range(len(ORD)))
    ax.set_yticklabels(ORD if ax is axA else [],fontsize=6.0)
    ax.invert_yaxis(); ax.set_xlim(33,47.2); ax.set_xticks([34,36,38,40,42,44])
    ax.set_xlabel("temperature of confirmed-growth loss (°C)",fontsize=7.2)
    ax.set_title(lab,fontsize=6.9,color=INK,pad=4,linespacing=1.3)
    for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y",length=0); ax.tick_params(axis="x",length=2.5)
    ax.text(44.15,len(ORD)-0.4,"assay\nceiling",fontsize=5.6,color=MUTED,va="top",linespacing=1.15)
axA.text(-0.62,1.06,"a",transform=axA.transAxes,fontweight="bold",fontsize=9)
axB.text(-0.10,1.06,"b",transform=axB.transAxes,fontweight="bold",fontsize=9)
h=[Line2D([],[],marker="o",color=BLUE,ls="",ms=5,label="C. auris (12 isolates)"),
   Line2D([],[],marker="o",color=ORANGE,ls="",ms=5,label="relatives (8 isolates)"),
   Line2D([],[],marker=">",color=MUTED,ls="",ms=5,label="right-censored: still growing at 44 °C")]
fig.legend(handles=h,frameon=False,fontsize=6.5,ncol=3,loc="upper center",bbox_to_anchor=(0.5,1.02))
fig.text(0.5,-0.02,"Shaded band 37-40 °C. Bars extend to the assay ceiling: respiration remained detectable "
         "in every isolate at every temperature,\nso the respiration-without-growth interval is right-censored "
         "throughout and can only be bounded below.",ha="center",fontsize=6.1,color=MUTED,linespacing=1.4)
fig.savefig(f"{FIG}/FIG_phase2_transition.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG_phase2_transition.pdf",bbox_inches="tight"); plt.close(fig)
# ================= FIGURE: PRODUCTIVE FRACTION ================================
fig,axes=plt.subplots(1,3,figsize=(7.6,3.4),gridspec_kw=dict(wspace=0.52,width_ratios=[1,1,0.95]))
Ts=[34,36,38,40,42,44]
for ax,colname,lab in ((axes[0],"frac_lo","lower bound\n(ambiguous non-productive)"),
                       (axes[1],"frac_hi","upper bound\n(ambiguous productive)")):
    for iso in ORD:
        s=p[p.isolate==iso].set_index("T").reindex(Ts)
        ax.plot(Ts,s[colname].values,color=col(iso),lw=1.0,alpha=0.55,
                marker="o",ms=2.6,zorder=2)
    for cl,c in (("C. auris",BLUE),("relatives",ORANGE)):
        m=p[p.clade==cl].pivot_table(index="T",values=colname,aggfunc="median").reindex(Ts)
        ax.plot(Ts,m[colname].values,color=c,lw=2.2,zorder=3)
    ax.axvspan(37,40,color=INK,alpha=0.08,lw=0,zorder=1)
    ax.set_ylim(-0.03,1.05); ax.set_xticks(Ts); ax.set_xlabel("temperature (°C)",fontsize=7.2)
    ax.set_title(lab,fontsize=6.9,pad=4,linespacing=1.3)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.5)
axes[0].set_ylabel("fraction of O₂ consumption in\nconfirmed-growth wells",fontsize=7.2)
axes[1].set_yticklabels([])
ax=axes[2]
r40=cons[(cons["T"]==40)&(cons.S=="resp_only")]
grp=["Clade1","Clade2","para","Duo","Hae"]
NMS={"Clade1":"C. auris I","Clade2":"C. auris II","para":"C. parapsilosis","Duo":"C. duobus.","Hae":"C. haemulonii"}
a40=cons[(cons["T"]==40)&(cons.S=="growth")&(cons.clade=="C. auris")].cons.median()
for yi,g in enumerate(grp):
    s=r40[r40.group==g]
    if not len(s): continue
    c=BLUE if g.startswith("Clade") else ORANGE
    ax.scatter(s.cons,[yi]*len(s),s=20,color=c,alpha=0.75,edgecolor="white",linewidth=0.4,zorder=3)
    ax.plot([s.cons.median()]*2,[yi-0.3,yi+0.3],color=c,lw=1.8,zorder=4)
ax.axvline(a40,color=MUTED,lw=1.0,ls=(0,(3,2)),zorder=1)
ax.text(a40+0.15,-0.55,"C. auris confirmed-growth at 40 °C",fontsize=5.6,color=MUTED,ha="right",va="center")
ax.axvline(1.0,color=ORANGE,lw=0.8,ls=":",zorder=1)
ax.text(1.2,-0.35,"respiration\ngate",fontsize=5.6,color=ORANGE,va="center",linespacing=1.2)
ax.set_yticks(range(len(grp))); ax.set_yticklabels([NMS[g] for g in grp],fontsize=6.3)
ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
ax.invert_yaxis(); ax.set_xlabel("O₂ consumed in the fixed window (mg L⁻¹)\nrespiration-only wells at 40 °C",fontsize=6.6)
ax.set_xlim(0,9.8)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y",length=0); ax.tick_params(axis="x",length=2.5)
for a,l in ((axes[0],"a"),(axes[1],"b"),(axes[2],"c")):
    a.text(-0.20 if l!="c" else -0.10,1.12,l,transform=a.transAxes,fontweight="bold",fontsize=9)
fig.savefig(f"{FIG}/FIG_phase2_productive.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG_phase2_productive.pdf",bbox_inches="tight")
print("saved FIG_phase2_transition and FIG_phase2_productive")
