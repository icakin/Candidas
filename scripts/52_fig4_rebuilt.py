"""FIGURE 4 (rebuilt) - upper thermal limits from the FROZEN automatic classifier,
aligned to phylogeny. Isolate-level. Censored isolates are shown as >44, never as 46."""
import os, re, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
BLUE,ORANGE,INK,MUTED,FAINT,AMB="#2a78d6","#eb6834","#0b0b0b","#52514e","#e8e8e6","#b9b7b2"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
t=pd.read_csv(f"{C}/results/tables/phase2a_transition.csv")
R=t[t.bound=="R"].set_index("isolate"); P=t[t.bound=="P"].set_index("isolate")
# phylogeny: depths from phylo/trees/pg_rooted.nwk, taxa in tree order
TREE_ORDER=["Clade1","Clade3","Clade2","Clade4","Duo","Hae","para"]
NAME={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
      "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
nwk=open(f"{C}/phylo/trees/pg_rooted.nwk").read()
KEY={"Clade1":"auris_cladeI","Clade2":"auris_cladeII","Clade3":"auris_cladeIII","Clade4":"auris_cladeIV",
     "Duo":"duobushaemulonii","Hae":"haemulonii","para":"parapsilosis"}
rows=[]
for g in TREE_ORDER:
    for iso in sorted(t[t.group==g].isolate.unique()): rows.append((g,iso))
fig,(axT,axL)=plt.subplots(1,2,figsize=(7.8,4.5),gridspec_kw=dict(width_ratios=[0.52,1],wspace=0.30))
# ---- left: schematic cladogram from the rooted topology ----------------------
# topology: (((I,III),II),IV) | ((Duo,pseudo),Hae) | para  (outgroup)
ypos={iso:i for i,(g,iso) in enumerate(rows)}
tax_y={g:np.mean([ypos[i] for gg,i in rows if gg==g]) for g in TREE_ORDER}
def hline(x0,x1,y,**kw): axT.plot([x0,x1],[y,y],color=MUTED,lw=0.9,**kw)
def vline(x,y0,y1,**kw): axT.plot([x,x],[y0,y1],color=MUTED,lw=0.9,**kw)
for g in TREE_ORDER:
    for gg,iso in rows:
        if gg==g: hline(0.86,1.0,ypos[iso])
    hline(0.80,0.86,tax_y[g])
    vline(0.86,min(ypos[i] for gg,i in rows if gg==g),max(ypos[i] for gg,i in rows if gg==g))
y_I_III=np.mean([tax_y["Clade1"],tax_y["Clade3"]]); vline(0.72,tax_y["Clade1"],tax_y["Clade3"]); hline(0.72,0.80,tax_y["Clade1"]); hline(0.72,0.80,tax_y["Clade3"])
y_iii=np.mean([y_I_III,tax_y["Clade2"]]); vline(0.64,y_I_III,tax_y["Clade2"]); hline(0.64,0.72,y_I_III); hline(0.64,0.80,tax_y["Clade2"])
y_auris=np.mean([y_iii,tax_y["Clade4"]]); vline(0.54,y_iii,tax_y["Clade4"]); hline(0.54,0.64,y_iii); hline(0.54,0.80,tax_y["Clade4"])
y_rel=np.mean([tax_y["Duo"],tax_y["Hae"]]); vline(0.54,tax_y["Duo"],tax_y["Hae"]); hline(0.54,0.80,tax_y["Duo"]); hline(0.54,0.80,tax_y["Hae"])
y_core=np.mean([y_auris,y_rel]); vline(0.30,y_auris,y_rel); hline(0.30,0.54,y_auris); hline(0.30,0.54,y_rel)
y_root=np.mean([y_core,tax_y["para"]]); vline(0.10,y_core,tax_y["para"]); hline(0.10,0.30,y_core); hline(0.10,0.80,tax_y["para"])
hline(0.02,0.10,y_root)
for g in TREE_ORDER:
    axT.text(1.0,tax_y[g],"  "+NAME[g],fontsize=6.2,style="italic",ha="left",va="center",color=INK)
axT.set_xlim(0,2.35); axT.set_ylim(len(rows)-0.4,-0.9); axT.axis("off")
axT.text(0.0,-0.55,"phylogenomic topology",fontsize=6.4,color=MUTED)
axT.text(-0.06,1.04,"a",transform=axT.transAxes,fontweight="bold",fontsize=9)
# ---- right: confirmed upper thermal limit per isolate ------------------------
for (g,iso) in rows:
    y=ypos[iso]; c=BLUE if g.startswith("Clade") else ORANGE
    r=R.loc[iso]; p=P.loc[iso]
    lo = 44.0 if r.growth_censored else r.T_growthloss     # conservative bound
    hi = 44.0 if p.growth_censored else p.T_growthloss     # permissive bound
    axL.plot([lo,max(hi,lo)],[y,y],color=AMB,lw=4.0,solid_capstyle="butt",zorder=2)
    if r.growth_censored:
        axL.scatter(44,y,s=30,marker=">",color=c,zorder=4)
        axL.plot([44,45.6],[y,y],color=c,lw=1.6,zorder=3)
    else:
        axL.scatter(lo,y,s=30,color=c,edgecolor="white",linewidth=0.6,zorder=4)
    if p.growth_censored and not r.growth_censored:
        axL.scatter(44,y,s=22,marker=">",facecolor="white",edgecolor=c,linewidth=0.8,zorder=4)
axL.axvspan(36.8,37.2,color=INK,alpha=0.10,lw=0,zorder=0)
axL.axvspan(39.8,40.2,color=INK,alpha=0.06,lw=0,zorder=0)
axL.axvline(44,color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1)
axL.set_yticks(range(len(rows)))
axL.set_yticklabels([iso+("  ←" if iso=="Hae_1724" else "") for g,iso in rows],fontsize=6.2)
axL.set_ylim(len(rows)-0.4,-0.9)
axL.set_xlim(34.5,46.2); axL.set_xticks([36,38,40,42,44])
axL.set_xlabel("temperature of confirmed-growth loss (°C)",fontsize=7.4)
axL.text(37,-0.85,"body",ha="center",fontsize=6.2,color=MUTED)
axL.text(40,-0.85,"fever",ha="center",fontsize=6.2,color=MUTED)
axL.text(44.15,-0.85,"assay ceiling",ha="left",fontsize=6.0,color=MUTED)
for sp in ("top","right","left"): axL.spines[sp].set_visible(False)
axL.tick_params(axis="y",length=0); axL.tick_params(axis="x",length=2.5)
axL.text(-0.22,1.04,"b",transform=axL.transAxes,fontweight="bold",fontsize=9)
h=[Line2D([],[],marker="o",color=BLUE,ls="",ms=5,label="C. auris isolate"),
   Line2D([],[],marker="o",color=ORANGE,ls="",ms=5,label="relative isolate"),
   Line2D([],[],color=AMB,lw=4,label="ambiguity interval (conservative → permissive bound)"),
   Line2D([],[],marker=">",color=MUTED,ls="",ms=5,label="right-censored: still growing at 44 °C")]
fig.legend(handles=h,frameon=False,fontsize=6.3,ncol=2,loc="upper center",bbox_to_anchor=(0.55,1.07))
fig.text(0.55,-0.04,"Limits are derived from the frozen blind classifier, not from the manual trimming decision.\n"
 "Censored isolates are shown as >44 °C and are never assigned a measured value.",
 ha="center",fontsize=6.2,color=MUTED,linespacing=1.4)
fig.savefig(f"{FIG}/FIG4_rebuilt_limits.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG4_rebuilt_limits.pdf",bbox_inches="tight")
print("saved FIG4_rebuilt_limits")
