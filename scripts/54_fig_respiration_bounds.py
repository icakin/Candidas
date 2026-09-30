"""Oxygen consumption after confirmed growth is lost, drawn as BOUNDS rather than as two
near-identical panels. Replaces FIG6_productive_respiration."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
BLUE,ORANGE,INK,MUTED,FAINT="#2a78d6","#eb6834","#0b0b0b","#52514e","#e8e8e6"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
p=pd.read_csv(f"{C}/results/tables/phase2c_productive.csv")
cons=pd.read_csv(f"{C}/results/tables/phase2b_consumption.csv")
# same order as the phylogeny figure
ORD=["Clade1_2068","Clade1_2069","Clade1_2070","Clade3_2074","Clade3_2075","Clade3_2076",
 "Clade2_2071","Clade2_2072","Clade2_2073","Clade4_2077","Clade4_2078","Clade4_2079",
 "Duo_1770","Duo_1771","Hae_1724","Hae_1768","Hae_1769","para_2051","para_2052","para_2053"]
col=lambda i:(BLUE if i.startswith("Clade") else ORANGE)

def pooled(df,c):
    """oxygen-weighted fraction, per temperature, for one clade label"""
    d=df[df.clade==c]
    g=d.groupby("T").apply(lambda x: pd.Series({
        "lo":np.average(x.frac_lo,weights=x.tot),"hi":np.average(x.frac_hi,weights=x.tot)}))
    return g.reset_index()

fig,axes=plt.subplots(1,3,figsize=(8.0,3.3),gridspec_kw=dict(wspace=0.46,width_ratios=[1.05,0.95,0.95]))

# ---- a: pooled fraction vs temperature, the two bounds as one band -----------
ax=axes[0]
for c,colr,lab in (("C. auris",BLUE,"C. auris"),("relatives",ORANGE,"relatives")):
    g=pooled(p,c)
    ax.fill_between(g["T"],g.lo,g.hi,color=colr,alpha=0.22,lw=0)
    ax.plot(g["T"],(g.lo+g.hi)/2,color=colr,lw=1.6,label=lab)
ax.axvspan(36.8,37.2,color=INK,alpha=0.10,lw=0,zorder=0)
ax.axvspan(39.8,40.2,color=INK,alpha=0.06,lw=0,zorder=0)
ax.set_xlabel("temperature (°C)",fontsize=7.4)
ax.set_ylabel("share of O₂ use in\nconfirmed-growth wells",fontsize=7.4)
ax.set_ylim(-0.03,1.05); ax.set_xlim(21.5,45); ax.set_xticks([22,26,30,34,38,42])
ax.legend(frameon=False,fontsize=6.4,loc="lower left")
ax.text(37,0.04,"body",rotation=90,ha="center",va="bottom",fontsize=5.8,color=MUTED)
ax.text(40,0.04,"fever",rotation=90,ha="center",va="bottom",fontsize=5.8,color=MUTED)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.text(-0.36,1.10,"a",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("band spans the two ambiguity bounds",fontsize=6.2,color=MUTED,pad=10,loc="center")

# ---- b: per isolate at 40 C, bound to bound --------------------------------
ax=axes[1]
q=p[p["T"]==40].set_index("isolate")
for yi,iso in enumerate(ORD):
    if iso not in q.index: continue
    r=q.loc[iso]; c=col(iso)
    if abs(r.frac_hi-r.frac_lo)<0.005:
        ax.scatter([r.frac_lo],[yi],s=20,color=c,edgecolor=INK,linewidth=0.5,zorder=4)
        continue
    ax.plot([r.frac_lo,r.frac_hi],[yi,yi],color=c,lw=3.0,alpha=0.35,solid_capstyle="butt",zorder=2)
    ax.scatter([r.frac_lo],[yi],s=16,color=c,zorder=3)
    ax.scatter([r.frac_hi],[yi],s=16,facecolor="white",edgecolor=c,linewidth=0.8,zorder=3)
ax.set_yticks(range(len(ORD)))
ax.set_yticklabels([i+("  ←" if i=="Hae_1724" else "") for i in ORD],fontsize=5.6)
ax.set_ylim(len(ORD)-0.4,-0.9); ax.set_xlim(-0.04,1.08); ax.set_xticks([0,0.5,1])
ax.set_xlabel("share at 40 °C",fontsize=7.4)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y",length=0)
ax.text(-0.62,1.10,"b",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("filled = lower bound, open = upper",fontsize=6.2,color=MUTED,pad=10,loc="center")

# ---- c: what the respiration-only wells actually consume at 40 C ------------
ax=axes[2]
r40=cons[(cons["T"]==40)&(cons.S=="resp_only")]
lab={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
     "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobus."}
GR=[g for g in ["Clade1","Clade2","Clade3","Clade4","para","Duo","Hae"] if (r40.group==g).any()]
for yi,g in enumerate(GR):
    v=r40[r40.group==g].cons.values; c=BLUE if g.startswith("Clade") else ORANGE
    ax.scatter(v,np.full(len(v),yi)+np.random.RandomState(0).uniform(-.12,.12,len(v)),
               s=12,color=c,alpha=0.65,zorder=3)
    ax.plot([np.median(v)]*2,[yi-0.28,yi+0.28],color=c,lw=1.8,zorder=4)
a40=cons[(cons["T"]==40)&(cons.S=="growth")].cons.median()
ax.axvline(a40,color=MUTED,lw=1.0,ls=(0,(3,2)),zorder=1)
ax.text(a40-0.2,-1.15,"median of growing\nwells at 40 °C",fontsize=5.6,color=MUTED,ha="right",linespacing=1.2)
ax.axvline(1.0,color=ORANGE,lw=0.8,ls=":",zorder=1)
ax.text(1.0,-1.15,"1 mg L⁻¹\ngate",fontsize=5.6,color=ORANGE,va="center",ha="center",linespacing=1.2)
ax.set_yticks(range(len(GR))); ax.set_yticklabels([lab[g] for g in GR],fontsize=6.0)
ax.set_ylim(len(GR)-0.5,-1.5)
ax.set_xlabel("O₂ used in the fixed window (mg L⁻¹)",fontsize=7.4)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y",length=0)
ax.text(-0.50,1.10,"c",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("wells respiring without confirmed growth",fontsize=6.2,color=MUTED,pad=10,loc="center")

fig.savefig(f"{FIG}/FIG_respiration_bounds.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG_respiration_bounds.pdf",bbox_inches="tight")
print("saved FIG_respiration_bounds")
for c in ["C. auris","relatives"]:
    g=pooled(p,c); r=g[g["T"]==40]
    print(f"  pooled at 40 C, {c:10s}: {float(r.lo):.2f} to {float(r.hi):.2f}")
q=p[p["T"]==40]
au=q[q.clade=="C. auris"]; print(f"  C. auris per isolate at 40 C: lo {au.frac_lo.min():.2f}-{au.frac_lo.max():.2f}, hi {au.frac_hi.min():.2f}-{au.frac_hi.max():.2f}")
rel=q[q.clade=="relatives"]
print("  relatives per isolate at 40 C:"); 
for _,r in rel.iterrows(): print(f"    {r.isolate:12s} {r.frac_lo:.2f} to {r.frac_hi:.2f}")
ex=rel[rel.isolate!="Hae_1724"]
print(f"  relatives excluding Hae_1724, pooled: {np.average(ex.frac_lo,weights=ex.tot):.2f} to {np.average(ex.frac_hi,weights=ex.tot):.2f}")
