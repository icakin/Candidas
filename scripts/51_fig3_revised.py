"""FIGURE 3 (revised) - thermal displacement as a RELATIVE result.
Leads with relative displacement and invariant ordering, not the sign of headroom.
Reports sensitivity across 1.5x, 2x, 3x and flags 3x as closest to observed growth loss.
Uses the ORIGINAL expert-selected quantitative fits (unchanged), per the plan."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
from matplotlib.lines import Line2D
RAMP=OI; NAME=SHORT
dr=pd.read_csv(f"{C}/results/tables/thermal_headroom_draws.csv")
MULT=[("T_1.5xcost","1.5×"),("T_2xcost","2×"),("T_3xcost","3×")]
ORD=["Clade4","Clade3","Clade1","Clade2","para","Hae","Duo"]
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.3,3.6),gridspec_kw=dict(width_ratios=[1.25,1],wspace=0.40))
# ---- (a) displacement RELATIVE to C. haemulonii, at all three multipliers ----
off={"1.5×":-0.26,"2×":0.0,"3×":0.26}; mk={"1.5×":"^","2×":"o","3×":"s"}
for col,lab in MULT:
    ref=dr[dr.Group=="Hae"][col].dropna().values
    for yi,g in enumerate(ORD):
        v=dr[dr.Group==g][col].dropna().values
        n=min(len(v),len(ref)); d=v[:n]-ref[:n]
        med=np.median(d); lo,hi=np.percentile(d,[2.5,97.5])
        y=yi+off[lab]
        axA.plot([lo,hi],[y,y],color=RAMP[g],lw=1.0,alpha=0.85,solid_capstyle="round",zorder=2)
        axA.scatter(med,y,s=20,color=RAMP[g],marker=mk[lab],edgecolor="white",linewidth=0.5,zorder=3)
axA.axvline(0,color=INK,lw=0.9,zorder=1)
axA.set_yticks(range(len(ORD))); axA.set_yticklabels([NAME[g] for g in ORD],fontsize=7)
axA.invert_yaxis(); axA.set_xlabel("thermal displacement relative to C. haemulonii (°C)",fontsize=7.2)
axA.set_xlim(-4,14)
for sp in ("top","right","left"): axA.spines[sp].set_visible(False)
axA.tick_params(axis="y",length=0); axA.tick_params(axis="x",length=2.5)
h=[Line2D([],[],marker=mk[l],color=MUTED,ls="",ms=4.5,label=f"{l} cost multiplier") for _,l in MULT]
axA.legend(handles=h,frameon=False,fontsize=6.2,loc="lower right")
axA.text(-0.42,1.06,"a",transform=axA.transAxes,fontweight="bold",fontsize=9)
# ---- (b) the ordering is invariant across multipliers ------------------------
for col,lab in MULT:
    xs=[dr[dr.Group==g][col].median() for g in ORD]
    axB.plot(xs,range(len(ORD)),color=MUTED,lw=0.7,alpha=0.5,zorder=1)
    for yi,g in enumerate(ORD):
        axB.scatter(xs[yi],yi,s=22,color=RAMP[g],marker=mk[lab],edgecolor="white",linewidth=0.5,zorder=3)
axB.axvspan(36.8,37.2,color=INK,alpha=0.10,lw=0,zorder=0)
axB.axvspan(39.8,40.2,color=INK,alpha=0.06,lw=0,zorder=0)
axB.axvline(44,color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1)
axB.set_yticks(range(len(ORD))); axB.set_yticklabels([])
axB.invert_yaxis(); axB.set_xlabel("temperature at which respiration per unit\ngrowth reaches the multiple (°C)",fontsize=7.2)
axB.set_xlim(33,49)
axB.text(37,-0.75,"body",ha="center",fontsize=6.2,color=MUTED)
axB.text(40,-0.75,"fever",ha="center",fontsize=6.2,color=MUTED)
axB.text(44.2,-0.75,"assay ceiling",ha="center",fontsize=6.2,color=MUTED)
for sp in ("top","right","left"): axB.spines[sp].set_visible(False)
axB.tick_params(axis="y",length=0); axB.tick_params(axis="x",length=2.5)
axB.text(-0.10,1.06,"b",transform=axB.transAxes,fontweight="bold",fontsize=9)
save(fig,"FIG3_revised_displacement")
# numbers for the caption
print("=== Clade IV minus C. haemulonii, by multiplier ===")
for col,lab in MULT:
    a=dr[dr.Group=="Clade4"][col].dropna().values; r=dr[dr.Group=="Hae"][col].dropna().values
    n=min(len(a),len(r)); d=a[:n]-r[:n]
    print(f"  {lab}: {np.median(d):.2f} °C  [{np.percentile(d,2.5):.2f}, {np.percentile(d,97.5):.2f}]  P(>0)={np.mean(d>0):.3f}")
print("\n=== ordering at each multiplier (warmest first) ===")
for col,lab in MULT:
    o=sorted(ORD,key=lambda g:-dr[dr.Group==g][col].median())
    print(f"  {lab}: "+" > ".join(NAME[g] for g in o))
