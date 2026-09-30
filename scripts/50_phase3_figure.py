"""PHASE 3 figure: does metabolism at 22-38 °C add predictive information beyond growth?"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
BLUE,ORANGE,INK,MUTED,FAINT,GREEN="#2a78d6","#eb6834","#0b0b0b","#52514e","#e8e8e6","#1baf7a"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
r=pd.read_csv(f"{C}/results/tables/phase3_loio_results.csv")
L=pd.read_csv(f"{C}/results/tables/phase3_limit_per_isolate.csv")
r=r[r.model.notna()&(r.model!="-")]
OUT=[("g40_lo","40 °C\nlower"),("g40_hi","40 °C\nupper"),("g42_lo","42 °C\nlower"),
     ("g42_hi","42 °C\nupper"),("g44_hi","44 °C\nupper")]
COL={"M1_growth":BLUE,"M2_metab":ORANGE,"M3_combined":GREEN}
LAB={"M1_growth":"growth only (r₃₈, Topt)","M2_metab":"metabolism only (log K/r at 38, slope)",
     "M3_combined":"combined (best of each)"}
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.4,3.5),gridspec_kw=dict(width_ratios=[1.15,1],wspace=0.34))
# (a) delta LPD vs growth-only
w=0.26
for k,m in enumerate(["M2_metab","M3_combined"]):
    d=[]
    for oc,_ in OUT:
        s=r[r.outcome==oc].set_index("model")
        d.append(s.lpd[m]-s.lpd["M1_growth"])
    axA.bar(np.arange(len(OUT))+(k-0.5)*w,d,width=w,color=COL[m],zorder=2,label=LAB[m])
axA.axhline(0,color=INK,lw=0.9,zorder=3)
axA.set_xticks(range(len(OUT))); axA.set_xticklabels([l for _,l in OUT],fontsize=6.6)
axA.set_ylabel("Δ held-out log predictive score\nvs growth-only  (>0 = metabolism helps)",fontsize=7)
axA.set_xlabel("outcome: confirmed growth at…",fontsize=7.2)
axA.set_ylim(-0.22,0.12)
axA.text(0.02,0.94,"metabolism HELPS",transform=axA.transAxes,fontsize=6.2,color=MUTED,va="top")
axA.text(0.02,0.03,"metabolism HURTS",transform=axA.transAxes,fontsize=6.2,color=MUTED,va="bottom")
axA.legend(frameon=False,fontsize=6.2,loc="lower center",bbox_to_anchor=(0.55,0.02))
for sp in ("top","right"): axA.spines[sp].set_visible(False)
axA.tick_params(length=2.5)
# (b) observed vs predicted upper thermal limit
L["cenb"]=L.cen==1
for m in ("M1_growth","M2_metab"):
    o=L[~L.cenb]
    axB.scatter(o.lim,o[f"pred_{m}"],s=26,color=COL[m],alpha=0.85,edgecolor="white",linewidth=0.5,
                zorder=3,label=LAB[m].split(" (")[0])
    c=L[L.cenb]
    axB.scatter([44.6]*len(c),c[f"pred_{m}"],s=26,marker=">",color=COL[m],alpha=0.85,zorder=3)
axB.plot([35,46],[35,46],color=MUTED,lw=0.9,ls=(0,(3,2)),zorder=1)
axB.axvline(44.2,color=FAINT,lw=6,zorder=0)
axB.text(44.35,45.6,">44 censored",fontsize=5.9,color=MUTED,ha="center",rotation=90,va="top")
axB.set_xlim(35,46); axB.set_ylim(35,46)
axB.set_xlabel("observed temperature of confirmed-growth loss (°C)",fontsize=7.2)
axB.set_ylabel("predicted from 22–38 °C only (°C)",fontsize=7.2)
axB.legend(frameon=False,fontsize=6.2,loc="upper left")
for sp in ("top","right"): axB.spines[sp].set_visible(False)
axB.tick_params(length=2.5)
mae=pd.read_csv(f"{C}/results/tables/phase3_limit_results.csv").set_index("model")
axB.text(0.03,0.70,f"MAE (uncensored, n=18)\n  growth only      {mae.MAE_uncensored['M1_growth']:.2f} °C\n"
                   f"  metabolism only {mae.MAE_uncensored['M2_metab']:.2f} °C\n"
                   f"  combined         {mae.MAE_uncensored['M3_combined']:.2f} °C",
         transform=axB.transAxes,fontsize=6.0,ha="left",color=INK,linespacing=1.45)
for a,l in ((axA,"a"),(axB,"b")): a.text(-0.17,1.05,l,transform=a.transAxes,fontweight="bold",fontsize=9)
fig.savefig(f"{FIG}/FIG_phase3_prediction.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG_phase3_prediction.pdf",bbox_inches="tight")
print("saved FIG_phase3_prediction")
