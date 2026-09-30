"""Supplementary: does the choice of fit window (manual vs automatic) move any conclusion?
Both arms are the same brms model, refit on the same traces; only the window rule differs."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
d=pd.read_csv(f"{C}/reprocess/refit_growth_comparison.csv")
NAME=SHORT
d["name"]=d.Group.map(NAME)
d=d.sort_values("Topt_legacy",ascending=False).reset_index(drop=True)
col=lambda g: OI[g]

fig,(axA,axB)=plt.subplots(1,2,figsize=(7.4,3.1),gridspec_kw=dict(wspace=0.60,width_ratios=[1,0.8]))
# ---- a: shift in the fitted growth optimum -----------------------------------
for i,r in d.iterrows():
    c=col(r.Group)
    axA.plot([r.diff_lo,r.diff_hi],[i,i],color=c,lw=1.1,alpha=0.85,solid_capstyle="round",zorder=2)
    axA.scatter(r.diff_med,i,s=22,color=c,edgecolor="white",linewidth=0.5,zorder=3)
axA.axvline(0,color=INK,lw=0.9,zorder=1)
axA.set_yticks(range(len(d))); axA.set_yticklabels(d.name,fontsize=6.6)
axA.set_ylim(len(d)-0.4,-0.8); axA.invert_yaxis()
axA.set_xlabel("shift in fitted growth optimum,\nautomatic − manual windows (°C)",fontsize=7.2)
for sp in ("top","right","left"): axA.spines[sp].set_visible(False)
axA.tick_params(axis="y",length=0)
axA.text(-0.55,1.06,"a",transform=axA.transAxes,fontweight="bold",fontsize=9)
axA.set_title("every interval includes zero",fontsize=6.3,color=MUTED,pad=8,loc="center")
# ---- b: the optima themselves, ordering preserved ---------------------------
for i,r in d.iterrows():
    c=col(r.Group)
    axB.plot([0,1],[r.Topt_legacy,r.Topt_auto],color=c,lw=1.0,alpha=0.8,zorder=2)
    axB.scatter([0,1],[r.Topt_legacy,r.Topt_auto],s=20,color=c,edgecolor="white",linewidth=0.5,zorder=3)
    axB.text(-0.06,r.Topt_legacy,r["name"],fontsize=6.0,ha="right",va="center",color=INK)
axB.set_xlim(-1.45,1.25); axB.set_xticks([0,1]); axB.set_xticklabels(["manual","automatic"],fontsize=7.0)
axB.set_ylabel("fitted growth optimum (°C)",fontsize=7.2)
for sp in ("top","right"): axB.spines[sp].set_visible(False)
axB.tick_params(axis="x",length=2.5)
axB.text(-0.72,1.06,"b",transform=axB.transAxes,fontweight="bold",fontsize=9)
axB.set_title("order preserved, except the two groups\nwhose optima are 0.2 °C apart",fontsize=6.3,color=MUTED,pad=8,loc="center")
save(fig,"FIG_window_sensitivity")
print("saved FIG_window_sensitivity")
print("  Topt shift range: %.2f to %.2f C" % (d.diff_med.min(),d.diff_med.max()))
print("  all intervals span zero:", bool(((d.diff_lo<0)&(d.diff_hi>0)).all()))
print("  P(shift>0) range: %.2f to %.2f" % (d.P_diff_gt0.min(),d.P_diff_gt0.max()))
print("  |dE| max: %.2f eV   |dTh| max: %.2f K (%s)" % (d.dE_med.abs().max(),d.dTh_med.abs().max(),
      d.loc[d.dTh_med.abs().idxmax(),"name"]))
print("  |dTh| excluding that group: %.2f K" % d.drop(d.dTh_med.abs().idxmax()).dTh_med.abs().max())
o1=list(d.sort_values("Topt_legacy",ascending=False).Group); o2=list(d.sort_values("Topt_auto",ascending=False).Group)
print("  order manual   :",o1); print("  order automatic:",o2); print("  identical:",o1==o2)
