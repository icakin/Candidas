"""Supplementary figure: plate-reader OD validation of the oxygen-derived growth calls.
Four isolates × 37, 40, 42 °C in 24-well plates (data/od_validation). Panel a: blank-corrected OD600
per well; panel b: fold change in OD over 16 h against the likelihood-model call for the same isolate
and temperature in the sealed PreSens vials. Writes results/tables/od_validation.csv."""
import os, numpy as np, pandas as pd
from style_common import *
from matplotlib.lines import Line2D
D=f"{C}/data/od_validation"; TEMPS=[37,40,42]
pm=pd.read_csv(f"{D}/plate_map.csv").set_index("well")["isolate"]
ISO=["Clade4_2077","Clade2_2073","para_2052","Hae_1768"]
lik=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv")
lik["grow"]=(lik.p_boot<0.01)&(lik.drawdown>=1)
rows=[]; curves={}
for T in TEMPS:
    d=pd.read_csv(f"{D}/od_{T}.csv",index_col=0); t=d.index.values
    blank=d[[w for w in d.columns if pm[w]=="blank"]].mean(axis=1)
    for w in d.columns:
        iso=pm[w]
        if iso in("blank","water"): continue
        y=(d[w]-blank).values; curves[(T,w)]=(t,y)
        # early rate: log-linear slope over the first 3 h; late rate: 8 to 16 h
        e=t<=3; l=t>=8
        rows.append(dict(T=T,well=w,isolate=iso,od0=y[0],od_end=y[-1],fold=y[-1]/y[0],
                         r_early=np.polyfit(t[e],np.log(np.clip(y[e],1e-3,None)),1)[0],
                         r_late=np.polyfit(t[l],np.log(np.clip(y[l],1e-3,None)),1)[0]))
res=pd.DataFrame(rows)
summ=res.groupby(["isolate","T"]).agg(od_end=("od_end","mean"),od_end_sd=("od_end","std"),fold=("fold","mean"),fold_sd=("fold","std"),r_early=("r_early","mean"),r_late=("r_late","mean"),n=("well","size")).reset_index()
ox=[]
for iso in ISO:
    for T in TEMPS:
        Tv=38 if T==37 else T   # vials were run on a 2 °C grid; 37 °C is compared with the 38 °C vials
        s=lik[(lik.otu_name==iso)&(lik["T"]==Tv)]
        ox.append(dict(isolate=iso,T=T,vial_T=Tv,vial_growing=int(s.grow.sum()),vial_n=len(s),vial_r=s[s.grow].r_hat.median() if s.grow.any() else np.nan))
summ=summ.merge(pd.DataFrame(ox),on=["isolate","T"]); summ.to_csv(f"{C}/results/tables/od_validation.csv",index=False)
print(summ.round(3).to_string(index=False))
# ---------- figure ----------
fig=plt.figure(figsize=(7.2,4.6)); gs=fig.add_gridspec(2,3,height_ratios=[1,0.75],hspace=0.5,wspace=0.28)
axs=[fig.add_subplot(gs[0,i]) for i in range(3)]
for ax,T in zip(axs,TEMPS):
    for (TT,w),(t,y) in curves.items():
        if TT!=T: continue
        ax.plot(t,y,color=col(pm[w]),lw=0.9,alpha=0.9)
    ax.set_title(f"{T} °C",fontsize=8,pad=4); ax.set_xlabel("time (h)"); ax.set_xticks([0,4,8,12,16]); ax.set_ylim(0,1.0)
    if T!=37: ax.tick_params(labelleft=False)
axs[0].set_ylabel("OD$_{600}$, blank-corrected")
h=[Line2D([],[],color=col(i),lw=1.6,label=iso_label(i)) for i in ISO]
axs[2].legend(handles=h,loc="upper left",fontsize=6,handlelength=1.4); panel_label(axs[0],"a",x=-0.28)
# b: fold change vs vial call
axb=fig.add_subplot(gs[1,:])
x=np.arange(len(ISO)); wdt=0.25
for k,T in enumerate(TEMPS):
    s=summ[summ["T"]==T].set_index("isolate").loc[ISO]
    xs=x+(k-1)*wdt
    axb.bar(xs,s.fold,wdt,color=[col(i) for i in ISO],alpha=[1.0,0.7,0.45][k],edgecolor="white",lw=0.5)
    axb.errorbar(xs,s.fold,yerr=s.fold_sd,fmt="none",ecolor=INK,elinewidth=0.6,capsize=1.5)
    for xi,(iso,r) in zip(xs,s.iterrows()):
        lab=f"{int(r.vial_growing)}/{int(r.vial_n)}"; axb.text(xi,r.fold+r.fold_sd+0.8,lab,ha="center",fontsize=5.8,color=INK if r.vial_growing>=3 else RED)
axb.axhline(1,color=AMB,lw=0.7); axb.set_yscale("log"); axb.set_ylim(0.8,60)
from matplotlib.ticker import FixedLocator,FixedFormatter,NullFormatter
axb.yaxis.set_major_locator(FixedLocator([1,2,5,10,20,50])); axb.yaxis.set_major_formatter(FixedFormatter(["1×","2×","5×","10×","20×","50×"])); axb.yaxis.set_minor_locator(FixedLocator([])); axb.yaxis.set_minor_formatter(NullFormatter())
axb.set_xticks(x); axb.set_xticklabels([iso_label(i) for i in ISO],fontsize=7,style="italic"); axb.set_ylabel("OD$_{600}$ fold change, 0 to 16 h")
panel_label(axb,"b",x=-0.085,y=1.08)
from matplotlib.patches import Patch
axb.legend(handles=[Patch(color=MUTED,alpha=a,label=f"{T} °C") for T,a in zip(TEMPS,[1.0,0.7,0.45])],loc="upper right",fontsize=6,ncol=3,handlelength=1.2,columnspacing=0.8)
save(fig,"FIG_SUPP_od_validation")
