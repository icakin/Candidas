"""FIGURE 4 (likelihood model). Every well scored by the per-well likelihood
model (reprocess/lik_model.py): growth when the bootstrap p of the exponential-vs-linear
likelihood ratio is < 0.01; respiration without growth otherwise; no detectable respiration
when total drawdown < 1.0 mg/L. Panel (b) shows the fitted growth rate r_hat itself.
Reads reprocess/lik_results_labelled.csv; nothing is recomputed here."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
from matplotlib.patches import Patch
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
P_GROW=0.01; DRAWDOWN_MIN=1.0
d=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv")
d["S"]=np.where(d.drawdown<DRAWDOWN_MIN,"no_resp",np.where(d.p_boot<P_GROW,"growth","resp_only"))
ORDER=["growth","resp_only","no_resp"]; COL={"growth":GROW,"resp_only":RESP,"no_resp":NORESP}
LAB={"growth":"growth detected (p < 0.01)","resp_only":"respiration without detectable growth","no_resp":"no detectable respiration"}
TAXA=GROUPS; TNAME=SHORT
Ts=sorted(d["T"].unique())
fig=plt.figure(figsize=(7.2,9.0))
gs=fig.add_gridspec(3,1,height_ratios=[1.05,1.62,0.70],hspace=0.46)
# ---- (a) state proportions ---------------------------------------------------
gsa=gs[0].subgridspec(1,7,wspace=0.12)
for i,t in enumerate(TAXA):
    ax=fig.add_subplot(gsa[0,i]); s=d[d.group==t]
    frac=(s.groupby(["T","S"]).size().unstack(fill_value=0).reindex(columns=ORDER,fill_value=0))
    frac=frac.div(frac.sum(1),axis=0); bot=np.zeros(len(frac))
    for st in ORDER:
        ax.bar(range(len(frac)),frac[st].values,bottom=bot,color=COL[st],width=0.92,lw=0.3,edgecolor="white",zorder=2)
        bot+=frac[st].values
    ax.axvspan(Ts.index(38)-0.5,Ts.index(40)+0.5,color=INK,alpha=0.07,lw=0,zorder=1)
    ax.set_xticks(range(len(Ts))); ax.set_xticklabels([str(x) if x in (22,30,38,44) else "" for x in Ts],fontsize=5.8)
    ax.set_ylim(0,1); ax.set_xlim(-0.6,len(Ts)-0.4); ax.set_title(TNAME[t],fontsize=6.3,style="italic",pad=3)
    if i==0: ax.set_ylabel("fraction of vials",fontsize=7)
    else: ax.set_yticklabels([])
    ax.tick_params(length=2)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
axA=fig.axes[0]; axA.text(-0.42,1.30,"a",transform=axA.transAxes,fontweight="bold",fontsize=9)
axA.text(0.0,1.16,"shaded band: 37-40 °C, host to fever",transform=axA.transAxes,fontsize=6.2,color=MUTED)
# ---- (b) per-isolate growth rate, every well -----------------------------------
axB=fig.add_subplot(gs[1])
isos=[i for t in TAXA for i in sorted(d[d.group==t].otu_name.unique())]
ypos={iso:k for k,iso in enumerate(isos)}
rmax=np.nanpercentile(d[d.S=="growth"].r_hat,99)
cmap=LinearSegmentedColormap.from_list("grow",["#dfe7f3",GROW]); norm=LogNorm(vmin=0.02,vmax=max(rmax,0.5))
for iso in isos:
    for T in Ts:
        w=d[(d.otu_name==iso)&(d["T"]==T)].sort_values("Replicate")
        x0=Ts.index(T)-0.46
        for j,(_,r) in enumerate(w.iterrows()):
            if r.S=="growth": c=cmap(norm(max(r.r_hat,0.02)))
            elif r.S=="resp_only": c=RESP
            else: c=NORESP
            axB.add_patch(plt.Rectangle((x0+j*0.184,ypos[iso]-0.42),0.184,0.84,facecolor=c,edgecolor="white",lw=0.3))
axB.axvspan(Ts.index(38)-0.5,Ts.index(40)+0.5,color=INK,alpha=0.07,lw=0,zorder=0)
for t in TAXA[:-1]:
    last=max(ypos[i] for i in isos if i.startswith(t) or (t=="para" and i.startswith("para")) or (t=="Hae" and i.startswith("Hae")))
    axB.axhline(last+0.5,color=MUTED,lw=0.5)
axB.set_xlim(-0.6,len(Ts)-0.4); axB.set_ylim(len(isos)-0.5,-0.5)
axB.set_yticks(range(len(isos))); axB.set_yticklabels([f"{SHORT[i.split('_')[0]]}  {i.split('_')[1]}" for i in isos],fontsize=6.2,style="italic")
for lab,iso in zip(axB.get_yticklabels(),isos): lab.set_color(col(iso))
axB.set_xticks(range(len(Ts))); axB.set_xticklabels([str(x) for x in Ts],fontsize=6.5)
axB.set_xlabel("temperature (°C)",fontsize=7.5); axB.tick_params(length=0)
for sp in ("top","right","left","bottom"): axB.spines[sp].set_visible(False)
axB.text(-0.175,1.04,"b",transform=axB.transAxes,fontweight="bold",fontsize=9)
axB.text(0.0,1.015,"each row is one isolate; each cell one vial; colour depth = fitted growth rate r (h⁻¹)",transform=axB.transAxes,fontsize=6.2,color=MUTED)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm); sm.set_array([])
cb=fig.colorbar(sm,ax=axB,fraction=0.018,pad=0.01,ticks=[0.02,0.05,0.1,0.2,0.5,1.0]); cb.ax.set_yticklabels(["<0.02","0.05","0.1","0.2","0.5","1.0"],fontsize=6); cb.ax.tick_params(length=2)
cb.set_label("r (h⁻¹)",fontsize=6.5)
# ---- (c) the 40 C result -----------------------------------------------------
axC=fig.add_subplot(gs[2]); s40=d[d["T"]==40].copy()
s40["clade"]=np.where(s40.group.str.startswith("Clade"),"C. auris","relatives")
rng=np.random.default_rng(3); ypos2={"C. auris":1,"relatives":0}; txt={}
for g,x in s40.groupby("clade"):
    p=(x.S=="growth").mean(); isol=x.otu_name.unique(); b=[]
    for _ in range(4000):
        pick=rng.choice(isol,len(isol),replace=True)
        b.append(pd.concat([x[x.otu_name==i] for i in pick]).S.eq("growth").mean())
    lo,hi=np.percentile(b,[2.5,97.5]); y=ypos2[g]; c=AURIS if g=="C. auris" else REL
    axC.plot([lo,hi],[y,y],color=c,lw=2.0,solid_capstyle="round",zorder=2)
    axC.scatter(p,y,s=46,color=c,edgecolor="white",linewidth=0.8,zorder=3)
    n=(x.S=="growth").sum(); txt[g]=(n,len(x),p,lo,hi)
    axC.text(hi+0.03,y,f"{n}/{len(x)} vials   {100*p:.1f}%  [{100*lo:.0f}, {100*hi:.0f}]",va="center",fontsize=6.6,color=INK)
iso_ok=s40.groupby(["clade","otu_name"]).apply(lambda x:(x.S=="growth").mean()>=0.5,include_groups=False).groupby("clade").agg(["sum","size"])
axC.set_yticks([0,1]); axC.set_yticklabels(["relatives (8 isolates)","C. auris (12 isolates)"],fontsize=7)
axC.set_xlim(-0.03,1.0); axC.set_xticks([0,0.25,0.5,0.75,1.0]); axC.set_xticklabels(["0","25","50","75","100%"])
axC.set_xlabel("vials with detectable growth at 40 °C (isolate-clustered 95% CI)",fontsize=7.5); axC.set_ylim(-0.7,1.7)
for sp in ("top","right","left"): axC.spines[sp].set_visible(False)
axC.tick_params(axis="y",length=0); axC.tick_params(axis="x",length=2)
axC.text(-0.175,1.10,"c",transform=axC.transAxes,fontweight="bold",fontsize=9)
axC.text(0.0,-0.62,f"Isolate level (at least half of the vials growing): C. auris {int(iso_ok.loc['C. auris','sum'])}/12, relatives {int(iso_ok.loc['relatives','sum'])}/8.",
         transform=axC.transAxes,fontsize=6.2,color=MUTED,linespacing=1.35)
h=[Patch(facecolor=COL[s],edgecolor="white",label=LAB[s]) for s in ORDER]
fig.legend(handles=h,frameon=False,fontsize=6.8,ncol=3,loc="upper center",bbox_to_anchor=(0.5,1.0),columnspacing=1.6,handlelength=1.5)
fig.subplots_adjust(top=0.905)
save(fig,"FIG4_growth_lik"); print(d.S.value_counts().to_dict()); print(txt); print(iso_ok)
