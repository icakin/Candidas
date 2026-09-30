"""FIGURE 5 (corrected) - states from the FROZEN v4 blind classifier.
Denominator 1,200 wells. Four observed categories, ambiguous drawn
explicitly as an UNCERTAINTY category (not a biological state).
States are read from reprocess/unblinded_results.csv and are NEVER recomputed here.
"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
BLUE,ORANGE,INK,MUTED,FAINT="#2a78d6","#eb6834","#0b0b0b","#52514e","#e8e8e6"
AMB,NORESP="#b9b7b2","#f2f1ef"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
d=pd.read_csv(f"{C}/reprocess/unblinded_results.csv"); ISO=set(pd.read_csv(f"{C}/results/tables/otu_names.csv").OTU); d=d[d.OTU.isin(ISO)].copy()  # study isolates only
NM={1:"growth",2:"resp_only",3:"ambig",4:"no_resp"}; d["S"]=d.state.map(NM)
ORDER=["growth","resp_only","ambig","no_resp"]
COL={"growth":BLUE,"resp_only":ORANGE,"ambig":AMB,"no_resp":NORESP}
LAB={"growth":"confirmed growth + respiration","resp_only":"respiration without confirmed growth",
     "ambig":"ambiguous (uncertainty)","no_resp":"no detectable respiration"}
TAXA=["Clade1","Clade2","Clade3","Clade4","para","Hae","Duo"]
TNAME={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
       "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
Ts=sorted(d["T"].unique())
fig=plt.figure(figsize=(7.2,9.0))
gs=fig.add_gridspec(3,1,height_ratios=[1.05,1.62,0.70],hspace=0.46)
# ---- (a) state proportions across temperature, per taxon ---------------------
gsa=gs[0].subgridspec(1,7,wspace=0.12)
for i,t in enumerate(TAXA):
    ax=fig.add_subplot(gsa[0,i]); s=d[d.group==t]
    frac=(s.groupby(["T","S"]).size().unstack(fill_value=0).reindex(columns=ORDER,fill_value=0))
    frac=frac.div(frac.sum(1),axis=0)
    bot=np.zeros(len(frac))
    for st in ORDER:
        ax.bar(range(len(frac)),frac[st].values,bottom=bot,color=COL[st],width=0.92,
               lw=0.3,edgecolor="white",zorder=2)
        bot+=frac[st].values
    ax.axvspan(Ts.index(38)-0.5,Ts.index(40)+0.5,color=INK,alpha=0.07,lw=0,zorder=1)
    ax.set_xticks(range(len(Ts))); ax.set_xticklabels([str(x) if x in (22,30,38,44) else "" for x in Ts],fontsize=5.8)
    ax.set_ylim(0,1); ax.set_xlim(-0.6,len(Ts)-0.4)
    ax.set_title(TNAME[t],fontsize=6.3,style="italic",pad=3)
    if i==0: ax.set_ylabel("fraction of wells",fontsize=7)
    else: ax.set_yticklabels([])
    ax.tick_params(length=2)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
axA=fig.axes[0]; axA.text(-0.42,1.30,"a",transform=axA.transAxes,fontweight="bold",fontsize=9)
axA.text(0.0,1.16,"shaded band: 37-40 °C, host to fever",transform=axA.transAxes,fontsize=6.2,color=MUTED)
# ---- (b) per-ISOLATE state, every well ---------------------------------------
axB=fig.add_subplot(gs[1])
isos=[]
for t in TAXA:
    for o in sorted(d[d.group==t].otu_name.unique()): isos.append((t,o))
for yi,(t,o) in enumerate(isos):
    s=d[d.otu_name==o]
    for xi,T in enumerate(Ts):
        w=s[s["T"]==T]
        if len(w)==0: continue
        cnt=w.S.value_counts(); x0=xi-0.45; tot=len(w)
        for st in ORDER:
            n=cnt.get(st,0)
            if n==0: continue
            wdt=0.9*n/tot
            axB.barh(yi,wdt,left=x0,height=0.74,color=COL[st],lw=0.25,edgecolor="white",zorder=2)
            x0+=wdt
axB.axvspan(Ts.index(38)-0.5,Ts.index(40)+0.5,color=INK,alpha=0.07,lw=0,zorder=1)
axB.set_yticks(range(len(isos)))
axB.set_yticklabels([f"{o}" + ("  ←" if o=="Hae_1724" else "") for _,o in isos],fontsize=6.1)
axB.invert_yaxis(); axB.set_xticks(range(len(Ts))); axB.set_xticklabels([str(x) for x in Ts],fontsize=6.4)
axB.set_xlabel("temperature (°C)",fontsize=7.5); axB.set_xlim(-0.6,len(Ts)-0.4)
prev=None
for yi,(t,o) in enumerate(isos):
    if prev is not None and t!=prev: axB.axhline(yi-0.5,color=MUTED,lw=0.5,alpha=0.5)
    prev=t
for sp in ("top","right","left"): axB.spines[sp].set_visible(False)
axB.tick_params(axis="y",length=0); axB.tick_params(axis="x",length=2)
axB.text(-0.175,1.02,"b",transform=axB.transAxes,fontweight="bold",fontsize=9)
axB.text(0.0,1.02,"each bar is one isolate at one temperature, split by the 5 wells",transform=axB.transAxes,fontsize=6.2,color=MUTED)
# Hae_1724 is flagged by the arrow on its y-label and named in the caption;
# no in-panel annotation, to avoid obscuring wells.
# ---- (c) the 40 C result -----------------------------------------------------
axC=fig.add_subplot(gs[2]); s40=d[d["T"]==40].copy()
s40["clade"]=np.where(s40.group.str.startswith("Clade"),"C. auris","relatives")
rng=np.random.default_rng(3); ypos={"C. auris":1,"relatives":0}
for g,x in s40.groupby("clade"):
    p=(x.S=="growth").mean(); isol=x.otu_name.unique(); b=[]
    for _ in range(4000):
        pick=rng.choice(isol,len(isol),replace=True)
        b.append(pd.concat([x[x.otu_name==i] for i in pick]).S.eq("growth").mean())
    lo,hi=np.percentile(b,[2.5,97.5]); y=ypos[g]
    c=BLUE if g=="C. auris" else ORANGE
    axC.plot([lo,hi],[y,y],color=c,lw=2.0,solid_capstyle="round",zorder=2)
    axC.scatter(p,y,s=46,color=c,edgecolor="white",linewidth=0.8,zorder=3)
    n=(x.S=="growth").sum()
    axC.text(hi+0.03,y,f"{n}/{len(x)} wells   {100*p:.1f}%  [{100*lo:.0f}, {100*hi:.0f}]",
             va="center",fontsize=6.6,color=INK)
axC.set_yticks([0,1]); axC.set_yticklabels(["relatives (8 isolates)","C. auris (12 isolates)"],fontsize=7)
axC.set_xlim(-0.03,1.0); axC.set_xticks([0,0.25,0.5,0.75,1.0]); axC.set_xticklabels(["0","25","50","75","100%"])
axC.set_xlabel("confirmed-growth wells at 40 °C (isolate-clustered 95% CI)",fontsize=7.5)
axC.set_ylim(-0.7,1.7)
for sp in ("top","right","left"): axC.spines[sp].set_visible(False)
axC.tick_params(axis="y",length=0); axC.tick_params(axis="x",length=2)
axC.text(-0.175,1.10,"c",transform=axC.transAxes,fontweight="bold",fontsize=9)
axC.text(0.0,-0.62,"Isolate level: C. auris 11/12. Relatives bounded 1/8 to 5/8 by four unresolved isolates\n"
                   "(risk difference +0.29 to +0.79). Direction holds under both allocations; magnitude does not.",
         transform=axC.transAxes,fontsize=6.2,color=MUTED,linespacing=1.35)
h=[Patch(facecolor=COL[s],edgecolor="white",label=LAB[s]) for s in ORDER]
fig.legend(handles=h,frameon=False,fontsize=6.8,ncol=2,loc="upper center",bbox_to_anchor=(0.5,1.005),columnspacing=1.6,handlelength=1.5)
fig.subplots_adjust(top=0.905)
fig.savefig(f"{FIG}/FIG5_corrected_states.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG}/FIG5_corrected_states.pdf",bbox_inches="tight")
print("saved FIG5_corrected_states")
print(f"total wells drawn: {len(d)} | ambiguous drawn: {(d.S=='ambig').sum()}")
