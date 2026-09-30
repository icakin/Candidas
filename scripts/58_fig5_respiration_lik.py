"""FIGURE 5 (likelihood model). Oxygen consumption after growth is lost. Well states come from
the per-well likelihood model (growth if bootstrap p < 0.01; no detectable respiration if
drawdown < 1 mg/L); consumption on the fixed, prespecified 90-390 min window is read from
results/tables/phase2b_consumption.csv (43_phase2b_consumption.py) and is state-independent.
Share = oxygen consumed in growing wells / oxygen consumed in all wells of that isolate at that
temperature. Uncertainty in (a) is an isolate-level bootstrap (isolates resampled within clade)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
P_GROW=0.01; DRAWDOWN_MIN=1.0
l=pd.read_csv(f"{C}/reprocess/lik_results_labelled.csv")
l["S"]=np.where(l.drawdown<DRAWDOWN_MIN,"no_resp",np.where(l.p_boot<P_GROW,"growth","resp_only"))
cons=pd.read_csv(f"{C}/results/tables/phase2b_consumption.csv")[["WELLID","cons","censored"]]
d=l.merge(cons,on="WELLID"); d=d[np.isfinite(d.cons)].copy()
d["clade"]=np.where(d.group.str.startswith("Clade"),"C. auris","relatives")
# per isolate x T share
rows=[]
for (g,iso,T),s in d.groupby(["group","otu_name","T"]):
    tot=s.cons.sum()
    if tot>0: rows.append(dict(group=g,isolate=iso,T=T,tot=tot,frac=s[s.S=="growth"].cons.sum()/tot,
                              clade="C. auris" if g.startswith("Clade") else "relatives"))
p=pd.DataFrame(rows); p.to_csv(f"{C}/results/tables/lik_productive_share.csv",index=False)
def pooled(df,c,rng=None,B=0):
    x=df[df.clade==c]; isol=sorted(x.isolate.unique()); out=[]
    for T,y in x.groupby("T"):
        est=np.average(y.frac,weights=y.tot)
        if B:
            b=[]
            for _ in range(B):
                pick=rng.choice(isol,len(isol),replace=True)
                z=pd.concat([y[y.isolate==i] for i in pick]); b.append(np.average(z.frac,weights=z.tot))
            lo,hi=np.percentile(b,[2.5,97.5])
        else: lo=hi=est
        out.append(dict(T=T,est=est,lo=lo,hi=hi))
    return pd.DataFrame(out)
rng=np.random.default_rng(5)
ORD=ISOLATES_TREE
fig,axes=plt.subplots(1,3,figsize=(8.0,3.3),gridspec_kw=dict(wspace=0.46,width_ratios=[1.05,0.95,0.95]))
ax=axes[0]; summ={}
for c,colr in (("C. auris",AURIS),("relatives",REL)):
    g=pooled(p,c,rng,2000); summ[c]=g
    ax.fill_between(g["T"],g.lo,g.hi,color=colr,alpha=0.20,lw=0)
    ax.plot(g["T"],g.est,color=colr,lw=1.6,label=c)
ax.axvspan(36.8,37.2,color=INK,alpha=0.10,lw=0,zorder=0); ax.axvspan(39.8,40.2,color=INK,alpha=0.06,lw=0,zorder=0)
ax.set_xlabel("temperature (°C)",fontsize=7.4); ax.set_ylabel("share of O₂ use in\nwells with detectable growth",fontsize=7.4)
ax.set_ylim(-0.03,1.05); ax.set_xlim(21.5,45); ax.set_xticks([22,26,30,34,38,42])
ax.legend(frameon=False,fontsize=6.4,loc="lower left")
ax.text(37,0.04,"body",rotation=90,ha="center",va="bottom",fontsize=5.8,color=MUTED)
ax.text(40,0.04,"fever",rotation=90,ha="center",va="bottom",fontsize=5.8,color=MUTED)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.text(-0.36,1.10,"a",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("band: isolate-level bootstrap 95% CI",fontsize=6.2,color=MUTED,pad=10)
ax=axes[1]; q=p[p["T"]==40].set_index("isolate")
for yi,iso in enumerate(ORD):
    if iso in q.index: ax.scatter([q.loc[iso].frac],[yi],s=20,color=col(iso),edgecolor=INK,linewidth=0.5,zorder=4)
ax.set_yticks(range(len(ORD))); ax.set_yticklabels([iso_label(i)+(" *" if i in ("Hae_1724","Clade2_2073") else "") for i in ORD],fontsize=5.4,style="italic")
for lab,i in zip(ax.get_yticklabels(),ORD): lab.set_color(col(i))
ax.set_ylim(len(ORD)-0.4,-0.9); ax.set_xlim(-0.04,1.08); ax.set_xticks([0,0.5,1]); ax.set_xlabel("share at 40 °C",fontsize=7.4)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y",length=0); ax.text(-0.62,1.10,"b",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("per isolate",fontsize=6.2,color=MUTED,pad=10)
ax=axes[2]; r40=d[(d["T"]==40)&(d.S=="resp_only")]
lab=dict(SHORT); lab["Duo"]="C. duobus."
GR=[g for g in ["Clade1","Clade2","Clade3","Clade4","para","Duo","Hae"] if (r40.group==g).any()]
for yi,g in enumerate(GR):
    v=r40[r40.group==g].cons.values; c=OI[g]
    ax.scatter(v,np.full(len(v),yi)+np.random.RandomState(0).uniform(-.12,.12,len(v)),s=12,color=c,alpha=0.65,zorder=3)
    ax.plot([np.median(v)]*2,[yi-0.28,yi+0.28],color=c,lw=1.8,zorder=4)
a40=d[(d["T"]==40)&(d.S=="growth")].cons.median()
ax.axvline(a40,color=MUTED,lw=1.0,ls=(0,(3,2)),zorder=1)
ax.text(a40-0.2,-1.15,"median of growing\nwells at 40 °C",fontsize=5.6,color=MUTED,ha="right",linespacing=1.2)
ax.axvline(1.0,color=RED,lw=0.8,ls=":",zorder=1); ax.text(1.0,-1.15,"1 mg L⁻¹\ngate",fontsize=5.6,color=RED,va="center",ha="center",linespacing=1.2)
ax.set_yticks(range(len(GR))); ax.set_yticklabels([lab[g] for g in GR],fontsize=6.0); ax.set_ylim(len(GR)-0.5,-1.5)
ax.set_xlabel("O₂ used in the fixed window (mg L⁻¹)",fontsize=7.4)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y",length=0); ax.text(-0.50,1.10,"c",transform=ax.transAxes,fontweight="bold",fontsize=9)
ax.set_title("vials respiring without detectable growth",fontsize=6.2,color=MUTED,pad=10)
save(fig,"FIG5_respiration_lik")
for c in summ: r=summ[c][summ[c]["T"]==40].iloc[0]; print(f"  pooled share at 40 C, {c}: {r.est:.2f} [{r.lo:.2f}, {r.hi:.2f}]")
q=p[p["T"]==40]; print("  per isolate at 40 C:"); print(q[["isolate","frac"]].round(2).to_string(index=False))
ex=q[(q.clade=="relatives")&(q.isolate!="Hae_1724")]; print(f"  relatives excl. Hae_1724 pooled: {np.average(ex.frac,weights=ex.tot):.2f}")
print("  resp-only vials at 40 C, consumption by group (median, n):"); print(r40.groupby("group").cons.agg(['median','size']).round(2).to_string())
print(f"  growing vials at 40 C median cons {a40:.2f}; censored share growing {d[(d['T']==40)&(d.S=='growth')].censored.mean():.2f}, resp-only {r40.censored.mean():.2f}")
# ratio resp-only / growing at 38 C within isolate (relatives)
g38=d[(d['T']==38)&(d.S=='growth')].groupby('otu_name').cons.median()
rr=r40.groupby('otu_name').cons.median(); ratio=(rr/g38.reindex(rr.index)).dropna(); print("  resp-only@40 / growing@38 same isolate:"); print(ratio.round(2).to_string())
