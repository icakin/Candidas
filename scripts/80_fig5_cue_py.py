"""FIGURE 5 and Supplementary Figs 1-4 (Python version of 17_fig2.R): apparent CUE curves,
CUE optima, the 37->40 C bill, relative R/G curves, T_opt vs fever cost, leave-one-out rho,
pairwise contrasts. Values with intervals are read from results/tables/fig_values.csv and
fig_contrasts.csv (written by 17_fig2.R) so the printed numbers are unchanged; curves and
isolate points are recomputed from the posterior draws dumped by 68_isolate_params_dump.R."""
from tpc_common import *
from scipy.stats import spearmanr
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
def logticks(ax,vals,labels,axis='x'):
    a=ax.xaxis if axis=='x' else ax.yaxis
    a.set_major_locator(FixedLocator(vals)); a.set_major_formatter(FixedFormatter(labels)); a.set_minor_locator(FixedLocator([])); a.set_minor_formatter(NullFormatter())
G,I,quota,vals=load(); P=points()
Tmax={g:P[P.group==g]["T"].max() for g in GROUPS}
# supported range: last temperature at which at least half of the group's wells (of 15, or 10 for Duo) carry a fitted growth rate
_nw={g:P[P.group==g].groupby("T").size() for g in GROUPS}
Tsup={g:max(t for t,k in _nw[g].items() if k>=0.5*_nw[g].max()) for g in GROUPS}
def split_curve(ax,g,fn,lw):
    grid=TG[TG<=Tmax[g]]; med,lo,hi=curve(G[g],fn,grid); m=grid<=Tsup[g]
    ax.fill_between(grid[m],lo[m],hi[m],color=OI[g],alpha=0.10,lw=0); ax.plot(grid[m],med[m],color=OI[g],lw=lw)
    if (~m).any():
        mx=grid>=Tsup[g]; ax.plot(grid[mx],med[mx],color=OI[g],lw=lw*0.7,ls=(0,(2,1.5)))
    return grid,med
wells40=pd.read_csv(f"{C}/results/tables/wells_by_isolate.csv").rename(columns={"otu_name":"Isolate"}).set_index("Isolate")["40"]
def cue_fn(t,x): Gc=np.exp(x.lnB0.values+lnG(t,x.E.values,x.Eh.values,x.Th.values)); R=np.exp(x.alpha.values+lnR(t,x.ER.values)); return Gc/(Gc+R)
# ================= FIGURE 5 =================
fig=plt.figure(figsize=(7.2,8.3)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=0.38,wspace=0.3)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); gsc=gs[1,:].subgridspec(1,2,wspace=0.08); axC1=fig.add_subplot(gsc[0]); axC2=fig.add_subplot(gsc[1],sharey=axC1)
# a: CUE curves
axA.axvspan(T_BODY,T_FEVER,color=RED,alpha=0.06,lw=0); axA.axvline(T_BODY,color=RED,lw=0.5,ls=(0,(2,2))); axA.axvline(T_FEVER,color=RED,lw=0.5,ls=(0,(2,2)))
for g in GROUPS:
    grid,med=split_curve(axA,g,cue_fn,1.6)
    k=int(np.argmax(med)); axA.scatter(grid[k],med[k],s=30,color=OI[g],edgecolor="white",linewidth=0.8,zorder=4)
axA.set_ylim(0,1); axA.set_xlabel("temperature (°C)"); axA.set_ylabel("apparent CUE, G/(G+R)"); axA.set_xticks([24,28,32,36,40,44])
axA.legend(handles=[Line2D([],[],marker="o",ls="",color=MUTED,markeredgecolor="white",ms=5,label="peak of each curve"),Line2D([],[],color=MUTED,lw=1.1,ls=(0,(2,1.5)),label="beyond supported growth range")],loc="lower left",fontsize=6.2); panel_label(axA,"a",x=-0.2)
# b: CUE optimum, growth optimum, isolates, P(<37)
axB.axvspan(T_BODY,T_FEVER,color=RED,alpha=0.06,lw=0); axB.axvline(T_BODY,color=RED,lw=0.5,ls=(0,(2,2))); axB.axvline(T_FEVER,color=RED,lw=0.5,ls=(0,(2,2)))
for i,g in enumerate(GROUPS):
    v=vals.loc[g]; y=len(GROUPS)-1-i
    axB.plot([v.Tcue_lo,v.Tcue_hi],[y,y],color=OI[g],lw=1.6); axB.scatter(v.Tcue,y,s=34,color=OI[g],zorder=4)
    axB.plot([v.Topt_lo,v.Topt_hi],[y,y],color=AMB,lw=0.9,zorder=2); axB.scatter(v.Topt,y,s=26,marker="s",facecolor="white",edgecolor=MUTED,linewidth=0.8,zorder=3)
    iso_t=[np.nanmedian(topt(I[k].E.values-I[k].ER.values,I[k].Eh.values,I[k].Th.values)) for k in I if k.startswith(g)]
    axB.scatter(iso_t,[y+0.32]*len(iso_t),s=6,color=OI[g],alpha=0.7,lw=0)
    axB.text(v.Tcue,y-0.42,f"{v.Tcue:.1f}",ha="center",fontsize=6.2,color=OI[g])
    axB.text(45.6,y,f"P = {v.P_below_body:.3f}" if v.P_below_body<0.9995 else "P = 1.000",fontsize=6.0,color=MUTED,va="center")
axB.set_yticks(range(len(GROUPS))); axB.set_yticklabels([SHORT[g] for g in GROUPS][::-1],style="italic",fontsize=6.6); axB.tick_params(axis="y",length=0); axB.spines["left"].set_visible(False)
axB.set_xlim(27,48.5); axB.set_xticks([28,32,36,40,44]); axB.set_xlabel("temperature (°C)"); axB.set_ylim(-0.8,len(GROUPS)-0.3)
axB.legend(handles=[Line2D([],[],marker="o",ls="-",color=MUTED,ms=5,label="CUE optimum (95% CrI)"),Line2D([],[],marker="s",ls="",markerfacecolor="white",markeredgecolor=MUTED,color=MUTED,ms=5,label="growth optimum"),Line2D([],[],marker="o",ls="",color=MUTED,ms=2.5,label="isolates")],loc="lower right",fontsize=5.8,ncol=3,bbox_to_anchor=(1.0,1.0),columnspacing=0.8,handletextpad=0.3); axB.text(45.6,len(GROUPS)-0.55,"P = P(T_opt(CUE) < 37 °C)",fontsize=5.8,color=MUTED)
panel_label(axB,"b",x=-0.3,y=1.07)
# c: the bill
for ax,(vk,lok,hik,lab,fmt,isofn) in zip((axC1,axC2),(("growth_kept","gk_lo","gk_hi","growth-rate retention, 37 to 40 °C","{:.2f}",lambda x: np.exp(lnG(T_FEVER,x.E.values,x.Eh.values,x.Th.values)-lnG(T_BODY,x.E.values,x.Eh.values,x.Th.values))),
        ("fever_cost","fc_lo","fc_hi","R/G change, 37 to 40 °C","{:.2f}×",lambda x: np.exp((lnR(T_FEVER,x.ER.values)-lnG(T_FEVER,x.E.values,x.Eh.values,x.Th.values))-(lnR(T_BODY,x.ER.values)-lnG(T_BODY,x.E.values,x.Eh.values,x.Th.values)))))):
    ax.axvline(1.0,color=AMB,lw=1.0)
    for i,g in enumerate(GROUPS):
        y=len(GROUPS)-1-i; v=vals.loc[g]
        if g=="Duo":
            ax.text(0.72 if vk=="growth_kept" else 1.55,y,"no detectable growth at 40 °C" if vk=="growth_kept" else "not data-supported at 40 °C",fontsize=6.0,color=MUTED,style="italic",ha="center",va="center"); continue
        ax.plot([v[lok],v[hik]],[y,y],color=OI[g],lw=1.8); ax.scatter(v[vk],y,s=40,color=OI[g],zorder=4)
        ax.text(v[vk],y+0.3,fmt.format(v[vk]),ha="center",fontsize=6.4,color=OI[g],fontweight="bold")
        iso=[(k,np.nanmedian(isofn(I[k]))) for k in I if k.startswith(g)]
        for k,val in iso:
            grew=wells40.get(k,1)>0
            ax.scatter(val,y-0.3,s=9,alpha=0.75 if grew else 1,facecolor=OI[g] if grew else "white",edgecolor=OI[g],linewidth=0.6,zorder=3)
        pos=[val for k,val in iso if wells40.get(k,1)>0]
        if 0<len(pos)<len(iso): ax.scatter(np.median(pos),y-0.15,marker="D",s=16,color=OI[g],zorder=4)
    ax.set_xscale("log"); ax.spines["left"].set_visible(False); ax.tick_params(axis="y",length=0)
    ax.set_title(lab,fontsize=7.2,pad=6)
axC1.set_yticks(range(len(GROUPS))); axC1.set_yticklabels([SHORT[g] for g in GROUPS][::-1],style="italic",fontsize=6.6); axC1.set_ylim(-0.8,len(GROUPS)-0.2)
axC1.set_xlim(0.42,1.08); logticks(axC1,[0.5,0.7,1.0],["0.5×","0.7×","1×"]); axC2.set_xlim(0.95,2.6); logticks(axC2,[1,1.5,2,2.5],["1×","1.5×","2×","2.5×"]); plt.setp(axC2.get_yticklabels(),visible=False)
axC1.legend(handles=[Line2D([],[],marker="o",ls="-",color=MUTED,ms=5.5,label="hierarchical estimate (95% CrI)"),Line2D([],[],marker="o",ls="",color=MUTED,ms=2.5,label="isolates"),Line2D([],[],marker="o",ls="",markerfacecolor="white",markeredgecolor=MUTED,ms=3,label="isolate without growth at 40 °C"),Line2D([],[],marker="D",ls="",color=MUTED,ms=3.5,label="median among growth-positive isolates")],loc="lower left",bbox_to_anchor=(0.0,1.06),ncol=4,fontsize=5.8,columnspacing=0.8,handletextpad=0.3)
panel_label(axC1,"c",x=-0.34,y=1.12)
save(fig,"FIG2_consequences_py")
# ================= SUPP 1: relative R/G =================
fig,ax=plt.subplots(figsize=(6.0,4.2))
ax.axhline(1,color=AMB,lw=0.8); ax.axvline(T_BODY,color=RED,lw=0.6,ls=(0,(2,2))); ax.axvline(T_FEVER,color=RED,lw=0.6,ls=(0,(2,2)))
ax.text(T_BODY-0.3,7.5,"body 37 °C",ha="right",fontsize=6.2,color=RED); ax.text(T_FEVER+0.3,7.5,"fever 40 °C",ha="left",fontsize=6.2,color=RED)
for g in GROUPS:
    split_curve(ax,g,lambda t,x: tax(t,x.E.values,x.ER.values,x.Eh.values,x.Th.values),1.5)
    ax.scatter(vals.loc[g].Tcue,1.0,s=30,color=OI[g],edgecolor="white",linewidth=0.8,zorder=4)
ax.set_yscale("log"); ax.set_ylim(0.95,8.5); logticks(ax,[1,1.5,2,3,5,8],["1×","1.5×","2×","3×","5×","8×"],"y"); ax.set_xticks([24,28,32,36,40,44])
ax.set_xlabel("temperature (°C)"); ax.set_ylabel("respiration per unit growth, relative to each group's minimum")
h=[Line2D([],[],color=OI[g],lw=2,label=SHORT[g]) for g in GROUPS]+[Line2D([],[],marker="o",ls="",color=MUTED,markeredgecolor="white",ms=5,label="minimum = T_opt(CUE)"),Line2D([],[],color=MUTED,lw=1.1,ls=(0,(2,1.5)),label="beyond supported growth range")]; ax.legend(handles=h,loc="upper left",fontsize=6.4)
save(fig,"FIG_SUPP_relative_RG_py")
# ================= SUPP 2: Topt vs fever cost =================
fig,ax=plt.subplots(figsize=(4.8,4.2))
for g in GROUPS:
    if g=="Duo": continue
    v=vals.loc[g]; ax.plot([v.Topt_lo,v.Topt_hi],[v.fever_cost]*2,color=OI[g],lw=0.9,alpha=0.6); ax.plot([v.Topt]*2,[v.fc_lo,v.fc_hi],color=OI[g],lw=0.9,alpha=0.6)
    ax.scatter(v.Topt,v.fever_cost,s=44,color=OI[g],edgecolor="white",linewidth=0.8,zorder=4)
ax.axvline(T_BODY,color=RED,lw=0.6,ls=(0,(2,2))); ax.text(T_BODY,ax.get_ylim()[1],"body 37 °C",fontsize=6,color=RED,ha="right",va="top")
ax.set_yscale("log"); logticks(ax,[1.2,1.5,2,3,5,8],["1.2×","1.5×","2×","3×","5×","8×"],"y"); ax.set_xlim(29.5,38.5)
ax.set_xlabel("fitted growth optimum T_opt (°C)"); ax.set_ylabel("respiration per unit growth, 40 vs 37 °C (×)"); group_legend(ax,loc="upper left",bbox_to_anchor=(1.01,1.0),fontsize=6.2,labelspacing=0.3,handletextpad=0.3)
save(fig,"FIG_SUPP_topt_vs_fever_cost_py")
# ================= SUPP 3: leave-one-out rho =================
n=min(len(G[g]) for g in GROUPS)
Tm=np.column_stack([topt(G[g].E.values[:n],G[g].Eh.values[:n],G[g].Th.values[:n]) for g in GROUPS])
Cm=np.column_stack([np.exp((lnR(T_FEVER,G[g].ER.values[:n])-lnG(T_FEVER,G[g].E.values[:n],G[g].Eh.values[:n],G[g].Th.values[:n]))-(lnR(T_BODY,G[g].ER.values[:n])-lnG(T_BODY,G[g].E.values[:n],G[g].Eh.values[:n],G[g].Th.values[:n]))) for g in GROUPS])
def rho(cols,collapse=False):
    out=[]
    for j in range(n):
        a=Tm[j,cols]; b=Cm[j,cols]
        if collapse: a=np.r_[np.mean(Tm[j,:4]),Tm[j,4:]]; b=np.r_[np.mean(Cm[j,:4]),Cm[j,4:]]
        if np.all(np.isfinite(a))&np.all(np.isfinite(b)): out.append(spearmanr(a,b).correlation)
    out=np.array(out); return np.median(out),np.percentile(out,2.5),np.percentile(out,97.5),np.mean(out<0)
rows=[("six groups with growth at 40 °C",rho([k for k in range(7) if GROUPS[k]!="Duo"])),("all seven groups (Duo extrapolated)",rho(list(range(7))))]+[(f"without {SHORT[g]}",rho([k for k in range(7) if k!=i])) for i,g in enumerate(GROUPS)]+[("C. auris clades collapsed (n = 4)",rho(None,True))]
fig,ax=plt.subplots(figsize=(5.6,3.6)); ax.axvline(0,color=AMB,lw=0.9)
for i,(lab,(m,lo,hi,pn)) in enumerate(rows):
    y=len(rows)-1-i; ax.plot([lo,hi],[y,y],color=INK,lw=1.2); ax.scatter(m,y,s=26,color=INK,zorder=3); ax.text(m,y+0.28,f"P(ρ<0) = {pn:.2f}",fontsize=6,color=MUTED,ha="center")
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1],fontsize=6.6); ax.tick_params(axis="y",length=0); ax.spines["left"].set_visible(False)
ax.set_xlim(-1.05,1.05); ax.set_xlabel("posterior Spearman ρ (growth T_opt vs 37 to 40 °C change in R/G)"); ax.set_ylim(-0.6,len(rows)-0.2)
save(fig,"FIG_SUPP_loo_correlation_py")
print("rho six groups: median %.2f [%.2f, %.2f], P(rho<0)=%.3f"%rows[0][1])
# ================= SUPP 4: pairwise contrasts =================
ct=pd.read_csv(f"{C}/results/tables/fig_contrasts.csv"); ab={"Clade1":"I","Clade2":"II","Clade3":"III","Clade4":"IV","para":"para","Hae":"hae","Duo":"duo"}
ct=ct[(ct.a!="Duo")&(ct.b!="Duo")]; ct["lab"]=ct.a.map(ab)+" vs "+ct.b.map(ab); ct=ct.sort_values("ratio",ascending=False).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(6.0,5.2)); ax.axvline(1,color=INK,lw=0.8)
for i,r in ct.iterrows():
    y=len(ct)-1-i; c=INK if r.credible else AMB; ax.plot([r.lo,r.hi],[y,y],color=c,lw=1.4); ax.scatter(r.ratio,y,s=22,color=c,zorder=3)
ax.set_yticks(range(len(ct))); ax.set_yticklabels(ct.lab[::-1],fontsize=6.4); ax.tick_params(axis="y",length=0); ax.spines["left"].set_visible(False)
ax.set_xscale("log"); logticks(ax,[0.1,0.2,0.3,0.5,0.7,1,1.5,2],["0.1","0.2","0.3","0.5","0.7","1.0","1.5","2.0"])
ax.set_xlabel("ratio of relative R/G at 40 °C  (< 1 = the first group pays less)")
ax.text(0.0,1.01,"dark, 95% CrI excludes 1; light, not separable.  I–IV = C. auris clades; hae = C. haemulonii. C. duobushaemulonii has no growth at 40 °C and is excluded",transform=ax.transAxes,fontsize=5.8,color=MUTED)
save(fig,"FIG_SUPP_pairwise_contrasts_py")
