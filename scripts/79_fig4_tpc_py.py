"""FIGURE 4 (Python version of 16_fig1.R): growth and respiration thermal performance
curves per group, per-isolate median curves, wells as points, activation energies."""
from tpc_common import *
G,I,quota,vals=load(); P=points()
fig=plt.figure(figsize=(7.2,7.6)); gs=fig.add_gridspec(2,2,height_ratios=[1,0.95],hspace=0.35,wspace=0.28)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,:])
for ax in (axA,axB):
    ax.axvspan(T_BODY,T_FEVER,color=RED,alpha=0.06,lw=0)
# ---- a: growth ----
for g in GROUPS:
    d=G[g]; q=quota[g]; Tmax=P[P.group==g]["T"].max(); grid=TG[TG<=Tmax]
    med,lo,hi=curve(d,lambda t,x: np.exp(x.lnB0.values+lnG(t,x.E.values,x.Eh.values,x.Th.values))/q,grid)
    axA.fill_between(grid,lo,hi,color=OI[g],alpha=0.12,lw=0); axA.plot(grid,med,color=OI[g],lw=1.6,zorder=3)
    for iso,di in I.items():
        if not iso.startswith(g): continue
        m,_,_=curve(di,lambda t,x: np.exp(x.lnB0.values+lnG(t,x.E.values,x.Eh.values,x.Th.values))/q,grid)
        axA.plot(grid,m,color=OI[g],lw=0.5,alpha=0.55,zorder=2)
    p=P[P.group==g]; axA.scatter(p["T"]+np.random.RandomState(0).uniform(-.25,.25,len(p)),p.r_h,s=4,color=OI[g],alpha=0.35,lw=0,zorder=1)
    nd=[t for t in np.arange(22,45,2) if t>Tmax]
    if nd: axA.scatter(nd,[0.05]*len(nd),marker="x",s=18,color=OI[g],lw=1.0,zorder=4)
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
axA.set_yscale("log"); axA.set_ylim(0.045,1.6); axA.yaxis.set_major_locator(FixedLocator([0.1,0.3,1.0])); axA.yaxis.set_major_formatter(FixedFormatter(["0.1","0.3","1.0"])); axA.yaxis.set_minor_formatter(NullFormatter())
axA.set_xlabel("temperature (°C)"); axA.set_ylabel("specific growth rate r (h⁻¹)"); axA.set_xticks([24,28,32,36,40,44])
axA.text(0.02,0.03,"× assayed, no detectable growth",transform=axA.transAxes,fontsize=6.2,color=MUTED)
panel_label(axA,"a",x=-0.2)
# ---- b: respiration ----
_nw={g:P[P.group==g].groupby("T").size() for g in GROUPS}
Tsup={g:max(t for t,k in _nw[g].items() if k>=0.5*_nw[g].max()) for g in GROUPS}
for g in GROUPS:
    d=G[g]; Tmax=P[P.group==g]["T"].max(); grid=TG[TG<=Tmax]
    med,lo,hi=curve(d,lambda t,x: np.exp(x.alpha.values+lnR(t,x.ER.values)),grid); m=grid<=Tsup[g]
    axB.fill_between(grid[m],lo[m],hi[m],color=OI[g],alpha=0.12,lw=0); axB.plot(grid[m],med[m],color=OI[g],lw=1.6,zorder=3)
    if (~m).any(): mx=grid>=Tsup[g]; axB.plot(grid[mx],med[mx],color=OI[g],lw=1.1,ls=(0,(2,1.5)),zorder=3)
    p=P[P.group==g]; axB.scatter(p["T"]+np.random.RandomState(0).uniform(-.25,.25,len(p)),p.respiration_fgC_h,s=4,color=OI[g],alpha=0.35,lw=0,zorder=1)
axB.set_yscale("log"); axB.yaxis.set_major_locator(FixedLocator([300,1000,3000])); axB.yaxis.set_major_formatter(FixedFormatter(["300","1000","3000"])); axB.yaxis.set_minor_formatter(NullFormatter()); axB.set_xlabel("temperature (°C)"); axB.set_ylabel("per-cell respiration (fg C cell⁻¹ h⁻¹)"); axB.set_xticks([24,28,32,36,40,44])
axB.text(0.02,0.95,"dashed: beyond the range in which at least half the group's vials grew",transform=axB.transAxes,fontsize=6.0,color=MUTED,va="top")
panel_label(axB,"b",x=-0.2)
# ---- c: activation energies ----
order=GROUPS
for i,g in enumerate(order):
    v=vals.loc[g]; y=len(order)-1-i
    axC.plot([v.ER_lo,v.ER_hi],[y,y],color=OI[g],lw=1.4); axC.scatter(v.ER,y,marker="^",s=34,color=OI[g],zorder=3)
    axC.plot([v.E_lo,v.E_hi],[y,y],color=OI[g],lw=1.4); axC.scatter(v.E,y,marker="o",s=30,color=OI[g],zorder=3)
axC.axvline(0.65,color=AMB,lw=0.8,ls=(0,(2,2))); axC.text(0.65,len(order)-0.35,"0.65 eV, MTE reference",ha="center",fontsize=6.2,color=MUTED)
axC.set_yticks(range(len(order))); axC.set_yticklabels([SHORT[g] for g in order][::-1],style="italic"); axC.tick_params(axis="y",length=0); axC.spines["left"].set_visible(False)
axC.set_xlabel("activation energy (eV)"); axC.set_xlim(0.15,1.4); axC.set_ylim(-0.6,len(order)-0.2)
from matplotlib.lines import Line2D
axC.legend(handles=[Line2D([],[],marker="o",ls="",color=INK,ms=5,label="growth"),Line2D([],[],marker="^",ls="",color=INK,ms=5,label="respiration")],loc="lower right",fontsize=6.4)
panel_label(axC,"c",x=-0.09)
fig.legend(handles=[Line2D([],[],color=OI[g],lw=2,label=SHORT[g]) for g in GROUPS],loc="lower center",ncol=4,bbox_to_anchor=(0.5,-0.02),fontsize=6.4,handlelength=1.6)
save(fig,"FIG1_decoupling_py")
