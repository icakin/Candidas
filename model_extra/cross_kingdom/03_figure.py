"""Cross-kingdom figure, model-free. (a) T_carbon vs T_opt per taxon; (b) respiration per
unit growth relative to its minimum against T - T_opt(growth); (c) state of wells above
each taxon's growth limit (respiring with growth below detection vs inert)."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BLUE, ORANGE, INK, MUTED, FAINT = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e8e8e6"
C=os.path.expanduser("~/mnt/Candidas"); D=f"{C}/model_extra/cross_kingdom"; FIG=f"{C}/results/figures/manuscript"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],"font.size":7.5,
  "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"axes.edgecolor":MUTED,"xtick.color":MUTED,
  "ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
t=pd.read_csv(f"{D}/cross_kingdom_modelfree.csv"); cv=pd.read_csv(f"{D}/cross_kingdom_curves.csv")
col={"Bacteria":BLUE,"Fungi":ORANGE}
# ---- states above the growth limit -------------------------------------------------
names={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV","para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
cs=pd.read_csv(f"{C}/results/tables/fig4_well_states.csv"); cs["taxon"]=cs.group.map(names); cs["kingdom"]="Fungi"
rs=pd.read_csv(f"{D}/r2a_wells.csv"); rs["taxon"]=np.where(rs.strain.str.startswith("E."),rs.strain,"Soil "+rs.strain); rs["kingdom"]="Bacteria"
st=pd.concat([cs[["taxon","kingdom","T","state"]],rs[["taxon","kingdom","T","state"]]])
rows=[]
for tax,s in st.groupby("taxon"):
    g=s[s.state=="growing"].groupby("T").size(); lim=g[g>=3].index.max()
    a=s[s["T"]>lim]
    if len(a)==0: rows.append(dict(taxon=tax,kingdom=s.kingdom.iloc[0],limit=lim,n=0,resp=0,inert=0,grow=0)); continue
    rows.append(dict(taxon=tax,kingdom=s.kingdom.iloc[0],limit=lim,n=len(a),resp=(a.state=="respiring").sum(),inert=(a.state=="inert").sum(),grow=(a.state=="growing").sum()))
sa=pd.DataFrame(rows); sa["f_resp"]=np.where(sa.n>0,(sa.resp+sa.grow)/sa.n.replace(0,np.nan),-1); sa.to_csv(f"{D}/states_above_limit.csv",index=False)
print(sa.to_string())
# ---- figure --------------------------------------------------------------------------
fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=0.42,wspace=0.32)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,:])
# (a)
for _,r in t.iterrows():
    c=col[r.kingdom]; mk="o" if r.dataset.startswith("Candidas") or r.kingdom=="Bacteria" else "s"
    j=np.random.default_rng(hash(r.taxon)%1000).uniform(-0.5,0.5)
    axA.plot([r.T_opt+j,r.T_opt+j],[r.T_carbon_lo,r.T_carbon_hi],color=c,lw=0.7,alpha=0.55,zorder=2)
    axA.plot([r.T_opt_lo+j,r.T_opt_hi+j],[r.T_carbon,r.T_carbon],color=c,lw=0.7,alpha=0.55,zorder=2)
    axA.scatter(r.T_opt+j,r.T_carbon,s=22,color=c,marker=mk,edgecolor="white",linewidth=0.6,zorder=3)
axA.plot([10,50],[10,50],color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1); axA.text(47.5,46.2,"1:1",fontsize=6.4,color=MUTED,ha="right",va="top")
axA.set_xlim(10,50); axA.set_ylim(10,50); axA.set_xlabel("Growth optimum (°C)"); axA.set_ylabel("Carbon-economy optimum (°C)")
axA.set_xticks([10,20,30,40,50]); axA.set_yticks([10,20,30,40,50])
nb=int((t.gap>0).sum()); axA.text(0.03,0.97,f"below the line: {nb} of {len(t)} taxa",transform=axA.transAxes,fontsize=6.8,va="top")
# (b)
for tax,s in cv.groupby("taxon"):
    s=s.sort_values("T"); axB.plot(s.dT,s.q_rel,color=col[s.kingdom.iloc[0]],lw=1.0,alpha=0.7,zorder=2)
axB.axvline(0,color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1); axB.set_yscale("log")
axB.set_xlabel("Temperature relative to growth optimum (°C)"); axB.set_ylabel("Respiration per unit growth\n(relative to its minimum)")
axB.set_yticks([1,2,5,10,20,50]); axB.set_yticklabels(["1","2","5","10","20","50"]); axB.set_xlim(-22,16)
axB.text(0.6,0.97,"growth\noptimum",transform=axB.get_xaxis_transform(),ha="left",fontsize=6.4,color=MUTED,va="top")
# (c)
sa=sa.sort_values(["kingdom","f_resp","taxon"],ascending=[False,False,True]).reset_index(drop=True)
y=np.arange(len(sa))
for i,r in sa.iterrows():
    c=col[r.kingdom]
    if r.n==0:
        axC.text(0.01,i,f"grows at every temperature measured (to {int(r.limit)} °C)",va="center",fontsize=6.2,color=MUTED,style="italic"); continue
    fr=(r.resp+r.grow)/r.n
    axC.barh(i,fr,color=c,height=0.66,zorder=2); axC.barh(i,1-fr,left=fr+0.004,color=FAINT,height=0.66,zorder=2)
    axC.text(1.02,i,f"{int(r.resp+r.grow)}/{int(r.n)}  (limit {int(r.limit)} °C)",va="center",fontsize=6.2,color=MUTED)
axC.set_yticks(y); axC.set_yticklabels(sa.taxon,fontsize=6.6); axC.invert_yaxis(); axC.set_xlim(0,1); axC.set_xticks([0,0.5,1]); axC.set_xticklabels(["0","50","100%"])
axC.set_xlabel("Wells above the growth limit that keep consuming oxygen (≥ 2 mg L⁻¹)")
for ax in (axA,axB,axC):
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(length=2.5)
axC.spines["left"].set_visible(False); axC.tick_params(axis="y",length=0)
h=[Line2D([],[],marker="o",color=BLUE,ls="",ms=5,label="Bacteria"),Line2D([],[],marker="o",color=ORANGE,ls="",ms=5,label="Fungi (Candida)"),Line2D([],[],marker="s",color=ORANGE,ls="",ms=5,label="Fungi (Zymoseptoria)")]
axA.legend(handles=h,frameon=False,fontsize=6.4,loc="upper left",bbox_to_anchor=(0.0,0.93),handletextpad=0.3)
for ax,l in ((axA,"a"),(axB,"b"),(axC,"c")): ax.text(-0.12 if ax is not axC else -0.16,1.04,l,transform=ax.transAxes,fontweight="bold",fontsize=9)
fig.savefig(f"{FIG}/FIG_cross_kingdom.png",dpi=300,bbox_inches="tight"); fig.savefig(f"{FIG}/FIG_cross_kingdom.pdf",bbox_inches="tight")
print("saved")
