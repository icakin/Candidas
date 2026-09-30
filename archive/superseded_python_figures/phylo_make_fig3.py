#!/usr/bin/env python3
"""Figure 3 - febrile persistence and conditional respiratory cost (final).

Run from anywhere:
    python3 ~/Desktop/Projects/Candidas/phylo/make_fig3.py

Needs: matplotlib, pandas, numpy, biopython  (all in the conda base env)
Reads : results/tables/, phylo/pg_noCl.nwk
Writes: results/figures/manuscript/FIG3_persistence_cost.png
        phylo/fig3_tree_layout.csv  (tree geometry, consumed by scripts/15_fig3.R)
"""
from pathlib import Path
PROJ   = Path(__file__).resolve().parent.parent
TABLES = PROJ / "results" / "tables"
FIGS   = PROJ / "results" / "figures" / "manuscript"
PHYLO  = PROJ / "phylo"
FIGS.mkdir(parents=True, exist_ok=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.patches import Rectangle, Polygon
from Bio import Phylo

W=pd.read_csv(TABLES/"wells_by_isolate.csv",index_col=0)
W.columns=[int(c) for c in W.columns]
CT=pd.read_csv(TABLES/"carbon_tax_isolate.csv")
w=CT.pivot_table(index="Isolate",columns="T_C",values="tax"); w["fever"]=w[40]/w[37]
TEMPS=[36,38,40,42,44]

ISO=(["Clade%d_%d"%(c,n) for c,ns in
      [(1,[2068,2069,2070]),(2,[2071,2072,2073]),(3,[2074,2075,2076]),(4,[2077,2078,2079])] for n in ns]
     +["Hae_1724","Hae_1768","Hae_1769","Duo_1770","Duo_1771",
       "para_2051","para_2052","para_2053"])
ROW={}
for i,s in enumerate(ISO[:12]):  ROW[s]=i
for i,s in enumerate(ISO[12:15]):ROW[s]=12+i
for i,s in enumerate(ISO[15:17]):ROW[s]=15+i
for i,s in enumerate(ISO[17:]):  ROW[s]=18+i
NROW=22
# final palette: Hae violet, Duo teal (matches OI in scripts/14_main_figures.R)
COL={"Clade1":"#0072B2","Clade2":"#0072B2","Clade3":"#0072B2","Clade4":"#0072B2",
     "Hae":"#984EA3","Duo":"#00A6A6","para":"#D55E00"}
SPBLK={"C. auris":(0,11,"#0072B2"),"C. haemulonii":(12,14,"#984EA3"),
       "C. duobushaemulonii":(15,16,"#00A6A6"),"C. parapsilosis":(18,20,"#D55E00")}
TIPY={"auris":5.5,"haemulonii":13,"duobushaemulonii":15.5,
      "pseudohaemulonii":17,"parapsilosis":19,"albicans":21}

# ---------------- tree geometry ----------------
t=Phylo.read(str(PHYLO/"pg_noCl.nwk"),"newick")
t.root_with_outgroup({"name":"albicans"},{"name":"parapsilosis"})
dep=t.depths()
aur=[x for x in t.get_terminals() if x.name.startswith("auris")]
AUR_X0=dep[t.common_ancestor(*aur)]
D={x.name:dep[x] for x in t.get_terminals()}
def node(*names):
    return dep[t.common_ancestor(*[x for x in t.get_terminals() if x.name in names])]
NODE_CAND=node("auris_cladeI","haemulonii")
NODE_HC  =node("haemulonii","duobushaemulonii","pseudohaemulonii")
NODE_DP  =node("duobushaemulonii","pseudohaemulonii")
NODE_DEB =node("albicans","parapsilosis")
XM=max(D.values()); AUR_TIP=AUR_X0+0.070*XM

fig=plt.figure(figsize=(13.4,6.9))
gs=fig.add_gridspec(1,3,width_ratios=[1.48,2.48,0.52],wspace=.02,top=.93)
axT=fig.add_subplot(gs[0]); axM=fig.add_subplot(gs[1]); axF=fig.add_subplot(gs[2])

# ---------------- a: tree ----------------
SEGS=[]                                   # collected for the R port
def hline(x0,x1,y,c="#333"): SEGS.append(("seg",x0,x1,y,y,c))
def vline(x,y0,y1,c="#333"): SEGS.append(("seg",x,x,y0,y1,c))
mid_c=(TIPY["auris"]+TIPY["pseudohaemulonii"])/2
mid_d=(TIPY["parapsilosis"]+TIPY["albicans"])/2
vline(0,mid_c,mid_d);            hline(0,NODE_CAND,mid_c);  hline(0,NODE_DEB,mid_d)
mid_h=(TIPY["haemulonii"]+TIPY["pseudohaemulonii"])/2
vline(NODE_CAND,TIPY["auris"],mid_h); hline(NODE_CAND,AUR_X0,TIPY["auris"])
hline(NODE_CAND,NODE_HC,mid_h)
mid_dp=(TIPY["duobushaemulonii"]+TIPY["pseudohaemulonii"])/2
vline(NODE_HC,TIPY["haemulonii"],mid_dp); hline(NODE_HC,D["haemulonii"],TIPY["haemulonii"])
hline(NODE_HC,NODE_DP,mid_dp)
vline(NODE_DP,TIPY["duobushaemulonii"],TIPY["pseudohaemulonii"])
hline(NODE_DP,D["duobushaemulonii"],TIPY["duobushaemulonii"])
hline(NODE_DP,D["pseudohaemulonii"],TIPY["pseudohaemulonii"],c="#b0b0b0")
vline(NODE_DEB,TIPY["parapsilosis"],TIPY["albicans"])
hline(NODE_DEB,D["parapsilosis"],TIPY["parapsilosis"])
hline(NODE_DEB,D["albicans"],TIPY["albicans"],c="#b0b0b0")
for k,x0,x1,y0,y1,c in SEGS:
    axT.plot([x0,x1],[y0,y1],color=c,lw=1.7,solid_capstyle="round",zorder=3)
axT.add_patch(Polygon([[AUR_X0,TIPY["auris"]],
                       [AUR_TIP,TIPY["auris"]-.62],[AUR_TIP,TIPY["auris"]+.62]],
              closed=True,facecolor="#0072B2",alpha=.32,edgecolor="#0072B2",lw=1.1,zorder=3))
LBL={"auris":("C. auris","#0072B2",1),"haemulonii":("C. haemulonii","#984EA3",1),
     "duobushaemulonii":("C. duobushaemulonii","#00A6A6",1),
     "pseudohaemulonii":("C. pseudohaemulonii","#b0b0b0",0),
     "parapsilosis":("C. parapsilosis","#D55E00",1),"albicans":("C. albicans","#b0b0b0",0)}
for k,(txt,c,meas) in LBL.items():
    x=AUR_TIP if k=="auris" else D[k]
    axT.plot([x,XM*1.03],[TIPY[k],TIPY[k]],ls=":",lw=.6,color="#d0d0d0",zorder=2)
    axT.text(XM*1.06,TIPY[k],txt,fontsize=9.0,va="center",style="italic",
             fontweight="bold" if meas else "normal",color=c)
    if not meas:
        axT.text(XM*1.06,TIPY[k]+.62,"not phenotyped",fontsize=6.3,va="center",color="#bdbdbd")
axT.text(AUR_TIP,TIPY["auris"]+1.05,"4 clade references collapsed",
         fontsize=6.2,color="#0072B2",ha="left",va="center")
# species-level persistence counts, aligned with the tips
GRP_ROWS={"auris":[i for i in ISO if i.startswith("Clade")],
          "Hae":[i for i in ISO if i.startswith("Hae")],
          "Duo":[i for i in ISO if i.startswith("Duo")],
          "para":[i for i in ISO if i.startswith("para")]}
CNT_Y={"auris":TIPY["auris"],"Hae":TIPY["haemulonii"],
       "Duo":TIPY["duobushaemulonii"],"para":TIPY["parapsilosis"]}
CNT_C={"auris":"#0072B2","Hae":"#984EA3","Duo":"#00A6A6","para":"#D55E00"}
XC=XM*3.02
for gk,rows_ in GRP_ROWS.items():
    n=sum(int(W.loc[i,40])>0 for i in rows_)
    axT.text(XC,CNT_Y[gk],f"{n}/{len(rows_)}",fontsize=8.6,ha="center",va="center",
             color=CNT_C[gk],fontweight="bold")
axT.text(XC,-1.45,"growing\nat 40 °C",fontsize=6.8,ha="center",va="center",
         color="#444",fontweight="bold")
axT.set_xlim(-.04,XM*3.45); axT.set_ylim(-1.9,NROW+.35); axT.invert_yaxis(); axT.axis("off")
axT.set_title("a   Phylogenomic context",loc="left",fontweight="bold",fontsize=10.2)
axT.plot([0,.3],[NROW+.02]*2,color="#333",lw=1.7)
axT.text(0,NROW+.30,"0.3 substitutions/site",fontsize=6.6,ha="left",va="top",color="#333")

# export tree geometry for the R port (scripts/15_fig3.R)
rows=[dict(kind=k,x=a,xend=b,y=c_,yend=d_,text="",color=cc,bold=0)
      for k,a,b,c_,d_,cc in SEGS]
rows.append(dict(kind="wedge",x=AUR_X0,xend=AUR_TIP,y=TIPY["auris"]-.62,
                 yend=TIPY["auris"]+.62,text="",color="#0072B2",bold=0))
for k,(txt,c,meas) in LBL.items():
    rows.append(dict(kind="label",x=XM*1.06,xend=XM*1.03,
                     y=TIPY[k],yend=AUR_TIP if k=="auris" else D[k],
                     text=txt,color=c,bold=meas))
pd.DataFrame(rows).to_csv(PHYLO/"fig3_tree_layout.csv",index=False)

# ---------------- b: matrix + fused cost column ----------------
SC=["#f7dcdc","#e8e8ea","#cfd0d4","#adaeb5","#84868f","#686c77"]
for iso,y in ROW.items():
    g=iso.split("_")[0]
    for j,T in enumerate(TEMPS):
        n=int(W.loc[iso,T])
        axM.add_patch(Rectangle((j-.46,y-.42),.92,.84,facecolor=SC[n],
                      edgecolor="white",lw=1.1,zorder=3))
        axM.text(j,y,str(n),ha="center",va="center",fontsize=7.9,zorder=4,
                 color="white" if n>=4 else "#444",fontweight="bold" if n==0 else "normal")
    axM.text(-.72,y,iso.split("_")[1],ha="right",va="center",fontsize=7.5,color=COL[g])
for sp,(a,b,c) in SPBLK.items():
    axM.plot([-1.62,-1.62],[a-.40,b+.40],color=c,lw=2.6,solid_capstyle="butt",clip_on=False)
axM.axvline(1.5,color="#B22222",lw=1.5,ls="--",zorder=5)
axM.text(1.5,-1.15,"fever threshold: 40 °C",fontsize=7.2,color="#B22222",ha="center",
         bbox=dict(facecolor="white",edgecolor="none",pad=1.2))
axM.set_xlim(-1.75,len(TEMPS)-.4); axM.set_ylim(-1.9,NROW+.35); axM.invert_yaxis()
axM.set_xticks(range(len(TEMPS))); axM.set_xticklabels([f"{t} °C" for t in TEMPS],fontsize=8.4)
axM.set_yticks([]); axM.xaxis.set_ticks_position("bottom")
for s in axM.spines.values(): s.set_visible(False)
axM.set_title("b   Persistence — growth-positive wells per isolate (of 5)",
              loc="left",fontweight="bold",fontsize=10.2)

# cost column: printed fold changes, an annotation not an analysis
for iso,y in ROW.items():
    g=iso.split("_")[0]
    if int(W.loc[iso,40])>0:
        axF.text(.5,y,f"{w.loc[iso,'fever']:.2f}×",ha="center",va="center",
                 fontsize=7.6,color=COL[g],fontweight="bold")
    else:
        axF.text(.5,y,"—",ha="center",va="center",fontsize=7.6,color="#c4c4c4")
axF.text(.5,-1.45,"R/G fold change\n37→40 °C",fontsize=6.8,color="#444",ha="center",
         va="center",fontweight="bold")
axF.set_xlim(0,1); axF.set_ylim(-1.9,NROW+.35); axF.invert_yaxis(); axF.axis("off")

fig.text(.008,.975,r"Consistent febrile growth in $\it{C.\ auris}$ carries an increased respiratory cost",
         fontsize=12.6,fontweight="bold",ha="left",va="bottom")
fig.text(.80,.055,"— = no detectable growth at 40 °C",
         fontsize=6.9,color="#777",ha="right",va="top")
plt.savefig(FIGS/"FIG3_persistence_cost.png",dpi=300,bbox_inches="tight",facecolor="white")
print("saved", FIGS/"FIG3_persistence_cost.png")
