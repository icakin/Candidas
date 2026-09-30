#!/usr/bin/env python3
"""Supplementary - phylogenetic distribution of thermal traits and homolog counts.

Run from anywhere:
    python3 ~/Desktop/Projects/Candidas/phylo/make_fig_supp_phylogenomics.py

Needs: matplotlib, pandas, numpy, biopython  (all in the conda base env)
Reads : results/tables/, phylo/pg_noCl.nwk
Writes: results/figures/manuscript/FIG_SUPP_phylogenomics.png
"""
from pathlib import Path
PROJ   = Path(__file__).resolve().parent.parent
TABLES = PROJ / "results" / "tables"
FIGS   = PROJ / "results" / "figures" / "manuscript"
PHYLO  = PROJ / "phylo"
FIGS.mkdir(parents=True, exist_ok=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from Bio import Phylo

fv=pd.read_csv(TABLES/"fig_values.csv").set_index("Group")
cen=pd.read_csv(PHYLO/"thermal_census.csv",index_col=0)
OI={"Clade1":"#0072B2","Clade2":"#CC79A7","Clade3":"#009E73","Clade4":"#E69F00",
    "para":"#D55E00","Hae":"#984EA3","Duo":"#00A6A6"}
G={"auris_cladeI":"Clade1","auris_cladeII":"Clade2","auris_cladeIII":"Clade3",
   "auris_cladeIV":"Clade4","haemulonii":"Hae","duobushaemulonii":"Duo","parapsilosis":"para"}
NICE={"auris_cladeI":"C. auris clade I","auris_cladeII":"C. auris clade II",
      "auris_cladeIII":"C. auris clade III","auris_cladeIV":"C. auris clade IV",
      "haemulonii":"C. haemulonii","duobushaemulonii":"C. duobushaemulonii",
      "pseudohaemulonii":"C. pseudohaemulonii","parapsilosis":"C. parapsilosis",
      "albicans":"C. albicans"}

# rooted 9-taxon tree (C. lusitaniae removed: its long branch attracted to the
# Debaryomycetaceae, see caption)
t=Phylo.read(str(PHYLO/"pg_noCl.nwk"),"newick")
t.root_with_outgroup({"name":"albicans"},{"name":"parapsilosis"})
t.ladderize(reverse=True)
tips=[x.name for x in t.get_terminals()]
Y={x.name:i for i,x in enumerate(t.get_terminals())}
def yof(c): return Y[c.name] if c.is_terminal() else float(np.mean([yof(k) for k in c.clades]))

fig=plt.figure(figsize=(16.0,6.4))
gs=fig.add_gridspec(1,4,width_ratios=[1.85,.74,.74,1.52],wspace=.34)
axT,axA,axB,axC=[fig.add_subplot(gs[i]) for i in range(4)]

def draw(cl,x0):
    x1=x0+(cl.branch_length or 0); y=yof(cl)
    axT.plot([x0,x1],[y,y],color="#333",lw=1.5,solid_capstyle="round")
    if not cl.is_terminal():
        ys=[yof(k) for k in cl.clades]
        axT.plot([x1,x1],[min(ys),max(ys)],color="#333",lw=1.5,solid_capstyle="round")
        for k in cl.clades: draw(k,x1)
draw(t.root,0.0)
dep=t.depths(); xmax=max(dep[x] for x in t.get_terminals())
for tip in t.get_terminals():
    nm=tip.name; y=Y[nm]; g=G.get(nm)
    axT.plot([dep[tip],xmax*1.02],[y,y],ls=":",lw=.6,color="#c9c9c9")
    axT.text(xmax*1.05,y,NICE[nm],fontsize=8.5,va="center",style="italic",
             fontweight="bold" if g else "normal",color=OI[g] if g else "#909090")
aur=[Y[k] for k in Y if k.startswith("auris")]
BX=xmax*1.92
axT.plot([BX,BX],[min(aur)-.22,max(aur)+.22],color="#555",lw=2.4,solid_capstyle="butt")
axT.text(BX+xmax*0.04,float(np.mean(aur)),
         "4 clades, ~25× closer to each other\nthan to any other species here",
         fontsize=7.1,va="center",color="#444")
axT.set_xlim(-.02,xmax*2.75); axT.set_ylim(-1.35,len(tips)-.15)
axT.invert_yaxis(); axT.axis("off")
axT.set_title("a   Phylogenomic tree",loc="left",fontweight="bold",fontsize=10.5)
axT.text(0,-.78,"600 single-copy orthologs · 289,787 aa · ML (LG+Γ)\nall nodes 1.00 · rooted on C. albicans + C. parapsilosis",
         fontsize=7,color="#555",va="top")
axT.plot([0,0.1],[len(tips)-.40]*2,color="#333",lw=1.6)
axT.text(.05,len(tips)-.58,"0.1 subs/site",fontsize=6.6,ha="center",color="#333")

def trait(ax,lo,mid,hi,title,xlab,log=False,ref=None,skip=()):
    for nm,y in Y.items():
        g=G.get(nm)
        if not g:
            ax.text(.5,y,"not phenotyped",transform=ax.get_yaxis_transform(),fontsize=6.3,
                    color="#d0d0d0",ha="center",va="center",style="italic"); continue
        if g in skip:
            ax.text(.5,y,"no detectable growth\nat 40 °C",transform=ax.get_yaxis_transform(),
                    fontsize=6.1,color="#666",ha="center",va="center",style="italic"); continue
        m,l,h=fv.loc[g,mid],fv.loc[g,lo],fv.loc[g,hi]
        ax.plot([l,h],[y,y],color=OI[g],lw=2.2,solid_capstyle="round",alpha=.92)
        ax.plot(m,y,"o",color=OI[g],ms=5.6,mec="white",mew=.8,zorder=3)
    if ref is not None: ax.axvline(ref,color="#B22222",ls="--",lw=.95)
    if log: ax.set_xscale("log",base=2)
    ax.set_ylim(-1.35,len(tips)-.15); ax.invert_yaxis(); ax.set_yticks([])
    ax.spines[["left","right","top"]].set_visible(False)
    ax.grid(axis="x",lw=.3,color="#eee"); ax.tick_params(labelsize=7.6)
    ax.set_xlabel(xlab,fontsize=8); ax.set_title(title,loc="left",fontweight="bold",fontsize=10.5)

trait(axA,"Tcue_lo","Tcue","Tcue_hi","b   CUE optimum","°C",ref=37)
trait(axB,"fc_lo","fever_cost","fc_hi","c   R/G, 37→40 °C","fold change",log=True,ref=1,skip=("Duo",))
axB.set_xlim(.95,3); axB.set_xticks([1,2]); axB.set_xticklabels(["1×","2×"])

cols=["HSP90","HSP70","sHSP","Trehalose","Calcineurin","Desaturase","Ergosterol","altNADH","AOX"]
M=cen.reindex(tips)[cols]
shade=np.zeros(M.shape); shade[:,cols.index("AOX")]=(M["AOX"].values-1)*0.85
axC.imshow(shade,cmap="Reds",vmin=0,vmax=1.5,aspect="auto")
for i in range(M.shape[0]):
    for j,c in enumerate(cols):
        v=int(M.values[i,j]); isaox=(c=="AOX")
        axC.text(j,i,str(v),ha="center",va="center",fontsize=8.3,
                 color="#B22222" if (isaox and v>1) else "#333",
                 fontweight="bold" if isaox else "normal")
axC.set_xticks(range(len(cols))); axC.set_xticklabels(cols,rotation=40,ha="right",fontsize=7.5)
axC.set_yticks([]); axC.set_ylim(len(tips)-.5,-.5)
for sp in axC.spines.values(): sp.set_visible(False)
axC.set_title("d   Candidate thermal-protein homolog counts",loc="left",fontweight="bold",fontsize=10.5)
axC.text(0,-.17,
  "Unique target proteins per family (DIAMOND, e<1e-20, \u226540 % id, \u226560 % cov).  Largely conserved;\n"
  "two AOX homologs in C. albicans and C. parapsilosis, one elsewhere \u2014 tracks family, not performance.",
  transform=axC.transAxes,fontsize=6.9,color="#666",va="top")

fig.suptitle("Phylogenetic distribution of thermal traits and candidate-protein homolog counts",
             fontsize=12.6,fontweight="bold",x=.008,ha="left",y=1.015)
fig.text(.008,-.045,
 "Genomes are lineage reference strains, not the assayed isolates.  C. lusitaniae excluded (long-branch attraction).",
 fontsize=7.2,color="#666",ha="left",va="top")
plt.savefig(str(FIGS/"FIG_SUPP_phylogenomics.png"),dpi=300,bbox_inches="tight",facecolor="white")
print("saved")
