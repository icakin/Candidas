"""One aesthetic for every Python figure in the manuscript, matched to the R figures
(fig_common.R): the same group palette (OI), the same short names, the same font stack
(Helvetica > Arial), the same ink/muted greys, the same panel-label style."""
import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=f"{C}/results/figures/manuscript"
# group palette, identical to OI in fig_common.R
OI={"Clade1":"#0072B2","Clade2":"#CC79A7","Clade3":"#009E73","Clade4":"#E69F00",
    "para":"#D55E00","Hae":"#984EA3","Duo":"#00A6A6"}
SHORT={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
       "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
GROUPS=["Clade1","Clade2","Clade3","Clade4","para","Hae","Duo"]          # figure order
TREE_ORDER=["Clade1","Clade3","Clade2","Clade4","Duo","Hae","para"]       # phylogeny order
ISOLATES_TREE=["Clade1_2068","Clade1_2069","Clade1_2070","Clade3_2074","Clade3_2075","Clade3_2076",
 "Clade2_2071","Clade2_2072","Clade2_2073","Clade4_2077","Clade4_2078","Clade4_2079",
 "Duo_1770","Duo_1771","Hae_1724","Hae_1768","Hae_1769","para_2051","para_2052","para_2053"]
INK,MUTED,FAINT,AMB="#0b0b0b","#52514e","#e8e8e6","#b9b7b2"
RED="#B22222"                       # body/fever markers, as in the R figures
GROW,RESP,NORESP="#2F5597","#B22222","#f2f1ef"   # well states (not taxa)
AURIS,REL="#0b0b0b","#8c8c8c"       # pooled species-level lines
def group_of(iso): return iso.split("_")[0]
def col(iso_or_group): return OI[iso_or_group if iso_or_group in OI else group_of(iso_or_group)]
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
 "font.size":7.5,"axes.titlesize":7.5,"axes.labelsize":7.5,"xtick.labelsize":6.8,"ytick.labelsize":6.8,
 "legend.fontsize":6.4,"axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
 "xtick.major.size":2.5,"ytick.major.size":2.5,"axes.edgecolor":"#4d4d4d","xtick.color":"#4d4d4d",
 "ytick.color":"#4d4d4d","axes.labelcolor":"#1a1a1a","text.color":INK,"axes.spines.top":False,
 "axes.spines.right":False,"legend.frameon":False,"savefig.facecolor":"white","pdf.fonttype":42})
def panel_label(ax,letter,x=-0.16,y=1.06): ax.text(x,y,letter,transform=ax.transAxes,fontweight="bold",fontsize=9,va="bottom")
def group_legend(ax_or_fig,groups=GROUPS,**kw):
    from matplotlib.lines import Line2D
    h=[Line2D([],[],marker="o",ls="",color=OI[g],ms=4.5,label=SHORT[g]) for g in groups]
    return ax_or_fig.legend(handles=h,**kw)
def save(fig,name,dpi=300):
    fig.savefig(f"{FIG}/{name}.png",dpi=dpi,bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.pdf",bbox_inches="tight")
    print("saved",name)
def iso_label(iso):
    """'Clade2_2073' -> 'C. auris II 2073', 'Hae_1724' -> 'C. haemulonii 1724' (for tick labels)."""
    g,n=iso.split("_"); return f"{SHORT.get(g,g)} {n}"
