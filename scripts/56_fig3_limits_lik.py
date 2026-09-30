"""Upper thermal limits from the per-well likelihood model (results/tables/lik_transition.csv,
written by 56_lik_transition.py), aligned to a to-scale phylogeny. Derived from 53_fig_limits_phylo.py;
the tree code is unchanged, only the right-hand panel's data source and encoding differ."""
import os, sys, re, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
from matplotlib.lines import Line2D

# ---------------- newick -> cumulative root distance per tip -------------------
def parse_newick(s):
    s=s.strip().rstrip(";"); pos=[0]
    def node():
        ch=[]
        if s[pos[0]]=="(":
            pos[0]+=1
            while True:
                ch.append(node())
                if s[pos[0]]==",": pos[0]+=1; continue
                if s[pos[0]]==")": pos[0]+=1; break
        m=re.match(r"[^,()\:]*",s[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        bl=0.0
        if pos[0]<len(s) and s[pos[0]]==":":
            pos[0]+=1; m=re.match(r"[0-9.eE+-]+",s[pos[0]:]); bl=float(m.group(0)); pos[0]+=len(m.group(0))
        return {"label":lab.strip(),"bl":bl,"children":ch}
    return node()

def depths(nd,acc=0.0,out=None,path=None):
    out={} if out is None else out; path=[] if path is None else path
    d=acc+nd["bl"]
    if not nd["children"]: out[nd["label"]]=(d,tuple(path)); return out
    for i,c in enumerate(nd["children"]): depths(c,d,out,path+[i])
    return out

root=parse_newick(open(f"{C}/phylo/trees/pg_rooted.nwk").read())
DEP=depths(root)   # tip -> (root distance, path)

KEY={"Clade1":"auris_cladeI","Clade2":"auris_cladeII","Clade3":"auris_cladeIII","Clade4":"auris_cladeIV",
     "Duo":"duobushaemulonii","Hae":"haemulonii","para":"parapsilosis"}
NAME=SHORT

# internal-node depth = root distance of the deepest shared prefix of two tips
def mrca_depth(tips):
    paths=[DEP[KEY[g]][1] for g in tips]
    k=0
    while all(len(p)>k for p in paths) and len({p[k] for p in paths})==1: k+=1
    nd=root; acc=0.0
    for i in paths[0][:k]:
        acc+=nd["bl"]; nd=nd["children"][i]
    return acc+nd["bl"]

t=pd.read_csv(f"{C}/results/tables/lik_transition.csv"); R=t.set_index("isolate")
rows=[]
for g in TREE_ORDER:
    for iso in sorted(t[t.group==g].isolate.unique()): rows.append((g,iso))
ypos={iso:i for i,(g,iso) in enumerate(rows)}
tax_y={g:np.mean([ypos[i] for gg,i in rows if gg==g]) for g in TREE_ORDER}
tip_x={g:DEP[KEY[g]][0] for g in TREE_ORDER}

fig,(axT,axL)=plt.subplots(1,2,figsize=(8.4,4.6),gridspec_kw=dict(width_ratios=[0.62,1],wspace=0.34))

def draw_tree(ax,xmin,xmax,label_tips=True,lw=0.9,fs=6.2):
    """rectangular phylogram, x = substitutions/site from the root"""
    def clip(x): return min(max(x,xmin),xmax)
    def rec(tips):
        if len(tips)==1:
            g=tips[0]
            ax.plot([clip(mrca_depth(tips)),clip(tip_x[g])],[tax_y[g]]*2,color=MUTED,lw=lw)
            return tax_y[g], mrca_depth(tips)
        return None
    # explicit clade structure, drawn from real depths
    groups=[(["Clade1"],),(["Clade3"],),(["Clade2"],),(["Clade4"],),(["Duo"],),(["Hae"],),(["para"],)]
    joins=[(["Clade1","Clade3"],),(["Clade1","Clade3","Clade2"],),(["Clade1","Clade3","Clade2","Clade4"],),
           (["Duo","Hae"],),(["Clade1","Duo"],),(["Clade1","para"],)]
    # terminal branches
    for g in TREE_ORDER:
        par = {"Clade1":["Clade1","Clade3"],"Clade3":["Clade1","Clade3"],
               "Clade2":["Clade1","Clade3","Clade2"],"Clade4":["Clade1","Clade3","Clade2","Clade4"],
               "Duo":["Duo","Hae"],"Hae":["Duo","Hae"],"para":["Clade1","para"]}[g]
        ax.plot([clip(mrca_depth(par)),clip(tip_x[g])],[tax_y[g]]*2,color=MUTED,lw=lw)
    def node(tips,parent_tips):
        x=clip(mrca_depth(tips)); xp=clip(mrca_depth(parent_tips))
        ys=[tax_y[g] for g in tips]
        ax.plot([xp,x],[np.mean(ys)]*2,color=MUTED,lw=lw)
        return x,np.mean(ys)
    def vbar(tips_a,tips_b,tips_join):
        x=clip(mrca_depth(tips_join))
        ya=np.mean([tax_y[g] for g in tips_a]); yb=np.mean([tax_y[g] for g in tips_b])
        ax.plot([x,x],[ya,yb],color=MUTED,lw=lw)
    vbar(["Clade1"],["Clade3"],["Clade1","Clade3"])
    node(["Clade1","Clade3"],["Clade1","Clade3","Clade2"])
    vbar(["Clade1","Clade3"],["Clade2"],["Clade1","Clade3","Clade2"])
    node(["Clade1","Clade3","Clade2"],["Clade1","Clade3","Clade2","Clade4"])
    vbar(["Clade1","Clade3","Clade2"],["Clade4"],["Clade1","Clade3","Clade2","Clade4"])
    node(["Clade1","Clade3","Clade2","Clade4"],["Clade1","Duo"])
    vbar(["Duo"],["Hae"],["Duo","Hae"])
    node(["Duo","Hae"],["Clade1","Duo"])
    vbar(["Clade1","Clade3","Clade2","Clade4"],["Duo","Hae"],["Clade1","Duo"])
    node(["Clade1","Duo"],["Clade1","para"])
    vbar(["Clade1","Clade2","Clade3","Clade4","Duo","Hae"],["para"],["Clade1","para"])
    ax.plot([clip(mrca_depth(["Clade1","para"])-0.03),clip(mrca_depth(["Clade1","para"]))],
            [np.mean([tax_y[g] for g in TREE_ORDER])]*2,color=MUTED,lw=lw)
    if label_tips:
        for g in TREE_ORDER:
            ax.plot([clip(tip_x[g]),0.985],[tax_y[g]]*2,color="#dedcd8",lw=0.5,ls=(0,(1,2)),zorder=0)
            ax.text(0.978,tax_y[g],NAME[g],fontsize=fs,style="italic",ha="right",va="center",color=INK)

draw_tree(axT,0.20,0.78)
# connect each group tip to its isolate rows
for g in TREE_ORDER:
    ys=[ypos[i] for gg,i in rows if gg==g]
    axT.plot([1.02,1.02],[min(ys),max(ys)],color=AMB,lw=0.7)
    axT.plot([tip_x[g]+0.115,1.02],[tax_y[g]]*2,color=AMB,lw=0.7)
    for y in ys: axT.plot([1.02,1.09],[y,y],color=AMB,lw=0.7)
axT.set_xlim(0.20,1.10); axT.set_ylim(len(rows)-0.4,-1.6); axT.axis("off")
# scale bar
axT.plot([0.22,0.32],[-1.25,-1.25],color=INK,lw=1.1)
axT.text(0.27,-1.42,"0.1 substitutions/site",fontsize=5.9,ha="center",color=MUTED)
axT.text(0.20,1.04,"a",transform=axT.transAxes,fontweight="bold",fontsize=9)

# ---------------- right: per-isolate limits ----------------------------------
for (g,iso) in rows:
    y=ypos[iso]; c=OI[g]
    r=R.loc[iso]
    if r.growth_censored:
        axL.plot([r.T_last_growth,44],[y,y],color=c,lw=1.6,zorder=3)
        axL.scatter(44,y,s=30,marker=">",color=c,zorder=4); axL.plot([44,45.6],[y,y],color=c,lw=1.6,zorder=3)
    else:
        # grid interval: growth still present at T_last_growth, lost by T_growthloss
        axL.plot([r.T_last_growth,r.T_growthloss],[y,y],color=AMB,lw=4.0,solid_capstyle="butt",zorder=2)
        axL.scatter(r.T_growthloss,y,s=30,color=c,edgecolor="white",linewidth=0.6,zorder=4)
axL.axvspan(36.8,37.2,color=INK,alpha=0.10,lw=0,zorder=0)
axL.axvspan(39.8,40.2,color=INK,alpha=0.06,lw=0,zorder=0)
axL.axvline(44,color=MUTED,lw=0.8,ls=(0,(3,2)),zorder=1)
axL.set_yticks(range(len(rows))); axL.set_yticklabels([iso_label(iso) for g,iso in rows],fontsize=6.0,style="italic")
for lab,(g,iso) in zip(axL.get_yticklabels(),rows): lab.set_color(OI[g])
axL.set_ylim(len(rows)-0.4,-1.6)
axL.set_xlim(34.5,46.2); axL.set_xticks([36,38,40,42,44])
axL.set_xlabel("temperature at which growth is lost (°C)",fontsize=7.4)
axL.text(37,-0.85,"body",ha="center",fontsize=6.2,color=MUTED)
axL.text(40,-0.85,"fever",ha="center",fontsize=6.2,color=MUTED)
axL.text(44.15,-0.85,"assay ceiling",ha="left",fontsize=6.0,color=MUTED)
for sp in ("top","right","left"): axL.spines[sp].set_visible(False)
axL.tick_params(axis="y",length=0); axL.tick_params(axis="x",length=2.5)
axL.text(-0.26,1.04,"b",transform=axL.transAxes,fontweight="bold",fontsize=9)
h=[Line2D([],[],color=AMB,lw=4,label="point = first temperature at which fewer than half the vials grow;\nbar = the 2 °C grid step over which growth was lost"),
   Line2D([],[],marker=">",color=MUTED,ls="",ms=5,label="still growing at 44 °C: limit is >44 °C, never a number")]
fig.legend(handles=h,frameon=False,fontsize=6.2,ncol=1,loc="upper center",bbox_to_anchor=(0.68,1.10))
save(fig,"FIG_thermal_limits_lik")
for g in TREE_ORDER: print(f"  {NAME[g]:24s} root distance {tip_x[g]:.4f}")
print("  max within C. auris (I-IV): %.4f" % (tip_x["Clade1"]+tip_x["Clade4"]-2*mrca_depth(["Clade1","Clade4"])))
print("  C. auris I to C. haemulonii: %.4f" % (tip_x["Clade1"]+tip_x["Hae"]-2*mrca_depth(["Clade1","Hae"])))
print("  C. auris I to C. parapsilosis: %.4f" % (tip_x["Clade1"]+tip_x["para"]-2*mrca_depth(["Clade1","para"])))
