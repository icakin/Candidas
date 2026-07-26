import numpy as np, matplotlib.pyplot as plt, textwrap
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

INK="#1a1a1a"; SUB="#5b5b5b"; GRND="#3d5a80"; GRNDF="#eaf0f6"; CORE="#2C6FB0"; COREF="#eaf2fb"
RED="#B22222"; GRN="#2e7d54"; GRNF="#e9f4ee"; OUTF="#fdf2e6"; OUT="#c9791f"
plt.rcParams.update({"font.family":"DejaVu Sans","text.color":INK,"mathtext.default":"it"})
fig,ax=plt.subplots(figsize=(8.7,11.7)); ax.set_xlim(0,100); ax.set_ylim(-3,100); ax.axis("off")

def box(x,y,w,h,ec,fc="white",lw=1.2,r=0.9):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.2,rounding_size={r}",fc=fc,ec=ec,lw=lw,zorder=2))
def T(x,y,s,size=8,c=INK,w="normal",ha="center",it=False,fam=None,va="center"):
    ax.text(x,y,s,fontsize=size,color=c,ha=ha,va=va,fontweight=w,style="italic" if it else "normal",
            family=fam,zorder=4,linespacing=1.25)
def arr(x0,y0,x1,y1,c="#555",lw=1.4,rad=0,style="-|>"):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle=style,mutation_scale=13,lw=lw,color=c,
                zorder=1,connectionstyle=f"arc3,rad={rad}"))
def stage(x,y,n,t,c):
    ax.add_patch(plt.matplotlib.patches.Circle((x,y),1.5,fc=c,ec="none",zorder=4))
    T(x,y,n,7.6,"white",w="bold"); T(x+3.2,y,t,10.5,c,w="bold",ha="left")
def tpc(cx,cy,w,h,pts=False,col=CORE):
    xs=np.linspace(0,1,100); ys=np.exp(-((xs-0.62)**2)/0.05)*(0.4+0.6*xs); ys/=ys.max()
    X=cx-w/2+xs*w; Y=cy-h/2+ys*h
    ax.plot(X,Y,color=col,lw=1.5,zorder=4); ax.plot([cx-w/2,cx+w/2],[cy-h/2,cy-h/2],color="#aaa",lw=0.7,zorder=3)
    ax.plot([cx-w/2,cx-w/2],[cy-h/2,cy+h/2],color="#aaa",lw=0.7,zorder=3)
    if pts:
        rng=np.random.default_rng(1); ix=np.linspace(8,95,11).astype(int)
        ax.scatter(X[ix],Y[ix]+rng.uniform(-.04,.04,11)*h,s=5,color="#333",zorder=5)

# ===== TITLE =====
T(50,98.4,"The sequence-grounded etc-GEM pipeline",13.5,w="bold",va="top")
T(50,95.6,"grounded inputs $\\rightarrow$ constrained simulation $\\rightarrow$ Bayesian calibration to measured growth",8.4,c=SUB,va="top")

# ===== STAGE 1: INPUTS =====
stage(4,91,"1","Grounded inputs   (fixed; from genome, structure & databases)",GRND)
xin=[3,22.4,41.8,61.2,80.6]; wob=17.4; yin=77.5; hin=11.0
inputs=[("Metabolic network","iRV973 GEM\n973 genes · 2,863 rxns"),
        ("Turnover  $k_{cat}$","per-EC medians\n(DLKcat / BRENDA)"),
        ("Enzyme MW + seq","per-clade proteome\n(.faa)"),
        ("Thermal envelope","$T_{opt}$ · $T_m$ · $\\Delta C_p$\n(meltome, TemStaPro)"),
        ("Proteome sectors","iBAQ $\\rightarrow$ $f_{metab}$,\n$P_{total}$, $\\sigma$")]
for (x,(t,s)) in zip(xin,inputs):
    box(x,yin,wob,hin,GRND,GRNDF)
    tt="\n".join(textwrap.wrap(t,16))
    T(x+wob/2,yin+hin-1.9,tt,7.1,GRND,w="bold",va="top")
    T(x+wob/2,yin+hin-4.9,s,6.5,SUB,va="top")
for x in xin: arr(x+wob/2,yin,50,70.5,c="#9bb0c6",lw=1.0,rad=0)

# ===== STAGE 2: etc-GEM CORE =====
stage(4,68.5,"2","etc-GEM   (three constraint layers, solved at each temperature)",CORE)
box(7,44,86,22,CORE,COREF,lw=1.6,r=1.2)
T(50,63.7,"Enzyme- & Temperature-Constrained genome-scale model",9.3,CORE,w="bold")
# three layers
ly=[("L1  Enzyme constraint (GECKO)","every reaction $j$ draws its enzyme from one shared protein pool"),
    ("L2  Temperature dependence","$k_{cat,j}(T)$ from MMRT: activity rises with $T$, then denatures near $T_{m,j}$"),
    ("L3  Proteome allocation","the pool is capped by the metabolic-sector budget")]
for i,(a,b) in enumerate(ly):
    yy=60.3-i*3.1
    T(10,yy,a,7.7,INK,w="bold",ha="left"); T(41,yy,b,7.2,SUB,ha="left")
# the constraint + objective
box(11,45.4,78,4.7,"#c9c9c9",fc="white",lw=1.0,r=0.7)
T(20,47.75,"solve  $\\max\\ \\mu(T)$   s.t.",8.6,INK,ha="center")
T(58,47.75,"$\\sum_j \\frac{v_j}{k_{cat,j}(T)}\\, MW_j \\ \\leq\\ P_{total}\\, f_{metab}\\, \\sigma$",9.2,INK,ha="center")

# ===== STAGE 3: SWEEP -> predicted TPC =====
arr(50,44,50,40.5,lw=1.5)
stage(4,38.5,"3","Temperature sweep",OUT)
T(30,34.5,"repeat over  $T = 20$–$46$ °C",8.2,ha="center")
T(30,31.8,"$\\rightarrow$  predicted growth curve  $\\mu(T)$",8.2,ha="center")
tpc(66,33,15,7,pts=False,col=CORE); T(66,27.8,"a-priori TPC (genome-only)",6.9,SUB)

# ===== STAGE 4: CALIBRATION =====
arr(50,26.5,50,23.0,lw=1.5)
stage(4,21,"4","Bayesian calibration to measured curves   (per clade I–IV)",GRN)
box(6,7.5,54,11.5,GRN,GRNF,lw=1.3)
T(33,16.6,"tune 3 parameters until predicted = measured",7.8,INK,w="bold")
T(10,13.6,"•  $\\Delta T_{opt}$  (optimum shift)",7.4,ha="left")
T(10,11.5,"•  $\\Delta C_{p}$-scale  (curvature)",7.4,ha="left")
ax.text(10,9.4,"•  ",fontsize=7.4,ha="left",va="center",zorder=4)
T(11.5,9.4,"$k_{cat}$-scale  (effective growth-capacity)",7.6,RED,w="bold",ha="left")
T(45,13.2,"differential evolution\n$\\rightarrow$ emcee MCMC\n(posterior)",7.0,SUB)
# measured-data box feeding in
box(64,11,30,7,GRN,"white",lw=1.2)
T(79,15.6,"measured growth TPCs",7.4,GRN,w="bold")
T(79,12.9,"O$_2$ respirometry · 12 isolates · 22–44 °C",6.7,SUB)
arr(64,14.5,60,14.5,c=GRN,lw=1.3)
# loop back to core (re-simulate)
arr(60,18.5,88,45.5,c="#9a9a9a",lw=1.1,rad=-0.35,style="-|>")
T(93,32,"re-simulate",6.6,SUB,it=True)

# ===== OUTPUT =====
CIA={"I":"#0072B2","II":"#CC79A7","III":"#009E73","IV":"#E69F00"}
def colored_center(cx,y,segs,size):
    fig.canvas.draw(); r=fig.canvas.get_renderer(); inv=ax.transData.inverted(); ws=[]
    for (t,c,wt,it) in segs:
        tmp=ax.text(0,-50,t,fontsize=size,fontweight=wt,style="italic" if it else "normal")
        bb=tmp.get_window_extent(renderer=r)
        ws.append(inv.transform((bb.width,0))[0]-inv.transform((0,0))[0]); tmp.remove()
    x=cx-sum(ws)/2
    for (t,c,wt,it),w in zip(segs,ws):
        ax.text(x,y,t,fontsize=size,color=c,fontweight=wt,style="italic" if it else "normal",ha="left",va="center",zorder=4); x+=w
arr(33,7.5,33,4.6,lw=1.5)
box(6,-1.0,88,5.4,OUT,OUTF,lw=1.4)
T(50,2.7,"OUTPUT:  calibrated per-clade TPC ($R^2$ = 0.94–0.97)  +  effective growth-capacity ($k_{cat}$-scale):",7.9,INK,w="bold")
colored_center(50,0.35,[("C. auris ",INK,"bold",True),("II 0.39",CIA["II"],"bold",False),
    ("   <   ",SUB,"normal",False),("I 0.60",CIA["I"],"bold",False),("   ≈   ",SUB,"normal",False),
    ("III 0.62",CIA["III"],"bold",False),("   <   ",SUB,"normal",False),("IV 0.71",CIA["IV"],"bold",False)],8.4)

import os as _os
# C1: the results tree follows CANDIDAS_RESULTS, matching config.R, so a
# non-destructive re-run writes to runs/C1_reproduction/ instead of results/.
_base = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_results = _os.environ.get("CANDIDAS_RESULTS") or _os.path.join(_base, "results")
_out = _os.path.join(_results, "figures", "manuscript")
_os.makedirs(_out, exist_ok=True)
_png = _os.path.join(_out, "FIG_model_schematic.png")
plt.savefig(_png, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", _png)
