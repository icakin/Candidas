"""Carbon-conversion tipping point (Robustness): how strongly the per-cell carbon quota q or the
oxygen-to-carbon conversion c (respiratory quotient) would have to depend on temperature for
the apparent-CUE optimum to reach 37 C, or for the 37->40 C rise in respiration per unit growth
to vanish. A log-linear factor alpha(T) = exp(beta (T - 37)) multiplies R/G (beta > 0 means
q/c falls, or c/q rises, with temperature). The optimum moves to 37 C when the cost curve is flat
there, i.e. beta* = -d ln(R/G)/dT at 37 C; the fever cost vanishes when beta = -[ln cost(40) -
ln cost(37)]/3. Both are computed per posterior draw and group. Writes
results/tables/conversion_tipping.csv and the supplementary figure FIG_SUPP_conversion_tipping."""
import numpy as np, pandas as pd
from tpc_common import *
from matplotlib.lines import Line2D
G,I,quota,vals=load()
rows=[]; draws={}
for g in GROUPS:
    d=G[g]; E,ER,Eh,Th=d.E.values,d.ER.values,d.Eh.values,d.Th.values
    h=0.05
    slope37=(ln_cost(37+h,E,ER,Eh,Th)-ln_cost(37-h,E,ER,Eh,Th))/(2*h)      # d ln(R/G)/dT at 37
    beta_opt=-slope37                                                      # per-degree factor that flattens the cost at 37
    beta_cost=-(ln_cost(40,E,ER,Eh,Th)-ln_cost(37,E,ER,Eh,Th))/3          # per-degree factor that erases the 37->40 rise
    # sign: beta<0 means q/c must RISE with temperature (quota up or RQ-conversion down) by |beta| per degree
    draws[g]=(beta_opt,beta_cost)
    for name,b in (("optimum_to_37",beta_opt),("fever_cost_erased",beta_cost)):
        if g=="Duo" and name=="fever_cost_erased": continue   # no measurable growth at 40 C
        pc=100*(np.exp(-b)-1)  # % change in q/c per degree needed (positive = q/c must increase with T)
        rows.append(dict(group=g,target=name,beta_median=np.median(b),beta_lo=np.percentile(b,2.5),beta_hi=np.percentile(b,97.5),
                         pct_per_degree_median=np.median(pc),pct_per_degree_lo=np.percentile(pc,2.5),pct_per_degree_hi=np.percentile(pc,97.5),
                         pct_over_37_to_40_median=np.median(100*(np.exp(-3*b)-1))))
R=pd.DataFrame(rows); R.to_csv(f"{C}/results/tables/conversion_tipping.csv",index=False)
pd.set_option("display.width",200); print(R.round(2).to_string(index=False))
# ---- figure: two panels, full range and 0 to 15% per degree ----
fig,(ax,ins)=plt.subplots(1,2,figsize=(7.2,3.6),gridspec_kw=dict(width_ratios=[1.35,1],wspace=0.08))
def draw(a,ms,lw):
    for i,g in enumerate(GROUPS):
        y=len(GROUPS)-1-i
        for name,mk,off in (("optimum_to_37","o",0.15),("fever_cost_erased","s",-0.15)):
            rr=R[(R.group==g)&(R.target==name)]
            if len(rr)==0: continue
            r=rr.iloc[0]
            a.plot([r.pct_per_degree_lo,r.pct_per_degree_hi],[y+off,y+off],color=OI[g],lw=lw,alpha=0.9)
            a.scatter(r.pct_per_degree_median,y+off,marker=mk,s=ms,color=OI[g],edgecolor="white",linewidth=0.6,zorder=3)
    a.axvline(0,color=AMB,lw=0.8); a.axvspan(-1,1,color=INK,alpha=0.08,lw=0); a.set_ylim(-0.7,len(GROUPS)-0.2)
draw(ax,26,1.3); draw(ins,26,1.3)
ax.set_xlim(-2,66); ins.set_xlim(-1.5,15); ins.set_xticks([0,5,10,15])
ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels([SHORT[g] for g in GROUPS][::-1],style="italic",fontsize=6.6); ax.tick_params(axis="y",length=0); ax.spines["left"].set_visible(False)
ins.set_yticks([]); ins.spines["left"].set_visible(False)
ax.set_title("full range",fontsize=7,loc="left",pad=3); ins.set_title("0 to 15% per °C",fontsize=7,loc="left",pad=3)
fig.supxlabel("required rise in q/c with temperature (% per °C); shaded, ±1% per °C illustrative reference",fontsize=7.2,y=0.0)
ax.legend(handles=[Line2D([],[],marker="o",ls="",color=INK,ms=5,label="CUE optimum to 37 °C"),Line2D([],[],marker="s",ls="",color=INK,ms=5,label="37 to 40 °C cost erased")],loc="lower right",fontsize=6)
panel_label(ax,"a",x=-0.3); panel_label(ins,"b",x=-0.04)
save(fig,"FIG_SUPP_conversion_tipping")
