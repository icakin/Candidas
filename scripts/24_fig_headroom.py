"""FIGURE: thermal headroom inside the host.
(a) the temperature at which respiration per unit growth reaches twice its own minimum,
    computed on every posterior draw with that draw's own minimum, so it is scale-free.
(b) how far that temperature sits above body temperature and above fever.
Input: results/tables/thermal_headroom_draws.csv (scripts/23_thermal_headroom.R)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE, ORANGE, INK, MUTED, FAINT = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e8e8e6"
LAB = {"Clade4":"C. auris IV","Clade3":"C. auris III","Clade1":"C. auris I","Clade2":"C. auris II",
       "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
ORDER = ["Clade4","Clade3","Clade1","Clade2","para","Hae","Duo"]

d = pd.read_csv("results/tables/thermal_headroom_draws.csv")
q = lambda x, p: np.nanpercentile(x, p)
S = {}
for g, s in d.groupby("Group"):
    t2 = s.T_2xcost.dropna(); hb = s.headroom_body.dropna(); hf = s.headroom_fever.dropna()
    S[g] = dict(t2=t2.median(), t2lo=q(t2,2.5), t2hi=q(t2,97.5),
                hb=hb.median(), hblo=q(hb,2.5), hbhi=q(hb,97.5), Pb=float((hb>0).mean()),
                hf=hf.median(), hflo=q(hf,2.5), hfhi=q(hf,97.5), Pf=float((hf>0).mean()))

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.05]))
Y = {g: -i for i, g in enumerate(ORDER)}
col = lambda g: BLUE if g.startswith("Clade") else ORANGE

# ---- (a) cost-doubling temperature -----------------------------------------
axA.axvspan(36.8, 37.2, color=INK, alpha=0.10, lw=0, zorder=0)
axA.axvspan(39.8, 40.2, color=INK, alpha=0.06, lw=0, zorder=0)
axA.axvline(44, color=MUTED, lw=0.8, ls=(0,(3,2)), zorder=1)
axA.text(37, 1.0, "body", ha="center", va="bottom", fontsize=6.4, color=MUTED)
axA.text(40, 1.0, "fever", ha="center", va="bottom", fontsize=6.4, color=MUTED)
axA.text(44.25, -6.4, "assay\nceiling", ha="left", va="bottom", fontsize=6.0, color=MUTED, linespacing=1.15)
for g in ORDER:
    v = S[g]; y = Y[g]; c = col(g)
    axA.plot([v["t2lo"], v["t2hi"]], [y, y], color=c, lw=1.1, solid_capstyle="round", zorder=2)
    for e in (v["t2lo"], v["t2hi"]): axA.plot([e, e], [y-.16, y+.16], color=c, lw=1.1, zorder=2)
    axA.plot(v["t2"], y, "o", ms=6, mfc=c, mec="white", mew=1.1, zorder=3)
axA.set_yticks(list(Y.values())); axA.set_yticklabels([LAB[g] for g in ORDER], fontsize=7, style="italic")
axA.set_xlim(33.5, 47.6); axA.set_ylim(-6.6, 1.5); axA.set_xticks([34,37,40,43,46])
axA.set_xlabel("temperature at which the carbon cost doubles (°C)")
axA.grid(axis="x", color=FAINT, lw=0.6, zorder=0); axA.set_axisbelow(True)
for sd in ("top","right","left"): axA.spines[sd].set_visible(False)
axA.tick_params(axis="y", length=0)


# ---- (b) headroom ------------------------------------------------------------
axB.axvline(0, color=INK, lw=0.9, zorder=1)
axB.text(0, 1.0, "no headroom left", ha="center", va="bottom", fontsize=6.4, color=INK)
for g in ORDER:
    v = S[g]; y = Y[g]; c = col(g)
    axB.plot([v["hf"], v["hb"]], [y, y], color=c, lw=1.6, alpha=0.30, solid_capstyle="round", zorder=2)
    for key, lo, hi, filled, dy in (("hb","hblo","hbhi",True, .12), ("hf","hflo","hfhi",False,-.12)):
        axB.plot([v[lo], v[hi]], [y+dy]*2, color=c, lw=0.9, alpha=0.85, zorder=3)
        axB.plot(v[key], y+dy, "o", ms=5.4, mfc=c if filled else "white", mec=c, mew=1.2, zorder=4)
    axB.text(9.9, y, f"{v['Pb']:.2f} / {v['Pf']:.2f}", fontsize=6.2, color=MUTED, va="center", ha="right")
axB.text(9.9, 0.62, "P(headroom > 0)\nbody / fever", fontsize=6.0, color=MUTED, va="bottom", ha="right", linespacing=1.2)
axB.set_yticks(list(Y.values())); axB.set_yticklabels([])
axB.set_xlim(-6.6, 10.0); axB.set_ylim(-6.6, 1.5); axB.set_xticks([-6,-4,-2,0,2,4,6,8])
axB.set_xlabel("headroom: °C between the host and cost doubling")
axB.grid(axis="x", color=FAINT, lw=0.6, zorder=0); axB.set_axisbelow(True)
for sd in ("top","right","left"): axB.spines[sd].set_visible(False)
axB.tick_params(axis="y", length=0)
axB.legend(handles=[Line2D([],[],marker="o",ls="",mfc=MUTED,mec=MUTED,ms=5.4,label="at 37 °C"),
                    Line2D([],[],marker="o",ls="",mfc="white",mec=MUTED,mew=1.2,ms=5.4,label="at 40 °C")],
           frameon=False, fontsize=6.6, loc="lower left", handletextpad=0.3, bbox_to_anchor=(-0.01,-0.02))


fig.subplots_adjust(left=0.175, right=0.985, top=0.86, bottom=0.155, wspace=0.10)
fig.text(0.012, 0.945, "a", fontweight="bold", fontsize=10, va="top")
fig.text(0.545, 0.945, "b", fontweight="bold", fontsize=10, va="top")
fig.savefig("results/figures/manuscript/FIG_headroom.png", dpi=300)
fig.savefig("results/figures/manuscript/FIG_headroom.pdf")
print("written results/figures/manuscript/FIG_headroom.png")
