"""Thermal position figure: growth optimum of every organism measured on this platform,
against the temperature of the environment it lives in."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# slot 1 blue, 2 orange, 3 aqua, 7 violet — validated all-pairs, light mode
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
INK, MUTED, FAINT = "#0b0b0b", "#52514e", "#e8e8e6"

D = [
 # name, Topt, lo, hi, group, host (None = no host), undetermined
 ("Zymoseptoria tritici", 23.2, 21.0, 26.4, "plant pathogen", 18.5, False),
 ("OTU13", 27.9, 26.6, 32.7, "environmental bacteria", None, False),
 ("OTU4",  29.2, 29.0, 29.5, "environmental bacteria", None, False),
 ("OTU9",  29.5, 28.6, 30.3, "environmental bacteria", None, False),
 ("OTU14", 30.5, 29.6, 33.5, "environmental bacteria", None, False),
 ("OTU7",  32.3, 31.9, 32.7, "environmental bacteria", None, False),
 ("OTU8",  32.3, 32.0, 32.7, "environmental bacteria", None, False),
 ("OTU20", 32.9, 32.6, 33.3, "environmental bacteria", None, False),
 ("OTU6",  35.9, 34.2, 41.2, "environmental bacteria", None, True),
 ("C. haemulonii",       31.9, 30.2, 34.0, "Candida", 37, False),
 ("C. parapsilosis",     33.4, 31.6, 35.3, "Candida", 37, False),
 ("C. duobushaemulonii", 33.5, 31.5, 35.9, "Candida", 37, False),
 ("C. auris II",         34.8, 33.0, 36.7, "Candida", 37, False),
 ("C. auris I",          35.5, 33.6, 37.4, "Candida", 37, False),
 ("C. auris III",        35.7, 33.9, 37.7, "Candida", 37, False),
 ("C. auris IV",         36.4, 34.5, 38.4, "Candida", 37, False),
 ("E. coli K-12 (B)",    39.4, 38.8, 40.2, "E. coli", 37, False),
 ("E. coli K-12 (A)",    40.4, 39.8, 40.9, "E. coli", 37, False),
]
COL = {"plant pathogen": VIOLET, "environmental bacteria": AQUA, "Candida": BLUE, "E. coli": ORANGE}
ORDER = ["plant pathogen", "environmental bacteria", "Candida", "E. coli"]

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":7,"axes.linewidth":0.6,"xtick.major.width":0.6,"axes.edgecolor":MUTED,
    "xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
fig, ax = plt.subplots(figsize=(7.2, 5.0))

rows, y, ticks, labs, seps, heads = [], 0, [], [], [], []
for g in ORDER:
    heads.append((y + 0.85, g))
    for nm, t, lo, hi, gg, host, und in D:
        if gg != g: continue
        rows.append((y, nm, t, lo, hi, g, host, und)); ticks.append(y); labs.append(nm); y -= 1
    seps.append(y + 0.5); y -= 1.35

ax.axvspan(36.6, 37.4, color=INK, alpha=0.07, lw=0, zorder=0)
ax.text(37, 1.05, "human body\n37 °C", ha="center", va="bottom", fontsize=7, color=MUTED, linespacing=1.15)
ax.axvspan(15, 22, color=VIOLET, alpha=0.07, lw=0, zorder=0)
ax.text(18.5, 1.05, "wheat canopy\n15–22 °C", ha="center", va="bottom", fontsize=7, color=MUTED, linespacing=1.15)

for yy, nm, t, lo, hi, g, host, und in rows:
    c = COL[g]
    if host is not None:                       # distance from optimum to the host it lives in
        ax.plot([t, host], [yy, yy], color=c, lw=2.0, alpha=0.35, solid_capstyle="round", zorder=2)
    ax.plot([lo, hi], [yy, yy], color=c, lw=1.1, solid_capstyle="round", zorder=3,
            ls=(0,(2.2,1.6)) if und else "-")
    if not und:
        for e in (lo, hi):
            ax.plot([e, e], [yy-0.17, yy+0.17], color=c, lw=1.1, zorder=3)
    ax.plot([t], [yy], "o", ms=6.5, mfc="white" if und else c, mec=c, mew=1.4, zorder=4)
    if und:
        ax.text(hi+0.6, yy, "optimum not determined\n(still rising at 35 °C)", fontsize=5.8,
                color=MUTED, va="center", ha="left", linespacing=1.1)


ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=7)
for lab, (_, nm, *_ ) in zip(ax.get_yticklabels(), rows):
    if nm.startswith(("C.","Zymo","E. coli")): lab.set_style("italic")
ax.set_xlim(14, 50); ax.set_ylim(y + 0.8, 1.9)
ax.set_xticks([15, 20, 25, 30, 35, 40, 45])
ax.set_xlabel("growth optimum (°C), 95% CI")
for sd in ("top","right","left"): ax.spines[sd].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color=FAINT, lw=0.6, zorder=0); ax.set_axisbelow(True)

for yy, g in heads:
    ax.text(0.0, yy, g.upper(), transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=6.8, color=COL[g], fontweight="bold")

ax.legend(handles=[Line2D([],[],color=MUTED,lw=2.4,alpha=0.35,
          label="gap between growth optimum and the temperature of the host")],
          frameon=False, fontsize=6.8, loc="upper left", bbox_to_anchor=(0.0, -0.085), handlelength=2.2)
fig.subplots_adjust(left=0.22, right=0.985, top=0.885, bottom=0.175)
fig.savefig("results/figures/manuscript/FIG_thermal_position.png", dpi=300)
fig.savefig("results/figures/manuscript/FIG_thermal_position.pdf")
print("written results/figures/manuscript/FIG_thermal_position.png")
