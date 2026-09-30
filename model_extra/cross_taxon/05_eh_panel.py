"""One curve parameter explains fever tolerance: growth kept at 40 C against Eh,
the deactivation energy (how steeply growth collapses past the optimum).
Medians and 95% CrI from results/tables/fig_values.csv (the manuscript's own posteriors)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, INK, MUTED, FAINT = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e8e8e6"
f = pd.read_csv("results/tables/fig_values.csv")
LAB = {"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
       "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
OFF = {"C. auris I":(8,-10),"C. auris II":(9,4),"C. auris III":(-6,10),"C. auris IV":(9,2),
       "C. parapsilosis":(9,-4),"C. haemulonii":(-8,-13),"C. duobushaemulonii":(-10,12)}
HA  = {"C. auris III":"right","C. haemulonii":"right","C. duobushaemulonii":"right"}
f["taxon"] = f.Group.map(LAB)
f["is_auris"] = f.Group.str.startswith("Clade")

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":7.5,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK,"text.color":INK})
fig, ax = plt.subplots(figsize=(4.4, 3.6))

for _, r in f.iterrows():
    c = BLUE if r.is_auris else ORANGE
    ax.plot([r.Eh_lo, r.Eh_hi], [r.growth_kept]*2, color=c, lw=0.9, alpha=0.75, zorder=2)
    ax.plot([r.Eh]*2, [r.gk_lo, r.gk_hi], color=c, lw=0.9, alpha=0.75, zorder=2)
    ax.plot(r.Eh, r.growth_kept, "o", ms=6, mfc=c, mec="white", mew=1.1, zorder=3)
    dx, dy = OFF[r.taxon]
    ax.annotate(r.taxon, (r.Eh, r.growth_kept), xytext=(dx, dy), textcoords="offset points",
                fontsize=6.4, style="italic", color=INK, ha=HA.get(r.taxon, "left"),
                va="center", zorder=4)

ax.set_xscale("log")
ax.set_xticks([2, 3, 4, 5, 6, 7]); ax.set_xticklabels(["2","3","4","5","6","7"])
ax.set_xlim(1.75, 8.6); ax.set_ylim(0, 1.0)
ax.set_xlabel("$E_h$  (eV) — steepness of the collapse past the optimum")
ax.set_ylabel("growth retained, 37 → 40 °C")
ax.grid(color=FAINT, lw=0.6, zorder=0); ax.set_axisbelow(True)
for sd in ("top","right"): ax.spines[sd].set_visible(False)
ax.plot([], [], "o", ms=5, mfc=BLUE, mec="white", label="C. auris clades")
ax.plot([], [], "o", ms=5, mfc=ORANGE, mec="white", label="relatives")
ax.legend(frameon=False, fontsize=6.6, loc="lower left", handletextpad=0.3)
ax.text(0.97, 0.96, "gentle collapse\n→ keeps growing at fever", transform=ax.transAxes,
        ha="right", va="top", fontsize=6.2, color=MUTED, linespacing=1.25)
fig.tight_layout()
fig.savefig("results/figures/manuscript/FIG_eh_explains.png", dpi=300)
fig.savefig("results/figures/manuscript/FIG_eh_explains.pdf")
print(f[["taxon","Eh","Eh_lo","Eh_hi","growth_kept","gk_lo","gk_hi"]].round(2).to_string(index=False))
