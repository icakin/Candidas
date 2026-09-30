"""SUPPLEMENTARY: draw-wise pairwise ordering probabilities among the C. auris
clades on the cost-multiplier axis, at 1.5x, 2x and 3x the draw's own minimum R/G.
The absolute crossing temperatures per taxon are already shown in the displacement
figure (FIG3_revised_displacement, panel b), so they are not repeated here."""
import os, sys, pandas as pd, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from style_common import *
LAB = {"Clade4":"IV","Clade3":"III","Clade1":"I","Clade2":"II"}
A = ["Clade4","Clade3","Clade1","Clade2"]
MULT = [("1.5×","T_1.5xcost"), ("2×","T_2xcost"), ("3×","T_3xcost")]

d = pd.read_csv(f"{C}/results/tables/thermal_headroom_draws.csv")

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
for ax, (lab, col) in zip(axes, MULT):
    W = d.pivot(index=".draw", columns="Group", values=col)
    M = np.full((4,4), np.nan)
    for i, gi in enumerate(A):
        for j, gj in enumerate(A):
            if i != j: M[i, j] = np.nanmean(W[gi] > W[gj])
    from matplotlib.colors import LinearSegmentedColormap
    ax.imshow(M, cmap=LinearSegmentedColormap.from_list("g",["#f4f6fa",GROW]), vmin=0.5, vmax=1.0)
    for i in range(4):
        for j in range(4):
            if i == j:
                ax.text(j, i, "-", ha="center", va="center", fontsize=8, color=MUTED); continue
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if M[i,j] > 0.85 else INK,
                    fontweight="bold" if M[i,j] >= 0.95 else "normal")
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels([LAB[g] for g in A], fontsize=7); ax.set_yticklabels([LAB[g] for g in A], fontsize=7)
    ax.set_xlabel("column clade (compared against)")
    ax.set_title(f"M = {lab}", fontsize=8, pad=8)
    ax.tick_params(length=0)
    for sd in ("top","right","bottom","left"): ax.spines[sd].set_visible(False)
axes[0].set_ylabel("row clade: P(reaches M × its minimum R/G\nat a higher temperature than the column clade)", fontsize=7)
fig.text(0.5, 0.965, "bold: 0.95 or more", ha="center", va="top", fontsize=6.3, color=MUTED)
fig.subplots_adjust(left=0.12, right=0.985, top=0.80, bottom=0.17, wspace=0.35)
save(fig,"SUPP_threshold_sensitivity")
