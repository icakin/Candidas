"""
c2_figures.py -- C2 PART D figures.

    cauris_etcgem/.venv/bin/python reports/tools/c2_figures.py

Writes reports/figures_N0/:
    FIG1_decoupling_<arm>.png     copied from each ARM'S OWN output tree
    FIG2_the_bill_<arm>.png       (never from results/figures/, which C1 found
    FIG_MODEL_<arm>.png            holds 7 stale duplicates)
    c2_arm_comparison.png          purpose-built: the four load-bearing
                                   quantities, all taxa, all three arms

Design notes for the comparison figure, following the dataviz method:
  * FORM. These are estimates with credible intervals, not magnitudes from a
    zero baseline, so it is a dot-and-interval plot faceted by quantity -- not
    bars. Taxa on the y axis; arms are the categorical series.
  * COLOUR. Three categorical slots from the validated default palette
    (#2a78d6 blue, #eb6834 orange, #1baf7a aqua). Validated with
    scripts/validate_palette.js --mode light: lightness band PASS, chroma floor
    PASS, CVD separation PASS (worst adjacent dE 9.2 deutan), normal-vision
    floor PASS (dE 27.6). Contrast WARNs on the aqua slot, so identity never
    rests on colour alone: each arm also has its own MARKER SHAPE and a legend,
    and every value appears numerically in reports/N0_SENSITIVITY.md (the table
    view the warning obliges).
  * One x axis per facet; no dual axes; recessive grid.
"""
import json, os, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "reports", "figures_N0")
os.makedirs(OUT, exist_ok=True)
ARMS = ["current", "ramp", "nobp"]
ARM_LABEL = {"current": "ARM 1  current  N0 = N_inoc·e^(r·δ)",
             "ramp":    "ARM 2  ramp-aware  (bound)",
             "nobp":    "ARM 3  no back-projection  (lower bound on N0)"}
COLOR = {"current": "#2a78d6", "ramp": "#eb6834", "nobp": "#1baf7a"}
MARKER = {"current": "o", "ramp": "s", "nobp": "^"}
TAXA = ["Clade1", "Clade2", "Clade3", "Clade4", "para"]
TAXLAB = {"Clade1": "C. auris I", "Clade2": "C. auris II", "Clade3": "C. auris III",
          "Clade4": "C. auris IV", "para": "C. parapsilosis"}
INK, SUB, GRID = "#0b0b0b", "#52514e", "#e6e5e1"

# ---------------------------------------------------------------------------
# 1. collect each arm's own figures
# ---------------------------------------------------------------------------
copied = []
for arm in ARMS:
    src = os.path.join(ROOT, "runs", f"C2_arm_{arm}", "figures", "manuscript")
    if not os.path.isdir(src):
        print(f"  no figures for arm {arm} at {src}")
        continue
    for name in ("FIG1_decoupling", "FIG2_the_bill", "FIG_MODEL"):
        p = os.path.join(src, name + ".png")
        if os.path.exists(p):
            q = os.path.join(OUT, f"{name}_{arm}.png")
            shutil.copyfile(p, q); copied.append(os.path.basename(q))
print(f"collected {len(copied)} per-arm figures into {os.path.relpath(OUT, ROOT)}")

# ---------------------------------------------------------------------------
# 2. the comparison figure
# ---------------------------------------------------------------------------
fv = {}
for arm in ARMS:
    p = os.path.join(ROOT, "runs", f"C2_arm_{arm}", "tables", "fig_values.csv")
    if os.path.exists(p):
        fv[arm] = pd.read_csv(p).set_index("Group")
have = [a for a in ARMS if a in fv]
print(f"arms with fig_values.csv: {have}")
if not have:
    raise SystemExit("no arm results yet")

PANELS = [("E_R  (respiration activation energy, eV)", "ER", "ER_lo", "ER_hi", None),
          ("ΔE = E_G − E_R  (eV)", "dE", "dE_lo", "dE_hi", 0.0),
          ("T_opt(CUE)  (°C)", "Tcue", "Tcue_lo", "Tcue_hi", 37.0),
          ("Fever tax at 40 °C  (×)", "tax40", "tax40_lo", "tax40_hi", 1.0)]

fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.6))
off = {a: (i - (len(have) - 1) / 2) * 0.26 for i, a in enumerate(have)}

for ax, (title, c, clo, chi, ref) in zip(axes, PANELS):
    for a in have:
        ys, xs, los, his = [], [], [], []
        for j, g in enumerate(TAXA):
            if g not in fv[a].index:
                continue
            ys.append(len(TAXA) - 1 - j + off[a])
            xs.append(float(fv[a].loc[g, c]))
            los.append(float(fv[a].loc[g, clo])); his.append(float(fv[a].loc[g, chi]))
        xs = np.array(xs); los = np.array(los); his = np.array(his); ys = np.array(ys)
        ax.hlines(ys, los, his, color=COLOR[a], lw=2.0, alpha=.9,
                  zorder=2, capstyle="round")
        ax.scatter(xs, ys, s=42, marker=MARKER[a], facecolor=COLOR[a],
                   edgecolor="white", linewidth=1.1, zorder=3,
                   label=ARM_LABEL[a] if ax is axes[0] else None)
    if ref is not None:
        ax.axvline(ref, color="#B22222", ls=(0, (2, 2)), lw=1.0, zorder=1)
        ax.annotate(f"{ref:g}", xy=(ref, len(TAXA) - 0.55), color="#B22222",
                    fontsize=7.5, ha="center", va="bottom",
                    bbox=dict(fc="white", ec="none", pad=1.0))
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=8)
    ax.set_yticks(range(len(TAXA)))
    ax.set_yticklabels([TAXLAB[g] for g in reversed(TAXA)], fontsize=8.5, color=INK)
    ax.set_ylim(-0.6, len(TAXA) - 0.4)
    ax.tick_params(axis="x", labelsize=8, colors=SUB)
    ax.grid(axis="x", color=GRID, lw=.7)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)

for ax in axes[1:]:
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)

h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8.5,
           bbox_to_anchor=(0.5, -0.015))
fig.suptitle("N0 sensitivity: the three treatments, with 95% credible intervals",
             fontsize=12, color=INK, x=0.008, ha="left", y=0.99)
fig.text(0.008, 0.925,
         "Growth-only quantities are bit-identical across arms; only respiration moves. "
         "Every value is tabulated in reports/N0_SENSITIVITY.md.",
         fontsize=8.5, color=SUB, ha="left")
fig.tight_layout(rect=[0, 0.06, 1, 0.90])
p = os.path.join(OUT, "c2_arm_comparison.png")
fig.savefig(p, dpi=220, facecolor="white")
print("wrote", p)

# ---------------------------------------------------------------------------
# 3. the PART B diagnostic panel -- transient + measured ramp
# ---------------------------------------------------------------------------
pb = pd.read_csv(os.path.join(ROOT, "reports", "tools", "c2_partB_per_series.csv"))
ti = pd.read_csv(os.path.join(ROOT, "reports", "tools", "c2_Tinternal.csv"))
F = pb[pb.status == "ok"]

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

ax = axes[0]
byT = F.groupby("T").A
ax.scatter(F["T"] + np.random.default_rng(0).uniform(-.25, .25, len(F)), F.A,
           s=6, color=COLOR["current"], alpha=.25, edgecolor="none")
ax.plot(sorted(byT.median().index), byT.median().values, color=COLOR["current"], lw=2)
ax.set_title("Transient amplitude A scales with the temperature step",
             fontsize=9.5, color=INK, loc="left")
ax.set_xlabel("Set point (°C)", fontsize=8.5, color=SUB)
ax.set_ylabel("A (mg/L)", fontsize=8.5, color=SUB)
ax.annotate("r(A, |Δsolubility|) = +0.996", xy=(23, 5.4), fontsize=8, color=INK)

ax = axes[1]
byT = F.groupby("T").tau
ax.scatter(F["T"] + np.random.default_rng(1).uniform(-.25, .25, len(F)), F.tau,
           s=6, color=COLOR["ramp"], alpha=.25, edgecolor="none")
ax.plot(sorted(byT.median().index), byT.median().values, color=COLOR["ramp"], lw=2)
ax.set_ylim(0, 60)
ax.set_title("…but the time constant τ does not", fontsize=9.5, color=INK, loc="left")
ax.set_xlabel("Set point (°C)", fontsize=8.5, color=SUB)
ax.set_ylabel("τ (min)", fontsize=8.5, color=SUB)
ax.annotate("17–20 min at every temperature", xy=(23, 50), fontsize=8, color=INK)

ax = axes[2]
_canon = {g.lower(): g for g in ["Clade1", "Clade2", "Clade3", "Clade4", "glab", "para"]}
ti["group"] = ti.group.str.strip().str.lower().map(_canon)
for Ts, col in ((22, "#8ab4e8"), (34, "#eb6834"), (44, "#7a1f00")):
    g = ti[ti.T_set == Ts]
    grid = np.arange(0, 200, 1.0)
    Y = np.mean([np.interp(grid, gg.Time.values, gg.T_internal.values)
                 for _, gg in g.groupby("File")], axis=0)
    ax.plot(grid, Y - Y[-1], color=col, lw=2, label=f"set point {Ts} °C")
ax.axhline(0, color=GRID, lw=1)
ax.axvline(75.2, color="#B22222", ls=(0, (2, 2)), lw=1.0)
ax.annotate("median δ\n75 min", xy=(78, -2.6), fontsize=7.5, color="#B22222")
ax.set_title("Measured thermal ramp (T_internal), relative to settled",
             fontsize=9.5, color=INK, loc="left")
ax.set_xlabel("Time (min)", fontsize=8.5, color=SUB)
ax.set_ylabel("T − T settled (K)", fontsize=8.5, color=SUB)
ax.legend(frameon=False, fontsize=8, loc="lower right")

for ax in axes:
    ax.tick_params(labelsize=8, colors=SUB)
    ax.grid(color=GRID, lw=.7); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(GRID)

fig.suptitle("The pre-window transient is instrumental, and the thermal ramp is measurable",
             fontsize=12, color=INK, x=0.006, ha="left", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.92])
p = os.path.join(OUT, "c2_transient_diagnostics.png")
fig.savefig(p, dpi=220, facecolor="white")
print("wrote", p)
