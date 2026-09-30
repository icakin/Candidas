# =============================================================================
# 20_fig4.py - FIGURE 4: fever switches the relatives into a respiring, non-growing state
#
# Model-free. Every one of the 1,380 oxygen series is classified from the raw trace and
# the pipeline's own fit tables; nothing is fitted here that is not already in the
# pipeline, except a straight line through the non-growing wells.
#
# STATES (per well)
#   growing   : the pipeline accepted a growth fit (07_oxygen_fits.R: curvature by AICc,
#               valid fit) - exactly the wells that feed Figs 1-2. At 40 C this rule and
#               the stricter Fig 3 rule (rt >= 1.05, drawdown >= 2 mg/L, QC) give
#               identical counts in every taxon.
#   respiring : no accepted growth fit (growth below detection, net rate < ~0.1 /h), but the well drew down >= MIN_O2_DRAWDOWN mg/L
#               of O2 from the pipeline's window start (after thermal equilibration,
#               04_trim_selector.R) to the end of the trace - the Fig 3 signal test.
#   inert     : drew down less than that.
#
# VOLUMETRIC RESPIRATION (panel d) compares like with like: for growing wells the
# pipeline's K (O2 consumption rate at the window start, mg/L/min -> /h); for
# non-growing wells the slope of a straight line through the first EARLY_H hours of the
# same window. Both start from the same inoculum.
#
# PER-CELL numbers for non-growing wells take N = N_inoc. Undetected growth (rt < ~1,
# i.e. < ~1.5 doublings over the window) could raise N by at most ~2x, so those numbers
# are UPPER bounds and are reported as such.
#
#
# Needs: pandas, numpy, matplotlib.   Run:  python3 scripts/20_fig4.py
# OUTPUTS results/figures/manuscript/FIG4_fever_state.png/.pdf
#         results/tables/fig4_well_states.csv, fig4_state_counts.csv, fig4_summary.csv
# =============================================================================
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

root = "." if os.path.exists("results/tables/fit_metrics.csv") else ".."
TAB = f"{root}/results/tables"; FIG = f"{root}/results/figures/manuscript"
os.makedirs(FIG, exist_ok=True)

MIN_O2_DRAWDOWN = 2.0      # mg/L, matches 04_trim_selector.R and Fig 3
EARLY_H         = 5.0      # h of the trimmed window used for the non-growing slope
FGC_PER_MG_O2   = 3.7536721e11   # pipeline conversion (12/32 g C per g O2, RQ = 1)
EXCLUDE         = set()
FEVER           = (37, 40)

# ---- 1. tables ---------------------------------------------------------------
otu  = pd.read_csv(f"{TAB}/otu_names.csv")
inoc = pd.read_csv(f"{TAB}/otu_inoc.csv")[["OTU", "N_inoc_cells_per_L"]]
der  = pd.read_csv(f"{TAB}/derived_N0_R_results_with_carbon.csv")[
        ["T", "OTU", "Replicate", "r", "K", "fit_start_time", "fit_end_time",
         "growth_fgC_h", "respiration_fgC_h", "CUE"]]
meta = pd.read_csv(f"{TAB}/Oxygen_Trimmed_Series_Metadata.csv").set_index(["T", "OTU", "Replicate"])
raw  = pd.read_csv(f"{TAB}/Oxygen_All_Long.csv")

rows = []
for (T, o, rep), s in raw.groupby(["T", "OTU", "Replicate"]):
    s = s.sort_values("Time"); t = s.Time.values / 60; y = s.Oxygen.values
    try:
        t0 = meta.loc[(T, o, rep), "main_run_start_time"] / 60
        t1 = meta.loc[(T, o, rep), "chosen_end_time"] / 60
    except KeyError:
        t0, t1 = 3.0, t[-1]
    sel = (t >= t0)                      # from the pipeline's window start to the END of the trace
    if sel.sum() < 20: sel = t >= 3.0
    ts, ys = t[sel], y[sel]
    e = ts <= ts[0] + EARLY_H
    rows.append(dict(T=T, OTU=o, Replicate=rep,
                     drawdown=ys[:5].mean() - ys[-5:].mean(),   # over the whole trace after the window start (Fig 3 rule)
                     slope_early=-np.polyfit(ts[e], ys[e], 1)[0],
                     window_h=ts[-1] - ts[0]))
x = (pd.DataFrame(rows)
     .merge(der, how="left", on=["T", "OTU", "Replicate"])
     .merge(otu[["OTU", "otu_name", "group", "species"]], on="OTU")
     .merge(inoc, on="OTU"))
x = x[~x.group.isin(EXCLUDE)].copy()
x["sp"] = np.where(x.group.str.startswith("Clade"), "auris",
          np.where(x.group == "Hae", "haemulonii",
          np.where(x.group == "Duo", "duobushaemulonii", "parapsilosis")))
x["state"] = np.where(x.respiration_fgC_h.notna(), "growing",
             np.where(x.drawdown >= MIN_O2_DRAWDOWN, "respiring", "inert"))
x["V_h"]   = np.where(x.state == "growing", x.K * 60, x.slope_early)        # mg O2 /L/h
x["R_ub_fgC_h"] = np.where(x.state == "growing", x.respiration_fgC_h,
                           x.slope_early / x.N_inoc_cells_per_L * FGC_PER_MG_O2)  # upper bound if not growing
x["CUE_well"] = np.where(x.state == "growing", x.CUE, np.where(x.state == "respiring", 0.0, np.nan))
x.to_csv(f"{TAB}/fig4_well_states.csv", index=False)

counts = (x.groupby(["group", "T", "state"]).size().unstack("state")
            .reindex(columns=["growing", "respiring", "inert"]).fillna(0).astype(int))
counts.to_csv(f"{TAB}/fig4_state_counts.csv")

# ---- 2. numbers for the caption ---------------------------------------------
summ = []
for sp, s in x.groupby("sp"):
    for T in (38, 40, 42, 44):
        u = s[s["T"] == T]
        if len(u) == 0: continue
        ng = u[u.state != "growing"]
        summ.append(dict(species=sp, T=T, n=len(u),
                         growing=(u.state == "growing").sum(), respiring=(u.state == "respiring").sum(),
                         inert=(u.state == "inert").sum(),
                         V_growing_med=u[u.state == "growing"].V_h.median(),
                         V_nongrowing_med=ng.V_h.median(),
                         Rcell_ub_nongrowing_med=ng.R_ub_fgC_h.median(),
                         Rcell_growing_med=u[u.state == "growing"].respiration_fgC_h.median(),
                         O2_consumed_nongrowing_med=ng.drawdown.median(),
                         CUE_med=u.CUE_well.median()))
summ = pd.DataFrame(summ); summ.to_csv(f"{TAB}/fig4_summary.csv", index=False)
iso = (x[x["T"] >= 38].groupby(["sp", "T", "otu_name"]).apply(lambda u: (u.state == "growing").sum() >= 3)
         .groupby(level=[0, 1]).agg(["sum", "size"]).rename(columns={"sum": "isolates_growing", "size": "isolates"}))
iso.to_csv(f"{TAB}/fig4_isolate_summary.csv"); print(iso.to_string())
pd.set_option("display.width", 250); print(summ.round(2).to_string())

# ---- 3. figure ----------------------------------------------------------------
INK, MUTED, FAINT = "#1a1a1a", "#4d4d4d", "#dedede"
NAVY, TEAL, RED, GREY = "#1d3f73", "#2f9aa6", "#a5453b", "#b8b8b8"
COL = {"growing": NAVY, "respiring": "#d9822b", "inert": GREY}
LAB = {"Clade1": "C. auris I", "Clade2": "C. auris II", "Clade3": "C. auris III", "Clade4": "C. auris IV",
       "Hae": "C. haemulonii", "Duo": "C. duobushaemulonii", "para": "C. parapsilosis"}
SPL = {"auris": "C. auris (12 isolates)", "haemulonii": "C. haemulonii (3)",
       "duobushaemulonii": "C. duobushaemulonii (2)", "parapsilosis": "C. parapsilosis (3)"}
ORDER = ["Clade1", "Clade2", "Clade3", "Clade4", "Hae", "Duo", "para"]
SPO   = ["auris", "haemulonii", "duobushaemulonii", "parapsilosis"]
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK})

def clean(ax):
    for sd in ("top", "right"): ax.spines[sd].set_visible(False)
def fever(ax): ax.axvspan(*FEVER, color="#000000", alpha=0.06, lw=0, zorder=0)
def letter(ax, L, title, y=1.06, italic=False):
    ax.text(-0.02, y, L, transform=ax.transAxes, fontweight="bold", fontsize=10, ha="right", va="bottom")
    ax.text(0.0, y, title, transform=ax.transAxes, fontsize=8, ha="left", va="bottom", color=INK,
            style="italic" if italic else "normal")

fig = plt.figure(figsize=(7.2, 4.0))
gs = gridspec.GridSpec(2, 4, figure=fig, width_ratios=[1.9, 1, 1, 0.0001], hspace=0.55, wspace=0.28,
                       left=0.08, right=0.99, top=0.9, bottom=0.13)

# (a) raw traces: one C. parapsilosis isolate at 38-44 C, one C. auris Clade IV well at 40 C as reference
TCOL = {38: "#7fb3c9", 40: "#d9822b", 42: "#a5453b", 44: "#5a1f1a"}
axA = fig.add_subplot(gs[:, 0])
axA.axvspan(0, 1.5, color="#000000", alpha=0.04, lw=0)
iso_auris = x[(x.group == "Clade4")].OTU.min(); iso_para = x[(x.group == "para")].OTU.min()
ref = raw[(raw.OTU == iso_auris) & (raw["T"] == 40)]
w = ref[ref.Replicate == ref.Replicate.iloc[0]].sort_values("Time")
axA.plot(w.Time / 60, w.Oxygen, color=GREY, lw=1.3, label="C. auris IV, 40 °C")
for T in (38, 40, 42, 44):
    s_ = raw[(raw.OTU == iso_para) & (raw["T"] == T)]
    for rep, w in s_.groupby("Replicate"):
        w = w.sort_values("Time"); axA.plot(w.Time / 60, w.Oxygen, color=TCOL[T], lw=0.6, alpha=0.9)
    axA.plot([], [], color=TCOL[T], lw=1.2, label=f"{T} °C")
axA.set_xlim(0, 21); axA.set_ylim(0, 11); clean(axA)
axA.set_xlabel("time (h)"); axA.set_ylabel("dissolved O$_2$ (mg L$^{-1}$)")
axA.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.4)
letter(axA, "a", "C. parapsilosis, one isolate", y=1.03, italic=True)

# (b) O2 consumption across the growth limit, per species (2 x 2)
rng = np.random.default_rng(1)
pos = {"auris": (0, 1), "haemulonii": (0, 2), "duobushaemulonii": (1, 1), "parapsilosis": (1, 2)}
axes = {}
for sp in SPO:
    r_, c_ = pos[sp]; ax = fig.add_subplot(gs[r_, c_]); axes[sp] = ax
    s_ = x[x.sp == sp]
    for st, m in (("growing", "o"), ("respiring", "o"), ("inert", "x")):
        u = s_[s_.state == st]
        ax.scatter(u["T"] + rng.uniform(-0.45, 0.45, len(u)), u.V_h.clip(lower=0.02), s=4, marker=m,
                   color=COL[st], alpha=0.5 if st == "growing" else 0.9, lw=0.5, zorder=2)
    med = s_[s_.state != "inert"].groupby("T").V_h.median()
    ax.plot(med.index, med.values, color=INK, lw=0.9, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(0.03, 3); ax.set_xlim(21, 45); ax.set_xticks([22, 30, 38, 44])
    fever(ax); clean(ax); ax.set_title(SPL[sp], style="italic", fontsize=7, pad=2)
    if r_ == 1: ax.set_xlabel("temperature (°C)")
    else: plt.setp(ax.get_xticklabels(), visible=False)
    if c_ == 1: ax.set_ylabel("O$_2$ consumption\n(mg L$^{-1}$ h$^{-1}$)", fontsize=7)
    else: plt.setp(ax.get_yticklabels(), visible=False)
letter(axes["auris"], "b", "Every well: blue growing, orange respiring with growth below detection, grey x inert", y=1.22)

fig.savefig(f"{FIG}/FIG4_fever_state.png", dpi=300); fig.savefig(f"{FIG}/FIG4_fever_state.pdf")
print("written", f"{FIG}/FIG4_fever_state.png")
