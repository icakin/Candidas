"""
c2_partB5_ramp_measured.py -- C2 PART B2, rebuilt on the MEASURED thermal ramp.

    cauris_etcgem/.venv/bin/python reports/tools/c2_partB5_ramp_measured.py

Supersedes the modelled ramp in c2_partB_transient.py (B5).

WHY. That version had to assume the sample started at ambient (20 C) and warmed
with the time constant of the OXYGEN signal (~19 min), because nothing better
was thought to be available. It is: every raw PreSens export carries a
`T_internal [degC]` column - a real measured temperature trace - which
01_convert_xlsx.R discards. reports/tools/c2_extract_Tinternal.R recovers it.

The measured trace tells a different story from the assumption:

  * the plate does NOT start at ambient. At a 44 C set point T_internal begins
    at 41.97 C, only ~1.6 C below its own asymptote;
  * it then DIPS (39.32 C at t = 20 min) as the load equilibrates, and recovers
    over roughly an hour: 41.45 at 40 min, 43.01 at 75 min, 43.39 at 120 min;
  * at a 22 C set point there is essentially no ramp at all (22.24 -> 22.67).

So the thermal excursion is a dip-and-recovery on a ~1 h timescale, not a
warm-up from ambient on a ~19 min one, and the ramp error is correspondingly
SMALLER than the modelled version implied. This script reports the measured
version and writes the ARM 2 factors from it.

CAVEAT, stated in the report: T_internal is the reader's internal sensor, not a
probe in the liquid, so it is still a proxy for vial temperature. It is a
MEASURED proxy, which the previous assumption was not.
"""
import json, os, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS = os.path.join(ROOT, "reports", "tools")
GROUPS = ["Clade1", "Clade2", "Clade3", "Clade4", "glab", "para"]
R = {}


def head(t):
    print("\n" + "=" * 88 + f"\n{t}\n" + "=" * 88)


def med(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else np.nan


def qs(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if x.size else np.nan


ti = pd.read_csv(os.path.join(TOOLS, "c2_Tinternal.csv"))
# The raw file names are inconsistently cased -- clade3/Clade3, clade4/Clade4,
# para/Para all occur in data/. config.R's canonical_group() normalises them for
# the pipeline; do the same here or a third of the traces go unmatched.
_canon = {g.lower(): g for g in GROUPS}
ti["group"] = ti.group.str.strip().str.lower().map(_canon)
if ti.group.isna().any():
    raise SystemExit("unrecognised group prefix in c2_Tinternal.csv")
der = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables",
                               "derived_N0_R_results_with_carbon.csv"))
der["group"] = der.otu_name.str.extract(r"^([A-Za-z0-9]+)_")[0]
ss = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables",
                              "sharpe_schoolfield_growth_fgC_h_coefs.csv"))
ssm = ss.pivot_table(index="OTU", columns="parameter", values="Estimate")
pb = pd.read_csv(os.path.join(TOOLS, "c2_partB_per_series.csv"))

# one measured trace per (group, set point)
traces = {}
for (grp, Ts), g in ti.groupby(["group", "T_set"]):
    g = g.sort_values("Time")
    # several plates can share a (group, set point); average them on a grid
    grid = np.arange(0.0, min(400.0, g.Time.max()), 0.5)
    ys = []
    for _, gg in g.groupby("File"):
        ys.append(np.interp(grid, gg.Time.values, gg.T_internal.values))
    traces[(grp, int(Ts))] = (grid, np.mean(ys, axis=0))
print(f"measured traces: {len(traces)} (group x set point) combinations")


# ---------------------------------------------------------------------------
# B5m-1. the measured ramp itself
# ---------------------------------------------------------------------------
head("B5m-1. the MEASURED thermal trajectory (T_internal), by set point")
print(f"  {'T_set':>5} {'plates':>6} | {'t=0':>7} {'t=20':>7} {'t=40':>7} {'t=75':>7} "
      f"{'t=120':>7} | {'asympt':>7} | {'deficit at':>11} | {'mean deficit':>13}")
print(f"  {'':>5} {'':>6} | {'':>7} {'':>7} {'':>7} {'':>7} {'':>7} | {'':>7} | "
      f"{'t=75 (K)':>11} | {'0-75 min (K)':>13}")
ramp_prof = {}
for Ts in sorted(ti.T_set.unique()):
    ks = [k for k in traces if k[1] == Ts]
    grid = traces[ks[0]][0]
    Y = np.mean([np.interp(grid, traces[k][0], traces[k][1]) for k in ks], axis=0)
    asy = float(Y[(grid >= 250)].mean()) if (grid >= 250).any() else float(Y[-50:].mean())
    at = lambda q: float(np.interp(q, grid, Y))
    m = grid <= 75
    mean_def = float(asy - Y[m].mean())
    ramp_prof[int(Ts)] = dict(n_plates=len(ks), t0=at(0), t20=at(20), t40=at(40),
                              t75=at(75), t120=at(120), asymptote=asy,
                              deficit_at_75=float(asy - at(75)), mean_deficit_0_75=mean_def)
    p = ramp_prof[int(Ts)]
    print(f"  {Ts:>5} {len(ks):>6} | {p['t0']:>7.2f} {p['t20']:>7.2f} {p['t40']:>7.2f} "
          f"{p['t75']:>7.2f} {p['t120']:>7.2f} | {asy:>7.2f} | {p['deficit_at_75']:>11.2f} | "
          f"{mean_def:>13.2f}")
R["measured_ramp"] = ramp_prof

s_, i_, rv_, p_, se_ = stats.linregress([t for t in ramp_prof],
                                        [ramp_prof[t]["mean_deficit_0_75"] for t in ramp_prof])
print(f"\n  mean deficit over 0-75 min ~ T_set: slope {s_:+.4f} K/degC "
      f"(p = {p_:.3g}), R2 = {rv_**2:.3f}")
print("  -> the excursion is a DIP AND RECOVERY on a ~1 h timescale, not a warm-up")
print("     from ambient; and it is larger at higher set points.")
R["deficit_vs_Tset"] = dict(slope=float(s_), p=float(p_), r2=float(rv_ ** 2))


# ---------------------------------------------------------------------------
# B5m-2. ramp-aware N0 from the measured trajectory
# ---------------------------------------------------------------------------
head("B5m-2. ramp-aware N0 using the measured trajectory and each isolate's own TPC")


def ss_pred_log(TK, lnB0, E, Eh, Th):
    """config.R pred_ss_log(), verbatim."""
    return (lnB0 - (E * 11604.51812) * ((1.0 / TK) - (1.0 / 293.15))
            - np.log1p(np.exp((Eh * 11604.51812) * ((1.0 / Th) - (1.0 / TK)))))


rows = []
for _, rr in der.iterrows():
    otu = int(rr.OTU); Tset = int(rr["T"]); grp = rr.group
    delta = float(rr.delta_Ninoc_to_N0_min); r_set = float(rr.r)
    o = dict(T=Tset, OTU=otu, Replicate=rr.Replicate, group=grp,
             delta=delta, r_set=r_set)
    key = (grp, Tset)
    if key not in traces or otu not in ssm.index or not np.isfinite(r_set) or r_set <= 0 \
       or not np.isfinite(delta) or delta <= 0:
        o.update(status="no_input", n0_factor=1.0); rows.append(o); continue
    c = ssm.loc[otu]
    p = (float(c.lnB0), float(c.E), float(c.Eh), float(c.Th))
    if not all(np.isfinite(p)):
        o.update(status="bad_tpc", n0_factor=1.0); rows.append(o); continue

    grid, Y = traces[key]
    tt = np.linspace(0.0, delta, 400)
    Tt = np.interp(tt, grid, Y)
    asy = float(Y[(grid >= 250)].mean()) if (grid >= 250).any() else float(Y[-50:].mean())
    # the isolate's TPC gives only the SHAPE; scale it through this series' own
    # fitted r at the temperature the plate actually settled at.
    lg_t = ss_pred_log(Tt + 273.15, *p)
    lg_set = float(ss_pred_log(np.array([asy + 273.15]), *p)[0])
    if not np.isfinite(lg_set):
        o.update(status="tpc_nonfinite", n0_factor=1.0); rows.append(o); continue
    r_t = r_set * np.exp(lg_t - lg_set)
    r_t = np.where(np.isfinite(r_t) & (r_t > 0), r_t, 0.0)
    integral = float(np.trapz(r_t, tt))
    o.update(status="ok", rdelta_current=r_set * delta, rdelta_ramp=integral,
             n0_factor=float(np.exp(integral - r_set * delta)),
             T_at_delta=float(Tt[-1]), T_mean_0_delta=float(Tt.mean()),
             T_asymptote=asy)
    rows.append(o)

RM = pd.DataFrame(rows)
ok = RM[RM.status == "ok"]
print(f"  computed for {len(ok)} of {len(RM)} series "
      f"({', '.join(f'{s}:{n}' for s, n in RM.status.value_counts().items() if s!='ok') or 'no failures'})")

print(f"\n  {'T':>4} {'n':>4} | {'mean T over':>12} | {'settle T':>9} | "
      f"{'N0_ramp / N0_current':>26} | {'r*delta cur -> ramp':>22}")
print(f"  {'':>4} {'':>4} | {'0..delta':>12} | {'':>9} | {'median [IQR]':>26} | {'medians':>22}")
byT = {}
for T, g in ok.groupby("T"):
    byT[int(T)] = dict(n=int(len(g)), Tmean=med(g.T_mean_0_delta), Tasy=med(g.T_asymptote),
                       ratio=med(g.n0_factor), q25=qs(g.n0_factor, .25),
                       q75=qs(g.n0_factor, .75),
                       rd_cur=med(g.rdelta_current), rd_ramp=med(g.rdelta_ramp))
    x = byT[int(T)]
    print(f"  {T:>4} {x['n']:>4} | {x['Tmean']:>12.2f} | {x['Tasy']:>9.2f} | "
          f"{x['ratio']:>10.4f} [{x['q25']:>6.4f}, {x['q75']:>6.4f}] | "
          f"{x['rd_cur']:>10.3f} -> {x['rd_ramp']:<9.3f}")
R["ramp_by_T_measured"] = byT

print(f"\n  by taxon (median N0_ramp / N0_current):")
for gn in GROUPS:
    g = ok[ok.group == gn]
    if len(g):
        print(f"      {gn:>8} n={len(g):>4}  {med(g.n0_factor):>8.4f} "
              f"[{qs(g.n0_factor,.25):.4f}, {qs(g.n0_factor,.75):.4f}]")
R["ramp_by_taxon_measured"] = {gn: dict(n=int((ok.group == gn).sum()),
                                        ratio=med(ok.n0_factor[ok.group == gn]))
                               for gn in GROUPS if (ok.group == gn).any()}

# how does this compare with the MODELLED ramp it replaces?
prev = os.path.join(TOOLS, "c2_ramp_n0_factors_modelled.csv")
old = os.path.join(TOOLS, "c2_ramp_n0_factors.csv")
# compare against the MODELLED file specifically -- once it has been archived,
# `old` already holds a measured version and comparing to it says nothing.
_cmp_src = prev if os.path.exists(prev) else old
if os.path.exists(_cmp_src):
    O = pd.read_csv(_cmp_src)[["T", "OTU", "Replicate", "n0_factor"]].rename(
        columns={"n0_factor": "n0_modelled"})
    cmp_ = ok.merge(O, on=["T", "OTU", "Replicate"], how="inner")
    if len(cmp_):
        head("B5m-3. measured vs the modelled ramp it replaces")
        print(f"  {'T':>4} | {'modelled (T0=ambient, tau=19min)':>33} | {'MEASURED':>10}")
        for T, g in cmp_.groupby("T"):
            print(f"  {T:>4} | {med(g.n0_modelled):>33.4f} | {med(g.n0_factor):>10.4f}")
        print(f"\n  the modelled version overstated the ramp: it assumed the plate started")
        print(f"  at 20 C, whereas it measurably starts within ~1.6 C of its asymptote.")
        R["modelled_vs_measured"] = {int(T): dict(modelled=med(g.n0_modelled),
                                                  measured=med(g.n0_factor))
                                     for T, g in cmp_.groupby("T")}


# ---------------------------------------------------------------------------
# B5m-4. second iteration
# ---------------------------------------------------------------------------
head("B5m-4. the circularity: size of the second iteration")
print("  The growth TPC is N0-INVARIANT (r comes from the oxygen curve shape; N0")
print("  enters only respiration = K/N0), so there is no feedback through growth.")
print("  The real coupling is series whose fit window opens before the plate has")
print("  settled, where the fitted r is itself a ramp average.")
tol = 0.25
ok2 = ok.copy()
ok2["deficit_at_delta"] = ok2.T_asymptote - ok2.T_at_delta
n_uns = int((ok2.deficit_at_delta > tol).sum())
print(f"\n  series still more than {tol} K below the settled temperature at t = delta: "
      f"{n_uns} of {len(ok2)} ({100*n_uns/len(ok2):.1f}%)")
ok2["r_iter2"] = ok2.r_set * (ok2.rdelta_current / ok2.rdelta_ramp).replace([np.inf, -np.inf], np.nan)
ok2["ratio_iter2"] = np.exp(ok2.rdelta_ramp * (ok2.r_iter2 / ok2.r_set)
                            - ok2.r_iter2 * ok2.delta)
chg = (ok2.ratio_iter2 - ok2.n0_factor).abs()
print(f"  |ratio(iter 2) - ratio(iter 1)|: median {med(chg):.5f}, p95 {qs(chg,.95):.5f}, "
      f"max {np.nanmax(chg):.5f}")
print("  -> ARM 2 is a BOUND on the size and direction of the ramp error, not an estimator.")
R["iteration_measured"] = dict(n_unsettled=n_uns, n=int(len(ok2)),
                               med=med(chg), p95=qs(chg, .95), max=float(np.nanmax(chg)))


# ---------------------------------------------------------------------------
# B5m-5. resulting E_R
# ---------------------------------------------------------------------------
head("B5m-5. what the measured ramp does to E_R")
K_B = 8.617e-5; T_REF = 293.15
D = der.merge(ok[["T", "OTU", "Replicate", "n0_factor"]], on=["T", "OTU", "Replicate"],
              how="left")
D["n0_factor"] = D.n0_factor.fillna(1.0)
D["boltz"] = 1 / (K_B * T_REF) - 1 / (K_B * (D["T"] + 273.15))
D["logR"] = np.log(D.respiration_fgC_h)
D["logR_ramp"] = D.logR - np.log(D.n0_factor)
D["logR_nobp"] = D.logR + D.r * D.delta_Ninoc_to_N0_min
print(f"  {'group':>8} {'n':>4} | {'ramp-aware':>11} | {'shipped':>9} | "
      f"{'no-backproj':>12} | {'shift (ramp)':>13}")
er = {}
for gn in GROUPS:
    g = D[(D.group == gn) & np.isfinite(D.logR)]
    if len(g) < 5:
        continue
    s0, *_ = stats.linregress(g.boltz, g.logR)
    s1, *_ = stats.linregress(g.boltz, g.logR_ramp)
    s2, *_ = stats.linregress(g.boltz, g.logR_nobp)
    er[gn] = dict(E_R_shipped=float(s0), E_R_ramp=float(s1), E_R_nobp=float(s2),
                  shift=float(s1 - s0))
    print(f"  {gn:>8} {len(g):>4} | {s1:>11.3f} | {s0:>9.3f} | {s2:>12.3f} | {s1-s0:>+13.3f}")
R["E_R_measured_ramp"] = er

# ---------------------------------------------------------------------------
# write the ARM 2 input (measured)
# ---------------------------------------------------------------------------
if os.path.exists(old) and not os.path.exists(prev):
    os.rename(old, prev)
    print(f"\n  kept the modelled version as {os.path.basename(prev)}")
arm2 = RM[["T", "OTU", "Replicate", "n0_factor", "status"]].copy()
arm2["n0_factor"] = arm2.n0_factor.where(np.isfinite(arm2.n0_factor), 1.0)
arm2.to_csv(old, index=False)
print(f"wrote {old}  ({len(arm2)} series; "
      f"{int((arm2.status!='ok').sum())} fall back to 1.0)")

with open(os.path.join(TOOLS, "c2_partB5_measured.json"), "w") as f:
    json.dump(R, f, indent=1, default=float)
print(f"wrote {os.path.join(TOOLS,'c2_partB5_measured.json')}")
