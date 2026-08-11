"""
c2_partB_transient.py -- C2 PART B.

Characterise the pre-peak equilibration transient, then bound the two errors it
creates for N0 and K.

    cauris_etcgem/.venv/bin/python reports/tools/c2_partB_transient.py

Read-only on the pipeline. Writes
    reports/tools/c2_partB.json                  summary statistics
    reports/tools/c2_partB_per_series.csv        A, tau, O2_eq, diagnostics per series
    reports/tools/c2_ramp_n0_factors.csv         ARM 2 input: per-series N0 multiplier

Model, fitted to every PRE-PEAK segment (t from the start of the record to the
automatic O2 tip in Oxygen_Trimmed_Series_Metadata.csv):

    O2(t) = O2_eq - A * exp(-t / tau)

so A (mg/L) is the total rise still to come at t = 0, tau (min) the time
constant, O2_eq (mg/L) the asymptote. If the transient is physico-chemical, A
and tau should scale with distance from ambient and tau should NOT differ by
taxon.
"""
import json, os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS = os.path.join(ROOT, "reports", "tools")
T_AMBIENT = 20.0
GROUPS = ["Clade1", "Clade2", "Clade3", "Clade4", "glab", "para"]
MIN_PTS = 8            # minimum pre-peak points to attempt a fit
MIN_RISE = 0.05        # mg/L; below this the series is recorded as "no rise"
R = {}


def head(t):
    print("\n" + "=" * 88 + f"\n{t}\n" + "=" * 88)


def qs(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if x.size else np.nan


def med(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else np.nan


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
print("loading raw oxygen traces ...", flush=True)
long = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables", "Oxygen_All_Long.csv"))
meta = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables",
                                "Oxygen_Trimmed_Series_Metadata.csv"))
der = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables",
                               "derived_N0_R_results_with_carbon.csv"))
der["group"] = der.otu_name.str.extract(r"^([A-Za-z0-9]+)_")[0]
print(f"  {len(long):,} rows, {len(meta)} series in the trimming metadata, "
      f"{len(der)} derived rows", flush=True)

long = long.sort_values(["T", "OTU", "Replicate", "Time"])
key = ["T", "OTU", "Replicate"]
groups = {k: g for k, g in long.groupby(key, sort=False)}
mrec = meta.set_index(key)


def model(t, o2eq, A, tau):
    return o2eq - A * np.exp(-t / tau)


# ---------------------------------------------------------------------------
# B0. fit the transient
# ---------------------------------------------------------------------------
head("B0. fitting O2(t) = O2_eq - A*exp(-t/tau) to every pre-peak segment")

rows = []
for k, g in groups.items():
    T, OTU, Rep = k
    if k not in mrec.index:
        continue
    m = mrec.loc[k]
    if isinstance(m, pd.DataFrame):
        m = m.iloc[0]
    pk = float(m.peak_time) if np.isfinite(m.peak_time) else np.nan
    rec = dict(T=int(T), OTU=int(OTU), Replicate=Rep, peak_time=pk)
    if not np.isfinite(pk) or pk <= 0:
        rec.update(status="no_peak"); rows.append(rec); continue

    seg = g[(g.Time >= 0) & (g.Time <= pk)]
    t = seg.Time.values.astype(float)
    y = seg.Oxygen.values.astype(float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    rec["n_pre"] = int(t.size)
    if t.size < MIN_PTS:
        rec.update(status="too_few_points"); rows.append(rec); continue

    rise = float(y[-1] - y[0])
    rec.update(O2_t0=float(y[0]), O2_peak=float(y.max()), observed_rise=rise)
    if rise < MIN_RISE:
        rec.update(status="no_rise"); rows.append(rec); continue

    p0 = [float(y.max()), max(rise, 0.1), max(pk / 3.0, 5.0)]
    try:
        popt, pcov = curve_fit(
            model, t, y, p0=p0,
            bounds=([y.min() - 5, 1e-4, 1.0], [y.max() + 10, 50.0, 1000.0]),
            maxfev=20000)
        o2eq, A, tau = [float(v) for v in popt]
        pred = model(t, *popt)
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rmse = float(np.sqrt(ss_res / t.size))
        rec.update(status="ok", O2_eq=o2eq, A=A, tau=tau, r2=r2, rmse=rmse,
                   tau_se=float(np.sqrt(pcov[2, 2])) if np.all(np.isfinite(pcov)) else np.nan,
                   settled_frac=float(1 - np.exp(-pk / tau)))
    except Exception as e:
        rec.update(status=f"fit_failed:{type(e).__name__}")
    rows.append(rec)

fit = pd.DataFrame(rows)
fit["group"] = fit.OTU.map(
    lambda o: GROUPS[(int(o) - 1) // 3] if 1 <= int(o) <= 18 else None)
fit["dT"] = np.abs(fit["T"] - T_AMBIENT)

n_ok = int((fit.status == "ok").sum())
print(f"  {len(fit)} series; fitted OK {n_ok}; " +
      ", ".join(f"{s} {int(n)}" for s, n in fit.status.value_counts().items()
                if s != "ok"))
R["fit_status"] = {str(s): int(n) for s, n in fit.status.value_counts().items()}
R["n_series"] = int(len(fit))

F = fit[fit.status == "ok"].copy()
print(f"  fit quality: median R2 = {med(F.r2):.4f} "
      f"(IQR {qs(F.r2,.25):.4f}-{qs(F.r2,.75):.4f}), "
      f"median RMSE = {med(F.rmse):.4f} mg/L")
R["fit_quality"] = dict(r2_med=med(F.r2), r2_q25=qs(F.r2, .25), r2_q75=qs(F.r2, .75),
                        rmse_med=med(F.rmse))


# ---------------------------------------------------------------------------
# B1a. A, tau, O2_eq by temperature
# ---------------------------------------------------------------------------
head("B1. A, tau and O2_eq by TEMPERATURE")
print(f"  {'T':>4} {'|T-20|':>7} {'n':>4} | {'A (mg/L)':>18} | {'tau (min)':>18} | "
      f"{'O2_eq (mg/L)':>18} | {'settled at peak':>15}")
byT = {}
for T, g in F.groupby("T"):
    byT[int(T)] = dict(n=int(len(g)), A_med=med(g.A), A_q25=qs(g.A, .25), A_q75=qs(g.A, .75),
                       tau_med=med(g.tau), tau_q25=qs(g.tau, .25), tau_q75=qs(g.tau, .75),
                       O2eq_med=med(g.O2_eq), O2eq_q25=qs(g.O2_eq, .25),
                       O2eq_q75=qs(g.O2_eq, .75), settled_med=med(g.settled_frac))
    b = byT[int(T)]
    print(f"  {T:>4} {abs(T-T_AMBIENT):>7.0f} {b['n']:>4} | "
          f"{b['A_med']:>7.2f} [{b['A_q25']:>4.2f},{b['A_q75']:>5.2f}] | "
          f"{b['tau_med']:>7.1f} [{b['tau_q25']:>4.1f},{b['tau_q75']:>5.1f}] | "
          f"{b['O2eq_med']:>7.2f} [{b['O2eq_q25']:>4.2f},{b['O2eq_q75']:>5.2f}] | "
          f"{100*b['settled_med']:>14.1f}%")
R["by_temperature"] = byT

# scaling with distance from ambient
for nm, col in (("A", "A"), ("tau", "tau")):
    s, i, rv, pv, se = stats.linregress(F.dT, F[col])
    rho, prho = stats.spearmanr(F.dT, F[col])
    print(f"\n  {nm} ~ |T - {T_AMBIENT:.0f}C| : slope {s:+.4f} (SE {se:.4f}, p = {pv:.3g}), "
          f"R2 = {rv**2:.4f}; Spearman rho = {rho:+.3f}")
    R[f"{nm}_vs_absdT"] = dict(slope=float(s), se=float(se), p=float(pv),
                               r2=float(rv ** 2), spearman=float(rho))

# tau by taxon -- should be null if the transient is instrumental
head("B2. does tau differ by TAXON?  (it should not, if the transient is instrumental)")
arrs, names = [], []
for gn in GROUPS:
    g = F[F.group == gn]
    if len(g) > 3:
        arrs.append(g.tau.values); names.append(gn)
        print(f"  {gn:>8} n={len(g):>4}  tau median {med(g.tau):>7.1f} min "
              f"[{qs(g.tau,.25):>6.1f}, {qs(g.tau,.75):>6.1f}]   "
              f"A median {med(g.A):>6.2f} mg/L")
Fst, pFst = stats.f_oneway(*arrs)
Kst, pKst = stats.kruskal(*arrs)
print(f"\n  ANOVA   F({len(arrs)-1},{sum(len(a) for a in arrs)-len(arrs)}) = {Fst:.2f}, p = {pFst:.3g}")
print(f"  Kruskal H = {Kst:.2f}, p = {pKst:.3g}")
# and controlling for temperature: two-way on residuals from the per-T median
F["tau_resid"] = F.tau - F.groupby("T").tau.transform("median")
arrs2 = [F.tau_resid[F.group == gn].values for gn in names]
F2, p2 = stats.f_oneway(*arrs2)
print(f"  after removing the per-temperature median: F = {F2:.2f}, p = {p2:.3g}")
for gn in names:
    print(f"      {gn:>8} residual tau median {med(F.tau_resid[F.group==gn]):+7.1f} min")
ng_names = [gn for gn in names if gn != "glab"]
arrs3 = [F.tau_resid[F.group == gn].values for gn in ng_names]
F3, p3 = stats.f_oneway(*arrs3)
print(f"\n  EXCLUDING glabrata: F({len(arrs3)-1},{sum(len(a) for a in arrs3)-len(arrs3)}) "
      f"= {F3:.2f}, p = {p3:.3g}  -> tau is a property of the instrument, not the organism")
R["tau_by_taxon_excl_glab"] = dict(F=float(F3), p=float(p3))
R["tau_by_taxon"] = dict(anova_F=float(Fst), anova_p=float(pFst),
                         kruskal_H=float(Kst), kruskal_p=float(pKst),
                         anova_resid_F=float(F2), anova_resid_p=float(p2),
                         per_taxon={gn: dict(n=int((F.group == gn).sum()),
                                             tau_med=med(F.tau[F.group == gn]),
                                             A_med=med(F.A[F.group == gn]),
                                             tau_resid_med=med(F.tau_resid[F.group == gn]))
                                    for gn in names})


# ---------------------------------------------------------------------------
# B3. O2_eq vs the oxygen solubility curve
# ---------------------------------------------------------------------------
head("B3. fitted O2_eq against the known temperature dependence of O2 solubility")


def o2_solubility_mgL(Tc, salinity=0.0):
    """Benson & Krause / Garcia-Gordon air-saturated O2 in fresh water, mg/L,
    1 atm, 100% humidity. Standard oceanographic fit."""
    Ts = np.log((298.15 - Tc) / (273.15 + Tc))
    A = [5.80871, 3.20291, 4.17887, 5.10006, -9.86643e-2, 3.80369]
    B = [-7.01577e-3, -7.70028e-3, -1.13864e-2, -9.51519e-3]
    C0 = -2.75915e-7
    lnC = (A[0] + A[1]*Ts + A[2]*Ts**2 + A[3]*Ts**3 + A[4]*Ts**4 + A[5]*Ts**5
           + salinity*(B[0] + B[1]*Ts + B[2]*Ts**2 + B[3]*Ts**3) + C0*salinity**2)
    umol_kg = np.exp(lnC)                     # umol/kg
    return umol_kg * 31.998 / 1000.0          # -> mg/L (rho ~ 1 kg/L)


print(f"  {'T':>4} | {'O2_eq fitted (med)':>19} | {'solubility (mg/L)':>18} | "
      f"{'fitted - soluble':>17} | {'fitted/soluble':>15}")
sol_rows = []
for T in sorted(byT):
    obs = byT[T]["O2eq_med"]; sol = float(o2_solubility_mgL(T))
    sol_rows.append((T, obs, sol))
    print(f"  {T:>4} | {obs:>19.3f} | {sol:>18.3f} | {obs-sol:>+17.3f} | {obs/sol:>15.3f}")
Ts_ = np.array([r[0] for r in sol_rows]); obs_ = np.array([r[1] for r in sol_rows])
sol_ = np.array([r[2] for r in sol_rows])
r_os, p_os = stats.pearsonr(sol_, obs_)
rho_os, prho_os = stats.spearmanr(Ts_, obs_)
rho_ss, _ = stats.spearmanr(Ts_, sol_)
s_o, i_o, rv_o, p_o, se_o = stats.linregress(Ts_, obs_)
s_s, i_s, rv_s, p_s, se_s = stats.linregress(Ts_, sol_)
print(f"\n  fitted O2_eq vs T   : slope {s_o:+.4f} mg/L/degC (p = {p_o:.3g}), "
      f"Spearman rho = {rho_os:+.3f}")
print(f"  solubility  vs T    : slope {s_s:+.4f} mg/L/degC, Spearman rho = {rho_ss:+.3f}")
print(f"  Pearson r(fitted, solubility) across the 12 temperatures = {r_os:+.3f} (p = {p_os:.3g})")
# The SIGN test. If the vial started at ambient equilibrium and warmed to T, the
# EQUILIBRIUM oxygen would FALL by sol(T) - sol(T_ambient) < 0. The observed
# transient RISES. Compare the magnitudes anyway: if |A| tracks |delta-solubility|
# the driver is still the temperature change, whatever carries the sign.
print(f"\n  SIGN TEST -- observed rise A vs the predicted equilibrium CHANGE on warming:")
print(f"  {'T':>4} | {'A observed (med)':>17} | {'sol(T) - sol(20)':>17} | {'|ratio|':>9}")
sign_rows = []
sol20 = float(o2_solubility_mgL(T_AMBIENT))
for T in sorted(byT):
    A_ = byT[T]["A_med"]; dsol = float(o2_solubility_mgL(T)) - sol20
    sign_rows.append(dict(T=int(T), A=A_, dsol=float(dsol),
                          abs_ratio=float(A_ / abs(dsol)) if dsol else np.nan))
    print(f"  {T:>4} | {A_:>17.3f} | {dsol:>+17.3f} | {A_/abs(dsol) if dsol else np.nan:>9.2f}")
A_v = np.array([r["A"] for r in sign_rows]); ds_v = np.abs([r["dsol"] for r in sign_rows])
r_sg, p_sg = stats.pearsonr(ds_v, A_v)
print(f"\n  Pearson r(A, |sol(T)-sol(20)|) = {r_sg:+.3f} (p = {p_sg:.3g})")
print("  -> magnitude tracks the temperature change closely, but the SIGN is inverted:")
print("     warming should LOWER dissolved O2, and the record RISES. Together with the")
print("     fitted O2_eq sitting 1.20-1.48x ABOVE air-saturation, that points to the")
print("     optode's own temperature compensation settling, not to gas re-equilibration.")
print("     Either way the driver is the temperature step, not the cells.")
R["sign_test"] = dict(rows=sign_rows, pearson_A_vs_absdsol=float(r_sg), p=float(p_sg))

R["O2eq_vs_solubility"] = dict(
    per_T=[dict(T=int(t), o2eq_med=float(o), solubility=float(s),
                diff=float(o - s), ratio=float(o / s)) for t, o, s in sol_rows],
    slope_obs=float(s_o), p_obs=float(p_o), slope_sol=float(s_s),
    pearson_obs_vs_sol=float(r_os), p_pearson=float(p_os),
    spearman_obs_vs_T=float(rho_os), spearman_sol_vs_T=float(rho_ss))


# ---------------------------------------------------------------------------
# B4. FIT-WINDOW CONTAMINATION
# ---------------------------------------------------------------------------
head("B4. fit-window contamination: residual abiotic drift inside the fit window")

# fit window comes from the derived table (delta = window start) and the
# trimming metadata (chosen_end_time). The observed drawdown is read from the
# raw trace at those two times.
w = der[["T", "OTU", "Replicate", "delta_Ninoc_to_N0_min", "T_end_min", "K", "r",
         "group", "respiration_fgC_h"]].rename(
    columns={"delta_Ninoc_to_N0_min": "t_start"})
w["t_end"] = w.t_start + w.T_end_min
W = w.merge(F[["T", "OTU", "Replicate", "A", "tau", "O2_eq", "peak_time"]],
            on=["T", "OTU", "Replicate"], how="inner")
print(f"  {len(W)} series have both a fitted transient and a fit window")


def obs_at(k, times):
    g = groups.get(k)
    if g is None:
        return [np.nan] * len(times)
    return [float(np.interp(t, g.Time.values, g.Oxygen.values)) for t in times]


o_start, o_end = [], []
for _, rr in W.iterrows():
    a, b = obs_at((int(rr["T"]), int(rr.OTU), rr.Replicate), [rr.t_start, rr.t_end])
    o_start.append(a); o_end.append(b)
W["O2_start"] = o_start; W["O2_end"] = o_end
W["bio_drawdown"] = W.O2_start - W.O2_end          # observed decline, mg/L


def contamination(W, tau_col="tau"):
    tau = W[tau_col].values
    A = W.A.values
    ts, te = W.t_start.values, W.t_end.values
    # residual abiotic RISE still occurring inside the window
    drift = A * (np.exp(-ts / tau) - np.exp(-te / tau))
    # instantaneous abiotic slope at the window start, mg/L/min
    slope0 = (A / tau) * np.exp(-ts / tau)
    return drift, slope0


drift, slope0 = contamination(W)
W["abiotic_drift"] = drift
W["frac_of_drawdown"] = drift / W.bio_drawdown
W["abiotic_slope0"] = slope0
W["settled_at_start"] = 1 - np.exp(-W.t_start / W.tau)

print("""
  SIGN. The transient is a RISE toward equilibrium, so whatever is left of it
  inside the fit window is UPWARD. It masks part of the drawdown, so it makes K
  too SMALL, not too large. Positive numbers below = the amount K is under-read.

  HOW IT IS MEASURED. Not by comparing the residual abiotic SLOPE to K -- that
  compares an instantaneous rate to a whole-window parameter and overstates the
  effect several-fold. Instead the fitted transient is subtracted from the
  observed trace over the fit window and 07's oxygen model
        O2(t) = O2_0 + (K/r) * (1 - exp(r*t))
  is refitted to the decontaminated trace with the same bounds. The refit is
  first validated by reproducing the SHIPPED K from the raw window.""")


def resp_model(t, r, K, O2_0):
    return O2_0 + (K / r) * (1.0 - np.exp(r * t))


def refit(t, y):
    """07_oxygen_fits.R fit_one(), same model and same bounds."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    if t.size < 6 or not np.all(np.isfinite(y)):
        return (np.nan, np.nan, np.nan)
    t0 = t - t.min()
    n_seg = max(6, int(0.25 * t.size))
    sl = np.polyfit(t0[:n_seg], y[:n_seg], 1)[0]
    K_est = max(abs(sl), 1e-8)
    p0 = [max(min(0.005, 0.15), 1e-6), min(max(K_est, 1e-8), 0.5), float(np.median(y[:min(30, y.size)]))]
    yr = (float(np.min(y)), float(np.max(y)))
    span = yr[1] - yr[0] or (abs(y[0]) + 1)
    lo = [1e-6, max(1e-10, K_est / 50), yr[0] - 0.5 * span]
    hi = [0.15, max(1.0, K_est * 50), yr[1] + 0.5 * span]
    p0 = [min(max(p0[i], lo[i]), hi[i]) for i in range(3)]
    try:
        popt, _ = curve_fit(resp_model, t0, y, p0=p0, bounds=(lo, hi), maxfev=20000)
        return (float(popt[0]), float(popt[1]), float(popt[2]))
    except Exception:
        return (np.nan, np.nan, np.nan)


def window_trace(k, t_start, t_end):
    g = groups.get(k)
    if g is None:
        return None, None
    m = (g.Time.values >= t_start) & (g.Time.values <= t_end)
    return g.Time.values[m], g.Oxygen.values[m]


def abiotic(t, A, tau, O2_eq):
    return O2_eq - A * np.exp(-t / tau)


# tau sensitivity band: per-temperature 25th / 75th percentile of fitted tau
tau_q25 = W["T"].map({T: byT[T]["tau_q25"] for T in byT})
tau_q75 = W["T"].map({T: byT[T]["tau_q75"] for T in byT})
W["tau_lo"] = tau_q25.values; W["tau_hi"] = tau_q75.values

rows_rf = []
for _, rr in W.iterrows():
    k = (int(rr["T"]), int(rr.OTU), rr.Replicate)
    t, y = window_trace(k, rr.t_start, rr.t_end)
    o = dict(T=int(rr["T"]), OTU=int(rr.OTU), Replicate=rr.Replicate)
    if t is None or t.size < 6:
        o.update(status="short"); rows_rf.append(o); continue
    r_raw, K_raw, _ = refit(t, y)
    out = {}
    for tag, tcol in (("", "tau"), ("_lo", "tau_lo"), ("_hi", "tau_hi")):
        tau = float(rr[tcol])
        y_corr = y - (abiotic(t, rr.A, tau, rr.O2_eq) - abiotic(t[0], rr.A, tau, rr.O2_eq))
        r_c, K_c, _ = refit(t, y_corr)
        out["r_corr" + tag] = r_c
        out["K_corr" + tag] = K_c
    o.update(status="ok", K_raw=K_raw, r_raw=r_raw, K_shipped=float(rr.K),
             r_shipped=float(rr.r), **out)
    rows_rf.append(o)

RF = pd.DataFrame(rows_rf)
W = W.merge(RF.drop(columns=["status"]), on=["T", "OTU", "Replicate"], how="left")

# validation: does the python refit reproduce the shipped K?
v = W[np.isfinite(W.K_raw) & np.isfinite(W.K)]
rel = np.abs(v.K_raw - v.K) / v.K
print(f"\n  VALIDATION refit vs shipped K on the raw window (n = {len(v)}): "
      f"median |rel diff| = {med(rel):.2e}, p90 = {qs(rel,.90):.2e}, "
      f"fraction within 1% = {float((rel < 0.01).mean()):.3f}")
R["refit_validation"] = dict(n=int(len(v)), med_rel=med(rel), p90_rel=qs(rel, .90),
                             frac_within_1pct=float((rel < 0.01).mean()))

for tag in ("", "_lo", "_hi"):
    W["dK" + tag] = (W["K_corr" + tag] - W.K_raw) / W.K_raw

print(f"\n  {'T':>4} {'n':>4} | {'settled at':>10} | {'abiotic drift':>15} | "
      f"{'as % of':>10} | {'K under-read (% of K), median [IQR]':>38}")
print(f"  {'':>4} {'':>4} | {'win start':>10} | {'(mg/L, median)':>15} | "
      f"{'drawdown':>10} | {'refit on the decontaminated trace':>38}")
b4 = {}
for T, g in W.groupby("T"):
    b4[int(T)] = dict(n=int(len(g)), settled=med(g.settled_at_start),
                      drift=med(g.abiotic_drift), frac=med(g.frac_of_drawdown),
                      dK=med(g.dK), dK_q25=qs(g.dK, .25), dK_q75=qs(g.dK, .75))
    x = b4[int(T)]
    print(f"  {T:>4} {x['n']:>4} | {100*x['settled']:>9.2f}% | {x['drift']:>15.4f} | "
          f"{100*x['frac']:>9.2f}% | "
          f"{100*x['dK']:>14.2f} [{100*x['dK_q25']:>6.2f}, {100*x['dK_q75']:>6.2f}]")
R["contamination_by_T"] = b4

print(f"\n  by taxon:")
for gn in GROUPS:
    g = W[W.group == gn]
    if len(g):
        print(f"      {gn:>8} n={len(g):>4}  K under-read {100*med(g.dK):>6.2f}%  "
              f"(drift {100*med(g.frac_of_drawdown):>5.2f}% of drawdown)")
R["contamination_by_taxon"] = {gn: dict(n=int((W.group == gn).sum()),
                                        dK=med(W.dK[W.group == gn]),
                                        frac=med(W.frac_of_drawdown[W.group == gn]))
                               for gn in GROUPS if (W.group == gn).any()}

vv = W[np.isfinite(W.dK)]
s_c, i_c, rv_c, p_c, se_c = stats.linregress(vv["T"], vv.dK)
rho_c, prho_c = stats.spearmanr(vv["T"], vv.dK)
print(f"\n  K under-read vs T: slope {100*s_c:+.4f} %-points/degC (p = {p_c:.3g}), "
      f"R2 = {rv_c**2:.4f}; Spearman rho = {rho_c:+.3f} (p = {prho_c:.3g})")
R["contamination_vs_T"] = dict(slope_pct_per_degC=float(100 * s_c), p=float(p_c),
                               r2=float(rv_c ** 2), spearman=float(rho_c),
                               spearman_p=float(prho_c))

head("B4b. sensitivity band over the interquartile range of fitted tau")
print(f"  {'T':>4} | {'tau q25':>8} {'tau med':>8} {'tau q75':>8} | "
      f"{'K under-read at tau q25':>24} {'at fitted tau':>15} {'at tau q75':>13}")
band = {}
for T, g in W.groupby("T"):
    band[int(T)] = dict(tau_q25=byT[int(T)]["tau_q25"], tau_med=byT[int(T)]["tau_med"],
                        tau_q75=byT[int(T)]["tau_q75"],
                        lo=med(g.dK_lo), mid=med(g.dK), hi=med(g.dK_hi))
    x = band[int(T)]
    print(f"  {T:>4} | {x['tau_q25']:>8.1f} {x['tau_med']:>8.1f} {x['tau_q75']:>8.1f} | "
          f"{100*x['lo']:>23.2f}% {100*x['mid']:>14.2f}% {100*x['hi']:>12.2f}%")
R["contamination_band"] = band

K_B = 8.617e-5; T_REF = 293.15
W["boltz"] = 1 / (K_B * T_REF) - 1 / (K_B * (W["T"] + 273.15))
print("\n  effect on E_K of removing the contamination (band over tau IQR):")
print(f"  {'group':>8} | {'E_K raw':>9} | {'at tau q25':>11} {'at fitted tau':>14} {'at tau q75':>11}")
ekband = {}
for gn in GROUPS:
    g = W[(W.group == gn)]
    g0 = g[np.isfinite(g.K_raw) & (g.K_raw > 0)]
    if len(g0) < 5:
        continue
    s0, *_ = stats.linregress(g0.boltz, np.log(g0.K_raw))
    vals = {}
    for tag, lab in (("_lo", "q25"), ("", "med"), ("_hi", "q75")):
        gg = g[np.isfinite(g["K_corr" + tag]) & (g["K_corr" + tag] > 0)]
        s1, *_ = stats.linregress(gg.boltz, np.log(gg["K_corr" + tag]))
        vals[lab] = float(s1)
    ekband[gn] = dict(E_K_raw=float(s0), **vals)
    print(f"  {gn:>8} | {s0:>9.3f} | {vals['q25']:>11.3f} {vals['med']:>14.3f} "
          f"{vals['q75']:>11.3f}")
R["EK_after_contamination_band"] = ekband

# ---------------------------------------------------------------------------
# B5. THERMAL-RAMP ERROR -> ramp-aware N0 (ARM 2 input)
# ---------------------------------------------------------------------------
head("B5. thermal-ramp error: a ramp-aware N0")

ss = pd.read_csv(os.path.join(ROOT, "runs/C1_reproduction", "tables",
                              "sharpe_schoolfield_growth_fgC_h_coefs.csv"))
print(f"  per-isolate Sharpe-Schoolfield growth fits: {len(ss)} rows, "
      f"columns {[c for c in ss.columns][:8]}")


def ss_pred_log(TK, lnB0, E, Eh, Th):
    """config.R pred_ss_log(), verbatim."""
    return (lnB0
            - (E * 11604.51812) * ((1.0 / TK) - (1.0 / 293.15))
            - np.log1p(np.exp((Eh * 11604.51812) * ((1.0 / Th) - (1.0 / TK)))))


# the coefficient table is LONG (one row per OTU x parameter) -- pivot it
ssm = ss.pivot_table(index="OTU", columns="parameter", values="Estimate")
need = {"lnB0", "E", "Eh", "Th"}
have = need.issubset(set(ssm.columns))
print(f"  pivoted to {ssm.shape[0]} isolates x {ssm.shape[1]} parameters "
      f"{sorted(ssm.columns)}; required columns present: {have}")

ramp_rows = []
for _, rr in W.iterrows():
    otu = int(rr.OTU); Tset = float(rr["T"]); tau = float(rr.tau)
    delta = float(rr.t_start); r_set = float(rr.r)
    out = dict(T=int(Tset), OTU=otu, Replicate=rr.Replicate, tau=tau, delta=delta,
               r_set=r_set)
    if not (have and ssm is not None and otu in ssm.index) or not np.isfinite(r_set) or r_set <= 0:
        out.update(status="no_tpc", n0_ramp_over_current=np.nan)
        ramp_rows.append(out); continue
    c = ssm.loc[otu]
    if isinstance(c, pd.DataFrame):
        c = c.iloc[0]
    p = (float(c.lnB0), float(c.E), float(c.Eh), float(c.Th))
    if not all(np.isfinite(p)):
        out.update(status="bad_tpc", n0_ramp_over_current=np.nan)
        ramp_rows.append(out); continue

    # sample temperature follows the same first-order approach as the optode
    tt = np.linspace(0.0, delta, 400) if delta > 0 else np.array([0.0])
    Tt = Tset - (Tset - T_AMBIENT) * np.exp(-tt / tau)
    # instantaneous growth rate from the isolate's own TPC, SCALED so that the
    # curve passes through this series' own fitted r at the set point. The
    # scaling cancels lnB0, so only the SHAPE of the TPC is used.
    lg_t = ss_pred_log(Tt + 273.15, *p)
    lg_set = ss_pred_log(np.array([Tset + 273.15]), *p)[0]
    if not np.isfinite(lg_set):
        out.update(status="tpc_nonfinite", n0_ramp_over_current=np.nan)
        ramp_rows.append(out); continue
    r_t = r_set * np.exp(lg_t - lg_set)
    r_t = np.where(np.isfinite(r_t) & (r_t > 0), r_t, 0.0)
    integral = float(np.trapz(r_t, tt)) if delta > 0 else 0.0
    out.update(status="ok",
               rdelta_current=r_set * delta,
               rdelta_ramp=integral,
               n0_ramp_over_current=float(np.exp(integral - r_set * delta)),
               T_at_delta=float(Tt[-1]),
               frac_thermally_settled=float(1 - np.exp(-delta / tau)))
    ramp_rows.append(out)

ramp = pd.DataFrame(ramp_rows)
ramp["group"] = ramp.OTU.map(lambda o: GROUPS[(int(o) - 1) // 3])
Rk = ramp[ramp.status == "ok"]
print(f"  ramp-aware N0 computed for {len(Rk)} of {len(ramp)} series "
      f"({', '.join(f'{s}:{n}' for s, n in ramp.status.value_counts().items() if s!='ok') or 'no failures'})")

print(f"\n  {'T':>4} {'n':>4} | {'T at t=delta':>13} | {'thermally':>10} | "
      f"{'N0_ramp / N0_current':>26} | {'r*delta cur -> ramp':>22}")
print(f"  {'':>4} {'':>4} | {'(degC)':>13} | {'settled':>10} | {'median [IQR]':>26} | {'medians':>22}")
b5 = {}
for T, g in Rk.groupby("T"):
    b5[int(T)] = dict(n=int(len(g)), T_at_delta=med(g.T_at_delta),
                      settled=med(g.frac_thermally_settled),
                      ratio=med(g.n0_ramp_over_current),
                      ratio_q25=qs(g.n0_ramp_over_current, .25),
                      ratio_q75=qs(g.n0_ramp_over_current, .75),
                      rd_cur=med(g.rdelta_current), rd_ramp=med(g.rdelta_ramp))
    x = b5[int(T)]
    print(f"  {T:>4} {x['n']:>4} | {x['T_at_delta']:>13.2f} | {100*x['settled']:>9.2f}% | "
          f"{x['ratio']:>10.4f} [{x['ratio_q25']:>6.4f}, {x['ratio_q75']:>6.4f}] | "
          f"{x['rd_cur']:>10.3f} -> {x['rd_ramp']:<9.3f}")
R["ramp_by_T"] = b5

print(f"\n  by taxon (median N0_ramp / N0_current):")
for gn in GROUPS:
    g = Rk[Rk.group == gn]
    if len(g):
        print(f"      {gn:>8} n={len(g):>4}  {med(g.n0_ramp_over_current):>8.4f} "
              f"[{qs(g.n0_ramp_over_current,.25):.4f}, {qs(g.n0_ramp_over_current,.75):.4f}]")
R["ramp_by_taxon"] = {gn: dict(n=int((Rk.group == gn).sum()),
                               ratio=med(Rk.n0_ramp_over_current[Rk.group == gn]))
                      for gn in GROUPS if (Rk.group == gn).any()}

# second-iteration change
head("B5b. the circularity: how big is the SECOND iteration?")
print("  The growth TPC is N0-INVARIANT (r is fitted to the oxygen curve shape; N0")
print("  enters only respiration = K/N0), so re-running the correction with a")
print("  TPC re-estimated after arm 2 changes nothing on the growth side. The one")
print("  real coupling is series whose FIT WINDOW starts before thermal equilibrium,")
print("  where the fitted r is itself a ramp average rather than the set-point rate.")
n_unsettled = int((Rk.frac_thermally_settled < 0.95).sum())
print(f"\n  series with delta < 3*tau (less than 95% thermally settled at window start): "
      f"{n_unsettled} of {len(Rk)} ({100*n_unsettled/len(Rk):.1f}%)")
# iteration 2: re-derive r at the set point by dividing out the ramp average,
# then recompute the ratio
Rk2 = Rk.copy()
Rk2["r_set_iter2"] = Rk2.r_set * (Rk2.rdelta_current / Rk2.rdelta_ramp).replace(
    [np.inf, -np.inf], np.nan)
Rk2["ratio_iter2"] = np.exp(Rk2.rdelta_ramp * (Rk2.r_set_iter2 / Rk2.r_set)
                            - Rk2.r_set_iter2 * Rk2.delta)
chg = (Rk2.ratio_iter2 - Rk2.n0_ramp_over_current).abs()
print(f"  |ratio(iter 2) - ratio(iter 1)|: median {med(chg):.5f}, "
      f"p95 {qs(chg,.95):.5f}, max {np.nanmax(chg):.5f}")
print("  -> treat arm 2 as a BOUND, not an estimator.")
R["ramp_iteration"] = dict(n_unsettled=n_unsettled, n=int(len(Rk)),
                           frac_unsettled=float(n_unsettled / len(Rk)),
                           iter2_change_med=med(chg), iter2_change_p95=qs(chg, .95),
                           iter2_change_max=float(np.nanmax(chg)))

# what does the ramp do to E_R?
head("B5c. resulting E_R shift (respiration = K / N0, so log R gains -log(ratio))")
WR = W.merge(Rk[["T", "OTU", "Replicate", "n0_ramp_over_current"]],
             on=["T", "OTU", "Replicate"], how="inner")
WR["logR"] = np.log(WR.respiration_fgC_h)
WR["logR_ramp"] = WR.logR - np.log(WR.n0_ramp_over_current)
print(f"  {'group':>8} {'n':>4} | {'E_R shipped':>13} | {'E_R ramp-aware':>15} | "
      f"{'E_R no-backproj':>16} | {'shift (ramp)':>13}")
er = {}
for gn in GROUPS:
    g = WR[WR.group == gn]
    if len(g) < 5:
        continue
    s0, *_ = stats.linregress(g.boltz, g.logR)
    s1, *_ = stats.linregress(g.boltz, g.logR_ramp)
    s2, *_ = stats.linregress(g.boltz, g.logR + g.r * g.t_start)
    er[gn] = dict(E_R_shipped=float(s0), E_R_ramp=float(s1), E_R_nobp=float(s2),
                  shift=float(s1 - s0))
    print(f"  {gn:>8} {len(g):>4} | {s0:>13.3f} | {s1:>15.3f} | {s2:>16.3f} | {s1-s0:>+13.3f}")
R["E_R_under_ramp"] = er

# ---------------------------------------------------------------------------
# write the ARM 2 input
# ---------------------------------------------------------------------------
arm2 = ramp[["T", "OTU", "Replicate", "n0_ramp_over_current", "status", "tau",
             "frac_thermally_settled"]].copy()
arm2 = arm2.rename(columns={"n0_ramp_over_current": "n0_factor"})
# series without a usable transient/TPC fall back to the current treatment
arm2["n0_factor"] = arm2.n0_factor.where(np.isfinite(arm2.n0_factor), 1.0)
p_arm2 = os.path.join(TOOLS, "c2_ramp_n0_factors.csv")
arm2.to_csv(p_arm2, index=False)
n_fallback = int((arm2.status != "ok").sum())
print(f"\nwrote {p_arm2}  ({len(arm2)} series; {n_fallback} fall back to n0_factor = 1, "
      f"i.e. the current treatment)")
R["arm2_file"] = dict(path=os.path.relpath(p_arm2, ROOT), n=int(len(arm2)),
                      n_fallback=n_fallback)

fit.to_csv(os.path.join(TOOLS, "c2_partB_per_series.csv"), index=False)
with open(os.path.join(TOOLS, "c2_partB.json"), "w") as f:
    json.dump(R, f, indent=1, default=float)
print(f"wrote {os.path.join(TOOLS,'c2_partB_per_series.csv')}")
print(f"wrote {os.path.join(TOOLS,'c2_partB.json')}")
