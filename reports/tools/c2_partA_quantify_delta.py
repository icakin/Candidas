"""
c2_partA_quantify_delta.py -- C2 PART A.

What is `delta`, how big is exp(r*delta), and how much of E_R is it?

    cauris_etcgem/.venv/bin/python reports/tools/c2_partA_quantify_delta.py

Read-only. Writes reports/tools/c2_partA.json and prints the tables.

Inputs (all from the C1 regeneration, so provenance is known):
    results_C1/tables/derived_N0_R_results_with_carbon.csv
    results_C1/tables/Oxygen_Trimmed_Series_Metadata.csv

Definitions, taken from the source rather than assumed:
    07_oxygen_fits.R:361   delta_Ninoc_to_N0_min = fit_start_time   (if N0_BACKPROJECT)
    07_oxygen_fits.R:922   N0 = N_inoc * exp(r * delta)
    07_oxygen_fits.R:928   C_tot_O2 = (K/r) * (exp(r*T_end) - 1)
    07_oxygen_fits.R:935   biomass_integral = N0 * (exp(r*T_end) - 1) / r
    07_oxygen_fits.R:941   R_O2_per_cell = C_tot_O2 / biomass_integral
  => R_O2_per_cell = K / N0 EXACTLY; T_end cancels. So
        log R = log K - log N_inoc - r*delta
  and respiration is an explicit DECREASING function of the fitted growth rate.
"""
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T_AMBIENT = 20.0          # assumed lab ambient; stated in the report
K_B = 8.617e-5            # eV/K
T_REF = 293.15
GROUPS = ["Clade1", "Clade2", "Clade3", "Clade4", "glab", "para"]
OUT = os.path.join(ROOT, "reports", "tools", "c2_partA.json")
R = {}


def head(t):
    print("\n" + "=" * 86 + f"\n{t}\n" + "=" * 86)


def q(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if x.size else np.nan


def desc(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not x.size:
        return dict(n=0)
    return dict(n=int(x.size), median=float(np.median(x)),
                q25=q(x, .25), q75=q(x, .75), p95=q(x, .95),
                min=float(x.min()), max=float(x.max()), mean=float(x.mean()))


d = pd.read_csv(os.path.join(ROOT, "results_C1", "tables",
                             "derived_N0_R_results_with_carbon.csv"))
meta = pd.read_csv(os.path.join(ROOT, "results_C1", "tables",
                                "Oxygen_Trimmed_Series_Metadata.csv"))

d["group"] = d.otu_name.str.extract(r"^([A-Za-z0-9]+)_")[0]
d["delta"] = d.delta_Ninoc_to_N0_min
d["rdelta"] = d.r * d.delta
d["expf"] = np.exp(d.rdelta)
d["TK"] = d["T"] + 273.15
d["boltz"] = 1.0 / (K_B * T_REF) - 1.0 / (K_B * d.TK)
# K is the volumetric consumption rate at the window start, mg O2 / L / min.
# It is the N0-FREE quantity: respiration = K / N0.
d["logK"] = np.log(d.K)
d["logR"] = np.log(d.respiration_fgC_h)
# respiration with the back-projection term removed, i.e. N0 = N_inoc
d["logR_nobp"] = d.logR + d.rdelta

use = d[d.keep & np.isfinite(d.r) & (d.r > 0) & np.isfinite(d.K) & (d.K > 0)].copy()
print(f"series used: {len(use)} of {len(d)} derived rows "
      f"({d.group.nunique()} groups, {use['T'].nunique()} temperatures)")
R["n_series"] = int(len(use))


# =============================================================================
# A1. distributions of delta, r, r*delta, exp(r*delta)
# =============================================================================
head("A1. delta, r, r*delta and exp(r*delta) -- overall")
for nm, col in [("delta (min)", "delta"), ("r (1/min)", "r"),
                ("r*delta", "rdelta"), ("exp(r*delta)", "expf")]:
    s = desc(use[col])
    print(f"  {nm:16} n={s['n']:4d}  median={s['median']:10.4f}  "
          f"IQR=[{s['q25']:.4f}, {s['q75']:.4f}]  p95={s['p95']:.4f}  max={s['max']:.4f}")
    R.setdefault("overall", {})[col] = s

head("A2. by TEMPERATURE")
print(f"  {'T':>4} {'n':>4} | {'delta med':>10} {'IQR':>18} | {'r*delta med':>12} "
      f"{'p95':>8} | {'exp med':>8} {'p95':>8} {'max':>9}")
byT = {}
for T, g in use.groupby("T"):
    dd, rd, ef = desc(g.delta), desc(g.rdelta), desc(g.expf)
    byT[int(T)] = dict(n=dd["n"], delta=dd, rdelta=rd, expf=ef)
    print(f"  {T:>4} {dd['n']:>4} | {dd['median']:>10.1f} "
          f"[{dd['q25']:>6.1f},{dd['q75']:>6.1f}] | {rd['median']:>12.3f} {rd['p95']:>8.3f} "
          f"| {ef['median']:>8.2f} {ef['p95']:>8.2f} {ef['max']:>9.2f}")
R["by_temperature"] = byT

head("A3. by TAXON")
print(f"  {'group':>8} {'n':>4} | {'delta med':>10} {'IQR':>18} | {'r*delta med':>12} "
      f"{'p95':>8} | {'exp med':>8} {'p95':>8} {'max':>9}")
byG = {}
for gname in GROUPS:
    g = use[use.group == gname]
    if not len(g):
        continue
    dd, rd, ef = desc(g.delta), desc(g.rdelta), desc(g.expf)
    byG[gname] = dict(n=dd["n"], delta=dd, rdelta=rd, expf=ef)
    print(f"  {gname:>8} {dd['n']:>4} | {dd['median']:>10.1f} "
          f"[{dd['q25']:>6.1f},{dd['q75']:>6.1f}] | {rd['median']:>12.3f} {rd['p95']:>8.3f} "
          f"| {ef['median']:>8.2f} {ef['p95']:>8.2f} {ef['max']:>9.2f}")
R["by_taxon"] = byG

head("A4. values at 22, 34 and 44 C, per taxon  (median [IQR])")
print(f"  {'group':>8} | " + " | ".join(f"{'T=%d' % T:^46}" for T in (22, 34, 44)))
print(f"  {'':>8} | " + " | ".join(f"{'delta':>12} {'r*delta':>14} {'exp(r*delta)':>16}"
                                   for _ in (22, 34, 44)))
key = {}
for gname in GROUPS:
    cells = []
    for T in (22, 34, 44):
        g = use[(use.group == gname) & (use["T"] == T)]
        if not len(g):
            cells.append(f"{'-':>12} {'-':>14} {'-':>16}"); continue
        cells.append(f"{np.median(g.delta):>12.1f} {np.median(g.rdelta):>14.3f} "
                     f"{np.median(g.expf):>16.2f}")
        key[f"{gname}@{T}"] = dict(n=int(len(g)),
                                   delta_med=float(np.median(g.delta)),
                                   rdelta_med=float(np.median(g.rdelta)),
                                   expf_med=float(np.median(g.expf)),
                                   expf_max=float(np.max(g.expf)))
    print(f"  {gname:>8} | " + " | ".join(cells))
R["at_22_34_44"] = key


# =============================================================================
# A5. does delta depend on temperature, taxon, or r?
# =============================================================================
head("A5. what does delta depend on?")

# (a) delta vs T -- linear, and Spearman (monotone but not necessarily linear)
sl, ic, rv, pv, se = stats.linregress(use["T"], use.delta)
rho, prho = stats.spearmanr(use["T"], use.delta)
print(f"  delta ~ T                       slope = {sl:+.4f} min/degC "
      f"(SE {se:.4f}, p = {pv:.3g}), R2 = {rv**2:.4f}; Spearman rho = {rho:+.3f} (p = {prho:.3g})")
R["delta_vs_T"] = dict(slope=float(sl), se=float(se), p=float(pv), r2=float(rv**2),
                       spearman_rho=float(rho), spearman_p=float(prho))

# (b) delta vs |T - T_ambient| -- the instrumental prediction
use["dT"] = np.abs(use["T"] - T_AMBIENT)
sl2, ic2, rv2, pv2, se2 = stats.linregress(use.dT, use.delta)
print(f"  delta ~ |T - {T_AMBIENT:.0f}C|                slope = {sl2:+.4f} min/degC "
      f"(SE {se2:.4f}, p = {pv2:.3g}), R2 = {rv2**2:.4f}, intercept = {ic2:.1f} min")
R["delta_vs_absdT"] = dict(T_ambient=T_AMBIENT, slope=float(sl2), se=float(se2),
                           p=float(pv2), r2=float(rv2**2), intercept=float(ic2))

# (c) delta vs r
slr, icr, rvr, pvr, ser = stats.linregress(use.r, use.delta)
rhor, prhor = stats.spearmanr(use.r, use.delta)
print(f"  delta ~ r                       slope = {slr:+.1f} min per (1/min) "
      f"(p = {pvr:.3g}), R2 = {rvr**2:.4f}; Spearman rho = {rhor:+.3f} (p = {prhor:.3g})")
R["delta_vs_r"] = dict(slope=float(slr), p=float(pvr), r2=float(rvr**2),
                       spearman_rho=float(rhor), spearman_p=float(prhor))

# (d) taxon effect -- one-way ANOVA and the per-group means
groups_present = [gn for gn in GROUPS if (use.group == gn).any()]
arrs = [use.delta[use.group == gn].values for gn in groups_present]
F, pF = stats.f_oneway(*arrs)
print(f"  delta ~ taxon                   ANOVA F({len(arrs)-1},{len(use)-len(arrs)}) "
      f"= {F:.2f}, p = {pF:.3g}")
grand = use.delta.mean()
print("     per-taxon mean delta (min), and residual after removing the grand mean:")
tax_eff = {}
for gn in groups_present:
    v = use.delta[use.group == gn]
    tax_eff[gn] = dict(n=int(len(v)), mean=float(v.mean()), median=float(v.median()),
                       residual=float(v.mean() - grand))
    print(f"       {gn:>8}  n={len(v):>4}  mean={v.mean():7.1f}  median={v.median():7.1f}  "
          f"residual={v.mean()-grand:+7.1f}")
R["delta_by_taxon_anova"] = dict(F=float(F), p=float(pF), grand_mean=float(grand),
                                 per_taxon=tax_eff)

# (e) non-glabrata only: is there still a temperature trend?
ng = use[use.group != "glab"]
sl3, ic3, rv3, pv3, se3 = stats.linregress(ng["T"], ng.delta)
print(f"\n  EXCLUDING glabrata (n={len(ng)}):")
print(f"    delta ~ T                     slope = {sl3:+.4f} min/degC (p = {pv3:.3g}), "
      f"R2 = {rv3**2:.4f}")
arrs2 = [use.delta[use.group == gn].values for gn in groups_present if gn != "glab"]
F2, pF2 = stats.f_oneway(*arrs2)
print(f"    delta ~ taxon                 ANOVA F = {F2:.2f}, p = {pF2:.3g}")
R["delta_excl_glab"] = dict(n=int(len(ng)), slope_vs_T=float(sl3), p_vs_T=float(pv3),
                            r2_vs_T=float(rv3**2), anova_F=float(F2), anova_p=float(pF2))

# (f) how does delta relate to the automatic O2 peak? delta is the FIT-WINDOW
#     start, which the trim-selector app may have moved by hand.
mm = meta[["T", "OTU", "Replicate", "peak_time", "main_run_start_time"]].copy()
j = use.merge(mm, on=["T", "OTU", "Replicate"], how="left")
same = np.isclose(j.delta, j.peak_time, atol=1e-6)
print(f"\n  delta vs the automatic O2 peak: identical in {same.sum()} of {len(j)} series "
      f"({100*same.mean():.1f}%)")
dif = (j.delta - j.peak_time)
print(f"    delta - peak_time: median {np.nanmedian(dif):+.1f} min, "
      f"IQR [{q(dif,.25):+.1f}, {q(dif,.75):+.1f}], "
      f"range [{np.nanmin(dif):+.1f}, {np.nanmax(dif):+.1f}]")
print("    -> delta is the HAND-CHOSEN fit-window start (04_trim_selector.R), not the")
print("       automatic peak. Both are 'time from liquid-in to the start of the fit'.")
R["delta_vs_peak"] = dict(n=int(len(j)), n_identical=int(same.sum()),
                          frac_identical=float(same.mean()),
                          med_diff=float(np.nanmedian(dif)),
                          iqr=[q(dif, .25), q(dif, .75)],
                          min=float(np.nanmin(dif)), max=float(np.nanmax(dif)))


# =============================================================================
# A6. how much of E_R is the correction?
# =============================================================================
head("A6. Boltzmann-Arrhenius slopes: E_K (no N0) vs E_R (shipped) vs E_R (term removed)")


def ols_slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return dict(E=np.nan, se=np.nan, n=int(m.sum()), r2=np.nan)
    s, i, r_, p_, se_ = stats.linregress(x[m], y[m])
    return dict(E=float(s), se=float(se_), n=int(m.sum()), r2=float(r_ ** 2), p=float(p_))


def slopes(df):
    return dict(
        E_K=ols_slope(df.boltz.values, df.logK.values),
        E_R_shipped=ols_slope(df.boltz.values, df.logR.values),
        E_R_nobp=ols_slope(df.boltz.values, df.logR_nobp.values),
    )


for label, sub in [("FULL RANGE 22-44 C", use),
                   ("RISING LIMB T <= 34 C", use[use["T"] <= 34])]:
    print(f"\n  --- {label} ---")
    print(f"  {'group':>8} {'n':>4} | {'E_K (eV)':>16} | {'E_R shipped':>16} | "
          f"{'E_R no-backproj':>17} | {'artefact E_K-E_R':>17}")
    tab = {}
    for gn in groups_present:
        g = sub[sub.group == gn]
        s = slopes(g)
        art = s["E_K"]["E"] - s["E_R_shipped"]["E"]
        tab[gn] = dict(**s, artefact_EK_minus_ER=float(art))
        print(f"  {gn:>8} {len(g):>4} | {s['E_K']['E']:>7.3f} +- {s['E_K']['se']:<6.3f} | "
              f"{s['E_R_shipped']['E']:>7.3f} +- {s['E_R_shipped']['se']:<6.3f} | "
              f"{s['E_R_nobp']['E']:>8.3f} +- {s['E_R_nobp']['se']:<6.3f} | {art:>17.3f}")
    # ordering
    ship = {gn: tab[gn]["E_R_shipped"]["E"] for gn in tab}
    nobp = {gn: tab[gn]["E_R_nobp"]["E"] for gn in tab}
    ok = {gn: tab[gn]["E_K"]["E"] for gn in tab}
    o_ship = sorted(ship, key=ship.get)
    o_nobp = sorted(nobp, key=nobp.get)
    o_k = sorted(ok, key=ok.get)
    print(f"\n    ordering of E_R, shipped        : {' < '.join(o_ship)}")
    print(f"    ordering of E_R, term removed   : {' < '.join(o_nobp)}")
    print(f"    ordering of E_K (no N0 at all)  : {' < '.join(o_k)}")
    print(f"    ORDERING SURVIVES REMOVAL?      : {o_ship == o_nobp}")
    tau_s, tau_p = stats.spearmanr([ship[g] for g in groups_present],
                                   [nobp[g] for g in groups_present])
    print(f"    Spearman(E_R shipped, E_R removed) across taxa = {tau_s:+.3f} (p = {tau_p:.3g})")
    sp_ship = max(ship.values()) - min(ship.values())
    sp_nobp = max(nobp.values()) - min(nobp.values())
    print(f"    between-taxon SPREAD of E_R: shipped {sp_ship:.3f} eV -> removed {sp_nobp:.3f} eV "
          f"({100*(1-sp_nobp/sp_ship):+.0f}% change)")
    R[f"slopes_{'full' if 'FULL' in label else 'rising'}"] = dict(
        per_taxon=tab, order_shipped=o_ship, order_nobp=o_nobp, order_EK=o_k,
        ordering_survives=bool(o_ship == o_nobp),
        spearman_shipped_vs_removed=float(tau_s),
        spread_shipped=float(sp_ship), spread_nobp=float(sp_nobp))

# the arithmetic identity, stated as a check rather than assumed
head("A7. identity check: is log R exactly log K - log N_inoc - r*delta?")
lhs = use.logR.values
rhs = (use.logK - np.log(use.N_inoculation_cells_per_L) - use.rdelta).values
# respiration_fgC_h carries the unit conversion MG_TO_FG * O2_TO_C_MASS * RQ * MIN_TO_H,
# a single constant, so the two differ by exactly that constant if the identity holds.
diff = lhs - rhs
print(f"  log R - [log K - log N_inoc - r*delta] : mean {diff.mean():.9f}, "
      f"sd {diff.std():.3e}, range [{diff.min():.9f}, {diff.max():.9f}]")
print(f"  expected constant = log(MG_TO_FG * O2_TO_C_MASS * RQ * MIN_TO_H) = "
      f"{np.log(1e12 * (12.011/31.998) * 1 * 60):.9f}")
print(f"  -> the identity holds to {diff.std():.1e}; respiration IS K/N0 and nothing else.")
R["identity_check"] = dict(mean=float(diff.mean()), sd=float(diff.std()),
                           expected=float(np.log(1e12 * (12.011 / 31.998) * 60)))

with open(OUT, "w") as f:
    json.dump(R, f, indent=1, default=float)
print(f"\nwrote {OUT}")
