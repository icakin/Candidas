"""
known_inconsistencies.py -- C1 PART F.4. Verify or refute each reported
internal inconsistency BY DIRECT TEST, and report the numbers either way.
Fixes nothing.

    cauris_etcgem/.venv/bin/python reports/tools/known_inconsistencies.py

Read-only. Writes reports/tools/known_inconsistencies.json and prints a report.

Tests here are the ones that need only committed CSVs and the R sources.
The two that need a BUILT MODEL (the anchor, and the Clade IV ceiling) live in
reports/tools/anchor_and_ceiling_test.py.
"""
import hashlib, json, os, re, sys
import numpy as np
import pandas as pd

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = os.path.dirname(os.path.dirname(HERE))
STRAIN = os.path.join(ROOT, "cauris_etcgem", "strains", "eci_cauris")
SUPP   = os.path.join(STRAIN, "outputs", "supp_data")
BACKUP = os.path.join(STRAIN, "outputs", "supp_data_original_backup")
FIG3   = os.path.join(STRAIN, "outputs", "figure3_data")
TABLES = os.path.join(ROOT, "results", "tables")
CLADES = ["I", "II", "III", "IV"]
OUT    = os.path.join(HERE, "known_inconsistencies.json")

R = {}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def head(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# =============================================================================
# 3. capacity_isolates.csv: identical to the backup, and pinned to the BACKUP
#    clade kcat_scale on the 120-point capacity grid
# =============================================================================
head("3. capacity_isolates.csv -- which build is it from?")

ci_supp, ci_back = os.path.join(SUPP, "capacity_isolates.csv"), os.path.join(BACKUP, "capacity_isolates.csv")
h1, h2 = sha(ci_supp), sha(ci_back)
print(f"  sha256 supp_data : {h1}")
print(f"  sha256 backup    : {h2}")
print(f"  BYTE-IDENTICAL   : {h1 == h2}")

cap = pd.read_csv(ci_supp)
r2_supp = pd.read_csv(os.path.join(SUPP, "calibration_r2.csv")).set_index("clade")
r2_back = pd.read_csv(os.path.join(BACKUP, "calibration_r2.csv")).set_index("clade")

# KG = np.linspace(0.20, 1.35, 120) in stage_boot
KG   = np.linspace(0.20, 1.35, 120)
STEP = float(KG[1] - KG[0])
print(f"\n  capacity grid: np.linspace(0.20, 1.35, 120) -> step = {STEP:.10f}")
print("  (the prompt quotes 0.0096639; that is this step)\n")

rows = []
print(f"  {'clade':5} {'isolate':14} {'capacity':>12} "
      f"{'steps vs BACKUP kcat':>22} {'steps vs CURRENT kcat':>23}")
for _, r in cap.iterrows():
    cl = r["clade"]
    d_back = (r["capacity"] - float(r2_back.loc[cl, "kcat_scale"])) / STEP
    d_curr = (r["capacity"] - float(r2_supp.loc[cl, "kcat_scale"])) / STEP
    rows.append(dict(clade=cl, isolate=r["isolate"], capacity=float(r["capacity"]),
                     steps_vs_backup=float(d_back), steps_vs_current=float(d_curr),
                     backup_is_integer=bool(abs(d_back - round(d_back)) < 1e-6),
                     current_is_integer=bool(abs(d_curr - round(d_curr)) < 1e-6)))
    print(f"  {cl:5} {r['isolate']:14} {r['capacity']:12.10f} "
          f"{d_back:>18.6f} {'INT' if abs(d_back-round(d_back))<1e-6 else '   '} "
          f"{d_curr:>19.6f} {'INT' if abs(d_curr-round(d_curr))<1e-6 else '   '}")

n_back = sum(x["backup_is_integer"] for x in rows)
n_curr = sum(x["current_is_integer"] for x in rows)
print(f"\n  integer multiples of the grid step away from BACKUP  kcat_scale: {n_back}/12")
print(f"  integer multiples of the grid step away from CURRENT kcat_scale: {n_curr}/12")

R["test3_capacity_isolates_provenance"] = dict(
    sha_supp=h1, sha_backup=h2, byte_identical=(h1 == h2),
    grid_step=STEP, per_isolate=rows,
    n_integer_vs_backup=n_back, n_integer_vs_current=n_curr)


# ---- the downstream numbers that depend on capacity_isolates.csv -------------
def isolate_stats(df):
    """ICC1, F(3,8), p from a one-way ANOVA of capacity by clade (n=3 each)."""
    from scipy import stats
    groups = [df.capacity[df.clade == c].values for c in CLADES]
    k, n = len(groups), len(groups[0])
    grand = np.concatenate(groups).mean()
    ssb = n * sum((g.mean() - grand) ** 2 for g in groups)
    ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
    dfb, dfw = k - 1, k * (n - 1)
    msb, msw = ssb / dfb, ssw / dfw
    F = msb / msw
    p = float(stats.f.sf(F, dfb, dfw))
    icc = (msb - msw) / (msb + (n - 1) * msw)
    return dict(ICC1=float(icc), F=float(F), df=[dfb, dfw], p=p,
                MSb=float(msb), MSw=float(msw))

st = isolate_stats(cap)
print(f"\n  from THIS capacity_isolates.csv:  ICC1 = {st['ICC1']:.4f} "
      f"({100*st['ICC1']:.0f}%)  F(3,8) = {st['F']:.2f}  p = {st['p']:.3e}")
R["test3_downstream_stats_committed"] = st


# =============================================================================
# 4. figure3_data/capacity.csv vs calibration_r2.csv kcat_scale
# =============================================================================
head("4. figure3_data/capacity.csv vs supp_data/calibration_r2.csv kcat_scale")

f3 = pd.read_csv(os.path.join(FIG3, "capacity.csv")).set_index("clade")
rows = []
print(f"  {'clade':5} {'fig3 capacity':>16} {'supp kcat_scale':>17} {'diff':>12} "
      f"{'backup kcat_scale':>19} {'diff vs backup':>16}")
for cl in CLADES:
    a = float(f3.loc[cl, "capacity"])
    b = float(r2_supp.loc[cl, "kcat_scale"])
    c = float(r2_back.loc[cl, "kcat_scale"])
    rows.append(dict(clade=cl, fig3_capacity=a, supp_kcat=b, backup_kcat=c,
                     diff_vs_supp=a - b, diff_vs_backup=a - c,
                     equals_backup_exactly=bool(a == c)))
    print(f"  {cl:5} {a:16.10f} {b:17.10f} {a-b:12.2e} {c:19.10f} {a-c:16.2e}")
print("\n  fig3 capacity == BACKUP kcat_scale exactly (bit-for-bit)? "
      f"{all(x['equals_backup_exactly'] for x in rows)}")
# the r2 column too
r2_match_backup = all(float(f3.loc[cl, "r2"]) == float(r2_back.loc[cl, "r2_calibrated"])
                      for cl in CLADES)
r2_match_supp = all(float(f3.loc[cl, "r2"]) == float(r2_supp.loc[cl, "r2_calibrated"])
                    for cl in CLADES)
print(f"  fig3 r2 column == BACKUP r2_calibrated exactly?  {r2_match_backup}")
print(f"  fig3 r2 column == CURRENT r2_calibrated exactly? {r2_match_supp}")
R["test4_fig3_vs_calibration_r2"] = dict(
    per_clade=rows, r2_matches_backup=r2_match_backup, r2_matches_current=r2_match_supp)


# =============================================================================
# 5. apriori_Topt differs across clades although the a-priori curve is
#    clade-independent by construction
# =============================================================================
head("5. apriori_Topt across clades vs the actual a-priori curve")

for tag, d in (("supp_data", SUPP), ("backup", BACKUP)):
    cur = pd.read_csv(os.path.join(d, "calibration_curves.csv"))
    r2t = pd.read_csv(os.path.join(d, "calibration_r2.csv")).set_index("clade")
    piv = cur.pivot(index="temp_C", columns="clade", values="apriori")
    spread = float((piv.max(axis=1) - piv.min(axis=1)).max())
    ap = piv["I"].to_numpy(float)
    Ts = piv.index.to_numpy(float)
    pk = float(ap.max())
    argmax_T = float(Ts[int(np.argmax(ap))])
    ties = int((ap == pk).sum())
    tie_range = [float(Ts[ap == pk].min()), float(Ts[ap == pk].max())]
    within = {}
    for frac in (1e-6, 1e-3, 1e-2, 5e-2):
        w = Ts[ap >= pk * (1 - frac)]
        within[f"{frac:g}"] = [float(w.min()), float(w.max()), int(w.size)]
    reported = [float(r2t.loc[cl, "apriori_Topt"]) for cl in CLADES]
    print(f"\n  --- {tag} ---")
    print(f"  a-priori curve, max across-clade |difference| at any T : {spread:.3e}")
    print(f"  a-priori peak                                          : {pk:.10f} at {argmax_T} C")
    print(f"  exact float ties at the peak                           : {ties} "
          f"(T = {tie_range[0]} .. {tie_range[1]} C)")
    for k, v in within.items():
        print(f"    within {k:>6} of the peak: {v[0]:5.1f} - {v[1]:5.1f} C ({v[2]} of {len(Ts)} points)")
    print(f"  np.argmax on that array MUST give                       : {argmax_T} C for EVERY clade")
    print(f"  calibration_r2.csv apriori_Topt actually reports        : {reported}")
    print(f"  consistent?                                            : "
          f"{all(x == argmax_T for x in reported)}")
    R[f"test5_apriori_Topt_{tag}"] = dict(
        across_clade_spread=spread, apriori_peak=pk, apriori_argmax_T=argmax_T,
        n_exact_ties=ties, tie_range=tie_range, plateau_within=within,
        reported_apriori_Topt=reported,
        consistent=bool(all(x == argmax_T for x in reported)))


# =============================================================================
# 6. drop_excluded() / EXCLUDE_GROUPS ("glab") is applied only in 08 and 09
# =============================================================================
head("6. EXCLUDE_GROUPS = glab -- where is drop_excluded() actually called?")

scripts_dir = os.path.join(ROOT, "scripts")
calls = {}
for fn in sorted(os.listdir(scripts_dir)):
    if not fn.endswith(".R"):
        continue
    txt = open(os.path.join(scripts_dir, fn), encoding="utf-8", errors="replace").read()
    n_def  = len(re.findall(r"drop_excluded\s*<-\s*function", txt))
    n_call = len(re.findall(r"drop_excluded\s*\(", txt)) - n_def
    if n_call > 0 or n_def > 0:
        calls[fn] = dict(defined=n_def, called=n_call)
        print(f"  {fn:32} defined={n_def}  called={n_call}")

# quantify: which committed tables still contain glabrata?
GLAB_OTUS = [13, 14, 15]
tabs = {}
print("\n  glabrata content of every committed results/tables/*.csv:")
print(f"  {'file':52} {'rows':>7} {'glab rows':>10} {'%':>7}")
for fn in sorted(os.listdir(TABLES)):
    if not fn.endswith(".csv"):
        continue
    p = os.path.join(TABLES, fn)
    try:
        df = pd.read_csv(p)
    except Exception:
        continue
    n = len(df)
    g = 0
    hit = None
    for col in ("Group", "group"):
        if col in df.columns:
            g = int((df[col].astype(str) == "glab").sum()); hit = col; break
    else:
        for col in ("OTU", "otu"):
            if col in df.columns:
                g = int(pd.to_numeric(df[col], errors="coerce").isin(GLAB_OTUS).sum())
                hit = col; break
        else:
            for col in ("otu_name", "Isolate", "isolate"):
                if col in df.columns:
                    g = int(df[col].astype(str).str.contains("glab", case=False).sum())
                    hit = col; break
    if hit is None:
        continue
    tabs[fn] = dict(rows=n, glab_rows=g, key_col=hit,
                    pct=round(100.0 * g / n, 2) if n else None)
    print(f"  {fn:52} {n:7d} {g:10d} {100.0*g/n if n else 0:6.1f}%")

R["test6_exclude_groups"] = dict(drop_excluded_calls=calls, tables=tabs)


# =============================================================================
# 7. capacity_bootstrap_draws.csv pairs capacity and peak by two independent
#    rng.choice calls -> within-clade corr(capacity, peak) should be ~0
# =============================================================================
head("7. capacity_bootstrap_draws.csv -- is (capacity, peak) a real pairing?")

dr = pd.read_csv(os.path.join(SUPP, "capacity_bootstrap_draws.csv"))
per = {}
print(f"  {'clade':6} {'n draws':>8} {'corr(capacity, peak)':>22}")
for cl in CLADES:
    d = dr[dr.clade == cl]
    r = float(np.corrcoef(d.capacity.values, d.peak.values)[0, 1])
    per[cl] = dict(n=int(len(d)), corr=r)
    print(f"  {cl:6} {len(d):8d} {r:22.6f}")
allr = float(np.corrcoef(dr.capacity.values, dr.peak.values)[0, 1])
print(f"\n  pooled across clades (this one is DRIVEN BY the between-clade spread): {allr:.6f}")
print("  The source is:")
print("    for c, pk in zip(*[caps[rng.choice(B, 1500, replace=False)],")
print("                       peaks[rng.choice(B, 1500, replace=False)]]):")
print("  i.e. two INDEPENDENT resamples of the same 4000 bootstrap replicates, so")
print("  draw i's capacity and draw i's peak come from different bootstrap samples.")
R["test7_bootstrap_pairing"] = dict(per_clade=per, pooled_corr=allr)


# =============================================================================
# summary of which files in supp_data/ match the backup byte-for-byte
# =============================================================================
head("cross-check: which supp_data/ files are byte-identical to the backup?")
prov = {}
for fn in sorted(os.listdir(SUPP)):
    pa, pb = os.path.join(SUPP, fn), os.path.join(BACKUP, fn)
    if not (os.path.isfile(pa) and os.path.isfile(pb)):
        continue
    same = sha(pa) == sha(pb)
    prov[fn] = same
    print(f"  {fn:36} {'IDENTICAL to backup' if same else 'differs from backup'}")
R["supp_vs_backup_identical"] = prov

with open(OUT, "w") as f:
    json.dump(R, f, indent=1, default=float)
print(f"\nwrote {OUT}")
