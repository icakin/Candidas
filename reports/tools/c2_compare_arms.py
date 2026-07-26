"""
c2_compare_arms.py -- C2 PARTs C0 and D.

    cauris_etcgem/.venv/bin/python reports/tools/c2_compare_arms.py

Read-only. Writes reports/tools/c2_arms.json and prints every table the
sensitivity report needs.

  C0  the reproducibility floor: results_C2_current vs results_C1 (same
      treatment, different run). Anything smaller than the floor is not
      evidence.
  D   the three arms side by side, plus the growth-side invariance check.
"""
import glob, json, os, re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS = os.path.join(ROOT, "reports", "tools")
ARMS = ["current", "ramp", "nobp"]
ARM_DIR = {a: os.path.join(ROOT, f"results_C2_{a}") for a in ARMS}
TAXA = ["Clade1", "Clade2", "Clade3", "Clade4", "para"]
CLADES = ["I", "II", "III", "IV"]
R = {}


def head(t):
    print("\n" + "=" * 92 + f"\n{t}\n" + "=" * 92)


def rd(arm, f):
    p = os.path.join(ARM_DIR[arm], "tables", f)
    return pd.read_csv(p) if os.path.exists(p) else None


# =============================================================================
# C0. the reproducibility floor
# =============================================================================
head("C0. reproducibility floor -- results_C2_current vs results_C1 (same treatment)")

a_dir = os.path.join(ROOT, "results_C1", "tables")
b_dir = os.path.join(ARM_DIR["current"], "tables")
floor = {}
print(f"  {'file':46} {'rows':>6} {'max abs diff':>14} {'max rel diff':>14}  verdict")
for fn in sorted(os.listdir(b_dir)):
    if not fn.endswith(".csv"):
        continue
    pa = os.path.join(a_dir, fn)
    if not os.path.exists(pa):
        continue
    try:
        A = pd.read_csv(pa); B = pd.read_csv(os.path.join(b_dir, fn))
    except Exception:
        continue
    if len(A) != len(B):
        print(f"  {fn:46} {'-':>6} {'-':>14} {'-':>14}  ROW COUNT {len(A)} vs {len(B)}")
        floor[fn] = dict(status="row_count", n_a=len(A), n_b=len(B)); continue
    kc = [c for c in A.columns if not pd.api.types.is_numeric_dtype(A[c]) and c in B.columns]
    if kc:
        A = A.sort_values(kc).reset_index(drop=True); B = B.sort_values(kc).reset_index(drop=True)
    nc = [c for c in A.columns if pd.api.types.is_numeric_dtype(A[c]) and c in B.columns
          and pd.api.types.is_numeric_dtype(B[c])]
    if not nc:
        continue
    ma = mr = 0.0
    for c in nc:
        x = pd.to_numeric(A[c], errors="coerce").to_numpy(float)
        y = pd.to_numeric(B[c], errors="coerce").to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        if not m.any():
            continue
        d = np.abs(x[m] - y[m])
        den = np.maximum.reduce([np.abs(x[m]), np.abs(y[m]), np.full(m.sum(), 1e-300)])
        ma = max(ma, float(d.max())); mr = max(mr, float((d / den).max()))
    v = "IDENTICAL" if ma < 1e-10 else ("EQUIVALENT" if ma < 1e-6 else "DIFFERS")
    floor[fn] = dict(max_abs=ma, max_rel=mr, verdict=v, rows=int(len(A)))
    if ma > 0:
        print(f"  {fn:46} {len(A):>6} {ma:>14.3e} {mr:>14.3e}  {v}")
n_ident = sum(1 for v in floor.values() if v.get("verdict") == "IDENTICAL")
print(f"\n  {n_ident} of {len(floor)} comparable files are IDENTICAL (max |diff| < 1e-10).")
print("  Files that differ define the floor: anything smaller than this is not evidence.")
R["floor"] = floor

# what did SEEDING do?  C1 measured the unseeded run-to-run floor for
# carbon_tax_curves.csv at 0.205 absolute / 2.21% relative.
print("\n  Seeding check. C1 measured the UNSEEDED run-to-run floor by running the")
print("  whole pipeline twice: 39 of 40 tables byte-identical, the sole exception")
print("  carbon_tax_curves.csv at 0.205 absolute / 2.21% relative. Here, with")
print("  CANDIDAS_SEED set:")
for fn in ("carbon_tax_curves.csv", "carbon_tax_group.csv", "carbon_tax_isolate.csv",
           "fig_values.csv", "fig_contrasts.csv"):
    v = floor.get(fn)
    if v and "max_abs" in v:
        print(f"      {fn:32} max abs {v['max_abs']:.3e}  max rel {v['max_rel']:.3e}  {v['verdict']}")


# =============================================================================
# D. the three arms
# =============================================================================
head("D1. per-taxon quantities, all three arms")

got = [a for a in ARMS if rd(a, "fig_values.csv") is not None]
print(f"  arms available: {', '.join(got)}")
fv = {a: rd(a, "fig_values.csv").set_index("Group") for a in got}

QUANT = [("E_G", "E", "E_lo", "E_hi"),
         ("E_R", "ER", "ER_lo", "ER_hi"),
         ("dE", "dE", "dE_lo", "dE_hi"),
         ("Topt_CUE", "Tcue", "Tcue_lo", "Tcue_hi"),
         ("tax37", "tax37", "tax37_lo", "tax37_hi"),
         ("tax40", "tax40", "tax40_lo", "tax40_hi")]

per_taxon = {}
for label, c, clo, chi in QUANT:
    print(f"\n  --- {label} ---")
    print(f"  {'taxon':>8} | " + " | ".join(f"{a:^30}" for a in got))
    for g in TAXA:
        cells = []
        for a in got:
            v, lo, hi = (float(fv[a].loc[g, c]), float(fv[a].loc[g, clo]),
                         float(fv[a].loc[g, chi]))
            cells.append(f"{v:>8.4f} [{lo:>8.4f},{hi:>8.4f}]")
            per_taxon.setdefault(g, {}).setdefault(label, {})[a] = dict(v=v, lo=lo, hi=hi)
        print(f"  {g:>8} | " + " | ".join(cells))
R["per_taxon"] = per_taxon

print("\n  dE credible (95% CrI excludes zero)?")
for g in TAXA:
    row = " | ".join(f"{a}: {str(fv[a].loc[g,'dE_credible'])}" for a in got)
    print(f"      {g:>8}  {row}")
R["dE_credible"] = {g: {a: str(fv[a].loc[g, "dE_credible"]) for a in got} for g in TAXA}

print("\n  P(T_opt(CUE) < 37 C):")
for g in TAXA:
    print(f"      {g:>8}  " + " | ".join(f"{a}: {float(fv[a].loc[g,'P_below_body']):.4f}"
                                         for a in got))
R["P_Topt_below_body"] = {g: {a: float(fv[a].loc[g, "P_below_body"]) for a in got}
                          for g in TAXA}

# ---- contrasts ------------------------------------------------------------
head("D2. the 10 pairwise fever contrasts, per arm")
fc = {a: rd(a, "fig_contrasts.csv") for a in got}
print(f"  {'a':>8} {'b':>8} | " + " | ".join(f"{a:^32}" for a in got))
cons = {}
base = fc[got[0]]
for _, r0 in base.iterrows():
    cells = []
    for a in got:
        m = fc[a][(fc[a].a == r0.a) & (fc[a].b == r0.b)]
        if not len(m):
            cells.append(f"{'-':^32}"); continue
        m = m.iloc[0]
        star = "*" if str(m.credible).upper() == "TRUE" else " "
        cells.append(f"{float(m.ratio):>7.4f} [{float(m.lo):>6.4f},{float(m.hi):>6.4f}]{star}")
        cons.setdefault(f"{r0.a}_vs_{r0.b}", {})[a] = dict(
            ratio=float(m.ratio), lo=float(m.lo), hi=float(m.hi),
            credible=str(m.credible))
    print(f"  {r0.a:>8} {r0.b:>8} | " + " | ".join(cells))
print("\n  resolved (95% CrI excludes 1):  " +
      "   ".join(f"{a}: {int((fc[a].credible.astype(str).str.upper()=='TRUE').sum())}/10"
                 for a in got))
R["contrasts"] = cons
R["contrasts_resolved"] = {a: int((fc[a].credible.astype(str).str.upper() == "TRUE").sum())
                           for a in got}

# ---- Spearman rho, from each arm's own log --------------------------------
head("D3. Spearman rho (environmental optimum vs fever cost) -- parsed from the arm logs")
rho = {}
for a in got:
    logs = sorted(glob.glob(os.path.join(ROOT, "logs", f"c2_arm_{a}_*.log")))
    val = None
    for lg in reversed(logs):
        txt = open(lg, errors="replace").read()
        m = re.findall(r"Spearman\(Topt, fever cost\): rho = ([-\d.]+) \[([-\d.]+), ([-\d.]+)\], "
                       r"P\(rho<0\) = ([\d.]+)", txt)
        if m:
            val = dict(rho=float(m[-1][0]), lo=float(m[-1][1]), hi=float(m[-1][2]),
                       P_neg=float(m[-1][3]), log=os.path.basename(lg))
            break
    rho[a] = val
    if val:
        print(f"  {a:>8}  rho = {val['rho']:+.2f} [{val['lo']:+.2f}, {val['hi']:+.2f}], "
              f"P(rho<0) = {val['P_neg']:.3f}")
    else:
        print(f"  {a:>8}  not found in logs")
print("\n  BASELINE: the committed brms fits give rho = -1.00 [-1.00, -0.70] (C1, twice).")
print("  The MANUSCRIPT states rho = -0.90 [-1.00, -0.30]; that mismatch is C1's, not C2's,")
print("  and no fourth number is introduced here.")
R["spearman"] = rho

# ---- CUE -------------------------------------------------------------------
head("D4. CUE by clade")
print(f"  {'taxon':>8} | " + " | ".join(f"{'%s: median  peak' % a:^28}" for a in got))
cue = {}
for g in TAXA:
    cells = []
    for a in got:
        d = rd(a, "bayes_cue_by_clade.csv")
        if d is None:
            cells.append(f"{'-':^28}"); continue
        gg = d[d.Group == g] if "Group" in d.columns else d[d.iloc[:, 0] == g]
        vc = [c for c in ("med", "median", "cue", "CUE") if c in gg.columns]
        col = vc[0] if vc else gg.select_dtypes("number").columns[-1]
        cells.append(f"{float(gg[col].median()):>13.5f} {float(gg[col].max()):>13.5f}")
        cue.setdefault(g, {})[a] = dict(median=float(gg[col].median()),
                                        peak=float(gg[col].max()), column=col)
    print(f"  {g:>8} | " + " | ".join(cells))
R["cue"] = cue

# ---- Fig 4d/4e capacity correlations --------------------------------------
head("D5. Fig 4d / 4e capacity correlations (capacity held FIXED at supp_data_C1)")
cap = pd.read_csv(os.path.join(ROOT, "cauris_etcgem", "strains", "eci_cauris",
                               "outputs", "supp_data_C1", "capacity_isolates.csv"))
print(f"  {'arm':>8} | {'r(capacity, fever tax 40C)':>28} | {'r(capacity, respiration)':>26}")
corr = {}
for a in got:
    tax = rd(a, "carbon_tax_isolate.csv"); der = rd(a, "derived_N0_R_results_with_carbon.csv")
    if tax is None or der is None:
        continue
    fev = tax[tax["T_C"] == 40][["Isolate", "tax"]].rename(
        columns={"Isolate": "otu_name", "tax": "fever_tax40"})
    db = der[(der["keep"] == True) & (der["T"] >= 32) & (der["T"] <= 36)]
    resp = (db.groupby("otu_name")["respiration_fgC_h"].mean().reset_index()
              .rename(columns={"respiration_fgC_h": "resp_opt"}))
    iso = cap.merge(fev, left_on="isolate", right_on="otu_name").merge(
        resp, left_on="isolate", right_on="otu_name")
    r1 = float(np.corrcoef(iso.capacity, iso.fever_tax40)[0, 1])
    p1 = float(stats.pearsonr(iso.capacity, iso.fever_tax40)[1])
    r2 = float(np.corrcoef(iso.capacity, iso.resp_opt)[0, 1])
    p2 = float(stats.pearsonr(iso.capacity, iso.resp_opt)[1])
    corr[a] = dict(n=int(len(iso)), r_tax=r1, p_tax=p1, r_resp=r2, p_resp=p2)
    print(f"  {a:>8} | {r1:>+13.4f} (p={p1:>8.2e}) | {r2:>+11.4f} (p={p2:>8.2e})")
print("\n  The capacity axis is identical in all three arms by construction; only the")
print("  economics axis moves. C1 established the etc-GEM itself is not reproducible,")
print("  so the capacity values remain unresolved pending that re-run.")
R["fig4de"] = corr


# =============================================================================
# D6. INVARIANCE CHECK -- growth must not move
# =============================================================================
head("D6. INVARIANCE CHECK: N0 touches only respiration, so growth must not move")

base_arm = "current"
inv = {}
print(f"  {'quantity':>34} | " + " | ".join(f"{a:^24}" for a in got if a != base_arm))
for col in ("r", "K", "growth_fgC_h", "growth_C_per_C_h"):
    A = rd(base_arm, "derived_N0_R_results_with_carbon.csv")
    cells = []
    for a in got:
        if a == base_arm:
            continue
        B = rd(a, "derived_N0_R_results_with_carbon.csv")
        m = A[["T", "OTU", "Replicate", col]].merge(
            B[["T", "OTU", "Replicate", col]], on=["T", "OTU", "Replicate"],
            suffixes=("_a", "_b"))
        d = np.abs(m[col + "_a"] - m[col + "_b"])
        mx = float(np.nanmax(d))
        cells.append(f"max |diff| = {mx:.3e}")
        inv.setdefault(col, {})[a] = mx
    print(f"  {col:>34} | " + " | ".join(cells))

for f, cols in (("bayes_growth_ss_summary.csv", ("mean", "sd", "q2.5", "q50", "q97.5")),
                ("sharpe_schoolfield_growth_fgC_h_coefs.csv", ("Estimate",)),
                ("arrhenius_growth_fgC_h_coefs.csv", ("alpha", "E"))):
    A = rd(base_arm, f)
    if A is None:
        continue
    cells = []
    for a in got:
        if a == base_arm:
            continue
        B = rd(a, f)
        mx = 0.0
        for c in cols:
            if c in A.columns and c in B.columns and len(A) == len(B):
                mx = max(mx, float(np.nanmax(np.abs(pd.to_numeric(A[c], errors="coerce")
                                                    - pd.to_numeric(B[c], errors="coerce")))))
        cells.append(f"max |diff| = {mx:.3e}")
        inv.setdefault(f, {})[a] = mx
    print(f"  {f:>34} | " + " | ".join(cells))

# E_G and the growth T_opt straight out of fig_values
for label, c in (("E_G (fig_values)", "E"), ("growth Topt (fig_values)", "Topt")):
    cells = []
    for a in got:
        if a == base_arm:
            continue
        mx = float(np.nanmax(np.abs(fv[base_arm][c].astype(float) - fv[a][c].astype(float))))
        cells.append(f"max |diff| = {mx:.3e}")
        inv.setdefault(label, {})[a] = mx
    print(f"  {label:>34} | " + " | ".join(cells))

worst = max((v for d in inv.values() for v in d.values()), default=0.0)
print(f"\n  worst growth-side disagreement across all arms and quantities: {worst:.3e}")
print("  PASS" if worst < 1e-8 else "  NOT BIT-IDENTICAL -- see the report")
R["invariance"] = dict(per_quantity=inv, worst=float(worst), pass_=bool(worst < 1e-8))

with open(os.path.join(TOOLS, "c2_arms.json"), "w") as f:
    json.dump(R, f, indent=1, default=float)
print(f"\nwrote {os.path.join(TOOLS, 'c2_arms.json')}")
