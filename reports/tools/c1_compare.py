"""
c1_compare.py -- row-by-row numeric comparison of two output trees.

    cauris_etcgem/.venv/bin/python reports/tools/c1_compare.py \
        --a results --b results_C1 --label "results vs results_C1" \
        --out reports/tools/compare_results.json

Read-only on the inputs. Writes one JSON (and a markdown table on stdout).

For every CSV present in BOTH trees it:
  * aligns rows on the non-numeric (key) columns when they exist, else on row order;
  * for every shared numeric column reports max |a-b| and max |a-b|/max(|a|,|b|,eps);
  * looks for a credible/confidence interval (lo/hi, l-95% CI/u-95% CI, Q2.5/Q97.5)
    and, when present, tests whether EVERY interval in A overlaps its partner in B;
  * classifies the file:
        IDENTICAL                 max abs diff < 1e-10
        NUMERICALLY EQUIVALENT    max abs diff < 1e-6
        STOCHASTIC-BUT-CONSISTENT differs, and every 95% interval overlaps
        DIVERGENT                 differs, and at least one interval does not overlap
        DIVERGENT (no interval)   differs, and the file carries no interval to test
        NOT REPRODUCED            present in A, absent (or unreadable) in B
        SHAPE MISMATCH            present in both but not row-alignable
"""
import argparse, json, math, os, sys
import numpy as np
import pandas as pd

IDENTICAL_TOL = 1e-10
EQUIV_TOL     = 1e-6

# interval column pairs, in priority order
INTERVAL_PAIRS = [
    ("lo", "hi"), ("cap_lo", "cap_hi"), ("peak_lo", "peak_hi"),
    ("r_lo", "r_hi"), ("q2.5", "q97.5"), ("Q2.5", "Q97.5"),
    ("l-95% CI", "u-95% CI"), ("lower", "upper"), ("conf.low", "conf.high"),
    ("ci_lo", "ci_hi"),
]


def numeric_cols(df):
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]


def key_cols(df):
    return [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]


def align(a, b):
    """Return (a, b) aligned row-for-row, plus how it was done."""
    ka, kb = key_cols(a), key_cols(b)
    shared_keys = [c for c in ka if c in kb]
    if shared_keys:
        try:
            aa = a.sort_values(shared_keys).reset_index(drop=True)
            bb = b.sort_values(shared_keys).reset_index(drop=True)
            if len(aa) == len(bb) and all(
                    (aa[c].astype(str).values == bb[c].astype(str).values).all()
                    for c in shared_keys):
                return aa, bb, "keys: " + ",".join(shared_keys)
        except Exception:
            pass
    if len(a) == len(b):
        return a.reset_index(drop=True), b.reset_index(drop=True), "row order"
    return None, None, f"row counts differ ({len(a)} vs {len(b)})"


def find_intervals(cols):
    """Every (lo, hi) column pair in the file, not just the first.

    Covers the fixed vocabulary in INTERVAL_PAIRS *and* the per-quantity
    convention this project uses in fig_values.csv / bayes_clade_params.csv,
    where each estimate X carries X_lo and X_hi (and sometimes Xlo/Xhi).
    Classifying a file as DIVERGENT because the detector only knew about a
    bare `lo`/`hi` pair would be wrong.
    """
    out, seen = [], set()
    low = {c.lower(): c for c in cols}
    for lo, hi in INTERVAL_PAIRS:
        if lo.lower() in low and hi.lower() in low:
            p = (low[lo.lower()], low[hi.lower()])
            if p not in seen:
                seen.add(p); out.append(p)
    for c in cols:
        for suf in ("_lo", "lo"):
            if c.lower().endswith(suf):
                stem = c[: len(c) - len(suf)]
                for hsuf in ("_hi", "hi"):
                    cand = stem + hsuf
                    match = next((x for x in cols if x.lower() == cand.lower()), None)
                    if match:
                        p = (c, match)
                        if p not in seen:
                            seen.add(p); out.append(p)
                        break
                break
    return out


def compare_file(pa, pb):
    rec = {"file": os.path.basename(pa)}
    try:
        a = pd.read_csv(pa)
    except Exception as e:
        return {**rec, "verdict": "UNREADABLE (A)", "detail": str(e)}
    if not os.path.exists(pb):
        return {**rec, "verdict": "NOT REPRODUCED",
                "detail": "absent from B", "rows_a": len(a)}
    try:
        b = pd.read_csv(pb)
    except Exception as e:
        return {**rec, "verdict": "NOT REPRODUCED", "detail": f"unreadable in B: {e}"}

    aa, bb, how = align(a, b)
    if aa is None:
        return {**rec, "verdict": "SHAPE MISMATCH", "detail": how,
                "rows_a": len(a), "rows_b": len(b)}

    ncols = [c for c in numeric_cols(aa) if c in numeric_cols(bb)]
    if not ncols:
        return {**rec, "verdict": "NO NUMERIC COLUMNS", "detail": how,
                "rows_a": len(aa), "rows_b": len(bb)}

    per_col, max_abs, max_rel, worst_col = {}, 0.0, 0.0, None
    for c in ncols:
        x = pd.to_numeric(aa[c], errors="coerce").to_numpy(float)
        y = pd.to_numeric(bb[c], errors="coerce").to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() == 0:
            per_col[c] = {"max_abs": None, "max_rel": None, "n": 0}
            continue
        d = np.abs(x[m] - y[m])
        den = np.maximum.reduce([np.abs(x[m]), np.abs(y[m]),
                                 np.full(m.sum(), 1e-300)])
        r = d / den
        per_col[c] = {"max_abs": float(d.max()), "max_rel": float(r.max()),
                      "n": int(m.sum()),
                      "n_finite_mismatch": int((np.isfinite(x) != np.isfinite(y)).sum())}
        if d.max() > max_abs:
            max_abs, worst_col = float(d.max()), c
        max_rel = max(max_rel, float(r.max()))

    pairs = find_intervals(ncols)
    interval = None
    if pairs:
        per_pair, tot, tot_ov = [], 0, 0
        for lo, hi in pairs:
            alo = pd.to_numeric(aa[lo], errors="coerce").to_numpy(float)
            ahi = pd.to_numeric(aa[hi], errors="coerce").to_numpy(float)
            blo = pd.to_numeric(bb[lo], errors="coerce").to_numpy(float)
            bhi = pd.to_numeric(bb[hi], errors="coerce").to_numpy(float)
            m = np.isfinite(alo) & np.isfinite(ahi) & np.isfinite(blo) & np.isfinite(bhi)
            ov = (alo[m] <= bhi[m]) & (blo[m] <= ahi[m])
            per_pair.append({"cols": [lo, hi], "n": int(m.sum()),
                             "n_overlap": int(ov.sum()),
                             "all_overlap": bool(ov.all()) if m.sum() else None})
            tot += int(m.sum()); tot_ov += int(ov.sum())
        interval = {"pairs": per_pair, "n_pairs": len(pairs),
                    "n": tot, "n_overlap": tot_ov,
                    "all_overlap": bool(tot_ov == tot) if tot else None}

    if max_abs < IDENTICAL_TOL:
        verdict = "IDENTICAL"
    elif max_abs < EQUIV_TOL:
        verdict = "NUMERICALLY EQUIVALENT"
    elif interval and interval["all_overlap"]:
        verdict = "STOCHASTIC-BUT-CONSISTENT"
    elif interval:
        verdict = "DIVERGENT"
    else:
        verdict = "DIVERGENT (no interval in file)"

    return {**rec, "verdict": verdict, "aligned_by": how,
            "rows": len(aa), "n_numeric_cols": len(ncols),
            "max_abs_diff": max_abs, "max_rel_diff": max_rel,
            "worst_column": worst_col, "interval": interval,
            "per_column": per_col}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="reference tree (committed)")
    ap.add_argument("--b", required=True, help="regenerated tree")
    ap.add_argument("--label", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--recurse", action="store_true",
                    help="walk subdirectories too")
    args = ap.parse_args()

    def csvs(root):
        out = {}
        if args.recurse:
            for dp, _, fns in os.walk(root):
                for fn in fns:
                    if fn.lower().endswith(".csv"):
                        p = os.path.join(dp, fn)
                        out[os.path.relpath(p, root)] = p
        else:
            for fn in sorted(os.listdir(root)):
                if fn.lower().endswith(".csv"):
                    out[fn] = os.path.join(root, fn)
        return out

    A, B = csvs(args.a), csvs(args.b)
    rows = []
    for rel in sorted(A):
        rows.append({**compare_file(A[rel], B.get(rel, os.path.join(args.b, rel))),
                     "rel": rel})
    extra = sorted(set(B) - set(A))

    res = {"label": args.label, "a": os.path.abspath(args.a),
           "b": os.path.abspath(args.b),
           "n_files_a": len(A), "n_files_b": len(B),
           "only_in_b": extra, "files": rows}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(res, f, indent=1)

    def fmt(v):
        if v is None:
            return "-"
        if v == 0:
            return "0"
        return f"{v:.3e}"

    print(f"\n### {args.label}")
    print(f"A = {args.a}   B = {args.b}\n")
    print("| file | rows | max abs diff | max rel diff | worst column | verdict |")
    print("|---|---:|---:|---:|---|---|")
    for r in rows:
        print(f"| `{r['rel']}` | {r.get('rows', r.get('rows_a','-'))} | "
              f"{fmt(r.get('max_abs_diff'))} | {fmt(r.get('max_rel_diff'))} | "
              f"{r.get('worst_column') or '-'} | **{r['verdict']}** |")
    if extra:
        print(f"\nOnly in B: {', '.join(extra)}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
