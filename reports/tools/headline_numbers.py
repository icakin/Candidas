"""
headline_numbers.py -- pull every load-bearing number of the manuscript out of
one (results tree, etc-GEM supp_data tree) pair, so the committed set and the
regenerated set can be put side by side.

    cauris_etcgem/.venv/bin/python reports/tools/headline_numbers.py \
        --tables results/tables \
        --supp   cauris_etcgem/strains/eci_cauris/outputs/supp_data \
        --out    reports/tools/headline_committed.json

Read-only.

Most numbers are read straight from the CSVs the pipeline writes. Three are
NOT written to any file by the pipeline and are recomputed here exactly as the
scripts compute them:

  * ICC1, F(3,8), p        -- 13_supplementary_figures.R panel a, from
                              aov(capacity ~ clade) on capacity_isolates.csv
  * r(capacity, fever tax at 40 C) and r(capacity, respiration at 32-36 C)
                           -- 12_main_figures.R panels d and e, quoted in the
                              manuscript as r = -0.86 and r = -0.64.
                              Only the POINT estimate is reproduced; the CI that
                              12 prints beside it comes from an UNSEEDED
                              sample.int() bootstrap and is not reproducible.
"""
import argparse, json, os
import numpy as np
import pandas as pd
from scipy import stats

CLADES = ["I", "II", "III", "IV"]
GROUP2CLADE = {"Clade1": "I", "Clade2": "II", "Clade3": "III", "Clade4": "IV"}


def rd(d, f):
    p = os.path.join(d, f)
    return pd.read_csv(p) if os.path.exists(p) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--supp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    T, S = a.tables, a.supp
    R = {"tables": os.path.abspath(T), "supp": os.path.abspath(S)}

    # ---- E_G, E_R, dE, Topt(CUE), tax37, tax40 (all in fig_values.csv) -------
    fv = rd(T, "fig_values.csv")
    if fv is not None:
        fv = fv.set_index("Group")
        R["per_taxon"] = {
            g: dict(
                E_G=float(fv.loc[g, "E"]), E_G_lo=float(fv.loc[g, "E_lo"]), E_G_hi=float(fv.loc[g, "E_hi"]),
                E_R=float(fv.loc[g, "ER"]), E_R_lo=float(fv.loc[g, "ER_lo"]), E_R_hi=float(fv.loc[g, "ER_hi"]),
                dE=float(fv.loc[g, "dE"]), dE_lo=float(fv.loc[g, "dE_lo"]), dE_hi=float(fv.loc[g, "dE_hi"]),
                dE_credible=str(fv.loc[g, "dE_credible"]),
                Topt_CUE=float(fv.loc[g, "Tcue"]), Topt_CUE_lo=float(fv.loc[g, "Tcue_lo"]),
                Topt_CUE_hi=float(fv.loc[g, "Tcue_hi"]),
                tax37=float(fv.loc[g, "tax37"]), tax37_lo=float(fv.loc[g, "tax37_lo"]), tax37_hi=float(fv.loc[g, "tax37_hi"]),
                tax40=float(fv.loc[g, "tax40"]), tax40_lo=float(fv.loc[g, "tax40_lo"]), tax40_hi=float(fv.loc[g, "tax40_hi"]),
                fever_cost=float(fv.loc[g, "fever_cost"]),
            ) for g in fv.index
        }
    else:
        R["per_taxon"] = None

    # ---- the 10 pairwise fever contrasts ------------------------------------
    fc = rd(T, "fig_contrasts.csv")
    if fc is not None:
        R["contrasts"] = [
            dict(a=r.a, b=r.b, ratio=float(r.ratio), lo=float(r.lo), hi=float(r.hi),
                 credible=str(r.credible)) for r in fc.itertuples()
        ]
        R["contrasts_n"] = int(len(fc))
        R["contrasts_n_resolved"] = int(sum(str(x).upper() == "TRUE" for x in fc.credible))
    else:
        R["contrasts"] = None

    # ---- Spearman rho (Topt vs fever cost across the 5 taxa) ----------------
    # 12 computes a POSTERIOR of rho draw-by-draw and only prints it; the point
    # estimate on the posterior medians is reproducible from fig_values.csv.
    if fv is not None:
        rho, p = stats.spearmanr(fv["Topt"].values, fv["fever_cost"].values)
        R["spearman_Topt_vs_fever_cost"] = dict(
            rho_on_posterior_medians=float(rho), p=float(p), n=int(len(fv)),
            note="12_main_figures.R reports the POSTERIOR MEDIAN of rho computed "
                 "draw-by-draw (printed, never written to a CSV). This is the "
                 "point estimate on the posterior medians -- the closest "
                 "reproducible proxy from the committed tables.")

    # ---- etc-GEM: R2 and kcat_scale per clade -------------------------------
    r2 = rd(S, "calibration_r2.csv")
    if r2 is not None:
        r2 = r2.set_index("clade")
        R["etcgem"] = {cl: dict(
            r2_calibrated=float(r2.loc[cl, "r2_calibrated"]),
            r2_apriori=float(r2.loc[cl, "r2_apriori"]),
            kcat_scale=float(r2.loc[cl, "kcat_scale"]),
            dTopt=float(r2.loc[cl, "dTopt"]),
            dCp_scale=float(r2.loc[cl, "dCp_scale"]),
            apriori_Topt=float(r2.loc[cl, "apriori_Topt"]),
        ) for cl in CLADES if cl in r2.index}

    # ---- the 12 isolate capacities + ICC1 / F(3,8) / p -----------------------
    cap = rd(S, "capacity_isolates.csv")
    if cap is not None:
        R["capacity_isolates"] = [
            dict(clade=r.clade, isolate=r.isolate, capacity=float(r.capacity),
                 peak=float(r.peak)) for r in cap.itertuples()]
        groups = [cap.capacity[cap.clade == c].values for c in CLADES]
        k, n = len(groups), len(groups[0])
        grand = np.concatenate(groups).mean()
        ssb = n * sum((g.mean() - grand) ** 2 for g in groups)
        ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
        dfb, dfw = k - 1, k * (n - 1)
        msb, msw = ssb / dfb, ssw / dfw
        F = msb / msw
        R["capacity_anova"] = dict(
            ICC1=float((msb - msw) / (msb + (n - 1) * msw)),
            F=float(F), df=[dfb, dfw], p=float(stats.f.sf(F, dfb, dfw)))

    # ---- r = -0.86 and r = -0.64 (12_main_figures.R panels d and e) ----------
    tax = rd(T, "carbon_tax_isolate.csv")
    der = rd(T, "derived_N0_R_results_with_carbon.csv")
    if cap is not None and tax is not None and der is not None:
        # fever tax at 40 C per isolate
        fev = tax[tax["T_C"] == 40][["Isolate", "tax"]].rename(
            columns={"Isolate": "otu_name", "tax": "fever_tax40"})
        # measured respiration in the growth-optimal band, exactly as 12 does
        db = der[(der["keep"] == True) & (der["T"] >= 32) & (der["T"] <= 36)]
        resp = (db.groupby("otu_name")["respiration_fgC_h"].mean()
                  .reset_index().rename(columns={"respiration_fgC_h": "resp_opt"}))
        iso = cap.merge(fev, left_on="isolate", right_on="otu_name")
        iso = iso.merge(resp, left_on="isolate", right_on="otu_name")
        out = {"n": int(len(iso))}
        for nm, col in (("r_capacity_vs_fever_tax40", "fever_tax40"),
                        ("r_capacity_vs_respiration", "resp_opt"),
                        ("r_capacity_vs_peak_growth", "peak")):
            if col in iso.columns and len(iso) > 2:
                r = float(np.corrcoef(iso["capacity"].values, iso[col].values)[0, 1])
                pv = float(stats.pearsonr(iso["capacity"].values, iso[col].values)[1])
                out[nm] = dict(r=r, p=pv)
        out["note"] = ("Point estimates only. 12_main_figures.R prints a 95% CI "
                       "from replicate(5000, sample.int(...)) with NO seed set, "
                       "so that CI is not reproducible run to run.")
        R["isolate_correlations"] = out

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(R, f, indent=1)
    print(json.dumps(R, indent=1)[:4000])
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
