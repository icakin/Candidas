#!/usr/bin/env python3
"""03_apply_medium_fba.py - apply a resolved per-model YMS medium and run the three
preregistered nutrient scenarios (min / proxy / rich) with a POOLED amino-acid
carbon budget (not 20 free caps). Verified working on iDC1003 (parapsilosis).

    python3 gem/03_apply_medium_fba.py \\
        gem/models/parapsilosis_iDC1003.xml \\
        gem/medium_iDC1003_parapsilosis.csv

Needs: cobra, pandas. Writes gem/medium_fba_<model>.csv.

Per-model medium maps resolve the abstract medium_YMS.csv onto each model's real
exchange ids + namespace (iDC1003 is KEGG; iRV973 may differ - make its own map
after 02_inspect_models.py). Columns: role, exchange_id, kegg, name, setting, note.
`setting`: OPEN | FIT (carbon, scaled) | AA_POOL (pooled budget) | ABSENT.
"""
import sys, re
from pathlib import Path
import pandas as pd, cobra
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

GEM = Path(__file__).resolve().parent
GLC_BOUND = 10.0                                   # FIT placeholder (calibration knob)
AA_C_BUDGET = {"scenario_min": 0.0, "scenario_proxy": 5.0, "scenario_rich": 20.0}

def carbon_atoms(note):
    m = re.search(r"nC=(\d+)", str(note)); return int(m.group(1)) if m else 6

def run(model, med, scenario):
    EX = {r.id for r in model.reactions
          if r.id.startswith("EX_") or r.id.startswith("Drain_")}
    aa = med[med.setting == "AA_POOL"]
    with model:
        for r in model.reactions:
            if r.id in EX: r.lower_bound = 0.0
        for _, row in med.iterrows():
            rid = str(row.exchange_id)
            if rid not in EX: continue
            if row.setting == "OPEN": model.reactions.get_by_id(rid).lower_bound = -1000
            elif row.setting == "FIT": model.reactions.get_by_id(rid).lower_bound = -GLC_BOUND
        budget = AA_C_BUDGET[scenario]
        aa_ex = [r for r in aa.exchange_id if r in EX]
        if aa_ex and budget > 0:
            for rid in aa_ex: model.reactions.get_by_id(rid).lower_bound = -1000
            expr = sum(carbon_atoms(med.loc[med.exchange_id == rid, "note"].iloc[0])
                       * model.reactions.get_by_id(rid).reverse_variable for rid in aa_ex)
            model.add_cons_vars(model.problem.Constraint(expr, ub=budget, name="aa_pool_C"))
            model.solver.update()
        try: g = model.slim_optimize()
        except Exception: g = float("nan")
    return g

def main():
    if len(sys.argv) < 3:
        print("usage: python3 03_apply_medium_fba.py <model.xml> <medium_map.csv>"); return
    model = cobra.io.read_sbml_model(sys.argv[1])
    med = pd.read_csv(sys.argv[2])
    rows = [dict(scenario=s, growth_per_h=round(run(model, med, s), 4))
            for s in ("scenario_min", "scenario_proxy", "scenario_rich")]
    df = pd.DataFrame(rows)
    out = TABLES / f"medium_fba_{Path(sys.argv[1]).stem}.csv"
    df.to_csv(out, index=False)
    print(df.to_string(index=False)); print("\nwrote", out)
    print("Sanity: growth must rise min->proxy->rich. If min is already high AND "
          "flat vs rich, an exchange or enzyme-free bypass is leaking.")

if __name__ == "__main__":
    main()
