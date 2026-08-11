"""
anchor_and_ceiling_test.py -- C1 known-inconsistency tests that need a BUILT MODEL.

Read-only. Writes one JSON to reports/tools/anchor_and_ceiling_test.json.
Run with the pinned venv, from anywhere:

    cauris_etcgem/.venv/bin/python reports/tools/anchor_and_ceiling_test.py

Two tests:

 (A) ANCHOR. _supp_common.build() calibrates the enzyme budget to ONE measured
     growth rate: ANCHOR_MU = 0.95 at ANCHOR_T = 307.15 K. The comment above it
     still names the previous anchor ("Was an arbitrary 1.1 at 40 C"). Evaluate
     the COMMITTED per-clade parameters under BOTH anchors and report the R2
     each gives against the measured clade means, next to the R2 recorded in
     outputs/supp_data/calibration_r2.csv. Whichever anchor reproduces the
     committed R2 is the anchor the committed file was actually built with.

 (B) CEILING. outputs/supp_data/calibration_curves.csv gives Clade IV a
     `calibrated` curve peaking at ~0.9314, above its own `apriori` plateau of
     ~0.7955 -- with kcat_scale = 0.710, i.e. LESS enzyme turnover than a priori.
     Rebuild both curves from one build and see whether that is reproducible.
"""
import json, os, sys, importlib
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = os.path.dirname(os.path.dirname(HERE))            # .../Candidas
STRAIN = os.path.join(ROOT, "cauris_etcgem", "strains", "eci_cauris")
sys.path.insert(0, os.path.join(STRAIN, "scripts"))

import pandas as pd
import _supp_common as SC

SUPP   = os.path.join(STRAIN, "outputs", "supp_data")
CLADES = ["I", "II", "III", "IV"]
OUT    = os.path.join(HERE, "anchor_and_ceiling_test.json")

r2tab = pd.read_csv(os.path.join(SUPP, "calibration_r2.csv")).set_index("clade")
curtab = pd.read_csv(os.path.join(SUPP, "calibration_curves.csv"))

iso = SC.isolate_frame()
iso["clade"] = iso.group.map(SC.GROUP2CLADE)
iso = iso.dropna(subset=["clade"])


def clade_mean(cl):
    s = iso[iso.clade == cl].groupby("temp_C")["rate_per_h"].mean()
    return s.index.values.astype(float), s.values.astype(float)


def build_with_anchor(mu, T):
    """Rebuild the model with a patched (ANCHOR_MU, ANCHOR_T).

    build() hard-codes the anchor inside its body, so the only honest way to
    swap it without editing the shipped file is to exec a patched copy of the
    source. Nothing on disk is modified.
    """
    src = open(os.path.join(STRAIN, "scripts", "_supp_common.py")).read()
    old = "ANCHOR_MU, ANCHOR_T = 0.95, 307.15"
    assert old in src, "anchor line not found -- _supp_common.py has changed"
    src = src.replace(old, f"ANCHOR_MU, ANCHOR_T = {mu!r}, {T!r}")
    ns = {"__name__": "_patched_supp_common",
          "__file__": os.path.join(STRAIN, "scripts", "_supp_common.py")}
    exec(compile(src, "_supp_common(patched)", "exec"), ns)
    return ns


results = {"anchors": {}, "ceiling": {}}

for label, (mu, T) in {
    "current_0.95@307.15K": (0.95, 307.15),
    "older_1.1@313.15K":    (1.1,  313.15),
}.items():
    print(f"\n=== building with anchor {label} ===", flush=True)
    ns = build_with_anchor(mu, T)
    pm, _ = ns["build"]()
    curve = ns["curve"]
    tf = np.arange(18, 46.5, 0.5)
    entry = {}
    for cl in CLADES:
        mt, mv = clade_mean(cl)
        p = [float(r2tab.loc[cl, "dTopt"]),
             float(r2tab.loc[cl, "dCp_scale"]),
             float(r2tab.loc[cl, "kcat_scale"])]
        cam = curve(pm, p, mt)
        r2 = 1 - np.sum((cam - mv) ** 2) / np.sum((mv - mv.mean()) ** 2)
        rmse = float(np.sqrt(np.mean((cam - mv) ** 2)))
        ap = curve(pm, [0, 1, 1], tf)
        ca = curve(pm, p, tf)
        entry[cl] = dict(
            params=p,
            r2_recomputed=float(r2),
            r2_committed=float(r2tab.loc[cl, "r2_calibrated"]),
            r2_abs_diff=float(abs(r2 - r2tab.loc[cl, "r2_calibrated"])),
            rmse_recomputed=rmse,
            rmse_committed=float(r2tab.loc[cl, "rmse_calibrated"]),
            apriori_peak=float(np.max(ap)),
            apriori_Topt=float(tf[int(np.argmax(ap))]),
            calibrated_peak=float(np.max(ca)),
            calibrated_Topt=float(tf[int(np.argmax(ca))]),
        )
        print(f"  {cl}: R2 recomputed={r2: .4f}  committed={r2tab.loc[cl,'r2_calibrated']: .4f}"
              f"  |diff|={abs(r2-r2tab.loc[cl,'r2_calibrated']):.4f}"
              f"  apriori_peak={np.max(ap):.4f}  calib_peak={np.max(ca):.4f}", flush=True)
    results["anchors"][label] = entry

# ---- (B) ceiling: committed columns, one build --------------------------------
for cl in CLADES:
    d = curtab[curtab.clade == cl]
    results["ceiling"][cl] = dict(
        committed_apriori_peak=float(d.apriori.max()),
        committed_apriori_Topt=float(d.temp_C.values[int(np.argmax(d.apriori.values))]),
        committed_calibrated_peak=float(d.calibrated.max()),
        committed_calibrated_Topt=float(d.temp_C.values[int(np.argmax(d.calibrated.values))]),
        committed_kcat_scale=float(r2tab.loc[cl, "kcat_scale"]),
        calibrated_exceeds_apriori=bool(d.calibrated.max() > d.apriori.max()),
        excess=float(d.calibrated.max() - d.apriori.max()),
    )

with open(OUT, "w") as f:
    json.dump(results, f, indent=1)
print("\nwrote", OUT)
