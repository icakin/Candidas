# High-oxygen refit of the per-well likelihood model (Supplementary Note 4)

Every trace was truncated where dissolved oxygen first fell below 50% of its value at the
window start (90 min), and the likelihood model (lik_model.py, 199 bootstrap draws) was rerun
on the truncated series; wells with fewer than 60 points or less than 2 h of truncated data
are reported as unidentifiable.

    python3 reprocess/hiO2/run_hiO2.py 0.5        # ~25 min on 2 cores; writes lik_hiO2_0.50.csv
    python3 reprocess/hiO2/compare_hiO2.py 0.5    # comparison with reprocess/lik_results_labelled.csv

Run 29 Sep 2026 (numpy 2.4.4, scipy 1.17.1). Outputs: lik_hiO2_0.50.csv (per well),
hiO2_comparison_0.50.csv (per group: model-free CUE optimum and 38->40 C change in K/(rN), both arms).
