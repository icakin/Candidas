# etcGEM reconstruction layer — run order

Runs on the Mac (cobra + diamond; a good MILP solver — CPLEX/Gurobi academic — for
gap-filling and later GECKO). Produces four comparable base GEMs on YMS; the
enzyme (GECKO) and temperature (etcGEM) layers come after these.

Setup:  `pip install cobra pandas`  and DIAMOND in PATH.

1. **`bash 10_fetch_curated_models.sh`** — get iRV973 (auris) and iDC1003
   (parapsilosis) SBML into `gem/models/`. Do NOT substitute CarveMe/gapseq drafts.

2. **`python3 11_inspect_models.py`** — reports each model's namespace, compartments,
   biomass id, exchanges, genes, default growth → `gem/model_report.md`. Read this
   before step 3: it tells you what has to be reconciled and what to put in the
   `ID_MAP` of step 12.

3. **Harmonise** (manual, guided by the report): one namespace, compartment codes,
   biomass/GAM/NGAM/P-O convention across auris + parapsilosis. This is a curation
   step, not a script — the report checklist drives it. Keep gene-supported vs
   gap-filled reactions tracked.

4. **`python3 12_apply_medium_fba.py gem/models/<model>.xml`** — maps `medium_YMS.csv`
   onto the model, runs the 3 scenarios (min/proxy/rich). Edit `ID_MAP` first so the
   BiGG hints resolve to real EX ids. Growth must rise min→proxy→rich; if `min`
   already grows fast, an exchange or enzyme-free bypass is leaking — fix before the
   enzyme layer. Replace the even-split AA bound with a proper carbon-weighted
   pooled constraint for the real run.

5. **`python3 13_transfer_orthologs.py`** (twice) — build haemulonii and
   duobushaemulonii drafts by RBH orthology transfer from the harmonised template.
   Then gap-fill each against YMS. Dropped reactions are candidates, NOT proven
   absences — no unique-capability claims without gene/orthology/localization evidence.

Output of this layer: four base GEMs (auris, parapsilosis curated; haemulonii,
duobushaemulonii draft) on a common namespace, all runnable on YMS.

## Then (separate, later)
- GECKO ecModel layer + kcat (DLKcat) — verify the enzyme pool BINDS (utilisation,
  shadow price, 10-20% relax) before trusting any separation.
- etcGEM temperature layer — calibrate ONLY low-dim shared thermal hyperparameters
  on auris+parapsilosis (frozen), predict haemulonii & duobushaemulonii TPCs.
  Calibration target: `results/tables/gem_growth_tpc_target.csv`.
- Benchmarks: taxon-blind mean-TPC + sequence-ablated etcGEM.
- Primary endpoint: normalized TPC shape (never held-out mu_max). See GEM_pipeline_plan.md.
