# Curated-model harmonisation — auris (iRV973) + parapsilosis (iDC1003)

## Result: the two curated models are already in a common framework
Both were built with the same pipeline (KEGG namespace, `__cyto/__extr/__mito/__pero`
suffixes, `EX_C#####__extr` exchanges, `e_Biomass__cyto` objective). The
namespace-reconciliation that was flagged as the hard part is largely unnecessary.

| | C. auris iRV973 | C. parapsilosis iDC1003 |
|---|---|---|
| model id | model_Cauris2020 | model_Cparapsilosis_Clean |
| reactions | 2863 | 2162 |
| metabolites | 2150 | 1637 |
| genes | 973 (CJI97_*) | 1003 (CPAR2_*) |
| compartments (met suffixes) | cyto, extr, mito, pero | cyto, extr, mito, pero |
| biomass objective | e_Biomass__cyto | e_Biomass__cyto |
| default growth | 0.49 /h (glc -7) | 0.17 /h (glc -2) |

**Shared: 1381 metabolite ids, 284 exchange reactions.** Same 4 compartments.

## Only per-model wrinkle: a few core-metabolite exchange-id conventions differ
| component | auris | parapsilosis |
|---|---|---|
| O2 | EX_C00007__extr | EX_Oxygen__extr |
| CO2 | EX_C00011__extr | EX_CO2__extr |
| ammonia | EX_Ammonia__extr | EX_C00014__extr |
| water | EX_C00001__extr | Drain_to_H2O__extr |
Handled by per-model medium maps (medium_iRV973_auris.csv, medium_iDC1003_parapsilosis.csv).

## Shared essentials / auxotrophy
Both require, beyond C source: O2, ammonia, sulfate, phosphate, and **biotin**
(biotin auxotrophy — neither synthesises it; YE supplies it). scenario_min must
include biotin or growth = 0.

## Amino-acid exchanges (for the pooled AA proxy)
- auris present 17/20 (missing Ala, Asn, Cys)
- parapsilosis present 17/20 (missing Ala, Cys, Trp)
- **cross-species intersection = 16** (both lack Ala, Cys; auris lacks Asn, parap lacks Trp).
  Use the 16-AA intersection for a fair cross-species proxy, or per-model 17 with a note.

## YMS scenario growth (pooled AA carbon budget; glucose FIT = 10 mmol/gDW/h)
| scenario | auris | parapsilosis |
|---|---|---|
| min (no AA) | 0.70 | 0.86 |
| proxy (AA C<=5) | 0.80 | 0.95 |
| rich (AA C<=20) | 1.00 | 1.21 |
Monotonic rise (no leak). Both reach the measured mu range (0.5-0.9/h) at plausible
glucose uptake — the calibration scale knob works, and absolute rates are NOT a
mismatch (the low default 0.17-0.49 was just glucose-limited default bounds).

## Status / next
- [x] both curated base models load, medium-mapped, scenarios validated
- [x] harmonisation confirmed trivial (shared framework)
- [ ] settle biomass/GAM/NGAM/P-O comparability between the two (check reaction ids)
- [ ] ortholog transfer -> haemulonii, duobushaemulonii drafts (13_transfer_orthologs.py; needs DIAMOND)
- [ ] GECKO ecModel layer; verify pool binds
- [ ] etcGEM temperature layer -> calibrate low-dim hyperparams on auris+parap, predict hae/duo
