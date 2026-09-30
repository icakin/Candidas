# Curated-model inspection report

## parapsilosis_iDC1003.xml
- id: `model_Cparapsilosis_Clean`
- reactions: 2162   metabolites: 1637   genes: 1003
- namespace: **KEGG/ModelSEED-like**   (e.g. `C00350__extr`, `C15777__cyto`)
- compartments: {'C_00001': 'extracellular', 'C_00002': 'cytoplasmic', 'C_00003': 'mitochondria', 'C_00004': 'endoplasmic_reticulum', 'C_00005': 'peroxisome'}
- objective / biomass: obj=['e_Biomass__cyto']  biomass-like=['EX_Biomass__cyto', 'e_Biomass__cyto']
- exchanges: 358  (e.g. EX_C00214__extr, EX_CO2__extr, EX_C01217__extr, EX_Oxygen__extr, EX_C00881__extr, EX_C02341__extr)
- gene id sample: ['CPAR2_806240', 'CPAR2_805390', 'CPAR2_804060', 'CPAR2_212320', 'CPAR2_210140']
- **default-medium growth: 0.1725 /h**

## Harmonisation checklist (compare the two curated models)
- [ ] reconcile namespace
- [ ] reconcile compartments
- [ ] reconcile biomass id
- [ ] reconcile exchange prefix
- [ ] common biomass composition / GAM / NGAM / P-O convention
- [ ] map medium_YMS.csv exchange hints onto each model's real EX ids
