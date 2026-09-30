#!/usr/bin/env python3
"""build_constraints.py - assemble gem_constraints.csv from the respirometry data.

Emits one row per SPECIES (auris clades pooled) with the measured quantities the
ecModel needs at 37 °C, plus the unit-conversion columns. Interspecies scope.

    python3 gem/audits/build_constraints.py

Writes results/tables/gem_constraints.csv.

TASK-0 WARNING (see GEM_pipeline_plan.md): the per-cell flux columns and the
fitted r do not reconcile under a naive reading. The C_FRACTION and the meaning
of growth_fgC_h MUST be confirmed before these numbers are used as hard bounds.
Everything below is emitted as measured values + explicit formulas, NOT as
finished mmol/gDW/h constraints.
"""
from pathlib import Path
import pandas as pd, numpy as np
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

T = RESULTS_TABLES

C_FRACTION = 0.48          # carbon mass fraction of dry weight - PLACEHOLDER, set from CHN
MW_C = 12.011              # g/mol

d = pd.read_csv(T / "derived_N0_R_results_with_carbon.csv")
d["Group"] = d.otu_name.str.split("_").str[0]
d["Species"] = d.Group.replace({"Clade1":"auris","Clade2":"auris",
                                "Clade3":"auris","Clade4":"auris"})
sizes = pd.read_csv(T / "otu_cell_sizes.csv").set_index("OTU")

# growth-positive curves near body temperature (36-38 C flank the un-assayed 37)
sub = d[(d["T"].isin([36,38])) & (d.keep==True) &
        (d.fit_valid==True) & (d.has_curvature==True)]

rows=[]
for sp, g in sub.groupby("Species"):
    otus = d[d.Species==sp].OTU.unique()
    cellC = sizes.loc[[o for o in otus if o in sizes.index],"cell_carbon_fg"].median()
    gdw_per_cell = cellC / C_FRACTION * 1e-15                       # g DW per cell
    growth_fgC = g.growth_fgC_h.median()      # fg C -> biomass  (UNITS: confirm task 0)
    resp_fgC   = g.respiration_fgC_h.median() # fg C -> CO2       (UNITS: confirm task 0)
    cue        = g.growth_C_per_C_h.median()
    r          = g.r.median()
    # per-cell fg C/h -> mmol C /gDW/h  =  (fgC/h /cell) / (gDW/cell) * 1e-15 g/fg / MW_C *1e3
    to_mmol = lambda fgc: (fgc * 1e-15) / MW_C * 1e3 / gdw_per_cell
    rows.append(dict(
        species=sp, n_isolates=d[d.Species==sp].OTU.nunique(),
        cell_carbon_fg=round(cellC,1), gDW_per_cell=f"{gdw_per_cell:.3e}",
        meas_r_per_h=round(r,4),
        meas_growth_fgC_cell_h=round(growth_fgC,1),
        meas_resp_fgC_cell_h=round(resp_fgC,1),
        meas_CUE=round(cue,4),
        pred_growth_mmolC_gDW_h=round(to_mmol(growth_fgC),3),
        pred_resp_mmolC_gDW_h=round(to_mmol(resp_fgC),3),
        carbon_source_uptake_bound="SET FROM YMS RECIPE",
    ))

out = pd.DataFrame(rows).set_index("species").reindex(
        ["auris","Hae","Duo","para"])
dest = T / "gem_constraints.csv"
out.to_csv(dest)
print("wrote", dest, "\n")
print(out.to_string())
print("\nTASK 0: confirm growth_fgC_h units and set C_FRACTION from CHN before "
      "treating pred_* columns as hard bounds. CUE is unit-free and usable now "
      "as the primary validation target.")
