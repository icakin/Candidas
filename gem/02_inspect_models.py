#!/usr/bin/env python3
"""02_inspect_models.py - report the structure of each curated GEM so harmonisation
decisions are made from fact, not assumption. Run BEFORE 03 (medium mapping).

    python3 gem/02_inspect_models.py

Needs: cobra (pip install cobra). Reads gem/models/*.xml, writes gem/model_report.md.
Harmonisation depends entirely on what this reveals: namespace prefix, compartment
codes, biomass reaction id, objective, exchange ids, gene ids, and whether each
model grows under its default medium.
"""
from pathlib import Path
import re, collections
import cobra
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

GEM = Path(__file__).resolve().parent
MODELS = GEM / "models"
OUT = NOTES / "model_report.md"

def guess_namespace(ids):
    # BiGG metabolites look like glc__D_c; MetaCyc/KEGG/other differ
    if any(re.search(r"__[A-Za-z]_[a-z]$", i) for i in ids): return "BiGG-like"
    if any(i.startswith("C0") or i.startswith("cpd") for i in ids): return "KEGG/ModelSEED-like"
    return "other/custom"

def report(path):
    m = cobra.io.read_sbml_model(str(path))
    mets = [x.id for x in m.metabolites]
    comps = dict(m.compartments)
    ex = [r.id for r in m.reactions if r.id.startswith("EX_") or r.boundary]
    bio = [r.id for r in m.reactions if "biomass" in r.id.lower()
           or "growth" in r.id.lower()]
    obj = [r.id for r in m.reactions if r.objective_coefficient != 0]
    try:
        sol = m.slim_optimize(); grow = f"{sol:.4f}"
    except Exception as e:
        grow = f"ERROR: {e}"
    L = []
    L.append(f"## {path.name}")
    L.append(f"- id: `{m.id}`")
    L.append(f"- reactions: {len(m.reactions)}   metabolites: {len(m.metabolites)}   genes: {len(m.genes)}")
    L.append(f"- namespace: **{guess_namespace(mets)}**   (e.g. `{mets[0]}`, `{mets[len(mets)//2]}`)")
    L.append(f"- compartments: {comps}")
    L.append(f"- objective / biomass: obj={obj}  biomass-like={bio[:5]}")
    L.append(f"- exchanges: {len(ex)}  (e.g. {', '.join(ex[:6])})")
    L.append(f"- gene id sample: {[g.id for g in list(m.genes)[:5]]}")
    L.append(f"- **default-medium growth: {grow} /h**")
    return "\n".join(L), m

def main():
    xmls = sorted(MODELS.glob("*.xml"))
    if not xmls:
        print("No .xml in", MODELS, "- run 01_fetch_curated_models.sh first."); return
    blocks=[]; models={}
    for p in xmls:
        try:
            b, m = report(p); blocks.append(b); models[p.stem]=m
        except Exception as e:
            blocks.append(f"## {p.name}\n- FAILED to load: {e}")
    # cross-model comparability check
    cmp=["\n## Harmonisation checklist (compare the two curated models)"]
    for k in ("namespace","compartments","biomass id","exchange prefix"):
        cmp.append(f"- [ ] reconcile {k}")
    cmp.append("- [ ] common biomass composition / GAM / NGAM / P-O convention")
    cmp.append("- [ ] map medium_YMS.csv exchange hints onto each model's real EX ids")
    OUT.write_text("# Curated-model inspection report\n\n" + "\n\n".join(blocks) +
                   "\n" + "\n".join(cmp) + "\n")
    print("wrote", OUT)
    for b in blocks: print("\n"+b)

if __name__ == "__main__":
    main()
