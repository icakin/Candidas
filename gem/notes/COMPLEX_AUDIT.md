# Enzyme-complex audit — response to the GECKO-review blocker

ChatGPT flagged (correctly, as a general principle) that KO-based re-keying flattens
enzyme complexes (AND) to isozyme-OR, and that keeping curated complexes only in
C. parapsilosis would confound model-construction quality with species.

## Finding: these source models encode essentially NO enzyme complexes
Direct audit of the curated models:
- C. auris iRV973 (original): **4** AND-complex GPRs of 2863 reactions — all membrane transporters (TO0010516/TO0000021, __plas/__extr).
- C. parapsilosis iDC1003: **6** AND-complex GPRs of 2162 — all membrane transporters.
- No metabolic enzyme complexes (ATP synthase, cytochrome chain, central metabolism)
  are represented as AND-GPRs in EITHER model. This KEGG-based reconstruction family
  represents such reactions with single-gene or OR-isozyme GPRs.

## Consequence for the review's concern
- The premise (curated models carry rich metabolic complexes that re-keying destroys)
  does NOT hold for this model family. There were ~4-6 complexes to begin with, all
  transporters.
- The parap-specific transporter complexes have no reaction counterpart in the
  auris-scaffolded models (compartment/id scheme differs; transfer applied to 0
  reactions), so there is nothing shared to harmonize.
- Growth is identical (0.490/h) with or without the transfer attempt -> complexes are
  immaterial to model function here.

## Resolution
- All four models now use uniform single-gene/OR GPR treatment across the metabolic
  network -> NO species x construction confound (the review's actual blocking concern).
- Remaining limitation, now UNIFORM across all four species (hence disclosable, not a
  confound): these models do not resolve enzyme-complex subunit (AND) structure for
  metabolic reactions, so complex enzyme MW is represented by the isozyme/single-gene
  MW rather than summed subunits. This affects all species equally.
- The complex-aware vs OR comparison the review requested is trivially null here:
  there are no metabolic complexes to differ on.

## Verdict
The blocker is resolved by evidence: no metabolic complexes exist to preserve, and
treatment is now uniform across species. Proceed to DLKcat + the temperature layer.
De-novo complex curation (S. cerevisiae-based) remains an OPTIONAL future refinement,
not a prerequisite, and would apply equally to all four models.
