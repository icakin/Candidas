# Phylogenomics + comparative genomics — methods and results

## Genomes (10 annotated proteomes, NCBI Datasets v2)

| Tag | Accession | Strain | Proteins | Phenotype measured |
|---|---|---|---|---|
| auris_cladeI | GCF_002759435.1 | B8441 | 5424 | Clade I |
| auris_cladeII | GCF_003013715.1 | B11220 | 5327 | Clade II |
| auris_cladeIII | GCF_002775015.1 | B11221 | 5521 | Clade III |
| auris_cladeIV | GCA_008275145.1 | B11245 | 5506 | Clade IV |
| haemulonii | GCF_002926055.2 | B11899 | 5249 | Hae |
| duobushaemulonii | GCF_002926085.2 | B09383 | 5173 | Duo |
| pseudohaemulonii | GCF_003013735.1 | B12108 | 5134 | — |
| parapsilosis | GCF_000182765.1 | (CDC317 ref) | 5830 | para |
| albicans | GCF_000182965.3 | SC5314 | 6030 | — |
| lusitaniae | GCF_014636115.1 | FDAARGOS_655 | 5485 | — |

Note: clade IV is a GenBank (GCA) annotation and is more sparsely named than the
RefSeq genomes (15% vs ~61–73% named products). **All analyses below are
sequence-homology based, not annotation-text based**, so this does not bias results.

## Phylogenomic tree
- Reciprocal best hits (DIAMOND blastp, e<1e-10) against *C. auris* clade I,
  requiring 1:1 orthology in every genome → **2,931 single-copy orthologs**.
- Random sample of 600 (seed 42), each aligned with MAFFT `--auto`, columns with
  >20% gaps trimmed, concatenated → **supermatrix 289,787 aa × 10 taxa**,
  occupancy 95.3–100%.
- ML tree, FastTree LG+Γ. **Every node 1.00 SH-like support; 0/7 bad splits.**
- Rooted on *C. albicans* / *C. parapsilosis* / *C. lusitaniae*.

### Topology (all nodes 1.00)
- **The four *C. auris* clades + the *haemulonii* complex form one clade** (Metschnikowiaceae).
- *C. auris* clades resolved: (I, III) < II < IV.
- *haemulonii* complex: (*haemulonii*, (*duobushaemulonii*, *pseudohaemulonii*)).
- *C. parapsilosis* + *C. albicans* (Debaryomycetaceae) fall outside that clade.

### Divergence (patristic, subs/site)
| comparison | distance |
|---|---|
| max within *C. auris* (4 clades) | 0.013 |
| *C. auris* ↔ *C. haemulonii* | 0.330 (**25×**) |
| *C. auris* ↔ *C. parapsilosis* | 1.146 (**87×**) |

## Thermal-machinery census (homology, DIAMOND e<1e-20, ≥40% id, ≥60% cov)
Query proteins extracted from the well-annotated *C. auris* B8441 and *C. albicans*
SC5314 proteomes (52 queries across 9 families), searched against all 10 genomes.

**Copy number is invariant across all ten genomes** for HSP90 (4), sHSP (1),
trehalose pathway (4), calcineurin (2), fatty-acid desaturase (1) and ergosterol
(2); HSP70 is 9–10 (within detection tolerance).

**The single exception is alternative oxidase (AOX): 2 copies in the two
Debaryomycetaceae (*C. parapsilosis*, *C. albicans*) vs 1 in all C. auris and haemulonii-complex
and *Clavispora*.** This is consistent with *C. parapsilosis*'s ~3× higher
measured per-cell respiration and its documented reliance on a parallel
non-phosphorylating respiratory chain — but the difference tracks **family**, so
it is confounded with phylogeny, and copy number is not flux.

## Proteome thermostability signatures — and why they do NOT support a correlation
Computed per proteome: IVYWREL fraction, charged−polar (CvP), (E+K)/(Q+H).

Naive correlation against measured phenotype (n = 7 genomes with data) appears
strong: IVYWREL vs CUE optimum rho = +0.89, p = 0.007; CvP vs CUE optimum
rho = +0.79, p = 0.036.

**This is pseudoreplication and must not be reported as a result.** Three tests:
1. **Collapse the four near-identical *auris* clades to one phylogenetic point**
   (effective n = 4): IVYWREL rho = +0.80, **p = 0.20**; CvP rho = +0.60, p = 0.40.
   With n = 4 even a perfect rank match gives p = 0.083 — significance is
   mathematically unreachable.
2. **Within *C. auris***: IVYWREL varies by only 0.00048 across the four clades.
   The genome signature is flat. We do NOT pair this with a within-*auris*
   phenotype spread: the fitted optima give 2.74 °C and the model-free optima
   4.06 °C, but only one of six pairwise contrasts survives correction for six
   comparisons, and that one has a 0.81 °C point estimate inside an interval
   running to -6.37 °C. The within-*auris* fine structure is not resolved by
   these data (see CLADE_SPREAD_VERDICT.md). The flat genome signature is
   reported on its own.
3. **CvP separates the two families completely, with no overlap**
   (Debaryomycetaceae −0.028 to −0.015; Metschnikowiaceae +0.004 to +0.022),
   independent of thermal optimum. The "correlation" is a family marker.

## Conclusion
Across ten genomes spanning ~4.5 °C in growth optimum and ~1.6-fold in fever cost,
**neither gene content nor coarse proteome composition predicts thermal
phenotype**. The between-species contrasts carry this: *C. auris* is markedly more
thermotolerant than the *haemulonii* complex despite being 25× closer to it than
to *C. parapsilosis*, while thermal-machinery copy number is identical across all
ten genomes and proteome composition is effectively identical. Conversely
*C. parapsilosis*, the most distant taxon measured, converges with the
*haemulonii* complex on a cool optimum and an expensive fever. We do not quote a
within-*C. auris* clade spread: it does not survive pairwise contrasts on the
model-free optima (CLADE_SPREAD_VERDICT.md), and the premise does not need it.

Thermal performance in this group is therefore **evolutionarily labile and not
genomically predictable**: it must arise from regulatory and quantitative
differences (expression, allelic variation, protein-level tuning) rather than
from the presence or absence of thermal-response genes. This is a positive
argument for direct phenotyping — the physiology cannot be read off the genome.

## Naming convention
Recent revisions place the *auris*/*haemulonii* group in *Candidozyma* (2024) and
*C. parapsilosis* in *Lodderomyces* (2026). **This project uses *Candida* throughout**,
matching config.R and the clinical literature, so every tip abbreviates unambiguously
to *C.* One sentence in the Methods should acknowledge the reclassifications and state
that the older combinations are retained for continuity with the clinical literature —
a referee will otherwise raise it.

## Files
- `pg_rooted.nwk` — final rooted phylogenomic tree
- `supermatrix.faa` — 600-gene concatenated alignment
- `core_orthologs.json` — the 2,931 single-copy ortholog map
- `thermal_census.csv`, `proteome_signatures.csv`, `signature_vs_phenotype.csv`
- `FIG_phylogenomics.png`
