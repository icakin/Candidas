# ITS phylogeny — methods and provenance

## Purpose
Place the taxa measured in this study on a molecular phylogeny, to test whether
thermal phenotype tracks relatedness.

## Sequences
Nine curated ITS reference sequences were downloaded from NCBI (E-utilities).
All are **type material**; RefSeq Targeted-Loci (NR_) records were preferred.

| Tip | Accession | Strain | In this study |
|---|---|---|---|
| *Candidozyma auris* | NR_154998.1 | CBS 10913 (T) | yes — clades I–IV |
| *Candidozyma haemulonii* | NR_130669.1 | CBS 5149 (T) | yes — 1724/1768/1769 |
| *C. haemulonii* var. *vulnera* | NR_202147.1 | CNM CL7239 (T) | context |
| *Candidozyma duobushaemulonii* | NR_130694.1 | CBS 7798 (T) | yes — 1770/1771 |
| *Candidozyma pseudohaemulonii* | NR_163771.1 | CBS 10004 (T) | context |
| *Lodderomyces parapsilosis* | NR_130673.1 | ATCC 22019 (T) | yes — 2051/2052/2053 |
| *Candida albicans* | NR_125332.1 | CBS 562 (T) | context |
| *Clavispora lusitaniae* | NR_130677.1 | CBS 6936 (T) | context |
| *Saccharomyces cerevisiae* | NR_132221.1 | S288C | outgroup (excluded, see below) |

## Alignment and tree
- Aligned with **MAFFT** (`--localpair --maxiterate 1000`, L-INS-i).
- Columns with >50% gaps trimmed (ragged SSU/LSU flanks) → **390 aligned columns**.
- **FastTree** (`-nt -gtr -gamma`), SH-like local support.
- Rooted on the Debaryomycetaceae clade (*C. albicans* + *L. parapsilosis*).

**S. cerevisiae was excluded from the final tree.** Its ITS identity to the ingroup
was ~32% (saturated), producing a branch length of 1.41 (vs 0.01–0.21 for all
other tips) and dragging *Clavispora* into a 0.00-support grouping — classic
long-branch attraction. Removing it stabilised the topology without changing any
ingroup relationship.

## Results (SH-like support)
- Metschnikowiaceae clade (*Candidozyma* spp. + *Clavispora*): **1.00**
- *C. haemulonii* complex (haemulonii, var. *vulnera*, duobushaemulonii, pseudohaemulonii): **0.91**
- *C. duobushaemulonii* + *C. pseudohaemulonii*: **0.98**
- *L. parapsilosis* + *C. albicans* (Debaryomycetaceae): **1.00**
- Position of *C. auris* within Metschnikowiaceae: **0.25** — unresolved by ITS.

Pairwise ITS identity (trimmed core): haemulonii/var. *vulnera* 95.7%;
duobushaemulonii/pseudohaemulonii 93.6%; within haemulonii complex 85–87%;
*auris* vs haemulonii complex 73–76%; *L. parapsilosis* vs *C. albicans* 74.7%;
*L. parapsilosis* vs *Candidozyma* 64.7–74.5%.

## Limitations (state in the paper)
1. **ITS cannot resolve the four *C. auris* clades** — they are near-identical at
   this locus. Clade structure is taken from published whole-genome phylogenies
   (Lockhart et al. 2017); clades are plotted separately in the trait panels only.
2. **The position of *C. auris*** relative to the haemulonii complex is not
   resolved by ITS (support 0.25); cite genome-based phylogenies for that node.
3. **No formal phylogenetic comparative test** (PGLS, Pagel's λ, Mantel) is
   reported: with 7 measured taxa — four of them clades of one species — such
   tests have effectively no power, and the distance matrix is bimodal
   (near-zero within *auris*, large between species). The tree is presented as
   descriptive context, not as an independent statistical test.

## Interpretation
Thermal phenotype cuts across the phylogeny. *L. parapsilosis* is the most
distantly related taxon measured (separate family) yet converges with the
*haemulonii* complex on a cool optimum (~31.5 °C) and an expensive fever
(~1.9×), whereas *C. auris* — the closest relative of the haemulonii complex —
has the warmest optima and the cheapest fever cost (1.24–1.68×). Within
*C. haemulonii*, isolates span 36–44 °C in growth ceiling. Thermal performance
in this group is therefore evolutionarily labile and not predictable from
taxonomy.

## Files
- `its_all.fasta` — concatenated raw sequences (relabelled)
- `trim_noSc.fasta` — final trimmed alignment (8 tips × 390 cols)
- `tree_noSc_rooted.nwk` — final rooted ML tree
- `FIG_phylo_traits.png` — tree + thermal traits figure
