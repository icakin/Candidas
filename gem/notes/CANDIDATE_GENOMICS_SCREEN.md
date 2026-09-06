# Quick comparative-genomics screen — candidate thermal pathways

**Question:** across the tested Candida species, do coded-gene content or copy numbers of
canonical thermal / stress pathways differ in a way that tracks high-temperature growth?

**Method (quick, sequence-based).** Because haemulonii / duobushaemulonii /
pseudohaemulonii proteomes are largely annotated as "uncharacterized protein," a
keyword search on product names would falsely score genes as absent in exactly those
relatives. So candidate genes were defined from the well-annotated *C. albicans*
proteome and searched by **sequence** (DIAMOND blastp, pident ≥30 %, ≥50 % query and
subject coverage, E ≤ 1e-20) against every proteome; each protein was assigned to its
best-hitting family and counted → per-species family size (copy-number proxy).
Genomes available: species/clade **reference proteomes** only (auris clades I–IV, +
haemulonii, duobushaemulonii, pseudohaemulonii, parapsilosis, lusitaniae, albicans).
No per-isolate genomes, no raw reads.

## Result: no family tracks the thermal phenotype

Thermal ranking (Fig 3): *C. auris* (all 12 isolates, ≥44 °C) ≫ *C. haemulonii* /
*C. parapsilosis* (≈40 °C, mixed) > *C. duobushaemulonii* (≤38 °C).

Copy numbers of the core heat / stress machinery are **essentially conserved** across
all species — HSP90 (3–4), HSP70 (7–8), calcineurin, adenylate cyclase, PKA, TPS2,
trehalase, ergosterol pathway (10–11), FKS, PKC1, SOD, catalase, glutaredoxin, ferric
reductase. None separates thermotolerant from heat-sensitive species. Small ±1
differences are within assembly/annotation noise.

The one marked expansion — **siderophore transporters** (auris 7–12, haemulonii 14,
duobushaemulonii 9, pseudohaemulonii 11 vs parapsilosis 5, lusitaniae 1, albicans 1) —
is a **Candidozyma-clade** feature shared by the heat-sensitive *C. duobushaemulonii*
as much as by *C. auris*, so it tracks **phylogeny, not thermotolerance**. It is also
already reported (Muñoz et al. 2018). It does not explain the phenotype.

## Interpretation

Consistent with the rest of the project: the thermal difference is **not** in gene
content or copy number of the canonical pathways. That points to a **regulatory /
sequence-level** basis (promoter/heat-shock-element architecture, fixed coding
substitutions, or pathway deployment), not gene gain/loss — matching the etcGEM result
(enzyme Topt/Tm were clade-conserved) and ChatGPT's prediction.

## Caveats (why this is a screen, not a verdict)

- Copy number here is inferred from proteome gene counts — the weakest evidence type;
  a true test needs read-depth CNV (no reads available).
- Family counts merge paralogs; this is family size, not per-locus resolution.
- Keyword-defined queries missed a few families in albicans annotation (TPS1, HSP12 rows
  uninformative); the **key regulators HSF1 and CRZ1** are generically annotated and were
  not assessed here — they are the main thing a fuller OrthoFinder + sequence-level pass
  should add.
- Single thermotolerant lineage (auris) → any correlate is descriptive, not causal.

## Follow-up: HSF1 / CRZ1 regulator sequence check

Seeded with *S. cerevisiae* Hsf1 (P10961) and Crz1 (P53968), orthologs pulled from the
precomputed 1:1 ortholog table and aligned (MAFFT).

**HSF1 — clean negative.** Single copy in all 10 species. The HSF DNA-binding domain is
**invariant across the four auris clades (100 % identical)**, ~90–95 % vs close relatives
(haemulonii/duobushaemulonii/pseudohaemulonii), 74–82 % vs parapsilosis/lusitaniae/
albicans (ordinary phylogenetic divergence). **No auris-specific substitution at any
conserved DBD residue**; the two auris-specific residues within the DBD window fall at
positions that are already variable among the relatives themselves (fast-evolving loop,
not a targeted lesion). The 113 whole-protein auris-fixed substitutions are single-lineage
drift in HSF1's disordered regions and cannot be tied to the phenotype (one transition).

**CRZ1 — not assignable in this quick pass.** Fast-evolving; best hit not in the 1:1
ortholog set; zinc-finger DBD shared with many C2H2 factors. A trustworthy ortholog call
needs the fuller OrthoFinder pipeline.

**Promoter / heat-shock-element architecture** across species could not be tested: only
auris clade genomes are present locally; the relatives' assemblies + coordinates are not,
and would need downloading.

## Bottom line
A gene-content / copy-number explanation for *C. auris* thermotolerance is **not**
supported by the genomes we have. This is a clean negative that strengthens the
"phenotype–phylogeny discordance + regulation, not gene content" framing, and argues for
keeping comparative genomics minimal/supplementary unless HSF1/CRZ1 sequence-level
analysis turns up a fixed lesion.
