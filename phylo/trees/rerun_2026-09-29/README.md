# Re-run of 02_phylogenomic_tree.sh, 29 Sep 2026

Versions: DIAMOND 2.1.9, MAFFT v7.505 (2022/Apr/10), FastTree 2.1.11 (double precision), Biopython 1.88.
Command: `CPU=2 bash phylo/02_phylogenomic_tree.sh` (N_GENES=600, SEED=42), Ubuntu 24.04.

Result: 3,658 single-copy orthologs present in all 10 genomes (2,931 in the original run of 7 Sep 2026,
whose DIAMOND version was not recorded); supermatrix 321,158 aa x 10 taxa from 600 genes;
FastTree LogLk -2787769.5, 0/7 bad splits, every node SH-like support 1.00.

Comparison with the committed reference tree (phylo/trees/pg_rooted.nwk): identical topology
(all 7 non-trivial splits), all supports 1.00, patristic distances within 5%:
  cladeI-cladeII 0.0051 (ref 0.0048); cladeI-cladeIII 0.0035 (0.0028); cladeI-cladeIV 0.0134 (0.0124);
  cladeII-cladeIV 0.0131 (0.0131); cladeI-haemulonii 0.346 (0.330); cladeI-parapsilosis 1.198 (1.146).
The committed tree remains the reference used by the figures and the phylogenetic-signal analysis.
