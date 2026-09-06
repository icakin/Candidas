#!/usr/bin/env bash
# =============================================================================
# 02_phylogenomic_tree.sh -- the ten-taxon phylogenomic tree Figure 3 is drawn on
# =============================================================================
#   bash phylo/02_phylogenomic_tree.sh          (~30-60 min on 8 cores)
#
# Procedure (the one recorded in notes/phylogenomics_results.md, from which this script
# was written down on 2026-09-06; the original run was interactive):
#   1. Reciprocal best hits (DIAMOND blastp, e < 1e-10, sensitive) of every proteome
#      against C. auris clade I. An auris protein is a single-copy ortholog when it has an
#      RBH partner in EVERY other genome -> 2,931 orthologs in the original run.
#   2. A random sample of 600 of them (seed 42).
#   3. Each ortholog group aligned with MAFFT --auto; columns with > 20% gaps removed;
#      alignments concatenated -> supermatrix (289,787 aa x 10 taxa originally).
#   4. FastTree -lg -gamma (SH-like support) -> unrooted tree.
#   5. Rooted on the Debaryomycetaceae side (C. albicans + C. parapsilosis + C. lusitaniae),
#      pseudohaemulonii and lusitaniae dropped for Fig 3's four-taxon panel (pg_noCl.nwk is
#      the ten-taxon tree without C. lusitaniae).
#
# Inputs   phylo/proteomes/<genome>.faa (01_fetch_proteomes.py)
# Outputs  phylo/trees/core_orthologs.tsv     auris id -> ortholog id per genome
#          phylo/trees/supermatrix.faa        concatenated trimmed alignment
#          phylo/trees/pg_unrooted.nwk        FastTree output
#          phylo/trees/pg_rooted.nwk          rooted (the file 18_fig3.R reads)
#
# The committed trees/pg_rooted.nwk and trees/supermatrix.faa.gz are the ORIGINAL run's.
# The 600-gene sample and the ortholog map of that run were not saved, so a re-run with
# this script gives the same topology (every node had 1.00 support) but not byte-identical
# branch lengths. Keep the committed tree as the reference; use this to check topology.
#
# Needs: diamond, mafft, FastTree (or fasttree), python3 with biopython.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P="$HERE/proteomes"; T="$HERE/trees"; W="$T/work"; mkdir -p "$W"
CPU="${CPU:-8}"; N_GENES="${N_GENES:-600}"; SEED="${SEED:-42}"
REF=auris_cladeI
GENOMES=(auris_cladeI auris_cladeII auris_cladeIII auris_cladeIV haemulonii duobushaemulonii pseudohaemulonii parapsilosis albicans lusitaniae)
FT=$(command -v FastTree || command -v fasttree || { echo "!! FastTree not in PATH"; exit 1; })
command -v mafft >/dev/null || { echo "!! mafft not in PATH"; exit 1; }

# 1. reciprocal best hits against the reference
for g in "${GENOMES[@]}"; do diamond makedb --in "$P/$g.faa" -d "$W/$g" --quiet; done
for g in "${GENOMES[@]}"; do
  [ "$g" = "$REF" ] && continue
  diamond blastp -q "$P/$REF.faa" -d "$W/$g"  -o "$W/${REF}_vs_$g.tsv" --sensitive -p "$CPU" --evalue 1e-10 --max-target-seqs 1 --quiet --outfmt 6 qseqid sseqid bitscore
  diamond blastp -q "$P/$g.faa"   -d "$W/$REF" -o "$W/${g}_vs_$REF.tsv" --sensitive -p "$CPU" --evalue 1e-10 --max-target-seqs 1 --quiet --outfmt 6 qseqid sseqid bitscore
done

# 2-3. single-copy orthologs, sample, align, trim, concatenate
python3 - "$W" "$P" "$T" "$REF" "$N_GENES" "$SEED" "${GENOMES[@]}" <<'PY'
import sys, random, subprocess, collections
from pathlib import Path
from Bio import SeqIO
W, P, T, REF, N, SEED = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
G = sys.argv[7:]
def best(f):
    d = {}
    for l in open(f):
        q, s, b = l.split('\t'); b = float(b)
        if q not in d or b > d[q][1]: d[q] = (s, b)
    return {k: v[0] for k, v in d.items()}
orth = collections.defaultdict(dict)                     # auris id -> {genome: id}
for g in G:
    if g == REF: continue
    ab, ba = best(W / f'{REF}_vs_{g}.tsv'), best(W / f'{g}_vs_{REF}.tsv')
    for a, b in ab.items():
        if ba.get(b) == a: orth[a][g] = b
core = sorted(a for a, m in orth.items() if len(m) == len(G) - 1)
print(f'{len(core)} single-copy orthologs present in all {len(G)} genomes')
with open(T / 'core_orthologs.tsv', 'w') as fh:
    fh.write('\t'.join(G) + '\n')
    for a in core: fh.write('\t'.join([a] + [orth[a][g] for g in G if g != REF]) + '\n')
random.seed(SEED); sample = sorted(random.sample(core, min(N, len(core))))
seqs = {g: {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(str(P / f'{g}.faa'), 'fasta')} for g in G}
concat = {g: [] for g in G}
for i, a in enumerate(sample):
    fa = W / f'og{i:04d}.faa'; aln = W / f'og{i:04d}.aln'
    with open(fa, 'w') as fh:
        for g in G:
            gid = a if g == REF else orth[a][g]
            fh.write(f'>{g}\n{seqs[g][gid]}\n')
    with open(aln, 'w') as out:
        subprocess.run(['mafft', '--auto', '--quiet', str(fa)], stdout=out, check=True)
    al = {r.id: str(r.seq) for r in SeqIO.parse(str(aln), 'fasta')}
    L = len(next(iter(al.values())))
    keep = [j for j in range(L) if sum(al[g][j] == '-' for g in G) / len(G) <= 0.20]
    for g in G: concat[g].append(''.join(al[g][j] for j in keep))
with open(T / 'supermatrix.faa', 'w') as fh:
    for g in G: fh.write(f'>{g}\n{"".join(concat[g])}\n')
print(f'supermatrix: {sum(len(x) for x in concat[G[0]])} aa x {len(G)} taxa from {len(sample)} genes')
PY

# 4. tree
"$FT" -lg -gamma < "$T/supermatrix.faa" > "$T/pg_unrooted.nwk"

# 5. root on the Debaryomycetaceae / Clavispora side
python3 - "$T" <<'PY'
import sys
from pathlib import Path
T = Path(sys.argv[1])
try:
    from Bio import Phylo
except ImportError:
    sys.exit('biopython needed for rooting')
t = Phylo.read(str(T / 'pg_unrooted.nwk'), 'newick')
out = [c for c in t.get_terminals() if c.name in ('albicans', 'parapsilosis', 'lusitaniae')]
t.root_with_outgroup(*out)
Phylo.write(t, str(T / 'pg_rooted.nwk'), 'newick')
print('rooted tree -> trees/pg_rooted.nwk')
PY
echo "compare topology with the committed trees/pg_rooted.nwk before replacing it"
