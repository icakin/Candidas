#!/usr/bin/env python3
"""lineage_specific.py — what is present in C. auris and absent from every relative?

The paper's negative result ("not metabolic enzyme thermostability") invites the obvious
question: then what? This is an unbiased screen for the genetic difference, using proteomes
already in the repository and no new data.

Design. A protein is called lineage-specific when it is present in ALL FOUR C. auris clades
(I-IV) and has no detectable homolog in ANY of C. haemulonii, C. duobushaemulonii,
C. pseudohaemulonii or C. parapsilosis. Requiring all four clades removes strain-level
assembly and annotation artefacts, which is the usual way a screen like this produces
nothing but noise. Requiring absence from four relatives, including the two closest, makes
the complement specific to the thermotolerant lineage rather than to one comparison.

Homology by DIAMOND blastp, which is what makes the call sequence-based rather than
annotation-based. That matters here: the C. auris proteome is better annotated than the
relatives' (27% hypothetical against 39%), so ANY count keyed on header text will find more
of everything in C. auris. That bias is precisely what this avoids.

The output is a gene list, not a mechanism. It says where to look next, and it is honest
about the fact that copy-number presence/absence is a weak proxy for a phenotype that is
much more likely to be regulatory.
"""
import subprocess, sys, os, re
from pathlib import Path
from collections import defaultdict
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

P = PROTEOMES
WORK = TABLES / 'lineage_work'; WORK.mkdir(exist_ok=True)   # DIAMOND dbs and hits, gitignored
AURIS = ['auris_cladeI', 'auris_cladeII', 'auris_cladeIII', 'auris_cladeIV']
RELS = ['haemulonii', 'duobushaemulonii', 'pseudohaemulonii', 'parapsilosis']
EVALUE = '1e-5'          # generous: we want ABSENCE to mean absence, so hits are easy to get
THREADS = str(os.cpu_count() or 2)


def headers(f):
    d = {}
    for line in open(P / f'{f}.faa'):
        if line.startswith('>'):
            acc = line[1:].split()[0]
            desc = line[1:].split(' ', 1)[1].strip() if ' ' in line else ''
            d[acc] = re.sub(r'\s*\[[^\]]*\]\s*$', '', desc)
    return d


def seqs(f):
    d, acc, buf = {}, None, []
    for line in open(P / f'{f}.faa'):
        if line.startswith('>'):
            if acc: d[acc] = ''.join(buf)
            acc = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if acc: d[acc] = ''.join(buf)
    return d


print('building DIAMOND databases for the relatives', flush=True)
for r in RELS:
    db = WORK / f'{r}.dmnd'
    if not db.exists():
        subprocess.run(['diamond', 'makedb', '--in', str(P / f'{r}.faa'),
                        '-d', str(WORK / r), '--quiet'], check=True)

# For each auris clade, which of its proteins hit NO relative?
orphan_by_clade = {}
for a in AURIS:
    hit = set()
    for r in RELS:
        out = WORK / f'{a}__{r}.tsv'
        if not out.exists():
            subprocess.run(['diamond', 'blastp', '-q', str(P / f'{a}.faa'),
                            '-d', str(WORK / r), '-o', str(out), '-e', EVALUE,
                            '--max-target-seqs', '1', '--outfmt', '6', 'qseqid',
                            '--threads', THREADS, '--quiet'], check=True)
        hit |= {l.split('\t')[0].strip() for l in open(out) if l.strip()}
    allp = set(headers(a))
    orphan_by_clade[a] = allp - hit
    print(f'  {a:<18} {len(allp)} proteins, {len(allp)-len(hit & allp)} with no relative homolog',
          flush=True)

# Conserved across all four clades: reciprocal presence within C. auris.
# Take clade I orphans and require each to have a homolog in the orphan set of the
# other three clades -- i.e. the gene is really in the lineage, not just in one assembly.
print('\nrequiring the orphan to be present in all four C. auris clades', flush=True)
c1 = orphan_by_clade['auris_cladeI']
s1 = seqs('auris_cladeI')
qf = WORK / 'c1_orphans.faa'
with open(qf, 'w') as fh:
    for acc in sorted(c1):
        fh.write(f'>{acc}\n{s1[acc]}\n')
keep = set(c1)
for a in AURIS[1:]:
    sub = WORK / f'{a}_orphans.faa'
    sa = seqs(a)
    with open(sub, 'w') as fh:
        for acc in sorted(orphan_by_clade[a]):
            fh.write(f'>{acc}\n{sa[acc]}\n')
    subprocess.run(['diamond', 'makedb', '--in', str(sub), '-d', str(WORK / f'{a}_orph'),
                    '--quiet'], check=True)
    out = WORK / f'c1__{a}_orph.tsv'
    subprocess.run(['diamond', 'blastp', '-q', str(qf), '-d', str(WORK / f'{a}_orph'),
                    '-o', str(out), '-e', EVALUE, '--max-target-seqs', '1',
                    '--outfmt', '6', 'qseqid', '--threads', THREADS, '--quiet'], check=True)
    keep &= {l.split('\t')[0].strip() for l in open(out) if l.strip()}
    print(f'  after {a}: {len(keep)}', flush=True)

h1 = headers('auris_cladeI')
print(f'\n{len(keep)} proteins present in all four C. auris clades and absent from all '
      f'{len(RELS)} relatives\n')
named = [(a, h1[a]) for a in sorted(keep)
         if not re.search(r'hypothetical|uncharacterized', h1[a], re.I)]
unnamed = len(keep) - len(named)
print(f'  {len(named)} have a functional annotation; {unnamed} are hypothetical/uncharacterized\n')
for a, d in named:
    print(f'  {a:<18} {d}')

with open(TABLES / 'lineage_specific_auris.tsv', 'w') as fh:
    fh.write('accession\tdescription\n')
    for a in sorted(keep):
        fh.write(f'{a}\t{h1[a]}\n')
print('\nwrote', TABLES / 'lineage_specific_auris.tsv')
