#!/usr/bin/env python3
"""07_rbh_orthologs.py -- one-to-one orthologs between C. auris and each relative, with identity.

Reciprocal best hits (DIAMOND blastp, e < 1e-10, sensitive mode) between the C. auris B8441
proteome and each of the other three proteomes. A pair (a, b) is kept when b is a's best
hit and a is b's best hit. Percent identity of the a->b alignment is recorded.

These tables are used by
  * gem/audits/common_network.py   run every species' enzyme parameters on the identical
                                   auris reaction scaffold (the network-free control);
  * the notes' identity-stratified paired Tm comparison (orthologs > 92% identity vs < 79%).
They are NOT how the draft models are built (that is 06_build_drafts.py, by KO evidence),
and Fig 4's paired ortholog statistic pairs enzymes by shared reaction, not by these RBH
tables, so that the pairing is by function inside the model rather than by sequence search.

Output  gem/tables/rbh/pid_<species>.tsv     auris_gene, ortholog, pident   (no header)
Run:    python3 gem/07_rbh_orthologs.py       (~2 min; needs diamond in PATH)
"""
import subprocess, tempfile
from pathlib import Path
from gempaths import *  # PROTEOMES, PROTEOME, TABLES, load_proteome

AURIS = PROTEOMES / 'auris_cladeI.faa'
TARGETS = {'haemulonii': PROTEOMES / 'haemulonii.faa',
           'duobushaemulonii': PROTEOMES / 'duobushaemulonii.faa',
           'parapsilosis': None}       # written from the UniProt TSV below, keyed by CPAR2_ tag

def write_fasta(seqs, path):
    with open(path, 'w') as fh:
        for k, v in seqs.items():
            fh.write(f'>{k}\n{v}\n')

def best_hits(q, s, td, tag):
    db = td / f'{tag}.dmnd'; out = td / f'{tag}.tsv'
    subprocess.run(['diamond', 'makedb', '--in', str(s), '-d', str(db), '--quiet'], check=True)
    subprocess.run(['diamond', 'blastp', '-q', str(q), '-d', str(db), '-o', str(out), '--sensitive',
                    '--max-target-seqs', '1', '--evalue', '1e-10', '--quiet',
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'bitscore'], check=True)
    d = {}
    for line in open(out):
        qid, sid, pid, bs = line.split('\t'); bs = float(bs)
        if qid not in d or bs > d[qid][2]:
            d[qid] = (sid, float(pid), bs)
    return d

def main():
    out_dir = TABLES / 'rbh'; out_dir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for sp, faa in TARGETS.items():
            if faa is None:
                faa = td / f'{sp}.faa'; write_fasta(load_proteome(sp), faa)
            ab = best_hits(AURIS, faa, td, f'a2{sp}'); ba = best_hits(faa, AURIS, td, f'{sp}2a')
            pairs = [(a, b, pid) for a, (b, pid, _) in ab.items() if ba.get(b, (None,))[0] == a]
            with open(out_dir / f'pid_{sp}.tsv', 'w') as fh:
                for a, b, pid in sorted(pairs):
                    fh.write(f'{a}\t{b}\t{pid:.1f}\n')
            import statistics
            print(f'{sp:18s} {len(pairs):5d} reciprocal best hits, median identity {statistics.median(p for _,_,p in pairs):.1f}%')

if __name__ == '__main__':
    main()
