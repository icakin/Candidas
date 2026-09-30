#!/usr/bin/env python3
"""03_thermal_census.py -- copy number of the canonical thermal-tolerance machinery in each proteome.

Nine gene families that the fungal heat-stress literature names first: HSP90, HSP70,
small HSPs, the trehalose pathway, calcineurin, the OLE1 fatty-acid desaturase, ergosterol
biosynthesis (ERG3/ERG11), alternative oxidase (AOX) and the alternative NADH
dehydrogenases. Query proteins (52) were taken from the two best-annotated proteomes,
C. auris B8441 and C. albicans SC5314, by product name, and are in inputs/thermal_queries.faa
with headers <family>__<source>_<accession>.

Each query is searched against all ten proteomes with DIAMOND blastp (sensitive mode,
e < 1e-20) and a subject protein counts as a family member when some query hits it at
>= 40% identity over >= 60% of the query's length. The census is the number of distinct
subject proteins per family per genome. Homology, not annotation text, so the sparsely
named clade IV GenBank annotation is not penalised.

Output  phylo/tables/thermal_census.csv     (rows: genomes; columns: families)
Run:    python3 phylo/03_thermal_census.py   (~1 min; needs diamond in PATH)
The committed table is reproduced exactly by this script (checked 2026-09-06).
"""
import subprocess, tempfile, glob, os
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
PROT = HERE / 'proteomes'; Q = HERE / 'inputs' / 'thermal_queries.faa'; OUT = HERE / 'tables' / 'thermal_census.csv'
GENOMES = ['auris_cladeI', 'auris_cladeII', 'auris_cladeIII', 'auris_cladeIV', 'haemulonii',
           'duobushaemulonii', 'pseudohaemulonii', 'parapsilosis', 'albicans', 'lusitaniae']
FAMS = ['HSP90', 'HSP70', 'sHSP', 'Trehalose', 'Calcineurin', 'Desaturase', 'Ergosterol', 'AOX', 'altNADH']
MIN_PID, MIN_QCOV, EVALUE = 40, 60, '1e-20'

def main():
    rows = {}
    with tempfile.TemporaryDirectory() as td:
        for g in GENOMES:
            db = f'{td}/{g}'; out = f'{td}/{g}.tsv'
            subprocess.run(['diamond', 'makedb', '--in', str(PROT / f'{g}.faa'), '-d', db, '--quiet'], check=True)
            subprocess.run(['diamond', 'blastp', '-q', str(Q), '-d', db, '-o', out, '--quiet', '--sensitive',
                            '--max-target-seqs', '50', '--evalue', EVALUE,
                            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qcovhsp'], check=True)
            h = pd.read_csv(out, sep='\t', names=['q', 's', 'pid', 'qcov'])
            h = h[(h.pid >= MIN_PID) & (h.qcov >= MIN_QCOV)]
            h['fam'] = h.q.str.split('__').str[0]
            rows[g] = {f: int(h[h.fam == f].s.nunique()) for f in FAMS}
    df = pd.DataFrame(rows).T.loc[GENOMES, FAMS]
    OUT.parent.mkdir(exist_ok=True); df.to_csv(OUT)
    print(df.to_string()); print('wrote', OUT)

if __name__ == '__main__':
    main()
