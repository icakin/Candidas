#!/usr/bin/env python3
"""04_proteome_signatures.py -- whole-proteome amino-acid composition indices of thermostability.

Three indices that distinguish thermophilic from mesophilic proteomes across the tree of
life, computed over every residue of every protein in each of the ten proteomes:
  IVYWREL   fraction of Ile, Val, Tyr, Trp, Arg, Glu, Leu (Zeldovich et al. 2007)
  CvP       charged (DEKR) minus polar (NQST) fraction (Suhre & Claverie 2003)
  EK_QH     (Glu + Lys) / (Gln + His) (Farias & Bonato 2003)
plus the charged, polar, acidic and basic fractions and the protein count. Terminal '*' is
dropped; only the 20 standard residues enter the denominators.

Output  phylo/tables/proteome_signatures.csv
Run:    python3 phylo/04_proteome_signatures.py   (~20 s; needs biopython)
The committed table is reproduced to 1e-15 by this script (checked 2026-09-06). Whether
these indices track the measured phenotypes is tested in notes/phylogenomics_results.md
(they do not, once phylogeny is accounted for).
"""
from pathlib import Path
import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
GENOMES = ['auris_cladeI', 'auris_cladeII', 'auris_cladeIII', 'auris_cladeIV', 'haemulonii',
           'duobushaemulonii', 'pseudohaemulonii', 'parapsilosis', 'albicans', 'lusitaniae']
AA = 'ACDEFGHIKLMNPQRSTVWY'

def main():
    rows = {}
    for g in GENOMES:
        c = {a: 0 for a in AA}; n = 0
        for r in SeqIO.parse(str(HERE / 'proteomes' / f'{g}.faa'), 'fasta'):
            n += 1
            for a in str(r.seq).rstrip('*'):
                if a in c: c[a] += 1
        tot = sum(c.values()); frac = lambda s: sum(c[a] for a in s) / tot
        rows[g] = dict(IVYWREL=frac('IVYWREL'), charged=frac('DEKR'), polar=frac('NQST'),
                       EK_QH=frac('EK') / frac('QH'), acidic=frac('DE'), basic=frac('KR'),
                       n_prot=float(n), CvP=frac('DEKR') - frac('NQST'))
    df = pd.DataFrame(rows).T
    (HERE / 'tables').mkdir(exist_ok=True); df.to_csv(HERE / 'tables' / 'proteome_signatures.csv')
    print(df.round(4).to_string()); print('wrote phylo/tables/proteome_signatures.csv')

if __name__ == '__main__':
    main()
