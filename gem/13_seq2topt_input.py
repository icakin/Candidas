#!/usr/bin/env python3
"""13_seq2topt_input.py -- the sequence table Seq2Topt and Seq2Tm predict from.

One row per enzyme in any of the four models: id = species|gene, sequence. Species in the
order auris, haemulonii, duobushaemulonii, parapsilosis; genes sorted within species.
Sequences longer than 2000 residues are left out (18 of 3060: the ESM2 embedding of a
2000-residue protein in a padded batch of four is the memory ceiling on a laptop, and
none of the 18 is a best-kcat enzyme for any reaction, so nothing in the pool constraint
depends on them; they take the species-median Topt and Tm downstream).

Output  gem/tables/seq2topt_input.csv        3042 rows (704 / 679 / 662 / 997 per species)
Run:    python3 gem/13_seq2topt_input.py     then bash gem/14_run_seq2topt_seq2tm.sh
"""
import pandas as pd
import cobra
from gempaths import *  # MODELS, TABLES, SP, load_proteome

MAX_LEN = 2000

def main():
    rows = []
    for sp, (xml, _) in SP.items():
        model = cobra.io.read_sbml_model(str(MODELS / xml))
        seqs = load_proteome(sp)
        kept = dropped = 0
        for g in sorted(x.id for x in model.genes):
            s = seqs.get(g)
            if s is None:
                continue
            if len(s) > MAX_LEN:
                dropped += 1; continue
            rows.append((f'{sp}|{g}', s)); kept += 1
        print(f'{sp:18s} {kept:4d} sequences kept, {dropped} longer than {MAX_LEN} aa dropped')
    out = pd.DataFrame(rows, columns=['id', 'sequence'])
    out.to_csv(TABLES / 'seq2topt_input.csv', index=False)
    print(f'{len(out)} rows -> tables/seq2topt_input.csv')

if __name__ == '__main__':
    main()
