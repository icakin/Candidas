#!/usr/bin/env python3
"""11_aggregate_kcat.py -- one kcat per reaction from DLKcat's per-(enzyme, substrate) predictions.

The sMOMENT constraint needs a single turnover number per reaction. DLKcat returns one
prediction per (enzyme, substrate) pair, so a reaction with two isozymes and three
substrates has up to six. The base model takes the MAXIMUM: the reaction is assumed to be
carried by its fastest enzyme acting on the substrate it handles best. This is the
optimistic choice and is deliberate: the model's job in Fig 4 is to show how much thermal
separation the enzyme layer can generate, and a generous kcat gives it the most room. The
median-over-substrates alternative is tested in gem/audits/sens_medium_kcat.py and does not
change the conclusion.

Input   gem/tables/kcat_<species>.tsv        DLKcat output (Substrate Name = gene|reaction|kegg)
Output  gem/tables/kcat_reaction_<species>.csv
            reaction, kcat_per_s, best_gene, best_substrate, n_substrates
        n_substrates counts the (gene, substrate) predictions the maximum was taken over.
Reactions with no successful prediction (DLKcat skips salts and fails on a few SMILES) get
no row; the modelling scripts then treat them as enzyme-free, i.e. unconstrained, which
again favours the model.

Run:  python3 gem/11_aggregate_kcat.py
"""
import sys
import pandas as pd
from gempaths import *  # TABLES, SP

def main():
    for sp in SP:
        f = TABLES / f'kcat_{sp}.tsv'
        if not f.exists():
            print(f'{sp:18s} {f.name} missing (run 10_run_dlkcat.sh)', file=sys.stderr); continue
        d = pd.read_csv(f, sep='\t')
        d = d.rename(columns={'Kcat value (1/s)': 'kcat'})
        d['kcat'] = pd.to_numeric(d['kcat'], errors='coerce')
        d = d.dropna(subset=['kcat'])
        parts = d['Substrate Name'].str.split('|', expand=True)
        d['gene'], d['reaction'], d['substrate_kegg'] = parts[0], parts[1], parts[2]
        best = d.sort_values('kcat', ascending=False).drop_duplicates('reaction')
        n = d.groupby('reaction').size().rename('n_substrates')
        out = (best.set_index('reaction')[['kcat', 'gene', 'substrate_kegg']]
                   .join(n).reset_index()
                   .rename(columns={'kcat': 'kcat_per_s', 'gene': 'best_gene', 'substrate_kegg': 'best_substrate'}))
        out['kcat_per_s'] = out['kcat_per_s'].round(4)
        # keep the model's reaction order, as in the committed tables
        order = {r: i for i, r in enumerate(pd.read_csv(TABLES / f'dlkcat_pairs_{sp}.tsv', sep='\t').reaction.unique())}
        out = out.sort_values('reaction', key=lambda c: c.map(order)).reset_index(drop=True)
        out.to_csv(TABLES / f'kcat_reaction_{sp}.csv', index=False)
        print(f'{sp:18s} {len(d):5d} predictions over {out.shape[0]:4d} reactions -> kcat_reaction_{sp}.csv  '
              f'(median kcat {out.kcat_per_s.median():.1f} /s)')

if __name__ == '__main__':
    main()
