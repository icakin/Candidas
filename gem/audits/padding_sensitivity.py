#!/usr/bin/env python3
"""padding_sensitivity.py -- does Seq2Tm's padding artefact bias Figure 4B?

THE PROBLEM
Upstream code/seq2tm.py batches sequences four at a time in file order and passes the
padded batch to MultiAttModel with no attention mask, so PAD tokens enter the attention
and pooling. A protein's predicted Tm therefore depends on which proteins happen to sit
beside it in the input file. Measured on 200 C. auris enzymes, the same protein moves by
up to 5.36 C (sd 1.27 C) between batch=4 and batch=1. gem/thermal_tm.csv was produced at
batch=4, and is reproduced exactly at that setting (max difference 0.0000 C, r = 1.000000).

WHY THAT MIGHT NOT MATTER
Figure 4B's headline is a MEAN paired difference over 1041 ortholog pairs. Unbiased noise
of sd 1.27 C contributes 1.27/sqrt(1041) ~ 0.04 C to the standard error of that mean,
well inside the reported interval, and paired proteins sit in different files so their
noise is independent.

WHY IT MIGHT
If the four proteome files have different length distributions they have systematically
different padding, and the artefact becomes a PER-SPECIES offset rather than random
noise. That is precisely the axis the comparison runs along. This script settles it by
recomputing the paired differences from batch=1 predictions, where there is no padding at
all, and putting them next to the published values.

The pairing is the one Figure 4B uses: two enzymes are orthologs if they catalyse the
same reaction in their species' model (kcat_reaction_<sp>.csv, best_gene per reaction).

    python3 padding_sensitivity.py --gem gem --nopad gem/thermal_tm_batch1.csv

WHERE thermal_tm_batch1.csv COMES FROM
    python3 gem/15_run_seq2tm.py --seqs-from gem/thermal_tm.csv \
            gem/thermal_tm_batch1.csv --batch-size 1
"""
import argparse
import numpy as np, pandas as pd
from scipy import stats as st
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

SP4 = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']
RELS = SP4[1:]


def paired(vals, kc):
    """vals: {species: {gene: Tm}} -> per-relative paired stats against C. auris."""
    out = {}
    for r in RELS:
        m = kc['auris'].merge(kc[r], on='reaction', suffixes=('_a', '_b'))
        d = np.array([vals['auris'][x] - vals[r][y]
                      for x, y in zip(m.best_gene_a, m.best_gene_b)
                      if x in vals['auris'] and y in vals[r]])
        se = float(d.std(ddof=1) / np.sqrt(len(d)))
        out[r] = dict(n=len(d), mean=float(d.mean()), se=se,
                      lo=float(d.mean() - 1.96 * se), hi=float(d.mean() + 1.96 * se),
                      p=float(st.ttest_1samp(d, 0)[1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gem', default=str(TABLES), help='directory holding kcat_reaction_*.csv and thermal_tm.csv')
    ap.add_argument('--nopad', required=True, help='batch=1 predictions, id,pred_tm')
    a = ap.parse_args()

    kc = {s: pd.read_csv(f'{a.gem}/kcat_reaction_{s}.csv')[['reaction', 'best_gene']].dropna()
          for s in SP4}

    def load(path, col):
        d = pd.read_csv(path)
        d['sp'] = d.id.str.split('|').str[0]
        d['g'] = d.id.str.split('|').str[1]
        return {s: dict(zip(d[d.sp == s].g, d[d.sp == s][col])) for s in SP4}

    pub = load(f'{a.gem}/thermal_tm.csv', 'pred_tm')       # batch=4, the published values
    npd = load(a.nopad, 'pred_tm')                          # batch=1, no padding

    common = {s: set(pub[s]) & set(npd[s]) for s in SP4}
    print('per-species offset introduced by padding (batch=4 minus batch=1):')
    print(f'{"":<20}{"n":>7}{"mean":>9}{"sd":>8}{"median":>9}')
    off = {}
    for s in SP4:
        d = np.array([pub[s][g] - npd[s][g] for g in common[s]])
        off[s] = d.mean()
        print(f'{s:<20}{len(d):>7}{d.mean():>+9.3f}{d.std():>8.3f}{np.median(d):>+9.3f}')
    spread = max(off.values()) - min(off.values())
    print(f'\nspread of the per-species offsets: {spread:.3f} C')
    print('  This is the quantity that matters. It enters the paired difference directly,')
    print('  so compare it with the 0.52 C the figure reports.')

    print(f'\n{"paired mean Tm difference (C. auris minus relative)":<52}')
    print(f'{"":<20}{"published":>22}{"no padding":>22}{"shift":>9}')
    P, N = paired(pub, kc), paired(npd, kc)
    for r in RELS:
        p, n = P[r], N[r]
        print(f'{r:<20}'
              f'{p["mean"]:>8.3f} [{p["lo"]:.2f},{p["hi"]:.2f}]'
              f'{n["mean"]:>10.3f} [{n["lo"]:.2f},{n["hi"]:.2f}]'
              f'{n["mean"]-p["mean"]:>+9.3f}')
    bp = max(P, key=lambda r: P[r]['mean']); bn = max(N, key=lambda r: N[r]['mean'])
    print(f'\nlargest paired difference, published : {P[bp]["mean"]:.3f} C ({bp}, n={P[bp]["n"]})')
    print(f'largest paired difference, no padding: {N[bn]["mean"]:.3f} C ({bn}, n={N[bn]["n"]})')
    for req, lab in ((32.54, 'fitted formulation'), (18.0, 'corrected unfolding width')):
        print(f'  fold gap against {req:.1f} C ({lab}): '
              f'published {req/P[bp]["mean"]:.0f}x, no padding {req/N[bn]["mean"]:.0f}x')

    print('\nREADING')
    print('If the no-padding difference is close to the published one, the artefact is')
    print('unbiased noise, Figure 4B stands as published, and the methods gain a sentence.')
    print('If it moves materially, the published 0.52 C was partly a padding artefact and')
    print('the panel must be rebuilt on batch=1 predictions before submission.')


if __name__ == '__main__':
    main()
